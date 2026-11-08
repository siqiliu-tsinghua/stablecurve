//! Explicit C1 experiment: pull completed rows without encoding or SQL.
//!
//! This module is dormant unless a diagnostic caller starts a reader. Ordinary
//! mathematical calls neither create this thread nor wait for its completion.

use super::{Registry, ResultSource, record_bytes};
use std::hint::black_box;
use std::io;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

const COPY_WINDOW_BYTES: usize = 64 * 1024;
const IDLE_PAUSE: Duration = Duration::from_millis(1);

#[derive(Default, Debug)]
pub(crate) struct ReaderReport {
    pub(crate) sources: usize,
    pub(crate) pulls: usize,
    pub(crate) batches: usize,
    pub(crate) rows: usize,
    /// Effective copied row bytes, excluding allocation spare capacity.
    pub(crate) bytes: usize,
    pub(crate) busy: usize,
    pub(crate) empty: usize,
    pub(crate) reached_target: bool,
    pub(crate) invalid_cursor: bool,
    pub(crate) worker_panicked: bool,
}

#[derive(Default)]
struct Control {
    stop: AtomicBool,
    /// Fixed completion watermarks captured by explicit experimental finish.
    finish: Mutex<Option<Vec<usize>>>,
}

struct SourceCursor {
    source: Arc<dyn ResultSource>,
    cursor: usize,
}

/// One independently owned background worker, never installed by attach_cache.
pub(crate) struct Reader {
    registry: Arc<Registry>,
    control: Arc<Control>,
    worker: Option<JoinHandle<ReaderReport>>,
}

impl Reader {
    pub(super) fn start(registry: Arc<Registry>) -> io::Result<Self> {
        let control = Arc::new(Control::default());
        let worker_registry = Arc::clone(&registry);
        let worker_control = Arc::clone(&control);
        let worker = thread::Builder::new()
            .name("sc-faber-c1-reader".into())
            .spawn(move || run(worker_registry, worker_control))?;
        Ok(Self {
            registry,
            control,
            worker: Some(worker),
        })
    }

    /// Experimental teardown barrier. Capture already published completions,
    /// drain through that finite target, and join. This is not durability and
    /// must not be inserted into an ordinary mathematical return path.
    pub(crate) fn finish(mut self) -> ReaderReport {
        let targets = self
            .registry
            .sources()
            .iter()
            .map(|source| source.ready())
            .collect();
        *self.control.finish.lock().unwrap() = Some(targets);
        let worker = self.worker.take().unwrap();
        worker.thread().unpark();
        worker.join().unwrap_or_else(|_| ReaderReport {
            worker_panicked: true,
            ..ReaderReport::default()
        })
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        self.control.stop.store(true, Ordering::Release);
        if let Some(worker) = self.worker.take() {
            worker.thread().unpark();
            // Drop may run during unwinding; never introduce a second panic.
            let _ = worker.join();
        }
    }
}

fn run(registry: Arc<Registry>, control: Arc<Control>) -> ReaderReport {
    let mut report = ReaderReport::default();
    let mut sources: Vec<SourceCursor> = Vec::new();
    loop {
        if control.stop.load(Ordering::Acquire) {
            return report;
        }
        // Only new Arc handles are copied. No source lock is acquired while
        // the registry list lock is held.
        sources.extend(
            registry
                .sources_since(sources.len())
                .into_iter()
                .map(|source| SourceCursor {
                    cursor: source.floor(),
                    source,
                }),
        );
        report.sources = sources.len();
        let targets = control.finish.lock().unwrap().clone();
        if let Some(targets) = &targets {
            // A source may have been registered immediately after discovery
            // but before finish captured its target. Discover it next pass.
            if sources.len() >= targets.len()
                && sources
                    .iter()
                    .zip(targets)
                    .all(|(source, &target)| source.cursor >= target)
            {
                report.reached_target = true;
                return report;
            }
        }

        let mut progressed = false;
        for (index, source) in sources.iter_mut().enumerate() {
            if control.stop.load(Ordering::Acquire) {
                return report;
            }
            // Another durable consumer may have reclaimed an unread prefix.
            // This diagnostic reader cannot claim to have copied those rows,
            // and retrying a cursor below floor can never make progress.
            if source.cursor < source.source.floor() {
                report.invalid_cursor = true;
                return report;
            }
            let target = match &targets {
                Some(targets) => match targets.get(index) {
                    Some(&target) => target,
                    None => continue,
                },
                None => source.source.ready(),
            };
            // A source may prepare compact identities while planning has not
            // yet produced any Q values. Only a fixed finish target permits
            // skipping an already covered source; empty ordinary pulls are
            // temporary observations, never end-of-source markers.
            if targets.is_some() && source.cursor >= target {
                continue;
            }
            report.pulls = report.pulls.saturating_add(1);
            let Some(batch) = source.source.pull(source.cursor, COPY_WINDOW_BYTES) else {
                report.busy = report.busy.saturating_add(1);
                continue;
            };
            if source.cursor.checked_add(batch.rows.len()) != Some(batch.next) {
                report.invalid_cursor = true;
                return report;
            }
            if batch.rows.is_empty() {
                report.empty = report.empty.saturating_add(1);
            } else {
                report.batches = report.batches.saturating_add(1);
                report.rows = report.rows.saturating_add(batch.rows.len());
                for row in &batch.rows {
                    report.bytes = report.bytes.saturating_add(record_bytes(row));
                }
                // Observe the actual owned row/key/value batch; counters alone
                // must not permit optimizing away the copied mathematical data.
                black_box(&batch.rows);
            }
            progressed |= batch.next > source.cursor;
            source.cursor = batch.next;
            // Source pull has already released its locks. Drop this one copy
            // window before visiting another source; no all-query row buffer.
            drop(batch);
        }
        if progressed {
            thread::yield_now();
        } else {
            thread::park_timeout(IDLE_PAUSE);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache_source::{CompletedIntegral, IntegralKey, Pulled};
    use sc_num::Q;

    struct Rows {
        values: Vec<Q>,
    }

    struct Planning {
        prepared: Mutex<Option<std::sync::mpsc::Sender<()>>>,
    }

    impl ResultSource for Planning {
        fn ready(&self) -> usize {
            0
        }

        fn pull(&self, cursor: usize, _byte_limit: usize) -> Option<Pulled> {
            if let Some(prepared) = self.prepared.lock().unwrap().take() {
                prepared.send(()).unwrap();
            }
            Some(Pulled {
                next: cursor,
                rows: Vec::new(),
            })
        }
    }

    impl ResultSource for Rows {
        fn ready(&self) -> usize {
            self.values.len()
        }

        fn pull(&self, cursor: usize, byte_limit: usize) -> Option<Pulled> {
            let mut rows = Vec::new();
            let mut bytes = 0usize;
            let mut next = cursor;
            for (index, value) in self.values.iter().enumerate().skip(cursor) {
                let row = CompletedIntegral {
                    key: IntegralKey::Descendant {
                        g: 1,
                        d: vec![index as u32],
                        ch: Vec::new(),
                    },
                    value: value.clone(),
                };
                bytes = bytes.saturating_add(record_bytes(&row));
                rows.push(row);
                next = index + 1;
                if bytes >= byte_limit {
                    break;
                }
            }
            Some(Pulled { next, rows })
        }
    }

    #[test]
    fn finish_discovers_new_sources_and_drains_published_rows() {
        let registry = Arc::new(Registry::new());
        let reader = registry.start_reader().unwrap();
        registry.register(Arc::new(Rows {
            values: vec![Q::from(1), Q::from_signeds(1, 24)],
        }));
        registry.register(Arc::new(Rows {
            values: vec![Q::from_signeds(1, 1152)],
        }));
        let report = reader.finish();
        assert_eq!(report.sources, 2);
        assert_eq!(report.rows, 3);
        assert!(report.bytes >= 3 * std::mem::size_of::<CompletedIntegral>());
        assert!(report.reached_target);
        assert!(!report.invalid_cursor);
        assert!(!report.worker_panicked);
    }

    #[test]
    fn empty_registry_finishes_without_future_source_wait() {
        let registry = Arc::new(Registry::new());
        let report = registry.start_reader().unwrap().finish();
        assert_eq!(report.rows, 0);
        assert!(report.reached_target);
    }

    #[test]
    fn reclaimed_unread_prefix_ends_the_measurement_instead_of_waiting_forever() {
        use std::sync::atomic::AtomicUsize;
        struct Reclaimed {
            floor: AtomicUsize,
            observed: Mutex<Option<std::sync::mpsc::Sender<()>>>,
        }
        impl ResultSource for Reclaimed {
            fn ready(&self) -> usize {
                1
            }
            fn floor(&self) -> usize {
                self.floor.load(Ordering::Acquire)
            }
            fn pull(&self, _: usize, _: usize) -> Option<Pulled> {
                if let Some(observed) = self.observed.lock().unwrap().take() {
                    observed.send(()).unwrap();
                }
                None
            }
        }
        let registry = Arc::new(Registry::new());
        let (observed, received) = std::sync::mpsc::channel();
        let source = Arc::new(Reclaimed {
            floor: AtomicUsize::new(0),
            observed: Mutex::new(Some(observed)),
        });
        registry.register(source.clone());
        let reader = registry.start_reader().unwrap();
        received.recv_timeout(Duration::from_secs(5)).unwrap();
        source.floor.store(1, Ordering::Release);
        let report = reader.finish();
        assert!(report.invalid_cursor);
        assert!(!report.reached_target);
        assert_eq!(report.rows, 0);
    }

    #[test]
    fn planning_source_is_polled_before_any_values_are_ready() {
        let registry = Arc::new(Registry::new());
        let (prepared, received) = std::sync::mpsc::channel();
        registry.register(Arc::new(Planning {
            prepared: Mutex::new(Some(prepared)),
        }));
        let reader = registry.start_reader().unwrap();
        received.recv_timeout(Duration::from_secs(5)).unwrap();
        let report = reader.finish();
        assert!(report.pulls > 0);
        assert_eq!(report.rows, 0);
        assert!(report.reached_target);
    }

    #[test]
    fn unwind_stops_and_joins_the_reader() {
        let result = std::panic::catch_unwind(|| {
            let registry = Arc::new(Registry::new());
            let _reader = registry.start_reader().unwrap();
            panic!("reader scope test");
        });
        assert!(result.is_err());
    }
}
