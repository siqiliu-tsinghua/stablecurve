//! One process-level native memory monitor shared by persistent cache owners.
//!
//! SQL and identity work cannot delay this clock. The monitor acquires no
//! mathematical/source lock and owns no cache value. Its last lease interrupts
//! the timer and joins outside the singleton lock. An in-progress native sample
//! must return before shutdown; the timer itself is always interruptible.

use super::budget::process_budget;
use std::sync::{Arc, Condvar, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

const INTERVAL: Duration = Duration::from_millis(100);

struct Control {
    stop: Mutex<bool>,
    changed: Condvar,
}

struct Running {
    users: usize,
    control: Arc<Control>,
    thread: JoinHandle<()>,
}

type State = Mutex<Option<Running>>;
static SAMPLER: State = Mutex::new(None);

/// Every cache service holds one lease through SQL/reader shutdown.
/// Leases are moved, not cloned: each successful `start` adds one owner.
pub(crate) struct SamplerLease {
    state: &'static State,
    control: Arc<Control>,
}

pub(crate) fn start() -> Result<SamplerLease, String> {
    start_at(&SAMPLER, || {
        process_budget().refresh_if_due();
    })
}

fn start_at(state: &'static State, sample: fn()) -> Result<SamplerLease, String> {
    let mut current = state.lock().unwrap_or_else(|error| error.into_inner());
    if let Some(running) = current.as_mut() {
        running.users = running
            .users
            .checked_add(1)
            .ok_or("memory sampler lease count exhausted")?;
        return Ok(SamplerLease {
            state,
            control: Arc::clone(&running.control),
        });
    }
    let control = Arc::new(Control {
        stop: Mutex::new(false),
        changed: Condvar::new(),
    });
    let monitor_control = Arc::clone(&control);
    let monitor = thread::Builder::new()
        .name("sc-faber-memory".into())
        .spawn(move || monitor(monitor_control, sample))
        .map_err(|error| format!("cannot start memory sampler: {error}"))?;
    *current = Some(Running {
        users: 1,
        control: Arc::clone(&control),
        thread: monitor,
    });
    Ok(SamplerLease { state, control })
}

fn monitor(control: Arc<Control>, sample: fn()) {
    loop {
        if *control
            .stop
            .lock()
            .unwrap_or_else(|error| error.into_inner())
        {
            return;
        }
        sample();
        let stop = control
            .stop
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if *stop {
            return;
        }
        let (stop, _) = control
            .changed
            .wait_timeout_while(stop, INTERVAL, |stop| !*stop)
            .unwrap_or_else(|error| error.into_inner());
        if *stop {
            return;
        }
    }
}

impl Drop for SamplerLease {
    fn drop(&mut self) {
        let monitor = {
            let mut current = self.state.lock().unwrap_or_else(|error| error.into_inner());
            let Some(running) = current.as_mut() else {
                return;
            };
            if !Arc::ptr_eq(&running.control, &self.control) {
                return;
            }
            debug_assert!(running.users > 0);
            running.users -= 1;
            if running.users != 0 {
                return;
            }
            let running = current.take().unwrap();
            *running
                .control
                .stop
                .lock()
                .unwrap_or_else(|error| error.into_inner()) = true;
            running.control.changed.notify_all();
            running.thread
        };
        // A new owner may already have started a monitor. The budget's sampling
        // atomic prevents duplicate native captures during this short overlap.
        // Production monitoring owns no lease, so self-join cannot arise there.
        if monitor.thread().id() != thread::current().id() {
            let _ = monitor.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::panic::{AssertUnwindSafe, catch_unwind};

    fn idle() {}

    #[test]
    fn two_leases_share_one_monitor_and_last_drop_joins_it() {
        static STATE: State = Mutex::new(None);
        let first = start_at(&STATE, idle).unwrap();
        let second = start_at(&STATE, idle).unwrap();
        assert!(Arc::ptr_eq(&first.control, &second.control));
        assert_eq!(STATE.lock().unwrap().as_ref().unwrap().users, 2);
        let weak = Arc::downgrade(&first.control);
        drop(first);
        assert_eq!(STATE.lock().unwrap().as_ref().unwrap().users, 1);
        assert!(weak.upgrade().is_some());
        drop(second);
        assert!(STATE.lock().unwrap().is_none());
        assert!(
            weak.upgrade().is_none(),
            "last Drop must join and release the monitor"
        );
    }

    #[test]
    fn caller_unwind_stops_the_monitor_and_allows_a_fresh_start() {
        static STATE: State = Mutex::new(None);
        let mut stopped = None;
        assert!(
            catch_unwind(AssertUnwindSafe(|| {
                let lease = start_at(&STATE, idle).unwrap();
                stopped = Some(Arc::downgrade(&lease.control));
                panic!("simulate cache service setup failure");
            }))
            .is_err()
        );
        assert!(STATE.lock().unwrap().is_none());
        assert!(stopped.unwrap().upgrade().is_none());
        let restarted = start_at(&STATE, idle).unwrap();
        let weak = Arc::downgrade(&restarted.control);
        drop(restarted);
        assert!(weak.upgrade().is_none());
    }

    #[test]
    fn repeated_concurrent_start_and_stop_leaves_no_monitor_owner() {
        static STATE: State = Mutex::new(None);
        for _ in 0..8 {
            std::thread::scope(|scope| {
                for _ in 0..4 {
                    scope.spawn(|| {
                        let lease = start_at(&STATE, idle).unwrap();
                        let second = start_at(&STATE, idle).unwrap();
                        assert!(Arc::ptr_eq(&lease.control, &second.control));
                        drop(lease);
                        drop(second);
                    });
                }
            });
            assert!(STATE.lock().unwrap().is_none());
        }
    }
}
