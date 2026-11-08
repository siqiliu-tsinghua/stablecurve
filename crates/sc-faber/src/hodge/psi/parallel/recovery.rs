//! Recovery requests reuse the planner's immutable frontier. Key construction
//! runs on the storage owner; responses move into ordinary constant equations.

use super::State;
use crate::cache_source::{IntegralKey, LookupBatch, Registry};
use sc_num::Q;
use std::ops::Range;
use std::sync::Arc;

const LOOKUP_ROWS: usize = 512;

struct Frontier {
    states: Arc<Vec<State>>,
    range: Range<usize>,
}

impl LookupBatch for Frontier {
    fn descendants_only(&self) -> bool {
        true
    }

    fn len(&self) -> usize {
        self.range.len()
    }

    fn key(&self, index: usize) -> IntegralKey {
        assert!(index < self.len());
        let state = &self.states[self.range.start + index];
        IntegralKey::Descendant {
            g: state.g,
            d: state.d.clone(),
            ch: state.ch.clone(),
        }
    }
}

#[cfg(test)]
pub(super) fn frontier(
    registry: &Arc<Registry>,
    states: &Arc<Vec<State>>,
) -> std::vec::IntoIter<Option<Q>> {
    frontier_range(registry, states, 0..states.len())
}

pub(super) fn frontier_range(
    registry: &Arc<Registry>,
    states: &Arc<Vec<State>>,
    range: Range<usize>,
) -> std::vec::IntoIter<Option<Q>> {
    assert!(range.start <= range.end && range.end <= states.len());
    let mut restored = Vec::with_capacity(range.len());
    for start in (range.start..range.end).step_by(LOOKUP_ROWS) {
        if registry.check_error().is_err() {
            break;
        }
        let end = (start + LOOKUP_ROWS).min(range.end);
        let batch = Arc::new(Frontier {
            states: Arc::clone(states),
            range: start..end,
        });
        match registry.lookup_batch(batch) {
            Ok(values) if values.len() == end - start => restored.extend(values),
            Ok(_) => {
                registry.fail("cache lookup returned an invalid frontier length".into());
                break;
            }
            Err(error) => {
                // An I/O failure is visible on the owner. Compute the missing
                // states exactly instead of manufacturing zero or retrying I/O.
                registry.fail(error);
                break;
            }
        }
    }
    restored.resize_with(range.len(), || None);
    restored.into_iter()
}

pub(super) fn take(
    values: &mut Option<std::vec::IntoIter<Option<Q>>>,
    count: usize,
) -> Option<Vec<Option<Q>>> {
    values.as_mut().map(|values| {
        let result: Vec<_> = values.by_ref().take(count).collect();
        assert_eq!(result.len(), count);
        result
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache_source::{CompletedIntegral, Pulled, ResultSource, storage::Storage};
    use sc_num::{One, Zero};

    struct Rows(Vec<CompletedIntegral>);
    impl ResultSource for Rows {
        fn ready(&self) -> usize {
            self.0.len()
        }
        fn pull(&self, cursor: usize, _: usize) -> Option<Pulled> {
            Some(Pulled {
                next: self.0.len(),
                rows: self.0[cursor..].to_vec(),
            })
        }
    }

    #[test]
    fn frontier_range_preserves_complete_normalized_c_identity() {
        let states = Arc::new(vec![
            State {
                g: 0,
                d: vec![0, 0, 0],
                ch: vec![],
                rank: (0, 1),
                index: 41,
            },
            State {
                g: 2,
                d: vec![0, 0],
                ch: vec![(5, 1)],
                rank: (5, 4),
                index: 7,
            },
        ]);
        let request = Frontier {
            states,
            range: 1..2,
        };
        assert_eq!(request.len(), 1);
        assert_eq!(
            request.key(0),
            IntegralKey::Descendant {
                g: 2,
                d: vec![0, 0],
                ch: vec![(5, 1)]
            }
        );
    }

    #[test]
    fn job_partition_preserves_legal_zero_and_miss_as_different_states() {
        let mut values = Some(vec![Some(Q::ZERO), None, Some(Q::ONE)].into_iter());
        assert_eq!(take(&mut values, 1), Some(vec![Some(Q::ZERO)]));
        assert_eq!(take(&mut values, 2), Some(vec![None, Some(Q::ONE)]));
        assert!(values.unwrap().next().is_none());
        assert!(take(&mut None, 2).is_none());
    }

    #[test]
    fn failed_owner_uses_misses_and_retains_the_original_error() {
        let registry = Arc::new(Registry::new());
        registry.fail("simulated SQLite failure".into());
        let states = Arc::new(vec![State {
            g: 2,
            d: vec![4],
            ch: vec![],
            rank: (0, 3),
            index: 0,
        }]);
        assert_eq!(frontier(&registry, &states).collect::<Vec<_>>(), vec![None]);
        assert_eq!(
            registry.check_error(),
            Err("simulated SQLite failure".into())
        );
    }

    #[test]
    fn backend_frontier_read_distinguishes_c_value_zero_and_absent_row() {
        let registry = Arc::new(Registry::new());
        let storage = Storage::open(":memory:", Arc::clone(&registry)).unwrap();
        registry.bind_storage(storage.handle());
        registry.register(Arc::new(Rows(vec![
            CompletedIntegral {
                key: IntegralKey::Descendant {
                    g: 2,
                    d: vec![],
                    ch: vec![(1, 3)],
                },
                value: Q::from_signeds(3, 5),
            },
            CompletedIntegral {
                key: IntegralKey::Descendant {
                    g: 2,
                    d: vec![0, 0],
                    ch: vec![(5, 1)],
                },
                value: Q::ZERO,
            },
        ])));
        storage.flush().unwrap();
        let states = Arc::new(vec![
            State {
                g: 2,
                d: vec![],
                ch: vec![(1, 3)],
                rank: (3, 2),
                index: 0,
            },
            State {
                g: 2,
                d: vec![0, 0],
                ch: vec![(5, 1)],
                rank: (5, 4),
                index: 1,
            },
            State {
                g: 3,
                d: vec![2, 3, 4],
                ch: vec![],
                rank: (0, 7),
                index: 2,
            },
        ]);
        assert_eq!(
            frontier(&registry, &states).collect::<Vec<_>>(),
            vec![Some(Q::from_signeds(3, 5)), Some(Q::ZERO), None]
        );
        storage.close().unwrap();
    }

    #[test]
    fn historical_subdescendants_restore_without_recomputation_or_dirty_rewrite() {
        use super::super::{Descendents, integrate, with_engine};
        let roots = [(vec![], Q::ONE)];
        let poly = [(vec![(1, 6)], Q::ONE)];
        let oracle = Arc::new(Registry::new());
        let expected = integrate(
            3,
            &roots,
            &poly,
            &mut Descendents::with_registry(Some(&oracle)),
            2,
            Some(&oracle),
        );
        for pure in [true, false] {
            // Select an actual non-root state visited by the original planner,
            // keeping all markings and the original normalized C value.
            let cached = with_engine(Some(&oracle), |engine| {
                let values = engine.source.as_ref().unwrap().read();
                for (weight, row) in engine.index.iter().enumerate() {
                    if weight >= 6 {
                        continue;
                    }
                    for shard in row {
                        for (ch, states) in shard.read().unwrap().iter() {
                            if ch.is_empty() != pure {
                                continue;
                            }
                            if let Some((d, &id)) = states.first_key_value() {
                                let triple_g = d.iter().map(|&power| i64::from(power)).sum::<i64>()
                                    + weight as i64
                                    - d.len() as i64
                                    + 3;
                                let g = u32::try_from(triple_g / 3).unwrap();
                                return Some(CompletedIntegral {
                                    key: IntegralKey::Descendant {
                                        g,
                                        d: d.clone(),
                                        ch: ch.clone(),
                                    },
                                    value: values[id].clone(),
                                });
                            }
                        }
                    }
                }
                None
            })
            .expect("oracle visits both pure and C subdescendants");
            let registry = Arc::new(Registry::new());
            let storage = Storage::open(":memory:", Arc::clone(&registry)).unwrap();
            registry.bind_storage(storage.handle());
            registry.register(Arc::new(Rows(vec![cached])));
            storage.flush().unwrap();
            registry.set_history(true);
            let value = integrate(
                3,
                &roots,
                &poly,
                &mut Descendents::with_registry(Some(&registry)),
                2,
                Some(&registry),
            );
            assert_eq!(value, expected);
            with_engine(Some(&registry), |engine| {
                let source = engine.source.as_ref().unwrap();
                assert_eq!(source.ready() + 1, source.read().len());
            });
            storage.close().unwrap();
        }
    }

    #[test]
    fn one_thread_history_recovers_c_with_the_original_external_grr_factor() {
        let registry = Arc::new(Registry::new());
        let storage = Storage::open(":memory:", Arc::clone(&registry)).unwrap();
        registry.bind_storage(storage.handle());
        registry.register(Arc::new(Rows(vec![CompletedIntegral {
            key: IntegralKey::Descendant {
                g: 2,
                d: vec![],
                ch: vec![(1, 3)],
            },
            value: Q::from_signeds(3, 5),
        }])));
        storage.flush().unwrap();
        registry.set_history(true);
        let value = crate::cache_source::with_registry(Arc::clone(&registry), || {
            super::super::super::integrate_with_threads(2, 0, &[], &[(1, 3)], 1)
        });
        assert_eq!(value, Q::from_signeds(1, 2880));
        super::super::with_engine(Some(&registry), |engine| {
            let source = engine.source.as_ref().unwrap();
            assert_eq!(source.read().len(), 1);
            assert_eq!(source.ready(), 0);
        });
        storage.close().unwrap();
    }

    #[test]
    fn formula_only_disk_history_does_not_probe_the_rebuilt_dag() {
        struct Database(std::path::PathBuf);
        impl Drop for Database {
            fn drop(&mut self) {
                for suffix in ["-wal", "-shm", ""] {
                    let mut path = self.0.as_os_str().to_os_string();
                    path.push(suffix);
                    let _ = std::fs::remove_file(std::path::PathBuf::from(path));
                }
            }
        }
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "sc-faber-formula-history-{}-{stamp}.sqlite",
            std::process::id()
        ));
        drop(
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .unwrap(),
        );
        let database = Database(path);
        {
            let registry = Arc::new(Registry::new());
            let storage = Storage::open(database.0.to_str().unwrap(), registry.clone()).unwrap();
            registry.bind_storage(storage.handle());
            registry.register(Arc::new(Rows(vec![
                CompletedIntegral {
                    key: IntegralKey::Descendant {
                        g: 1,
                        d: vec![1],
                        ch: vec![],
                    },
                    value: Q::from_signeds(1, 24),
                },
                CompletedIntegral {
                    key: IntegralKey::Descendant {
                        g: 2,
                        d: vec![4],
                        ch: vec![],
                    },
                    value: Q::from_signeds(1, 1152),
                },
                CompletedIntegral {
                    key: IntegralKey::Descendant {
                        g: 2,
                        d: vec![0],
                        ch: vec![(1, 4)],
                    },
                    value: Q::ZERO,
                },
            ])));
            storage.close().unwrap();
        }
        let registry = Arc::new(Registry::new());
        let storage = Storage::open(database.0.to_str().unwrap(), registry.clone()).unwrap();
        registry.bind_storage(storage.handle());
        assert!(registry.has_history());
        assert!(!registry.has_descendant_history());
        let value = crate::cache_source::with_registry(registry.clone(), || {
            super::super::super::integrate_with_threads(3, 0, &[], &[(1, 6)], 4)
        });
        assert_eq!(value, Q::from_signeds(1, 90720));
        let report = storage.close().unwrap();
        assert_eq!(
            report.read_keys, 0,
            "formula history cannot serve planner states"
        );
        assert_eq!(report.read_hits, 0);
    }
}
