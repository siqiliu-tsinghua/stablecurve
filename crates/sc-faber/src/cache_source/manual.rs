//! Mathematical L1 for the explicit legacy cache_get/cache_put API.
//! Engine memos publish their own sources instead of copying values here.

use super::{CompletedIntegral, IntegralKey, Pulled, Registry, ResultSource, record_bytes};
use sc_num::Q;
use std::collections::BTreeMap;
use std::sync::RwLock;
use std::sync::atomic::{AtomicUsize, Ordering};

type Key = (String, Vec<u32>);

#[derive(Default)]
struct Values {
    values: BTreeMap<Key, Q>,
    completed: Vec<Key>,
    base_completion: usize,
}

impl Values {
    fn ready(&self) -> usize {
        self.base_completion
            .checked_add(self.completed.len())
            .expect("manual completion cursor exhausted")
    }
}

#[derive(Default)]
pub(super) struct ManualSource {
    values: RwLock<Values>,
    durable: AtomicUsize,
}

impl ManualSource {
    pub(super) fn get(&self, kind: &str, nums: &[u32]) -> Option<Q> {
        // Explicit legacy API calls are uncommon; the mathematical engines
        // retain their original, borrowed-key indexes in their own sources.
        self.values
            .read()
            .unwrap()
            .values
            .get(&(kind.into(), nums.to_vec()))
            .cloned()
    }

    pub(super) fn put(&self, kind: &str, nums: &[u32], value: &Q) {
        let key = (kind.to_owned(), nums.to_vec());
        let mut values = self.values.write().unwrap();
        if values.values.get(&key).is_some_and(|old| old == value) {
            return;
        }
        values.values.insert(key.clone(), value.clone());
        values.completed.push(key);
    }
}

impl ManualSource {
    fn reclaim_with_history(&self, registry: Option<&Registry>) -> usize {
        let retired = {
            let Ok(mut values) = self.values.try_write() else {
                return 0;
            };
            let ready = values.ready();
            if self.durable.load(Ordering::Acquire) != ready {
                return 0;
            }
            if let Some(registry) = registry.filter(|_| !values.values.is_empty()) {
                registry.set_history(true);
            }
            let retired = (
                std::mem::take(&mut values.values),
                std::mem::take(&mut values.completed),
            );
            values.base_completion = ready;
            retired
        };
        let key_heap = |key: &Key| {
            key.0
                .capacity()
                .saturating_add(key.1.capacity().saturating_mul(std::mem::size_of::<u32>()))
        };
        let mut bytes = retired
            .0
            .len()
            .saturating_mul(std::mem::size_of::<(Key, Q)>())
            .saturating_add(
                retired
                    .1
                    .capacity()
                    .saturating_mul(std::mem::size_of::<Key>()),
            );
        for (key, value) in &retired.0 {
            bytes = bytes
                .saturating_add(key_heap(key))
                .saturating_add(super::value_bytes(value).saturating_sub(std::mem::size_of::<Q>()));
        }
        for key in &retired.1 {
            bytes = bytes.saturating_add(key_heap(key));
        }
        drop(retired);
        bytes
    }
}

impl ResultSource for ManualSource {
    fn floor(&self) -> usize {
        self.values.read().unwrap().base_completion
    }

    fn ready(&self) -> usize {
        self.values.read().unwrap().ready()
    }

    fn durable(&self, cursor: usize) {
        self.durable.fetch_max(cursor, Ordering::Release);
    }

    fn reclaim(&self) -> usize {
        self.reclaim_with_history(None)
    }

    fn reclaim_during_query(&self, registry: &Registry) -> usize {
        self.reclaim_with_history(Some(registry))
    }

    fn pull(&self, cursor: usize, byte_limit: usize) -> Option<Pulled> {
        let values = self.values.try_read().ok()?;
        if cursor < values.base_completion {
            return None;
        }
        assert!(cursor <= values.ready(), "invalid manual completion cursor");
        let mut rows = Vec::new();
        let mut bytes = 0usize;
        for key in &values.completed[cursor - values.base_completion..] {
            let row = CompletedIntegral {
                key: IntegralKey::Legacy {
                    kind: key.0.clone(),
                    nums: key.1.clone(),
                },
                value: values.values[key].clone(),
            };
            bytes = bytes.saturating_add(record_bytes(&row));
            rows.push(row);
            if bytes >= byte_limit.max(1) {
                break;
            }
        }
        Some(Pulled {
            next: cursor + rows.len(),
            rows,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn active_reclaim_skips_dirty_or_borrowed_values_and_publishes_history() {
        let registry = std::sync::Arc::new(Registry::new());
        let source = ManualSource::default();
        source.put("test", &[1], &Q::from(7));
        super::super::with_registry(registry.clone(), || {
            assert!(registry.try_reclaim().is_none());
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            source.durable(1);
            let held = source.values.read().unwrap();
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            drop(held);
            let owned = source.get("test", &[1]).unwrap();
            assert!(source.reclaim_during_query(&registry) > 0);
            assert!(registry.has_history());
            assert_eq!(owned, Q::from(7));
            assert_eq!(source.floor(), 1);
            assert!(source.get("test", &[1]).is_none());
        });
    }

    #[test]
    fn an_old_commit_cannot_reclaim_a_new_manual_overwrite() {
        let source = ManualSource::default();
        source.put("example", &[1], &Q::from(1));
        source.durable(1);
        source.put("example", &[1], &Q::from(2));
        assert_eq!(source.reclaim(), 0);
        assert_eq!(source.get("example", &[1]), Some(Q::from(2)));
        source.durable(2);
        let held = source.values.read().unwrap();
        assert_eq!(source.reclaim(), 0);
        drop(held);
        assert!(source.reclaim() > 0);
        assert_eq!(source.floor(), 2);
        assert_eq!(source.ready(), 2);
        assert!(source.get("example", &[1]).is_none());
        assert!(source.pull(0, 1024).is_none());
        source.put("example", &[2], &Q::from(3));
        let next = source.pull(2, 1024).unwrap();
        assert_eq!(next.next, 3);
        assert_eq!(next.rows[0].value, Q::from(3));
    }
}
