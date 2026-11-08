//! Readable ownership for the existing ST/IC mathematical memo tables.
//!
//! The HashMap is the sole mathematical value table. Completion order stores
//! identities only; raw persistent keys and result copies are built by pull.

use super::{CompletedIntegral, IntegralKey, Pulled, Registry, ResultSource, record_bytes};
use sc_num::Q;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::hash::Hash;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{RwLock, TryLockError};

struct Store<K> {
    values: HashMap<K, Q>,
    completed: Vec<K>,
    base_completion: usize,
}

impl<K> Store<K> {
    fn ready(&self) -> usize {
        self.base_completion
            .checked_add(self.completed.len())
            .expect("legacy completion cursor exhausted")
    }
}

pub(crate) struct LegacySource<K> {
    kind: &'static str,
    nums: fn(&K) -> Vec<u32>,
    store: RwLock<Store<K>>,
    writers: AtomicUsize,
    durable: AtomicUsize,
}

struct WriteIntent<'a>(&'a AtomicUsize);

impl Drop for WriteIntent<'_> {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::Release);
    }
}

impl<K: Eq + Hash + Clone + Send + Sync> LegacySource<K> {
    pub(crate) fn new(kind: &'static str, nums: fn(&K) -> Vec<u32>) -> Self {
        Self {
            kind,
            nums,
            store: RwLock::new(Store {
                values: HashMap::new(),
                completed: Vec::new(),
                base_completion: 0,
            }),
            writers: AtomicUsize::new(0),
            durable: AtomicUsize::new(0),
        }
    }

    pub(crate) fn get(&self, key: &K) -> Option<Q> {
        self.store.read().unwrap().values.get(key).cloned()
    }

    /// Perform the original mathematical memo insertion. The single Q clone
    /// becomes the memo's value; the completion directory contains no Q.
    /// Duplicate keys retain their original value and completion position.
    pub(crate) fn insert(&self, key: K, value: &Q) -> bool {
        self.writers.fetch_add(1, Ordering::AcqRel);
        let _intent = WriteIntent(&self.writers);
        let mut store = self.store.write().unwrap();
        match store.values.entry(key) {
            Entry::Occupied(_) => false,
            Entry::Vacant(entry) => {
                let identity = entry.key().clone();
                entry.insert(value.clone());
                // Both structures are published when this write guard ends.
                store.completed.push(identity);
                true
            }
        }
    }
}

impl LegacySource<Vec<u32>> {
    pub(crate) fn get_slice(&self, key: &[u32]) -> Option<Q> {
        self.store.read().unwrap().values.get(key).cloned()
    }
}

impl<K: Eq + Hash + Clone + Send + Sync> LegacySource<K> {
    fn reclaim_with_history(&self, registry: Option<&Registry>) -> usize {
        let retired = {
            let Ok(mut store) = self.store.try_write() else {
                return 0;
            };
            let ready = store.ready();
            if self.durable.load(Ordering::Acquire) != ready {
                return 0;
            }
            if let Some(registry) = registry.filter(|_| !store.values.is_empty()) {
                registry.set_history(true);
            }
            let retired = (
                std::mem::take(&mut store.values),
                std::mem::take(&mut store.completed),
            );
            store.base_completion = ready;
            retired
        };
        // K is generic: count known table/identity capacity and Q heap only.
        // Unknown key-owned heaps, allocator metadata and retention are not
        // represented by this lower estimate; native sampling reports reality.
        let mut bytes = retired
            .0
            .capacity()
            .saturating_mul(std::mem::size_of::<(K, Q)>())
            .saturating_add(
                retired
                    .1
                    .capacity()
                    .saturating_mul(std::mem::size_of::<K>()),
            );
        for value in retired.0.values() {
            bytes = bytes
                .saturating_add(super::value_bytes(value).saturating_sub(std::mem::size_of::<Q>()));
        }
        drop(retired);
        bytes
    }
}

impl<K: Eq + Hash + Clone + Send + Sync> ResultSource for LegacySource<K> {
    fn floor(&self) -> usize {
        self.store.read().unwrap().base_completion
    }

    fn ready(&self) -> usize {
        self.store.read().unwrap().ready()
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
        if self.writers.load(Ordering::Acquire) != 0 {
            return None;
        }
        let store = match self.store.try_read() {
            Ok(store) => store,
            Err(TryLockError::WouldBlock) => return None,
            // A poisoned mathematical table must not turn an explicit finish
            // into an endless busy retry. Propagate failure to its consumer.
            Err(TryLockError::Poisoned(_)) => panic!("legacy mathematical memo is poisoned"),
        };
        if self.writers.load(Ordering::Acquire) != 0 {
            return None;
        }
        if cursor < store.base_completion {
            return None;
        }
        assert!(cursor <= store.ready(), "invalid legacy completion cursor");
        let mut rows = Vec::new();
        let mut bytes = 0usize;
        for identity in &store.completed[cursor - store.base_completion..] {
            let value = store
                .values
                .get(identity)
                .expect("completed legacy memo entry");
            let row = CompletedIntegral {
                key: IntegralKey::Legacy {
                    kind: self.kind.to_owned(),
                    nums: (self.nums)(identity),
                },
                value: value.clone(),
            };
            bytes = bytes.saturating_add(record_bytes(&row));
            rows.push(row);
            // A soft window always admits at least one complete value, even
            // when it is larger than the requested target.
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
    use sc_num::{One, Zero};
    use std::sync::Arc;

    fn st_nums<K: AsRef<[u32]>>(key: &K) -> Vec<u32> {
        key.as_ref().to_vec()
    }

    #[test]
    fn active_reclaim_never_clears_dirty_or_borrowed_legacy_memos() {
        let registry = Arc::new(Registry::new());
        let source = LegacySource::new("st", st_nums);
        source.insert(vec![1], &Q::ONE);
        super::super::with_registry(registry.clone(), || {
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            source.durable(1);
            let held = source.store.read().unwrap();
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            drop(held);
            assert!(source.reclaim_during_query(&registry) > 0);
            assert!(registry.has_history());
            assert_eq!(source.floor(), 1);
            assert!(source.get_slice(&[1]).is_none());
        });
    }

    #[test]
    fn reclaim_requires_all_commits_and_keeps_completion_positions_absolute() {
        let source = LegacySource::new("st", st_nums);
        source.insert(vec![1], &Q::from_signeds(1, 24));
        source.insert(vec![4], &Q::from_signeds(1, 1152));
        source.durable(1);
        assert_eq!(source.reclaim(), 0);
        source.durable(2);
        let held = source.store.read().unwrap();
        assert_eq!(source.reclaim(), 0);
        drop(held);
        assert!(source.reclaim() > 0);
        assert_eq!(source.floor(), 2);
        assert_eq!(source.ready(), 2);
        assert!(source.get_slice(&[1]).is_none());
        assert!(source.pull(0, 1024).is_none());
        assert!(source.pull(2, 1024).unwrap().rows.is_empty());
        source.insert(vec![7], &Q::from_signeds(1, 82944));
        let next = source.pull(2, 1024).unwrap();
        assert_eq!(next.next, 3);
        assert_eq!(next.rows[0].value, Q::from_signeds(1, 82944));
        assert_eq!(source.reclaim(), 0);
    }

    type IcKey = (u32, u32, Vec<(u32, u32)>, Vec<(u32, u32)>, Vec<(u32, u32)>);

    fn ic_nums(key: &IcKey) -> Vec<u32> {
        let (g, n, kappa, psi, ch) = key;
        [
            *g,
            *n,
            u32::try_from(kappa.len() * 2).unwrap(),
            u32::try_from(psi.len() * 2).unwrap(),
            u32::try_from(ch.len() * 2).unwrap(),
        ]
        .into_iter()
        .chain(
            kappa
                .iter()
                .chain(psi)
                .chain(ch)
                .flat_map(|&(index, power)| [index, power]),
        )
        .collect()
    }

    #[test]
    fn one_mathematical_value_and_identity_survive_duplicate_insertions() {
        let source = LegacySource::new("st", st_nums);
        let key = vec![0, 0, 3];
        assert!(source.insert(key.clone(), &Q::ZERO));
        assert!(!source.insert(key.clone(), &Q::from(7)));
        assert_eq!(source.get(&key), Some(Q::ZERO));
        assert_eq!(source.ready(), 1);
        {
            let store = source.store.read().unwrap();
            assert_eq!(store.values.len(), 1);
            assert_eq!(store.completed.as_slice(), std::slice::from_ref(&key));
        }
        let pulled = source.pull(0, 1).unwrap();
        assert_eq!(pulled.next, 1);
        assert_eq!(pulled.rows.len(), 1);
        assert_eq!(pulled.rows[0].value, Q::ZERO);
        assert_eq!(
            pulled.rows[0].key,
            IntegralKey::Legacy {
                kind: "st".into(),
                nums: key,
            }
        );
    }

    #[test]
    fn completion_order_handles_smaller_middle_and_zero_keys() {
        let source = LegacySource::new("st", st_nums);
        assert!(source.pull(0, 1024).unwrap().rows.is_empty());
        source.insert(vec![4], &Q::from(4));
        let first = source.pull(0, 0).unwrap();
        assert_eq!(first.next, 1);
        assert_eq!(first.rows[0].value, Q::from(4));
        source.insert(vec![0, 0, 3], &Q::ZERO);
        source.insert(vec![2, 2, 2], &Q::from_signeds(7, 240));
        assert_eq!(source.ready(), 3);
        let second = source.pull(first.next, 1).unwrap();
        assert_eq!(second.next, 2);
        assert_eq!(second.rows.len(), 1);
        assert_eq!(second.rows[0].value, Q::ZERO);
        let third = source.pull(second.next, usize::MAX).unwrap();
        assert_eq!(third.next, 3);
        assert_eq!(third.rows[0].value, Q::from_signeds(7, 240));
        let empty = source.pull(third.next, usize::MAX).unwrap();
        assert_eq!(empty.next, 3);
        assert!(empty.rows.is_empty());
    }

    #[test]
    fn raw_key_callback_runs_only_when_the_consumer_pulls() {
        static CALLS: AtomicUsize = AtomicUsize::new(0);
        fn nums<K: AsRef<[u32]>>(key: &K) -> Vec<u32> {
            CALLS.fetch_add(1, Ordering::Relaxed);
            key.as_ref().to_vec()
        }
        let source = LegacySource::new("st", nums);
        source.insert(vec![1], &Q::from_signeds(1, 24));
        assert_eq!(source.get(&vec![1]), Some(Q::from_signeds(1, 24)));
        assert!(!source.insert(vec![1], &Q::ONE));
        assert_eq!(source.ready(), 1);
        assert_eq!(CALLS.load(Ordering::Relaxed), 0);
        assert_eq!(source.pull(0, 1024).unwrap().rows.len(), 1);
        assert_eq!(CALLS.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn original_ic_tuple_keeps_component_boundaries_and_actual_ch() {
        let source = LegacySource::new("ic", ic_nums);
        let key = (2, 1, vec![(1, 2)], vec![(1, 1)], vec![(1, 1)]);
        let value = Q::from_signeds(1, 144);
        assert!(source.insert(key.clone(), &value));
        assert_eq!(source.get(&key), Some(value.clone()));
        let pulled = source.pull(0, 1024).unwrap();
        assert_eq!(pulled.rows[0].value, value);
        assert_eq!(
            pulled.rows[0].key,
            IntegralKey::Legacy {
                kind: "ic".into(),
                nums: vec![2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1],
            }
        );
    }

    #[test]
    fn reader_ownership_survives_the_producer_thread() {
        let source = Arc::new(LegacySource::new("st", st_nums));
        let producer = Arc::clone(&source);
        std::thread::spawn(move || {
            producer.insert(vec![2, 2, 2], &Q::from_signeds(7, 240));
            producer.insert(vec![0, 0, 3], &Q::ZERO);
        })
        .join()
        .unwrap();
        let reader: Arc<dyn ResultSource> = source;
        assert_eq!(reader.ready(), 2);
        let rows = reader.pull(0, usize::MAX).unwrap();
        assert_eq!(rows.next, 2);
        assert_eq!(rows.rows[0].value, Q::from_signeds(7, 240));
        assert_eq!(rows.rows[1].value, Q::ZERO);
    }

    #[test]
    fn write_intent_is_released_on_unwind() {
        let source = LegacySource::new("st", st_nums);
        let result = std::panic::catch_unwind(|| {
            source.writers.fetch_add(1, Ordering::AcqRel);
            let _intent = WriteIntent(&source.writers);
            assert!(source.pull(0, 1024).is_none());
            panic!("interrupted legacy write intent");
        });
        assert!(result.is_err());
        assert_eq!(source.writers.load(Ordering::Acquire), 0);
        assert!(source.insert(vec![1], &Q::ONE));
        assert_eq!(source.pull(0, 1024).unwrap().next, 1);
    }
}
