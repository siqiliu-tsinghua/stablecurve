//! Shared access to the original serial mathematical memo.
//!
//! Values live only in `SerialStore::values`. The append-only directory contains
//! identities, never copies of Q. Its insertion order lets a reader find new
//! results even when a mathematical key sorts before keys it has already read.
//! Identity copies and directory capacity are additional cache costs; unlike a
//! BTree key cursor, they provide a complete incremental completion history.

use crate::cache_source::{CompletedIntegral, IntegralKey, Pulled, Registry, ResultSource};
use sc_num::Q;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{RwLock, RwLockReadGuard};

pub(super) type SerialValues = BTreeMap<Vec<(u32, u32)>, BTreeMap<Vec<u32>, Q>>;

struct CompletedIdentity {
    g: u32,
    d: Vec<u32>,
    ch: Vec<(u32, u32)>,
}

#[derive(Default)]
pub(super) struct SerialStore {
    values: SerialValues,
    completed: Vec<CompletedIdentity>,
    base_completion: usize,
}

impl SerialStore {
    pub(super) fn values(&self) -> &SerialValues {
        &self.values
    }

    fn ready(&self) -> usize {
        self.base_completion
            .checked_add(self.completed.len())
            .expect("serial completion cursor exhausted")
    }
}

#[derive(Default)]
pub(super) struct SerialSource {
    store: RwLock<SerialStore>,
    writers: AtomicUsize,
    durable: AtomicUsize,
}

struct WriteIntent<'a>(&'a AtomicUsize);

impl Drop for WriteIntent<'_> {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::Release);
    }
}

impl SerialSource {
    pub(super) fn read(&self) -> RwLockReadGuard<'_, SerialStore> {
        self.store.read().unwrap()
    }

    pub(super) fn get(&self, d: &[u32], ch: &[(u32, u32)]) -> Option<Q> {
        self.read()
            .values
            .get(ch)
            .and_then(|row| row.get(d))
            .cloned()
    }

    fn write_intent(&self) -> WriteIntent<'_> {
        self.writers.fetch_add(1, Ordering::AcqRel);
        WriteIntent(&self.writers)
    }

    /// Insert a normally memoized value and publish its complete identity.
    ///
    /// The Q clone is the original serial memo insertion, not a persistence
    /// copy. The extra d/ch copies form completion metadata and must be budgeted
    /// and measured. A repeated key does not create another completion record.
    pub(super) fn insert(&self, g: u32, d: Vec<u32>, ch: &[(u32, u32)], value: &Q) -> bool {
        let _intent = self.write_intent();
        let mut store = self.store.write().unwrap();
        if store.values.get(ch).is_some_and(|row| row.contains_key(&d)) {
            return false;
        }
        let identity = CompletedIdentity {
            g,
            d: d.clone(),
            ch: ch.to_vec(),
        };
        store
            .values
            .entry(ch.to_vec())
            .or_default()
            .insert(d, value.clone());
        // The table and completion directory become visible under one lock.
        // No recursive mathematics, budget wait or I/O runs under this guard.
        store.completed.push(identity);
        true
    }
}

impl SerialSource {
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
        // These are known payload/capacity bytes, excluding BTree node metadata
        // and allocator retention. Large scans and destructors run unlocked.
        let mut bytes = retired
            .1
            .capacity()
            .saturating_mul(std::mem::size_of::<CompletedIdentity>());
        for identity in &retired.1 {
            bytes = bytes
                .saturating_add(
                    identity
                        .d
                        .capacity()
                        .saturating_mul(std::mem::size_of::<u32>()),
                )
                .saturating_add(
                    identity
                        .ch
                        .capacity()
                        .saturating_mul(std::mem::size_of::<(u32, u32)>()),
                );
        }
        bytes = bytes.saturating_add(retired.0.len().saturating_mul(std::mem::size_of::<(
            Vec<(u32, u32)>,
            std::collections::BTreeMap<Vec<u32>, Q>,
        )>()));
        for (ch, row) in &retired.0 {
            bytes = bytes
                .saturating_add(
                    ch.capacity()
                        .saturating_mul(std::mem::size_of::<(u32, u32)>()),
                )
                .saturating_add(
                    row.len()
                        .saturating_mul(std::mem::size_of::<(Vec<u32>, Q)>()),
                );
            for (d, value) in row {
                bytes = bytes
                    .saturating_add(d.capacity().saturating_mul(std::mem::size_of::<u32>()))
                    .saturating_add(heap_min(value));
            }
        }
        drop(retired);
        bytes
    }
}

impl ResultSource for SerialSource {
    fn floor(&self) -> usize {
        self.read().base_completion
    }

    fn ready(&self) -> usize {
        self.read().ready()
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
        let store = self.store.try_read().ok()?;
        if self.writers.load(Ordering::Acquire) != 0 {
            return None;
        }
        if cursor < store.base_completion {
            return None;
        }
        assert!(cursor <= store.ready(), "invalid serial completion cursor");
        let mut next = cursor;
        let mut rows = Vec::new();
        let mut bytes = 0usize;
        for identity in &store.completed[cursor - store.base_completion..] {
            let value = &store.values[identity.ch.as_slice()][identity.d.as_slice()];
            let row_bytes = std::mem::size_of::<CompletedIntegral>()
                .saturating_add(identity.d.len().saturating_mul(std::mem::size_of::<u32>()))
                .saturating_add(
                    identity
                        .ch
                        .len()
                        .saturating_mul(std::mem::size_of::<(u32, u32)>()),
                )
                .saturating_add(heap_min(value));
            // One large result may exceed the soft window so the source can
            // make progress. A hard-budget caller must reserve that headroom.
            if !rows.is_empty() && bytes.saturating_add(row_bytes) > byte_limit {
                break;
            }
            rows.push(CompletedIntegral {
                key: IntegralKey::Descendant {
                    g: identity.g,
                    d: identity.d.clone(),
                    ch: identity.ch.clone(),
                },
                value: value.clone(),
            });
            next += 1;
            bytes = bytes.saturating_add(row_bytes);
            if bytes >= byte_limit {
                break;
            }
        }
        Some(Pulled { next, rows })
    }
}

fn heap_min(value: &Q) -> usize {
    let limb_bytes = value
        .denominator_ref()
        .limbs()
        .next()
        .map_or(8, |limb| std::mem::size_of_val(&limb));
    [value.numerator_ref(), value.denominator_ref()]
        .into_iter()
        .map(|integer| match integer.limb_count() as usize {
            0 | 1 => 0,
            count => count.saturating_mul(limb_bytes),
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::Zero;

    #[test]
    fn only_fully_durable_unlocked_memos_reclaim_with_monotone_cursors() {
        let source = SerialSource::default();
        source.insert(1, vec![1], &[], &Q::from_signeds(1, 24));
        source.insert(2, vec![4], &[], &Q::from_signeds(1, 1152));
        assert_eq!(source.reclaim(), 0);
        source.durable(1);
        assert_eq!(source.reclaim(), 0);
        source.durable(2);
        let held = source.read();
        assert_eq!(source.reclaim(), 0);
        drop(held);
        assert!(source.reclaim() > 0);
        assert!(source.read().values.is_empty());
        assert_eq!(source.floor(), 2);
        assert_eq!(source.ready(), 2);
        assert!(source.pull(0, 1024).is_none());
        assert!(source.pull(2, 1024).unwrap().rows.is_empty());
        source.insert(3, vec![7], &[], &Q::from_signeds(1, 82944));
        assert_eq!(source.ready(), 3);
        let next = source.pull(2, 1024).unwrap();
        assert_eq!(next.next, 3);
        assert_eq!(next.rows[0].value, Q::from_signeds(1, 82944));
        assert_eq!(source.reclaim(), 0);
    }

    #[test]
    fn active_reclaim_respects_serial_memo_borrows_and_pending_completions() {
        let registry = std::sync::Arc::new(Registry::new());
        let source = SerialSource::default();
        source.insert(1, vec![1], &[], &Q::from_signeds(1, 24));
        crate::cache_source::with_registry(registry.clone(), || {
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            source.durable(1);
            // with_memo uses this same guard for its externally borrowed table.
            let held = source.read();
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            drop(held);
            assert!(source.reclaim_during_query(&registry) > 0);
            assert!(registry.has_history());
            assert_eq!(source.floor(), 1);
            assert!(source.get(&[1], &[]).is_none());
        });
    }

    #[test]
    fn insertion_order_directory_does_not_skip_smaller_or_middle_keys() {
        let source = SerialSource::default();
        assert_eq!(source.ready(), 0);
        assert!(source.insert(1, vec![1, 1, 1], &[], &Q::from(3)));
        let first = source.pull(0, 1).unwrap();
        assert_eq!(first.next, 1);
        assert_eq!(first.rows.len(), 1);
        assert!(source.insert(1, vec![0, 0, 3], &[], &Q::ZERO));
        assert!(source.insert(1, vec![0, 1, 2], &[], &Q::from(2)));
        let second = source.pull(first.next, usize::MAX).unwrap();
        assert_eq!(second.next, 3);
        assert_eq!(second.rows.len(), 2);
        assert_eq!(second.rows[0].value, Q::ZERO);
        assert_eq!(second.rows[1].value, Q::from(2));
        assert_eq!(source.get(&[0, 0, 3], &[]), Some(Q::ZERO));
        assert!(source.pull(second.next, 1024).unwrap().rows.is_empty());
        assert!(!source.insert(1, vec![0, 0, 3], &[], &Q::from(7)));
        assert_eq!(source.get(&[0, 0, 3], &[]), Some(Q::ZERO));
        assert_eq!(source.read().completed.len(), 3);
        assert_eq!(source.ready(), 3);
    }

    #[test]
    fn pull_yields_to_write_intent_and_intent_is_scoped() {
        let source = SerialSource::default();
        let intent = source.write_intent();
        assert!(source.pull(0, 1024).is_none());
        drop(intent);
        let pulled = source.pull(0, 1024).unwrap();
        assert_eq!(pulled.next, 0);
        assert!(pulled.rows.is_empty());
    }

    #[test]
    fn bounded_pull_keeps_all_markings_and_c_powers() {
        let source = SerialSource::default();
        source.insert(2, vec![0, 0, 3], &[(1, 2)], &Q::from(5));
        source.insert(2, vec![0, 3], &[(1, 2)], &Q::from(6));
        let first = source.pull(0, 1).unwrap();
        assert_eq!(first.next, 1);
        assert_eq!(first.rows.len(), 1);
        match &first.rows[0].key {
            IntegralKey::Descendant { g, d, ch } => {
                assert_eq!(*g, 2);
                assert_eq!(d, &[0, 0, 3]);
                assert_eq!(ch, &[(1, 2)]);
            }
            _ => panic!("serial memo emitted a non-descendant key"),
        }
        let second = source.pull(first.next, 1).unwrap();
        assert_eq!(second.next, 2);
        assert_eq!(second.rows.len(), 1);
    }

    #[test]
    fn formula_and_legal_zero_results_have_a_real_memo_source() {
        let registry = std::sync::Arc::new(crate::cache_source::Registry::default());
        let mut cached = super::super::Descendents::with_registry(Some(&registry));
        let mut original = super::super::Descendents::default();
        let cases = [
            (0, vec![0, 0, 0], vec![]),
            (2, vec![4], vec![]),
            (1, vec![0], vec![(1, 1)]),
            (0, vec![0, 0, 0, 0], vec![(1, 1)]),
            (2, vec![0, 0], vec![(5, 1)]),
        ];
        for (g, d, ch) in &cases {
            let expected = original.eval(*g, d.clone(), ch);
            assert_eq!(cached.eval(*g, d.clone(), ch), expected);
            assert_eq!(cached.eval(*g, d.clone(), ch), expected);
            assert_eq!(cached.source.as_ref().unwrap().get(d, ch), Some(expected));
        }
        // Dimension errors and unstable conventions are not completed stable
        // integrals and must not enter this source.
        assert_eq!(cached.eval(0, vec![], &[]), Q::ZERO);
        assert_eq!(cached.eval(2, vec![0], &[]), Q::ZERO);
        assert!(cached.values.is_empty());
        drop(cached);
        let sources = registry.sources();
        assert_eq!(sources.len(), 1);
        let rows = sources[0].pull(0, usize::MAX).unwrap();
        assert_eq!(rows.rows.len(), cases.len());
    }

    #[test]
    fn recursive_mathematics_uses_only_the_shared_memo() {
        let registry = std::sync::Arc::new(crate::cache_source::Registry::default());
        let mut cached = super::super::Descendents::with_registry(Some(&registry));
        let mut original = super::super::Descendents::default();
        for (g, d, ch) in [
            (3, vec![2, 3, 4], vec![]),
            (2, vec![], vec![(1, 3)]),
            (2, vec![0, 3], vec![(1, 2)]),
        ] {
            let expected = original.eval(g, d.clone(), &ch);
            assert_eq!(cached.eval(g, d.clone(), &ch), expected);
            let source = cached.source.as_ref().unwrap().clone();
            let cursor = source.pull(0, usize::MAX).unwrap().next;
            assert_eq!(cached.eval(g, d, &ch), expected);
            assert!(source.pull(cursor, usize::MAX).unwrap().rows.is_empty());
        }
        assert!(cached.values.is_empty());
        assert!(!original.values.is_empty());
        assert!(cached.with_memo(|values| values.values().any(|row| !row.is_empty())));
    }

    #[test]
    fn registry_contexts_reuse_their_own_memo_across_switches() {
        use super::super::with_descendents;
        use std::sync::Arc;

        let first = Arc::new(crate::cache_source::Registry::default());
        let second = Arc::new(crate::cache_source::Registry::default());
        let calculate = |evaluator: &mut super::super::Descendents| evaluator.eval(2, vec![4], &[]);
        let expected = with_descendents(Some(&first), calculate);
        assert_eq!(with_descendents(Some(&second), calculate), expected);
        assert_eq!(with_descendents(Some(&first), calculate), expected);
        assert_eq!(first.sources().len(), 1);
        assert_eq!(second.sources().len(), 1);
        assert_eq!(
            first.sources()[0].pull(0, usize::MAX).unwrap().rows.len(),
            1
        );
        assert_eq!(
            second.sources()[0].pull(0, usize::MAX).unwrap().rows.len(),
            1
        );
    }
}
