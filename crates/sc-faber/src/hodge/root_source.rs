//! 公共 κ/ψ/λ/ch 根积分的数学 L1。仅在启用结果源时建立。
//!
//! 插入顺序是完成游标，数学键的排序顺序只用于 L1 查找。值仅存一份；
//! 后台通过原 L1 自行复制，不维护另一张待写结果表。

use crate::cache_source::reuse::ReuseStats;
use crate::cache_source::{
    CompletedIntegral, IntegralKey, Pulled, Registry, ResultSource, value_bytes,
};
use sc_num::{Q, Zero};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, RwLock, Weak};

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct RootKey {
    g: u32,
    n: u32,
    kappa: Vec<(u32, u32)>,
    psi: Vec<(u32, u32)>,
    lambda: Vec<(u32, u32)>,
    ch: Vec<(u32, u32)>,
}

trait RootIdentity: Clone + Ord + Send + Sync + 'static {
    fn integral_key(&self) -> IntegralKey;
    fn payload_bytes(&self) -> usize;
    fn heap_capacity_bytes(&self) -> usize;
}

impl RootIdentity for RootKey {
    fn integral_key(&self) -> IntegralKey {
        IntegralKey::Integral {
            g: self.g,
            n: self.n,
            kappa: self.kappa.clone(),
            psi: self.psi.clone(),
            lambda: self.lambda.clone(),
            ch: self.ch.clone(),
        }
    }

    fn payload_bytes(&self) -> usize {
        (self.kappa.len() + self.psi.len() + self.lambda.len() + self.ch.len())
            .saturating_mul(std::mem::size_of::<(u32, u32)>())
    }

    fn heap_capacity_bytes(&self) -> usize {
        [&self.kappa, &self.psi, &self.lambda, &self.ch]
            .into_iter()
            .fold(0usize, |bytes, entries| {
                bytes.saturating_add(
                    entries
                        .capacity()
                        .saturating_mul(std::mem::size_of::<(u32, u32)>()),
                )
            })
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct NamedKey {
    kind: &'static str,
    nums: Vec<u32>,
}

impl RootIdentity for NamedKey {
    fn integral_key(&self) -> IntegralKey {
        IntegralKey::Legacy {
            kind: self.kind.to_owned(),
            nums: self.nums.clone(),
        }
    }

    fn payload_bytes(&self) -> usize {
        self.kind
            .len()
            .saturating_add(self.nums.len().saturating_mul(std::mem::size_of::<u32>()))
    }

    fn heap_capacity_bytes(&self) -> usize {
        self.nums
            .capacity()
            .saturating_mul(std::mem::size_of::<u32>())
    }
}

struct RootValues<K: RootIdentity> {
    index: BTreeMap<K, usize>,
    entries: Vec<(K, Q)>,
    completed: Vec<usize>,
    base_completion: usize,
}

impl<K: RootIdentity> Default for RootValues<K> {
    fn default() -> Self {
        Self {
            index: BTreeMap::new(),
            entries: Vec::new(),
            completed: Vec::new(),
            base_completion: 0,
        }
    }
}

impl<K: RootIdentity> RootValues<K> {
    fn ready(&self) -> usize {
        self.base_completion
            .checked_add(self.completed.len())
            .expect("root completion cursor exhausted")
    }
}

struct RootSource<K: RootIdentity = RootKey> {
    values: RwLock<RootValues<K>>,
    writing: AtomicBool,
    durable: AtomicUsize,
    reuse: ReuseStats,
}

impl<K: RootIdentity> Default for RootSource<K> {
    fn default() -> Self {
        Self {
            values: RwLock::new(RootValues::default()),
            writing: AtomicBool::new(false),
            durable: AtomicUsize::new(0),
            reuse: ReuseStats::default(),
        }
    }
}

// 源的强引用由 Registry 保留；TLS 不延长已关闭会话的生命周期。
thread_local! {
    static SOURCES: RefCell<BTreeMap<u64, Weak<RootSource>>> =
        const { RefCell::new(BTreeMap::new()) };
    static NAMED_SOURCES: RefCell<BTreeMap<(u64, &'static str), Weak<RootSource<NamedKey>>>> =
        const { RefCell::new(BTreeMap::new()) };
}

fn named_source(registry: &Registry, kind: &'static str) -> Arc<RootSource<NamedKey>> {
    NAMED_SOURCES.with(|cell| {
        let mut sources = cell.borrow_mut();
        if let Some(source) = sources.get(&(registry.id(), kind)).and_then(Weak::upgrade) {
            return source;
        }
        sources.retain(|_, source| source.strong_count() > 0);
        let source = Arc::new(RootSource::default());
        registry.register(source.clone());
        sources.insert((registry.id(), kind), Arc::downgrade(&source));
        source
    })
}

fn source(registry: &Registry) -> Arc<RootSource> {
    SOURCES.with(|cell| {
        let mut sources = cell.borrow_mut();
        if let Some(source) = sources.get(&registry.id()).and_then(Weak::upgrade) {
            return source;
        }
        sources.retain(|_, source| source.strong_count() > 0);
        let source = Arc::new(RootSource::default());
        registry.register(source.clone());
        sources.insert(registry.id(), Arc::downgrade(&source));
        source
    })
}

struct WriteIntent<'a>(&'a AtomicBool);

impl Drop for WriteIntent<'_> {
    fn drop(&mut self) {
        self.0.store(false, Ordering::Release);
    }
}

impl<K: RootIdentity> RootSource<K> {
    fn get(&self, key: &K) -> Option<Q> {
        let value = {
            let values = self.values.read().unwrap();
            values
                .index
                .get(key)
                .map(|&slot| values.entries[slot].1.clone())
        };
        self.reuse.record(usize::from(value.is_some()), 1);
        value
    }

    fn insert(&self, key: K, value: Q, publish: bool) -> Q {
        self.writing.store(true, Ordering::Release);
        let _intent = WriteIntent(&self.writing);
        let mut values = self.values.write().unwrap();
        if let Some(&slot) = values.index.get(&key) {
            return values.entries[slot].1.clone();
        }
        let slot = values.entries.len();
        // 索引和反向读取目录仅复制数学身份；结果移动到唯一的数学 L1 槽位。
        values.index.insert(key.clone(), slot);
        values.entries.push((key, value));
        if publish {
            values.completed.push(slot);
        }
        values.entries[slot].1.clone()
    }
}

impl<K: RootIdentity> RootSource<K> {
    fn reclaim_with_history(&self, registry: Option<&Registry>) -> usize {
        let retired = {
            let Ok(mut values) = self.values.try_write() else {
                return 0;
            };
            let ready = values.ready();
            if self.durable.load(Ordering::Acquire) != ready {
                return 0;
            }
            if let Some(registry) = registry.filter(|_| !values.entries.is_empty()) {
                registry.set_history(true);
            }
            let retired = (
                std::mem::take(&mut values.index),
                std::mem::take(&mut values.entries),
                std::mem::take(&mut values.completed),
            );
            values.base_completion = ready;
            self.reuse.reset();
            retired
        };
        let mut bytes = retired
            .0
            .len()
            .saturating_mul(std::mem::size_of::<(K, usize)>())
            .saturating_add(
                retired
                    .1
                    .capacity()
                    .saturating_mul(std::mem::size_of::<(K, Q)>()),
            )
            .saturating_add(
                retired
                    .2
                    .capacity()
                    .saturating_mul(std::mem::size_of::<usize>()),
            );
        for key in retired.0.keys() {
            bytes = bytes.saturating_add(key.heap_capacity_bytes());
        }
        for (key, value) in &retired.1 {
            bytes = bytes
                .saturating_add(key.heap_capacity_bytes())
                .saturating_add(value_bytes(value).saturating_sub(std::mem::size_of::<Q>()));
        }
        drop(retired);
        bytes
    }
}

impl<K: RootIdentity> ResultSource for RootSource<K> {
    fn floor(&self) -> usize {
        self.values.read().unwrap().base_completion
    }

    fn ready(&self) -> usize {
        self.values.read().unwrap().ready()
    }

    fn durable(&self, cursor: usize) {
        self.durable.fetch_max(cursor, Ordering::Release);
    }

    fn heat(&self) -> u64 {
        self.reuse.score()
    }

    fn reclaim(&self) -> usize {
        self.reclaim_with_history(None)
    }

    fn reclaim_during_query(&self, registry: &Registry) -> usize {
        self.reclaim_with_history(Some(registry))
    }

    fn pull(&self, cursor: usize, byte_limit: usize) -> Option<Pulled> {
        if self.writing.load(Ordering::Acquire) {
            return None;
        }
        let values = self.values.try_read().ok()?;
        if cursor < values.base_completion {
            return None;
        }
        if cursor >= values.ready() {
            return Some(Pulled {
                next: cursor,
                rows: Vec::new(),
            });
        }
        let mut rows = Vec::new();
        let mut bytes = 0usize;
        for &slot in &values.completed[cursor - values.base_completion..] {
            let (key, value) = &values.entries[slot];
            bytes = bytes
                .saturating_add(std::mem::size_of::<IntegralKey>())
                .saturating_add(key.payload_bytes())
                .saturating_add(value_bytes(value));
            rows.push(CompletedIntegral {
                key: key.integral_key(),
                value: value.clone(),
            });
            // 软窗口至少容纳一条，不能让单个较大的合法结果永远无法被读走。
            if bytes >= byte_limit {
                break;
            }
        }
        Some(Pulled {
            next: cursor + rows.len(),
            rows,
        })
    }
}

pub(crate) fn with_registry(
    registry: Arc<Registry>,
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
    calculate: impl FnOnce() -> Result<Q, String>,
) -> Result<Q, String> {
    with_registry_policy(registry, g, n, kappa, psi, lambda, ch, true, calculate)
}

/// An already validated, non-top fundamental monomial has the exact value
/// zero. Keep it in the same primary mathematical root table as Socle/KLP;
/// coefficient arithmetic and polynomial objects remain outside its identity.
pub(crate) fn complete_fund_zero(
    g: u32,
    n: u32,
    kappalist: &[u32],
    psi: &[(u32, u32)],
) -> Result<(), String> {
    crate::cache_source::with_current_registry(|registry| {
        let Some(registry) = registry else {
            return Ok(());
        };
        let mut indices = kappalist.to_vec();
        indices.sort_unstable();
        let mut kappa: Vec<(u32, u32)> = Vec::new();
        for index in indices {
            if let Some(last) = kappa.last_mut().filter(|last| last.0 == index) {
                last.1 += 1;
            } else {
                kappa.push((index, 1));
            }
        }
        let psi = crate::kppoly::merge_psi(psi);
        with_registry_policy(
            registry.clone(),
            g,
            n,
            &kappa,
            &psi,
            &[],
            &[],
            false,
            || Ok(Q::ZERO),
        )
        .map(|_| ())
    })
}

pub(crate) fn with_registry_policy(
    registry: Arc<Registry>,
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
    read_history: bool,
    calculate: impl FnOnce() -> Result<Q, String>,
) -> Result<Q, String> {
    registry.check_error()?;
    crate::cache_source::try_with_registry(registry.clone(), || {
        let source = source(&registry);
        let key = RootKey {
            g,
            n,
            kappa: kappa.to_vec(),
            psi: psi.to_vec(),
            lambda: lambda.to_vec(),
            ch: ch.to_vec(),
        };
        evaluate_root(&registry, &source, key, read_history, calculate)
    })?
}

fn evaluate_root<K: RootIdentity>(
    registry: &Registry,
    source: &RootSource<K>,
    key: K,
    read_history: bool,
    calculate: impl FnOnce() -> Result<Q, String>,
) -> Result<Q, String> {
    if let Some(value) = source.get(&key) {
        return Ok(value);
    }
    if read_history && registry.has_history() {
        if let Some(value) = registry.lookup(key.integral_key())? {
            // Loaded values are already durable; retain only the mathematical
            // L1 value, without publishing another completion for storage.
            return Ok(source.insert(key, value, false));
        }
    }
    let value = crate::cache_source::catch_resource(calculate)?;
    registry.check_error()?;
    Ok(source.insert(key, value, true))
}

pub(crate) fn with_named(
    kind: &'static str,
    key_builder: impl FnOnce() -> Vec<u32>,
    calculate: impl FnOnce() -> Result<Q, String>,
) -> Result<Q, String> {
    with_named_policy(kind, key_builder, true, calculate)
}

pub(crate) fn with_named_policy(
    kind: &'static str,
    key_builder: impl FnOnce() -> Vec<u32>,
    read_history: bool,
    calculate: impl FnOnce() -> Result<Q, String>,
) -> Result<Q, String> {
    let Some(registry) = crate::cache_source::current_registry() else {
        return crate::cache_source::with_optional_registry(None, || {
            crate::cache_source::catch_resource(calculate)
        });
    };
    registry.check_error()?;
    crate::cache_source::try_with_registry(registry.clone(), || {
        let source = named_source(&registry, kind);
        let key = NamedKey {
            kind,
            nums: key_builder(),
        };
        evaluate_root(&registry, &source, key, read_history, calculate)
    })?
}

pub(crate) fn named_i64(output: &mut Vec<u32>, value: i64) {
    let bits = value as u64;
    output.push(bits as u32);
    output.push((bits >> 32) as u32);
}

pub(crate) fn named_len(output: &mut Vec<u32>, length: usize) {
    let length = u64::try_from(length).expect("named integral slice exceeds u64");
    output.push(length as u32);
    output.push((length >> 32) as u32);
}

pub(crate) fn named_i64s(output: &mut Vec<u32>, values: &[i64]) {
    named_len(output, values.len());
    for &value in values {
        named_i64(output, value);
    }
}

pub(crate) fn named_powers(output: &mut Vec<u32>, values: &[(u32, u32)]) {
    named_len(output, values.len());
    for &(index, power) in values {
        output.extend_from_slice(&[index, power]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::{One, Zero};

    #[test]
    fn active_root_reclaim_requires_durability_and_an_unborrowed_table() {
        let registry = Arc::new(Registry::new());
        let source = RootSource::<NamedKey>::default();
        let key = NamedKey {
            kind: "test-active-reclaim",
            nums: vec![1],
        };
        source.insert(key.clone(), Q::ONE, true);
        crate::cache_source::with_registry(registry.clone(), || {
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            source.durable(1);
            let held = source.values.read().unwrap();
            assert_eq!(source.reclaim_during_query(&registry), 0);
            assert!(!registry.has_history());
            drop(held);
            assert!(source.reclaim_during_query(&registry) > 0);
            assert!(registry.has_history());
            assert_eq!(source.floor(), 1);
            assert!(source.get(&key).is_none());
        });
    }

    #[test]
    fn public_root_reuse_counts_hits_and_misses_until_whole_source_reclaim() {
        let source = RootSource::<NamedKey>::default();
        let stored = NamedKey {
            kind: "test-reuse-v1",
            nums: vec![1],
        };
        let missing = NamedKey {
            kind: "test-reuse-v1",
            nums: vec![2],
        };
        source.insert(stored.clone(), Q::ZERO, true);
        assert_eq!(source.get(&stored), Some(Q::ZERO));
        assert_eq!(source.get(&missing), None);
        assert_eq!(source.heat(), 500_000);
        assert_eq!(source.reclaim(), 0);
        assert_eq!(source.heat(), 500_000);
        source.durable(1);
        assert!(source.reclaim() > 0);
        assert_eq!(source.heat(), 0);
        assert_eq!(source.floor(), 1);
        assert_eq!(source.ready(), 1);
    }

    #[test]
    fn clean_history_and_fully_committed_roots_release_all_primary_storage() {
        let source = RootSource::default();
        let key = |g| RootKey {
            g,
            n: 0,
            kappa: Vec::new(),
            psi: Vec::new(),
            lambda: vec![(1, 3 * g - 3)],
            ch: Vec::new(),
        };
        source.insert(key(2), Q::from_signeds(1, 2880), false);
        assert_eq!(source.ready(), 0);
        assert!(source.reclaim() > 0);
        assert!(source.get(&key(2)).is_none());
        source.insert(key(3), Q::from_signeds(1, 90720), true);
        assert_eq!(source.reclaim(), 0);
        source.durable(1);
        let held = source.values.read().unwrap();
        assert_eq!(source.reclaim(), 0);
        drop(held);
        assert!(source.reclaim() > 0);
        assert_eq!(source.floor(), 1);
        assert_eq!(source.ready(), 1);
        assert!(source.pull(0, 1024).is_none());
        source.insert(key(2), Q::from_signeds(1, 2880), true);
        let next = source.pull(1, 1024).unwrap();
        assert_eq!(next.next, 2);
        assert_eq!(next.rows[0].value, Q::from_signeds(1, 2880));
        assert_eq!(source.reclaim(), 0);
    }

    fn root_rows(registry: &Registry) -> Vec<CompletedIntegral> {
        registry
            .sources()
            .into_iter()
            .flat_map(|source| {
                let mut rows = Vec::new();
                let mut cursor = 0;
                while let Some(pulled) = source.pull(cursor, usize::MAX) {
                    if pulled.rows.is_empty() {
                        assert_eq!(pulled.next, cursor);
                        break;
                    }
                    assert!(pulled.next > cursor);
                    cursor = pulled.next;
                    rows.extend(
                        pulled
                            .rows
                            .into_iter()
                            .filter(|row| matches!(&row.key, IntegralKey::Integral { .. })),
                    );
                }
                rows
            })
            .collect()
    }

    #[test]
    fn formulas_mixed_integrals_and_legal_zeros_are_completed_roots() {
        let registry = Arc::new(Registry::new());
        crate::cache_source::with_registry(registry.clone(), || {
            assert_eq!(
                super::super::klp(2, 1, &[], &[(1, 2)], &[(2, 1)], &[]).unwrap(),
                Q::from_signeds(7, 5760)
            );
            assert_eq!(
                super::super::klp(3, 2, &[(1, 1)], &[(1, 2)], &[(1, 2)], &[(1, 3)]).unwrap(),
                Q::from_signeds(1, 756)
            );
            assert_eq!(
                super::super::klp(3, 0, &[], &[], &[(3, 2)], &[]).unwrap(),
                Q::ZERO
            );
            assert_eq!(
                super::super::klp(2, 0, &[], &[], &[(3, 1)], &[]).unwrap(),
                Q::ZERO
            );
            assert_eq!(
                super::super::klp(2, 0, &[], &[], &[(1, 1)], &[]).unwrap(),
                Q::ZERO
            );
        });
        let rows = root_rows(&registry);
        assert_eq!(rows.len(), 5);
        assert_eq!(rows.iter().filter(|row| row.value == Q::ZERO).count(), 3);
    }

    #[test]
    fn normalized_root_hit_does_not_publish_another_completion() {
        let registry = Arc::new(Registry::new());
        crate::cache_source::with_registry(registry.clone(), || {
            let expected = super::super::klp(2, 2, &[], &[(1, 2)], &[(1, 3)], &[]).unwrap();
            assert_eq!(expected, Q::from_signeds(1, 1440));
            assert_eq!(
                super::super::klp(
                    2,
                    2,
                    &[(2, 0)],
                    &[(2, 0), (1, 1), (1, 1)],
                    &[(1, 1), (1, 0), (1, 2)],
                    &[(1, 0)],
                )
                .unwrap(),
                expected
            );
        });
        assert_eq!(root_rows(&registry).len(), 1);
    }

    #[test]
    fn root_l1_hit_skips_calculation_and_keeps_the_original_value() {
        let registry = Arc::new(Registry::new());
        assert_eq!(
            with_registry(registry.clone(), 0, 3, &[], &[], &[], &[], || Ok(Q::ONE)).unwrap(),
            Q::ONE
        );
        assert_eq!(
            with_registry(registry.clone(), 0, 3, &[], &[], &[], &[], || {
                panic!("completed root must be read from mathematical L1")
            })
            .unwrap(),
            Q::ONE
        );
        assert_eq!(root_rows(&registry).len(), 1);
    }

    #[test]
    fn disabled_cache_does_not_create_a_root_source() {
        let before = SOURCES.with(|sources| sources.borrow().len());
        crate::cache_source::with_optional_registry(None, || {
            assert_eq!(super::super::klp(0, 3, &[], &[], &[], &[]).unwrap(), Q::ONE);
        });
        assert_eq!(SOURCES.with(|sources| sources.borrow().len()), before);
    }

    #[test]
    fn invalid_public_inputs_never_become_completed_roots() {
        let registry = Arc::new(Registry::new());
        crate::cache_source::with_registry(registry.clone(), || {
            assert!(super::super::klp(0, 2, &[], &[], &[], &[]).is_err());
            assert!(super::super::klp(2, 1, &[], &[(2, 4)], &[], &[]).is_err());
            assert!(super::super::klp(2, 0, &[(0, 1)], &[], &[], &[]).is_err());
            assert!(super::super::klp(2, 0, &[], &[], &[(1, 1)], &[(2, 1)]).is_err());
        });
        assert!(root_rows(&registry).is_empty());
    }

    #[test]
    fn root_source_survives_query_stack_and_calling_thread() {
        let registry = Arc::new(Registry::new());
        let child_registry = registry.clone();
        std::thread::spawn(move || {
            crate::cache_source::with_registry(child_registry, || {
                assert_eq!(super::super::klp(0, 3, &[], &[], &[], &[]).unwrap(), Q::ONE);
            });
        })
        .join()
        .unwrap();
        let rows = root_rows(&registry);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].value, Q::ONE);
    }

    #[test]
    fn registry_changes_keep_root_sources_separate() {
        let first = Arc::new(Registry::new());
        let second = Arc::new(Registry::new());
        for registry in [&first, &second, &first] {
            crate::cache_source::with_registry(registry.clone(), || {
                assert_eq!(super::super::klp(0, 3, &[], &[], &[], &[]).unwrap(), Q::ONE);
            });
        }
        assert_eq!(root_rows(&first).len(), 1);
        assert_eq!(root_rows(&second).len(), 1);
    }

    #[test]
    fn insertion_cursor_includes_new_keys_before_the_previous_key() {
        let source = RootSource::default();
        let key = |g| RootKey {
            g,
            n: 0,
            kappa: vec![],
            psi: vec![],
            lambda: vec![],
            ch: vec![],
        };
        source.insert(key(4), Q::from(4), true);
        let first = source.pull(0, 0).unwrap();
        assert_eq!(first.next, 1);
        assert_eq!(first.rows.len(), 1);
        assert_eq!(first.rows[0].value, Q::from(4));
        source.insert(key(2), Q::ZERO, true);
        source.insert(key(3), Q::from(3), true);
        let second = source.pull(first.next, 0).unwrap();
        assert_eq!(second.next, 2);
        assert_eq!(second.rows.len(), 1);
        assert_eq!(second.rows[0].value, Q::ZERO);
        let third = source.pull(second.next, usize::MAX).unwrap();
        assert_eq!(third.next, 3);
        assert_eq!(third.rows[0].value, Q::from(3));
        let empty = source.pull(third.next, usize::MAX).unwrap();
        assert_eq!(empty.next, third.next);
        assert!(empty.rows.is_empty());
    }

    #[test]
    fn named_keys_preserve_signed_width_slice_boundaries_and_namespaces() {
        let mut signed = Vec::new();
        for value in [i64::MIN, i64::MAX, -1, 0, 1i64 << 32] {
            named_i64(&mut signed, value);
        }
        assert_eq!(
            signed,
            vec![
                0,
                0x8000_0000,
                u32::MAX,
                0x7fff_ffff,
                u32::MAX,
                u32::MAX,
                0,
                0,
                0,
                1
            ]
        );
        let mut first = Vec::new();
        named_i64s(&mut first, &[1, 2]);
        named_powers(&mut first, &[]);
        let mut second = Vec::new();
        named_i64s(&mut second, &[1]);
        named_powers(&mut second, &[(2, 0)]);
        assert_ne!(first, second);
        let mut empty = Vec::new();
        named_powers(&mut empty, &[]);
        let mut zero = Vec::new();
        named_powers(&mut zero, &[(1, 0)]);
        assert_ne!(empty, zero);

        let registry = Arc::new(Registry::new());
        crate::cache_source::with_registry(registry.clone(), || {
            for (kind, value) in [("test-named-a-v1", Q::ONE), ("test-named-b-v1", Q::from(2))] {
                assert_eq!(
                    with_named(kind, || signed.clone(), || Ok(value.clone())).unwrap(),
                    value
                );
                assert_eq!(
                    with_named(kind, || signed.clone(), || panic!("named L1 hit")).unwrap(),
                    value
                );
                assert_eq!(named_source(&registry, kind).ready(), 1);
            }
        });
        assert_eq!(registry.sources().len(), 2);
        crate::cache_source::with_optional_registry(None, || {
            assert_eq!(
                with_named(
                    "test-no-cache-v1",
                    || panic!("disabled key builder"),
                    || Ok(Q::ONE)
                )
                .unwrap(),
                Q::ONE
            );
        });
    }

    fn named_scalar_rows(registry: &Registry) -> Vec<CompletedIntegral> {
        registry.sources().into_iter().flat_map(|source| {
            let mut cursor = source.floor();
            let mut rows = Vec::new();
            while let Some(pulled) = source.pull(cursor, usize::MAX) {
                if pulled.rows.is_empty() { break; }
                cursor = pulled.next;
                rows.extend(pulled.rows.into_iter().filter(|row| {
                    matches!(&row.key, IntegralKey::Legacy { kind, .. } if kind.starts_with("scalar-"))
                }));
            }
            rows
        }).collect()
    }

    #[test]
    fn named_scalar_validation_and_original_early_zero_priority_are_preserved() {
        let registry = Arc::new(Registry::new());
        crate::cache_source::with_registry(registry.clone(), || {
            assert!(crate::dr::int_dr(1, &[1, 0], None, &[], &[]).is_err());
            assert!(crate::dr::int_dr(1, &[1, -1], None, &[], &[(3, 1)]).is_err());
            assert!(crate::witten::int_witten(1, &[(1, 1, -1)], &[], &[]).is_err());
            assert!(crate::witten::int_witten(1, &[(0, 0, 0)], &[], &[]).is_err());
            assert!(crate::classes::int_dr_chiodo_poly(0, &[], 0, 0, false, &[], &[]).is_err());
            assert!(crate::classes::int_dr_chiodo_poly(1, &[2, 0], 1, 1, true, &[], &[]).is_err());
            assert!(crate::classes::hyperell_int(2, 7, 0, &[], &[]).is_err());
            assert!(crate::classes::hyperell_int(2, 1, 0, &[], &[(2, 1)]).is_err());
            assert!(crate::classes::hyperell_pair(2, 1, 0, &[], &[(0, 0)]).is_err());
            assert!(crate::classes::biell_int(1, 0, 0, &[], &[]).is_err());
            assert!(crate::classes::biell_pair(2, 3, 0, &[], &[]).is_err());
        });
        assert!(registry.sources().is_empty());
        crate::cache_source::with_registry(registry.clone(), || {
            assert_eq!(
                crate::dr::int_dr(1, &[1, -1], Some(2), &[], &[]).unwrap(),
                Q::ZERO
            );
            // Historically this Witten zero precedes denominator validation.
            let early_zero = [(0, 0, -1), (0, -1, 1)];
            assert_eq!(
                crate::witten::int_witten(1, &early_zero, &[], &[]).unwrap(),
                Q::ZERO
            );
            assert!(crate::witten::int_witten(1, &early_zero, &[], &[(3, 1)]).is_err());
            assert_eq!(
                crate::classes::hyperell_pair(1, 1, 1, &[], &[(5, 1)]).unwrap(),
                Q::ZERO
            );
        });
        let rows = named_scalar_rows(&registry);
        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|row| row.value == Q::ZERO));
    }

    #[test]
    fn named_scalar_roots_reopen_with_exact_values_without_republication() {
        use crate::cache_source::storage::Storage;
        use std::path::PathBuf;
        static NEXT_DATABASE: AtomicUsize = AtomicUsize::new(0);
        struct Database(PathBuf);
        impl Drop for Database {
            fn drop(&mut self) {
                for suffix in ["", "-wal", "-shm"] {
                    let mut path = self.0.as_os_str().to_os_string();
                    path.push(suffix);
                    let _ = std::fs::remove_file(PathBuf::from(path));
                }
            }
        }
        let database = Database(std::env::temp_dir().join(format!(
            "sc-faber-named-roots-{}-{}.sqlite",
            std::process::id(),
            NEXT_DATABASE.fetch_add(1, Ordering::Relaxed)
        )));
        let calculate = || {
            assert_eq!(
                crate::dr::int_dr(2, &[], None, &[(1, 1)], &[]).unwrap(),
                Q::from_signeds(7, 5760)
            );
            assert_eq!(
                crate::witten::int_witten(1, &[(0, 1, 2), (1, 1, -2)], &[], &[(1, 1)]).unwrap(),
                Q::from_signeds(-1, 8)
            );
            assert_eq!(
                crate::classes::int_dr_chiodo_poly(
                    1,
                    &[2, 3, -5],
                    1,
                    0,
                    false,
                    &[],
                    &[(1, 1), (2, 1)]
                )
                .unwrap(),
                Q::from_signeds(3, 2)
            );
            assert_eq!(
                crate::classes::int_dr_chiodo_poly(1, &[1, 1], 1, 1, true, &[(1, 1)], &[]).unwrap(),
                Q::from_signeds(1, 24)
            );
            assert_eq!(
                crate::classes::hyperell_int(1, 1, 1, &[], &[(1, 2)]).unwrap(),
                Q::from_signeds(1, 24)
            );
            assert_eq!(
                crate::classes::biell_int(2, 0, 0, &[(2, 1)], &[]).unwrap(),
                Q::from_signeds(1, 32)
            );
            assert_eq!(
                crate::classes::hyperell_pair(2, 6, 0, &[], &[(1, 3)]).unwrap(),
                Q::from_signeds(1, 16)
            );
            assert_eq!(
                crate::classes::biell_pair(2, 1, 0, &[], &[(1, 2)]).unwrap(),
                Q::from_signeds(1, 48)
            );
            assert_eq!(
                crate::dr::int_dr(1, &[1, -1], Some(2), &[], &[]).unwrap(),
                Q::ZERO
            );
        };
        let first = Arc::new(Registry::new());
        let storage = Storage::open(database.0.to_str().unwrap(), first.clone()).unwrap();
        first.bind_storage(storage.handle());
        crate::cache_source::with_registry(first.clone(), calculate);
        assert_eq!(named_scalar_rows(&first).len(), 9);
        storage.close().unwrap();

        let second = Arc::new(Registry::new());
        let reopened = Storage::open(database.0.to_str().unwrap(), second.clone()).unwrap();
        second.bind_storage(reopened.handle());
        assert!(second.has_history());
        // Every nontrivial named result loads into its sole mathematical L1.
        // The cheap DR zero is intentionally recomputed and published.
        crate::cache_source::with_registry(second.clone(), calculate);
        let rows = named_scalar_rows(&second);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].value, Q::ZERO);
        for kind in [
            "scalar-witten-v1",
            "scalar-chiodo-v1",
            "scalar-hyperell-v1",
            "scalar-biell-v1",
            "scalar-hyperell-pair-v1",
            "scalar-biell-pair-v1",
        ] {
            assert_eq!(named_source(&second, kind).ready(), 0);
        }
        reopened.close().unwrap();
    }
}
