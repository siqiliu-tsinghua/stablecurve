//! Long-lived access to the original contiguous result table. Only the reader
//! constructs persistence rows; the planner publishes existing state batches.

use super::{Index, State};
use crate::cache_source::reuse::ReuseStats;
use crate::cache_source::{
    CompletedIntegral, IntegralKey, Pulled, Registry, ResultSource, encode_descendant_prefix_into,
    value_bytes,
};
use sc_num::Q;
use std::collections::VecDeque;
use std::ops::{Deref, DerefMut};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock, RwLockReadGuard, RwLockWriteGuard, TryLockError};

#[derive(Default)]
pub(super) struct ParallelSource {
    // A query's numeric IDs can outlive any one value/index guard. Hold this
    // read permit for their complete mathematical use, including worker joins.
    usage: RwLock<()>,
    values: RwLock<Vec<Q>>,
    writing: AtomicBool,
    states: Mutex<VecDeque<Arc<Vec<State>>>>,
    directory: Mutex<Directory>,
    batches: Mutex<Vec<Arc<Batch>>>,
    ready: AtomicUsize,
    durable: AtomicUsize,
    floor: AtomicUsize,
    reuse: ReuseStats,
    // These are the original mathematical index locks, not another key table.
    shards: Mutex<Vec<Arc<RwLock<Index>>>>,
}

pub(super) struct Batch {
    base: usize,
    start: usize,
    pub(super) layers: Vec<Vec<usize>>,
    dirty_layers: Option<Vec<Vec<usize>>>,
    offsets: Vec<usize>,
    published: AtomicUsize,
}

pub(super) struct WriteAccess<'a> {
    guard: Option<RwLockWriteGuard<'a, Vec<Q>>>,
    writing: &'a AtomicBool,
}

impl Deref for WriteAccess<'_> {
    type Target = Vec<Q>;

    fn deref(&self) -> &Self::Target {
        self.guard.as_ref().unwrap()
    }
}

impl DerefMut for WriteAccess<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.guard.as_mut().unwrap()
    }
}

impl Drop for WriteAccess<'_> {
    fn drop(&mut self) {
        drop(self.guard.take());
        self.writing.store(false, Ordering::Release);
    }
}

impl ParallelSource {
    pub(super) fn begin_use(&self) -> RwLockReadGuard<'_, ()> {
        self.usage.read().unwrap_or_else(|error| error.into_inner())
    }

    fn try_idle(&self) -> Option<RwLockWriteGuard<'_, ()>> {
        match self.usage.try_write() {
            Ok(guard) => Some(guard),
            Err(TryLockError::Poisoned(error)) => {
                self.usage.clear_poison();
                Some(error.into_inner())
            }
            Err(TryLockError::WouldBlock) => None,
        }
    }

    pub(super) fn observe_shards(&self, shards: Vec<Arc<RwLock<Index>>>) {
        if !shards.is_empty() {
            self.shards.lock().unwrap().extend(shards);
        }
    }

    /// Observe resolved query roots that reuse a preexisting parallel L1 ID.
    /// Formula constants and recursive descendant accesses are outside this
    /// coarse observation; no per-node or arithmetic-term counter is updated.
    pub(super) fn note_query(&self, reused_roots: usize, total_roots: usize) {
        self.reuse.record(reused_roots, total_roots);
    }

    pub(super) fn read(&self) -> RwLockReadGuard<'_, Vec<Q>> {
        self.values
            .read()
            .unwrap_or_else(|error| error.into_inner())
    }

    pub(super) fn write(&self) -> WriteAccess<'_> {
        self.writing.store(true, Ordering::Release);
        let guard = self
            .values
            .write()
            .unwrap_or_else(|error| error.into_inner());
        // Publication, rather than placeholder contents, defines completion.
        // An interrupted write cannot invalidate earlier immutable results.
        self.values.clear_poison();
        WriteAccess {
            guard: Some(guard),
            writing: &self.writing,
        }
    }

    pub(super) fn observe(&self, states: &Arc<Vec<State>>) {
        self.states.lock().unwrap().push_back(Arc::clone(states));
    }

    pub(super) fn begin(&self, base: usize, layers: Vec<Vec<usize>>) -> Arc<Batch> {
        self.begin_with(base, layers, None)
    }

    pub(super) fn begin_restored(
        &self,
        base: usize,
        layers: Vec<Vec<usize>>,
        mut restored: Vec<usize>,
    ) -> Arc<Batch> {
        restored.sort_unstable();
        let dirty = layers
            .iter()
            .map(|layer| {
                layer
                    .iter()
                    .copied()
                    .filter(|id| restored.binary_search(&(base + *id)).is_err())
                    .collect()
            })
            .collect();
        self.begin_with(base, layers, Some(dirty))
    }

    fn begin_with(
        &self,
        base: usize,
        layers: Vec<Vec<usize>>,
        dirty_layers: Option<Vec<Vec<usize>>>,
    ) -> Arc<Batch> {
        let mut offsets = Vec::with_capacity(layers.len() + 1);
        offsets.push(0);
        for layer in dirty_layers.as_deref().unwrap_or(&layers) {
            offsets.push(offsets.last().unwrap() + layer.len());
        }
        let batch = Arc::new(Batch {
            base,
            start: self.ready.load(Ordering::Acquire),
            layers,
            dirty_layers,
            offsets,
            published: AtomicUsize::new(0),
        });
        self.batches.lock().unwrap().push(Arc::clone(&batch));
        batch
    }

    pub(super) fn publish(&self, batch: &Batch, layer: usize) {
        // The coordinator has already released the ordinary table write access.
        batch.published.store(layer + 1, Ordering::Release);
        self.ready
            .store(batch.start + batch.offsets[layer + 1], Ordering::Release);
    }

    pub(super) fn completed_since(&self, base: usize, assigned: usize) -> Vec<bool> {
        let mut complete = vec![false; assigned.saturating_sub(base)];
        for batch in self
            .batches
            .lock()
            .unwrap()
            .iter()
            .filter(|batch| batch.base >= base)
        {
            let published = batch.published.load(Ordering::Acquire);
            for layer in &batch.layers[..published] {
                for &id in layer {
                    complete[batch.base + id - base] = true;
                }
            }
        }
        complete
    }

    /// The caller holds the exclusive source-use permit. Numeric IDs cannot
    /// remain live while this clears their original index and value table.
    /// Try every access before changing anything so a concurrent consumer's
    /// short read window can simply postpone this pressure-only operation.
    fn try_reclaim(&self, count_payload: bool) -> Option<usize> {
        if self.durable.load(Ordering::Acquire) < self.ready.load(Ordering::Acquire) {
            return None;
        }
        let mut directory = self.directory.try_lock().ok()?;
        let mut states = self.states.try_lock().ok()?;
        let mut batches = self.batches.try_lock().ok()?;
        let shards = self.shards.try_lock().ok()?;
        let mut indexes = Vec::with_capacity(shards.len());
        for shard in shards.iter() {
            indexes.push(shard.try_write().ok()?);
        }
        let mut values = self.values.try_write().ok()?;
        let ready = self.ready.load(Ordering::Acquire);
        if self.durable.load(Ordering::Acquire) < ready || self.writing.load(Ordering::Acquire) {
            return None;
        }

        let bytes = if count_payload {
            // Effective payload/capacity lower bound. BTree node allocation and
            // allocator spare space are intentionally not presented as exact RSS.
            // Successful close only needs to release these same allocations.
            let mut bytes = values.capacity().saturating_mul(size_of::<Q>());
            for value in values.iter() {
                bytes = bytes.saturating_add(value_bytes(value).saturating_sub(size_of::<Q>()));
            }
            for index in &indexes {
                for (ch, powers) in index.iter() {
                    bytes =
                        bytes.saturating_add(ch.capacity().saturating_mul(size_of::<(u32, u32)>()));
                    for d in powers.keys() {
                        bytes = bytes.saturating_add(d.capacity().saturating_mul(size_of::<u32>()));
                    }
                }
            }
            bytes = bytes
                .saturating_add(directory.bytes.capacity())
                .saturating_add(directory.encoded_scratch.capacity())
                .saturating_add(
                    directory
                        .offsets
                        .capacity()
                        .saturating_mul(size_of::<usize>()),
                )
                .saturating_add(
                    states
                        .capacity()
                        .saturating_mul(size_of::<Arc<Vec<State>>>()),
                )
                .saturating_add(batches.capacity().saturating_mul(size_of::<Arc<Batch>>()));
            for state in states.iter() {
                bytes = bytes.saturating_add(state_bytes(state));
            }
            if let Some((state, _)) = &directory.active {
                bytes = bytes.saturating_add(state_bytes(state));
            }
            for batch in batches.iter() {
                if Arc::strong_count(batch) == 1 {
                    bytes = bytes
                        .saturating_add(size_of::<Batch>())
                        .saturating_add(layer_bytes(&batch.layers))
                        .saturating_add(
                            batch.offsets.capacity().saturating_mul(size_of::<usize>()),
                        );
                    if let Some(layers) = &batch.dirty_layers {
                        bytes = bytes.saturating_add(layer_bytes(layers));
                    }
                }
            }
            bytes
        } else {
            0
        };

        for index in &mut indexes {
            **index = Index::new();
        }
        *values = Vec::new();
        *directory = Directory::default();
        *states = VecDeque::new();
        *batches = Vec::new();
        self.reuse.reset();
        // Keep ready/durable monotone. The live storage cursor is already at
        // this committed prefix; a new consumer starts from floor instead.
        self.floor.store(ready, Ordering::Release);
        Some(bytes)
    }
}

fn layer_bytes(layers: &Vec<Vec<usize>>) -> usize {
    layers.iter().fold(
        layers.capacity().saturating_mul(size_of::<Vec<usize>>()),
        |bytes, layer| bytes.saturating_add(layer.capacity().saturating_mul(size_of::<usize>())),
    )
}

fn state_bytes(states: &Arc<Vec<State>>) -> usize {
    if Arc::strong_count(states) != 1 {
        return 0;
    }
    states.iter().fold(
        states.capacity().saturating_mul(size_of::<State>()),
        |bytes, state| {
            bytes
                .saturating_add(state.d.capacity().saturating_mul(size_of::<u32>()))
                .saturating_add(state.ch.capacity().saturating_mul(size_of::<(u32, u32)>()))
        },
    )
}

#[derive(Default)]
struct Directory {
    offsets: Vec<usize>,
    bytes: Vec<u8>,
    encoded_scratch: Vec<u8>,
    active: Option<(Arc<Vec<State>>, usize)>,
}

fn put(mut value: u32, bytes: &mut Vec<u8>) {
    while value >= 128 {
        bytes.push((value as u8 & 127) | 128);
        value >>= 7;
    }
    bytes.push(value as u8);
}

fn take(bytes: &mut &[u8]) -> u32 {
    let mut value = 0;
    let mut shift = 0;
    loop {
        let byte = bytes[0];
        *bytes = &bytes[1..];
        value |= u32::from(byte & 127) << shift;
        if byte < 128 {
            return value;
        }
        shift += 7;
    }
}

impl Directory {
    fn prepare(&mut self, source: &ParallelSource, byte_limit: usize) {
        let start = self.bytes.len();
        while self.bytes.len() - start < byte_limit.max(1) {
            if self.active.is_none() {
                let Some(states) = source.states.lock().unwrap().pop_front() else {
                    break;
                };
                self.active = Some((states, 0));
            }
            let (states, position) = self.active.as_mut().unwrap();
            if *position == states.len() {
                self.active = None;
                continue;
            }
            let state = &states[*position];
            if self.offsets.len() <= state.index {
                self.offsets.resize(state.index + 1, usize::MAX);
            }
            self.offsets[state.index] = self.bytes.len();
            put(state.g, &mut self.bytes);
            put(state.d.len().try_into().unwrap(), &mut self.bytes);
            for &degree in &state.d {
                put(degree, &mut self.bytes);
            }
            put(state.ch.len().try_into().unwrap(), &mut self.bytes);
            for &(index, power) in &state.ch {
                put(index, &mut self.bytes);
                put(power, &mut self.bytes);
            }
            *position += 1;
        }
        if self
            .active
            .as_ref()
            .is_some_and(|(states, position)| *position == states.len())
        {
            self.active = None;
        }
    }

    fn key(&self, id: usize) -> Option<IntegralKey> {
        let &offset = self.offsets.get(id)?;
        if offset == usize::MAX {
            return None;
        }
        let mut bytes = &self.bytes[offset..];
        let g = take(&mut bytes);
        let n = take(&mut bytes);
        let d = (0..n).map(|_| take(&mut bytes)).collect();
        let count = take(&mut bytes);
        let ch = (0..count)
            .map(|_| (take(&mut bytes), take(&mut bytes)))
            .collect();
        Some(IntegralKey::Descendant { g, d, ch })
    }

    /// Reuse the raw mathematical identity verbatim, adding only the same
    /// version/rank prefix as the full IntegralKey codec. Return the original
    /// identity's in-memory payload size so both consumers use the same soft
    /// read window; no d/ch vector is reconstructed for storage.
    fn encoded_key(&self, id: usize, output: &mut Vec<u8>) -> Option<usize> {
        let &offset = self.offsets.get(id)?;
        if offset == usize::MAX {
            return None;
        }
        let raw = &self.bytes[offset..];
        let mut remaining = raw;
        let g = take(&mut remaining);
        let n = take(&mut remaining) as usize;
        for _ in 0..n {
            take(&mut remaining);
        }
        let count = take(&mut remaining) as usize;
        let mut weight = 0u128;
        for _ in 0..count {
            let index = take(&mut remaining);
            let power = take(&mut remaining);
            weight += u128::from(index) * u128::from(power);
        }
        encode_descendant_prefix_into(g, n, weight, output);
        output.extend_from_slice(&raw[..raw.len() - remaining.len()]);
        Some(
            n.saturating_mul(size_of::<u32>())
                .saturating_add(count.saturating_mul(size_of::<(u32, u32)>())),
        )
    }
}

impl ResultSource for ParallelSource {
    fn floor(&self) -> usize {
        self.floor.load(Ordering::Acquire)
    }

    fn durable(&self, cursor: usize) {
        self.durable.fetch_max(cursor, Ordering::Release);
    }

    fn reclaim(&self) -> usize {
        let Some(_idle) = self.try_idle() else {
            return 0;
        };
        self.try_reclaim(true).unwrap_or(0)
    }

    fn reclaim_during_query(&self, registry: &Registry) -> usize {
        let Some(_idle) = self.try_idle() else {
            return 0;
        };
        let Some(bytes) = self.try_reclaim(true) else {
            return 0;
        };
        if bytes != 0 {
            // A new mathematical user cannot acquire its read permit until
            // the cleared source is known to require historical recovery.
            registry.set_history(true);
        }
        bytes
    }

    fn release(&self) {
        let Some(_idle) = self.try_idle() else {
            return;
        };
        let _ = self.try_reclaim(false);
    }

    fn heat(&self) -> u64 {
        self.reuse.score()
    }

    fn ready(&self) -> usize {
        self.ready.load(Ordering::Acquire)
    }

    fn pull(&self, cursor: usize, byte_limit: usize) -> Option<Pulled> {
        let mut directory = self.directory.try_lock().ok()?;
        if cursor < self.floor.load(Ordering::Acquire) {
            return None;
        }
        // Identity discovery also progresses while planning has not yet yielded
        // a result. No source values or planner-index lock is held here.
        directory.prepare(self, byte_limit);
        if self.writing.load(Ordering::Acquire) {
            return None;
        }
        let ready = self.ready.load(Ordering::Acquire);
        if cursor >= ready {
            return Some(Pulled {
                next: cursor,
                rows: Vec::new(),
            });
        }
        let batch = {
            let batches = self.batches.try_lock().ok()?;
            let position = batches.partition_point(|batch| batch.start <= cursor);
            Arc::clone(&batches[position - 1])
        };
        let ready = ready.min(batch.start + batch.offsets[batch.published.load(Ordering::Acquire)]);
        let values = match self.values.try_read() {
            Ok(values) => values,
            Err(TryLockError::Poisoned(error)) => error.into_inner(),
            Err(TryLockError::WouldBlock) => return None,
        };
        if self.writing.load(Ordering::Acquire) {
            return None;
        }
        let relative = cursor - batch.start;
        let layers = batch.dirty_layers.as_deref().unwrap_or(&batch.layers);
        let mut layer = batch.offsets.partition_point(|&offset| offset <= relative) - 1;
        let mut position = relative - batch.offsets[layer];
        let mut next = cursor;
        let mut bytes = 0usize;
        let mut rows = Vec::new();
        while next < ready && layer < layers.len() {
            while position < layers[layer].len() && next < ready {
                let id = batch.base + layers[layer][position];
                let Some(key) = directory.key(id) else {
                    return if rows.is_empty() {
                        None
                    } else {
                        Some(Pulled { next, rows })
                    };
                };
                let value = &values[id];
                let IntegralKey::Descendant { d, ch, .. } = &key else {
                    unreachable!()
                };
                bytes = bytes.saturating_add(value_bytes(value)).saturating_add(
                    d.len() * size_of::<u32>() + ch.len() * size_of::<(u32, u32)>(),
                );
                rows.push(CompletedIntegral {
                    key,
                    value: value.clone(),
                });
                next += 1;
                position += 1;
                if bytes >= byte_limit.max(1) {
                    return Some(Pulled { next, rows });
                }
            }
            layer += 1;
            position = 0;
        }
        Some(Pulled { next, rows })
    }

    fn pull_encoded(
        &self,
        cursor: usize,
        byte_limit: usize,
        sink: &mut dyn FnMut(&[u8], &Q),
    ) -> Option<usize> {
        let mut directory = self.directory.try_lock().ok()?;
        if cursor < self.floor.load(Ordering::Acquire) {
            return None;
        }
        directory.prepare(self, byte_limit);
        if self.writing.load(Ordering::Acquire) {
            return None;
        }
        let ready = self.ready.load(Ordering::Acquire);
        if cursor >= ready {
            return Some(cursor);
        }
        let batch = {
            let batches = self.batches.try_lock().ok()?;
            let position = batches.partition_point(|batch| batch.start <= cursor);
            Arc::clone(&batches[position - 1])
        };
        let ready = ready.min(batch.start + batch.offsets[batch.published.load(Ordering::Acquire)]);
        let values = match self.values.try_read() {
            Ok(values) => values,
            Err(TryLockError::Poisoned(error)) => error.into_inner(),
            Err(TryLockError::WouldBlock) => return None,
        };
        if self.writing.load(Ordering::Acquire) {
            return None;
        }
        let relative = cursor - batch.start;
        let layers = batch.dirty_layers.as_deref().unwrap_or(&batch.layers);
        let mut layer = batch.offsets.partition_point(|&offset| offset <= relative) - 1;
        let mut position = relative - batch.offsets[layer];
        let mut next = cursor;
        let mut bytes = 0usize;
        // One backend scratch allocation is reused across windows. Moving the
        // Vec header lets raw identities stay immutably borrowed while it grows.
        let mut encoded = std::mem::take(&mut directory.encoded_scratch);
        let result = (|| {
            while next < ready && layer < layers.len() {
                while position < layers[layer].len() && next < ready {
                    let id = batch.base + layers[layer][position];
                    let Some(identity_bytes) = directory.encoded_key(id, &mut encoded) else {
                        return if next == cursor { None } else { Some(next) };
                    };
                    let value = &values[id];
                    bytes = bytes
                        .saturating_add(value_bytes(value))
                        .saturating_add(identity_bytes);
                    // Storage encodes/copies bytes here, before the next key
                    // reuses this scratch and before this bounded read ends.
                    sink(&encoded, value);
                    next += 1;
                    position += 1;
                    if bytes >= byte_limit.max(1) {
                        return Some(next);
                    }
                }
                layer += 1;
                position = 0;
            }
            Some(next)
        })();
        directory.encoded_scratch = encoded;
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::Zero;

    fn durable_source(value: u32) -> (Arc<ParallelSource>, Arc<RwLock<Index>>) {
        let source = Arc::new(ParallelSource::default());
        let shard = Arc::new(RwLock::new(Index::new()));
        shard
            .write()
            .unwrap()
            .entry(vec![])
            .or_default()
            .insert(vec![4], 0);
        source.observe_shards(vec![Arc::clone(&shard)]);
        source.write().push(Q::from(value));
        source.observe(&Arc::new(vec![state(0, vec![4])]));
        let batch = source.begin(0, vec![vec![0]]);
        source.publish(&batch, 0);
        source.durable(1);
        (source, shard)
    }

    #[test]
    fn live_numeric_ids_exclude_reclaim_between_table_borrows() {
        let registry = Registry::new();
        let (source, shard) = durable_source(7);
        let usage = source.begin_use();
        let id = shard.read().unwrap()[&vec![]][&vec![4]];
        // No index or value guard survives here. Only the whole-use read
        // permit protects this copied numeric ID across mathematical stages.
        let reclaimed = std::thread::scope(|scope| {
            scope
                .spawn(|| source.reclaim_during_query(&registry))
                .join()
                .unwrap()
        });
        assert_eq!(reclaimed, 0);
        assert_eq!(source.reclaim(), 0);
        assert_eq!(source.read()[id], Q::from(7));
        assert!(!registry.has_history());
        drop(usage);
        assert!(source.reclaim_during_query(&registry) > 0);
        let _next_use = source.begin_use();
        assert!(registry.has_history());
        assert!(shard.read().unwrap().is_empty());
        assert!(source.read().is_empty());
        assert_eq!(source.floor(), 1);
    }

    #[test]
    fn active_registry_can_reclaim_another_idle_parallel_source() {
        let registry = Arc::new(Registry::new());
        let (active, _) = durable_source(7);
        let (idle, idle_index) = durable_source(11);
        registry.register(active.clone());
        registry.register(idle.clone());
        crate::cache_source::with_registry(registry.clone(), || {
            let _use = active.begin_use();
            assert!(registry.try_reclaim().is_none());
            let reclaimed = std::thread::scope(|scope| {
                scope
                    .spawn(|| idle.reclaim_during_query(&registry))
                    .join()
                    .unwrap()
            });
            assert!(reclaimed > 0);
            let _next_use = idle.begin_use();
            assert!(registry.has_history());
            assert!(idle_index.read().unwrap().is_empty());
            assert!(idle.read().is_empty());
            assert_eq!(active.reclaim_during_query(&registry), 0);
            assert_eq!(active.read()[0], Q::from(7));
        });
    }

    #[test]
    fn unsuccessful_reclaim_does_not_publish_a_false_history_hint() {
        let registry = Registry::new();
        let source = ParallelSource::default();
        assert_eq!(source.reclaim_during_query(&registry), 0);
        assert!(!registry.has_history());
        source.write().push(Q::from(7));
        source.observe(&Arc::new(vec![state(0, vec![4])]));
        let batch = source.begin(0, vec![vec![0]]);
        source.publish(&batch, 0);
        assert_eq!(source.reclaim_during_query(&registry), 0);
        assert_eq!(source.read()[0], Q::from(7));
        assert!(!registry.has_history());
    }

    fn encoded_rows(
        source: &ParallelSource,
        cursor: usize,
        byte_limit: usize,
    ) -> Option<(usize, Vec<(Vec<u8>, Q)>)> {
        let mut rows = Vec::new();
        let next = source.pull_encoded(cursor, byte_limit, &mut |key, value| {
            rows.push((key.to_vec(), value.clone()));
        });
        match next {
            Some(next) => {
                assert_eq!(next.checked_sub(cursor), Some(rows.len()));
                Some((next, rows))
            }
            None => {
                assert!(rows.is_empty());
                None
            }
        }
    }

    fn encoded(key: &IntegralKey) -> Vec<u8> {
        let mut bytes = Vec::new();
        crate::cache_source::encode_key_into(key, &mut bytes);
        bytes
    }

    #[test]
    fn raw_directory_encoding_matches_full_mathematical_keys() {
        let keys = [
            IntegralKey::Descendant {
                g: 2,
                d: vec![0, 4],
                ch: vec![(1, 3), (3, 1)],
            },
            IntegralKey::Descendant {
                g: 0,
                d: vec![0],
                ch: vec![],
            },
            IntegralKey::Descendant {
                g: 1,
                d: vec![],
                ch: vec![(1, 0)],
            },
            IntegralKey::Descendant {
                g: u32::MAX,
                d: vec![u32::MAX, 0, 128],
                ch: vec![(u32::MAX, u32::MAX), (u32::MAX, u32::MAX), (1, 0)],
            },
            IntegralKey::Descendant {
                g: 2,
                d: vec![0; 128],
                ch: (0..128).map(|m| (m, m % 3)).collect(),
            },
        ];
        let source = ParallelSource::default();
        source.observe(&Arc::new(
            keys.iter()
                .enumerate()
                .map(|(index, key)| {
                    let IntegralKey::Descendant { g, d, ch } = key else {
                        unreachable!()
                    };
                    State {
                        g: *g,
                        d: d.clone(),
                        ch: ch.clone(),
                        index,
                        rank: (0, 0),
                    }
                })
                .collect(),
        ));
        let mut directory = source.directory.lock().unwrap();
        directory.prepare(&source, usize::MAX);
        let mut buffer = vec![255; 128];
        for (index, key) in keys.iter().enumerate() {
            let IntegralKey::Descendant { d, ch, .. } = key else {
                unreachable!()
            };
            assert_eq!(
                directory.encoded_key(index, &mut buffer),
                Some(d.len() * size_of::<u32>() + ch.len() * size_of::<(u32, u32)>())
            );
            assert_eq!(buffer, encoded(key));
        }
        assert!(directory.encoded_key(keys.len(), &mut buffer).is_none());
    }

    #[test]
    fn encoded_window_borrows_original_values_and_preserves_old_pull_order() {
        let source = ParallelSource::default();
        source
            .write()
            .extend([Q::ZERO, Q::from_signeds(-7, 11), Q::from(13)]);
        source.observe(&Arc::new(vec![
            state(0, vec![4]),
            state(1, vec![0, 5]),
            state(2, vec![2, 3]),
        ]));
        let batch = source.begin(0, vec![vec![2], vec![0, 1]]);
        assert_eq!(encoded_rows(&source, 0, 1), Some((0, vec![])));
        source
            .directory
            .lock()
            .unwrap()
            .prepare(&source, usize::MAX);
        source.publish(&batch, 0);
        let original = source.read();
        let mut calls = 0;
        assert_eq!(
            source.pull_encoded(0, 1, &mut |key, value| {
                assert!(
                    std::ptr::eq(value, &original[2]),
                    "borrow the mathematical Q itself"
                );
                assert_eq!(
                    key,
                    encoded(&IntegralKey::Descendant {
                        g: 2,
                        d: vec![2, 3],
                        ch: vec![]
                    })
                );
                calls += 1;
            }),
            Some(1)
        );
        assert_eq!(calls, 1);
        drop(original);
        source.publish(&batch, 1);
        let mut cursor = 0;
        while cursor < source.ready() {
            let old = source.pull(cursor, 1).unwrap();
            let (next, rows) = encoded_rows(&source, cursor, 1).unwrap();
            assert_eq!(next, old.next);
            let expected = old
                .rows
                .iter()
                .map(|row| (encoded(&row.key), row.value.clone()))
                .collect::<Vec<_>>();
            assert_eq!(rows, expected);
            cursor = next;
        }
        let capacity = source.directory.lock().unwrap().encoded_scratch.capacity();
        assert!(capacity > 0);
        encoded_rows(&source, 0, usize::MAX).unwrap();
        assert_eq!(
            source.directory.lock().unwrap().encoded_scratch.capacity(),
            capacity
        );
    }

    #[test]
    fn encoded_pull_is_nonblocking_and_can_resume_missing_identity() {
        let source = ParallelSource::default();
        source.write().extend([Q::from(3), Q::from(7)]);
        source.observe(&Arc::new(vec![state(0, vec![4])]));
        let batch = source.begin(0, vec![vec![0, 1]]);
        source.publish(&batch, 0);
        let write = source.write();
        assert!(encoded_rows(&source, 0, 1024).is_none());
        drop(write);
        let directory = source.directory.lock().unwrap();
        assert!(encoded_rows(&source, 0, 1024).is_none());
        drop(directory);
        let first = encoded_rows(&source, 0, usize::MAX).unwrap();
        assert_eq!(first.0, 1);
        assert_eq!(first.1[0].1, Q::from(3));
        assert!(encoded_rows(&source, 1, 1024).is_none());
        source.observe(&Arc::new(vec![state(1, vec![0, 5])]));
        let second = encoded_rows(&source, 1, usize::MAX).unwrap();
        assert_eq!(second.0, 2);
        assert_eq!(second.1[0].1, Q::from(7));
    }

    #[test]
    fn encoded_pull_skips_restored_values_and_unpublished_tails() {
        let source = ParallelSource::default();
        source
            .write()
            .extend([Q::from(3), Q::ZERO, Q::from(7), Q::ZERO]);
        source.observe(&Arc::new(vec![
            state(0, vec![4]),
            state(1, vec![0, 5]),
            state(2, vec![2, 3]),
            state(3, vec![1, 4]),
        ]));
        let batch = source.begin_restored(0, vec![vec![2], vec![0, 1]], vec![1, 2]);
        source.publish(&batch, 0);
        assert_eq!(encoded_rows(&source, 0, usize::MAX), Some((0, vec![])));
        source.publish(&batch, 1);
        let first = encoded_rows(&source, 0, usize::MAX).unwrap();
        assert_eq!(first.0, 1);
        assert_eq!(first.1.len(), 1);
        assert_eq!(first.1[0].1, Q::from(3));
        let _interrupted = source.begin(3, vec![vec![0]]);
        assert_eq!(encoded_rows(&source, 1, usize::MAX), Some((1, vec![])));
        source.write().push(Q::from(11));
        source.observe(&Arc::new(vec![state(4, vec![0, 0, 6])]));
        let resumed = source.begin(4, vec![vec![0]]);
        source.publish(&resumed, 0);
        let next = encoded_rows(&source, 1, usize::MAX).unwrap();
        assert_eq!(next.0, 2);
        assert_eq!(next.1.len(), 1);
        assert_eq!(next.1[0].1, Q::from(11));
        assert_eq!(
            source.completed_since(0, 5),
            vec![true, true, true, false, true]
        );
    }

    #[test]
    fn encoded_pull_observes_reclaimed_floor_before_reusing_runtime_ids() {
        let source = ParallelSource::default();
        source.write().push(Q::from(7));
        source.observe(&Arc::new(vec![state(0, vec![4])]));
        let first = source.begin(0, vec![vec![0]]);
        source.publish(&first, 0);
        encoded_rows(&source, 0, usize::MAX).unwrap();
        source.durable(1);
        assert!(source.reclaim() > 0);
        assert!(encoded_rows(&source, 0, 1024).is_none());
        assert_eq!(encoded_rows(&source, 1, 1024), Some((1, vec![])));
        assert_eq!(
            source.directory.lock().unwrap().encoded_scratch.capacity(),
            0
        );
        source.write().push(Q::ZERO);
        source.observe(&Arc::new(vec![state(0, vec![0, 5])]));
        let second = source.begin(0, vec![vec![0]]);
        source.publish(&second, 0);
        let next = encoded_rows(&source, 1, 1024).unwrap();
        assert_eq!(
            next,
            (
                2,
                vec![(
                    encoded(&IntegralKey::Descendant {
                        g: 2,
                        d: vec![0, 5],
                        ch: vec![]
                    }),
                    Q::ZERO
                )]
            )
        );
    }

    fn state(index: usize, d: Vec<u32>) -> State {
        State {
            g: 2,
            d,
            ch: Vec::new(),
            rank: (0, 0),
            index,
        }
    }

    #[test]
    fn completion_preserves_zero_and_out_of_order_ids_across_queries() {
        let source = ParallelSource::default();
        source.write().resize(3, Q::ZERO);
        source.observe(&Arc::new(vec![
            state(0, vec![4]),
            state(1, vec![0, 5]),
            state(2, vec![2, 3]),
        ]));
        let first = source.begin(0, vec![vec![2], vec![0, 1]]);
        assert!(source.pull(0, 65536).unwrap().rows.is_empty());
        source.write()[2] = Q::from(7);
        source.publish(&first, 0);
        let pulled = source.pull(0, 65536).unwrap();
        assert_eq!(pulled.rows.len(), 1);
        assert_eq!(pulled.rows[0].value, Q::from(7));
        assert_eq!(
            pulled.rows[0].key,
            IntegralKey::Descendant {
                g: 2,
                d: vec![2, 3],
                ch: vec![]
            }
        );
        source.write()[0] = Q::from(3);
        source.publish(&first, 1);
        let pulled = source.pull(1, 65536).unwrap();
        assert_eq!(pulled.rows.len(), 2);
        assert_eq!(pulled.rows[1].value, Q::ZERO);
        assert_eq!(
            pulled.rows[1].key,
            IntegralKey::Descendant {
                g: 2,
                d: vec![0, 5],
                ch: vec![]
            }
        );
        source.write().push(Q::from(11));
        source.observe(&Arc::new(vec![state(3, vec![1, 4])]));
        let second = source.begin(3, vec![vec![0]]);
        source.publish(&second, 0);
        assert_eq!(source.pull(3, 1).unwrap().rows[0].value, Q::from(11));
        assert!(source.pull(4, 65536).unwrap().rows.is_empty());
    }

    #[test]
    fn write_intent_is_released_after_unwind() {
        let source = ParallelSource::default();
        let _ = std::panic::catch_unwind(|| {
            let _write = source.write();
            panic!("simulate interrupted writeback");
        });
        assert!(!source.writing.load(Ordering::Acquire));
        assert!(source.read().is_empty());
    }

    #[test]
    fn interrupted_batch_does_not_publish_its_unfinished_tail() {
        let source = ParallelSource::default();
        source.write().resize(4, Q::ZERO);
        source.observe(&Arc::new(vec![
            state(0, vec![4]),
            state(1, vec![0, 5]),
            state(2, vec![2, 3]),
            state(3, vec![1, 4]),
        ]));
        let interrupted = source.begin(0, vec![vec![2], vec![0, 1]]);
        source.write()[2] = Q::from(7);
        source.publish(&interrupted, 0);
        let resumed = source.begin(3, vec![vec![0]]);
        source.write()[3] = Q::from(11);
        source.publish(&resumed, 0);
        let first = source.pull(0, 65536).unwrap();
        assert_eq!(first.next, 1);
        assert_eq!(first.rows.len(), 1);
        assert_eq!(first.rows[0].value, Q::from(7));
        let second = source.pull(first.next, 65536).unwrap();
        assert_eq!(second.rows[0].value, Q::from(11));
        assert_eq!(second.next, 2);
        assert_eq!(source.completed_since(0, 4), vec![false, false, true, true]);
    }

    #[test]
    fn restored_values_complete_mathematics_without_becoming_dirty_again() {
        let source = ParallelSource::default();
        source.write().resize(3, Q::ZERO);
        source.observe(&Arc::new(vec![
            state(0, vec![4]),
            state(1, vec![0, 5]),
            state(2, vec![2, 3]),
        ]));
        let batch = source.begin_restored(0, vec![vec![2], vec![0, 1]], vec![1, 2]);
        source.write()[2] = Q::from(7);
        source.publish(&batch, 0);
        assert_eq!(source.ready(), 0);
        assert_eq!(source.completed_since(0, 3), vec![false, false, true]);
        source.write()[0] = Q::from(3);
        source.publish(&batch, 1);
        assert_eq!(source.ready(), 1);
        assert_eq!(source.completed_since(0, 3), vec![true, true, true]);
        let dirty = source.pull(0, usize::MAX).unwrap();
        assert_eq!(dirty.next, 1);
        assert_eq!(dirty.rows.len(), 1);
        assert_eq!(dirty.rows[0].value, Q::from(3));
        assert_eq!(source.read()[1], Q::ZERO);
        assert!(source.pull(1, usize::MAX).unwrap().rows.is_empty());
    }

    #[test]
    fn reclamation_clears_original_shards_and_keeps_the_completion_cursor_monotone() {
        let source = ParallelSource::default();
        let shard = Arc::new(RwLock::new(Index::new()));
        shard
            .write()
            .unwrap()
            .entry(vec![])
            .or_default()
            .insert(vec![4], 0);
        source.observe_shards(vec![Arc::clone(&shard)]);
        source.write().push(Q::from(7));
        source.observe(&Arc::new(vec![state(0, vec![4])]));
        let first = source.begin(0, vec![vec![0]]);
        source.publish(&first, 0);
        source.note_query(1, 4);
        assert_eq!(source.heat(), 250_000);
        assert_eq!(source.reclaim(), 0); // Dirty values cannot be reclaimed.
        assert_eq!(source.heat(), 250_000);
        assert_eq!(source.read()[0], Q::from(7));
        source.durable(1);
        let held = shard.read().unwrap();
        assert_eq!(source.reclaim(), 0); // Failed try-lock does not clear Q.
        assert_eq!(source.read()[0], Q::from(7));
        drop(held);
        assert!(source.reclaim() > 0);
        assert_eq!(source.heat(), 0);
        assert!(shard.read().unwrap().is_empty());
        assert_eq!(source.read().capacity(), 0);
        assert_eq!(source.ready(), 1);
        assert_eq!(source.floor(), 1);
        assert!(source.pull(0, 1024).is_none());
        assert!(source.pull(1, 1024).unwrap().rows.is_empty());

        // Runtime ID zero is now available again, but its old key is gone.
        shard
            .write()
            .unwrap()
            .entry(vec![])
            .or_default()
            .insert(vec![0, 5], 0);
        source.write().push(Q::ZERO);
        source.observe(&Arc::new(vec![state(0, vec![0, 5])]));
        let second = source.begin(0, vec![vec![0]]);
        source.publish(&second, 0);
        let next = source.pull(1, 1024).unwrap();
        assert_eq!(next.next, 2);
        assert_eq!(next.rows[0].value, Q::ZERO);
        assert_eq!(
            next.rows[0].key,
            IntegralKey::Descendant {
                g: 2,
                d: vec![0, 5],
                ch: vec![]
            }
        );
        assert_eq!(source.reclaim(), 0);
        source.release();
        assert_eq!(source.floor(), 1);
        assert_eq!(source.read().len(), 1);
        source.durable(2);
        source.release();
        assert_eq!(source.floor(), 2);
        assert_eq!(source.ready(), 2);
        assert!(shard.read().unwrap().is_empty());
        assert_eq!(source.read().capacity(), 0);
    }

    #[test]
    fn unpublished_panic_holes_do_not_become_zero_cache_entries_after_reclaim() {
        let source = ParallelSource::default();
        let shard = Arc::new(RwLock::new(Index::new()));
        shard
            .write()
            .unwrap()
            .entry(vec![])
            .or_default()
            .extend([(vec![4], 0), (vec![0, 5], 1)]);
        source.observe_shards(vec![Arc::clone(&shard)]);
        source.write().extend([Q::from(3), Q::ZERO]);
        source.observe(&Arc::new(vec![state(0, vec![4]), state(1, vec![0, 5])]));
        let batch = source.begin(0, vec![vec![0], vec![1]]);
        source.publish(&batch, 0);
        source.durable(1);
        assert!(source.reclaim() > 0);
        assert_eq!(source.floor(), 1);
        assert_eq!(source.ready(), 1);
        assert!(shard.read().unwrap().is_empty());
        assert!(source.read().is_empty());
    }
}
