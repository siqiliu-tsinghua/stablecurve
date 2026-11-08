//! Completed mathematical values exposed to an optional storage consumer.
//!
//! Sources retain their original values. Only the consumer calls `pull`, so
//! copying, key construction, and storage batching stay off the compute path.

use sc_num::Q;
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock, Weak};

pub(crate) mod budget;
mod codec;
mod compat;
pub(crate) mod legacy_source;
mod manual;
mod reader;
pub(crate) mod reuse;
pub(crate) mod sampler;
pub(crate) mod storage;
pub(crate) mod work;
mod work_storage;

pub(crate) use codec::encode_descendant_prefix_into;
#[cfg(test)]
pub(crate) use codec::encode_key_into;
pub(crate) use reader::Reader;

/// An allocation-boundary failure crossing the original infallible internal
/// mathematics. Public Result entries catch only this payload; genuine bugs
/// retain their ordinary panic behavior.
#[derive(Debug)]
pub(crate) struct ResourceFailure(pub(crate) String);

pub(crate) fn resource_failure(message: impl Into<String>) -> ! {
    // resume_unwind avoids printing a panic diagnostic for an expected resource
    // error which the public Result boundary reports to its caller.
    std::panic::resume_unwind(Box::new(ResourceFailure(message.into())))
}

pub(crate) fn catch_resource<T>(
    calculate: impl FnOnce() -> Result<T, String>,
) -> Result<T, String> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(calculate)) {
        Ok(result) => result,
        Err(payload) => match payload.downcast::<ResourceFailure>() {
            Ok(error) => Err(error.0),
            Err(payload) => std::panic::resume_unwind(payload),
        },
    }
}

/// Exact mathematical identity; no hash or query-local node ID is persisted.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(crate) enum IntegralKey {
    Descendant {
        g: u32,
        d: Vec<u32>,
        ch: Vec<(u32, u32)>,
    },
    Integral {
        g: u32,
        n: u32,
        kappa: Vec<(u32, u32)>,
        psi: Vec<(u32, u32)>,
        lambda: Vec<(u32, u32)>,
        ch: Vec<(u32, u32)>,
    },
    /// Legacy ST includes auxiliary recursive conventions. Keep its namespace
    /// separate from normalized C and actual stable-space integrals.
    Legacy { kind: String, nums: Vec<u32> },
}

#[derive(Clone, Debug)]
pub(crate) struct CompletedIntegral {
    pub(crate) key: IntegralKey,
    pub(crate) value: Q,
}

#[derive(Debug)]
pub(crate) struct Pulled {
    pub(crate) next: usize,
    pub(crate) rows: Vec<CompletedIntegral>,
}

/// Incremental, nonblocking access to already completed values.
///
/// `None` means the source is busy; retry the same cursor later. An empty
/// batch means no new values are currently visible. Cursors belong to one
/// source and advance only over completed values. `byte_limit` is a soft
/// copy-window target: a source may return one larger row to make progress.
/// Returned values are consumer-owned copies, independent of source locks.
pub(crate) trait ResultSource: Send + Sync {
    /// Published, monotonically increasing completion cursor. Read this once
    /// to establish a finite target; future completions do not extend it.
    fn ready(&self) -> usize;

    /// Already durable and reclaimed prefix; new consumers begin here.
    fn floor(&self) -> usize {
        0
    }

    fn pull(&self, cursor: usize, byte_limit: usize) -> Option<Pulled>;

    /// Storage-only access. Each callback must consume the encoded key and
    /// borrowed value before returning; it must not wait for SQL, budget
    /// progress, or another source. One callback advances one completion.
    /// Overrides may keep the original value table borrowed for this bounded
    /// window, avoiding a temporary mathematical-value copy.
    fn pull_encoded(
        &self,
        cursor: usize,
        byte_limit: usize,
        sink: &mut dyn FnMut(&[u8], &Q),
    ) -> Option<usize> {
        let batch = self.pull(cursor, byte_limit)?;
        let mut key = Vec::new();
        for row in &batch.rows {
            codec::encode_key_into(&row.key, &mut key);
            sink(&key, &row.value);
        }
        Some(batch.next)
    }

    /// Called by storage only after the corresponding transaction commits.
    fn durable(&self, _cursor: usize) {}

    /// Pressure-only, nonblocking reclamation of an entirely durable source.
    /// The caller must hold its registry's quiescent reclamation guard.
    fn reclaim(&self) -> usize {
        0
    }

    /// Pressure-only opt-in while other queries remain active. Implementations
    /// must exclude their own borrowers and recheck durability under the source
    /// lock. Publish registry history before cleared storage becomes observable.
    /// Numeric-ID sources also need a whole-computation use permit.
    fn reclaim_during_query(&self, _registry: &Registry) -> usize {
        0
    }

    /// Release an entirely durable source when payload-size statistics are
    /// unnecessary. The same quiescence and durability requirements apply.
    fn release(&self) {
        let _ = self.reclaim();
    }

    /// Coarse query/batch reuse information; never a per-term counter.
    fn heat(&self) -> u64 {
        0
    }
}

pub(crate) trait LookupBatch: Send + Sync {
    fn len(&self) -> usize;
    /// Internal frontier contract: every key is a Descendant. This permits a
    /// read snapshot without entering Integral's lazy legacy-import writes.
    fn descendants_only(&self) -> bool {
        false
    }
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn key(&self, index: usize) -> IntegralKey;
}

const CLOSING: usize = 1usize << (usize::BITS - 1);
const RECLAIMING: usize = CLOSING >> 1;
const ACTIVE: usize = RECLAIMING - 1;

static NEXT_REGISTRY_ID: AtomicU64 = AtomicU64::new(1);

/// Strong ownership keeps completed values readable after their worker exits.
pub(crate) struct Registry {
    id: u64,
    sources: Mutex<Vec<Arc<dyn ResultSource>>>,
    storage: OnceLock<storage::StorageHandle>,
    history: AtomicBool,
    descendant_history: AtomicBool,
    failed: AtomicBool,
    error: Mutex<Option<String>>,
    manual: OnceLock<Weak<manual::ManualSource>>,
    activity: AtomicUsize,
    idle: Mutex<()>,
    idle_changed: Condvar,
}

impl Registry {
    pub(crate) fn new() -> Self {
        let mut id = NEXT_REGISTRY_ID.load(Ordering::Relaxed);
        loop {
            let next = id
                .checked_add(1)
                .expect("result source registry IDs exhausted");
            match NEXT_REGISTRY_ID.compare_exchange_weak(
                id,
                next,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(current) => id = current,
            }
        }
        Self {
            id,
            sources: Mutex::new(Vec::new()),
            storage: OnceLock::new(),
            history: AtomicBool::new(false),
            descendant_history: AtomicBool::new(false),
            failed: AtomicBool::new(false),
            error: Mutex::new(None),
            manual: OnceLock::new(),
            activity: AtomicUsize::new(0),
            idle: Mutex::new(()),
            idle_changed: Condvar::new(),
        }
    }

    pub(crate) fn id(&self) -> u64 {
        self.id
    }

    pub(crate) fn register(&self, source: Arc<dyn ResultSource>) {
        let activity = self.activity.load(Ordering::Acquire);
        assert!(
            activity & CLOSING == 0 || activity & ACTIVE != 0,
            "cache session is closed"
        );
        self.sources.lock().unwrap().push(source);
    }

    pub(crate) fn bind_storage(&self, storage: storage::StorageHandle) {
        assert!(
            self.storage.set(storage).is_ok(),
            "storage already attached to registry"
        );
    }

    pub(crate) fn set_history(&self, present: bool) {
        // Reclamation and callers without a namespace summary conservatively
        // permit descendant recovery. The SQLite owner may narrow this once,
        // at open, after proving the saved rows cannot serve planner states.
        self.descendant_history.store(present, Ordering::Release);
        self.history.store(present, Ordering::Release);
    }

    pub(crate) fn has_history(&self) -> bool {
        self.history.load(Ordering::Acquire)
    }

    pub(crate) fn set_descendant_history(&self, present: bool) {
        self.descendant_history.store(present, Ordering::Release);
    }

    pub(crate) fn has_descendant_history(&self) -> bool {
        self.descendant_history.load(Ordering::Acquire)
    }

    /// Called only before building a new DAG, after mathematical L1 misses.
    /// A narrowed opening snapshot needs one owner-side external-write check;
    /// ordinary resident queries and known relevant history require no request.
    pub(crate) fn query_descendant_history(&self) -> Result<bool, String> {
        let relevant = self.has_descendant_history();
        if !self.has_history() || relevant {
            return Ok(relevant);
        }
        self.check_error()?;
        match self.storage.get() {
            Some(storage) => storage.descendant_history(),
            None => Ok(false),
        }
    }

    pub(crate) fn fail(&self, message: String) {
        let mut error = self.error.lock().unwrap_or_else(|error| error.into_inner());
        if error.is_none() {
            *error = Some(message);
        }
        self.failed.store(true, Ordering::Release);
    }

    pub(crate) fn check_error(&self) -> Result<(), String> {
        if self.failed.load(Ordering::Acquire) {
            return Err(self
                .error
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .clone()
                .unwrap_or_else(|| "cache storage failed".into()));
        }
        if let Some(error) = self.storage.get().and_then(storage::StorageHandle::error) {
            return Err(error);
        }
        Ok(())
    }

    pub(crate) fn lookup(&self, key: IntegralKey) -> Result<Option<Q>, String> {
        self.check_error()?;
        match self.storage.get() {
            Some(storage) => storage.get(key),
            None => Ok(None),
        }
    }

    pub(crate) fn manual_get(&self, kind: &str, nums: &[u32]) -> Option<Q> {
        self.manual
            .get()
            .and_then(Weak::upgrade)
            .and_then(|source| source.get(kind, nums))
    }

    pub(crate) fn manual_put(&self, kind: &str, nums: &[u32], value: &Q) {
        self.manual
            .get_or_init(|| {
                let source = Arc::new(manual::ManualSource::default());
                self.register(source.clone());
                Arc::downgrade(&source)
            })
            .upgrade()
            .expect("manual source owner is closed")
            .put(kind, nums, value);
    }

    /// The list lock is released before a consumer accesses any source.
    pub(crate) fn sources(&self) -> Vec<Arc<dyn ResultSource>> {
        self.sources.lock().unwrap().clone()
    }

    /// Append-only discovery: an existing consumer need not clone old sources.
    fn sources_since(&self, cursor: usize) -> Vec<Arc<dyn ResultSource>> {
        self.sources.lock().unwrap()[cursor..].to_vec()
    }

    /// Experimental C1 consumer; attachments never start it implicitly.
    pub(crate) fn start_reader(self: &Arc<Self>) -> std::io::Result<Reader> {
        Reader::start(Arc::clone(self))
    }

    pub(crate) fn lookup_batch(
        &self,
        batch: Arc<dyn LookupBatch>,
    ) -> Result<Vec<Option<Q>>, String> {
        self.check_error()?;
        if batch.is_empty() {
            return Ok(Vec::new());
        }
        match self.storage.get() {
            Some(storage) => storage.lookup_batch(batch),
            None => Ok(vec![None; batch.len()]),
        }
    }

    fn enter(self: &Arc<Self>) -> Result<QueryPermit, String> {
        let mut activity = self.activity.load(Ordering::Acquire);
        loop {
            if activity & CLOSING != 0 {
                return Err("cache session is closing or closed".to_owned());
            }
            if activity & RECLAIMING != 0 {
                self.wait_reclaim();
                activity = self.activity.load(Ordering::Acquire);
                continue;
            }
            if activity & ACTIVE == ACTIVE {
                return Err("too many active cache queries".to_owned());
            }
            match self.activity.compare_exchange_weak(
                activity,
                activity + 1,
                Ordering::AcqRel,
                Ordering::Acquire,
            ) {
                Ok(_) => {
                    ACTIVE_OWNERS.with(|owners| owners.borrow_mut().push(self.id()));
                    return Ok(QueryPermit(Arc::clone(self)));
                }
                Err(current) => activity = current,
            }
        }
    }

    pub(crate) fn begin_close(&self) {
        self.activity.fetch_or(CLOSING, Ordering::AcqRel);
    }

    pub(crate) fn is_closing(&self) -> bool {
        self.activity.load(Ordering::Acquire) & CLOSING != 0
    }

    pub(crate) fn wait_idle(&self) {
        let mut idle = self.idle.lock().unwrap_or_else(|error| error.into_inner());
        while self.activity.load(Ordering::Acquire) & (ACTIVE | RECLAIMING) != 0 {
            idle = self
                .idle_changed
                .wait(idle)
                .unwrap_or_else(|error| error.into_inner());
        }
    }

    fn wait_reclaim(&self) {
        let mut idle = self.idle.lock().unwrap_or_else(|error| error.into_inner());
        while self.activity.load(Ordering::Acquire) & RECLAIMING != 0 {
            idle = self
                .idle_changed
                .wait(idle)
                .unwrap_or_else(|error| error.into_inner());
        }
    }

    pub(crate) fn try_reclaim(&self) -> Option<ReclaimGuard<'_>> {
        self.activity
            .compare_exchange(0, RECLAIMING, Ordering::AcqRel, Ordering::Acquire)
            .ok()
            .map(|_| ReclaimGuard(self))
    }

    pub(crate) fn register_workspace(
        &self,
        source: Arc<dyn work::WorkSource>,
    ) -> Result<(), String> {
        self.check_error()?;
        self.storage
            .get()
            .ok_or("plan workspace requires an attached SQLite service")?
            .register_workspace(source)
    }

    pub(crate) fn release_sources(&self) {
        self.sources
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clear();
    }
}

pub(crate) struct ReclaimGuard<'a>(&'a Registry);

impl Drop for ReclaimGuard<'_> {
    fn drop(&mut self) {
        self.0.activity.fetch_and(!RECLAIMING, Ordering::Release);
        let _idle = self
            .0
            .idle
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        self.0.idle_changed.notify_all();
    }
}

struct QueryPermit(Arc<Registry>);

impl Drop for QueryPermit {
    fn drop(&mut self) {
        ACTIVE_OWNERS.with(|owners| {
            let _ = owners.borrow_mut().pop();
        });
        let previous = self.0.activity.fetch_sub(1, Ordering::AcqRel);
        if previous & ACTIVE == 1 {
            let _idle = self
                .0
                .idle
                .lock()
                .unwrap_or_else(|error| error.into_inner());
            self.0.idle_changed.notify_all();
        }
    }
}

impl Default for Registry {
    fn default() -> Self {
        Self::new()
    }
}

thread_local! {
    // Outer None: inherit the legacy process-global attachment at query entry.
    // Some(None): explicitly disable registration in this scope.
    static SCOPED_REGISTRY: RefCell<Option<Option<Arc<Registry>>>> = const {
        RefCell::new(None)
    };
    // Scope bindings can temporarily mask an outer owner. These are the
    // actual per-thread permits, needed to avoid closing our own outer query.
    static ACTIVE_OWNERS: RefCell<Vec<u64>> = const { RefCell::new(Vec::new()) };
}

/// Resolve once at a query entry, then retain/bind that owner for the query.
/// The legacy fallback is process-global; this does not provide instance
/// isolation for callers that have not explicitly bound a registry.
pub(crate) fn current_registry() -> Option<Arc<Registry>> {
    match SCOPED_REGISTRY.with(|slot| slot.borrow().clone()) {
        Some(registry) => registry,
        None => crate::socle::result_registry(),
    }
}

/// Borrow the current owner for an internal mathematical leaf entry.
///
/// A scoped owner needs no Arc increment or decrement. The callback may nest
/// other read-only owner lookups or scopes with the same owner. Changing the
/// owner is forbidden while its TLS read borrow lasts through the callback.
/// Public entries that establish or replace a binding must use the owned
/// `current_registry` interface instead. An unbound legacy fallback is captured
/// as an owned Arc and kept alive until the callback returns.
pub(crate) fn with_current_registry<T>(f: impl FnOnce(Option<&Arc<Registry>>) -> T) -> T {
    SCOPED_REGISTRY.with(|slot| {
        let binding = slot.borrow();
        match binding.as_ref() {
            Some(registry) => f(registry.as_ref()),
            None => {
                drop(binding);
                let registry = crate::socle::result_registry();
                // Fix the legacy owner once for the entire mathematical call,
                // including nested public entries and an explicit no-cache call.
                with_optional_registry(registry, || with_current_registry(f))
            }
        }
    })
}

struct Scope {
    previous: Option<Option<Arc<Registry>>>,
    _permit: Option<QueryPermit>,
}

impl Drop for Scope {
    fn drop(&mut self) {
        SCOPED_REGISTRY.with(|slot| {
            *slot.borrow_mut() = self.previous.take();
        });
    }
}

/// Bind an owner without serializing independent calling threads.
#[cfg(test)]
pub(crate) fn with_registry<T>(registry: Arc<Registry>, f: impl FnOnce() -> T) -> T {
    with_optional_registry(Some(registry), f)
}

pub(crate) fn in_scope(id: u64) -> bool {
    SCOPED_REGISTRY.with(|slot| {
        slot.borrow()
            .as_ref()
            .and_then(Option::as_ref)
            .is_some_and(|owner| owner.id() == id)
    })
}

pub(crate) fn has_active_scope(id: u64) -> bool {
    ACTIVE_OWNERS.with(|owners| owners.borrow().contains(&id))
}

pub(crate) fn try_with_registry<T>(
    registry: Arc<Registry>,
    f: impl FnOnce() -> T,
) -> Result<T, String> {
    if in_scope(registry.id()) {
        return Ok(f());
    }
    let permit = if has_active_scope(registry.id()) {
        None
    } else {
        Some(registry.enter()?)
    };
    let previous = SCOPED_REGISTRY.with(|slot| slot.replace(Some(Some(registry))));
    let _scope = Scope {
        previous,
        _permit: permit,
    };
    Ok(f())
}

/// Borrowed permission to inherit an already accepted query in a synchronous
/// scoped pool. Its lifetime cannot outlive with_inherited_registry's callback.
pub(crate) struct InheritedRegistry<'scope> {
    registry: Option<&'scope Arc<Registry>>,
}

impl InheritedRegistry<'_> {
    pub(crate) fn run<T>(&self, calculate: impl FnOnce() -> T) -> T {
        let bound = SCOPED_REGISTRY.with(|slot| {
            let current = slot.borrow();
            match (current.as_ref(), self.registry) {
                (Some(None), None) => true,
                (Some(Some(current)), Some(owner)) => current.id() == owner.id(),
                _ => false,
            }
        });
        if bound {
            return calculate();
        }

        struct InheritedScope {
            previous: Option<Option<Arc<Registry>>>,
            owner: Option<u64>,
        }
        impl Drop for InheritedScope {
            fn drop(&mut self) {
                SCOPED_REGISTRY.with(|slot| {
                    *slot.borrow_mut() = self.previous.take();
                });
                if let Some(owner) = self.owner {
                    ACTIVE_OWNERS.with(|owners| {
                        let previous = owners.borrow_mut().pop();
                        debug_assert_eq!(previous, Some(owner));
                    });
                }
            }
        }
        let owner = self.registry.map(|owner| owner.id());
        let previous = SCOPED_REGISTRY.with(|slot| slot.replace(Some(self.registry.cloned())));
        let mut scope = InheritedScope {
            previous,
            owner: None,
        };
        if let Some(owner) = owner {
            // The caller's permit covers these scoped children, including
            // children starting after CLOSING. No new admission or count is
            // needed; the marker also prevents a child closing its own query.
            ACTIVE_OWNERS.with(|owners| owners.borrow_mut().push(owner));
            scope.owner = Some(owner);
        }
        let _scope = scope;
        calculate()
    }
}

/// Capture and accept one owner before starting a synchronous scoped pool.
/// The higher-ranked callback prevents its borrowed token escaping after the
/// parent permit is released. Explicit None is inherited as an actual binding.
pub(crate) fn with_inherited_registry<T>(
    calculate: impl for<'scope> FnOnce(InheritedRegistry<'scope>) -> Result<T, String>,
) -> Result<T, String> {
    let registry = current_registry();
    match registry.as_ref() {
        Some(owner) => try_with_registry(owner.clone(), || {
            calculate(InheritedRegistry {
                registry: Some(owner),
            })
        })?,
        None => with_optional_registry(None, || calculate(InheritedRegistry { registry: None })),
    }
}

/// Explicit `None` masks the process-global fallback, including in nested calls.
pub(crate) fn with_optional_registry<T>(
    registry: Option<Arc<Registry>>,
    f: impl FnOnce() -> T,
) -> T {
    let already_bound = SCOPED_REGISTRY.with(|slot| {
        let binding = slot.borrow();
        match (binding.as_ref(), registry.as_ref()) {
            (Some(None), None) => true,
            (Some(Some(current)), Some(registry)) => current.id() == registry.id(),
            _ => false,
        }
    });
    if already_bound {
        return f();
    }
    let (registry, permit) = match registry {
        Some(registry) if has_active_scope(registry.id()) => (Some(registry), None),
        Some(registry) => match registry.enter() {
            Ok(permit) => (Some(registry), Some(permit)),
            // Infallible internal helpers retain their mathematical behavior
            // after closure. Result-based external entries use try_with_registry.
            Err(_) => (None, None),
        },
        None => (None, None),
    };
    let previous = SCOPED_REGISTRY.with(|slot| slot.replace(Some(registry)));
    let _scope = Scope {
        previous,
        _permit: permit,
    };
    f()
}

/// Inline rational storage plus effective heap limbs, excluding allocation
/// spare capacity and allocator metadata. This sizes soft pull windows, not
/// a hard memory budget or a process RSS guarantee.
pub(crate) fn value_bytes(value: &Q) -> usize {
    let limb_size = value
        .denominator_ref()
        .limbs()
        .next()
        .map_or(8, |limb| std::mem::size_of_val(&limb));
    [value.numerator_ref(), value.denominator_ref()]
        .into_iter()
        .fold(std::mem::size_of::<Q>(), |bytes, natural| {
            let count = usize::try_from(natural.limb_count()).unwrap_or(usize::MAX);
            let heap = if count > 1 {
                count.saturating_mul(limb_size)
            } else {
                0
            };
            bytes.saturating_add(heap)
        })
}

/// Effective row storage, excluding Vec spare capacity and allocator metadata.
pub(crate) fn record_bytes(row: &CompletedIntegral) -> usize {
    let payload = match &row.key {
        IntegralKey::Descendant { d, ch, .. } => d
            .len()
            .saturating_mul(std::mem::size_of::<u32>())
            .saturating_add(ch.len().saturating_mul(std::mem::size_of::<(u32, u32)>())),
        IntegralKey::Integral {
            kappa,
            psi,
            lambda,
            ch,
            ..
        } => [kappa, psi, lambda, ch]
            .into_iter()
            .fold(0usize, |bytes, entries| {
                bytes.saturating_add(
                    entries
                        .len()
                        .saturating_mul(std::mem::size_of::<(u32, u32)>()),
                )
            }),
        IntegralKey::Legacy { kind, nums } => kind
            .len()
            .saturating_add(nums.len().saturating_mul(size_of::<u32>())),
    };
    std::mem::size_of::<CompletedIntegral>()
        .saturating_add(payload)
        .saturating_add(value_bytes(&row.value).saturating_sub(std::mem::size_of::<Q>()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    #[test]
    fn default_encoded_pull_preserves_namespaces_and_completion_count() {
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
        let rows = Rows(vec![
            CompletedIntegral {
                key: IntegralKey::Descendant {
                    g: 2,
                    d: vec![0, 5],
                    ch: vec![],
                },
                value: Q::from(0),
            },
            CompletedIntegral {
                key: IntegralKey::Integral {
                    g: 2,
                    n: 0,
                    kappa: vec![],
                    psi: vec![],
                    lambda: vec![],
                    ch: vec![(1, 3)],
                },
                value: Q::from_signeds(-7, 11),
            },
            CompletedIntegral {
                key: IntegralKey::Legacy {
                    kind: "st".into(),
                    nums: vec![],
                },
                value: Q::from(1),
            },
        ]);
        let mut encoded = Vec::new();
        let next = rows.pull_encoded(0, usize::MAX, &mut |key, value| {
            encoded.push((key.to_vec(), value.clone()));
        });
        assert_eq!(next, Some(encoded.len()));
        for (actual, original) in encoded.iter().zip(&rows.0) {
            let mut expected = Vec::new();
            codec::encode_key_into(&original.key, &mut expected);
            assert_eq!(actual, &(expected, original.value.clone()));
        }
        let next = rows.pull_encoded(rows.ready(), 1, &mut |_, _| panic!("caught-up callback"));
        assert_eq!(next, Some(rows.ready()));
    }

    struct Source {
        dropped: Arc<AtomicUsize>,
    }

    impl ResultSource for Source {
        fn ready(&self) -> usize {
            0
        }

        fn pull(&self, cursor: usize, _byte_limit: usize) -> Option<Pulled> {
            Some(Pulled {
                next: cursor,
                rows: Vec::new(),
            })
        }
    }

    impl Drop for Source {
        fn drop(&mut self) {
            self.dropped.fetch_add(1, Ordering::Relaxed);
        }
    }

    #[test]
    fn registry_owns_sources_after_producer_exit() {
        let registry = Registry::new();
        let dropped = Arc::new(AtomicUsize::new(0));
        let worker_registry = Arc::new(registry);
        let worker = Arc::clone(&worker_registry);
        let worker_dropped = Arc::clone(&dropped);
        std::thread::spawn(move || {
            worker.register(Arc::new(Source {
                dropped: worker_dropped,
            }));
        })
        .join()
        .unwrap();
        assert_eq!(dropped.load(Ordering::Relaxed), 0);
        let sources = worker_registry.sources();
        assert_eq!(sources.len(), 1);
        assert_eq!(sources[0].pull(7, 1024).unwrap().next, 7);
        drop(worker_registry);
        assert_eq!(dropped.load(Ordering::Relaxed), 0);
        drop(sources);
        assert_eq!(dropped.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn nested_scope_restores_owner_and_explicit_none() {
        let first = Arc::new(Registry::new());
        let second = Arc::new(Registry::default());
        assert_ne!(first.id(), second.id());
        with_registry(Arc::clone(&first), || {
            assert_eq!(current_registry().unwrap().id(), first.id());
            with_optional_registry(None, || {
                assert!(current_registry().is_none());
                with_registry(Arc::clone(&second), || {
                    assert_eq!(current_registry().unwrap().id(), second.id());
                });
                assert!(current_registry().is_none());
            });
            assert_eq!(current_registry().unwrap().id(), first.id());
        });
    }

    #[test]
    fn unwinding_restores_outer_scope() {
        let registry = Arc::new(Registry::new());
        with_registry(Arc::clone(&registry), || {
            let result = std::panic::catch_unwind(|| {
                with_optional_registry(None, || panic!("scope test"));
            });
            assert!(result.is_err());
            assert_eq!(current_registry().unwrap().id(), registry.id());
        });
    }

    #[test]
    fn borrowed_owner_respects_explicit_none_and_nested_read_scopes() {
        let registry = Arc::new(Registry::new());
        with_optional_registry(None, || {
            with_current_registry(|owner| {
                assert!(owner.is_none());
                with_current_registry(|nested| assert!(nested.is_none()));
            });
        });
        with_registry(Arc::clone(&registry), || {
            let strong_count = Arc::strong_count(&registry);
            with_current_registry(|owner| {
                let owner = owner.unwrap();
                assert!(Arc::ptr_eq(owner, &registry));
                assert_eq!(Arc::strong_count(&registry), strong_count);
                with_current_registry(|nested| {
                    assert!(Arc::ptr_eq(nested.unwrap(), owner));
                    assert_eq!(Arc::strong_count(&registry), strong_count);
                });
                assert_eq!(Arc::strong_count(&registry), strong_count);
            });
            assert_eq!(Arc::strong_count(&registry), strong_count);
        });
    }

    #[test]
    fn borrowed_owner_releases_tls_read_access_after_unwind() {
        let registry = Arc::new(Registry::new());
        with_registry(Arc::clone(&registry), || {
            let result = std::panic::catch_unwind(|| {
                with_current_registry(|_| panic!("borrowed owner scope test"));
            });
            assert!(result.is_err());
            with_optional_registry(None, || {
                with_current_registry(|owner| assert!(owner.is_none()));
            });
            assert_eq!(current_registry().unwrap().id(), registry.id());
        });
    }

    #[test]
    fn inherited_pool_keeps_main_thread_read_borrows_and_explicit_none() {
        for registry in [None, Some(Arc::new(Registry::new()))] {
            with_optional_registry(registry.clone(), || {
                with_current_registry(|borrowed| {
                    let expected = borrowed.map(|owner| owner.id());
                    with_inherited_registry(|context| {
                        let barrier = std::sync::Barrier::new(4);
                        let parts = sc_num::parallel::chunks(
                            4,
                            1,
                            4,
                            || (),
                            |_, _| {
                                barrier.wait();
                                context.run(|| {
                                    assert_eq!(
                                        current_registry().map(|owner| owner.id()),
                                        expected
                                    );
                                    with_current_registry(|nested| {
                                        assert_eq!(nested.map(|owner| owner.id()), expected);
                                    });
                                })
                            },
                        );
                        assert_eq!(parts.len(), 4);
                        Ok(())
                    })
                    .unwrap();
                });
            });
        }
    }

    #[test]
    fn inherited_workers_keep_owner_none_and_closing_query_with_global_attachment() {
        // Changing the actual compatibility fallback must not affect unrelated
        // unit tests. Run this one test in its own existing Rust test binary.
        const CHILD: &str = "SC_FABER_TEST_INHERITED_REGISTRY_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "cache_source::tests::inherited_workers_keep_owner_none_and_closing_query_with_global_attachment",
                    "--test-threads=1",
                ])
                .env(CHILD, "1")
                .env_remove("SC_FABER_CACHE_MEASURE_PULL")
                .env_remove("SC_FABER_CACHE_MEASURE_SOURCE_ONLY")
                .status().unwrap();
            assert!(
                status.success(),
                "isolated registry inheritance test failed: {status}"
            );
            return;
        }

        crate::socle::attach_cache(":memory:").unwrap();
        let global = current_registry().unwrap();
        let owner = Arc::new(Registry::new());
        let storage = storage::Storage::open(":memory:", owner.clone()).unwrap();
        owner.bind_storage(storage.handle());
        let expected = Q::from_signeds(1, 1152);
        with_registry(owner.clone(), || {
            with_inherited_registry(|context| {
                let barrier = std::sync::Barrier::new(4);
                let results = sc_num::parallel::chunks(
                    4,
                    1,
                    4,
                    || (),
                    |_, _| {
                        context.run(|| {
                            barrier.wait();
                            assert_eq!(current_registry().unwrap().id(), owner.id());
                            assert!(has_active_scope(owner.id()));
                            assert_eq!(owner.activity.load(Ordering::Acquire) & ACTIVE, 1);
                            (
                                std::thread::current().id(),
                                crate::socle::socle_formula(2, &[4], &[]).unwrap(),
                            )
                        })
                    },
                );
                let mut threads = Vec::new();
                for (thread, value) in results {
                    assert!(
                        !threads.contains(&thread),
                        "test must execute four real workers"
                    );
                    threads.push(thread);
                    assert_eq!(value, expected);
                }
                Ok(())
            })
            .unwrap();
        });
        let root = IntegralKey::Integral {
            g: 2,
            n: 1,
            kappa: vec![],
            psi: vec![(1, 4)],
            lambda: vec![],
            ch: vec![],
        };
        assert!(storage.flush().unwrap().rows > 0);
        assert_eq!(owner.lookup(root.clone()).unwrap(), Some(expected.clone()));
        assert!(global.sources().is_empty());
        assert_eq!(global.lookup(root.clone()).unwrap(), None);

        with_registry(owner.clone(), || {
            with_inherited_registry(|context| {
                std::thread::scope(|scope| {
                    scope
                        .spawn(|| {
                            assert_eq!(current_registry().unwrap().id(), global.id());
                            let panic =
                                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                    context.run(|| panic!("test inherited worker unwind"));
                                }));
                            assert!(panic.is_err());
                            assert!(!has_active_scope(owner.id()));
                            assert_eq!(current_registry().unwrap().id(), global.id());
                        })
                        .join()
                        .unwrap();
                });
                Ok(())
            })
            .unwrap();
        });

        let graph = crate::graph::StableGraph::new(&[2], &[vec![1]], &[]).unwrap();
        let class = crate::taut::TautClass {
            terms: (0..64)
                .map(|_| crate::taut::DecStratum {
                    graph: graph.clone(),
                    poly: crate::kppoly::KappaPsiPolynomial::one(),
                })
                .collect(),
        };
        with_registry(owner.clone(), || {
            assert_eq!(
                crate::classes::int_against(&class, &[], &[(1, 4)]).unwrap(),
                Q::from_signeds(1, 18)
            );
        });
        assert!(global.sources().is_empty());

        with_optional_registry(None, || {
            with_inherited_registry(|context| {
                let barrier = std::sync::Barrier::new(4);
                let results = sc_num::parallel::chunks(
                    4,
                    1,
                    4,
                    || (),
                    |_, _| {
                        context.run(|| {
                            barrier.wait();
                            assert!(current_registry().is_none());
                            crate::socle::socle_formula(2, &[4], &[]).unwrap()
                        })
                    },
                );
                assert!(results.into_iter().all(|value| value == expected));
                Ok(())
            })
            .unwrap();
        });
        assert!(global.sources().is_empty());

        // Close starts after parent admission but before any worker binding.
        // The closer cannot pass wait_idle until the synchronous pool and its
        // parent's original permit have both finished.
        let (started, closing) = std::sync::mpsc::channel();
        let (done, closed) = std::sync::mpsc::channel();
        let closer = with_registry(owner.clone(), || {
            with_inherited_registry(|context| {
                let closing_owner = owner.clone();
                let closer = std::thread::spawn(move || {
                    closing_owner.begin_close();
                    started.send(()).unwrap();
                    closing_owner.wait_idle();
                    done.send(()).unwrap();
                });
                closing.recv().unwrap();
                assert!(owner.is_closing());
                let barrier = std::sync::Barrier::new(4);
                let values = sc_num::parallel::chunks(
                    4,
                    1,
                    4,
                    || (),
                    |_, _| {
                        context.run(|| {
                            barrier.wait();
                            assert_eq!(current_registry().unwrap().id(), owner.id());
                            assert!(has_active_scope(owner.id()));
                            crate::socle::socle_formula(2, &[4], &[])
                        })
                    },
                );
                for value in values {
                    assert_eq!(value.unwrap(), expected);
                }
                assert!(closed.try_recv().is_err());
                Ok(closer)
            })
            .unwrap()
        });
        closed
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap();
        closer.join().unwrap();
        assert_eq!(owner.activity.load(Ordering::Acquire) & ACTIVE, 0);
        storage.flush().unwrap();
        assert_eq!(owner.lookup(root.clone()).unwrap(), Some(expected));
        assert!(global.sources().is_empty());
        assert_eq!(global.lookup(root).unwrap(), None);
        storage.close().unwrap();
        crate::socle::close_cache().unwrap();
    }

    #[test]
    fn captured_query_owner_stays_fixed_when_scope_changes() {
        let first = Arc::new(Registry::new());
        let second = Arc::new(Registry::new());
        with_registry(Arc::clone(&first), || {
            let query_owner = current_registry().unwrap();
            with_registry(Arc::clone(&second), || {
                let dropped = Arc::new(AtomicUsize::new(0));
                query_owner.register(Arc::new(Source { dropped }));
                assert_eq!(current_registry().unwrap().id(), second.id());
            });
            assert_eq!(first.sources().len(), 1);
            assert!(second.sources().is_empty());
        });
    }

    #[test]
    fn soft_value_size_includes_inline_storage_and_large_limbs() {
        assert_eq!(value_bytes(&Q::from(1)), std::mem::size_of::<Q>());
        let mut large = Q::from(1);
        for _ in 0..200 {
            large *= Q::from(2);
        }
        assert!(value_bytes(&large) > value_bytes(&Q::from(1)));
    }

    #[test]
    fn reclaim_gate_waits_for_queries_and_releases_waiting_entries() {
        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            assert!(registry.try_reclaim().is_none())
        });
        let reclaim = registry.try_reclaim().unwrap();
        let (started, start) = std::sync::mpsc::channel();
        let (entered, enter) = std::sync::mpsc::channel();
        let worker_registry = registry.clone();
        let worker = std::thread::spawn(move || {
            started.send(()).unwrap();
            try_with_registry(worker_registry, || entered.send(()).unwrap()).unwrap();
        });
        start.recv().unwrap();
        assert!(
            enter
                .recv_timeout(std::time::Duration::from_millis(10))
                .is_err()
        );
        drop(reclaim);
        enter
            .recv_timeout(std::time::Duration::from_secs(3))
            .unwrap();
        worker.join().unwrap();
        assert!(registry.try_reclaim().is_some());
    }

    #[test]
    fn resource_errors_return_without_swallowing_mathematical_panics() {
        let registry = Arc::new(Registry::new());
        let result = catch_resource(|| -> Result<Q, String> {
            with_registry(registry.clone(), || resource_failure("test memory floor"))
        });
        assert_eq!(result.unwrap_err(), "test memory floor");
        assert!(registry.try_reclaim().is_some());
        assert!(
            std::panic::catch_unwind(|| {
                let _ = catch_resource::<Q>(|| panic!("ordinary mathematical bug"));
            })
            .is_err()
        );
    }
}
