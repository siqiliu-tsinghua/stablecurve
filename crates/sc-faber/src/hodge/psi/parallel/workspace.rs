//! Whole-block pressure backing, separate from the ordinary Plan.terms path.
//!
//! No term is split, recollected, or evaluated here. Jobs obtain a bounded set
//! of Arc<TermBlock> before entering arithmetic. The storage worker alone pulls,
//! encodes, and restores blocks; no mathematical value is copied for handoff.

use super::{Coefficient, CompactTerm, Term, TermBlock, WIDE};
use crate::cache_source::budget::Reservation;
use crate::cache_source::value_bytes;
use crate::cache_source::work::{WorkBlock, WorkLoad, WorkPayload, WorkSource};
use sc_num::{Q, decode_rational, encode_rational_into};
use std::collections::VecDeque;
use std::io::{ErrorKind, Read, Write};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, OnceLock, TryLockError, Weak};
use std::thread::Thread;
use std::time::{Duration, Instant};

const VERSION: u8 = 1;
const COMPACT_BYTES: usize = 16;
const WINDOW_BYTES: usize = 64 * 1024;
const READ_GROUP_TIMEOUT: Duration = Duration::from_secs(30);
const EVICTION_VISITS: usize = 64;

static NEXT_WORKSPACE_ID: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy)]
struct Shape {
    compact: usize,
    wide: usize,
}

struct Slot {
    resident: Option<Arc<TermBlock>>,
    weak: Weak<TermBlock>,
    resident_bytes: usize,
    shape: Shape,
    stored: bool,
    claimed: bool,
    requested: bool,
    loading: bool,
    load_queued: bool,
    eviction_queued: bool,
    waiters: Arc<AtomicUsize>,
    request_credit: Option<Reservation>,
}

impl Slot {
    fn new(block: TermBlock) -> Self {
        let shape = Shape {
            compact: block.compact.len(),
            wide: block.wide.len(),
        };
        let resident_bytes = block.resident_bytes();
        let block = Arc::new(block);
        Self {
            weak: Arc::downgrade(&block),
            resident: Some(block),
            resident_bytes,
            shape,
            stored: false,
            claimed: false,
            requested: false,
            loading: false,
            load_queued: false,
            eviction_queued: false,
            waiters: Arc::new(AtomicUsize::new(0)),
            request_credit: None,
        }
    }
}

#[derive(Default)]
struct State {
    slots: Vec<Slot>,
    save_cursor: usize,
    loads: VecDeque<usize>,
    eviction: VecDeque<usize>,
    reclaim_requested: usize,
    error: Option<String>,
    retired: bool,
    progress_epoch: u64,
    reclaimed_bytes: u128,
    last_reclaimed_at: Option<Instant>,
    releasing_bytes: usize,
}

pub(super) struct PlanWorkspace {
    id: u64,
    state: Mutex<State>,
    changed: Condvar,
    worker: OnceLock<Thread>,
}

/// The bounded group's existing Arc pins and genuinely missing block IDs.
/// These reuse the original block ownership; no mathematical payload is copied.
pub(super) struct BlockGroup {
    blocks: Vec<(usize, Arc<TermBlock>)>,
    missing: Vec<usize>,
    missing_bytes: u64,
    largest_missing_bytes: u64,
}

impl BlockGroup {
    pub(super) fn decode_bytes(&self, window: u64) -> Result<u64, String> {
        if self.missing.is_empty() {
            return Ok(0);
        }
        self.missing_bytes
            .checked_add(self.largest_missing_bytes)
            .and_then(|bytes| bytes.checked_add(window))
            .ok_or_else(|| "pressure block admission overflow".into())
    }
}

/// Constant-time state for ordinary pressure retries. Full slot diagnostics
/// belong to `WorkspaceSnapshot`, outside the repeated waiting path.
#[derive(Clone, Copy, Debug)]
pub(super) struct WorkspaceProgress {
    pub(super) progress_epoch: u64,
    pub(super) last_reclaimed_at: Option<Instant>,
}

#[derive(Default)]
pub(super) struct WorkspaceSnapshot {
    pub(super) progress_epoch: u64,
    pub(super) blocks: usize,
    pub(super) stored_blocks: usize,
    pub(super) unsaved_blocks: usize,
    pub(super) claimed_blocks: usize,
    pub(super) resident_blocks: usize,
    pub(super) unsaved_resident_bytes: usize,
    /// Live block payload estimates, including pins held outside this source.
    pub(super) resident_bytes: usize,
    pub(super) source_owned_bytes: usize,
    /// Subset of resident_bytes, not additional memory.
    pub(super) pinned_bytes: usize,
    pub(super) pinned_blocks: usize,
    pub(super) waiting_clients: usize,
    pub(super) requested_blocks: usize,
    pub(super) loading_blocks: usize,
    /// Payload estimates counted only after the last owned Arc was dropped.
    /// Allocator retention means this is not a claim that RSS fell by as much.
    pub(super) reclaimed_bytes: u128,
    pub(super) last_reclaimed_at: Option<Instant>,
    pub(super) releasing_bytes: usize,
    pub(super) retired: bool,
    pub(super) failed: bool,
}

impl std::fmt::Debug for WorkspaceSnapshot {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("WorkspaceSnapshot")
            .field("progress_epoch", &self.progress_epoch)
            .field("blocks", &self.blocks)
            .field("stored_blocks", &self.stored_blocks)
            .field("unsaved_blocks", &self.unsaved_blocks)
            .field("claimed_blocks", &self.claimed_blocks)
            .field("resident_blocks", &self.resident_blocks)
            .field("unsaved_resident_bytes", &self.unsaved_resident_bytes)
            .field("resident_bytes", &self.resident_bytes)
            .field("source_owned_bytes", &self.source_owned_bytes)
            .field("pinned_bytes", &self.pinned_bytes)
            .field("pinned_blocks", &self.pinned_blocks)
            .field("waiting_clients", &self.waiting_clients)
            .field("requested_blocks", &self.requested_blocks)
            .field("loading_blocks", &self.loading_blocks)
            .field("reclaimed_bytes", &self.reclaimed_bytes)
            .field("last_reclaimed_at", &self.last_reclaimed_at)
            .field("releasing_bytes", &self.releasing_bytes)
            .field("retired", &self.retired)
            .field("failed", &self.failed)
            .finish()
    }
}

impl PlanWorkspace {
    /// Move already sealed blocks only when entering the pressure path.
    pub(super) fn new(blocks: Vec<TermBlock>) -> Arc<Self> {
        let mut id = NEXT_WORKSPACE_ID.load(Ordering::Relaxed);
        loop {
            let next = id.checked_add(1).expect("workspace IDs exhausted");
            match NEXT_WORKSPACE_ID.compare_exchange_weak(
                id,
                next,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(current) => id = current,
            }
        }
        Arc::new(Self {
            id,
            state: Mutex::new(State {
                slots: blocks.into_iter().map(Slot::new).collect(),
                ..State::default()
            }),
            changed: Condvar::new(),
            worker: OnceLock::new(),
        })
    }

    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(|error| error.into_inner())
    }

    fn wake_worker(&self) {
        if let Some(worker) = self.worker.get() {
            worker.unpark();
        }
    }

    fn try_lock(&self) -> Option<MutexGuard<'_, State>> {
        match self.state.try_lock() {
            Ok(state) => Some(state),
            Err(TryLockError::Poisoned(error)) => Some(error.into_inner()),
            Err(TryLockError::WouldBlock) => None,
        }
    }

    /// Admit the complete replacement allocation before growing the slot Vec.
    /// The planner is its single producer: storage changes slot contents but
    /// never appends. Budget callbacks run without the workspace mutex held.
    pub(super) fn reserve_blocks(
        &self,
        additional: usize,
        mut account: impl FnMut(usize),
    ) -> Result<(), String> {
        let (length, capacity, target) = {
            let state = self.lock();
            Self::check_state(&state)?;
            let length = state.slots.len();
            let capacity = state.slots.capacity();
            let needed = length
                .checked_add(additional)
                .ok_or("plan slot count overflow")?;
            if needed <= capacity {
                return Ok(());
            }
            let target = needed.max(capacity.saturating_mul(2)).max(4);
            (length, capacity, target)
        };
        let allocated = target
            .checked_mul(std::mem::size_of::<Slot>())
            .ok_or("plan slot capacity byte count overflow")?;
        // The old Vec is already part of native usage. Its entire replacement,
        // rather than only the added slots, coexists with it during a realloc.
        account(allocated);
        let actual = {
            let mut state = self.lock();
            Self::check_state(&state)?;
            if state.slots.len() != length || state.slots.capacity() != capacity {
                return Err("plan slot producer changed during capacity admission".into());
            }
            state
                .slots
                .try_reserve_exact(target - length)
                .map_err(|error| format!("reserve plan workspace slots: {error}"))?;
            state.slots.capacity()
        };
        if actual > target {
            let extra = (actual - target)
                .checked_mul(std::mem::size_of::<Slot>())
                .ok_or("plan slot allocator excess byte count overflow")?;
            // An allocator returning excess capacity must be accounted for
            // explicitly, before any new block is published into that buffer.
            account(extra);
        }
        Ok(())
    }

    pub(super) fn push_block(&self, block: TermBlock) -> Result<usize, String> {
        // Account/Arc creation happens before acquiring the workspace lock.
        let slot = Slot::new(block);
        let mut state = self.lock();
        Self::check_state(&state)?;
        let index = state.slots.len();
        if index == state.slots.capacity() {
            return Err("plan workspace slots must be reserved before publishing a block".into());
        }
        state.slots.push(slot);
        state.progress_epoch = state.progress_epoch.saturating_add(1);
        drop(state);
        self.changed.notify_all();
        Ok(index)
    }

    fn check_state(state: &State) -> Result<(), String> {
        if let Some(error) = &state.error {
            return Err(error.clone());
        }
        if state.retired {
            return Err("plan workspace has retired".into());
        }
        Ok(())
    }

    /// Pin live payload under the same mutex used by reclamation. Only missing
    /// blocks need fresh decode credit; a later successful read may make that
    /// credit conservative, while these existing pins cannot become missing.
    pub(super) fn pin_blocks(&self, mut indices: Vec<usize>) -> Result<BlockGroup, String> {
        indices.sort_unstable();
        indices.dedup();
        let mut group = BlockGroup {
            blocks: Vec::with_capacity(indices.len()),
            missing: Vec::new(),
            missing_bytes: 0,
            largest_missing_bytes: 0,
        };
        let mut missing = 0;
        {
            let state = self.lock();
            Self::check_state(&state)?;
            if indices
                .last()
                .is_some_and(|&index| index >= state.slots.len())
            {
                return Err("unknown plan block".into());
            }
            for position in 0..indices.len() {
                let index = indices[position];
                let slot = &state.slots[index];
                if let Some(block) = slot
                    .resident
                    .as_ref()
                    .map(Arc::clone)
                    .or_else(|| slot.weak.upgrade())
                {
                    group.blocks.push((index, block));
                } else {
                    if !slot.stored {
                        return Err("plan block lost before successful storage".into());
                    }
                    let bytes = slot.resident_bytes as u64;
                    group.missing_bytes = group
                        .missing_bytes
                        .checked_add(bytes)
                        .ok_or("pressure block admission overflow")?;
                    group.largest_missing_bytes = group.largest_missing_bytes.max(bytes);
                    indices[missing] = index;
                    missing += 1;
                }
            }
        }
        indices.truncate(missing);
        group.missing = indices;
        Ok(group)
    }

    /// Request a bounded drain at a pressure boundary. The consumer alone
    /// evicts blocks; this call changes only small metadata.
    pub(super) fn request_reclaim(&self, bytes: usize) {
        if bytes == 0 {
            return;
        }
        let mut state = self.lock();
        if state.retired || state.error.is_some() {
            return;
        }
        state.reclaim_requested = state.reclaim_requested.max(bytes);
        drop(state);
        self.wake_worker();
    }

    /// Pressure-only block preparation. Caller must hold no result-table,
    /// planner-index, budget, or SQL lock while waiting for backend readback.
    #[cfg(test)]
    pub(super) fn block(&self, index: usize) -> Result<Arc<TermBlock>, String> {
        self.block_inner(index, None, Instant::now() + READ_GROUP_TIMEOUT)
    }

    /// Carry the coordinator's existing whole-group credit to the backend.
    /// This does not start, enlarge, or settle that reservation. The caller's
    /// group deadline remains fixed across all blocks in that group.
    pub(super) fn block_with_reservation(
        &self,
        index: usize,
        reservation: &Reservation,
        deadline: Instant,
    ) -> Result<Arc<TermBlock>, String> {
        self.block_inner(index, Some(reservation), deadline)
    }

    fn block_inner(
        &self,
        index: usize,
        reservation: Option<&Reservation>,
        deadline: Instant,
    ) -> Result<Arc<TermBlock>, String> {
        let waiters = {
            let state = self.lock();
            Self::check_state(&state)?;
            let slot = state.slots.get(index).ok_or("unknown plan block")?;
            let mut count = slot.waiters.load(Ordering::Relaxed);
            loop {
                let next = count.checked_add(1).ok_or("too many plan block waiters")?;
                match slot.waiters.compare_exchange_weak(
                    count,
                    next,
                    Ordering::AcqRel,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(current) => count = current,
                }
            }
            Arc::clone(&slot.waiters)
        };
        struct Waiting<'a> {
            workspace: &'a PlanWorkspace,
            index: usize,
            waiters: Arc<AtomicUsize>,
        }
        impl Drop for Waiting<'_> {
            fn drop(&mut self) {
                if self.waiters.fetch_sub(1, Ordering::AcqRel) != 1 {
                    return;
                }
                let (credit, changed) = {
                    let mut state = self.workspace.lock();
                    let mut changed = false;
                    let credit = state.slots.get_mut(self.index).and_then(|slot| {
                        // Another client can arrive between the decrement and
                        // this lock. Its request and credit must remain live.
                        if slot.waiters.load(Ordering::Acquire) != 0 || slot.loading {
                            return None;
                        }
                        changed = slot.requested;
                        slot.requested = false;
                        // Keep its sole queue entry until next_load discards
                        // it. A new client can reactivate this same entry,
                        // without leaving duplicate IDs behind on cancellation.
                        slot.request_credit.take()
                    });
                    if changed {
                        state.progress_epoch = state.progress_epoch.saturating_add(1);
                    }
                    (credit, changed)
                };
                // A final credit drop may acquire the budget ledger. Never
                // hold the workspace mutex while releasing that ownership.
                drop(credit);
                if changed {
                    self.workspace.changed.notify_all();
                }
            }
        }
        // A pending reader prevents reclaiming a just-restored block before
        // the notified client has cloned its Arc. The state guard below drops
        // before this unwind-safe request cleanup can reacquire the mutex.
        let _waiting = Waiting {
            workspace: self,
            index,
            waiters,
        };
        let mut state = self.lock();
        loop {
            Self::check_state(&state)?;
            let remaining = deadline.saturating_duration_since(Instant::now());
            if remaining.is_zero() {
                return Err("plan workspace resource deadline exceeded while loading a block group (30-second production limit)".into());
            }
            let slot = state.slots.get_mut(index).ok_or("unknown plan block")?;
            if let Some(block) = &slot.resident {
                return Ok(Arc::clone(block));
            }
            if let Some(block) = slot.weak.upgrade() {
                return Ok(block);
            }
            if !slot.stored {
                return Err("plan block lost before successful storage".into());
            }
            let mut changed = !slot.requested;
            slot.requested = true;
            if !slot.loading && slot.request_credit.is_none() {
                if let Some(reservation) = reservation {
                    slot.request_credit = Some(reservation.clone());
                    changed = true;
                }
            }
            let enqueue = !slot.loading && !slot.load_queued;
            if enqueue {
                slot.load_queued = true;
            }
            if enqueue {
                state.loads.push_back(index);
            }
            if changed || enqueue {
                state.progress_epoch = state.progress_epoch.saturating_add(1);
                self.changed.notify_all();
                self.wake_worker();
            }
            let (next, _) = self
                .changed
                .wait_timeout(state, remaining)
                .unwrap_or_else(|error| error.into_inner());
            state = next;
        }
    }

    /// Only the prepared group's missing blocks can request readback. Its
    /// existing pins remain live through loading and the caller's arithmetic.
    pub(super) fn load_pinned_blocks(
        &self,
        group: BlockGroup,
        reservation: &Reservation,
    ) -> Result<Vec<(usize, Arc<TermBlock>)>, String> {
        self.load_pinned_blocks_inner(
            group,
            Some(reservation),
            Instant::now() + READ_GROUP_TIMEOUT,
        )
    }

    #[cfg(test)]
    fn load_blocks_inner(
        &self,
        indices: &[usize],
        reservation: Option<&Reservation>,
        deadline: Instant,
    ) -> Result<Vec<(usize, Arc<TermBlock>)>, String> {
        self.load_pinned_blocks_inner(self.pin_blocks(indices.to_vec())?, reservation, deadline)
    }

    fn load_pinned_blocks_inner(
        &self,
        mut group: BlockGroup,
        reservation: Option<&Reservation>,
        deadline: Instant,
    ) -> Result<Vec<(usize, Arc<TermBlock>)>, String> {
        {
            let state = self.lock();
            Self::check_state(&state)?;
        }
        // The deadline belongs to the complete group. A successful first load
        // never grants another thirty seconds to the next missing block.
        for index in group.missing {
            let block = match reservation {
                Some(reservation) => self.block_with_reservation(index, reservation, deadline),
                None => self.block_inner(index, None, deadline),
            }?;
            group.blocks.push((index, block));
        }
        group.blocks.sort_unstable_by_key(|(index, _)| *index);
        // Already resident pins bypass block_inner, but must still obey the
        // same original whole-group deadline as blocks restored from disk.
        if Instant::now() >= deadline {
            return Err("plan workspace resource deadline exceeded while loading a block group (30-second production limit)".into());
        }
        Ok(group.blocks)
    }

    pub(super) fn snapshot(&self) -> WorkspaceSnapshot {
        let state = self.lock();
        Self::snapshot_of(&state)
    }

    pub(super) fn progress(&self) -> Result<WorkspaceProgress, String> {
        let state = self.lock();
        Self::check_state(&state)?;
        Ok(Self::progress_of(&state))
    }

    fn progress_of(state: &State) -> WorkspaceProgress {
        WorkspaceProgress {
            progress_epoch: state.progress_epoch,
            last_reclaimed_at: state.last_reclaimed_at,
        }
    }

    /// Recheck a prospective stall against the latest published release clock.
    /// This short critical section never traverses the slot directory.
    pub(super) fn reclamation_stalled(
        &self,
        started: Instant,
        timeout: Duration,
    ) -> Result<bool, String> {
        self.reclamation_stalled_with_clock(started, timeout, Instant::now)
    }

    fn reclamation_stalled_with_clock(
        &self,
        started: Instant,
        timeout: Duration,
        now: impl FnOnce() -> Instant,
    ) -> Result<bool, String> {
        let state = self.lock();
        Self::check_state(&state)?;
        let since = state
            .last_reclaimed_at
            .map_or(started, |last| started.max(last));
        Ok(now().saturating_duration_since(since) >= timeout)
    }

    fn snapshot_of(state: &State) -> WorkspaceSnapshot {
        let mut result = WorkspaceSnapshot {
            progress_epoch: state.progress_epoch,
            blocks: state.slots.len(),
            retired: state.retired,
            failed: state.error.is_some(),
            reclaimed_bytes: state.reclaimed_bytes,
            last_reclaimed_at: state.last_reclaimed_at,
            releasing_bytes: state.releasing_bytes,
            ..WorkspaceSnapshot::default()
        };
        for slot in &state.slots {
            result.stored_blocks += usize::from(slot.stored);
            result.unsaved_blocks += usize::from(!slot.stored);
            result.claimed_blocks += usize::from(slot.claimed);
            result.requested_blocks += usize::from(slot.requested);
            result.loading_blocks += usize::from(slot.loading);
            result.waiting_clients = result
                .waiting_clients
                .saturating_add(slot.waiters.load(Ordering::Acquire));
            let owners = slot.weak.strong_count();
            let source_owns = usize::from(slot.resident.is_some());
            if owners != 0 {
                result.resident_blocks += 1;
                result.resident_bytes = result.resident_bytes.saturating_add(slot.resident_bytes);
                if !slot.stored {
                    result.unsaved_resident_bytes = result
                        .unsaved_resident_bytes
                        .saturating_add(slot.resident_bytes);
                }
            }
            if source_owns != 0 {
                result.source_owned_bytes = result
                    .source_owned_bytes
                    .saturating_add(slot.resident_bytes);
            }
            if owners > source_owns {
                result.pinned_blocks += 1;
                result.pinned_bytes = result.pinned_bytes.saturating_add(slot.resident_bytes);
            }
        }
        result
    }

    /// Bounded pressure retry point. No external lock may be held. The epoch
    /// records source transitions, not an RSS guarantee: a backend/job Arc may
    /// still pin a stored block, and its eventual drop needs a native resample.
    pub(super) fn wait_for_progress(
        &self,
        observed_epoch: u64,
        timeout: Duration,
    ) -> Result<WorkspaceProgress, String> {
        let state = self.lock();
        let (state, _) = self
            .changed
            .wait_timeout_while(state, timeout, |state| {
                state.progress_epoch == observed_epoch && state.error.is_none() && !state.retired
            })
            .unwrap_or_else(|error| error.into_inner());
        Self::check_state(&state)?;
        Ok(Self::progress_of(&state))
    }

    /// Call only after the plan's jobs have joined. In-flight backend Arcs remain
    /// safe; callbacks become no-ops and blocked preparations are awakened.
    pub(super) fn retire(&self) {
        let retired = {
            let mut state = self.lock();
            state.retired = true;
            state.progress_epoch = state.progress_epoch.saturating_add(1);
            state.save_cursor = 0;
            state.reclaim_requested = 0;
            (
                std::mem::take(&mut state.slots),
                std::mem::take(&mut state.loads),
                std::mem::take(&mut state.eviction),
            )
        };
        self.changed.notify_all();
        drop(retired);
    }
}

impl WorkSource for PlanWorkspace {
    fn bind_worker(&self, worker: Thread) {
        let _ = self.worker.set(worker);
    }

    fn workspace_id(&self) -> u64 {
        self.id
    }

    fn pull(&self) -> Option<WorkBlock> {
        let mut state = self.try_lock()?;
        Self::check_state(&state).ok()?;
        while state.save_cursor < state.slots.len() {
            let index = state.save_cursor;
            state.save_cursor += 1;
            let slot = &mut state.slots[index];
            if slot.stored || slot.claimed {
                continue;
            }
            let payload = Arc::clone(
                slot.resident
                    .as_ref()
                    .expect("unpersisted plan block must be resident"),
            );
            let payload: Arc<dyn WorkPayload> = payload;
            slot.claimed = true;
            return Some(WorkBlock {
                block: index,
                payload,
            });
        }
        None
    }

    fn next_load(&self) -> Option<WorkLoad> {
        let mut state = self.try_lock()?;
        Self::check_state(&state).ok()?;
        while let Some(index) = state.loads.pop_front() {
            let Some(slot) = state.slots.get_mut(index) else {
                continue;
            };
            slot.load_queued = false;
            if !slot.requested || slot.loading {
                continue;
            }
            slot.loading = true;
            let load = WorkLoad {
                block: index,
                resident_bytes: slot.resident_bytes,
                reservation: slot.request_credit.take(),
            };
            state.progress_epoch = state.progress_epoch.saturating_add(1);
            return Some(load);
        }
        None
    }

    fn stored(&self, index: usize) -> Result<(), String> {
        {
            let mut state = self.lock();
            if state.retired {
                return Ok(());
            }
            Self::check_state(&state)?;
            let slot = state
                .slots
                .get_mut(index)
                .ok_or("unknown stored plan block")?;
            if slot.stored || !slot.claimed {
                return Err("plan block storage acknowledgement is out of order".into());
            }
            slot.stored = true;
            slot.claimed = false;
            // A successful save permits eviction; it does not require it.
            // Retaining this owner avoids rereading a block at the next rank.
            if !slot.eviction_queued {
                slot.eviction_queued = true;
                state.eviction.push_back(index);
            }
            state.progress_epoch = state.progress_epoch.saturating_add(1);
        }
        self.changed.notify_all();
        Ok(())
    }

    fn restore(&self, index: usize, input: &mut dyn Read) -> Result<(), String> {
        let shape = {
            let state = self.lock();
            if state.retired {
                return Ok(());
            }
            Self::check_state(&state)?;
            let slot = state
                .slots
                .get(index)
                .ok_or("unknown restored plan block")?;
            if !slot.stored || !slot.loading {
                return Err("plan block readback was not requested".into());
            }
            slot.shape
        };
        let block = match decode_block(input, shape) {
            Ok(block) => block,
            Err(error) => {
                self.fail(error.clone());
                return Err(error);
            }
        };
        let resident_bytes = block.resident_bytes();
        let block = Arc::new(block);
        let mut state = self.lock();
        if state.retired {
            drop(state);
            drop(block);
            return Ok(());
        }
        Self::check_state(&state)?;
        let slot = state
            .slots
            .get_mut(index)
            .ok_or("unknown restored plan block")?;
        slot.weak = Arc::downgrade(&block);
        slot.resident = Some(block);
        slot.resident_bytes = resident_bytes;
        slot.requested = false;
        slot.loading = false;
        let credit = slot.request_credit.take();
        let enqueue = !slot.eviction_queued;
        if enqueue {
            slot.eviction_queued = true;
            state.eviction.push_back(index);
        }
        state.progress_epoch = state.progress_epoch.saturating_add(1);
        drop(state);
        drop(credit);
        self.changed.notify_all();
        Ok(())
    }

    fn fail(&self, error: String) {
        let mut state = self.lock();
        if state.error.is_none() {
            state.error = Some(error);
            state.progress_epoch = state.progress_epoch.saturating_add(1);
        }
        let credits: Vec<_> = state
            .slots
            .iter_mut()
            .filter_map(|slot| {
                slot.requested = false;
                slot.loading = false;
                slot.load_queued = false;
                slot.eviction_queued = false;
                slot.request_credit.take()
            })
            .collect();
        let queues = (
            std::mem::take(&mut state.loads),
            std::mem::take(&mut state.eviction),
        );
        state.reclaim_requested = 0;
        drop(state);
        self.changed.notify_all();
        drop(credits);
        drop(queues);
    }

    fn take_reclaim_request(&self) -> usize {
        std::mem::take(&mut self.lock().reclaim_requested)
    }

    fn reclaim_pending(&self) -> bool {
        let state = self.lock();
        !state.retired
            && state.error.is_none()
            && (!state.eviction.is_empty() || state.save_cursor < state.slots.len())
    }

    fn reclaim_stored(&self, byte_limit: usize) -> usize {
        if byte_limit == 0 {
            return 0;
        }
        let mut retired = Vec::new();
        let mut bytes = 0usize;
        {
            let mut state = self.lock();
            // Only saved/restored, source-owned blocks enter this rotation.
            // Old saved slots without an owner are never reconsidered here.
            let visits = state.eviction.len().min(EVICTION_VISITS);
            for _ in 0..visits {
                let index = state.eviction.pop_front().unwrap();
                let Some(slot) = state.slots.get_mut(index) else {
                    continue;
                };
                slot.eviction_queued = false;
                if !slot.stored || slot.resident.is_none() {
                    continue;
                }
                if slot.waiters.load(Ordering::Acquire) != 0
                    || Arc::strong_count(slot.resident.as_ref().unwrap()) != 1
                {
                    slot.eviction_queued = true;
                    state.eviction.push_back(index);
                    continue;
                }
                retired.push(slot.resident.take().unwrap());
                slot.weak = Weak::new();
                bytes = bytes.saturating_add(slot.resident_bytes);
                if bytes >= byte_limit {
                    break;
                }
            }
            if bytes != 0 {
                state.progress_epoch = state.progress_epoch.saturating_add(1);
                state.releasing_bytes = state.releasing_bytes.saturating_add(bytes);
            }
        }
        drop(retired);
        if bytes != 0 {
            let reclaimed_at = Instant::now();
            // Every selected Arc had no external pins under the source lock.
            // Only its completed destructor, not storage/claim/epoch activity,
            // is evidence that this payload allocation was actually released.
            let mut state = self.lock();
            state.releasing_bytes = state.releasing_bytes.saturating_sub(bytes);
            state.reclaimed_bytes = state.reclaimed_bytes.saturating_add(bytes as u128);
            state.last_reclaimed_at = Some(
                state
                    .last_reclaimed_at
                    .map_or(reclaimed_at, |last| last.max(reclaimed_at)),
            );
            state.progress_epoch = state.progress_epoch.saturating_add(1);
            drop(state);
            self.changed.notify_all();
        }
        bytes
    }

    fn is_retired(&self) -> bool {
        self.lock().retired
    }
}

fn add_size(left: u64, right: u64) -> Result<u64, String> {
    left.checked_add(right)
        .ok_or_else(|| "plan block encoding exceeds u64".into())
}

fn count(value: usize) -> Result<u64, String> {
    u64::try_from(value).map_err(|_| "plan block count exceeds u64".into())
}

fn rational_size(value: &Q) -> Result<u64, String> {
    let mut size = 2;
    for natural in [value.numerator_ref(), value.denominator_ref()] {
        let length = match natural.limbs().next_back() {
            None => 0,
            Some(high) => {
                let limb_bytes = count(std::mem::size_of_val(&high))?;
                let lower = natural
                    .limb_count()
                    .checked_sub(1)
                    .and_then(|limbs| limbs.checked_mul(limb_bytes))
                    .ok_or("rational limb size overflow")?;
                let bits = limb_bytes
                    .checked_mul(8)
                    .and_then(|bits| bits.checked_sub(u64::from(high.leading_zeros())))
                    .ok_or("rational high limb size overflow")?;
                add_size(lower, bits.div_ceil(8))?
            }
        };
        let mut prefix = 1;
        let mut remaining = length;
        while remaining >= 128 {
            prefix += 1;
            remaining >>= 7;
        }
        size = add_size(add_size(size, prefix)?, length)?;
    }
    Ok(size)
}

impl WorkPayload for TermBlock {
    fn resident_bytes(&self) -> usize {
        self.wide.iter().fold(
            self.compact
                .capacity()
                .saturating_mul(std::mem::size_of::<CompactTerm>())
                .saturating_add(
                    self.wide
                        .capacity()
                        .saturating_mul(std::mem::size_of::<Term>()),
                ),
            |bytes, term| match &term.factor {
                Coefficient::Rational(value) => bytes.saturating_add(value_bytes(value)),
                _ => bytes,
            },
        )
    }

    fn encoded_len(&self) -> Result<u64, String> {
        let mut size = add_size(
            17,
            count(self.compact.len())?
                .checked_mul(16)
                .ok_or("compact plan block size overflow")?,
        )?;
        for term in &self.wide {
            size = add_size(size, 25)?;
            if let Coefficient::Rational(value) = &term.factor {
                size = add_size(size, rational_size(value)?)?;
            }
        }
        Ok(size)
    }

    fn encode(&self, output: &mut dyn Write) -> Result<(), String> {
        let write = |output: &mut dyn Write, bytes: &[u8]| {
            output.write_all(bytes).map_err(|error| error.to_string())
        };
        write(output, &[VERSION])?;
        write(output, &count(self.compact.len())?.to_le_bytes())?;
        write(output, &count(self.wide.len())?.to_le_bytes())?;
        let mut window = [0u8; WINDOW_BYTES];
        for chunk in self.compact.chunks(WINDOW_BYTES / COMPACT_BYTES) {
            for (term, bytes) in chunk.iter().zip(window.chunks_exact_mut(COMPACT_BYTES)) {
                bytes[..8].copy_from_slice(&term.numerator.to_le_bytes());
                bytes[8..12].copy_from_slice(&term.left.to_le_bytes());
                bytes[12..16].copy_from_slice(&term.right.to_le_bytes());
            }
            write(output, &window[..chunk.len() * COMPACT_BYTES])?;
        }
        let mut rational = Vec::new();
        for term in &self.wide {
            match &term.factor {
                Coefficient::Integer(value) => {
                    write(output, &[0])?;
                    write(output, &value.to_le_bytes())?;
                }
                Coefficient::Half(value) => {
                    write(output, &[1])?;
                    write(output, &value.to_le_bytes())?;
                }
                Coefficient::Rational(value) => {
                    rational.clear();
                    encode_rational_into(value, &mut rational);
                    write(output, &[2])?;
                    write(output, &count(rational.len())?.to_le_bytes())?;
                    write(output, &rational)?;
                }
            }
            write(output, &count(term.left)?.to_le_bytes())?;
            let right = if term.right == usize::MAX {
                u64::MAX
            } else {
                count(term.right)?
            };
            write(output, &right.to_le_bytes())?;
        }
        Ok(())
    }
}

fn read_array<const N: usize>(input: &mut dyn Read) -> Result<[u8; N], String> {
    let mut bytes = [0; N];
    input
        .read_exact(&mut bytes)
        .map_err(|error| error.to_string())?;
    Ok(bytes)
}

fn decode_block(input: &mut dyn Read, expected: Shape) -> Result<TermBlock, String> {
    if read_array::<1>(input)?[0] != VERSION {
        return Err("unsupported plan block encoding version".into());
    }
    let compact = u64::from_le_bytes(read_array(input)?);
    let wide = u64::from_le_bytes(read_array(input)?);
    if compact != count(expected.compact)? || wide != count(expected.wide)? {
        return Err("plan block shape differs from the original sealed block".into());
    }
    let mut block = TermBlock::default();
    block
        .compact
        .try_reserve_exact(expected.compact)
        .map_err(|error| error.to_string())?;
    block
        .wide
        .try_reserve_exact(expected.wide)
        .map_err(|error| error.to_string())?;
    let mut window = [0u8; WINDOW_BYTES];
    let mut remaining = expected.compact;
    while remaining != 0 {
        let terms = remaining.min(WINDOW_BYTES / COMPACT_BYTES);
        let bytes = &mut window[..terms * COMPACT_BYTES];
        input.read_exact(bytes).map_err(|error| error.to_string())?;
        for bytes in bytes.chunks_exact(COMPACT_BYTES) {
            let numerator = i64::from_le_bytes(bytes[..8].try_into().unwrap());
            let left = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
            let right = u32::from_le_bytes(bytes[12..16].try_into().unwrap());
            if right == WIDE
                && usize::try_from(numerator)
                    .ok()
                    .is_none_or(|index| index >= expected.wide)
            {
                return Err("invalid wide coefficient reference in plan block".into());
            }
            block.compact.push(CompactTerm {
                numerator,
                left,
                right,
            });
        }
        remaining -= terms;
    }
    let mut rational = Vec::new();
    for _ in 0..expected.wide {
        let factor = match read_array::<1>(input)?[0] {
            0 => Coefficient::Integer(i64::from_le_bytes(read_array(input)?)),
            1 => Coefficient::Half(i64::from_le_bytes(read_array(input)?)),
            2 => {
                let length = u64::from_le_bytes(read_array(input)?);
                rational.clear();
                (&mut *input)
                    .take(length)
                    .read_to_end(&mut rational)
                    .map_err(|error| error.to_string())?;
                if count(rational.len())? != length {
                    return Err("truncated rational coefficient in plan block".into());
                }
                Coefficient::Rational(Box::new(decode_rational(&rational)?))
            }
            _ => return Err("invalid plan block coefficient tag".into()),
        };
        let left = usize::try_from(u64::from_le_bytes(read_array(input)?))
            .map_err(|_| "left plan index exceeds this platform")?;
        let right = u64::from_le_bytes(read_array(input)?);
        let right = if right == u64::MAX {
            usize::MAX
        } else {
            let index =
                usize::try_from(right).map_err(|_| "right plan index exceeds this platform")?;
            if index == usize::MAX {
                return Err(
                    "right plan index conflicts with this platform's linear sentinel".into(),
                );
            }
            index
        };
        block.wide.push(Term {
            factor,
            left,
            right,
        });
    }
    let mut extra = [0u8; 1];
    loop {
        match input.read(&mut extra) {
            Ok(0) => return Ok(block),
            Ok(_) => return Err("trailing data after plan block".into()),
            Err(error) if error.kind() == ErrorKind::Interrupted => continue,
            Err(error) => return Err(error.to_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::{One, Zero, pow_unsigned};
    use std::io::Cursor;

    fn sample_block() -> TermBlock {
        let mut block = TermBlock::default();
        for (factor, left, right) in [
            (Coefficient::Integer(-3), 0, usize::MAX),
            (Coefficient::Half(7), 1, 2),
            (
                Coefficient::Rational(Box::new(Q::from_signeds(-5, 7))),
                2,
                usize::MAX,
            ),
        ] {
            block.push(Term {
                factor,
                left,
                right,
            });
        }
        block
    }

    fn bytes(block: &TermBlock) -> Vec<u8> {
        let mut bytes = Vec::new();
        block.encode(&mut bytes).unwrap();
        assert_eq!(block.encoded_len().unwrap(), bytes.len() as u64);
        bytes
    }

    fn next_readback(workspace: &PlanWorkspace) -> WorkLoad {
        let deadline = Instant::now() + Duration::from_secs(5);
        loop {
            if let Some(load) = workspace.next_load() {
                return load;
            }
            assert!(Instant::now() < deadline, "missing queued plan readback");
            std::thread::yield_now();
        }
    }

    fn evict_unpinned_fixture(workspace: &PlanWorkspace) {
        loop {
            let before = workspace.lock().eviction.len();
            if before == 0 {
                return;
            }
            workspace.reclaim_stored(usize::MAX);
            assert!(
                workspace.lock().eviction.len() < before,
                "fixture unexpectedly pinned"
            );
        }
    }

    #[test]
    fn resident_group_needs_no_decode_credit_and_reuses_original_arcs() {
        let workspace = PlanWorkspace::new(vec![sample_block(), sample_block()]);
        let originals = [workspace.block(0).unwrap(), workspace.block(1).unwrap()];
        let group = workspace.pin_blocks(vec![1, 0, 1]).unwrap();
        assert!(group.missing.is_empty());
        assert_eq!(group.decode_bytes(1024 * 1024).unwrap(), 0);
        let blocks = workspace
            .load_pinned_blocks_inner(group, None, Instant::now() + Duration::from_secs(1))
            .unwrap();
        assert_eq!(blocks.len(), 2);
        for (index, block) in blocks {
            assert!(Arc::ptr_eq(&block, &originals[index]));
        }
        assert!(workspace.next_load().is_none());
    }

    #[test]
    fn weak_only_live_payload_is_pinned_without_decode_credit() {
        let workspace = PlanWorkspace::new(vec![sample_block()]);
        let original = workspace.block(0).unwrap();
        let owner = workspace.lock().slots[0].resident.take();
        drop(owner);
        let group = workspace.pin_blocks(vec![0]).unwrap();
        assert_eq!(group.decode_bytes(4096).unwrap(), 0);
        assert!(Arc::ptr_eq(&group.blocks[0].1, &original));
        let weak = Arc::downgrade(&original);
        drop(original);
        assert!(weak.upgrade().is_some());
        drop(group);
        assert!(weak.upgrade().is_none());
    }

    #[test]
    fn mixed_group_charges_only_missing_payload_and_handles_a_later_fill() {
        let block = sample_block();
        let encoded = bytes(&block);
        let missing_bytes = block.resident_bytes() as u64;
        let workspace = PlanWorkspace::new(vec![block, sample_block()]);
        let work = workspace.pull().unwrap();
        workspace.stored(0).unwrap();
        drop(work);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let group = workspace.pin_blocks(vec![0, 1]).unwrap();
        assert_eq!(group.missing, vec![0]);
        assert_eq!(group.blocks[0].0, 1);
        assert_eq!(group.decode_bytes(4096).unwrap(), 2 * missing_bytes + 4096);
        // This second preparation sees the same miss. A different admitted
        // group may fill it before the second group actually starts readback.
        let later = workspace.pin_blocks(vec![0]).unwrap();
        assert!(later.decode_bytes(4096).unwrap() > 0);
        let blocks = std::thread::scope(|scope| {
            let backend = scope.spawn(|| {
                let load = next_readback(&workspace);
                assert_eq!(load.block, 0);
                workspace.restore(0, &mut Cursor::new(&encoded)).unwrap();
            });
            let blocks = workspace
                .load_pinned_blocks_inner(group, None, Instant::now() + Duration::from_secs(5))
                .unwrap();
            backend.join().unwrap();
            blocks
        });
        assert_eq!(
            blocks.iter().map(|(index, _)| *index).collect::<Vec<_>>(),
            [0, 1]
        );
        for (_, block) in &blocks {
            assert_eq!(bytes(block), encoded);
        }
        let later = workspace
            .load_pinned_blocks_inner(later, None, Instant::now() + Duration::from_secs(1))
            .unwrap();
        assert!(Arc::ptr_eq(&later[0].1, &blocks[0].1));
        assert!(workspace.next_load().is_none());
    }

    #[test]
    fn group_pins_prevent_concurrent_reclaim_until_the_group_is_released() {
        let workspace = PlanWorkspace::new(vec![sample_block(), sample_block()]);
        for index in 0..2 {
            let work = workspace.pull().unwrap();
            workspace.stored(index).unwrap();
            drop(work);
        }
        let group = workspace.pin_blocks(vec![0, 1]).unwrap();
        let freed = std::thread::scope(|scope| {
            scope
                .spawn(|| workspace.reclaim_stored(usize::MAX))
                .join()
                .unwrap()
        });
        assert_eq!(freed, 0);
        assert_eq!(group.decode_bytes(4096).unwrap(), 0);
        drop(group);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        assert_eq!(workspace.snapshot().resident_bytes, 0);
    }

    #[test]
    fn failed_group_load_releases_existing_pins_before_reclaim() {
        let workspace = PlanWorkspace::new(vec![sample_block(), sample_block()]);
        let original = workspace.block(0).unwrap();
        for index in 0..2 {
            let work = workspace.pull().unwrap();
            workspace.stored(index).unwrap();
            drop(work);
        }
        // The external pin keeps the first block resident while the second
        // becomes missing. No backend services the ensuing group load.
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let group = workspace.pin_blocks(vec![0, 1]).unwrap();
        assert_eq!(group.missing, vec![1]);
        drop(original);
        let error = workspace
            .load_pinned_blocks_inner(group, None, Instant::now() + Duration::from_millis(10))
            .err()
            .unwrap();
        assert!(error.contains("resource deadline exceeded"));
        {
            let state = workspace.lock();
            assert_eq!(
                Arc::strong_count(state.slots[0].resident.as_ref().unwrap()),
                1
            );
            assert_eq!(state.slots[1].waiters.load(Ordering::Acquire), 0);
            assert!(!state.slots[1].requested);
        }
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        assert_eq!(workspace.snapshot().resident_bytes, 0);
    }

    #[test]
    fn slot_growth_from_8192_requires_the_complete_new_buffer_before_allocation() {
        const ORIGINAL: usize = 8192;
        let workspace = PlanWorkspace::new((0..ORIGINAL).map(|_| TermBlock::default()).collect());
        let (old_capacity, old_pointer) = {
            let state = workspace.lock();
            assert_eq!(state.slots.len(), ORIGINAL);
            assert_eq!(state.slots.capacity(), ORIGINAL);
            (state.slots.capacity(), state.slots.as_ptr())
        };
        let required = 2 * ORIGINAL * std::mem::size_of::<Slot>();
        let denied = crate::cache_source::catch_resource(|| {
            workspace.reserve_blocks(1, |bytes| {
                assert!(
                    workspace.state.try_lock().is_ok(),
                    "budget callback held source lock"
                );
                assert_eq!(bytes, required);
                let allowed = required - 1;
                if bytes > allowed {
                    crate::cache_source::resource_failure("test slot growth admission denied");
                }
            })
        });
        assert_eq!(denied.unwrap_err(), "test slot growth admission denied");
        {
            let state = workspace.lock();
            assert_eq!(state.slots.len(), ORIGINAL);
            assert_eq!(state.slots.capacity(), old_capacity);
            assert_eq!(state.slots.as_ptr(), old_pointer);
        }
        let mut admitted = 0usize;
        workspace
            .reserve_blocks(1, |bytes| {
                assert!(workspace.state.try_lock().is_ok());
                admitted = admitted.checked_add(bytes).unwrap();
            })
            .unwrap();
        let capacity = workspace.lock().slots.capacity();
        assert!(capacity >= 2 * ORIGINAL);
        assert_eq!(admitted, capacity * std::mem::size_of::<Slot>());
        assert_eq!(
            workspace.push_block(TermBlock::default()).unwrap(),
            ORIGINAL
        );
        assert_eq!(
            workspace.lock().slots.capacity(),
            capacity,
            "push must not silently grow"
        );
        workspace
            .reserve_blocks(0, |_| panic!("existing capacity needs no new credit"))
            .unwrap();
    }

    #[test]
    fn publishing_a_block_without_preflight_capacity_is_rejected() {
        let workspace = PlanWorkspace::new(Vec::new());
        let error = workspace.push_block(TermBlock::default()).unwrap_err();
        assert!(error.contains("must be reserved"));
        assert_eq!(workspace.lock().slots.len(), 0);
        assert_eq!(workspace.lock().slots.capacity(), 0);
    }

    fn decode(bytes: &[u8], original: &TermBlock) -> Result<TermBlock, String> {
        decode_block(
            &mut Cursor::new(bytes),
            Shape {
                compact: original.compact.len(),
                wide: original.wide.len(),
            },
        )
    }

    #[test]
    fn pressure_read_and_reclaim_requests_wake_the_bound_backend() {
        for read in [false, true] {
            let workspace = PlanWorkspace::new(vec![sample_block()]);
            if read {
                let payload = workspace.pull().unwrap();
                workspace.stored(payload.block).unwrap();
                drop(payload);
                evict_unpinned_fixture(&workspace);
            }
            let (bound, binding) = std::sync::mpsc::channel();
            let (woke, awake) = std::sync::mpsc::channel();
            let backend_workspace = Arc::clone(&workspace);
            let backend = std::thread::spawn(move || {
                backend_workspace.bind_worker(std::thread::current());
                bound.send(()).unwrap();
                std::thread::park_timeout(Duration::from_secs(5));
                if read {
                    assert_eq!(backend_workspace.lock().loads.front(), Some(&0));
                    backend_workspace.fail("test backend observed queued read".into());
                } else {
                    assert_eq!(backend_workspace.take_reclaim_request(), 123);
                }
                woke.send(()).unwrap();
            });
            binding.recv_timeout(Duration::from_secs(5)).unwrap();
            let reader = if read {
                let workspace = Arc::clone(&workspace);
                Some(std::thread::spawn(move || {
                    workspace.block_inner(0, None, Instant::now() + Duration::from_secs(3))
                }))
            } else {
                workspace.request_reclaim(123);
                None
            };
            let notified = awake.recv_timeout(Duration::from_secs(1));
            backend.thread().unpark();
            backend.join().unwrap();
            if let Some(reader) = reader {
                assert!(
                    reader
                        .join()
                        .unwrap()
                        .err()
                        .unwrap()
                        .contains("test backend observed queued read")
                );
            }
            assert!(
                notified.is_ok(),
                "pressure request did not wake the parked backend"
            );
        }
    }

    #[test]
    fn whole_block_codec_preserves_exact_arithmetic_and_compact_layout() {
        let block = sample_block();
        let encoded = bytes(&block);
        assert_eq!(
            &encoded[17..33],
            &[
                (-6i64).to_le_bytes().as_slice(),
                0u32.to_le_bytes().as_slice(),
                u32::MAX.to_le_bytes().as_slice()
            ]
            .concat()
        );
        let restored = decode(&encoded, &block).unwrap();
        let values = [Q::from(2), Q::from(3), Q::from(5)];
        assert_eq!(
            block.evaluate(0..block.len(), &Q::ONE, Some(&Q::from(3)), &values),
            restored.evaluate(0..restored.len(), &Q::ONE, Some(&Q::from(3)), &values)
        );
        assert_eq!(bytes(&restored), encoded);
    }

    #[test]
    fn wide_codec_preserves_large_signed_rationals_and_linear_indices() {
        let large = (pow_unsigned(&Q::from(2), 521) + Q::ONE) / pow_unsigned(&Q::from(3), 181);
        let mut block = TermBlock::default();
        for factor in [
            Coefficient::Integer(i64::MAX),
            Coefficient::Half(i64::MIN),
            Coefficient::Rational(Box::new(-large)),
        ] {
            let id = block.wide.len();
            block.compact.push(CompactTerm {
                numerator: id as i64,
                left: 0,
                right: WIDE,
            });
            block.wide.push(Term {
                factor,
                left: usize::MAX - 2,
                right: usize::MAX,
            });
        }
        let encoded = bytes(&block);
        assert_eq!(bytes(&decode(&encoded, &block).unwrap()), encoded);
    }

    #[test]
    fn malformed_or_incomplete_blocks_are_not_published() {
        let block = sample_block();
        let encoded = bytes(&block);
        for length in 0..encoded.len() {
            assert!(decode(&encoded[..length], &block).is_err());
        }
        let mut bad = encoded.clone();
        bad[0] = 2;
        assert!(decode(&bad, &block).is_err());
        let mut bad = encoded.clone();
        bad.push(0);
        assert!(decode(&bad, &block).is_err());
        let mut bad = encoded.clone();
        bad[1..9].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(decode(&bad, &block).is_err());
        let mut bad = encoded;
        bad[17 + 2 * 16..17 + 2 * 16 + 8].copy_from_slice(&(-1i64).to_le_bytes());
        assert!(decode(&bad, &block).is_err());
    }

    #[test]
    fn stored_ack_does_not_release_backend_or_job_pins() {
        let workspace = PlanWorkspace::new(vec![sample_block()]);
        let job = workspace.block(0).unwrap();
        let work = workspace.pull().unwrap();
        workspace.stored(0).unwrap();
        assert_eq!(workspace.snapshot().pinned_blocks, 1);
        drop(work);
        assert_eq!(workspace.snapshot().pinned_blocks, 1);
        assert!(Arc::ptr_eq(&job, &workspace.block(0).unwrap()));
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
        drop(job);
        assert!(workspace.snapshot().source_owned_bytes > 0);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        assert_eq!(workspace.snapshot().resident_bytes, 0);
    }

    #[test]
    fn saved_blocks_remain_resident_until_storage_handles_a_reclaim_request() {
        let workspace = PlanWorkspace::new(vec![sample_block()]);
        let original = workspace.block(0).unwrap();
        let work = workspace.pull().unwrap();
        workspace.stored(0).unwrap();
        drop(work);
        for _ in 0..3 {
            let reused = workspace.block(0).unwrap();
            assert!(Arc::ptr_eq(&original, &reused));
            assert!(workspace.next_load().is_none());
        }
        drop(original);
        assert!(workspace.snapshot().source_owned_bytes > 0);
        workspace.request_reclaim(usize::MAX);
        assert_eq!(workspace.take_reclaim_request(), usize::MAX);
        assert_eq!(workspace.take_reclaim_request(), 0);
        assert!(workspace.reclaim_pending());
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        assert!(!workspace.reclaim_pending());
        assert_eq!(workspace.snapshot().resident_bytes, 0);
    }

    #[test]
    fn reclamation_clock_moves_only_after_the_last_payload_arc_is_dropped() {
        let workspace = PlanWorkspace::new(vec![sample_block(), sample_block()]);
        let job = workspace.block(0).unwrap();
        let weak = Arc::downgrade(&job);
        let work = workspace.pull().unwrap();
        workspace.stored(0).unwrap();
        let stored = workspace.snapshot();
        assert!(stored.progress_epoch > 0);
        assert_eq!(stored.reclaimed_bytes, 0);
        assert!(stored.last_reclaimed_at.is_none());
        assert_eq!(stored.unsaved_blocks, 1);
        assert!(stored.unsaved_resident_bytes > 0);
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
        drop(work);
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
        assert!(workspace.snapshot().last_reclaimed_at.is_none());
        drop(job);
        let freed = workspace.reclaim_stored(usize::MAX);
        assert!(freed > 0);
        assert!(
            weak.upgrade().is_none(),
            "recording progress must follow the final Arc drop"
        );
        let first = workspace.snapshot();
        assert_eq!(first.reclaimed_bytes, freed as u128);
        assert!(first.last_reclaimed_at.is_some());
        assert_eq!(first.releasing_bytes, 0);

        workspace.request_reclaim(usize::MAX);
        assert_eq!(workspace.take_reclaim_request(), usize::MAX);
        let work = workspace.pull().unwrap();
        workspace.stored(1).unwrap();
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
        let pinned = workspace.snapshot();
        assert_eq!(pinned.last_reclaimed_at, first.last_reclaimed_at);
        assert_eq!(pinned.reclaimed_bytes, first.reclaimed_bytes);
        drop(work);
        let freed_again = workspace.reclaim_stored(usize::MAX);
        assert!(freed_again > 0);
        let finished = workspace.snapshot();
        assert_eq!(
            finished.reclaimed_bytes,
            first.reclaimed_bytes + freed_again as u128
        );
        assert!(finished.last_reclaimed_at.unwrap() >= first.last_reclaimed_at.unwrap());
        assert_eq!(finished.resident_bytes, 0);
        assert_eq!(finished.unsaved_resident_bytes, 0);
        assert_eq!(finished.releasing_bytes, 0);
    }

    #[test]
    fn stale_snapshot_cannot_override_a_release_recorded_during_native_sampling() {
        let workspace = PlanWorkspace::new(Vec::new());
        let started = Instant::now();
        let timeout = Duration::from_secs(30);
        let stale = workspace.snapshot();
        let now = started + Duration::from_secs(31);
        assert!(stale.last_reclaimed_at.is_none());
        {
            // Simulate the release-clock publication after the old snapshot,
            // while the caller was collecting native memory counters.
            let mut state = workspace.lock();
            state.last_reclaimed_at = Some(now - Duration::from_secs(1));
            state.reclaimed_bytes = 17;
        }
        assert!(
            !workspace
                .reclamation_stalled_with_clock(started, timeout, || now)
                .unwrap()
        );
        {
            let mut state = workspace.lock();
            state.last_reclaimed_at = Some(started);
            state.progress_epoch += 100;
        }
        assert!(
            workspace
                .reclamation_stalled_with_clock(started, timeout, || now)
                .unwrap()
        );
        workspace.fail("storage failed before final stall check".into());
        assert_eq!(
            workspace
                .reclamation_stalled_with_clock(started, timeout, || now)
                .unwrap_err(),
            "storage failed before final stall check"
        );
    }

    #[test]
    fn progress_wait_is_bounded_and_reports_failure_instead_of_retrying_forever() {
        let workspace = PlanWorkspace::new(vec![sample_block()]);
        let first = workspace.progress().unwrap();
        let unchanged = workspace
            .wait_for_progress(first.progress_epoch, Duration::ZERO)
            .unwrap();
        assert_eq!(unchanged.progress_epoch, first.progress_epoch);
        assert!(unchanged.last_reclaimed_at.is_none());
        let work = workspace.pull().unwrap();
        workspace.stored(0).unwrap();
        let stored = workspace
            .wait_for_progress(first.progress_epoch, Duration::ZERO)
            .unwrap();
        assert!(stored.progress_epoch > first.progress_epoch);
        // The spill epoch alone is not a declaration that the backend pin died.
        assert!(stored.last_reclaimed_at.is_none());
        let detail = workspace.snapshot();
        assert_eq!(detail.stored_blocks, 1);
        assert_eq!(detail.pinned_blocks, 1);
        drop(work);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let reclaimed = workspace.progress().unwrap();
        assert!(reclaimed.last_reclaimed_at.is_some());
        workspace.fail("workspace stopped".into());
        assert_eq!(workspace.progress().unwrap_err(), "workspace stopped");
        assert!(
            workspace
                .wait_for_progress(stored.progress_epoch, Duration::ZERO)
                .is_err()
        );
    }

    #[test]
    fn retired_workspace_rejects_progress_reads_and_wakes_waits() {
        let workspace = PlanWorkspace::new(Vec::new());
        let first = workspace.progress().unwrap();
        workspace.retire();
        assert!(workspace.progress().is_err());
        assert!(
            workspace
                .wait_for_progress(first.progress_epoch, Duration::ZERO)
                .is_err()
        );
    }

    #[test]
    fn absent_readback_worker_reaches_group_deadline_and_cancels_last_request() {
        let workspace = PlanWorkspace::new(vec![sample_block(), sample_block()]);
        // The first block remains pinned and can be acquired immediately; the
        // second is durable but has no backend to restore it for this group.
        let first = workspace.block(0).unwrap();
        for index in 0..2 {
            let work = workspace.pull().unwrap();
            assert_eq!(work.block, index);
            workspace.stored(index).unwrap();
            drop(work);
        }
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let waiters: Vec<_> = workspace
            .lock()
            .slots
            .iter()
            .map(|slot| Arc::clone(&slot.waiters))
            .collect();
        let deadline = Instant::now() + Duration::from_millis(10);
        let error = workspace
            .load_blocks_inner(&[0, 1], None, deadline)
            .err()
            .unwrap();
        assert!(error.contains("resource deadline exceeded"));
        assert!(Instant::now() >= deadline);
        assert!(
            waiters
                .iter()
                .all(|count| count.load(Ordering::Acquire) == 0)
        );
        let snapshot = workspace.snapshot();
        assert_eq!(snapshot.waiting_clients, 0);
        assert_eq!(snapshot.requested_blocks, 0);
        assert_eq!(snapshot.loading_blocks, 0);
        assert!(
            workspace
                .lock()
                .slots
                .iter()
                .all(|slot| slot.request_credit.is_none())
        );
        assert!(
            workspace.next_load().is_none(),
            "abandoned read must not reach backend"
        );
        assert!(Arc::ptr_eq(&first, &workspace.block(0).unwrap()));
    }

    #[test]
    fn incremental_directories_find_tail_blocks_after_a_large_saved_history() {
        const HISTORY: usize = 8192;
        let mut blocks: Vec<_> = (0..HISTORY).map(|_| TermBlock::default()).collect();
        blocks.push(sample_block());
        let encoded = bytes(blocks.last().unwrap());
        let workspace = PlanWorkspace::new(blocks);
        for index in 0..=HISTORY {
            let work = workspace.pull().expect("next unsaved block");
            assert_eq!(work.block, index);
            workspace.stored(index).unwrap();
            drop(work);
        }
        assert!(workspace.pull().is_none());
        assert_eq!(workspace.lock().save_cursor, HISTORY + 1);
        assert!(workspace.next_load().is_none());
        evict_unpinned_fixture(&workspace);
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);

        workspace.reserve_blocks(1, |_| {}).unwrap();
        let tail = workspace.push_block(sample_block()).unwrap();
        assert_eq!(tail, HISTORY + 1);
        let work = workspace.pull().unwrap();
        assert_eq!(work.block, tail);
        workspace.stored(tail).unwrap();
        drop(work);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        assert!(workspace.pull().is_none());
        assert_eq!(workspace.lock().save_cursor, HISTORY + 2);

        let waiting = workspace.clone();
        let client = std::thread::spawn(move || {
            waiting.load_blocks_inner(
                &[tail, HISTORY, tail],
                None,
                Instant::now() + Duration::from_secs(5),
            )
        });
        for index in [HISTORY, tail] {
            let load = next_readback(&workspace);
            assert_eq!(load.block, index);
            assert!(workspace.next_load().is_none(), "loading ID queued twice");
            workspace
                .restore(index, &mut Cursor::new(&encoded))
                .unwrap();
        }
        let loaded = client.join().unwrap().unwrap();
        assert_eq!(loaded.len(), 2);
        assert_eq!(
            loaded.iter().map(|(index, _)| *index).collect::<Vec<_>>(),
            [HISTORY, tail]
        );
        for (_, block) in &loaded {
            assert_eq!(bytes(block), encoded);
        }
        for _ in 0..10 {
            assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
            let state = workspace.lock();
            assert_eq!(state.eviction.len(), 2, "pinned candidates duplicated");
            assert!(state.loads.is_empty());
            assert!(state.slots.iter().all(|slot| !slot.load_queued));
        }
        drop(loaded);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        assert!(workspace.lock().eviction.is_empty());
        assert_eq!(workspace.snapshot().resident_bytes, 0);
    }

    #[test]
    fn cancelled_read_reuses_its_single_queue_entry_then_clears_all_flags() {
        let block = sample_block();
        let encoded = bytes(&block);
        let workspace = PlanWorkspace::new(vec![block]);
        let work = workspace.pull().unwrap();
        workspace.stored(0).unwrap();
        drop(work);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let cancel_queued_read = || {
            let retry_deadline = Instant::now() + Duration::from_secs(5);
            loop {
                let before = workspace.lock().progress_epoch;
                let error = workspace
                    .block_inner(0, None, Instant::now() + Duration::from_millis(1))
                    .err()
                    .unwrap();
                assert!(error.contains("resource deadline"));
                let state = workspace.lock();
                assert!(!state.slots[0].requested);
                assert_eq!(state.slots[0].waiters.load(Ordering::Acquire), 0);
                if state.progress_epoch == before {
                    // Scheduling can consume 1 ms before the request is even
                    // enqueued. Exercise an actual request/cancellation cycle
                    // before asserting its retained queue entry.
                    assert!(Instant::now() < retry_deadline, "request never enqueued");
                    continue;
                }
                assert_eq!(state.loads.len(), 1);
                assert!(state.slots[0].load_queued);
                break;
            }
        };
        for _ in 0..3 {
            cancel_queued_read();
        }
        let waiting = workspace.clone();
        let client = std::thread::spawn(move || waiting.block(0));
        let load = next_readback(&workspace);
        assert_eq!(load.block, 0);
        assert!(workspace.next_load().is_none());
        workspace.restore(0, &mut Cursor::new(&encoded)).unwrap();
        let loaded = client.join().unwrap().unwrap();
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
        drop(loaded);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        cancel_queued_read();
        assert!(workspace.next_load().is_none());
        {
            let state = workspace.lock();
            assert!(state.loads.is_empty());
            assert!(state.eviction.is_empty());
            assert!(!state.slots[0].load_queued);
            assert!(!state.slots[0].eviction_queued);
        }
        cancel_queued_read();
        workspace.fail("stop queued load".into());
        {
            let state = workspace.lock();
            assert!(state.loads.is_empty());
            assert!(!state.slots[0].load_queued);
            assert!(!state.slots[0].requested);
        }
        workspace.retire();
        let state = workspace.lock();
        assert!(state.slots.is_empty());
        assert!(state.loads.is_empty());
        assert!(state.eviction.is_empty());
    }

    #[test]
    fn eviction_rotates_only_a_bounded_candidate_window_and_keeps_waiters() {
        let workspace =
            PlanWorkspace::new((0..EVICTION_VISITS + 3).map(|_| sample_block()).collect());
        // Establish already-restored owners independently of read timing; the
        // preceding tests exercise the actual restore-to-candidate transition.
        let waiter = {
            let mut state = workspace.lock();
            for index in 0..state.slots.len() {
                state.slots[index].stored = true;
                state.slots[index].eviction_queued = true;
                state.eviction.push_back(index);
            }
            state.slots[0].waiters.store(1, Ordering::Release);
            Arc::clone(&state.slots[0].waiters)
        };
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        {
            let state = workspace.lock();
            assert!(state.slots[0].resident.is_some());
            assert_eq!(
                state
                    .slots
                    .iter()
                    .filter(|slot| slot.resident.is_some())
                    .count(),
                4
            );
            assert_eq!(state.eviction.len(), 4);
        }
        waiter.store(0, Ordering::Release);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let state = workspace.lock();
        assert!(state.eviction.is_empty());
        assert!(
            state
                .slots
                .iter()
                .all(|slot| slot.resident.is_none() && !slot.eviction_queued)
        );
    }

    #[test]
    fn expired_group_deadline_is_not_renewed_for_a_resident_block() {
        let workspace = PlanWorkspace::new(vec![sample_block(), sample_block()]);
        let deadline = Instant::now();
        let error = workspace
            .load_blocks_inner(&[0, 1], None, deadline)
            .err()
            .unwrap();
        assert!(error.contains("resource deadline exceeded"));
        assert_eq!(workspace.snapshot().waiting_clients, 0);
        assert_eq!(workspace.snapshot().requested_blocks, 0);
        assert_eq!(
            workspace.snapshot().resident_bytes,
            workspace.snapshot().source_owned_bytes
        );
        workspace.fail("readback stopped".into());
        assert_eq!(
            workspace.block_inner(0, None, deadline).err().unwrap(),
            "readback stopped"
        );
        workspace.retire();
        assert_eq!(workspace.snapshot().waiting_clients, 0);
    }

    #[test]
    fn requested_readback_restores_before_reclaim_and_failure_wakes_waiters() {
        let workspace = PlanWorkspace::new(vec![sample_block()]);
        let work = workspace.pull().unwrap();
        let mut encoded = Vec::new();
        work.payload.encode(&mut encoded).unwrap();
        workspace.stored(0).unwrap();
        drop(work);
        assert!(workspace.reclaim_stored(usize::MAX) > 0);
        let waiting = Arc::clone(&workspace);
        let thread = std::thread::spawn(move || waiting.block(0));
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        let request = loop {
            if let Some(request) = workspace.next_load() {
                break request;
            }
            assert!(std::time::Instant::now() < deadline);
            std::thread::yield_now();
        };
        assert_eq!(request.block, 0);
        workspace.restore(0, &mut Cursor::new(encoded)).unwrap();
        let restored = thread.join().unwrap().unwrap();
        assert_eq!(workspace.reclaim_stored(usize::MAX), 0);
        drop(restored);
        assert!(workspace.reclaim_stored(1) > 0);
        let waiting = Arc::clone(&workspace);
        let thread = std::thread::spawn(move || waiting.block(0));
        workspace.fail("readback failed".into());
        assert_eq!(thread.join().unwrap().err().unwrap(), "readback failed");
    }

    #[test]
    fn retirement_releases_owner_data_and_preserves_already_pinned_blocks() {
        let workspace = PlanWorkspace::new(vec![sample_block()]);
        let block = workspace.block(0).unwrap();
        workspace.retire();
        assert!(workspace.is_retired());
        assert!(workspace.block(0).is_err());
        assert!(workspace.pull().is_none());
        let values = [Q::from(2), Q::from(3), Q::from(5)];
        assert_ne!(
            block.evaluate(0..block.len(), &Q::ZERO, None, &values),
            Q::ZERO
        );
    }
}
