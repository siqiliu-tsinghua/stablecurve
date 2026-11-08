//! Process-wide memory observations and pressure policy.
//!
//! Call `check` only at existing planning/expansion/allocation boundaries. It
//! reads a cached observation and never calls the OS. A storage/management
//! worker calls `refresh_if_due`; initialization and explicit diagnostic slow
//! paths may use `force_check`. This module starts no thread, performs no SQL,
//! and never acquires a mathematical result-source lock.
//!
//! Thresholds are estimates, not RSS guarantees. Shared reservation groups
//! serialize declared batch envelopes; callers must ensure growth before
//! allocating and include all queued jobs, replies and arithmetic scratch.
//! Native observations are not atomic, Linux accounting differs from Windows
//! private commit, and allocations may grow between observations. In particular,
//! a plain check does not reserve headroom. Only try_reserve and group growth
//! do so; settlement keeps the full envelope charged until native sampling.
//! Environment observations cover the visible Linux memory-controller tree or
//! Windows global memory/commit and the visible current job. Hidden cgroup/job
//! ancestors and address-space rlimits are not represented by these counters.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant};

const MIB: u64 = 1024 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PrivateMemoryKind {
    /// Windows PROCESS_MEMORY_COUNTERS_EX::PrivateUsage.
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    Committed,
    /// Linux RssAnon + VmSwap; excludes unpopulated address-space reservations.
    #[cfg(target_os = "linux")]
    AnonymousResidentAndSwap,
    /// Apple's physical-footprint ledger, including compressed accounting.
    #[cfg(target_os = "macos")]
    PhysicalFootprint,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct MemorySnapshot {
    pub(crate) captured_at: Instant,
    pub(crate) process_private_bytes: Option<u64>,
    /// Native peak private commit where available (Windows), not sampled RSS.
    pub(crate) process_peak_commit_bytes: Option<u64>,
    pub(crate) process_private_kind: Option<PrivateMemoryKind>,
    pub(crate) process_resident_bytes: Option<u64>,
    pub(crate) total_physical_bytes: Option<u64>,
    pub(crate) available_physical_bytes: Option<u64>,
    pub(crate) commit_headroom_bytes: Option<u64>,
    /// Linux CommitLimit is enforced only with vm.overcommit_memory=2.
    pub(crate) commit_limit_enforced: bool,
    pub(crate) cgroup_headroom_bytes: Option<u64>,
    pub(crate) cgroup_limit_bytes: Option<u64>,
    pub(crate) job_process_limit_bytes: Option<u64>,
    #[cfg(any(target_os = "windows", test))]
    pub(crate) job_memory_limit_bytes: Option<u64>,
    pub(crate) job_headroom_bytes: Option<u64>,
    /// Job peak usage is an upper bound on current usage. The resulting
    /// headroom is conservative and may remain small after an earlier peak.
    #[cfg(any(target_os = "windows", test))]
    pub(crate) job_headroom_peak_based: bool,
    /// Whether this backend's documented environment checks completed; this is
    /// not a claim that every OS quota is observable (see module limitations).
    pub(crate) environment_limits_checked: bool,
}

impl MemorySnapshot {
    fn unknown() -> Self {
        Self {
            captured_at: Instant::now(),
            process_private_bytes: None,
            process_peak_commit_bytes: None,
            process_private_kind: None,
            process_resident_bytes: None,
            total_physical_bytes: None,
            available_physical_bytes: None,
            commit_headroom_bytes: None,
            commit_limit_enforced: false,
            cgroup_headroom_bytes: None,
            cgroup_limit_bytes: None,
            job_process_limit_bytes: None,
            #[cfg(any(target_os = "windows", test))]
            job_memory_limit_bytes: None,
            job_headroom_bytes: None,
            #[cfg(any(target_os = "windows", test))]
            job_headroom_peak_based: false,
            environment_limits_checked: false,
        }
    }

    pub(crate) fn capture() -> Self {
        let mut snapshot = platform::capture();
        snapshot.captured_at = Instant::now();
        snapshot
    }

    fn available_growth(&self) -> Option<u64> {
        let commit = self
            .commit_limit_enforced
            .then_some(self.commit_headroom_bytes)
            .flatten();
        minimum(
            minimum(
                minimum(self.available_physical_bytes, self.cgroup_headroom_bytes),
                self.job_headroom_bytes,
            ),
            commit,
        )
    }

    fn complete(&self) -> bool {
        self.process_private_bytes.is_some()
            && self.available_physical_bytes.is_some()
            && (!self.commit_limit_enforced || self.commit_headroom_bytes.is_some())
            && self.environment_limits_checked
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BudgetConfig {
    /// Maximum total process private usage, not an extra per-instance allowance.
    pub(crate) max_bytes: Option<u64>,
    pub(crate) warning_percent: u8,
    pub(crate) reclaim_percent: u8,
    pub(crate) system_headroom_bytes: u64,
    /// Shared space for minimum compute output, pull, encode, commit and readback.
    pub(crate) progress_reserve_bytes: u64,
    pub(crate) sample_interval: Duration,
}

impl Default for BudgetConfig {
    fn default() -> Self {
        Self {
            max_bytes: None,
            warning_percent: 85,
            reclaim_percent: 70,
            system_headroom_bytes: 1024 * MIB,
            progress_reserve_bytes: 256 * MIB,
            sample_interval: Duration::from_millis(100),
        }
    }
}

impl BudgetConfig {
    fn validate(&self) -> Result<(), String> {
        if self.warning_percent == 0 || self.warning_percent >= 100 {
            return Err("memory warning_percent must be between 1 and 99".into());
        }
        if self.reclaim_percent == 0 || self.reclaim_percent >= self.warning_percent {
            return Err("memory reclaim_percent must be positive and below warning_percent".into());
        }
        if self.progress_reserve_bytes == 0 {
            return Err("memory progress_reserve_bytes must be positive".into());
        }
        if self
            .max_bytes
            .is_some_and(|maximum| maximum <= self.progress_reserve_bytes)
        {
            return Err("memory max_bytes must exceed the shared progress reserve".into());
        }
        if self.sample_interval.is_zero() {
            return Err("memory sample_interval must be positive".into());
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pressure {
    Normal,
    Warning,
    Hard,
    /// Missing observations prevent a normal-pressure conclusion. Known limits
    /// and known overshoots are still returned, rather than discarding a cap.
    Unavailable,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Limits {
    pub(crate) hard_bytes: u64,
    pub(crate) computation_max_bytes: u64,
    pub(crate) warning_bytes: u64,
    pub(crate) reclaim_to_bytes: u64,
    pub(crate) progress_reserve_bytes: u64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct BudgetDecision {
    pub(crate) pressure: Pressure,
    pub(crate) memory: MemorySnapshot,
    pub(crate) sample_age: Duration,
    pub(crate) limits: Option<Limits>,
    /// May use the last known private counter if a refresh failed.
    pub(crate) known_private_bytes: Option<u64>,
    pub(crate) projected_private_bytes: Option<u64>,
    pub(crate) reclaim_bytes: u64,
    pub(crate) progress_at_risk: bool,
    pub(crate) observations_complete: bool,
    pub(crate) arithmetic_overflow: bool,
    pub(crate) active_reserved_bytes: u64,
    pub(crate) settled_reserved_bytes: u64,
    pub(crate) sample_generation: u64,
    pub(crate) stale: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReservationKind {
    Compute,
    Progress,
}

#[derive(Debug)]
pub(crate) enum ReservationError {
    Busy,
    Finished,
    AlreadyStarted,
    WouldGrow {
        reserved: u64,
        required: u64,
    },
    Stale(Duration),
    Unavailable,
    Pressure {
        pressure: Pressure,
        projected: Option<u64>,
        limits: Option<Limits>,
    },
    Underestimated {
        reserved: u64,
        measured: u64,
    },
    Overflow,
}

impl fmt::Display for ReservationError {
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Busy => output.write_str("memory reservation ledger is busy"),
            Self::Finished => output.write_str("memory reservation has finished"),
            Self::AlreadyStarted => {
                output.write_str("allocated memory cannot be cancelled as unused")
            }
            Self::WouldGrow { reserved, required } => write!(
                output,
                "shrinking a reservation cannot increase credit: reserved={reserved}, required={required}"
            ),
            Self::Stale(age) => write!(output, "memory observation is stale ({age:?})"),
            Self::Unavailable => output.write_str("memory observations are incomplete"),
            Self::Pressure {
                pressure,
                projected,
                limits,
            } => write!(
                output,
                "insufficient memory reservation: {pressure:?}, projected={projected:?}, limits={limits:?}"
            ),
            Self::Underestimated { reserved, measured } => write!(
                output,
                "memory batch exceeded reservation: reserved={reserved}, measured={measured}"
            ),
            Self::Overflow => output.write_str("memory reservation arithmetic overflow"),
        }
    }
}

#[derive(Clone, Copy)]
enum ReservationPhase {
    Active,
    Settled(u64),
}

struct Entry {
    bytes: u64,
    phase: ReservationPhase,
}

#[derive(Clone, Copy)]
struct SampleCutoff {
    generation: u64,
    event: u64,
}

struct State {
    config: BudgetConfig,
    snapshot: MemorySnapshot,
    #[cfg(test)]
    test_sample_age: Option<Duration>,
    last_private: Option<u64>,
    last_limits: Option<Limits>,
    pressure: Pressure,
    reservations: BTreeMap<u64, Entry>,
    next_reservation: u64,
    event: u64,
    next_generation: u64,
    sample_generation: u64,
    active_reserved: u64,
    settled_reserved: u64,
}

pub(crate) struct Budget {
    state: Arc<Mutex<State>>,
    sampling: AtomicBool,
}

static PROCESS_BUDGET: OnceLock<Budget> = OnceLock::new();

/// All cache instances must use this same process-level policy object.
pub(crate) fn process_budget() -> &'static Budget {
    PROCESS_BUDGET.get_or_init(|| {
        Budget::new(BudgetConfig::default(), MemorySnapshot::capture())
            .expect("default memory budget configuration")
    })
}

impl Budget {
    #[cfg(test)]
    pub(crate) fn for_test(config: BudgetConfig, snapshot: MemorySnapshot) -> Result<Self, String> {
        Self::new(config, snapshot)
    }

    fn new(config: BudgetConfig, snapshot: MemorySnapshot) -> Result<Self, String> {
        config.validate()?;
        Ok(Self {
            state: Arc::new(Mutex::new(State {
                config,
                last_private: snapshot.process_private_bytes,
                snapshot,
                #[cfg(test)]
                test_sample_age: None,
                last_limits: None,
                pressure: Pressure::Normal,
                reservations: BTreeMap::new(),
                next_reservation: 0,
                event: 0,
                next_generation: 0,
                sample_generation: 0,
                active_reserved: 0,
                settled_reserved: 0,
            })),
            sampling: AtomicBool::new(false),
        })
    }

    /// Explicit process-wide configuration; no environment is read by checks.
    pub(crate) fn configure(&self, config: BudgetConfig) -> Result<(), String> {
        config.validate()?;
        let mut state = self.state.lock().unwrap();
        state.config = config;
        state.last_limits = None;
        state.pressure = Pressure::Normal;
        Ok(())
    }

    /// Boundary-only cached estimate. Does not sample, reclaim, reserve, or wait.
    pub(crate) fn check(&self, additional_bytes: u64) -> BudgetDecision {
        self.state.lock().unwrap().decision(additional_bytes)
    }

    /// One token covers an entire dispatched wave, including queued jobs and
    /// completed replies. This never samples, reclaims, or waits for memory.
    /// Progress callers also avoid waiting for the ledger mutex.
    pub(crate) fn try_reserve(
        &self,
        bytes: u64,
        kind: ReservationKind,
    ) -> Result<Reservation, ReservationError> {
        let mut state = match kind {
            ReservationKind::Compute => {
                self.state.lock().unwrap_or_else(|error| error.into_inner())
            }
            ReservationKind::Progress => self.state.try_lock().map_err(|error| match error {
                std::sync::TryLockError::WouldBlock => ReservationError::Busy,
                std::sync::TryLockError::Poisoned(_) => ReservationError::Unavailable,
            })?,
        };
        state.approve(bytes, kind)?;
        let id = state.next_reservation;
        state.next_reservation = id.checked_add(1).ok_or(ReservationError::Overflow)?;
        state.active_reserved = state
            .active_reserved
            .checked_add(bytes)
            .ok_or(ReservationError::Overflow)?;
        state.reservations.insert(
            id,
            Entry {
                bytes,
                phase: ReservationPhase::Active,
            },
        );
        drop(state);
        Ok(Reservation(Arc::new(ReservationGroup {
            ledger: Arc::clone(&self.state),
            id,
            kind,
            state: Mutex::new(GroupState {
                reserved: bytes,
                started: false,
                finished: false,
            }),
        })))
    }

    /// Called by a management/storage worker without source/DB locks held.
    /// A concurrent sampler leaves the previous observation available to callers.
    pub(crate) fn refresh_if_due(&self) -> BudgetDecision {
        self.refresh(false);
        self.check(0)
    }

    /// Native sampling is allowed only at initialization or an explicit slow
    /// diagnostic/pressure boundary, never in a mathematical item loop.
    pub(crate) fn force_check(&self, additional_bytes: u64) -> BudgetDecision {
        self.refresh(true);
        self.check(additional_bytes)
    }

    /// Publish a background sample, also allowing deterministic boundary tests.
    #[cfg(test)]
    pub(crate) fn publish_snapshot(&self, snapshot: MemorySnapshot) {
        let mut state = self.state.lock().unwrap();
        if let Some(private) = snapshot.process_private_bytes {
            state.last_private = Some(private);
        }
        state.snapshot = snapshot;
        state.test_sample_age = None;
        // An externally captured sample has no proven capture-start cutoff.
        // Retaining liabilities is conservative; only refresh can retire them.
    }

    /// Synthetic age keeps tests independent of the platform clock's lower bound.
    #[cfg(test)]
    pub(crate) fn set_sample_age_for_test(&self, age: Duration) {
        self.state.lock().unwrap().test_sample_age = Some(age);
    }

    fn refresh(&self, force: bool) {
        if !force {
            let state = self.state.lock().unwrap();
            if state.sample_age() < state.config.sample_interval {
                return;
            }
        }
        if self
            .sampling
            .compare_exchange(false, true, Ordering::Acquire, Ordering::Relaxed)
            .is_err()
        {
            return;
        }
        struct Sampling<'a>(&'a AtomicBool);
        impl Drop for Sampling<'_> {
            fn drop(&mut self) {
                self.0.store(false, Ordering::Release);
            }
        }
        let _sampling = Sampling(&self.sampling);
        if !force {
            let state = self.state.lock().unwrap();
            if state.sample_age() < state.config.sample_interval {
                return;
            }
        }
        let cutoff = self.state.lock().unwrap().begin_sample();
        // Native calls and /proc reads happen outside the policy mutex.
        let snapshot = MemorySnapshot::capture();
        self.state.lock().unwrap().publish_sample(snapshot, cutoff);
    }
}

struct GroupState {
    reserved: u64,
    started: bool,
    finished: bool,
}

struct ReservationGroup {
    ledger: Arc<Mutex<State>>,
    id: u64,
    kind: ReservationKind,
    state: Mutex<GroupState>,
}

/// Clone once per participating worker/job, not once per mathematical item.
/// All growth requests describe the total envelope of this same wave.
#[derive(Clone)]
pub(crate) struct Reservation(Arc<ReservationGroup>);

impl fmt::Debug for Reservation {
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        let group = self
            .0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        output
            .debug_struct("Reservation")
            .field("id", &self.0.id)
            .field("kind", &self.0.kind)
            .field("reserved", &group.reserved)
            .field("started", &group.started)
            .field("finished", &group.finished)
            .finish()
    }
}

impl Reservation {
    pub(crate) fn start(&self) -> Result<(), ReservationError> {
        let mut group = self
            .0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if group.finished {
            return Err(ReservationError::Finished);
        }
        group.started = true;
        Ok(())
    }

    pub(crate) fn reserved_bytes(&self) -> u64 {
        self.0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .reserved
    }

    /// Call before allocating, at a container growth or batch boundary. Keeping
    /// within the existing envelope only takes this group's short mutex. A
    /// larger envelope is approved atomically against every other group's use.
    /// Credit alone does not dispatch work or allocate: callers mark that
    /// transition explicitly with `start` before touching the reserved memory.
    pub(crate) fn ensure(&self, required_total: u64) -> Result<(), ReservationError> {
        let mut group = self
            .0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if group.finished {
            return Err(ReservationError::Finished);
        }
        if required_total <= group.reserved {
            return Ok(());
        }
        let additional = required_total - group.reserved;
        let mut ledger = match self.0.kind {
            ReservationKind::Compute => self
                .0
                .ledger
                .lock()
                .unwrap_or_else(|error| error.into_inner()),
            ReservationKind::Progress => self.0.ledger.try_lock().map_err(|error| match error {
                std::sync::TryLockError::WouldBlock => ReservationError::Busy,
                std::sync::TryLockError::Poisoned(_) => ReservationError::Unavailable,
            })?,
        };
        ledger.approve(additional, self.0.kind)?;
        let active = ledger
            .active_reserved
            .checked_add(additional)
            .ok_or(ReservationError::Overflow)?;
        let entry = ledger
            .reservations
            .get_mut(&self.0.id)
            .ok_or(ReservationError::Finished)?;
        entry.bytes = required_total;
        ledger.active_reserved = active;
        group.reserved = required_total;
        Ok(())
    }

    /// Return only demonstrably unused future credit. The entry stays active
    /// under the same identity, with no release-and-reacquire admission race.
    /// Unlike growth this is valid even when the native observation is stale.
    pub(crate) fn shrink_unstarted(&self, required_total: u64) -> Result<(), ReservationError> {
        let mut group = self
            .0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if group.finished {
            return Err(ReservationError::Finished);
        }
        if group.started {
            return Err(ReservationError::AlreadyStarted);
        }
        if required_total > group.reserved {
            return Err(ReservationError::WouldGrow {
                reserved: group.reserved,
                required: required_total,
            });
        }
        if required_total == group.reserved {
            return Ok(());
        }
        let mut ledger = match self.0.kind {
            ReservationKind::Compute => self
                .0
                .ledger
                .lock()
                .unwrap_or_else(|error| error.into_inner()),
            ReservationKind::Progress => self.0.ledger.try_lock().map_err(|error| match error {
                std::sync::TryLockError::WouldBlock => ReservationError::Busy,
                std::sync::TryLockError::Poisoned(_) => ReservationError::Unavailable,
            })?,
        };
        let active = ledger
            .active_reserved
            .checked_sub(group.reserved - required_total)
            .ok_or(ReservationError::Overflow)?;
        let entry = ledger
            .reservations
            .get_mut(&self.0.id)
            .ok_or(ReservationError::Finished)?;
        if !matches!(entry.phase, ReservationPhase::Active) {
            return Err(ReservationError::Finished);
        }
        entry.bytes = required_total;
        ledger.active_reserved = active;
        group.reserved = required_total;
        Ok(())
    }

    /// All jobs/replies/store allocations must have finished before settlement.
    /// Actual measurement detects a broken estimate but cannot make a prior
    /// unreserved allocation safe. The full envelope stays charged until a new
    /// sample covers it, including temporary allocator pages retained after use.
    pub(crate) fn settle(&self, measured_bytes: u64) -> Result<(), ReservationError> {
        let mut group = self
            .0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if group.finished {
            return Err(ReservationError::Finished);
        }
        let reserved = group.reserved;
        self.0
            .ledger
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .settle_entry(self.0.id, measured_bytes.max(reserved))?;
        group.finished = true;
        if measured_bytes > reserved {
            return Err(ReservationError::Underestimated {
                reserved,
                measured: measured_bytes,
            });
        }
        Ok(())
    }

    pub(crate) fn cancel_unused(&self) -> Result<(), ReservationError> {
        let mut group = self
            .0
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if group.finished {
            return Err(ReservationError::Finished);
        }
        if group.started {
            return Err(ReservationError::AlreadyStarted);
        }
        self.0
            .ledger
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .cancel_entry(self.0.id)?;
        group.finished = true;
        Ok(())
    }
}

impl Drop for ReservationGroup {
    fn drop(&mut self) {
        let group = self
            .state
            .get_mut()
            .unwrap_or_else(|error| error.into_inner());
        if group.finished {
            return;
        }
        let mut ledger = self
            .ledger
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        // Failure/unwind after dispatch cannot return potentially materialized
        // bytes before native sampling. Never wait for budget progress in Drop.
        let _ = if group.started {
            ledger.settle_entry(self.id, group.reserved)
        } else {
            ledger.cancel_entry(self.id)
        };
    }
}

impl State {
    fn sample_age(&self) -> Duration {
        #[cfg(test)]
        if let Some(age) = self.test_sample_age {
            return age;
        }
        self.snapshot.captured_at.elapsed()
    }

    fn begin_sample(&mut self) -> SampleCutoff {
        self.next_generation = self
            .next_generation
            .checked_add(1)
            .expect("memory sampling generation exhausted");
        SampleCutoff {
            generation: self.next_generation,
            event: self.event,
        }
    }

    fn publish_sample(&mut self, snapshot: MemorySnapshot, cutoff: SampleCutoff) {
        if cutoff.generation < self.sample_generation {
            return;
        }
        if let Some(private) = snapshot.process_private_bytes {
            self.last_private = Some(private);
            let mut absorbed = 0u64;
            self.reservations.retain(|_, entry| {
                if matches!(entry.phase, ReservationPhase::Settled(event) if event <= cutoff.event)
                {
                    absorbed = absorbed.saturating_add(entry.bytes);
                    false
                } else {
                    true
                }
            });
            self.settled_reserved = self.settled_reserved.saturating_sub(absorbed);
        }
        self.sample_generation = cutoff.generation;
        self.snapshot = snapshot;
        #[cfg(test)]
        {
            self.test_sample_age = None;
        }
    }

    fn settle_entry(&mut self, id: u64, bytes: u64) -> Result<(), ReservationError> {
        let entry = self
            .reservations
            .get(&id)
            .ok_or(ReservationError::Finished)?;
        if !matches!(entry.phase, ReservationPhase::Active) {
            return Err(ReservationError::Finished);
        }
        let active = self
            .active_reserved
            .checked_sub(entry.bytes)
            .ok_or(ReservationError::Overflow)?;
        let settled = self
            .settled_reserved
            .checked_add(bytes)
            .ok_or(ReservationError::Overflow)?;
        let event = self
            .event
            .checked_add(1)
            .ok_or(ReservationError::Overflow)?;
        let entry = self.reservations.get_mut(&id).unwrap();
        entry.bytes = bytes;
        entry.phase = ReservationPhase::Settled(event);
        self.event = event;
        self.active_reserved = active;
        self.settled_reserved = settled;
        Ok(())
    }

    fn cancel_entry(&mut self, id: u64) -> Result<(), ReservationError> {
        let entry = self
            .reservations
            .get(&id)
            .ok_or(ReservationError::Finished)?;
        if !matches!(entry.phase, ReservationPhase::Active) {
            return Err(ReservationError::Finished);
        }
        self.active_reserved = self
            .active_reserved
            .checked_sub(entry.bytes)
            .ok_or(ReservationError::Overflow)?;
        self.reservations.remove(&id);
        Ok(())
    }

    fn approve(&mut self, additional: u64, kind: ReservationKind) -> Result<(), ReservationError> {
        let decision = self.decision(additional);
        if decision.stale {
            return Err(ReservationError::Stale(decision.sample_age));
        }
        if decision.arithmetic_overflow {
            return Err(ReservationError::Overflow);
        }
        if !decision.observations_complete || decision.limits.is_none() {
            return Err(ReservationError::Unavailable);
        }
        let permitted = match kind {
            // Warning and its reclamation hysteresis guide the controller's
            // spill/batch policy; an irreducible floor need not reach the
            // reclaim target before the next bounded wave can make progress.
            ReservationKind::Compute => decision
                .projected_private_bytes
                .zip(decision.limits)
                .is_some_and(|(projected, limits)| projected <= limits.computation_max_bytes),
            // Critical progress may consume the reserve excluded from compute.
            ReservationKind::Progress => decision
                .projected_private_bytes
                .zip(decision.limits)
                .is_some_and(|(projected, limits)| projected < limits.hard_bytes),
        };
        if permitted {
            Ok(())
        } else {
            Err(ReservationError::Pressure {
                pressure: decision.pressure,
                projected: decision.projected_private_bytes,
                limits: decision.limits,
            })
        }
    }

    fn decision(&mut self, additional: u64) -> BudgetDecision {
        let age = self.sample_age();
        let stale = age
            > self
                .config
                .sample_interval
                .checked_mul(2)
                .unwrap_or(Duration::MAX);
        let complete = self.snapshot.complete() && !stale;
        let available = self.snapshot.available_growth();
        let mut overflow = false;
        let dynamic_max = self
            .snapshot
            .process_private_bytes
            .zip(available)
            .map(|(used, free)| {
                used.checked_add(free.saturating_sub(self.config.system_headroom_bytes))
                    .unwrap_or_else(|| {
                        overflow = true;
                        u64::MAX
                    })
            });
        let mut maximum = minimum(
            minimum(self.config.max_bytes, dynamic_max),
            self.snapshot.job_process_limit_bytes,
        );
        if !complete {
            // Failed/partial samples must not silently raise a previous cap.
            maximum = minimum(maximum, self.last_limits.map(|limits| limits.hard_bytes));
        }
        let limits = maximum.map(|hard| {
            let computation = hard.saturating_sub(self.config.progress_reserve_bytes);
            Limits {
                hard_bytes: hard,
                computation_max_bytes: computation,
                warning_bytes: percent(computation, self.config.warning_percent),
                reclaim_to_bytes: percent(computation, self.config.reclaim_percent),
                progress_reserve_bytes: self.config.progress_reserve_bytes,
            }
        });
        if let Some(limits) = limits {
            self.last_limits = Some(limits);
        }
        let private = self.snapshot.process_private_bytes.or(self.last_private);
        let liabilities = self.active_reserved.checked_add(self.settled_reserved);
        let growth = liabilities.and_then(|reserved| reserved.checked_add(additional));
        let projected = private
            .zip(growth)
            .and_then(|(private, growth)| private.checked_add(growth));
        if growth.is_none() {
            overflow = true;
        }
        if private.is_some() && projected.is_none() {
            overflow = true;
        }
        let required = self
            .config
            .system_headroom_bytes
            .checked_add(self.config.progress_reserve_bytes)
            .and_then(|reserved| growth.and_then(|growth| reserved.checked_add(growth)));
        if required.is_none() {
            overflow = true;
        }
        let growth_threat = available
            .zip(required)
            .is_some_and(|(available, required)| available <= required);
        // With no process counter, requested growth is still a lower bound;
        // an explicit small cap therefore remains actionable.
        let projected_lower_bound = projected.unwrap_or(growth.unwrap_or(u64::MAX));
        let cap_threat =
            limits.is_some_and(|limits| projected_lower_bound >= limits.computation_max_bytes);
        let progress_at_risk = overflow || growth_threat || cap_threat;
        let pressure = if progress_at_risk {
            Pressure::Hard
        } else if limits.is_some_and(|limits| projected_lower_bound >= limits.warning_bytes) {
            Pressure::Warning
        } else if !complete || limits.is_none() {
            Pressure::Unavailable
        } else if matches!(self.pressure, Pressure::Warning | Pressure::Hard)
            && limits.is_some_and(|limits| projected_lower_bound > limits.reclaim_to_bytes)
        {
            Pressure::Warning
        } else {
            Pressure::Normal
        };
        if pressure != Pressure::Unavailable {
            self.pressure = pressure;
        }
        let reclaim_bytes = if matches!(pressure, Pressure::Warning | Pressure::Hard) {
            limits.map_or(0, |limits| {
                projected_lower_bound.saturating_sub(limits.reclaim_to_bytes)
            })
        } else {
            0
        };
        BudgetDecision {
            pressure,
            memory: self.snapshot,
            sample_age: age,
            limits,
            known_private_bytes: private,
            projected_private_bytes: projected,
            reclaim_bytes,
            progress_at_risk,
            observations_complete: complete,
            arithmetic_overflow: overflow,
            active_reserved_bytes: self.active_reserved,
            settled_reserved_bytes: self.settled_reserved,
            sample_generation: self.sample_generation,
            stale,
        }
    }
}

fn minimum(left: Option<u64>, right: Option<u64>) -> Option<u64> {
    match (left, right) {
        (Some(left), Some(right)) => Some(left.min(right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn percent(value: u64, percentage: u8) -> u64 {
    // Splitting the quotient and remainder avoids value * percentage overflow.
    (value / 100)
        .checked_mul(u64::from(percentage))
        .and_then(|whole| {
            (value % 100)
                .checked_mul(u64::from(percentage))
                .map(|fraction| (whole, fraction / 100))
        })
        .and_then(|(whole, fraction)| whole.checked_add(fraction))
        .expect("percentage between zero and one hundred")
}

#[cfg(target_os = "windows")]
mod platform {
    use super::{MemorySnapshot, PrivateMemoryKind};
    use std::ffi::c_void;

    #[repr(C)]
    #[derive(Default)]
    struct MemoryStatus {
        length: u32,
        load: u32,
        total_physical: u64,
        available_physical: u64,
        total_page_file: u64,
        available_page_file: u64,
        total_virtual: u64,
        available_virtual: u64,
        available_extended_virtual: u64,
    }

    #[repr(C)]
    #[derive(Default)]
    struct ProcessCounters {
        cb: u32,
        page_fault_count: u32,
        peak_working_set_size: usize,
        working_set_size: usize,
        quota_peak_paged_pool_usage: usize,
        quota_paged_pool_usage: usize,
        quota_peak_non_paged_pool_usage: usize,
        quota_non_paged_pool_usage: usize,
        pagefile_usage: usize,
        peak_pagefile_usage: usize,
        private_usage: usize,
    }

    #[repr(C)]
    #[derive(Default)]
    struct PerformanceInfo {
        cb: u32,
        commit_total: usize,
        commit_limit: usize,
        commit_peak: usize,
        physical_total: usize,
        physical_available: usize,
        system_cache: usize,
        kernel_total: usize,
        kernel_paged: usize,
        kernel_non_paged: usize,
        page_size: usize,
        handle_count: u32,
        process_count: u32,
        thread_count: u32,
    }

    #[repr(C)]
    #[derive(Default)]
    struct JobBasicLimits {
        process_user_time: i64,
        job_user_time: i64,
        flags: u32,
        minimum_working_set: usize,
        maximum_working_set: usize,
        active_processes: u32,
        affinity: usize,
        priority: u32,
        scheduling: u32,
    }

    #[repr(C)]
    #[derive(Default)]
    struct JobExtendedLimits {
        basic: JobBasicLimits,
        io_counters: [u64; 6],
        process_memory: usize,
        job_memory: usize,
        peak_process_memory: usize,
        peak_job_memory: usize,
    }

    #[repr(C)]
    #[derive(Default)]
    struct JobMemoryUsage {
        job_memory: u64,
        peak_job_memory: u64,
    }

    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetCurrentProcess() -> *mut c_void;
        fn GlobalMemoryStatusEx(status: *mut MemoryStatus) -> i32;
        fn IsProcessInJob(process: *mut c_void, job: *mut c_void, result: *mut i32) -> i32;
        fn QueryInformationJobObject(
            job: *mut c_void,
            class: i32,
            information: *mut c_void,
            length: u32,
            returned: *mut u32,
        ) -> i32;
    }

    #[link(name = "psapi")]
    unsafe extern "system" {
        fn GetProcessMemoryInfo(
            process: *mut c_void,
            counters: *mut ProcessCounters,
            cb: u32,
        ) -> i32;
        fn GetPerformanceInfo(info: *mut PerformanceInfo, cb: u32) -> i32;
    }

    pub(super) fn capture() -> MemorySnapshot {
        let mut snapshot = MemorySnapshot::unknown();
        snapshot.environment_limits_checked = true;
        snapshot.commit_limit_enforced = true;
        let mut status = MemoryStatus {
            length: std::mem::size_of::<MemoryStatus>() as u32,
            ..MemoryStatus::default()
        };
        // SAFETY: MemoryStatus has the MEMORYSTATUSEX C ABI and initialized size.
        if unsafe { GlobalMemoryStatusEx(&mut status) } != 0 {
            snapshot.total_physical_bytes = Some(status.total_physical);
            snapshot.available_physical_bytes = Some(status.available_physical);
        }
        let mut counters = ProcessCounters {
            cb: std::mem::size_of::<ProcessCounters>() as u32,
            ..ProcessCounters::default()
        };
        // SAFETY: ProcessCounters is PROCESS_MEMORY_COUNTERS_EX; the current
        // process pseudo-handle is valid and must not be closed.
        if unsafe { GetProcessMemoryInfo(GetCurrentProcess(), &mut counters, counters.cb) } != 0 {
            snapshot.process_private_bytes = Some(counters.private_usage as u64);
            snapshot.process_peak_commit_bytes = Some(counters.peak_pagefile_usage as u64);
            snapshot.process_private_kind = Some(PrivateMemoryKind::Committed);
            snapshot.process_resident_bytes = Some(counters.working_set_size as u64);
        }
        let mut performance = PerformanceInfo {
            cb: std::mem::size_of::<PerformanceInfo>() as u32,
            ..PerformanceInfo::default()
        };
        // SAFETY: PerformanceInfo has PERFORMANCE_INFORMATION's C ABI and size.
        if unsafe { GetPerformanceInfo(&mut performance, performance.cb) } != 0
            && performance.page_size != 0
        {
            let page_size = performance.page_size as u64;
            snapshot.commit_headroom_bytes = (performance
                .commit_limit
                .saturating_sub(performance.commit_total)
                as u64)
                .checked_mul(page_size);
            if snapshot.total_physical_bytes.is_none() {
                snapshot.total_physical_bytes =
                    (performance.physical_total as u64).checked_mul(page_size);
            }
            if snapshot.available_physical_bytes.is_none() {
                snapshot.available_physical_bytes =
                    (performance.physical_available as u64).checked_mul(page_size);
            }
        }
        job_limits(&mut snapshot);
        snapshot
    }

    fn job_limits(snapshot: &mut MemorySnapshot) {
        let mut in_job = 0;
        // SAFETY: the pseudo-handle and BOOL output follow IsProcessInJob's ABI;
        // null selects any job associated with this process, without opening it.
        if unsafe { IsProcessInJob(GetCurrentProcess(), std::ptr::null_mut(), &mut in_job) } == 0 {
            snapshot.environment_limits_checked = false;
            return;
        }
        if in_job == 0 {
            return;
        }
        let mut limits = JobExtendedLimits::default();
        // SAFETY: this is JOBOBJECT_EXTENDED_LIMIT_INFORMATION, information
        // class 9. Null queries the caller's job; the full output is initialized.
        if unsafe {
            QueryInformationJobObject(
                std::ptr::null_mut(),
                9,
                (&mut limits as *mut JobExtendedLimits).cast(),
                std::mem::size_of::<JobExtendedLimits>() as u32,
                std::ptr::null_mut(),
            )
        } == 0
        {
            snapshot.environment_limits_checked = false;
            return;
        }
        if limits.basic.flags & 0x100 != 0 {
            snapshot.job_process_limit_bytes = Some(limits.process_memory as u64);
        }
        if limits.basic.flags & 0x200 != 0 {
            snapshot.job_memory_limit_bytes = Some(limits.job_memory as u64);
            let mut usage = JobMemoryUsage::default();
            // SAFETY: JOBOBJECT_MEMORY_USAGE_INFORMATION has two ULONG64 fields;
            // class 28 reports current usage. Microsoft uses this layout in
            // hcsshim/internal/winapi/jobobject.go. Older systems may reject it.
            let current = unsafe {
                QueryInformationJobObject(
                    std::ptr::null_mut(),
                    28,
                    (&mut usage as *mut JobMemoryUsage).cast(),
                    std::mem::size_of::<JobMemoryUsage>() as u32,
                    std::ptr::null_mut(),
                )
            } != 0;
            let used = if current {
                usage.job_memory
            } else {
                limits.peak_job_memory as u64
            };
            snapshot.job_headroom_bytes = Some((limits.job_memory as u64).saturating_sub(used));
            snapshot.job_headroom_peak_based = !current;
        }
    }
}

#[cfg(target_os = "linux")]
mod platform {
    use super::{MemorySnapshot, PrivateMemoryKind, minimum};
    use std::ffi::OsString;
    use std::fs;
    use std::io::ErrorKind;
    use std::os::unix::ffi::OsStringExt;
    use std::path::{Component, Path, PathBuf};

    fn kibibytes(text: &str, key: &str) -> Option<u64> {
        text.lines()
            .find_map(|line| line.strip_prefix(key))?
            .split_whitespace()
            .next()?
            .parse::<u64>()
            .ok()?
            .checked_mul(1024)
    }

    pub(super) fn capture() -> MemorySnapshot {
        let mut snapshot = MemorySnapshot::unknown();
        if let Ok(status) = fs::read_to_string("/proc/self/status") {
            snapshot.process_resident_bytes = kibibytes(&status, "VmRSS:");
            snapshot.process_private_bytes = kibibytes(&status, "RssAnon:")
                .zip(kibibytes(&status, "VmSwap:"))
                .and_then(|(anonymous, swapped)| anonymous.checked_add(swapped));
            if snapshot.process_private_bytes.is_some() {
                snapshot.process_private_kind = Some(PrivateMemoryKind::AnonymousResidentAndSwap);
            }
        }
        if let Ok(memory) = fs::read_to_string("/proc/meminfo") {
            snapshot.total_physical_bytes = kibibytes(&memory, "MemTotal:");
            snapshot.available_physical_bytes = kibibytes(&memory, "MemAvailable:");
            snapshot.commit_headroom_bytes = kibibytes(&memory, "CommitLimit:")
                .zip(kibibytes(&memory, "Committed_AS:"))
                .map(|(limit, used)| limit.saturating_sub(used));
        }
        // In modes 0/1 this is advisory and may already be over CommitLimit.
        // Treating that as a hard limit would incorrectly stop normal Linux jobs.
        let overcommit = fs::read_to_string("/proc/sys/vm/overcommit_memory")
            .ok()
            .and_then(|value| value.trim().parse::<u8>().ok());
        snapshot.commit_limit_enforced = overcommit == Some(2);
        let cgroup = cgroup_memory();
        snapshot.cgroup_headroom_bytes = cgroup.headroom;
        snapshot.cgroup_limit_bytes = cgroup.limit;
        snapshot.environment_limits_checked = cgroup.checked && overcommit.is_some();
        snapshot
    }

    #[derive(Default)]
    struct CgroupMemory {
        headroom: Option<u64>,
        limit: Option<u64>,
        checked: bool,
    }

    // proc mount paths encode whitespace/backslashes as three octal digits.
    fn proc_path(text: &str) -> PathBuf {
        let input = text.as_bytes();
        let mut bytes = Vec::with_capacity(input.len());
        let mut index = 0;
        while index < input.len() {
            if input[index] == b'\\' && index + 3 < input.len() {
                let digits = &input[index + 1..index + 4];
                if digits.iter().all(|digit| (b'0'..=b'7').contains(digit)) {
                    let value = u16::from(digits[0] - b'0') * 64
                        + u16::from(digits[1] - b'0') * 8
                        + u16::from(digits[2] - b'0');
                    if let Ok(value) = u8::try_from(value) {
                        bytes.push(value);
                        index += 4;
                        continue;
                    }
                }
            }
            bytes.push(input[index]);
            index += 1;
        }
        PathBuf::from(OsString::from_vec(bytes))
    }

    fn group_directory(group: &Path, root: &Path, mount: &Path) -> Option<PathBuf> {
        if !group.is_absolute()
            || !root.is_absolute()
            || !mount.is_absolute()
            || mount.components().any(|part| part == Component::ParentDir)
        {
            return None;
        }
        if let Ok(relative) = group.strip_prefix(root) {
            // Both proc paths are relative to the cgroup namespace root, so
            // leading `..` components are valid (e.g. Colab's /../../...).
            // Remove the mount's root first; only the remaining suffix must
            // stay inside the mounted hierarchy.
            if relative
                .components()
                .any(|part| part == Component::ParentDir)
            {
                return None;
            }
            Some(mount.join(relative))
        } else if group == Path::new("/")
            && !root.components().any(|part| part == Component::ParentDir)
        {
            // A cgroup namespace can expose its current group as root while
            // mountinfo still names the host-side subtree root.
            Some(mount.to_path_buf())
        } else {
            None
        }
    }

    fn inspect_hierarchy(directory: &Path, mount: &Path, unified: bool) -> CgroupMemory {
        let (limit_file, usage_file) = if unified {
            ("memory.max", "memory.current")
        } else {
            ("memory.limit_in_bytes", "memory.usage_in_bytes")
        };
        let mut result = CgroupMemory {
            checked: true,
            ..CgroupMemory::default()
        };
        if !directory.is_dir() {
            return CgroupMemory::default();
        }
        for ancestor in directory.ancestors() {
            if !ancestor.starts_with(mount) {
                break;
            }
            let limit = fs::read_to_string(ancestor.join(limit_file));
            let usage = fs::read_to_string(ancestor.join(usage_file));
            match (limit, usage) {
                (Ok(limit), Ok(usage)) => {
                    if unified && limit.trim() == "max" {
                        continue;
                    }
                    match (limit.trim().parse::<u64>(), usage.trim().parse::<u64>()) {
                        (Ok(limit), Ok(used)) => {
                            result.limit = minimum(result.limit, Some(limit));
                            result.headroom =
                                minimum(result.headroom, Some(limit.saturating_sub(used)));
                        }
                        _ => result.checked = false,
                    }
                }
                // The controller can be absent at a node, and root cgroups
                // can omit their unlimited max file. Ancestors still matter.
                (Err(error), Ok(_)) if error.kind() == ErrorKind::NotFound => {}
                (Err(limit), Err(usage))
                    if limit.kind() == ErrorKind::NotFound
                        && usage.kind() == ErrorKind::NotFound => {}
                _ => result.checked = false,
            }
            if ancestor == mount {
                break;
            }
        }
        result
    }

    fn cgroup_memory() -> CgroupMemory {
        let Ok(membership) = fs::read_to_string("/proc/self/cgroup") else {
            return CgroupMemory::default();
        };
        let Ok(mountinfo) = fs::read_to_string("/proc/self/mountinfo") else {
            return CgroupMemory::default();
        };
        cgroup_memory_from_proc(&membership, &mountinfo)
    }

    fn cgroup_memory_from_proc(membership: &str, mountinfo: &str) -> CgroupMemory {
        let mut result = CgroupMemory {
            checked: true,
            ..CgroupMemory::default()
        };
        for member in membership.lines() {
            let Some((_, rest)) = member.split_once(':') else {
                continue;
            };
            let Some((controllers, group)) = rest.split_once(':') else {
                continue;
            };
            let unified = controllers.is_empty();
            if !unified && !controllers.split(',').any(|name| name == "memory") {
                continue;
            }
            let group = proc_path(group);
            let mut matched = false;
            for line in mountinfo.lines() {
                let Some((before, after)) = line.split_once(" - ") else {
                    continue;
                };
                let mut fields = after.split_whitespace();
                let Some(filesystem) = fields.next() else {
                    continue;
                };
                let _source = fields.next();
                let options = fields.next().unwrap_or("");
                if (unified && filesystem != "cgroup2")
                    || (!unified
                        && (filesystem != "cgroup"
                            || !options.split(',').any(|name| name == "memory")))
                {
                    continue;
                }
                let mut fields = before.split_whitespace().skip(3);
                let Some(root) = fields.next() else { continue };
                let Some(mount) = fields.next() else { continue };
                let root = proc_path(root);
                let mount = proc_path(mount);
                let Some(directory) = group_directory(&group, &root, &mount) else {
                    continue;
                };
                let observed = inspect_hierarchy(&directory, &mount, unified);
                result.limit = minimum(result.limit, observed.limit);
                result.headroom = minimum(result.headroom, observed.headroom);
                result.checked &= observed.checked;
                matched = true;
            }
            result.checked &= matched;
        }
        result
    }

    #[cfg(test)]
    #[test]
    fn proc_values_and_mount_roots_are_parsed_without_host_root_assumptions() {
        assert_eq!(
            kibibytes("VmRSS: 12 kB\nRssAnon: 7 kB\n", "RssAnon:"),
            Some(7168)
        );
        assert_eq!(
            proc_path("/sys/fs/cgroup\\040space"),
            Path::new("/sys/fs/cgroup space")
        );
        assert_eq!(
            group_directory(
                Path::new("/parent/job"),
                Path::new("/parent"),
                Path::new("/cg")
            ),
            Some(PathBuf::from("/cg/job"))
        );
        assert!(group_directory(Path::new("/../job"), Path::new("/"), Path::new("/cg")).is_none());
    }

    #[cfg(test)]
    #[test]
    fn cgroup_namespace_parents_are_resolved_before_mount_containment_checks() {
        let mount = Path::new("/cg");
        assert_eq!(
            group_directory(
                Path::new("/../../jupyter-children"),
                Path::new("/../.."),
                mount
            ),
            Some(mount.join("jupyter-children"))
        );
        assert_eq!(
            group_directory(Path::new("/../.."), Path::new("/../.."), mount),
            Some(mount.to_path_buf())
        );
        for (group, root) in [
            ("/../outside", "/"),
            ("/../../parent2/job", "/../../parent"),
            ("/../../parent/../outside", "/../../parent"),
            ("/", "/.."),
            ("relative/job", "/"),
            ("/job", "relative"),
        ] {
            assert!(group_directory(Path::new(group), Path::new(root), mount).is_none());
        }
        assert!(group_directory(Path::new("/job"), Path::new("/"), Path::new("/cg/..")).is_none());
    }

    #[cfg(test)]
    struct TestHierarchy(PathBuf);

    #[cfg(test)]
    impl Drop for TestHierarchy {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[cfg(test)]
    fn test_hierarchy() -> TestHierarchy {
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        TestHierarchy(
            std::env::temp_dir().join(format!("sc-faber-cgroup-{}-{stamp}", std::process::id())),
        )
    }

    #[cfg(test)]
    fn mount_field(path: &Path) -> String {
        path.to_str()
            .unwrap()
            .replace('\\', "\\134")
            .replace(' ', "\\040")
            .replace('\t', "\\011")
            .replace('\n', "\\012")
    }

    #[cfg(test)]
    #[test]
    fn cgroup_namespace_parent_paths_still_enforce_visible_ancestor_limits() {
        let hierarchy = test_hierarchy();
        let mount = &hierarchy.0;
        let child = mount.join("jupyter-children");
        fs::create_dir_all(&child).unwrap();
        let mount_field = mount_field(mount);
        for (membership, filesystem, options, limit_file, usage_file) in [
            (
                "0::/../../jupyter-children\n",
                "cgroup2",
                "rw",
                "memory.max",
                "memory.current",
            ),
            (
                "5:cpu,memory:/../../jupyter-children\n",
                "cgroup",
                "rw,cpu,memory",
                "memory.limit_in_bytes",
                "memory.usage_in_bytes",
            ),
        ] {
            fs::write(mount.join(limit_file), "90\n").unwrap();
            fs::write(mount.join(usage_file), "80\n").unwrap();
            fs::write(child.join(limit_file), "100\n").unwrap();
            fs::write(child.join(usage_file), "20\n").unwrap();
            let mountinfo =
                format!("42 21 0:42 /../.. {mount_field} rw - {filesystem} cgroup {options}\n");
            let observed = cgroup_memory_from_proc(membership, &mountinfo);
            assert!(observed.checked);
            assert_eq!(observed.limit, Some(90));
            assert_eq!(observed.headroom, Some(10));

            fs::write(child.join(limit_file), "invalid\n").unwrap();
            assert!(!cgroup_memory_from_proc(membership, &mountinfo).checked);
            let outside = membership.replace("/../../jupyter-children", "/../../../outside");
            assert!(!cgroup_memory_from_proc(&outside, &mountinfo).checked);
        }
    }

    #[cfg(test)]
    #[test]
    fn colab_dual_mounts_resolve_the_current_group_without_unrelated_limits() {
        let hierarchy = test_hierarchy();
        let unrelated = hierarchy.0.join("sys/fs/cgroup");
        let current = hierarchy.0.join("var/colab/cgroup/jupyter-children");
        for (directory, limit, used) in [(&unrelated, "10\n", "9\n"), (&current, "100\n", "30\n")] {
            fs::create_dir_all(directory).unwrap();
            fs::write(directory.join("memory.max"), limit).unwrap();
            fs::write(directory.join("memory.current"), used).unwrap();
        }
        let unrelated_line = format!(
            "936 935 0:27 / {} ro,nosuid,nodev,noexec,relatime - cgroup2 cgroup rw,nsdelegate,memory_recursiveprot\n",
            mount_field(&unrelated)
        );
        let current_line = format!(
            "945 940 0:27 /../../jupyter-children {} ro,nosuid,nodev,noexec,relatime - cgroup2 cgroup2 rw,nsdelegate,memory_recursiveprot\n",
            mount_field(&current)
        );
        let membership = "0::/../../jupyter-children\n";
        assert!(!cgroup_memory_from_proc(membership, &unrelated_line).checked);
        for mountinfo in [
            format!("{unrelated_line}{current_line}"),
            format!("{current_line}{unrelated_line}"),
        ] {
            let observed = cgroup_memory_from_proc(membership, &mountinfo);
            assert!(observed.checked);
            assert_eq!(observed.limit, Some(100));
            assert_eq!(observed.headroom, Some(70));
        }

        fs::write(current.join("memory.max"), "invalid\n").unwrap();
        let mountinfo = format!("{unrelated_line}{current_line}");
        assert!(!cgroup_memory_from_proc(membership, &mountinfo).checked);
    }
}

#[cfg(target_os = "macos")]
#[path = "budget/macos.rs"]
mod platform;

#[cfg(not(any(target_os = "windows", target_os = "linux", target_os = "macos")))]
mod platform {
    use super::MemorySnapshot;

    pub(super) fn capture() -> MemorySnapshot {
        // Unsupported platforms do not substitute virtual size or RSS for a
        // meaningful private/footprint counter or claim available growth.
        MemorySnapshot::unknown()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reservation_budget() -> Budget {
        let mut policy = config();
        policy.sample_interval = Duration::from_secs(60);
        Budget::new(policy, sample(100, 10_000)).unwrap()
    }

    fn publish_cutoff(budget: &Budget, observation: MemorySnapshot, cutoff: SampleCutoff) {
        budget
            .state
            .lock()
            .unwrap()
            .publish_sample(observation, cutoff);
    }

    #[test]
    fn reservations_share_headroom_and_only_unused_cancellation_refunds_it() {
        let budget = reservation_budget();
        let first = budget.try_reserve(600, ReservationKind::Compute).unwrap();
        assert_eq!(budget.check(0).active_reserved_bytes, 600);
        assert!(matches!(
            budget.try_reserve(301, ReservationKind::Compute),
            Err(ReservationError::Pressure { .. })
        ));
        first.cancel_unused().unwrap();
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        let first = budget.try_reserve(600, ReservationKind::Compute).unwrap();
        first.start().unwrap();
        assert!(matches!(
            first.cancel_unused(),
            Err(ReservationError::AlreadyStarted)
        ));
        first.settle(400).unwrap();
        assert_eq!(budget.check(0).settled_reserved_bytes, 600);
        assert!(matches!(first.ensure(601), Err(ReservationError::Finished)));
        assert!(matches!(
            budget.try_reserve(301, ReservationKind::Compute),
            Err(ReservationError::Pressure { .. })
        ));
    }

    #[test]
    fn cutoff_absorbs_only_preexisting_settlements_not_active_or_later_work() {
        let budget = reservation_budget();
        let first = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        first.start().unwrap();
        first.settle(80).unwrap();
        let active = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        active.start().unwrap();
        let cutoff = budget.state.lock().unwrap().begin_sample();
        // This group finishes during the OS call; the sample may have observed
        // it before allocation. Retiring it now would spend its bytes twice.
        active.settle(90).unwrap();
        publish_cutoff(&budget, sample(270, 10_000), cutoff);
        let decision = budget.check(0);
        assert_eq!(decision.active_reserved_bytes, 0);
        assert_eq!(decision.settled_reserved_bytes, 100);
        assert_eq!(decision.projected_private_bytes, Some(370));
        assert_eq!(decision.sample_generation, cutoff.generation);
        let next = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(270, 10_000), next);
        assert_eq!(budget.check(0).settled_reserved_bytes, 0);
        assert_eq!(budget.check(0).projected_private_bytes, Some(270));
        // An obsolete publication cannot overwrite a newer accounting base.
        publish_cutoff(&budget, sample(0, 10_000), cutoff);
        assert_eq!(budget.check(0).known_private_bytes, Some(270));
    }

    #[test]
    fn active_envelopes_survive_sampling_and_started_drop_keeps_liability() {
        let budget = reservation_budget();
        let group = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        let other = group.clone();
        group.start().unwrap();
        let cutoff = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(150, 10_000), cutoff);
        assert_eq!(budget.check(0).active_reserved_bytes, 100);
        drop(group);
        assert_eq!(budget.check(0).active_reserved_bytes, 100);
        drop(other);
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        assert_eq!(budget.check(0).settled_reserved_bytes, 100);
        // An unattributed external sample has no safe capture-start cutoff.
        budget.publish_snapshot(sample(200, 10_000));
        assert_eq!(budget.check(0).settled_reserved_bytes, 100);
        let cutoff = budget.state.lock().unwrap().begin_sample();
        let mut failed = MemorySnapshot::unknown();
        failed.available_physical_bytes = Some(10_000);
        publish_cutoff(&budget, failed, cutoff);
        assert_eq!(budget.check(0).settled_reserved_bytes, 100);
        let cutoff = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(200, 10_000), cutoff);
        assert_eq!(budget.check(0).settled_reserved_bytes, 0);
    }

    #[test]
    fn twenty_four_workers_expand_one_wave_envelope_without_duplicate_charges() {
        let budget = reservation_budget();
        // The initial envelope includes all threads*4 queued/replied jobs.
        let group = budget.try_reserve(384, ReservationKind::Compute).unwrap();
        group.start().unwrap();
        std::thread::scope(|scope| {
            for worker in 0..24u64 {
                let group = group.clone();
                scope.spawn(move || {
                    for job in 0..4 {
                        group.ensure(384 + worker * 4 + job).unwrap();
                    }
                });
            }
        });
        assert_eq!(group.reserved_bytes(), 479);
        assert_eq!(budget.check(0).active_reserved_bytes, 479);
        group.settle(450).unwrap();
        assert_eq!(budget.check(0).settled_reserved_bytes, 479);
    }

    #[test]
    fn growth_failure_and_underestimation_do_not_silently_return_credit() {
        let budget = reservation_budget();
        let group = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        assert!(matches!(
            group.ensure(901),
            Err(ReservationError::Pressure { .. })
        ));
        assert_eq!(group.reserved_bytes(), 100);
        assert_eq!(budget.check(0).active_reserved_bytes, 100);
        assert!(matches!(
            group.settle(120),
            Err(ReservationError::Underestimated {
                reserved: 100,
                measured: 120
            })
        ));
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        assert_eq!(budget.check(0).settled_reserved_bytes, 120);
        assert!(matches!(
            budget.try_reserve(u64::MAX, ReservationKind::Compute),
            Err(ReservationError::Overflow)
        ));
    }

    #[test]
    fn warning_floor_and_reclaim_hysteresis_do_not_deny_bounded_compute() {
        let mut policy = config();
        policy.sample_interval = Duration::from_secs(60);
        let budget = Budget::new(policy, sample(800, 10_000)).unwrap();
        assert_eq!(budget.check(0).pressure, Pressure::Warning);
        let wave = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        // A floor at 80% of the computation cap cannot reach the 60% reclaim
        // target, but this wave and its full growth envelope still fit.
        wave.ensure(200).unwrap();
        assert_eq!(budget.check(0).projected_private_bytes, Some(1000));
        assert!(matches!(
            wave.ensure(201),
            Err(ReservationError::Pressure { .. })
        ));
        assert!(matches!(
            budget.try_reserve(1, ReservationKind::Compute),
            Err(ReservationError::Pressure { .. })
        ));
        assert_eq!(wave.reserved_bytes(), 200);
        wave.settle(200).unwrap();
        let cutoff = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(800, 10_000), cutoff);
        let next = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        next.cancel_unused().unwrap();

        // Hysteresis remains a reclamation signal below the warning line;
        // it does not turn the reclaim target into a second admission cap.
        budget.publish_snapshot(sample(700, 10_000));
        assert_eq!(budget.check(0).pressure, Pressure::Warning);
        budget
            .try_reserve(200, ReservationKind::Compute)
            .unwrap()
            .cancel_unused()
            .unwrap();
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        assert_eq!(budget.check(0).settled_reserved_bytes, 0);
    }

    #[test]
    fn stale_samples_refuse_new_grants_and_critical_progress_can_use_the_reserve() {
        let budget = reservation_budget();
        let group = budget.try_reserve(600, ReservationKind::Compute).unwrap();
        let progress = budget.try_reserve(350, ReservationKind::Progress).unwrap();
        assert_eq!(budget.check(0).projected_private_bytes, Some(1050));
        assert!(matches!(
            budget.try_reserve(50, ReservationKind::Progress),
            Err(ReservationError::Pressure { .. })
        ));
        progress.cancel_unused().unwrap();
        group.cancel_unused().unwrap();
        budget.publish_snapshot(sample(100, 10_000));
        budget.set_sample_age_for_test(Duration::from_secs(121));
        assert_eq!(budget.check(0).sample_age, Duration::from_secs(121));
        assert!(budget.check(0).stale);
        assert_ne!(budget.check(0).pressure, Pressure::Normal);
        assert!(matches!(
            budget.try_reserve(1, ReservationKind::Compute),
            Err(ReservationError::Stale(_))
        ));
        assert!(matches!(
            budget.try_reserve(1, ReservationKind::Progress),
            Err(ReservationError::Stale(_))
        ));
        let cutoff = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(100, 10_000), cutoff);
        assert!(!budget.check(0).stale);
        budget
            .try_reserve(1, ReservationKind::Compute)
            .unwrap()
            .cancel_unused()
            .unwrap();
    }

    #[test]
    fn progress_reservation_does_not_wait_for_a_busy_ledger() {
        let budget = reservation_budget();
        let _locked = budget.state.lock().unwrap();
        assert!(matches!(
            budget.try_reserve(1, ReservationKind::Progress),
            Err(ReservationError::Busy)
        ));
    }

    #[test]
    fn future_credit_can_grow_shrink_and_cancel_until_explicit_start() {
        let budget = reservation_budget();
        let future = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        future.ensure(400).unwrap();
        assert_eq!(budget.check(0).active_reserved_bytes, 400);
        budget.publish_snapshot(sample(100, 10_000));
        budget.set_sample_age_for_test(Duration::from_secs(121));
        assert!(budget.check(0).stale);
        // Shrinking proved-unused memory does not compete for fresh admission.
        future.shrink_unstarted(100).unwrap();
        assert_eq!(future.reserved_bytes(), 100);
        assert_eq!(budget.check(0).active_reserved_bytes, 100);
        assert!(matches!(
            future.shrink_unstarted(101),
            Err(ReservationError::WouldGrow { .. })
        ));
        future.cancel_unused().unwrap();
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        assert_eq!(budget.check(0).settled_reserved_bytes, 0);
        assert!(matches!(
            future.shrink_unstarted(0),
            Err(ReservationError::Finished)
        ));

        let cutoff = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(100, 10_000), cutoff);
        let allocated = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        allocated.ensure(200).unwrap();
        allocated.start().unwrap();
        assert!(matches!(
            allocated.shrink_unstarted(100),
            Err(ReservationError::AlreadyStarted)
        ));
        assert_eq!(allocated.reserved_bytes(), 200);
        drop(allocated);
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        assert_eq!(budget.check(0).settled_reserved_bytes, 200);
    }

    #[test]
    fn emergency_future_credit_survives_materialized_value_sampling() {
        let budget = reservation_budget();
        let values = budget.try_reserve(600, ReservationKind::Compute).unwrap();
        let cleanup = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        values.start().unwrap();
        values.settle(600).unwrap();
        let cutoff = budget.state.lock().unwrap().begin_sample();
        publish_cutoff(&budget, sample(700, 10_000), cutoff);
        let decision = budget.check(0);
        assert_eq!(decision.settled_reserved_bytes, 0);
        assert_eq!(decision.active_reserved_bytes, 100);
        assert_eq!(decision.projected_private_bytes, Some(800));
        // A later cleanup allocation uses its original credit, without trying
        // to obtain new space at an already failing mathematical boundary.
        cleanup.start().unwrap();
        cleanup.settle(100).unwrap();
        assert_eq!(budget.check(0).settled_reserved_bytes, 100);

        let unused = budget.try_reserve(100, ReservationKind::Compute).unwrap();
        unused.ensure(150).unwrap();
        drop(unused);
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
        assert_eq!(budget.check(0).settled_reserved_bytes, 100);
    }

    #[test]
    fn visible_job_limits_constrain_process_and_group_headroom() {
        let mut observation = sample(100, 10_000);
        observation.job_process_limit_bytes = Some(700);
        let budget = Budget::new(config(), observation).unwrap();
        assert_eq!(budget.check(0).limits.unwrap().hard_bytes, 700);
        observation.job_memory_limit_bytes = Some(5000);
        observation.job_headroom_bytes = Some(400);
        observation.job_headroom_peak_based = true;
        budget.publish_snapshot(observation);
        assert_eq!(budget.check(0).limits.unwrap().hard_bytes, 400);
        assert_eq!(budget.check(0).memory.job_memory_limit_bytes, Some(5000));
        assert!(budget.check(0).memory.job_headroom_peak_based);
    }

    fn config() -> BudgetConfig {
        BudgetConfig {
            max_bytes: Some(1100),
            warning_percent: 80,
            reclaim_percent: 60,
            system_headroom_bytes: 100,
            progress_reserve_bytes: 100,
            sample_interval: Duration::from_millis(100),
        }
    }

    fn sample(private: u64, available: u64) -> MemorySnapshot {
        MemorySnapshot {
            process_private_bytes: Some(private),
            process_private_kind: Some({
                #[cfg(not(any(target_os = "linux", target_os = "macos")))]
                {
                    PrivateMemoryKind::Committed
                }
                #[cfg(target_os = "linux")]
                {
                    PrivateMemoryKind::AnonymousResidentAndSwap
                }
                #[cfg(target_os = "macos")]
                {
                    PrivateMemoryKind::PhysicalFootprint
                }
            }),
            available_physical_bytes: Some(available),
            environment_limits_checked: true,
            ..MemorySnapshot::unknown()
        }
    }

    #[test]
    fn warning_and_hard_pressure_preserve_progress_and_reclaim_hysteresis() {
        let budget = Budget::new(config(), sample(700, 10_000)).unwrap();
        assert_eq!(budget.check(0).pressure, Pressure::Normal);
        assert_eq!(budget.check(100).pressure, Pressure::Warning);
        budget.publish_snapshot(sample(700, 10_000));
        assert_eq!(budget.check(0).pressure, Pressure::Warning);
        budget.publish_snapshot(sample(600, 10_000));
        assert_eq!(budget.check(0).pressure, Pressure::Normal);
        let hard = budget.check(400);
        assert_eq!(hard.pressure, Pressure::Hard);
        assert!(hard.progress_at_risk);
        assert_eq!(hard.limits.unwrap().computation_max_bytes, 1000);
        assert_eq!(hard.reclaim_bytes, 400);
    }

    #[test]
    fn physical_commit_and_cgroup_headroom_tighten_one_shared_pool() {
        let mut observation = sample(300, 10_000);
        observation.commit_limit_enforced = true;
        observation.commit_headroom_bytes = Some(900);
        observation.cgroup_headroom_bytes = Some(600);
        let budget = Budget::new(config(), observation).unwrap();
        let decision = budget.check(0);
        assert_eq!(decision.limits.unwrap().hard_bytes, 800);
        assert_eq!(decision.limits.unwrap().computation_max_bytes, 700);
        assert_eq!(budget.check(400).pressure, Pressure::Hard);
        let mut tighter = observation;
        tighter.commit_headroom_bytes = Some(150);
        budget.publish_snapshot(tighter);
        assert_eq!(budget.check(0).pressure, Pressure::Hard);
    }

    #[test]
    fn advisory_linux_commit_limit_does_not_override_physical_headroom() {
        let mut observation = sample(300, 10_000);
        observation.commit_headroom_bytes = Some(0);
        let budget = Budget::new(config(), observation).unwrap();
        assert_eq!(budget.check(0).pressure, Pressure::Normal);
        assert_eq!(budget.check(0).limits.unwrap().hard_bytes, 1100);
    }

    #[test]
    fn missing_metrics_keep_explicit_and_previous_limits_actionable() {
        let budget = Budget::new(config(), MemorySnapshot::unknown()).unwrap();
        assert_eq!(budget.check(0).pressure, Pressure::Unavailable);
        assert_eq!(budget.check(0).limits.unwrap().hard_bytes, 1100);
        assert_eq!(budget.check(1000).pressure, Pressure::Hard);
        budget.publish_snapshot(sample(900, 10_000));
        assert_eq!(budget.check(0).pressure, Pressure::Warning);
        budget.publish_snapshot(MemorySnapshot::unknown());
        let decision = budget.check(100);
        assert_eq!(decision.pressure, Pressure::Hard);
        assert_eq!(decision.known_private_bytes, Some(900));
        assert!(!decision.observations_complete);
    }

    #[test]
    fn invalid_configuration_and_projected_overflow_are_explicit() {
        let mut invalid = config();
        invalid.max_bytes = Some(99);
        assert!(Budget::new(invalid, sample(1, 10_000)).is_err());
        let budget = Budget::new(config(), sample(u64::MAX - 1, 0)).unwrap();
        let decision = budget.check(2);
        assert_eq!(decision.pressure, Pressure::Hard);
        assert!(decision.arithmetic_overflow);
        assert_eq!(decision.projected_private_bytes, None);
        assert_eq!(
            percent(u64::MAX, 99),
            u64::MAX / 100 * 99 + u64::MAX % 100 * 99 / 100
        );
    }

    #[test]
    fn missing_sample_keeps_an_automatically_derived_ceiling() {
        let mut policy = config();
        policy.max_bytes = None;
        let budget = Budget::new(policy, sample(100, 1000)).unwrap();
        let ceiling = budget.check(0).limits.unwrap();
        assert_eq!(ceiling.hard_bytes, 1000);
        budget.publish_snapshot(MemorySnapshot::unknown());
        let unavailable = budget.check(0);
        assert_eq!(unavailable.pressure, Pressure::Unavailable);
        assert_eq!(unavailable.limits, Some(ceiling));
        assert_eq!(budget.check(800).pressure, Pressure::Hard);
    }

    #[test]
    fn cached_checks_do_not_refresh_the_native_observation() {
        let observation = sample(100, 10_000);
        let budget = Budget::new(config(), observation).unwrap();
        for _ in 0..4 {
            assert_eq!(budget.check(0).memory.captured_at, observation.captured_at);
        }
    }

    #[test]
    fn process_policy_is_shared_and_native_counters_have_explicit_availability() {
        assert!(std::ptr::eq(process_budget(), process_budget()));
        let observation = MemorySnapshot::capture();
        #[cfg(any(target_os = "windows", target_os = "linux", target_os = "macos"))]
        {
            assert!(observation.process_private_bytes.is_some());
            assert!(observation.available_physical_bytes.is_some());
            assert!(observation.process_resident_bytes.is_some());
        }
        #[cfg(not(any(target_os = "windows", target_os = "linux", target_os = "macos")))]
        assert!(observation.process_private_bytes.is_none());
    }
}
