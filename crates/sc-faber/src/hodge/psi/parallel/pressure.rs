//! Allocation boundaries and the bounded, disk-backed Plan slow path.
//!
//! The ordinary arithmetic kernels never enter this module. A planning wave
//! owns one credit shared by its jobs/replies; only Vec capacity changes touch
//! that credit. Native sampling and workspace waits happen after the borrowed
//! mathematical tables have been released.

use super::workspace::PlanWorkspace;
use super::{Batch, PackedEquation, ParallelSource, Plan, State, TermBlock};
use crate::cache_source::budget::{
    Budget, BudgetDecision, Pressure, Reservation, ReservationError, ReservationKind,
    process_budget,
};
use crate::cache_source::{Registry, resource_failure};
use sc_num::{Q, Zero};
use std::cell::RefCell;
use std::mem::size_of;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

const WAVE_STATES: usize = 4096;
const GROUP_EQUATIONS: usize = 1024;
const DECODE_WINDOW: u64 = 1024 * 1024;
const PROGRESS_WAIT: Duration = Duration::from_millis(100);
const PROGRESS_TIMEOUT: Duration = Duration::from_secs(30);

fn bytes(count: usize, width: usize) -> u64 {
    count
        .checked_mul(width)
        .and_then(|value| u64::try_from(value).ok())
        .unwrap_or_else(|| resource_failure("plan allocation size overflow"))
}

fn failure(context: &str, error: impl std::fmt::Display) -> ! {
    resource_failure(format!("{context}: {error}"))
}

fn retry_initial_reservation<T>(
    mut attempt: impl FnMut() -> Result<T, ReservationError>,
    mut refresh: impl FnMut(),
    mut pause: impl FnMut(),
    mut elapsed: impl FnMut() -> Duration,
    timeout: Duration,
) -> Result<T, ReservationError> {
    let mut error = match attempt() {
        Ok(credit) => return Ok(credit),
        Err(error @ (ReservationError::Busy | ReservationError::Stale(_))) => error,
        Err(error) => return Err(error),
    };
    // Only transient failures start the clock. An uncontended fresh admission
    // performs its original one attempt, without a timer, sample or wait.
    let _ = elapsed();
    loop {
        if matches!(error, ReservationError::Stale(_)) {
            refresh();
        }
        pause();
        if elapsed() >= timeout {
            return Err(error);
        }
        match attempt() {
            Ok(credit) => return Ok(credit),
            Err(next @ (ReservationError::Busy | ReservationError::Stale(_))) => error = next,
            Err(error) => return Err(error),
        }
    }
}

/// Coordinator-only admission. Both call sites are outside source/table guards
/// and worker scopes; storage's own progress path keeps its nonblocking policy.
fn initial_reservation(
    budget: &Budget,
    bytes: u64,
    kind: ReservationKind,
    context: &str,
) -> Reservation {
    let mut started: Option<Instant> = None;
    retry_initial_reservation(
        || budget.try_reserve(bytes, kind),
        || {
            budget.force_check(0);
        },
        || std::thread::park_timeout(Duration::from_millis(1)),
        || started.get_or_insert_with(Instant::now).elapsed(),
        PROGRESS_TIMEOUT,
    )
    .unwrap_or_else(|error| match error {
        ReservationError::Busy | ReservationError::Stale(_) => failure(
            context,
            format!("transient admission did not recover within 30 seconds: {error}"),
        ),
        error => failure(context, error),
    })
}

// One group shares the refresh across shrinking candidates. Stale observations
// may still request a newer sample; a genuine fresh denial alone can shrink it.
fn evaluation_admission<T>(
    refreshed: &mut bool,
    mut attempt: impl FnMut() -> Result<T, ReservationError>,
    mut refresh: impl FnMut(),
    mut pause: impl FnMut(),
    mut elapsed: impl FnMut() -> Duration,
) -> Result<T, ReservationError> {
    let result = retry_initial_reservation(
        &mut attempt,
        || {
            *refreshed = true;
            refresh();
        },
        &mut pause,
        &mut elapsed,
        PROGRESS_TIMEOUT,
    );
    if !*refreshed
        && matches!(
            &result,
            Err(ReservationError::Pressure { .. } | ReservationError::Unavailable)
        )
    {
        *refreshed = true;
        refresh();
        return retry_initial_reservation(attempt, refresh, pause, elapsed, PROGRESS_TIMEOUT);
    }
    result
}

fn evaluation_reserve<T>(
    budget: &Budget,
    refreshed: &mut bool,
    attempt: impl FnMut() -> Result<T, ReservationError>,
) -> Result<T, ReservationError> {
    // A group may have spent longer than 30 seconds waiting for real spill
    // progress. Only a new transient failure starts this independent clock.
    let mut transient_since: Option<Instant> = None;
    evaluation_admission(
        refreshed,
        attempt,
        || {
            budget.force_check(0);
        },
        || std::thread::park_timeout(Duration::from_millis(1)),
        || transient_since.get_or_insert_with(Instant::now).elapsed(),
    )
}

fn capacity_denial(error: &ReservationError) -> bool {
    matches!(error, ReservationError::Pressure { .. })
}

// Match Budget's Compute approval, including the inclusive computation limit.
// This observation grants no credit: the caller must retry its atomic reserve.
fn compute_fits(decision: &BudgetDecision) -> bool {
    !decision.stale
        && !decision.arithmetic_overflow
        && decision.observations_complete
        && decision
            .projected_private_bytes
            .zip(decision.limits)
            .is_some_and(|(projected, limits)| projected <= limits.computation_max_bytes)
}

#[cfg(test)]
mod evaluation_admission_tests {
    use super::*;
    use crate::cache_source::budget::{BudgetConfig, MemorySnapshot};
    use std::cell::Cell;
    use std::collections::VecDeque;

    fn pressure() -> ReservationError {
        ReservationError::Pressure {
            pressure: Pressure::Hard,
            projected: Some(1024),
            limits: None,
        }
    }

    fn budget() -> Budget {
        Budget::for_test(
            BudgetConfig {
                max_bytes: Some(1000),
                system_headroom_bytes: 0,
                progress_reserve_bytes: 32,
                sample_interval: Duration::from_secs(60),
                ..BudgetConfig::default()
            },
            MemorySnapshot {
                captured_at: Instant::now(),
                process_private_bytes: Some(64),
                process_peak_commit_bytes: None,
                process_private_kind: None,
                process_resident_bytes: None,
                total_physical_bytes: Some(20_000),
                available_physical_bytes: Some(10_000),
                commit_headroom_bytes: None,
                commit_limit_enforced: false,
                cgroup_headroom_bytes: None,
                cgroup_limit_bytes: None,
                cgroup_observation: None,
                job_process_limit_bytes: None,
                job_memory_limit_bytes: None,
                job_headroom_bytes: None,
                job_headroom_peak_based: false,
                environment_limits_checked: true,
            },
        )
        .unwrap()
    }

    #[test]
    fn a_fresh_fitting_request_can_skip_reclaim_even_at_the_compute_limit() {
        let budget = budget();
        let mut decision = budget.check(904);
        // Pressure is Hard at equality, but Compute's actual atomic approval
        // accepts equality. Warning/hysteresis alone must not force a wait.
        assert_eq!(decision.pressure, Pressure::Hard);
        assert!(compute_fits(&decision));
        let credit = budget.try_reserve(904, ReservationKind::Compute).unwrap();
        drop(credit);
        assert!(!compute_fits(&budget.check(905)));
        decision.stale = true;
        assert!(!compute_fits(&decision));
        decision.stale = false;
        decision.observations_complete = false;
        assert!(!compute_fits(&decision));
        decision.observations_complete = true;
        decision.arithmetic_overflow = true;
        assert!(!compute_fits(&decision));
    }

    #[test]
    fn rejected_future_outputs_leave_only_allocated_metadata_liability() {
        let budget = budget();
        let metadata = budget.try_reserve(200, ReservationKind::Compute).unwrap();
        metadata.start().unwrap();
        assert!(capacity_denial(&metadata.ensure(1200).unwrap_err()));
        assert_eq!(metadata.reserved_bytes(), 200);
        drop(metadata);
        let decision = budget.check(0);
        assert_eq!(decision.active_reserved_bytes, 0);
        assert_eq!(decision.settled_reserved_bytes, 200);
    }

    #[test]
    fn fresh_success_and_fatal_errors_do_not_refresh_or_wait() {
        assert_eq!(
            evaluation_admission(
                &mut false,
                || Ok(17usize),
                || panic!("fresh success sampled memory"),
                || panic!("fresh success waited"),
                || panic!("fresh success started a clock"),
            )
            .unwrap(),
            17
        );
        for error in [
            ReservationError::Overflow,
            ReservationError::Finished,
            ReservationError::AlreadyStarted,
            ReservationError::WouldGrow {
                reserved: 1,
                required: 2,
            },
            ReservationError::Underestimated {
                reserved: 1,
                measured: 2,
            },
        ] {
            let mut error = Some(error);
            let result = evaluation_admission::<usize>(
                &mut false,
                || Err(error.take().expect("fatal error was retried")),
                || panic!("fatal error sampled memory"),
                || panic!("fatal error waited"),
                || panic!("fatal error started a clock"),
            );
            assert!(!capacity_denial(&result.unwrap_err()));
        }
    }

    #[test]
    fn a_fresh_observation_retries_the_original_request_before_shrinking() {
        let mut attempts = VecDeque::from([Err(pressure()), Ok(23usize)]);
        let refreshes = Cell::new(0);
        let mut refreshed = false;
        assert_eq!(
            evaluation_admission(
                &mut refreshed,
                || attempts.pop_front().expect("extra attempt"),
                || refreshes.set(refreshes.get() + 1),
                || panic!("successful fresh retry waited"),
                || panic!("successful fresh retry started a clock"),
            )
            .unwrap(),
            23
        );
        assert!(attempts.is_empty());
        assert!(refreshed);
        assert_eq!(refreshes.get(), 1);
    }

    #[test]
    fn genuine_pressure_shares_one_refresh_across_smaller_candidates() {
        let refreshes = Cell::new(0);
        let mut refreshed = false;
        for _ in 0..4 {
            assert!(capacity_denial(
                &evaluation_admission::<usize>(
                    &mut refreshed,
                    || Err(pressure()),
                    || refreshes.set(refreshes.get() + 1),
                    || panic!("pressure classification waited"),
                    || panic!("pressure classification started a clock"),
                )
                .unwrap_err()
            ));
        }
        assert_eq!(refreshes.get(), 1);
    }

    #[test]
    fn missing_observations_are_refreshed_once_and_never_shrink() {
        let refreshes = Cell::new(0);
        let result = evaluation_admission::<usize>(
            &mut false,
            || Err(ReservationError::Unavailable),
            || refreshes.set(refreshes.get() + 1),
            || panic!("missing observations waited for block reclamation"),
            || panic!("missing observations started a clock"),
        );
        assert!(matches!(result, Err(ReservationError::Unavailable)));
        assert_eq!(refreshes.get(), 1);
    }

    #[test]
    fn busy_and_stale_retries_keep_the_original_request_and_are_bounded() {
        let mut attempts = VecDeque::from([
            Err(ReservationError::Busy),
            Err(ReservationError::Stale(Duration::from_secs(1))),
            Ok(29usize),
        ]);
        let refreshes = Cell::new(0);
        let pauses = Cell::new(0);
        let mut refreshed = false;
        assert_eq!(
            evaluation_admission(
                &mut refreshed,
                || attempts.pop_front().expect("extra transient attempt"),
                || refreshes.set(refreshes.get() + 1),
                || pauses.set(pauses.get() + 1),
                || Duration::ZERO,
            )
            .unwrap(),
            29
        );
        assert!(attempts.is_empty());
        assert!(refreshed);
        assert_eq!(refreshes.get(), 1);
        assert_eq!(pauses.get(), 2);

        let elapsed = Cell::new(0);
        let result = evaluation_admission::<usize>(
            &mut false,
            || Err(ReservationError::Busy),
            || panic!("busy ledger needs no native sample"),
            || {},
            || {
                let seconds = elapsed.get();
                elapsed.set(seconds + 15);
                Duration::from_secs(seconds)
            },
        );
        assert!(matches!(result, Err(ReservationError::Busy)));
        assert_eq!(elapsed.get(), 45);
    }
}

#[cfg(test)]
mod initial_admission_tests {
    use super::*;
    use std::cell::Cell;
    use std::collections::VecDeque;

    #[test]
    fn fresh_success_and_real_denials_never_time_sample_or_wait() {
        assert_eq!(
            retry_initial_reservation(
                || Ok(17usize),
                || panic!("unexpected refresh"),
                || panic!("unexpected wait"),
                || panic!("unexpected clock"),
                PROGRESS_TIMEOUT,
            )
            .unwrap(),
            17
        );
        for error in [
            ReservationError::Overflow,
            ReservationError::Unavailable,
            ReservationError::Pressure {
                pressure: Pressure::Hard,
                projected: Some(10),
                limits: None,
            },
        ] {
            let mut error = Some(error);
            let result = retry_initial_reservation::<usize>(
                || Err(error.take().expect("real denial retried")),
                || panic!("unexpected refresh"),
                || panic!("unexpected wait"),
                || panic!("unexpected clock"),
                PROGRESS_TIMEOUT,
            );
            assert!(result.is_err());
        }
    }

    #[test]
    fn busy_and_still_stale_sampler_retries_preserve_the_real_admission() {
        let mut attempts = VecDeque::from([
            Err(ReservationError::Busy),
            Err(ReservationError::Stale(Duration::from_secs(1))),
            Err(ReservationError::Stale(Duration::from_secs(1))),
            Ok(23usize),
        ]);
        let refreshes = Cell::new(0);
        let pauses = Cell::new(0);
        let result = retry_initial_reservation(
            || attempts.pop_front().expect("extra admission attempt"),
            || refreshes.set(refreshes.get() + 1),
            || pauses.set(pauses.get() + 1),
            || Duration::ZERO,
            PROGRESS_TIMEOUT,
        )
        .unwrap();
        assert_eq!(result, 23);
        assert!(attempts.is_empty());
        assert_eq!(refreshes.get(), 2);
        assert_eq!(pauses.get(), 3);
    }

    #[test]
    fn repeated_contention_stops_at_a_bounded_deadline_without_real_sleep() {
        let elapsed = Cell::new(0);
        let attempts = Cell::new(0);
        let result = retry_initial_reservation::<usize>(
            || {
                attempts.set(attempts.get() + 1);
                Err(ReservationError::Busy)
            },
            || panic!("busy ledger alone needs no native resample"),
            || {},
            || {
                let now = elapsed.get();
                elapsed.set(now + 1);
                Duration::from_millis(now)
            },
            Duration::from_millis(3),
        );
        assert!(matches!(result, Err(ReservationError::Busy)));
        assert_eq!(attempts.get(), 3);
    }
}

struct FutureValues {
    credit: Reservation,
    cleanup: Reservation,
    slots: usize,
    base: usize,
    old_capacity: usize,
    materialized: bool,
    cleanup_started: bool,
    finished: bool,
}

impl FutureValues {
    fn value_credit(&self, slots: usize) -> u64 {
        let end = self
            .base
            .checked_add(slots)
            .unwrap_or_else(|| resource_failure("future result count overflow"));
        if end > self.old_capacity {
            bytes(end, size_of::<Q>())
        } else {
            0
        }
    }
}

/// Kept outside integrate's catch_unwind, so cleanup has credit for every
/// assigned ID even if planning fails before the value vector is allocated.
pub(super) struct Control {
    future: Mutex<FutureValues>,
    state_bytes: usize,
}

impl Control {
    pub(super) fn new(base: usize, old_capacity: usize, max_key_bytes: usize) -> Self {
        let budget = process_budget();
        let reserve = |context| initial_reservation(budget, 0, ReservationKind::Compute, context);
        Self {
            future: Mutex::new(FutureValues {
                credit: reserve("future result slots"),
                cleanup: reserve("future cleanup bitmap"),
                slots: 0,
                base,
                old_capacity,
                materialized: false,
                cleanup_started: false,
                finished: false,
            }),
            // The original BTree keys and frontier state both own their
            // vectors. Node/allocator overhead remains a conservative estimate
            // and is subsequently covered by process-private sampling.
            state_bytes: size_of::<State>()
                .saturating_add(max_key_bytes.saturating_mul(2))
                .saturating_add(192),
        }
    }

    pub(super) fn state_bytes(&self) -> usize {
        self.state_bytes
    }

    /// Called only with no wave jobs in flight. Discovery may reserve spare
    /// local capacity, so a checkpoint records credits rather than assigned IDs.
    pub(super) fn checkpoint(&self) -> usize {
        self.future
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .slots
    }

    /// Undo only future allocations after all failed-wave jobs have joined.
    /// The surviving IDs keep their original credits throughout this operation.
    pub(super) fn rollback(&self, saved: usize) -> Result<(), String> {
        let mut future = self
            .future
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if future.materialized || future.cleanup_started || future.finished {
            return Err("cannot roll back materialized result or cleanup storage".into());
        }
        if saved > future.slots {
            return Err("invalid future-slot rollback checkpoint".into());
        }
        future
            .credit
            .shrink_unstarted(future.value_credit(saved))
            .map_err(|error| error.to_string())?;
        future
            .cleanup
            .shrink_unstarted(bytes(saved, size_of::<bool>()))
            .map_err(|error| error.to_string())?;
        future.slots = saved;
        Ok(())
    }

    /// Called before growing a worker's discovery Vec, outside its index lock.
    /// Unused local capacity is intentionally conservative; no per-node atomic
    /// or budget lookup is needed after this chunk has been granted.
    pub(super) fn claim_slots(&self, additional: usize) {
        let mut future = self
            .future
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if future.materialized || future.finished {
            resource_failure("cannot discover IDs after result allocation");
        }
        let slots = future
            .slots
            .checked_add(additional)
            .unwrap_or_else(|| resource_failure("future result count overflow"));
        // A realloc may hold the old allocation and the complete new buffer at
        // once. The old buffer is already in the native observation; reserving
        // only its incremental length would understate that temporary peak.
        let previous = future.credit.reserved_bytes();
        future
            .credit
            .ensure(future.value_credit(slots))
            .unwrap_or_else(|error| exhaustion(error, "future result slots"));
        // Keep this future bitmap active even after the value vector is
        // materialized and its own credit has been absorbed by native sampling.
        if let Err(error) = future.cleanup.ensure(bytes(slots, size_of::<bool>())) {
            future
                .credit
                .shrink_unstarted(previous)
                .unwrap_or_else(|error| failure("undo future result credit", error));
            exhaustion(error, "future cleanup bitmap");
        }
        future.slots = slots;
    }

    /// No new admission is needed when failure cleanup actually allocates its
    /// bitmap. Its separate active credit has survived every native cutoff.
    pub(super) fn start_cleanup(&self) {
        let mut future = self
            .future
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if future.cleanup_started {
            return;
        }
        future
            .cleanup
            .start()
            .unwrap_or_else(|error| failure("allocate cleanup bitmap", error));
        future.cleanup_started = true;
    }

    /// A successful query never allocates the emergency bitmap. Return that
    /// unused credit immediately; used storage retains normal sample liability.
    pub(super) fn finish(&self) {
        let mut future = self
            .future
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if future.finished {
            return;
        }
        if !future.materialized {
            future
                .credit
                .cancel_unused()
                .unwrap_or_else(|error| failure("finish future result credit", error));
        }
        if future.cleanup_started {
            future
                .cleanup
                .settle(future.cleanup.reserved_bytes())
                .unwrap_or_else(|error| failure("finish cleanup bitmap credit", error));
        } else {
            future
                .cleanup
                .cancel_unused()
                .unwrap_or_else(|error| failure("finish unused cleanup credit", error));
        }
        future.finished = true;
    }

    pub(super) fn resize_values(&self, source: &ParallelSource, end: usize) {
        let mut future = self
            .future
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if future.materialized {
            // Evaluation failures reach the same cleanup path after the slots
            // have already been allocated. It must not acquire a second credit.
            debug_assert!(source.read().len() >= end);
            return;
        }
        debug_assert!(end <= future.base.saturating_add(future.slots));
        future
            .credit
            .start()
            .unwrap_or_else(|error| failure("allocate result slots", error));
        let allocated = {
            let mut values = source.write();
            let old_capacity = values.capacity();
            if end > old_capacity {
                let additional = end.saturating_sub(values.len());
                values.reserve_exact(additional);
            }
            values.resize(end, Q::ZERO);
            if values.capacity() > old_capacity {
                bytes(values.capacity(), size_of::<Q>())
            } else {
                0
            }
        };
        // Publish completion before any diagnostic error can unwind back into
        // cleanup. All assigned slots now exist and remain mathematically owned.
        future.materialized = true;
        future
            .credit
            .settle(allocated)
            .unwrap_or_else(|error| failure("allocate result slots", error));
    }
}

pub(super) struct WaveExhausted {
    pub(super) message: String,
    pub(super) stale: bool,
}

pub(super) fn exhaustion(error: ReservationError, context: &str) -> ! {
    let stale = matches!(&error, ReservationError::Stale(_));
    match error {
        ReservationError::Pressure { .. } | ReservationError::Stale(_) | ReservationError::Busy => {
            std::panic::resume_unwind(Box::new(WaveExhausted {
                message: format!("{context}: {error}"),
                stale,
            }))
        }
        error => failure(context, error),
    }
}

pub(super) fn propagate(payload: Box<dyn std::any::Any + Send>) -> ! {
    match payload.downcast::<WaveExhausted>() {
        Ok(error) => resource_failure(error.message),
        Err(payload) => std::panic::resume_unwind(payload),
    }
}

struct WaveInner {
    credit: Reservation,
    allocated: Mutex<u64>,
    backed: bool,
    #[cfg(test)]
    fail_after: std::sync::atomic::AtomicUsize,
    #[cfg(test)]
    stale_injection: Option<Arc<StaleInjection>>,
    #[cfg(test)]
    normal_boundary: bool,
}

#[cfg(test)]
struct StaleInjection {
    at_capacity: bool,
    remaining: std::sync::atomic::AtomicUsize,
    armed: std::sync::atomic::AtomicBool,
    injected: std::sync::atomic::AtomicUsize,
}

#[derive(Clone)]
pub(super) struct Wave(Arc<WaveInner>);

thread_local! {
    static WAVE: RefCell<Option<Wave>> = const { RefCell::new(None) };
}

#[cfg(test)]
thread_local! {
    static FORCE_BACKING: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static BACKED_PLANS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static FAIL_NEXT_WAVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static STALE_WAVES: RefCell<Option<Arc<StaleInjection>>> = const { RefCell::new(None) };
}

#[cfg(test)]
pub(super) fn stale_across_waves<R>(
    at_capacity: bool,
    run: impl FnOnce() -> R,
) -> (R, usize, usize) {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    struct Restore(Option<Arc<StaleInjection>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            STALE_WAVES.with(|cell| {
                *cell.borrow_mut() = self.0.take();
            });
        }
    }
    let injection = Arc::new(StaleInjection {
        at_capacity,
        remaining: AtomicUsize::new(2),
        armed: AtomicBool::new(true),
        injected: AtomicUsize::new(0),
    });
    let _restore = Restore(STALE_WAVES.with(|cell| cell.replace(Some(injection.clone()))));
    let before = BACKED_PLANS.with(|cell| cell.get());
    let result = run();
    (
        result,
        injection.injected.load(Ordering::Relaxed),
        BACKED_PLANS.with(|cell| cell.get()) - before,
    )
}

#[cfg(test)]
fn next_stale() -> Option<Arc<StaleInjection>> {
    use std::sync::atomic::Ordering;
    STALE_WAVES.with(|cell| {
        let binding = cell.borrow();
        let injection = binding.as_ref()?;
        if injection.remaining.load(Ordering::Relaxed) == 0
            || !injection.armed.swap(false, Ordering::Relaxed)
        {
            return None;
        }
        injection.remaining.fetch_sub(1, Ordering::Relaxed);
        Some(injection.clone())
    })
}

#[cfg(test)]
pub(super) fn force_backing<R>(run: impl FnOnce() -> R) -> (R, usize) {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            FORCE_BACKING.with(|cell| cell.set(self.0));
        }
    }
    let _restore = Restore(FORCE_BACKING.with(|cell| cell.replace(true)));
    let before = BACKED_PLANS.with(|cell| cell.get());
    let result = run();
    (result, BACKED_PLANS.with(|cell| cell.get()) - before)
}

#[cfg(test)]
pub(super) fn exhaust_next_wave<R>(run: impl FnOnce() -> R) -> (R, usize) {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            FAIL_NEXT_WAVE.with(|cell| cell.set(self.0));
        }
    }
    let _restore = Restore(FAIL_NEXT_WAVE.with(|cell| cell.replace(true)));
    let before = BACKED_PLANS.with(|cell| cell.get());
    let result = run();
    (result, BACKED_PLANS.with(|cell| cell.get()) - before)
}

impl Wave {
    pub(super) fn begin(count: usize, jobs: usize, pressure_view: bool) -> Result<Self, String> {
        Self::begin_with(count, jobs, pressure_view, pressure_view)
    }

    pub(super) fn begin_root(count: usize, jobs: usize) -> Result<Self, String> {
        Self::begin_with(count, jobs, true, false)
    }

    fn begin_with(
        count: usize,
        jobs: usize,
        pressure_view: bool,
        backed: bool,
    ) -> Result<Self, String> {
        #[cfg(test)]
        if !pressure_view && FORCE_BACKING.with(|cell| cell.get()) {
            return Err("test pressure boundary".into());
        }
        #[cfg(test)]
        let stale_injection = if pressure_view { None } else { next_stale() };
        #[cfg(test)]
        if let Some(injection) = stale_injection
            .as_ref()
            .filter(|injection| !injection.at_capacity)
        {
            injection
                .injected
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return Err("test stale plan observation".into());
        }
        let budget = process_budget();
        let estimate = bytes(
            count,
            size_of::<PackedEquation>() + 12 * size_of::<usize>() + 128,
        )
        .saturating_add(bytes(jobs.max(1), 64 * 1024));
        let decision = budget.check(estimate);
        // Warning is a request to leave the borrowed fast path. Once blocks
        // are backed, admission is determined by the computation cap, rather
        // than requiring an irreducible floor to fall below the warning line.
        if !pressure_view && (decision.pressure != Pressure::Normal || decision.stale) {
            return Err(format!(
                "plan boundary: {:?}, stale={}",
                decision.pressure, decision.stale
            ));
        }
        let credit = budget
            .try_reserve(estimate, ReservationKind::Compute)
            .map_err(|error| error.to_string())?;
        credit.start().map_err(|error| error.to_string())?;
        Ok(Self(Arc::new(WaveInner {
            credit,
            allocated: Mutex::new(0),
            backed,
            #[cfg(test)]
            fail_after: std::sync::atomic::AtomicUsize::new(
                // Equation-vector reserve, discovery-vector reserve, then the
                // first term reserve: the child ID exists before this third
                // capacity boundary rejects the wave.
                if stale_injection.is_some()
                    || (!pressure_view && FAIL_NEXT_WAVE.with(|cell| cell.replace(false)))
                {
                    3
                } else {
                    0
                },
            ),
            #[cfg(test)]
            stale_injection,
            #[cfg(test)]
            normal_boundary: !pressure_view,
        })))
    }

    #[inline]
    fn check_injection(&self) {
        #[cfg(test)]
        {
            use std::sync::atomic::Ordering;
            let mut remaining = self.0.fail_after.load(Ordering::Relaxed);
            while remaining != 0 {
                match self.0.fail_after.compare_exchange_weak(
                    remaining,
                    remaining - 1,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => {
                        if remaining == 1 {
                            let stale = self.0.stale_injection.is_some();
                            if let Some(injection) = &self.0.stale_injection {
                                injection.injected.fetch_add(1, Ordering::Relaxed);
                            }
                            std::panic::resume_unwind(Box::new(WaveExhausted {
                                message: if stale {
                                    "test wave stale rejection"
                                } else {
                                    "test wave capacity rejection"
                                }
                                .into(),
                                stale,
                            }));
                        }
                        break;
                    }
                    Err(current) => remaining = current,
                }
            }
        }
    }

    /// A refused candidate must not consume any of this shared wave's usage
    /// or reservation. Other jobs can still use the previously admitted credit.
    fn try_allocate(&self, additional: u64) -> Result<(), ReservationError> {
        let mut allocated = self
            .0
            .allocated
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        let total = allocated
            .checked_add(additional)
            .ok_or(ReservationError::Overflow)?;
        self.0.credit.ensure(total)?;
        *allocated = total;
        Ok(())
    }

    fn allocate(&self, additional: u64) {
        self.check_injection();
        self.try_allocate(additional)
            .unwrap_or_else(|error| exhaustion(error, "plan wave capacity"));
    }

    fn capacity_target(&self, capacity: usize, needed: usize, width: usize) -> usize {
        // A forced rejection belongs to this logical growth operation, not to
        // each candidate. It must propagate before any adaptive admission.
        self.check_injection();
        let preferred = needed.max(capacity.saturating_mul(2)).max(4);
        let mut error = match self.try_allocate(bytes(preferred, width)) {
            Ok(()) => return preferred,
            Err(error) => error,
        };
        if self.0.backed && matches!(&error, ReservationError::Pressure { .. }) {
            let reduced = needed.max(capacity.saturating_add(capacity / 8));
            let mut previous = preferred;
            for target in [reduced, needed] {
                if target >= previous {
                    continue;
                }
                previous = target;
                // Reserve the complete replacement buffer each time. The old
                // buffer remains live and is already covered by native/credit.
                match self.try_allocate(bytes(target, width)) {
                    Ok(()) => return target,
                    Err(next) => error = next,
                }
                if !matches!(&error, ReservationError::Pressure { .. }) {
                    break;
                }
            }
        }
        exhaustion(error, "plan wave capacity")
    }

    pub(super) fn with<R>(wave: Option<&Self>, run: impl FnOnce() -> R) -> R {
        let Some(wave) = wave else {
            return run();
        };
        struct Restore(Option<Wave>);
        impl Drop for Restore {
            fn drop(&mut self) {
                WAVE.with(|cell| {
                    *cell.borrow_mut() = self.0.take();
                });
            }
        }
        let previous = WAVE.with(|cell| cell.replace(Some(wave.clone())));
        let _restore = Restore(previous);
        run()
    }

    pub(super) fn finish(self) {
        let allocated = *self
            .0
            .allocated
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        self.0
            .credit
            .settle(allocated)
            .unwrap_or_else(|error| failure("complete plan wave", error));
        #[cfg(test)]
        if self.0.normal_boundary {
            use std::sync::atomic::Ordering;
            // A rejected boundary stays unarmed while it is retried. Only
            // genuine committed-wave progress permits the next stale event.
            STALE_WAVES.with(|cell| {
                if let Some(injection) = cell.borrow().as_ref() {
                    if injection.remaining.load(Ordering::Relaxed) != 0 {
                        injection.armed.store(true, Ordering::Relaxed);
                    }
                }
            });
        }
    }
}

/// Call only on Vec capacity slow branches. No item-loop TLS/budget lookup.
pub(super) fn capacity<T>(values: &mut Vec<T>, needed: usize) {
    if needed <= values.capacity() {
        return;
    }
    WAVE.with(|cell| {
        let borrowed = cell.borrow();
        let Some(wave) = borrowed.as_ref() else {
            return;
        };
        let target = wave.capacity_target(values.capacity(), needed, size_of::<T>());
        values.reserve_exact(target - values.len());
        if values.capacity() > target {
            wave.allocate(bytes(values.capacity() - target, size_of::<T>()));
        }
    });
}

#[cfg(test)]
mod adaptive_capacity_tests {
    use super::*;
    use crate::cache_source::budget::{Budget, BudgetConfig, MemorySnapshot};
    use std::sync::atomic::AtomicUsize;

    fn local_wave(maximum: u64, backed: bool) -> (Budget, Wave, MemorySnapshot) {
        // The old 64-byte allocation is part of the synthetic native floor.
        // This policy object never reads or configures the process budget.
        let snapshot = MemorySnapshot {
            captured_at: Instant::now(),
            process_private_bytes: Some(64),
            process_peak_commit_bytes: None,
            process_private_kind: None,
            process_resident_bytes: None,
            total_physical_bytes: Some(20_000),
            available_physical_bytes: Some(10_000),
            commit_headroom_bytes: None,
            commit_limit_enforced: false,
            cgroup_headroom_bytes: None,
            cgroup_limit_bytes: None,
            cgroup_observation: None,
            job_process_limit_bytes: None,
            job_memory_limit_bytes: None,
            job_headroom_bytes: None,
            job_headroom_peak_based: false,
            environment_limits_checked: true,
        };
        let budget = Budget::for_test(
            BudgetConfig {
                max_bytes: Some(maximum),
                warning_percent: 85,
                reclaim_percent: 70,
                system_headroom_bytes: 0,
                progress_reserve_bytes: 32,
                sample_interval: Duration::from_secs(60),
            },
            snapshot,
        )
        .unwrap();
        let credit = budget.try_reserve(0, ReservationKind::Compute).unwrap();
        credit.start().unwrap();
        let wave = Wave(Arc::new(WaveInner {
            credit,
            allocated: Mutex::new(0),
            backed,
            fail_after: AtomicUsize::new(0),
            stale_injection: None,
            normal_boundary: !backed,
        }));
        (budget, wave, snapshot)
    }

    fn data() -> Vec<u8> {
        let mut values = Vec::with_capacity(64);
        values.extend(0u8..64);
        assert_eq!(values.capacity(), 64);
        values
    }

    fn unchanged_credit(budget: &Budget, wave: &Wave) {
        assert_eq!(*wave.0.allocated.lock().unwrap(), 0);
        assert_eq!(wave.0.credit.reserved_bytes(), 0);
        assert_eq!(budget.check(0).active_reserved_bytes, 0);
    }

    #[test]
    fn backed_growth_admits_full_smaller_buffer_after_doubling_is_refused() {
        // First case fits 1/8 growth; the second fits only the exact minimum.
        for (maximum, target) in [(192, 72), (164, 65)] {
            let (budget, wave, _) = local_wave(maximum, true);
            let mut values = data();
            assert!(matches!(
                wave.try_allocate(128),
                Err(ReservationError::Pressure { .. })
            ));
            unchanged_credit(&budget, &wave);
            Wave::with(Some(&wave), || capacity(&mut values, 65));
            assert_eq!(values.capacity(), target);
            assert_eq!(values, data());
            assert_eq!(*wave.0.allocated.lock().unwrap(), target as u64);
            let decision = budget.check(0);
            assert_eq!(decision.active_reserved_bytes, target as u64);
            assert_eq!(decision.projected_private_bytes, Some(64 + target as u64));
            assert!(
                decision.projected_private_bytes.unwrap()
                    <= decision.limits.unwrap().computation_max_bytes
            );
            values.push(99);
            assert_eq!(values.last(), Some(&99));
        }
    }

    #[test]
    fn minimum_replacement_refusal_keeps_vector_and_credit_unchanged() {
        // Compute cap is 128: old 64 plus the minimum new 65 does not fit.
        let (budget, wave, _) = local_wave(160, true);
        let mut values = data();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            Wave::with(Some(&wave), || capacity(&mut values, 65));
        }));
        let payload = result.unwrap_err();
        let error = payload
            .downcast_ref::<WaveExhausted>()
            .expect("capacity rejection");
        assert!(!error.stale);
        assert_eq!(values.capacity(), 64);
        assert_eq!(values, data());
        unchanged_credit(&budget, &wave);
    }

    #[test]
    fn resident_growth_keeps_doubling_and_does_not_use_pressure_fallback() {
        let (budget, wave, _) = local_wave(192, false);
        let mut values = data();
        let rejected = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            Wave::with(Some(&wave), || capacity(&mut values, 65));
        }));
        assert!(rejected.unwrap_err().is::<WaveExhausted>());
        assert_eq!(values.capacity(), 64);
        unchanged_credit(&budget, &wave);
        let (_, wave, _) = local_wave(544, false);
        Wave::with(Some(&wave), || capacity(&mut values, 65));
        assert_eq!(values.capacity(), 128);
    }

    #[test]
    fn forced_rejection_is_not_absorbed_by_adaptive_growth_and_fires_once() {
        let (budget, wave, _) = local_wave(192, true);
        wave.0
            .fail_after
            .store(1, std::sync::atomic::Ordering::Relaxed);
        let mut values = data();
        let rejected = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            Wave::with(Some(&wave), || capacity(&mut values, 65));
        }));
        let payload = rejected.unwrap_err();
        let error = payload.downcast_ref::<WaveExhausted>().unwrap();
        assert_eq!(error.message, "test wave capacity rejection");
        assert_eq!(values.capacity(), 64);
        unchanged_credit(&budget, &wave);
        Wave::with(Some(&wave), || capacity(&mut values, 65));
        assert_eq!(values.capacity(), 72);
    }

    #[test]
    fn stale_observation_and_size_overflow_do_not_try_smaller_candidates() {
        let (budget, wave, _) = local_wave(192, true);
        budget.set_sample_age_for_test(Duration::from_secs(300));
        assert_eq!(budget.check(0).sample_age, Duration::from_secs(300));
        assert!(budget.check(0).stale);
        let mut values = data();
        let rejected = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            Wave::with(Some(&wave), || capacity(&mut values, 65));
        }));
        assert!(
            rejected
                .unwrap_err()
                .downcast_ref::<WaveExhausted>()
                .unwrap()
                .stale
        );
        assert_eq!(values.capacity(), 64);
        unchanged_credit(&budget, &wave);
        let (budget, wave, _) = local_wave(192, true);
        let mut values = Vec::<u64>::new();
        let result = crate::cache_source::catch_resource(|| {
            Wave::with(Some(&wave), || capacity(&mut values, usize::MAX));
            Ok(())
        });
        assert!(result.unwrap_err().contains("allocation size overflow"));
        assert_eq!(values.capacity(), 0);
        unchanged_credit(&budget, &wave);
    }
}

pub(super) fn reallocation(bytes: usize) {
    WAVE.with(|cell| {
        if let Some(wave) = cell.borrow().as_ref() {
            wave.allocate(bytes as u64);
        }
    });
}

pub(super) fn discovery(values: &mut Vec<State>, control: Option<&Control>) {
    if values.len() < values.capacity() {
        return;
    }
    let Some(control) = control else {
        return;
    };
    let target = values.capacity().saturating_mul(2).max(4);
    // This must precede assigning an ID, and must not hold an index lock.
    control.claim_slots(target - values.capacity());
    WAVE.with(|cell| {
        if let Some(wave) = cell.borrow().as_ref() {
            wave.allocate(bytes(target, control.state_bytes()));
        }
    });
    values.reserve_exact(target - values.len());
    if values.capacity() > target {
        control.claim_slots(values.capacity() - target);
    }
}

pub(super) fn wave_limit(pressure_view: bool, remaining: usize) -> usize {
    if pressure_view {
        remaining.min(WAVE_STATES)
    } else {
        remaining
    }
}

pub(super) fn refresh_normal() -> bool {
    let decision = process_budget().force_check(0);
    decision.pressure == Pressure::Normal && !decision.stale
}

pub(super) fn normal_at_boundary() -> bool {
    let decision = process_budget().check(0);
    if decision.pressure == Pressure::Normal && !decision.stale {
        true
    } else {
        refresh_normal()
    }
}

pub(super) fn back_plan(plan: &mut Plan, registry: &Registry) {
    if plan.workspace.is_some() {
        return;
    }
    let budget = process_budget();
    budget.force_check(0);
    let bytes = bytes(plan.terms.len(), 512 + size_of::<TermBlock>());
    let credit = initial_reservation(
        budget,
        bytes,
        ReservationKind::Progress,
        "plan workspace metadata",
    );
    credit
        .start()
        .unwrap_or_else(|error| failure("plan workspace metadata", error));
    let blocks = std::mem::take(&mut plan.terms);
    let count = blocks.len();
    let workspace = PlanWorkspace::new(blocks);
    // Retire even on failed registration; no hidden backend is created for a
    // source-only diagnostic Registry that cannot actually spill to disk.
    if let Err(error) = registry.register_workspace(workspace.clone()) {
        workspace.retire();
        failure("plan spill unavailable", error);
    }
    plan.terms.resize_with(count, TermBlock::default);
    plan.workspace = Some(workspace);
    #[cfg(test)]
    BACKED_PLANS.with(|cell| cell.set(cell.get() + 1));
    credit
        .settle(bytes)
        .unwrap_or_else(|error| failure("plan workspace metadata", error));
}

fn reclaim_stalled_at(
    started: Instant,
    last_reclaimed: Option<Instant>,
    now: Instant,
    timeout: Duration,
) -> bool {
    let since = last_reclaimed.map_or(started, |last| started.max(last));
    now.saturating_duration_since(since) >= timeout
}

#[cfg(test)]
mod reclaim_stall_tests {
    use super::*;

    #[test]
    fn sustained_real_reclamation_extends_the_stall_window_without_a_total_deadline() {
        let started = Instant::now();
        let timeout = Duration::from_secs(30);
        for elapsed in [45, 90, 150, 300] {
            let now = started + Duration::from_secs(elapsed);
            let reclaimed = Some(now - Duration::from_secs(4));
            assert!(!reclaim_stalled_at(started, reclaimed, now, timeout));
        }
        assert!(reclaim_stalled_at(
            started,
            None,
            started + timeout,
            timeout
        ));
        let last = started + Duration::from_secs(120);
        assert!(!reclaim_stalled_at(
            started,
            Some(last),
            last + timeout - Duration::from_nanos(1),
            timeout
        ));
        assert!(reclaim_stalled_at(
            started,
            Some(last),
            last + timeout,
            timeout
        ));
    }

    #[test]
    fn old_reclamation_does_not_shorten_a_new_wait_and_metadata_has_no_clock() {
        let first = Instant::now();
        let started = first + Duration::from_secs(60);
        let timeout = Duration::from_secs(30);
        assert!(!reclaim_stalled_at(
            started,
            Some(first),
            started + Duration::from_secs(29),
            timeout
        ));
        // Stored/claimed/requested/epoch changes intentionally supply no new
        // last_reclaimed timestamp, even after arbitrarily many such changes.
        assert!(reclaim_stalled_at(
            started,
            Some(first),
            started + timeout,
            timeout
        ));
    }
}

/// Wait while real block reclamation continues, rather than imposing a total
/// I/O deadline on a long spill backlog. Epoch/claim/store activity alone never
/// extends the stall window. Admission still needs a fresh native observation.
pub(super) fn retry(workspace: &PlanWorkspace, started: Instant) -> bool {
    retry_with_request(workspace, started, None)
}

fn retry_for(workspace: &PlanWorkspace, started: Instant, required: u64) -> bool {
    retry_with_request(workspace, started, Some(required))
}

fn retry_with_request(workspace: &PlanWorkspace, started: Instant, required: Option<u64>) -> bool {
    let budget = process_budget();
    let decision = budget.force_check(required.unwrap_or(0));
    if required.is_some() && compute_fits(&decision) {
        // A release or settled-credit absorption may already have made the
        // original request fit. Atomic admission is still retried by the caller.
        return true;
    }
    // These observations must not scan the historical slot directory while
    // holding the source mutex. Detailed payload/pin counts are failure-only.
    let progress = workspace.progress().unwrap_or_else(|error| {
        failure(
            "plan workspace unavailable",
            format!(
                "{error}; pressure={decision:?}; workspace={:?}",
                workspace.snapshot()
            ),
        )
    });
    // Observe the epoch before requesting work. A fast reclaimer must not
    // finish its turn before we capture the epoch and leave us waiting for an
    // unrelated later change. Only planning lacks an exact request estimate.
    let reclaim = decision.reclaim_bytes.max(32 * 1024 * 1024);
    workspace.request_reclaim(usize::try_from(reclaim).unwrap_or(usize::MAX));
    if reclaim_stalled_at(
        started,
        progress.last_reclaimed_at,
        Instant::now(),
        PROGRESS_TIMEOUT,
    ) && workspace
        .reclamation_stalled(started, PROGRESS_TIMEOUT)
        .unwrap_or_else(|error| {
            failure(
                "check plan spill progress",
                format!(
                    "{error}; pressure={decision:?}; workspace={:?}",
                    workspace.snapshot()
                ),
            )
        })
    {
        fail_stalled(workspace, required);
    }
    let after = workspace
        .wait_for_progress(progress.progress_epoch, PROGRESS_WAIT)
        .unwrap_or_else(|error| {
            failure(
                "wait for plan spill",
                format!(
                    "{error}; pressure={decision:?}; workspace={:?}",
                    workspace.snapshot()
                ),
            )
        });
    budget.force_check(required.unwrap_or(0));
    // This sample grants no credit by itself. The caller retries its original
    // reservation, against native usage and all process-wide liabilities.
    if reclaim_stalled_at(
        started,
        after.last_reclaimed_at,
        Instant::now(),
        PROGRESS_TIMEOUT,
    ) && workspace
        .reclamation_stalled(started, PROGRESS_TIMEOUT)
        .unwrap_or_else(|error| {
            failure(
                "check plan spill progress",
                format!(
                    "{error}; pressure={:?}; workspace={:?}",
                    budget.check(0),
                    workspace.snapshot()
                ),
            )
        })
    {
        fail_stalled(workspace, required);
    }
    true
}

fn fail_stalled(workspace: &PlanWorkspace, required: Option<u64>) -> ! {
    let decision = process_budget().check(required.unwrap_or(0));
    // Count historical slots/pins only when an actual stall is about to fail.
    // A transient Hard observation alone does not justify this full traversal.
    let snapshot = workspace.snapshot();
    if snapshot.resident_bytes == 0
        && snapshot.releasing_bytes == 0
        && decision.pressure == Pressure::Hard
    {
        failure(
            "no resident workspace payload remains to relieve pressure",
            format!(
                "blocks={}, reclaimed_payload_bytes={}, required_bytes={required:?}; pressure={decision:?}; workspace={snapshot:?}",
                snapshot.blocks, snapshot.reclaimed_bytes
            ),
        );
    }
    failure(
        "workspace made no actual payload reclamation progress for 30 seconds",
        format!("required_bytes={required:?}; pressure={decision:?}; workspace={snapshot:?}"),
    );
}

/// A deliberately conservative working estimate, not a formal upper bound for
/// arbitrary rational expressions. No arithmetic-item probe is needed: completed
/// Q payloads are observed once at the existing layer/group barrier.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct OutputMemory {
    peak_heap: u64,
}

impl OutputMemory {
    fn per_output(self) -> u64 {
        self.peak_heap.saturating_mul(8).max(4096)
    }

    fn containers(count: usize, buffered: bool) -> u64 {
        if buffered {
            // Jobs' reply buffers coexist with the collecting Vec and its old
            // allocation during growth. The original scheduler is unchanged.
            bytes(count, 4 * size_of::<(usize, Q)>())
        } else {
            0
        }
    }

    fn scratch(self, workers: usize) -> u64 {
        self.per_output()
            .checked_mul(workers.max(1) as u64)
            .and_then(|value| value.checked_mul(16))
            .unwrap_or_else(|| resource_failure("evaluation scratch estimate overflow"))
    }

    pub(super) fn estimate(self, count: usize, workers: usize, buffered: bool) -> u64 {
        self.per_output()
            .checked_mul(count as u64)
            .and_then(|value| value.checked_add(Self::containers(count, buffered)))
            .and_then(|value| value.checked_add(self.scratch(workers)))
            .unwrap_or_else(|| resource_failure("evaluation output estimate overflow"))
    }

    pub(super) fn heap(value: &Q) -> u64 {
        crate::cache_source::value_bytes(value).saturating_sub(size_of::<Q>()) as u64
    }

    fn observe(
        &mut self,
        count: usize,
        workers: usize,
        buffered: bool,
        heap_bytes: u64,
        peak_heap: u64,
        other_bytes: u64,
    ) -> u64 {
        self.peak_heap = self.peak_heap.max(peak_heap);
        heap_bytes
            .checked_add(Self::containers(count, buffered))
            .and_then(|value| value.checked_add(self.scratch(workers)))
            .and_then(|value| value.checked_add(other_bytes))
            .unwrap_or_else(|| resource_failure("observed evaluation memory overflow"))
    }

    pub(super) fn finish(
        &mut self,
        credit: Reservation,
        count: usize,
        workers: usize,
        buffered: bool,
        heap_bytes: u64,
        peak_heap: u64,
        other_bytes: u64,
    ) {
        let observed = self.observe(count, workers, buffered, heap_bytes, peak_heap, other_bytes);
        // A broken estimate is reported even though the allocation has already
        // happened. Settlement retains the complete declared envelope until a
        // later native sample covers it, including allocator-retained pages.
        credit
            .settle(observed)
            .unwrap_or_else(|error| failure("evaluation output exceeded its reservation", error));
    }
}

#[cfg(test)]
mod output_memory_tests {
    use super::*;

    #[test]
    fn an_unexpected_large_exact_result_exceeds_the_old_envelope_and_updates_the_next() {
        let large = sc_num::pow_unsigned(&Q::from(2), 65_536);
        let heap = OutputMemory::heap(&large);
        let mut outputs = OutputMemory::default();
        let admitted = outputs.estimate(1, 1, true);
        let observed = outputs.observe(1, 1, true, heap, heap, 0);
        // finish forwards this full observation to Reservation::settle, whose
        // underestimation path must report the already-allocated overshoot.
        assert!(observed > admitted);
        assert!(outputs.estimate(1, 1, true) >= observed);
        assert_eq!(OutputMemory::heap(&Q::ZERO), 0);
    }

    #[test]
    fn an_unrepresentable_output_envelope_is_a_resource_error() {
        let result = crate::cache_source::catch_resource(|| {
            OutputMemory {
                peak_heap: u64::MAX,
            }
            .estimate(2, 1, true);
            Ok(())
        });
        assert!(result.unwrap_err().contains("estimate overflow"));
    }
}

// The previous layer's settled envelope can make a fresh-enough cached sample
// conservative. Refresh once before committing to the disk-backed evaluator.
// A Normal first admission retains its original one check and one attempt.
fn resident_layer_admission<T>(
    initial: (Pressure, bool),
    mut refresh: impl FnMut() -> (Pressure, bool),
    mut attempt: impl FnMut() -> Result<T, ReservationError>,
) -> Option<T> {
    let normal = |(pressure, stale): (Pressure, bool)| pressure == Pressure::Normal && !stale;
    let mut refreshed = false;
    if !normal(initial) {
        refreshed = true;
        if !normal(refresh()) {
            return None;
        }
    }
    match attempt() {
        Ok(credit) => Some(credit),
        Err(error)
            if !refreshed
                && matches!(
                    &error,
                    ReservationError::Pressure { .. }
                        | ReservationError::Stale(_)
                        | ReservationError::Unavailable
                ) =>
        {
            if !normal(refresh()) {
                return None;
            }
            attempt().ok()
        }
        Err(_) => None,
    }
}

#[cfg(test)]
mod resident_layer_admission_tests {
    use super::*;
    use std::cell::Cell;
    use std::collections::VecDeque;

    fn denial() -> ReservationError {
        ReservationError::Pressure {
            pressure: Pressure::Hard,
            projected: Some(1024),
            limits: None,
        }
    }

    #[test]
    fn normal_first_admission_has_no_extra_sample_or_attempt() {
        let attempts = Cell::new(0);
        assert_eq!(
            resident_layer_admission(
                (Pressure::Normal, false),
                || panic!("ordinary resident admission sampled memory"),
                || {
                    attempts.set(attempts.get() + 1);
                    Ok(17usize)
                },
            ),
            Some(17)
        );
        assert_eq!(attempts.get(), 1);
    }

    #[test]
    fn settled_warning_refreshes_before_unnecessary_backing() {
        let refreshes = Cell::new(0);
        let attempts = Cell::new(0);
        assert_eq!(
            resident_layer_admission(
                (Pressure::Warning, false),
                || {
                    refreshes.set(refreshes.get() + 1);
                    (Pressure::Normal, false)
                },
                || {
                    attempts.set(attempts.get() + 1);
                    Ok(23usize)
                },
            ),
            Some(23)
        );
        assert_eq!(refreshes.get(), 1);
        assert_eq!(attempts.get(), 1);
    }

    #[test]
    fn full_output_denial_retries_the_same_request_after_one_sample() {
        let mut replies = VecDeque::from([Err(denial()), Ok(29usize)]);
        let mut requests = Vec::new();
        let refreshes = Cell::new(0);
        let required = 4096u64;
        assert_eq!(
            resident_layer_admission(
                (Pressure::Normal, false),
                || {
                    refreshes.set(refreshes.get() + 1);
                    (Pressure::Normal, false)
                },
                || {
                    requests.push(required);
                    replies.pop_front().expect("extra full-layer attempt")
                },
            ),
            Some(29)
        );
        assert_eq!(requests, [4096, 4096]);
        assert_eq!(refreshes.get(), 1);
        assert!(replies.is_empty());
    }

    #[test]
    fn genuine_fresh_pressure_backs_without_repeated_sampling() {
        for pressure in [Pressure::Warning, Pressure::Hard] {
            let refreshes = Cell::new(0);
            assert_eq!(
                resident_layer_admission::<usize>(
                    (pressure, false),
                    || {
                        refreshes.set(refreshes.get() + 1);
                        (pressure, false)
                    },
                    || panic!("fresh pressure must pause before allocation"),
                ),
                None
            );
            assert_eq!(refreshes.get(), 1);
        }

        // Even after a Warning is cleared, an oversized full layer may still
        // be refused. That is one real fresh denial, not another sample loop.
        let refreshes = Cell::new(0);
        let attempts = Cell::new(0);
        assert_eq!(
            resident_layer_admission::<usize>(
                (Pressure::Warning, false),
                || {
                    refreshes.set(refreshes.get() + 1);
                    (Pressure::Normal, false)
                },
                || {
                    attempts.set(attempts.get() + 1);
                    Err(denial())
                },
            ),
            None
        );
        assert_eq!(refreshes.get(), 1);
        assert_eq!(attempts.get(), 1);
    }

    #[test]
    fn a_still_stale_sample_never_admits_the_resident_layer() {
        let refreshes = Cell::new(0);
        assert_eq!(
            resident_layer_admission::<usize>(
                (Pressure::Normal, true),
                || {
                    refreshes.set(refreshes.get() + 1);
                    (Pressure::Normal, true)
                },
                || panic!("stale sample must not allocate"),
            ),
            None
        );
        assert_eq!(refreshes.get(), 1);
    }
}

pub(super) fn admit_layer(
    outputs: OutputMemory,
    count: usize,
    workers: usize,
    buffered: bool,
) -> Option<Reservation> {
    let budget = process_budget();
    let decision = budget.check(0);
    let credit = resident_layer_admission(
        (decision.pressure, decision.stale),
        || {
            let decision = budget.force_check(0);
            (decision.pressure, decision.stale)
        },
        || {
            budget.try_reserve(
                outputs.estimate(count, workers, buffered),
                ReservationKind::Compute,
            )
        },
    )?;
    credit.start().ok()?;
    Some(credit)
}

pub(super) fn evaluate(
    plan: &Plan,
    source: &ParallelSource,
    batch: &Batch,
    base: usize,
    threads: usize,
) {
    evaluate_from(
        plan,
        source,
        batch,
        base,
        threads,
        0,
        OutputMemory::default(),
    );
}

#[derive(Clone, Copy)]
enum EvaluationPhase {
    Admission,
    Load,
    Compute,
    Writeback,
}

#[derive(Clone, Copy, Default)]
struct PartCounts {
    groups: usize,
    equations: usize,
}

struct EvaluationProbe {
    started: Instant,
    last_live_report: Instant,
    active: Option<(EvaluationPhase, Instant)>,
    phase_ns: [u128; 4],
    by_parts: Vec<PartCounts>,
    group_limit: usize,
    next_layer: usize,
    groups: usize,
    equations: usize,
    attempts: usize,
    shrinks: usize,
    retry_waits: usize,
    min_group: usize,
    max_group: usize,
    group_ns: u128,
    max_group_ns: u128,
    current_parts: usize,
    unfinished_group: bool,
}

impl EvaluationProbe {
    fn new(threads: usize, group_limit: usize, next_layer: usize) -> Self {
        let started = Instant::now();
        Self {
            started,
            last_live_report: started,
            active: None,
            phase_ns: [0; 4],
            by_parts: vec![PartCounts::default(); threads.max(1) + 1],
            group_limit,
            next_layer,
            groups: 0,
            equations: 0,
            attempts: 0,
            shrinks: 0,
            retry_waits: 0,
            min_group: usize::MAX,
            max_group: 0,
            group_ns: 0,
            max_group_ns: 0,
            current_parts: 0,
            unfinished_group: false,
        }
    }

    fn phase(&mut self, phase: Option<EvaluationPhase>) {
        let now = Instant::now();
        if let Some((previous, start)) = self.active.take() {
            self.phase_ns[previous as usize] += now.duration_since(start).as_nanos();
        }
        self.active = phase.map(|phase| (phase, now));
    }

    fn begin_group(&mut self) {
        self.unfinished_group = true;
        self.current_parts = 0;
        self.phase(Some(EvaluationPhase::Admission));
    }

    fn complete_group(&mut self, count: usize, started: Instant) {
        let elapsed = started.elapsed().as_nanos();
        self.groups = self.groups.saturating_add(1);
        self.equations = self.equations.saturating_add(count);
        self.min_group = self.min_group.min(count);
        self.max_group = self.max_group.max(count);
        self.group_ns += elapsed;
        self.max_group_ns = self.max_group_ns.max(elapsed);
        let parts = &mut self.by_parts[self.current_parts];
        parts.groups = parts.groups.saturating_add(1);
        parts.equations = parts.equations.saturating_add(count);
        self.unfinished_group = false;
        self.current_parts = 0;
    }

    fn report_if_due(
        &mut self,
        layer: usize,
        start: usize,
        count: usize,
        event: &str,
        required_bytes: u64,
    ) {
        let now = Instant::now();
        if now.duration_since(self.last_live_report) < Duration::from_secs(5) {
            return;
        }
        self.last_live_report = now;
        let mut phase_ns = self.phase_ns;
        let phase = match self.active {
            Some((phase, began)) => {
                phase_ns[phase as usize] += now.duration_since(began).as_nanos();
                match phase {
                    EvaluationPhase::Admission => "admission",
                    EvaluationPhase::Load => "load",
                    EvaluationPhase::Compute => "compute",
                    EvaluationPhase::Writeback => "writeback",
                }
            }
            None => "between_phases",
        };
        // Read the existing cached ledger only for a live diagnostic line.
        // Reporting must not sample memory or request reclamation itself.
        let decision = process_budget().check(0);
        let seconds = |ns: u128| ns as f64 / 1e9;
        eprintln!(
            "sc-faber: pressure evaluation live wall_s={:.6} groups={} equations={} layer={} start={} count={} phase={} event={} admission_attempts={} shrinks={} retry_waits={} required_bytes={} admission_s={:.6} load_s={:.6} compute_scope_s={:.6} writeback_heap_s={:.6} private_bytes={:?} projected_bytes={:?} active_bytes={} settled_bytes={} pressure={:?} limits={:?} sample_age_ms={} sample_generation={} stale={}",
            now.duration_since(self.started).as_secs_f64(),
            self.groups,
            self.equations,
            layer,
            start,
            count,
            phase,
            event,
            self.attempts,
            self.shrinks,
            self.retry_waits,
            required_bytes,
            seconds(phase_ns[EvaluationPhase::Admission as usize]),
            seconds(phase_ns[EvaluationPhase::Load as usize]),
            seconds(phase_ns[EvaluationPhase::Compute as usize]),
            seconds(phase_ns[EvaluationPhase::Writeback as usize]),
            decision.known_private_bytes,
            decision.projected_private_bytes,
            decision.active_reserved_bytes,
            decision.settled_reserved_bytes,
            decision.pressure,
            decision.limits,
            decision.sample_age.as_millis(),
            decision.sample_generation,
            decision.stale,
        );
        #[cfg(target_os = "linux")]
        if let Some(observation) = decision.memory.cgroup_observation {
            eprintln!("sc-faber: pressure cgroup {observation:?}");
        }
    }
}

impl Drop for EvaluationProbe {
    fn drop(&mut self) {
        self.phase(None);
        let wall = self.started.elapsed().as_nanos();
        let measured: u128 = self.phase_ns.iter().sum();
        let seconds = |ns: u128| ns as f64 / 1e9;
        eprintln!(
            "sc-faber: pressure evaluation summary from_layer={} thread_limit={} group_limit={} groups={} equations={} min_group={} max_group={} mean_group_eq={:.3} admission_attempts={} shrinks={} retry_waits={} admission_s={:.6} load_s={:.6} compute_scope_s={:.6} writeback_heap_s={:.6} other_s={:.6} wall_s={:.6} mean_group_ms={:.6} max_group_ms={:.6} unfinished_group={} unfinished_parts={}",
            self.next_layer,
            self.by_parts.len() - 1,
            self.group_limit,
            self.groups,
            self.equations,
            if self.groups == 0 { 0 } else { self.min_group },
            self.max_group,
            self.equations as f64 / self.groups.max(1) as f64,
            self.attempts,
            self.shrinks,
            self.retry_waits,
            seconds(self.phase_ns[EvaluationPhase::Admission as usize]),
            seconds(self.phase_ns[EvaluationPhase::Load as usize]),
            seconds(self.phase_ns[EvaluationPhase::Compute as usize]),
            seconds(self.phase_ns[EvaluationPhase::Writeback as usize]),
            seconds(wall.saturating_sub(measured)),
            seconds(wall),
            self.group_ns as f64 / self.groups.max(1) as f64 / 1e6,
            self.max_group_ns as f64 / 1e6,
            self.unfinished_group,
            self.current_parts,
        );
        // Parts count actual mathematical chunks, including the caller's chunk;
        // they are not a measurement of simultaneous OS CPU activity.
        for (parts, counts) in self
            .by_parts
            .iter()
            .enumerate()
            .filter(|(_, counts)| counts.groups != 0)
        {
            eprintln!(
                "sc-faber: pressure evaluation parts={} groups={} equations={}",
                parts, counts.groups, counts.equations
            );
        }
    }
}

pub(super) fn evaluate_from(
    plan: &Plan,
    source: &ParallelSource,
    batch: &Batch,
    base: usize,
    threads: usize,
    next_layer: usize,
    mut outputs: OutputMemory,
) {
    let workspace = plan.workspace.as_ref().expect("pressure plan workspace");
    let group_equations = match std::env::var("SC_FABER_CACHE_MEASURE_GROUP_EQUATIONS") {
        Ok(value) => value
            .parse::<usize>()
            .ok()
            .filter(|&value| value > 0)
            .unwrap_or_else(|| {
                resource_failure("SC_FABER_CACHE_MEASURE_GROUP_EQUATIONS must be positive")
            }),
        Err(std::env::VarError::NotPresent) => GROUP_EQUATIONS,
        Err(error) => resource_failure(format!("cannot read pressure group size: {error}")),
    };
    let mut probe = super::REPORT_STAGES
        .get()
        .copied()
        .unwrap_or(false)
        .then(|| EvaluationProbe::new(threads, group_equations, next_layer));
    for (layer_index, layer) in batch.layers.iter().enumerate().skip(next_layer) {
        let mut start = 0;
        while start < layer.len() {
            let started = Instant::now();
            if let Some(probe) = &mut probe {
                probe.begin_group();
            }
            let mut count = group_equations.min(layer.len() - start);
            let mut refreshed = false;
            let (credit, blocks, mut ids, other_bytes, workers) = loop {
                if let Some(probe) = &mut probe {
                    probe.attempts = probe.attempts.saturating_add(1);
                }
                // Only metadata is allocated before ensure(total). Reserving
                // future output/scratch here would charge that unused envelope
                // as settled each time block admission fails and shrinks.
                let metadata = bytes(count, 6 * size_of::<usize>() + 1024);
                let workers = threads.max(1).min(count.div_ceil(32).max(1));
                let output_bytes = outputs.estimate(count, workers, true);
                let budget = process_budget();
                let credit = match evaluation_reserve(budget, &mut refreshed, || {
                    budget.try_reserve(metadata, ReservationKind::Compute)
                }) {
                    Ok(credit) => credit,
                    Err(error) => {
                        if !capacity_denial(&error) {
                            failure("pressure evaluation metadata", error);
                        }
                        if count > 1 {
                            if let Some(probe) = &mut probe {
                                probe.shrinks = probe.shrinks.saturating_add(1);
                            }
                            count = count.div_ceil(2);
                            continue;
                        }
                        if let Some(probe) = &mut probe {
                            probe.retry_waits = probe.retry_waits.saturating_add(1);
                            probe.report_if_due(
                                layer_index,
                                start,
                                count,
                                "metadata_retry",
                                metadata,
                            );
                        }
                        if retry_for(workspace, started, metadata) {
                            continue;
                        }
                        failure("pressure evaluation metadata", error);
                    }
                };
                credit
                    .start()
                    .unwrap_or_else(|error| failure("pressure evaluation metadata", error));
                let mut ids = layer[start..start + count].to_vec();
                ids.sort_unstable_by_key(|&id| (plan.equations[id].block, id));
                let mut indices = ids
                    .iter()
                    .map(|&id| plan.equations[id].block)
                    .collect::<Vec<_>>();
                indices.dedup();
                let group = workspace
                    .pin_blocks(indices)
                    .unwrap_or_else(|error| failure("pin plan block group", error));
                // Existing Arc payload is already covered by native usage or
                // its original active/settled credit. Only missing blocks can
                // allocate new payload; their pins survive the entire group.
                let decode_bytes = group
                    .decode_bytes(DECODE_WINDOW)
                    .unwrap_or_else(|error| failure("pressure block admission", error));
                let other_bytes = metadata
                    .checked_add(decode_bytes)
                    .unwrap_or_else(|| resource_failure("pressure block admission overflow"));
                let total = other_bytes
                    .checked_add(output_bytes)
                    .unwrap_or_else(|| resource_failure("pressure evaluation admission overflow"));
                if let Err(error) =
                    evaluation_reserve(budget, &mut refreshed, || credit.ensure(total))
                {
                    // Pins must not prevent the independent reclaimer from
                    // making progress while admission shrinks or waits.
                    drop(group);
                    drop(ids);
                    drop(credit);
                    if !capacity_denial(&error) {
                        failure("pressure evaluation block group", error);
                    }
                    if count > 1 {
                        if let Some(probe) = &mut probe {
                            probe.shrinks = probe.shrinks.saturating_add(1);
                        }
                        count = count.div_ceil(2);
                        continue;
                    }
                    if let Some(probe) = &mut probe {
                        probe.retry_waits = probe.retry_waits.saturating_add(1);
                        probe.report_if_due(layer_index, start, count, "block_retry", total);
                    }
                    if retry_for(workspace, started, total) {
                        continue;
                    }
                    failure("pressure evaluation block group", error);
                }
                if let Some(probe) = &mut probe {
                    probe.phase(Some(EvaluationPhase::Load));
                    probe.report_if_due(layer_index, start, count, "load_group", total);
                }
                let blocks = workspace
                    .load_pinned_blocks(group, &credit)
                    .unwrap_or_else(|error| failure("restore plan blocks", error));
                if let Some(probe) = &mut probe {
                    probe.phase(None);
                }
                break (credit, blocks, ids, other_bytes, workers);
            };
            let completed = {
                let values = source.read();
                let chunk = ids.len().div_ceil(threads.max(1)).max(32);
                if let Some(probe) = &mut probe {
                    probe.current_parts = ids.len().div_ceil(chunk);
                    probe.phase(Some(EvaluationPhase::Compute));
                }
                let completed = std::thread::scope(|scope| {
                    let mut jobs = Vec::new();
                    for ids in ids.chunks(chunk).skip(1) {
                        let blocks = &blocks;
                        let values = &*values;
                        jobs.push(scope.spawn(move || calculate(plan, blocks, ids, values)));
                    }
                    let mut completed =
                        calculate(plan, &blocks, &ids[..chunk.min(ids.len())], &values);
                    for job in jobs {
                        completed.extend(
                            job.join()
                                .unwrap_or_else(|payload| std::panic::resume_unwind(payload)),
                        );
                    }
                    completed
                });
                if let Some(probe) = &mut probe {
                    probe.phase(None);
                }
                completed
            };
            if let Some(probe) = &mut probe {
                probe.phase(Some(EvaluationPhase::Writeback));
            }
            let mut heap_bytes = 0u64;
            let mut peak_heap = 0u64;
            for (_, value) in &completed {
                let heap = OutputMemory::heap(value);
                heap_bytes = heap_bytes
                    .checked_add(heap)
                    .unwrap_or_else(|| resource_failure("pressure output heap overflow"));
                peak_heap = peak_heap.max(heap);
            }
            {
                let mut values = source.write();
                for (id, value) in completed {
                    values[base + id] = value;
                }
            }
            if let Some(probe) = &mut probe {
                probe.phase(None);
            }
            ids.clear();
            drop(ids);
            drop(blocks);
            outputs.finish(
                credit,
                count,
                workers,
                true,
                heap_bytes,
                peak_heap,
                other_bytes,
            );
            start += count;
            if let Some(probe) = &mut probe {
                probe.complete_group(count, started);
                probe.report_if_due(layer_index, start, count, "group_complete", 0);
            }
        }
        // Same-rank groups cannot observe one another. Only a fully completed
        // rank is exposed to the asynchronous integral reader.
        source.publish(batch, layer_index);
    }
}

fn calculate(
    plan: &Plan,
    blocks: &[(usize, Arc<TermBlock>)],
    ids: &[usize],
    values: &[Q],
) -> Vec<(usize, Q)> {
    ids.iter()
        .map(|&id| {
            let equation = &plan.equations[id];
            let slot = blocks
                .binary_search_by_key(&equation.block, |(index, _)| *index)
                .expect("loaded pressure equation block");
            (
                id,
                blocks[slot].1.evaluate(
                    equation.terms.clone(),
                    equation.constant(),
                    equation.divisor(),
                    values,
                ),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::super::{Coefficient, Term};
    use super::*;

    #[test]
    fn block_group_uses_exact_factors_and_original_equation_ranges() {
        let mut left = TermBlock::default();
        left.push(Term {
            factor: Coefficient::Half(3),
            left: 0,
            right: 1,
        });
        left.push(Term {
            factor: Coefficient::Integer(-2),
            left: 1,
            right: usize::MAX,
        });
        let mut right = TermBlock::default();
        right.push(Term {
            factor: Coefficient::Rational(Box::new(Q::from_signeds(7, 13))),
            left: 0,
            right: usize::MAX,
        });
        let plan = Plan {
            equations: vec![
                PackedEquation::new(Q::from(5), Some(Q::from(3)), 1, 0..1),
                PackedEquation::new(Q::ZERO, None, 0, 1..2),
                PackedEquation::new(Q::from(1), None, 0, 0..1),
            ],
            terms: Vec::new(),
            layers: Vec::new(),
            restored: Vec::new(),
            workspace: None,
            ..Plan::default()
        };
        let values = vec![Q::from(2), Q::from(3)];
        let blocks = vec![(0, Arc::new(left)), (1, Arc::new(right))];
        let actual = calculate(&plan, &blocks, &[1, 2, 0], &values);
        assert_eq!(
            actual,
            vec![
                (1, Q::from(-6)),
                (2, Q::from(10)),
                (0, Q::from_signeds(79, 39))
            ]
        );
    }
}
