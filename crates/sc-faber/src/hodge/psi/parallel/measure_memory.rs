//! Native memory observations for manual diagnostics, never production timing.
//!
//! Windows counters include allocations outside Rust (including SQLite and the
//! allocator). They do not attribute memory to individual Rust containers.
//! Native peak counters cover the entire process lifetime; phase summaries below
//! instead take maxima of the current counters sampled during that phase.
#![cfg(test)]
#![allow(dead_code)]

use std::io;
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Phase {
    Idle,
    Planning,
    Evaluation,
    Release,
}

impl Phase {
    pub(super) const ALL: [Self; 4] = [Self::Idle, Self::Planning, Self::Evaluation, Self::Release];
}

#[derive(Clone, Copy, Debug)]
pub(super) struct ProcessMemory {
    pub working_set_bytes: usize,
    pub private_commit_bytes: usize,
    pub lifetime_peak_working_set_bytes: usize,
    /// The native PeakPagefileUsage counter measures peak process commit, not I/O.
    pub lifetime_peak_commit_bytes: usize,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct SystemMemory {
    pub total_physical_bytes: u64,
    pub available_physical_bytes: u64,
    pub commit_total_bytes: u64,
    pub commit_limit_bytes: u64,
    pub available_commit_bytes: u64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MemorySnapshot {
    pub process: Option<ProcessMemory>,
    pub system: Option<SystemMemory>,
}

/// Calls the operating system immediately. Individual unavailable counters are None.
pub(super) fn snapshot() -> MemorySnapshot {
    MemorySnapshot {
        process: platform::process_memory(),
        system: platform::system_memory(),
    }
}

/// Accumulated kernel plus user CPU time of the calling thread, in nanoseconds.
/// Intended for diagnostic worker start/end samples, not mathematical hot loops.
pub(super) fn current_thread_cpu_ns() -> Option<u128> {
    platform::current_thread_cpu_ns()
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MemorySample {
    pub elapsed: Duration,
    pub phase: Phase,
    pub label: Option<&'static str>,
    pub memory: MemorySnapshot,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct PhaseSummary {
    pub phase: Phase,
    pub samples: usize,
    pub sampled_peak_working_set_bytes: Option<usize>,
    pub sampled_peak_private_commit_bytes: Option<usize>,
    pub minimum_available_physical_bytes: Option<u64>,
    pub minimum_available_commit_bytes: Option<u64>,
}

#[derive(Debug)]
pub(super) struct SamplingReport {
    pub interval: Duration,
    pub worker_panicked: bool,
    pub samples: Vec<MemorySample>,
}

impl SamplingReport {
    pub(super) fn phase_summaries(&self) -> [PhaseSummary; 4] {
        Phase::ALL.map(|phase| {
            let mut summary = PhaseSummary {
                phase,
                samples: 0,
                sampled_peak_working_set_bytes: None,
                sampled_peak_private_commit_bytes: None,
                minimum_available_physical_bytes: None,
                minimum_available_commit_bytes: None,
            };
            for sample in self.samples.iter().filter(|sample| sample.phase == phase) {
                summary.samples += 1;
                if let Some(process) = sample.memory.process {
                    summary.sampled_peak_working_set_bytes = Some(
                        summary
                            .sampled_peak_working_set_bytes
                            .unwrap_or(0)
                            .max(process.working_set_bytes),
                    );
                    summary.sampled_peak_private_commit_bytes = Some(
                        summary
                            .sampled_peak_private_commit_bytes
                            .unwrap_or(0)
                            .max(process.private_commit_bytes),
                    );
                }
                if let Some(system) = sample.memory.system {
                    summary.minimum_available_physical_bytes = Some(
                        summary
                            .minimum_available_physical_bytes
                            .unwrap_or(u64::MAX)
                            .min(system.available_physical_bytes),
                    );
                    summary.minimum_available_commit_bytes = Some(
                        summary
                            .minimum_available_commit_bytes
                            .unwrap_or(u64::MAX)
                            .min(system.available_commit_bytes),
                    );
                }
            }
            summary
        })
    }
}

struct State {
    phase: Phase,
    stop: bool,
    samples: Vec<MemorySample>,
}

struct Shared {
    started: Instant,
    state: Mutex<State>,
    changed: Condvar,
}

impl Shared {
    fn sample(&self, state: &mut State, label: Option<&'static str>) -> MemorySample {
        let sample = MemorySample {
            elapsed: self.started.elapsed(),
            phase: state.phase,
            label,
            memory: snapshot(),
        };
        state.samples.push(sample);
        sample
    }
}

/// A process sampler for one isolated diagnostic run. It performs no file I/O.
/// Its own thread and sample storage are included in the native process counters.
pub(super) struct Sampler {
    shared: Arc<Shared>,
    interval: Duration,
    thread: Option<JoinHandle<()>>,
    worker_panicked: bool,
}

impl Sampler {
    /// Accept 1..=20 ms. Polling is diagnostic overhead and must not be timed as production.
    pub(super) fn start(interval: Duration) -> io::Result<Self> {
        if interval < Duration::from_millis(1) || interval > Duration::from_millis(20) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "memory sampling interval must be between 1 and 20 milliseconds",
            ));
        }
        let shared = Arc::new(Shared {
            started: Instant::now(),
            state: Mutex::new(State {
                phase: Phase::Idle,
                stop: false,
                samples: Vec::new(),
            }),
            changed: Condvar::new(),
        });
        {
            let mut state = shared
                .state
                .lock()
                .unwrap_or_else(|error| error.into_inner());
            shared.sample(&mut state, Some("sampler_start"));
        }
        let worker = shared.clone();
        let thread = std::thread::Builder::new()
            .name("faber-memory-sampler".into())
            .spawn(move || {
                let mut state = worker
                    .state
                    .lock()
                    .unwrap_or_else(|error| error.into_inner());
                while !state.stop {
                    let result = worker
                        .changed
                        .wait_timeout(state, interval)
                        .unwrap_or_else(|error| error.into_inner());
                    state = result.0;
                    if !state.stop {
                        worker.sample(&mut state, None);
                    }
                }
            })?;
        Ok(Self {
            shared,
            interval,
            thread: Some(thread),
            worker_panicked: false,
        })
    }

    /// Changing the phase also records an immediate sample, so short phases remain visible.
    /// Phase changes and samples serialize only with this diagnostic sampler.
    pub(super) fn set_phase(&self, phase: Phase) {
        let mut state = self
            .shared
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        state.phase = phase;
        self.shared.sample(&mut state, None);
    }

    /// Record a named boundary immediately and return the same structured observation.
    pub(super) fn sample(&self, label: &'static str) -> MemorySample {
        let mut state = self
            .shared
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        self.shared.sample(&mut state, Some(label))
    }

    fn stop_worker(&mut self) {
        {
            let mut state = self
                .shared
                .state
                .lock()
                .unwrap_or_else(|error| error.into_inner());
            if !state.stop {
                self.shared.sample(&mut state, Some("sampler_stop"));
                state.stop = true;
            }
            self.shared.changed.notify_all();
        }
        if let Some(thread) = self.thread.take() {
            self.worker_panicked = thread.join().is_err();
        }
    }

    pub(super) fn stop(mut self) -> SamplingReport {
        self.stop_worker();
        let mut state = self
            .shared
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        SamplingReport {
            interval: self.interval,
            worker_panicked: self.worker_panicked,
            samples: std::mem::take(&mut state.samples),
        }
    }
}

impl Drop for Sampler {
    fn drop(&mut self) {
        self.stop_worker();
    }
}

#[cfg(target_os = "windows")]
mod platform {
    use super::{ProcessMemory, SystemMemory};
    use std::ffi::c_void;

    #[repr(C)]
    #[derive(Default)]
    struct FileTime {
        low_date_time: u32,
        high_date_time: u32,
    }

    impl FileTime {
        fn ticks(&self) -> u128 {
            (u128::from(self.high_date_time) << 32) | u128::from(self.low_date_time)
        }
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
        kernel_nonpaged: usize,
        page_size: usize,
        handle_count: u32,
        process_count: u32,
        thread_count: u32,
    }

    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetCurrentProcess() -> *mut c_void;
        fn GetCurrentThread() -> *mut c_void;
        fn GetThreadTimes(
            thread: *mut c_void,
            creation: *mut FileTime,
            exit: *mut FileTime,
            kernel: *mut FileTime,
            user: *mut FileTime,
        ) -> i32;
    }

    pub(super) fn current_thread_cpu_ns() -> Option<u128> {
        let mut creation = FileTime::default();
        let mut exit = FileTime::default();
        let mut kernel = FileTime::default();
        let mut user = FileTime::default();
        // SAFETY: FileTime has the Windows FILETIME layout. The current-thread
        // pseudo-handle remains valid for this call and must not be closed.
        if unsafe {
            GetThreadTimes(
                GetCurrentThread(),
                &mut creation,
                &mut exit,
                &mut kernel,
                &mut user,
            )
        } == 0
        {
            return None;
        }
        // GetThreadTimes reports CPU time in 100-nanosecond units.
        Some((kernel.ticks() + user.ticks()) * 100)
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

    pub(super) fn process_memory() -> Option<ProcessMemory> {
        let mut counters = ProcessCounters {
            cb: std::mem::size_of::<ProcessCounters>() as u32,
            ..ProcessCounters::default()
        };
        // SAFETY: both structures use the documented Windows C ABI. The current
        // process pseudo-handle is valid without acquiring or closing a handle.
        if unsafe { GetProcessMemoryInfo(GetCurrentProcess(), &mut counters, counters.cb) } == 0 {
            return None;
        }
        Some(ProcessMemory {
            working_set_bytes: counters.working_set_size,
            private_commit_bytes: counters.private_usage,
            lifetime_peak_working_set_bytes: counters.peak_working_set_size,
            lifetime_peak_commit_bytes: counters.peak_pagefile_usage,
        })
    }

    pub(super) fn system_memory() -> Option<SystemMemory> {
        let mut info = PerformanceInfo {
            cb: std::mem::size_of::<PerformanceInfo>() as u32,
            ..PerformanceInfo::default()
        };
        // SAFETY: info has the documented PERFORMANCE_INFORMATION layout and size.
        if unsafe { GetPerformanceInfo(&mut info, info.cb) } == 0 || info.page_size == 0 {
            return None;
        }
        let page = info.page_size as u64;
        Some(SystemMemory {
            total_physical_bytes: (info.physical_total as u64).saturating_mul(page),
            available_physical_bytes: (info.physical_available as u64).saturating_mul(page),
            commit_total_bytes: (info.commit_total as u64).saturating_mul(page),
            commit_limit_bytes: (info.commit_limit as u64).saturating_mul(page),
            available_commit_bytes: (info.commit_limit.saturating_sub(info.commit_total) as u64)
                .saturating_mul(page),
        })
    }
}

#[cfg(not(target_os = "windows"))]
mod platform {
    use super::{ProcessMemory, SystemMemory};

    pub(super) fn current_thread_cpu_ns() -> Option<u128> {
        None
    }

    pub(super) fn process_memory() -> Option<ProcessMemory> {
        None
    }
    pub(super) fn system_memory() -> Option<SystemMemory> {
        None
    }
}

#[test]
fn named_boundaries_and_phase_summaries_survive_sampler_stop() {
    let sampler = Sampler::start(Duration::from_millis(20)).unwrap();
    sampler.set_phase(Phase::Planning);
    let boundary = sampler.sample("plan_ready");
    assert_eq!(boundary.phase, Phase::Planning);
    assert_eq!(boundary.label, Some("plan_ready"));
    sampler.set_phase(Phase::Release);
    let report = sampler.stop();
    assert!(!report.worker_panicked);
    assert!(
        report
            .samples
            .windows(2)
            .all(|pair| pair[0].elapsed <= pair[1].elapsed)
    );
    assert!(
        report
            .samples
            .iter()
            .any(|sample| sample.label == Some("plan_ready"))
    );
    assert!(
        report
            .phase_summaries()
            .iter()
            .any(|phase| phase.phase == Phase::Release && phase.samples > 0)
    );
}

#[cfg(target_os = "windows")]
#[test]
fn native_windows_counters_are_available_for_the_current_process() {
    let memory = snapshot();
    let process = memory.process.expect("GetProcessMemoryInfo failed");
    let system = memory.system.expect("GetPerformanceInfo failed");
    assert!(process.working_set_bytes > 0);
    assert!(process.lifetime_peak_working_set_bytes >= process.working_set_bytes);
    assert!(system.total_physical_bytes >= system.available_physical_bytes);
    assert!(system.commit_limit_bytes >= system.commit_total_bytes);
    assert_eq!(
        system.available_commit_bytes,
        system.commit_limit_bytes - system.commit_total_bytes
    );
}
