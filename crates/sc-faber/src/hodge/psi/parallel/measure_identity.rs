//! Test-only identity-source experiment. No values or SQL are involved.
//!
//! Hooks observe an existing `Arc<Vec<State>>` once per planning frontier. The
//! calling thread retains one batch reference and records an event, without
//! scanning the states. Payload estimates are collected by the consumer or
//! after the measured calculation. They describe retained source allocations,
//! not additional physical process memory while mathematics also owns a batch.
//! The catalog contains only existing identity batches, never computed results.
//! This diagnostic has no production budget or persistence. `finish` drains and
//! releases the reverse directory; it does not implement a long-lived source.
//! Busy/wait/drain durations are wall time with overlapping scopes. Consumer
//! CPU time covers the whole consumer thread, including directory release.

use super::State;
use std::cell::RefCell;
use std::collections::VecDeque;
use std::mem::size_of;
use std::rc::Rc;
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Mode {
    Off,
    Retain,
    Consume,
    Packed,
}

#[derive(Clone, Debug)]
pub(super) struct Report {
    pub mode: Mode,
    pub observed_batches: usize,
    pub observed_states: usize,
    pub scanned_states: usize,
    pub duplicate_ids: usize,
    pub peak_retained_payload_bytes: usize,
    pub directory_slots: usize,
    pub directory_capacity: usize,
    pub directory_payload_bytes: usize,
    pub packed_key_bytes: usize,
    pub packed_key_capacity: usize,
    pub observe_elapsed: Duration,
    pub observe_lock_wait: Duration,
    pub consumer_busy: Duration,
    pub consumer_wait: Duration,
    pub consumer_cpu_ns: Option<u128>,
    pub drain_elapsed: Duration,
    pub statistics_elapsed: Duration,
    pub release_elapsed: Duration,
    pub finish_elapsed: Duration,
}

impl Report {
    pub(super) fn print(&self) {
        eprintln!(
            "identity mode={:?} observed_batches={} observed_states={} scanned_states={} duplicate_ids={} peak_retained_payload_bytes={} directory_slots={} directory_capacity={} directory_payload_bytes={} packed_key_bytes={} packed_key_capacity={} observe={:?} observe_lock_wait={:?} consumer_busy={:?} consumer_wait={:?} consumer_cpu_ns={:?} drain={:?} statistics={:?} release={:?} finish={:?}",
            self.mode,
            self.observed_batches,
            self.observed_states,
            self.scanned_states,
            self.duplicate_ids,
            self.peak_retained_payload_bytes,
            self.directory_slots,
            self.directory_capacity,
            self.directory_payload_bytes,
            self.packed_key_bytes,
            self.packed_key_capacity,
            self.observe_elapsed,
            self.observe_lock_wait,
            self.consumer_busy,
            self.consumer_wait,
            self.consumer_cpu_ns,
            self.drain_elapsed,
            self.statistics_elapsed,
            self.release_elapsed,
            self.finish_elapsed,
        );
    }
}

thread_local! {
    static ACTIVE: RefCell<Option<Rc<RefCell<IdentityProbe>>>> = const { RefCell::new(None) };
}

/// Nested runs must be finished or dropped in stack order.
pub(super) struct Run {
    probe: Rc<RefCell<IdentityProbe>>,
    previous: Option<Rc<RefCell<IdentityProbe>>>,
    active: bool,
}

impl Run {
    pub(super) fn start(mode: Mode) -> Self {
        let probe = Rc::new(RefCell::new(IdentityProbe::new(mode)));
        let previous =
            ACTIVE.with(|slot| slot.replace((mode != Mode::Off).then(|| Rc::clone(&probe))));
        Self {
            probe,
            previous,
            active: true,
        }
    }

    pub(super) fn finish(mut self) -> Report {
        self.deactivate();
        self.probe.borrow_mut().finish()
    }

    fn deactivate(&mut self) {
        if self.active {
            ACTIVE.with(|slot| {
                slot.replace(self.previous.take());
            });
            self.active = false;
        }
    }
}

impl Drop for Run {
    fn drop(&mut self) {
        self.deactivate();
    }
}

/// Called only where the original planner already constructs a state batch.
pub(super) fn observe_states(states: &Arc<Vec<State>>) {
    ACTIVE.with(|slot| {
        if let Some(probe) = slot.borrow().as_ref() {
            probe.borrow_mut().observe_states(states);
        }
    });
}

struct Batch {
    states: Arc<Vec<State>>,
    observed: u64,
}

#[derive(Default)]
struct Catalog {
    batches: VecDeque<Batch>,
    event: u64,
    stop: bool,
}

#[derive(Default)]
struct Shared {
    catalog: Mutex<Catalog>,
    changed: Condvar,
}

struct Lifetime {
    observed: u64,
    released: u64,
    bytes: usize,
}

struct Identity {
    _g: u32,
    d: Vec<u32>,
    ch: Vec<(u32, u32)>,
}

struct DirectorySummary {
    slots: usize,
    capacity: usize,
    payload_bytes: usize,
    packed_key_bytes: usize,
    packed_key_capacity: usize,
}

trait Directory: Default {
    /// False means the runtime ID was already present.
    fn insert(&mut self, state: &State) -> bool;
    fn summary(&self) -> DirectorySummary;
    fn observe(&self);
}

#[derive(Default)]
struct ExpandedDirectory {
    entries: Vec<Option<Identity>>,
    heap_bytes: usize,
}

impl Directory for ExpandedDirectory {
    fn insert(&mut self, state: &State) -> bool {
        if self.entries.len() <= state.index {
            self.entries.resize_with(state.index + 1, || None);
        }
        if self.entries[state.index].is_some() {
            return false;
        }
        let identity = Identity {
            _g: state.g,
            d: state.d.clone(),
            ch: state.ch.clone(),
        };
        self.heap_bytes += identity.d.capacity() * size_of::<u32>()
            + identity.ch.capacity() * size_of::<(u32, u32)>();
        self.entries[state.index] = Some(identity);
        true
    }

    fn summary(&self) -> DirectorySummary {
        DirectorySummary {
            slots: self.entries.len(),
            capacity: self.entries.capacity(),
            payload_bytes: self.entries.capacity() * size_of::<Option<Identity>>()
                + self.heap_bytes,
            packed_key_bytes: 0,
            packed_key_capacity: 0,
        }
    }

    fn observe(&self) {
        std::hint::black_box(self.entries.as_slice());
    }
}

/// One arena of mathematical identities and one offset per runtime ID. Each
/// record has a four-byte little-endian body length, then ULEB128 fields:
/// genus, marking count, all psi exponents, C-entry count, and (index, power).
/// The runtime ID and redundant rank are not part of mathematical identity.
#[derive(Default)]
struct PackedDirectory {
    offsets: Vec<usize>,
    keys: Vec<u8>,
}

fn append_u32(bytes: &mut Vec<u8>, mut value: u32) {
    while value >= 128 {
        bytes.push(value as u8 | 128);
        value >>= 7;
    }
    bytes.push(value as u8);
}

impl Directory for PackedDirectory {
    fn insert(&mut self, state: &State) -> bool {
        if self.offsets.len() <= state.index {
            self.offsets.resize(state.index + 1, usize::MAX);
        }
        if self.offsets[state.index] != usize::MAX {
            return false;
        }
        let start = self.keys.len();
        self.keys.extend_from_slice(&[0; 4]);
        append_u32(&mut self.keys, state.g);
        append_u32(
            &mut self.keys,
            u32::try_from(state.d.len()).expect("identity marking count exceeds u32"),
        );
        for &power in &state.d {
            append_u32(&mut self.keys, power);
        }
        append_u32(
            &mut self.keys,
            u32::try_from(state.ch.len()).expect("identity C-entry count exceeds u32"),
        );
        for &(index, power) in &state.ch {
            append_u32(&mut self.keys, index);
            append_u32(&mut self.keys, power);
        }
        let length = u32::try_from(self.keys.len() - start - 4)
            .expect("one identity record exceeds u32 byte length");
        self.keys[start..start + 4].copy_from_slice(&length.to_le_bytes());
        self.offsets[state.index] = start;
        true
    }

    fn summary(&self) -> DirectorySummary {
        DirectorySummary {
            slots: self.offsets.len(),
            capacity: self.offsets.capacity(),
            payload_bytes: self.offsets.capacity() * size_of::<usize>() + self.keys.capacity(),
            packed_key_bytes: self.keys.len(),
            packed_key_capacity: self.keys.capacity(),
        }
    }

    fn observe(&self) {
        std::hint::black_box(self.offsets.as_slice());
        std::hint::black_box(self.keys.as_slice());
    }
}

impl PackedDirectory {
    fn record(&self, id: usize) -> Option<&[u8]> {
        let start = *self.offsets.get(id)?;
        let length_end = start.checked_add(4)?;
        let length =
            u32::from_le_bytes(self.keys.get(start..length_end)?.try_into().ok()?) as usize;
        self.keys.get(length_end..length_end.checked_add(length)?)
    }
}

#[derive(Default)]
struct ConsumerReport {
    lifetimes: Vec<Lifetime>,
    scanned_states: usize,
    duplicate_ids: usize,
    directory_slots: usize,
    directory_capacity: usize,
    directory_payload_bytes: usize,
    packed_key_bytes: usize,
    packed_key_capacity: usize,
    busy: Duration,
    wait: Duration,
    release: Duration,
    cpu_ns: Option<u128>,
}

struct IdentityProbe {
    mode: Mode,
    retained: Vec<Batch>,
    shared: Option<Arc<Shared>>,
    consumer: Option<JoinHandle<ConsumerReport>>,
    observed_batches: usize,
    observed_states: usize,
    observe_elapsed: Duration,
    observe_lock_wait: Duration,
}

impl IdentityProbe {
    fn new(mode: Mode) -> Self {
        let shared =
            matches!(mode, Mode::Consume | Mode::Packed).then(|| Arc::new(Shared::default()));
        let consumer = shared.as_ref().map(|shared| {
            let shared = Arc::clone(shared);
            std::thread::Builder::new()
                .name("faber-identity-probe".into())
                .spawn(move || consume(&shared, mode))
                .expect("start identity consumer")
        });
        Self {
            mode,
            retained: Vec::new(),
            shared,
            consumer,
            observed_batches: 0,
            observed_states: 0,
            observe_elapsed: Duration::ZERO,
            observe_lock_wait: Duration::ZERO,
        }
    }

    fn observe_states(&mut self, states: &Arc<Vec<State>>) {
        if self.mode == Mode::Off {
            return;
        }
        let started = Instant::now();
        self.observed_batches += 1;
        self.observed_states += states.len();
        match self.mode {
            Mode::Off => unreachable!(),
            Mode::Retain => self.retained.push(Batch {
                states: Arc::clone(states),
                observed: self.observed_batches as u64,
            }),
            Mode::Consume | Mode::Packed => {
                let shared = self.shared.as_ref().unwrap();
                let waiting = Instant::now();
                let mut catalog = shared.catalog.lock().unwrap();
                self.observe_lock_wait += waiting.elapsed();
                catalog.event += 1;
                let observed = catalog.event;
                catalog.batches.push_back(Batch {
                    states: Arc::clone(states),
                    observed,
                });
                shared.changed.notify_one();
            }
        }
        self.observe_elapsed += started.elapsed();
    }

    fn finish(&mut self) -> Report {
        let started = Instant::now();
        let draining = Instant::now();
        let consumer = self.consumer.take().map(|consumer| {
            self.stop();
            consumer.join().expect("identity consumer panicked")
        });
        let drain_elapsed = draining.elapsed();
        let statistics = Instant::now();
        let peak_retained_payload_bytes = if let Some(consumer) = &consumer {
            retained_peak(&consumer.lifetimes)
        } else {
            self.retained
                .iter()
                .map(|batch| source_payload_bytes(&batch.states))
                .sum()
        };
        let statistics_elapsed = statistics.elapsed();
        let releasing = Instant::now();
        self.retained.clear();
        let release_elapsed = releasing.elapsed();
        let consumer = consumer.unwrap_or_default();
        Report {
            mode: self.mode,
            observed_batches: self.observed_batches,
            observed_states: self.observed_states,
            scanned_states: consumer.scanned_states,
            duplicate_ids: consumer.duplicate_ids,
            peak_retained_payload_bytes,
            directory_slots: consumer.directory_slots,
            directory_capacity: consumer.directory_capacity,
            directory_payload_bytes: consumer.directory_payload_bytes,
            packed_key_bytes: consumer.packed_key_bytes,
            packed_key_capacity: consumer.packed_key_capacity,
            observe_elapsed: self.observe_elapsed,
            observe_lock_wait: self.observe_lock_wait,
            consumer_busy: consumer.busy,
            consumer_wait: consumer.wait,
            consumer_cpu_ns: consumer.cpu_ns,
            drain_elapsed,
            statistics_elapsed,
            release_elapsed: release_elapsed + consumer.release,
            finish_elapsed: started.elapsed(),
        }
    }

    fn stop(&self) {
        if let Some(shared) = &self.shared {
            shared.catalog.lock().unwrap().stop = true;
            shared.changed.notify_one();
        }
    }
}

impl Drop for IdentityProbe {
    fn drop(&mut self) {
        if let Some(consumer) = self.consumer.take() {
            self.stop();
            let _ = consumer.join();
        }
    }
}

fn source_inline_bytes(states: &Vec<State>) -> usize {
    2 * size_of::<usize>() + size_of::<Vec<State>>() + states.capacity() * size_of::<State>()
}

fn state_heap_bytes(state: &State) -> usize {
    state.d.capacity() * size_of::<u32>() + state.ch.capacity() * size_of::<(u32, u32)>()
}

fn source_payload_bytes(states: &Vec<State>) -> usize {
    source_inline_bytes(states) + states.iter().map(state_heap_bytes).sum::<usize>()
}

fn retained_peak(lifetimes: &[Lifetime]) -> usize {
    let mut events = Vec::with_capacity(lifetimes.len() * 2);
    for lifetime in lifetimes {
        events.push((lifetime.observed, true, lifetime.bytes));
        events.push((lifetime.released, false, lifetime.bytes));
    }
    events.sort_unstable_by_key(|event| event.0);
    let (mut held, mut peak) = (0, 0);
    for (_, acquired, bytes) in events {
        if acquired {
            held += bytes;
            peak = peak.max(held);
        } else {
            held -= bytes;
        }
    }
    peak
}

fn consume(shared: &Shared, mode: Mode) -> ConsumerReport {
    let cpu_started = super::measure_memory::current_thread_cpu_ns();
    let mut report = match mode {
        Mode::Consume => consume_into::<ExpandedDirectory>(shared),
        Mode::Packed => consume_into::<PackedDirectory>(shared),
        _ => unreachable!(),
    };
    report.cpu_ns = super::measure_memory::current_thread_cpu_ns()
        .zip(cpu_started)
        .map(|(end, start)| end.saturating_sub(start));
    report
}

fn consume_into<D: Directory>(shared: &Shared) -> ConsumerReport {
    let mut report = ConsumerReport::default();
    let mut directory = D::default();
    loop {
        let waiting = Instant::now();
        let batch = {
            let mut catalog = shared.catalog.lock().unwrap();
            while catalog.batches.is_empty() && !catalog.stop {
                catalog = shared.changed.wait(catalog).unwrap();
            }
            catalog.batches.pop_front()
        };
        report.wait += waiting.elapsed();
        let Some(batch) = batch else {
            break;
        };
        let started = Instant::now();
        let mut source_bytes = source_inline_bytes(&batch.states);
        for state in batch.states.iter() {
            source_bytes += state_heap_bytes(state);
            report.scanned_states += 1;
            if !directory.insert(state) {
                report.duplicate_ids += 1;
            }
        }
        drop(batch.states);
        let released = {
            let mut catalog = shared.catalog.lock().unwrap();
            catalog.event += 1;
            catalog.event
        };
        report.lifetimes.push(Lifetime {
            observed: batch.observed,
            released,
            bytes: source_bytes,
        });
        report.busy += started.elapsed();
    }
    // Keep the actual identity contents observable, not just their lengths.
    directory.observe();
    let summary = directory.summary();
    report.directory_slots = summary.slots;
    report.directory_capacity = summary.capacity;
    report.directory_payload_bytes = summary.payload_bytes;
    report.packed_key_bytes = summary.packed_key_bytes;
    report.packed_key_capacity = summary.packed_key_capacity;
    let releasing = Instant::now();
    drop(directory);
    report.release = releasing.elapsed();
    report
}

#[test]
fn consumer_handles_sparse_ids_and_releases_observed_batches() {
    for mode in [Mode::Consume, Mode::Packed] {
        let states = Arc::new(vec![
            State {
                g: 1,
                d: vec![0],
                ch: vec![(1, 1)],
                rank: (1, 1),
                index: 2,
            },
            State {
                g: 2,
                d: vec![1, 2, 2],
                ch: vec![(1, 1)],
                rank: (1, 5),
                index: 7,
            },
        ]);
        let weak = Arc::downgrade(&states);
        let run = Run::start(mode);
        observe_states(&states);
        drop(states);
        let report = run.finish();
        report.print();
        assert!(weak.upgrade().is_none());
        assert_eq!(report.observed_states, 2);
        assert_eq!(report.scanned_states, 2);
        assert_eq!(report.directory_slots, 8);
        assert_eq!(report.duplicate_ids, 0);
        assert!(report.peak_retained_payload_bytes > 0);
        assert!(report.directory_payload_bytes > 0);
        if mode == Mode::Packed {
            assert!(report.packed_key_bytes > 0);
        }
    }
}

#[test]
fn packed_identity_roundtrip_preserves_zero_markings_and_all_c_powers() {
    fn take_u32(bytes: &mut &[u8]) -> u32 {
        let mut value = 0;
        for shift in (0..=28).step_by(7) {
            let (byte, rest) = bytes.split_first().expect("incomplete identity field");
            *bytes = rest;
            assert!(shift < 28 || *byte <= 15);
            value |= u32::from(byte & 127) << shift;
            if byte & 128 == 0 {
                return value;
            }
        }
        panic!("unterminated identity field");
    }

    fn decode(mut bytes: &[u8]) -> (u32, Vec<u32>, Vec<(u32, u32)>) {
        let g = take_u32(&mut bytes);
        let n = take_u32(&mut bytes);
        let d = (0..n).map(|_| take_u32(&mut bytes)).collect();
        let count = take_u32(&mut bytes);
        let ch = (0..count)
            .map(|_| (take_u32(&mut bytes), take_u32(&mut bytes)))
            .collect();
        assert!(bytes.is_empty());
        (g, d, ch)
    }

    let first = State {
        g: 130,
        d: vec![0, 130, 16_384],
        ch: vec![(1, 255), (129, 65_536)],
        rank: (0, 0),
        index: 2,
    };
    let extra_zero = State {
        g: first.g,
        d: vec![0, 0, 130, 16_384],
        ch: first.ch.clone(),
        rank: first.rank,
        index: 7,
    };
    let same_math = State {
        g: first.g,
        d: first.d.clone(),
        ch: first.ch.clone(),
        rank: first.rank,
        index: 9,
    };
    let unmarked = State {
        g: 2,
        d: vec![],
        ch: vec![(1, 3)],
        rank: (3, 2),
        index: 14,
    };
    let mut directory = PackedDirectory::default();
    for state in [&first, &extra_zero, &same_math, &unmarked] {
        assert!(directory.insert(state));
        assert_eq!(
            decode(directory.record(state.index).unwrap()),
            (state.g, state.d.clone(), state.ch.clone())
        );
    }
    assert_ne!(
        directory.record(first.index),
        directory.record(extra_zero.index)
    );
    assert_eq!(
        directory.record(first.index),
        directory.record(same_math.index)
    );
    assert!(directory.record(0).is_none());
    let bytes = directory.keys.len();
    assert!(!directory.insert(&first));
    assert_eq!(directory.keys.len(), bytes);
}
