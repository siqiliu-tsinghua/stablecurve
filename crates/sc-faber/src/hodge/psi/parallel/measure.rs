//! Opt-in, test-only experiments. No persistence or production execution path.

use super::*;
use std::sync::atomic::{AtomicBool, AtomicU64};
use std::time::{Duration, Instant};

pub(super) struct Run {
    sampler: Option<super::measure_memory::Sampler>,
    identity: Option<super::measure_identity::Run>,
}

impl Run {
    pub(super) fn start() -> Self {
        let sampler = std::env::var_os("SC_FABER_MEASURE_MEMORY").map(|_| {
            let sampler = super::measure_memory::Sampler::start(Duration::from_millis(10)).unwrap();
            sampler.set_phase(super::measure_memory::Phase::Planning);
            println!("memory {:?}", sampler.sample("query_start"));
            sampler
        });
        let identity = match std::env::var("SC_FABER_MEASURE_IDENTITY").as_deref() {
            Ok("retain") => Some(super::measure_identity::Run::start(
                super::measure_identity::Mode::Retain,
            )),
            Ok("consume") => Some(super::measure_identity::Run::start(
                super::measure_identity::Mode::Consume,
            )),
            Ok("packed") => Some(super::measure_identity::Run::start(
                super::measure_identity::Mode::Packed,
            )),
            _ => None,
        };
        Self { sampler, identity }
    }

    pub(super) fn planned(&self, engine: &Engine, plan: &Plan) {
        if let Some(sampler) = &self.sampler {
            println!("memory {:?}", sampler.sample("plan_ready"));
            inspect_plan(engine, plan);
            sampler.set_phase(super::measure_memory::Phase::Evaluation);
        }
    }

    pub(super) fn evaluated(&self, engine: &Engine) {
        if let Some(sampler) = &self.sampler {
            println!("memory {:?}", sampler.sample("values_ready"));
            inspect_values(engine);
            sampler.set_phase(super::measure_memory::Phase::Release);
        }
    }

    pub(super) fn resized(&self) {
        if let Some(sampler) = &self.sampler {
            println!("memory {:?}", sampler.sample("value_slots_allocated"));
        }
    }

    pub(super) fn terms_released(&self) {
        if let Some(sampler) = &self.sampler {
            println!("memory {:?}", sampler.sample("terms_released"));
        }
    }

    pub(super) fn finish(&mut self, engine: &Engine) {
        if let Some(identity) = self.identity.take() {
            let report = identity.finish();
            if matches!(
                report.mode,
                super::measure_identity::Mode::Consume | super::measure_identity::Mode::Packed
            ) {
                assert_eq!(report.scanned_states, report.observed_states);
                assert_eq!(report.duplicate_ids, 0);
            }
            report.print();
        }
        if let Some(sampler) = self.sampler.take() {
            println!("memory {:?}", sampler.sample("plan_released"));
            let report = sampler.stop();
            for phase in report.phase_summaries() {
                println!("memory_summary {phase:?}");
            }
            println!(
                "memory_samples count={} worker_panicked={} residual_rows={}",
                report.samples.len(),
                report.worker_panicked,
                engine.values.len(),
            );
        }
    }
}

pub(super) fn check_memory_limit() {
    let limit_mib = env_usize("SC_FABER_MEASURE_LIMIT_MIB", 0);
    if limit_mib == 0 {
        return;
    }
    let snapshot = super::measure_memory::snapshot();
    if let Some(process) = snapshot.process {
        assert!(
            process.private_commit_bytes < limit_mib * 1024 * 1024,
            "controlled diagnostic stop: private commit reached {limit_mib} MiB"
        );
    }
    if let Some(system) = snapshot.system {
        assert!(
            system.available_physical_bytes > 4 * 1024 * 1024 * 1024,
            "controlled diagnostic stop: less than 4 GiB physical headroom"
        );
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Mode {
    Baseline,
    Owned,
    Pull,
}

impl Mode {
    fn from_env() -> Self {
        match std::env::var("SC_FABER_MEASURE_SOURCE").as_deref() {
            Ok("owned") => Self::Owned,
            Ok("pull") => Self::Pull,
            _ => Self::Baseline,
        }
    }
}

pub(super) fn heap_min(value: &Q) -> usize {
    let limb_size = value
        .denominator_ref()
        .limbs()
        .next()
        .map_or(8, |limb| std::mem::size_of_val(&limb));
    [value.numerator_ref(), value.denominator_ref()]
        .into_iter()
        .map(|natural| match natural.limb_count() as usize {
            0 | 1 => 0,
            count => count * limb_size,
        })
        .sum()
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

struct Source {
    values: RwLock<Vec<Q>>,
    layers: Vec<Vec<usize>>,
    ready: AtomicUsize,
    writing: AtomicBool,
    finished: AtomicBool,
}

#[derive(Default, Debug)]
struct PullReport {
    rows: usize,
    bytes_min: usize,
    windows: usize,
    read_ns: u128,
    max_read_ns: u128,
    retries: usize,
    cpu_ns: Option<u128>,
}

fn pull(source: &Source, base: usize, window_bytes: usize) -> PullReport {
    let cpu_start = super::measure_memory::current_thread_cpu_ns();
    let backoff_us = env_usize("SC_FABER_MEASURE_BACKOFF_US", 50);
    let backoff = || {
        if backoff_us == 0 {
            std::thread::yield_now();
        } else {
            std::thread::sleep(Duration::from_micros(backoff_us as u64));
        }
    };
    let mut report = PullReport::default();
    let mut layer = 0;
    let mut offset = 0;
    let mut copied = Vec::new();
    loop {
        let ready = source.ready.load(Ordering::Acquire);
        if layer == ready {
            if source.finished.load(Ordering::Acquire)
                && layer == source.ready.load(Ordering::Acquire)
            {
                break;
            }
            std::thread::sleep(Duration::from_millis(1));
            continue;
        }
        if source.writing.load(Ordering::Acquire) {
            report.retries += 1;
            backoff();
            continue;
        }
        let Ok(values) = source.values.try_read() else {
            report.retries += 1;
            backoff();
            continue;
        };
        let start = Instant::now();
        let mut bytes = 0;
        while layer < ready {
            let ids = &source.layers[layer];
            while offset < ids.len() {
                let value = &values[base + ids[offset]];
                bytes += std::mem::size_of::<Q>() + heap_min(value);
                copied.push(value.clone());
                offset += 1;
                if bytes >= window_bytes {
                    break;
                }
            }
            if offset == ids.len() {
                layer += 1;
                offset = 0;
            }
            if bytes >= window_bytes {
                break;
            }
        }
        drop(values);
        let elapsed = start.elapsed().as_nanos();
        report.rows += copied.len();
        report.bytes_min += bytes;
        report.windows += 1;
        report.read_ns += elapsed;
        report.max_read_ns = report.max_read_ns.max(elapsed);
        // Deliberately outside the source access: C1 measures copying, not SQL.
        std::hint::black_box(copied.as_slice());
        copied.clear();
    }
    report.cpu_ns = cpu_start.and_then(|start| {
        super::measure_memory::current_thread_cpu_ns().map(|end| end.saturating_sub(start))
    });
    report
}

fn write_values<'a>(
    source: &'a Source,
    wait_ns: &AtomicU64,
    details: bool,
) -> std::sync::RwLockWriteGuard<'a, Vec<Q>> {
    source.writing.store(true, Ordering::Release);
    let start = details.then(Instant::now);
    let guard = source.values.write().unwrap();
    if let Some(start) = start {
        wait_ns.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
    guard
}

fn publish(source: &Source, layer: usize) {
    source.ready.store(layer + 1, Ordering::Release);
    source.writing.store(false, Ordering::Release);
}

fn evaluate_owned(
    plan: &mut Plan,
    values: &mut Vec<Q>,
    base: usize,
    requested_threads: usize,
    mode: Mode,
) {
    let start = Instant::now();
    let source = Arc::new(Source {
        values: RwLock::new(std::mem::take(values)),
        layers: std::mem::take(&mut plan.layers),
        ready: AtomicUsize::new(0),
        writing: AtomicBool::new(false),
        finished: AtomicBool::new(false),
    });
    let window_bytes = env_usize("SC_FABER_MEASURE_WINDOW_KIB", 64).max(1) * 1024;
    let reader = (mode == Mode::Pull).then(|| {
        let source = Arc::clone(&source);
        std::thread::spawn(move || pull(&source, base, window_bytes))
    });
    let details = std::env::var_os("SC_FABER_MEASURE_DETAILS").is_some();
    let wait_ns = AtomicU64::new(0);
    let layers = &source.layers;
    let width = layers.iter().map(Vec::len).max().unwrap_or(0);
    let threads = if width < 128 {
        1
    } else {
        requested_threads.min(width.div_ceil(32)).max(1)
    };
    if threads == 1 {
        for (layer_index, layer) in layers.iter().enumerate() {
            let mut values = write_values(&source, &wait_ns, details);
            for &id in layer {
                values[base + id] = plan.evaluate(id, &values);
            }
            drop(values);
            publish(&source, layer_index);
        }
    } else {
        let (jobs_tx, jobs_rx) = mpsc::channel::<Job>();
        let jobs_rx = Mutex::new(jobs_rx);
        let (reply_tx, reply_rx) = mpsc::channel();
        let calculate = |job: Job| {
            let values = source.values.read().unwrap();
            layers[job.layer][job.start..job.end]
                .iter()
                .map(|&id| (id, plan.evaluate(id, &values)))
                .collect::<Vec<_>>()
        };
        let outcome = std::thread::scope(|scope| {
            for _ in 1..threads {
                let reply_tx = reply_tx.clone();
                let jobs_rx = &jobs_rx;
                let calculate = &calculate;
                scope.spawn(move || {
                    loop {
                        let Ok(job) = jobs_rx.lock().unwrap().recv() else {
                            break;
                        };
                        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            calculate(job)
                        }));
                        if reply_tx.send(result).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(reply_tx);
            let result = (|| {
                for (layer_index, layer) in layers.iter().enumerate() {
                    if layer.len() < 128 {
                        let mut values = write_values(&source, &wait_ns, details);
                        for &id in layer {
                            let result = plan.evaluate(id, &values);
                            values[base + id] = result;
                        }
                        drop(values);
                        publish(&source, layer_index);
                        continue;
                    }
                    let chunk = layer.len().div_ceil(threads * 4).max(32);
                    let mut jobs = 0;
                    for start in (chunk..layer.len()).step_by(chunk) {
                        jobs_tx
                            .send(Job {
                                layer: layer_index,
                                start,
                                end: (start + chunk).min(layer.len()),
                            })
                            .unwrap();
                        jobs += 1;
                    }
                    let mut completed = calculate(Job {
                        layer: layer_index,
                        start: 0,
                        end: chunk.min(layer.len()),
                    });
                    for _ in 0..jobs {
                        completed.extend(reply_rx.recv().unwrap()?);
                    }
                    let mut values = write_values(&source, &wait_ns, details);
                    for (id, value) in completed {
                        values[base + id] = value;
                    }
                    drop(values);
                    publish(&source, layer_index);
                }
                Ok(())
            })();
            drop(jobs_tx);
            result
        });
        if let Err(payload) = outcome {
            source.finished.store(true, Ordering::Release);
            // The experiment runs in an isolated child; do not hide a math panic.
            std::panic::resume_unwind(payload);
        }
    }
    let math_seconds = start.elapsed().as_secs_f64();
    source.finished.store(true, Ordering::Release);
    let report = reader.map(|reader| reader.join().unwrap());
    if let Some(report) = &report {
        assert_eq!(
            report.rows,
            layers.iter().map(Vec::len).sum::<usize>(),
            "background experiment missed completed results",
        );
    }
    let finish_seconds = start.elapsed().as_secs_f64();
    // The harness restores the original Engine for oracle comparisons. This is
    // not a production ownership design, nor an implicit persistence barrier.
    let source = match Arc::try_unwrap(source) {
        Ok(source) => source,
        Err(_) => panic!("measurement retained a result source"),
    };
    *values = source.values.into_inner().unwrap();
    plan.layers = source.layers;
    println!(
        "source mode={mode:?} math_seconds={math_seconds:.6} copy_finished_seconds={finish_seconds:.6} write_wait_seconds={:?} window_kib={} report={report:?}",
        details.then(|| wait_ns.load(Ordering::Relaxed) as f64 / 1e9),
        window_bytes / 1024,
    );
}

pub(super) fn evaluate(plan: &mut Plan, values: &mut Vec<Q>, base: usize, threads: usize) {
    match Mode::from_env() {
        Mode::Baseline => super::evaluate_layers(plan, values, base, threads),
        mode => evaluate_owned(plan, values, base, threads, mode),
    }
}

pub(super) fn reevaluate_same_plan(
    plan: &mut Plan,
    values: &mut Vec<Q>,
    base: usize,
    threads: usize,
) {
    if std::env::var_os("SC_FABER_MEASURE_REEVALUATE").is_none() {
        return;
    }
    // Reuse mathematical identities, equations and their exact physical layout;
    // recompute every new result, rather than time a memo hit or clone a plan.
    for (round, mode) in [
        Mode::Baseline,
        Mode::Owned,
        Mode::Pull,
        Mode::Pull,
        Mode::Baseline,
        Mode::Owned,
        Mode::Owned,
        Mode::Pull,
        Mode::Baseline,
    ]
    .into_iter()
    .enumerate()
    {
        for value in &mut values[base..] {
            *value = Q::ZERO;
        }
        let started = Instant::now();
        match mode {
            Mode::Baseline => super::evaluate_layers(plan, values, base, threads),
            mode => evaluate_owned(plan, values, base, threads, mode),
        }
        println!(
            "same_plan round={} mode={mode:?} evaluation_seconds={:.6}",
            round + 1,
            started.elapsed().as_secs_f64(),
        );
    }
}

fn print_values(label: &str, values: &[Q], capacity: usize) {
    let heap_min: usize = values.iter().map(heap_min).sum();
    println!(
        "layout {label} rows={} slot_capacity_bytes={} rational_heap_min_bytes={heap_min}",
        values.len(),
        capacity * std::mem::size_of::<Q>(),
    );
}

pub(super) fn inspect_plan(engine: &Engine, plan: &Plan) {
    let compact: usize = plan
        .terms
        .iter()
        .map(|block| block.compact.capacity())
        .sum();
    let wide: usize = plan.terms.iter().map(|block| block.wide.capacity()).sum();
    let layer_slots: usize = plan.layers.iter().map(Vec::capacity).sum();
    let largest_block_bytes = plan
        .terms
        .iter()
        .map(|block| {
            block.compact.capacity() * std::mem::size_of::<CompactTerm>()
                + block.wide.capacity() * std::mem::size_of::<Term>()
        })
        .max()
        .unwrap_or(0);
    let largest_equation_terms = plan
        .equations
        .iter()
        .map(|equation| equation.terms.len())
        .max()
        .unwrap_or(0);
    let wide_heap_min: usize = plan
        .terms
        .iter()
        .flat_map(|block| &block.wide)
        .map(|term| match &term.factor {
            Coefficient::Rational(value) => std::mem::size_of::<Q>() + heap_min(value),
            _ => 0,
        })
        .sum();
    let equation_heap_min: usize = plan
        .equations
        .iter()
        .map(|equation| {
            heap_min(&equation.constant) + equation.divisor.as_ref().map_or(0, heap_min)
        })
        .sum();
    let mut keys = 0;
    let mut key_bytes = 0;
    for row in &engine.index {
        for shard in row {
            let index = shard.read().unwrap();
            for (ch, entries) in index.iter() {
                key_bytes += ch.capacity() * std::mem::size_of::<(u32, u32)>();
                for d in entries.keys() {
                    keys += 1;
                    key_bytes += d.capacity() * std::mem::size_of::<u32>();
                }
            }
        }
    }
    println!(
        "layout plan equations={} equation_capacity_bytes={} compact_capacity_bytes={} wide_capacity_bytes={} layer_capacity_bytes={} wide_extra_heap_min_bytes={wide_heap_min} equation_heap_min_bytes={equation_heap_min} blocks={} largest_block_capacity_bytes={largest_block_bytes} largest_equation_terms={largest_equation_terms} index_keys={keys} index_key_payload_bytes={key_bytes} btree_nodes_not_included=true",
        plan.equations.len(),
        plan.equations.capacity() * std::mem::size_of::<PackedEquation>(),
        compact * std::mem::size_of::<CompactTerm>(),
        wide * std::mem::size_of::<Term>(),
        layer_slots * std::mem::size_of::<usize>(),
        plan.terms.len(),
    );
    print_values(
        "before_evaluation",
        &engine.values,
        engine.values.capacity(),
    );
}

pub(super) fn inspect_values(engine: &Engine) {
    print_values("after_evaluation", &engine.values, engine.values.capacity());
}

#[test]
#[ignore = "opt-in cache design experiment; isolated process, no database"]
fn cache_source_experiment() {
    let genus = env_usize("SC_FABER_MEASURE_GENUS", 14) as u32;
    let threads = env_usize("SC_FABER_THREADS", 24);
    let start = Instant::now();
    let value = super::super::integrate_with_threads(genus, 0, &[], &[(1, 3 * genus - 3)], threads);
    let elapsed = start.elapsed();
    let (planning, evaluation, equations, terms) = super::profile();
    println!(
        "query genus={genus} mode={:?} threads={threads} total_seconds={:.6} planning_seconds={:.6} evaluation_seconds={:.6} equations={equations} terms={terms} exact={value}",
        Mode::from_env(),
        elapsed.as_secs_f64(),
        planning.as_secs_f64(),
        evaluation.as_secs_f64(),
    );
    if std::env::var_os("SC_FABER_MEASURE_MEMORY").is_some() {
        ENGINE.with(|cell| {
            let mut engine = cell.borrow_mut();
            engine.values.clear();
            println!(
                "memory values_cleared_keep_slots {:?}",
                super::measure_memory::snapshot(),
            );
            engine.index.clear();
            println!(
                "memory index_cleared {:?}",
                super::measure_memory::snapshot(),
            );
            engine.values.shrink_to_fit();
            println!(
                "memory slots_released {:?}",
                super::measure_memory::snapshot(),
            );
        });
    }
}
