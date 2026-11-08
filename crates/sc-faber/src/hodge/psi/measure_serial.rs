//! Test-only comparisons of the original and production serial memo sources.
//!
//! `owned` stores values solely in the shared copy of the original BTree shape;
//! it does not mirror them into a second result table. Lookup and insertion
//! guards end before recursive mathematics runs. `source` uses the production
//! SerialSource and completion directory; `pull` also starts the single C1
//! reader. All cases use private mathematical entries, without a RootSource.

use super::{DESCENDENTS, Descendents, integrate_mixed_with_threads, integrate_with_threads};
use crate::cache_source::{Reader, Registry, with_optional_registry};
use sc_num::Q;
use std::collections::BTreeMap;
use std::hint::black_box;
use std::sync::{Arc, RwLock};
use std::time::Instant;

pub(super) type SerialValues = BTreeMap<Vec<(u32, u32)>, BTreeMap<Vec<u32>, Q>>;

impl Descendents {
    #[inline]
    pub(super) fn measurement_memo_get(&self, d: &[u32], ch: &[(u32, u32)]) -> Option<Q> {
        if let Some(shared) = &self.shared_values {
            let values = shared.read().unwrap();
            values.get(ch).and_then(|row| row.get(d)).cloned()
        } else {
            self.values.get(ch).and_then(|row| row.get(d)).cloned()
        }
    }

    #[inline]
    pub(super) fn measurement_memo_insert(&mut self, d: Vec<u32>, ch: &[(u32, u32)], value: &Q) {
        if let Some(shared) = &self.shared_values {
            shared
                .write()
                .unwrap()
                .entry(ch.to_vec())
                .or_default()
                .insert(d, value.clone());
        } else {
            self.values
                .entry(ch.to_vec())
                .or_default()
                .insert(d, value.clone());
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Baseline,
    Owned,
    Source,
    Pull,
}

impl Mode {
    fn from_env() -> Self {
        match std::env::var("SC_FABER_MEASURE_SERIAL_SOURCE")
            .unwrap_or_else(|_| "baseline".into())
            .as_str()
        {
            "baseline" => Self::Baseline,
            "owned" => Self::Owned,
            "source" => Self::Source,
            "pull" => Self::Pull,
            value => panic!("unknown SC_FABER_MEASURE_SERIAL_SOURCE: {value}"),
        }
    }
}

fn reset(mode: Mode) {
    DESCENDENTS.with(|cell| {
        *cell.borrow_mut() = Descendents {
            shared_values: match mode {
                Mode::Baseline | Mode::Source | Mode::Pull => None,
                Mode::Owned => Some(Arc::new(RwLock::new(SerialValues::new()))),
            },
            ..Descendents::default()
        };
    });
    // Reclaim a previous case's expired context during setup, rather than
    // attributing its release to the next case's cold mathematical call.
    super::CACHED_DESCENDENTS.with(|cell| {
        cell.borrow_mut()
            .retain(|_, context| context.registry.strong_count() != 0);
    });
}

fn setup(mode: Mode) -> (Option<Arc<Registry>>, Option<Reader>) {
    reset(mode);
    let registry = matches!(mode, Mode::Source | Mode::Pull).then(|| Arc::new(Registry::new()));
    let reader = (mode == Mode::Pull).then(|| registry.as_ref().unwrap().start_reader().unwrap());
    (registry, reader)
}

fn source_sizes(registry: Option<&Arc<Registry>>) -> (usize, usize) {
    let Some(registry) = registry else {
        return (0, 0);
    };
    let sources = registry.sources();
    let rows = sources.iter().map(|source| source.ready()).sum();
    (sources.len(), rows)
}

fn finish_reader(mode: Mode, case: &str, reader: Option<Reader>) {
    if let Some(reader) = reader {
        let started = Instant::now();
        let report = reader.finish();
        let close = started.elapsed();
        assert!(report.reached_target);
        assert!(!report.invalid_cursor);
        assert!(!report.worker_panicked);
        eprintln!(
            "serial_source_reader mode={mode:?} case={case} finish_ns={} report={report:?}",
            close.as_nanos(),
        );
    }
}

fn memo_sizes() -> (usize, usize) {
    DESCENDENTS.with(|cell| {
        let evaluator = cell.borrow();
        let original = evaluator.values.values().map(BTreeMap::len).sum();
        let shared = evaluator.shared_values.as_ref().map_or(0, |values| {
            values.read().unwrap().values().map(BTreeMap::len).sum()
        });
        (original, shared)
    })
}

struct Case {
    name: &'static str,
    g: u32,
    n: u32,
    kappa: &'static [(u32, u32)],
    psi: &'static [(u32, u32)],
    lambda: &'static [(u32, u32)],
    ch: &'static [(u32, u32)],
}

const CASES: &[Case] = &[
    Case {
        name: "psi-g5-n3",
        g: 5,
        n: 3,
        kappa: &[],
        psi: &[(1, 3), (2, 4), (3, 8)],
        lambda: &[],
        ch: &[],
    },
    Case {
        name: "lambda1-g5",
        g: 5,
        n: 0,
        kappa: &[],
        psi: &[],
        lambda: &[(1, 12)],
        ch: &[],
    },
    Case {
        name: "kappa-psi-lambda-ch-g4",
        g: 4,
        n: 1,
        kappa: &[(1, 2)],
        psi: &[(1, 2)],
        lambda: &[(1, 3)],
        ch: &[(3, 1)],
    },
];

fn calculate(case: &Case) -> Q {
    integrate_mixed_with_threads(
        case.g,
        case.n,
        case.kappa,
        case.psi,
        case.lambda,
        case.ch,
        1,
    )
}

#[test]
fn shared_serial_memo_matches_original_and_is_the_only_value_source() {
    with_optional_registry(None, || {
        for case in CASES {
            reset(Mode::Baseline);
            let expected = calculate(case);
            let (original, shared) = memo_sizes();
            assert!(original > 0);
            assert_eq!(shared, 0);
            reset(Mode::Owned);
            assert_eq!(calculate(case), expected);
            assert_eq!(calculate(case), expected);
            let (original, shared) = memo_sizes();
            assert_eq!(original, 0);
            assert!(shared > 0);
        }
    });
    reset(Mode::Baseline);
}

#[test]
#[ignore = "manual serial source/C1 pull measurement; no SQL or persistence"]
fn measure_serial_source() {
    let warm_repeats: usize = std::env::var("SC_FABER_MEASURE_SERIAL_REPEATS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(10_000);
    let mode = Mode::from_env();
    for case in CASES {
        let setup = Instant::now();
        let (registry, reader) = self::setup(mode);
        let setup = setup.elapsed();
        let mathematical = Instant::now();
        let (value, cold, warm) = with_optional_registry(registry.clone(), || {
            let started = Instant::now();
            let value = black_box(calculate(black_box(case)));
            let cold = started.elapsed();
            let started = Instant::now();
            let mut last = value.clone();
            for _ in 0..warm_repeats {
                last = black_box(calculate(black_box(case)));
            }
            let warm = started.elapsed();
            assert_eq!(last, value);
            (value, cold, warm)
        });
        let mathematical = mathematical.elapsed();
        finish_reader(mode, case.name, reader);
        let (original_entries, shared_entries) = memo_sizes();
        let (sources, actual_source_rows) = source_sizes(registry.as_ref());
        if matches!(mode, Mode::Source | Mode::Pull) {
            assert_eq!(sources, 1);
            assert!(actual_source_rows > 0);
            assert_eq!(original_entries, 0);
            assert_eq!(shared_entries, 0);
        }
        eprintln!(
            "serial_source mode={mode:?} case={} threads=1 reset_ns={} cold_ns={} warm_repeats={warm_repeats} warm_batch_ns={} math_return_ns={} original_entries={original_entries} owned_prototype_entries={shared_entries} sources={sources} actual_source_rows={actual_source_rows} value={value}",
            case.name,
            setup.as_nanos(),
            cold.as_nanos(),
            warm.as_nanos(),
            mathematical.as_nanos(),
        );
    }

    measure_formula(mode, "lambda-g-theorem", || {
        integrate_with_threads(4, 2, &[(1, 3), (2, 4)], &[(4, 1)], 1)
    });
    measure_formula(mode, "one-point-psi", || {
        integrate_with_threads(6, 1, &[(1, 16)], &[], 1)
    });
    reset(Mode::Baseline);
}

/// Same-process diagnostic for small hot-path differences. Every measured
/// invocation still enters the real private integral and acquires its own
/// evaluator/memo access. No evaluator or memo guard spans the repeat loop.
/// The prewarmed original/source tables coexist; this is not a production
/// cold-query acceptance measurement and never exercises the public root L1.
#[test]
#[ignore = "manual interleaved hot serial baseline/source/C1 diagnostic; no SQL"]
fn measure_serial_source_paired() {
    let warm_repeats: usize = match std::env::var("SC_FABER_MEASURE_SERIAL_REPEATS") {
        Ok(value) => value
            .parse()
            .expect("serial repeats must be a positive integer"),
        Err(std::env::VarError::NotPresent) => 200_000,
        Err(error) => panic!("cannot read SC_FABER_MEASURE_SERIAL_REPEATS: {error}"),
    };
    assert!(warm_repeats > 0, "serial repeats must be positive");
    let rounds: usize = match std::env::var("SC_FABER_MEASURE_SERIAL_ROUNDS") {
        Ok(value) => value
            .parse()
            .expect("serial rounds must be a positive integer"),
        Err(std::env::VarError::NotPresent) => 9,
        Err(error) => panic!("cannot read SC_FABER_MEASURE_SERIAL_ROUNDS: {error}"),
    };
    assert!(rounds > 0, "serial rounds must be positive");
    let selected = match std::env::var("SC_FABER_MEASURE_SERIAL_CASE") {
        Ok(value) => Some(value),
        Err(std::env::VarError::NotPresent) => None,
        Err(error) => panic!("cannot read SC_FABER_MEASURE_SERIAL_CASE: {error}"),
    };
    if let Some(name) = &selected {
        assert!(
            CASES.iter().any(|case| case.name == name.as_str()),
            "unknown SC_FABER_MEASURE_SERIAL_CASE: {name}"
        );
    }
    let orders = [
        [Mode::Baseline, Mode::Source, Mode::Pull],
        [Mode::Pull, Mode::Baseline, Mode::Source],
        [Mode::Source, Mode::Pull, Mode::Baseline],
    ];

    for case in CASES
        .iter()
        .filter(|case| selected.as_deref().is_none_or(|name| case.name == name))
    {
        reset(Mode::Baseline);
        let registry = Arc::new(Registry::new());
        let expected = with_optional_registry(None, || black_box(calculate(black_box(case))));
        let source_value = with_optional_registry(Some(Arc::clone(&registry)), || {
            black_box(calculate(black_box(case)))
        });
        assert_eq!(source_value, expected, "{} source warmup", case.name);
        let original_sizes = memo_sizes();
        assert!(original_sizes.0 > 0);
        assert_eq!(original_sizes.1, 0);
        let targets: Vec<_> = registry
            .sources()
            .iter()
            .map(|source| source.ready())
            .collect();
        let ready_rows: usize = targets.iter().sum();
        assert!(!targets.is_empty());
        assert!(ready_rows > 0);
        eprintln!(
            "serial_source_paired_setup case={} warm_repeats={warm_repeats} rounds={rounds} original_entries={} sources={} source_rows={ready_rows} value={expected} note=prewarmed-private-entries-no-root-cache",
            case.name,
            original_sizes.0,
            targets.len(),
        );

        for round in 0..rounds {
            for (position, mode) in orders[round % orders.len()].into_iter().enumerate() {
                let round = round + 1;
                let position = position + 1;
                let before: Vec<_> = registry
                    .sources()
                    .iter()
                    .map(|source| source.ready())
                    .collect();
                assert_eq!(before, targets, "{} before round {round}", case.name);
                assert_eq!(memo_sizes(), original_sizes);
                let startup = Instant::now();
                let reader = (mode == Mode::Pull)
                    .then(|| registry.start_reader().expect("start paired C1 reader"));
                let startup = startup.elapsed();
                let owner = match mode {
                    Mode::Baseline => None,
                    Mode::Source | Mode::Pull => Some(Arc::clone(&registry)),
                    Mode::Owned => unreachable!(),
                };
                // Bind once outside timing. Each calculate still obtains its
                // own production context and memo access, exactly as before.
                let (last, elapsed) = with_optional_registry(owner, || {
                    let mut last = expected.clone();
                    let started = Instant::now();
                    for _ in 0..warm_repeats {
                        last = black_box(calculate(black_box(case)));
                    }
                    (last, started.elapsed())
                });
                assert_eq!(last, expected, "{} round {round} {mode:?}", case.name);
                if let Some(reader) = reader {
                    let finish = Instant::now();
                    let report = reader.finish();
                    let finish = finish.elapsed();
                    assert!(report.reached_target);
                    assert!(!report.invalid_cursor);
                    assert!(!report.worker_panicked);
                    assert_eq!(report.rows, ready_rows);
                    eprintln!(
                        "serial_source_paired_reader case={} round={round} position={position} start_ns={} finish_ns={} report={report:?} note=reader-copies-prewarmed-rows-from-zero",
                        case.name,
                        startup.as_nanos(),
                        finish.as_nanos(),
                    );
                }
                let after: Vec<_> = registry
                    .sources()
                    .iter()
                    .map(|source| source.ready())
                    .collect();
                assert_eq!(after, targets, "{} after round {round}", case.name);
                assert_eq!(memo_sizes(), original_sizes);
                eprintln!(
                    "serial_source_paired case={} round={round} position={position} mode={mode:?} threads=1 warm_repeats={warm_repeats} warm_batch_ns={} warm_mean_ns={:.3} source_rows={ready_rows} value={last}",
                    case.name,
                    elapsed.as_nanos(),
                    elapsed.as_nanos() as f64 / warm_repeats as f64,
                );
            }
        }
    }
    reset(Mode::Baseline);
}

fn measure_formula(mode: Mode, case: &str, calculate: impl FnOnce() -> Q) {
    let (registry, reader) = setup(mode);
    let started = Instant::now();
    let value = with_optional_registry(registry.clone(), || black_box(calculate()));
    let mathematical = started.elapsed();
    finish_reader(mode, case, reader);
    let (original_entries, shared_entries) = memo_sizes();
    let (sources, actual_source_rows) = source_sizes(registry.as_ref());
    eprintln!(
        "serial_source_formula mode={mode:?} case={case} math_return_ns={} original_entries={original_entries} owned_prototype_entries={shared_entries} sources={sources} actual_source_rows={actual_source_rows} value={value} note=private-entry-lambda-g-is-covered-by-public-root-source",
        mathematical.as_nanos(),
    );
}
