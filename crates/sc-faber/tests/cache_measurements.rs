//! Manual measurements through the production library's public entry points.
//!
//! Run one ignored test at a time in release mode, for example:
//! `cargo test -p sc-faber --release --test cache_measurements measure_lambda1_top -- --ignored --nocapture --test-threads=1`
//!
//! Each outer repetition starts a fresh process. SC_FABER_MEASURE_REPEATS
//! defaults to 1; SC_FABER_THREADS defaults to 24 and remains an upper bound.
//! SC_FABER_MEASURE_SOURCE_MODE defaults to `off`. `owned` enables result
//! sources without opening SQLite; `pull` adds the C1 background consumer.
//! `persist` uses a fresh temporary disk SQLite database unless
//! SC_FABER_MEASURE_CACHE_PATH supplies a filesystem path. Provided paths are
//! reported as fresh only when absent and are never deleted. Automatic SQLite/WAL/SHM files
//! are removed after a successful close; no other measurement files are made.
//! Public-call times exclude attachment and printing. `durable_total_s` spans
//! work (including its result printing) and flush (including its reporting).
//! Process times include the harness, attachment, flush, close, and cleanup.
//! Repeated small cases with sources enabled measure public root L1 hits,
//! rather than the cost of hitting the internal descendant memo.
//! `measure_lambda1_recovery` instead reserves its own database, interrupts
//! after committed nonzero C DAG values but before the root, then reports
//! restart and complete-root-hit math, flush, process and storage-read counts.
//! SC_FABER_MEASURE_RECOVERY_GENUS defaults to 15 and also accepts 12; both
//! workloads have known exact reference values. Genus 12 can finish before
//! the observer sees committed descendants, so it is only a short diagnostic.
//! `measure_lambda_top_triple` computes lambda_g lambda_(g-1) lambda_(g-2)
//! through public klp, including its algebraic shortcuts. It now measures the
//! closed-form path, not general-GRR memory pressure. A separate exact-value
//! check runs outside the timed call. Every supplied SQLite path must
//! be absent before launch, so this workload cannot become a cached root hit.
//! For general workloads, work_written_blocks and work_loaded_blocks confirm
//! actual spill/readback; a formula-only triple should not create work blocks.
//! SC_FABER_MEASURE_TIMEOUT_SECS optionally bounds each owned child, including
//! flush/close. On timeout the child is terminated and its database retained.

use rusqlite::{Connection, ErrorCode, OpenFlags, params};
use sc_faber::hodge::klp;
use sc_num::Q;
use std::fs::{self, OpenOptions};
use std::hint::black_box;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

const CHILD: &str = "SC_FABER_MEASURE_CHILD";
const PULL: &str = "SC_FABER_CACHE_MEASURE_PULL";
const SOURCE_ONLY: &str = "SC_FABER_CACHE_MEASURE_SOURCE_ONLY";
const REPORT: &str = "SC_FABER_CACHE_MEASURE_REPORT";
const CACHE_PATH: &str = "SC_FABER_MEASURE_CACHE_PATH";
const RECOVERY_PHASE: &str = "SC_FABER_MEASURE_RECOVERY_PHASE";
const RECOVERY_GENUS: &str = "SC_FABER_MEASURE_RECOVERY_GENUS";

fn source_mode() -> &'static str {
    match std::env::var("SC_FABER_MEASURE_SOURCE_MODE").as_deref() {
        Ok("off") | Err(std::env::VarError::NotPresent) => "off",
        Ok("owned") => "owned",
        Ok("pull") => "pull",
        Ok("persist") => "persist",
        Ok(value) => panic!("unknown SC_FABER_MEASURE_SOURCE_MODE: {value}"),
        Err(error) => panic!("cannot read SC_FABER_MEASURE_SOURCE_MODE: {error}"),
    }
}

struct DiskCache {
    path: PathBuf,
    automatic: bool,
}

impl DiskCache {
    fn prepare() -> Self {
        if let Some(path) = std::env::var_os(CACHE_PATH) {
            let path = PathBuf::from(path);
            assert!(
                !path.as_os_str().is_empty() && path != Path::new(":memory:"),
                "{CACHE_PATH} must name a disk database"
            );
            // Treat the option as a filesystem path, including relative names,
            // rather than permitting an SQLite URI to select memory storage.
            let path = if path.is_absolute() {
                path
            } else {
                std::env::current_dir()
                    .expect("measurement working directory")
                    .join(path)
            };
            return Self {
                path,
                automatic: false,
            };
        }
        Self::fresh()
    }

    fn fresh() -> Self {
        let directory = std::env::temp_dir();
        let process = std::process::id();
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("measurement clock before Unix epoch")
            .as_nanos();
        for attempt in 0..100 {
            let path = directory.join(format!(
                "sc-faber-c2-{process}-{timestamp}-{attempt}.sqlite"
            ));
            match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(file) => {
                    drop(file);
                    return Self {
                        path,
                        automatic: true,
                    };
                }
                Err(error) if error.kind() == ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("create fresh SQLite database {}: {error}", path.display()),
            }
        }
        panic!("could not reserve a unique temporary SQLite database");
    }

    fn remove_files(&self) -> std::io::Result<()> {
        if !self.automatic {
            return Ok(());
        }
        // Only these exact automatically reserved paths are eligible, and only
        // after close_cache or the owned child's wait has released all handles.
        for suffix in ["-wal", "-shm", ""] {
            let mut name = self.path.as_os_str().to_os_string();
            name.push(suffix);
            let path = PathBuf::from(name);
            match fs::remove_file(&path) {
                Ok(()) => {}
                Err(error) if error.kind() == ErrorKind::NotFound => {}
                Err(error) => return Err(error),
            }
        }
        Ok(())
    }

    fn cleanup(self) {
        self.remove_files().unwrap_or_else(|error| {
            panic!(
                "remove temporary SQLite files {}: {error}",
                self.path.display()
            )
        });
    }
}

fn setting(name: &str, default: u32) -> u32 {
    match std::env::var(name) {
        Ok(value) => value
            .parse()
            .ok()
            .filter(|&value| value > 0)
            .unwrap_or_else(|| panic!("{name} must be a positive integer")),
        Err(std::env::VarError::NotPresent) => default,
        Err(error) => panic!("cannot read {name}: {error}"),
    }
}

fn isolated(name: &str, threads: Option<u32>, work: impl FnOnce()) {
    if cfg!(debug_assertions) {
        panic!("run cache_measurements with cargo test --release");
    }
    if std::env::var(CHILD).is_ok_and(|child| child == name) {
        let mode = source_mode();
        if std::env::var_os("SC_FABER_MEASURE_MEMORY_MIB").is_some() {
            let mib = setting("SC_FABER_MEASURE_MEMORY_MIB", 1);
            sc_faber::socle::set_cache_memory_limit(Some(u64::from(mib) * 1024 * 1024))
                .expect("set process-wide measurement memory budget");
            println!("measurement_memory_limit_mib={mib}");
        }
        let disk_cache = (mode == "persist").then(DiskCache::prepare);
        if let Some(cache) = &disk_cache {
            println!(
                "measurement_cache workload={name} source_mode={mode} path={} cache_state={}",
                cache.path.display(),
                if cache.automatic || !cache.path.exists() {
                    "fresh"
                } else {
                    "provided_non_cold"
                }
            );
        }
        match mode {
            "off" => sc_faber::socle::close_cache().expect("start without a cache attachment"),
            "owned" | "pull" => sc_faber::socle::attach_cache(":memory:")
                .expect("attach diagnostic result source registry without SQLite"),
            "persist" => sc_faber::socle::attach_cache(
                disk_cache
                    .as_ref()
                    .unwrap()
                    .path
                    .to_str()
                    .expect("SQLite measurement path must be valid UTF-8"),
            )
            .expect("attach persistent disk SQLite cache"),
            _ => unreachable!(),
        }
        println!(
            "measurement workload={name} source_mode={mode} requested_threads={} effective_thread_limit={}",
            std::env::var("SC_FABER_THREADS").unwrap(),
            sc_num::parallel::thread_count("SC_FABER_THREADS", 24)
        );
        let work_start = Instant::now();
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(work));
        let flush_start = Instant::now();
        let flushed = sc_faber::socle::flush_cache();
        let flush_elapsed = flush_start.elapsed();
        let durable_total = work_start.elapsed();
        println!(
            "measurement_flush workload={name} source_mode={mode} elapsed_s={:.9} durable_total_s={:.9} persistent_storage={} includes_work_and_flush_reporting=true",
            flush_elapsed.as_secs_f64(),
            durable_total.as_secs_f64(),
            mode == "persist"
        );
        if let Ok(phase) = std::env::var(RECOVERY_PHASE) {
            println!(
                "measurement_recovery_flush phase={phase} elapsed_s={:.9} durable_total_s={:.9} includes_work_and_flush_reporting=true",
                flush_elapsed.as_secs_f64(),
                durable_total.as_secs_f64()
            );
        }
        let start = Instant::now();
        let closed = sc_faber::socle::close_cache();
        println!(
            "measurement_source_close workload={name} source_mode={mode} elapsed_s={:.9}",
            start.elapsed().as_secs_f64()
        );
        if let Err(payload) = outcome {
            if let Some(cache) = &disk_cache {
                eprintln!(
                    "measurement_failed_cache_retained path={} flush={flushed:?} close={closed:?}",
                    cache.path.display()
                );
            }
            std::panic::resume_unwind(payload);
        }
        flushed.expect("flush completed measurement results");
        closed.expect("close measurement cache normally");
        if let Some(cache) = disk_cache {
            cache.cleanup();
        }
        return;
    }

    let threads = threads.unwrap_or_else(|| setting("SC_FABER_THREADS", 24));
    let repetitions = setting("SC_FABER_MEASURE_REPEATS", 1);
    let executable = std::env::current_exe().expect("measurement executable path");
    let mode = source_mode();
    for run in 1..=repetitions {
        println!(
            "measurement_start workload={name} run={run} source_mode={mode} requested_threads={threads}",
        );
        let start = Instant::now();
        let mut command = Command::new(&executable);
        command
            .args([
                "--ignored",
                "--exact",
                name,
                "--nocapture",
                "--test-threads=1",
            ])
            .env(CHILD, name)
            .env("SC_FABER_MEASURE_SOURCE_MODE", mode)
            .env("SC_FABER_THREADS", threads.to_string());
        // Configure the isolated child's environment before its threads or
        // cache reader exist, without mutating the running test process.
        match mode {
            "owned" => {
                command
                    .env(SOURCE_ONLY, "1")
                    .env_remove(PULL)
                    .env_remove(REPORT);
            }
            "pull" => {
                command
                    .env(PULL, "1")
                    .env_remove(SOURCE_ONLY)
                    .env_remove(REPORT);
            }
            "off" => {
                command
                    .env_remove(SOURCE_ONLY)
                    .env_remove(PULL)
                    .env_remove(REPORT);
            }
            "persist" => {
                command
                    .env_remove(SOURCE_ONLY)
                    .env_remove(PULL)
                    .env(REPORT, "1");
            }
            _ => unreachable!(),
        }
        let status = if std::env::var_os("SC_FABER_MEASURE_TIMEOUT_SECS").is_some() {
            let timeout =
                Duration::from_secs(u64::from(setting("SC_FABER_MEASURE_TIMEOUT_SECS", 1)));
            let mut child = RecoveryChild {
                child: command.spawn().expect("start bounded measurement process"),
            };
            loop {
                if let Some(status) = child.child.try_wait().expect("poll bounded measurement") {
                    break status;
                }
                assert!(
                    start.elapsed() < timeout,
                    "measurement child exceeded {timeout:?}; terminating owned child and retaining its cache"
                );
                thread::sleep(Duration::from_millis(10));
            }
        } else {
            command.status().expect("start fresh measurement process")
        };
        println!(
            "measurement_process workload={name} run={run} elapsed_s={:.9}",
            start.elapsed().as_secs_f64()
        );
        assert!(status.success(), "measurement child failed: {status}");
    }
}

fn lambda1_top(g: u32) -> Q {
    assert!(g >= 2, "lambda1 measurements require genus at least 2");
    let exponent = g
        .checked_mul(3)
        .and_then(|value| value.checked_sub(3))
        .expect("lambda1 exponent exceeds u32");
    klp(g, 0, &[], &[], &[(1, exponent)], &[]).expect("lambda1 top integral")
}

fn check_known_top(g: u32, value: &Q) {
    let expected = match g {
        2 => "1/2880",
        3 => "1/90720",
        12 => "1787081855869384294007862200032/1331337837123669375",
        15 => "2738843642342159071177877514044195843584/41976078300646875",
        // Independent production run of the original v0.3.0 public entry.
        18 => {
            "20689758456866608474576165063672737895330310255093176113198592/248288903138376392489071875"
        }
        _ => return,
    };
    assert_eq!(value.to_string(), expected, "lambda1 genus {g}");
}

#[test]
#[ignore = "manual fresh-process lambda1 measurement through the production public entry"]
fn measure_lambda1_top() {
    isolated("measure_lambda1_top", None, || {
        let g = setting("SC_FABER_MEASURE_GENUS", 12);
        let start = Instant::now();
        let value = lambda1_top(g);
        let elapsed = start.elapsed();
        check_known_top(g, &value);
        println!(
            "measurement_query workload=lambda1_top genus={g} elapsed_s={:.9} value={value}",
            elapsed.as_secs_f64()
        );
        if let Ok(phase) = std::env::var(RECOVERY_PHASE) {
            println!(
                "measurement_recovery_math phase={phase} genus={g} elapsed_s={:.9} exact_golden_checked=true",
                elapsed.as_secs_f64()
            );
        }
    });
}

fn lambda_top_triple_reference(g: u32) -> Q {
    assert!(g >= 2, "top lambda triples require genus at least two");
    let twice_g = g.checked_mul(2).expect("Bernoulli index exceeds u32");
    // Faber--Pandharipande, Hodge integrals and Gromov--Witten theory,
    // Theorem 4 and Section 4.3: https://arxiv.org/abs/math/9810173.
    // Adjacent even Bernoulli numbers have opposite signs.
    -(sc_faber::socle::bernoulli(u64::from(twice_g))
        * sc_faber::socle::bernoulli(u64::from(twice_g - 2)))
        / (Q::from(2) * Q::from(twice_g) * Q::from(twice_g - 2) * sc_num::factorial(twice_g - 2))
}

#[test]
fn lambda_top_triple_reference_has_independent_small_genus_anchors() {
    for (g, expected) in [(2, "1/5760"), (3, "1/1451520"), (4, "1/87091200")] {
        assert_eq!(lambda_top_triple_reference(g).to_string(), expected);
    }
}

#[test]
#[ignore = "manual cold public lambda triple, including the closed-form shortcut"]
fn measure_lambda_top_triple() {
    if source_mode() == "persist" {
        if let Some(path) = std::env::var_os(CACHE_PATH) {
            assert!(
                !Path::new(&path).exists(),
                "hard lambda triple requires a fresh SQLite path; refusing an existing cache"
            );
        }
    }
    isolated("measure_lambda_top_triple", None, || {
        let g = setting("SC_FABER_MEASURE_GENUS", 12);
        assert!(g >= 2, "top lambda triples require genus at least two");
        // lambda_0 is the unit in genus two. No hodge_data formula/reduction
        // participates in this computation; mathematical shortcuts belong to klp.
        let lambda: Vec<_> = (g - 2..=g)
            .filter(|&index| index != 0)
            .map(|index| (index, 1))
            .collect();
        let start = Instant::now();
        let value = klp(g, 0, &[], &[], &lambda, &[]).expect("top lambda triple");
        let elapsed = start.elapsed();
        assert_eq!(
            value,
            lambda_top_triple_reference(g),
            "lambda triple genus {g}"
        );
        println!(
            "measurement_query workload=lambda_top_triple genus={g} elapsed_s={:.9} value={value} direct_klp=true exact_closed_formula_checked=true",
            elapsed.as_secs_f64()
        );
    });
}

// Only this owned handle is ever terminated, including when its parent test
// panics. Inherited output cannot fill an undrained child pipe during polling.
struct RecoveryChild {
    child: Child,
}

impl RecoveryChild {
    fn spawn(path: &Path, phase: &str, g: u32) -> Self {
        let child = Command::new(std::env::current_exe().expect("measurement executable path"))
            .args([
                "--ignored",
                "--exact",
                "measure_lambda1_top",
                "--nocapture",
                "--test-threads=1",
            ])
            .env(CHILD, "measure_lambda1_top")
            .env("SC_FABER_MEASURE_SOURCE_MODE", "persist")
            .env("SC_FABER_MEASURE_GENUS", g.to_string())
            .env("SC_FABER_THREADS", "24")
            .env(CACHE_PATH, path)
            .env(RECOVERY_PHASE, phase)
            .env(REPORT, "1")
            .env_remove(SOURCE_ONLY)
            .env_remove(PULL)
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("start owned recovery measurement child");
        Self { child }
    }

    fn stop(&mut self) -> bool {
        let sent_kill = match self.child.kill() {
            Ok(()) => true,
            Err(error) => {
                let status = self
                    .child
                    .try_wait()
                    .expect("check child after termination race")
                    .unwrap_or_else(|| panic!("terminate owned recovery child: {error}"));
                assert!(
                    status.success(),
                    "recovery child failed before kill: {status}"
                );
                false
            }
        };
        let status = self
            .child
            .wait()
            .expect("wait for terminated recovery child");
        sent_kill && !status.success()
    }
}

impl Drop for RecoveryChild {
    fn drop(&mut self) {
        if !matches!(self.child.try_wait(), Ok(Some(_))) {
            // Drop also runs on observer errors and deadlines. Never leave a
            // measurement process consuming memory after a failed measurement.
            let _ = self.child.kill();
            let _ = self.child.wait();
        }
    }
}

// Stable version-1 keys: tag 4 orders descendants by C-weight, so this range
// excludes weight-zero pure-psi results. In the genus-12/15 workloads all
// W/E ranks are below 128: [1,4,W,E,g,...] has genus at SQL byte position 5.
// Restrict to genus 2..g and exclude canonical Q zero (01 00 00 01 01), so early
// genus-zero zeros and genus-one closed formulas cannot trigger interruption.
// The root is Integral(g,n=0,kappa=[],psi=[],lambda=[(1,3g-3)],ch=[]), with
// one-byte ULEB128 fields. These byte offsets are not a general key decoder.
const RECOVERY_DESC_LOWER: &[u8] = &[1, 4, 1];
const RECOVERY_DESC_UPPER: &[u8] = &[1, 5];

fn recovery_observer(path: &Path) -> Connection {
    let connection = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .expect("open recovery database read-only");
    connection
        .busy_timeout(Duration::from_millis(1))
        .expect("bound read-only recovery observation waits");
    connection
}

fn recovery_state(connection: &Connection, g: u32) -> rusqlite::Result<(bool, bool)> {
    let genus = [u8::try_from(g).expect("one-byte recovery genus")];
    let root = [
        1,
        2,
        genus[0],
        0,
        0,
        0,
        1,
        1,
        u8::try_from(3 * g - 3).expect("one-byte recovery lambda exponent"),
        0,
    ];
    // One statement observes both conditions in the same committed snapshot.
    // No read transaction is held between polls, and both probes use the key.
    connection.query_row(
        "SELECT EXISTS(SELECT 1 FROM faber_integrals WHERE k >= ?1 AND k < ?2
                       AND substr(k,5,1) >= x'02' AND substr(k,5,1) <= ?4
                       AND v != x'0100000101'),
                EXISTS(SELECT 1 FROM faber_integrals WHERE k = ?3)",
        params![
            RECOVERY_DESC_LOWER,
            RECOVERY_DESC_UPPER,
            root.as_slice(),
            genus.as_slice()
        ],
        |row| Ok((row.get(0)?, row.get(1)?)),
    )
}

fn observer_busy(error: &rusqlite::Error) -> bool {
    matches!(
        error.sqlite_error_code(),
        Some(ErrorCode::DatabaseBusy | ErrorCode::DatabaseLocked)
    )
}

struct InterruptedCache {
    rows: u64,
    expensive_descendant_rows: u64,
    elapsed: Duration,
}

fn interrupt_before_root(path: &Path, g: u32) -> Option<InterruptedCache> {
    let start = Instant::now();
    let mut child = RecoveryChild::spawn(path, "interrupted", g);
    let observer = recovery_observer(path);
    let mut schema_ready = false;
    let deadline = Duration::from_secs(300);
    loop {
        assert!(
            start.elapsed() < deadline,
            "no committed nonzero genus >= 2 C-descendant before root within five minutes"
        );
        if let Some(status) = child.child.try_wait().expect("poll owned recovery child") {
            assert!(
                status.success(),
                "recovery child failed before interruption: {status}"
            );
            println!("measurement_recovery_race genus={g} reason=child_completed_before_kill");
            return None;
        }
        let observed = if schema_ready {
            recovery_state(&observer, g)
        } else {
            observer
                .query_row(
                    "SELECT EXISTS(SELECT 1 FROM sqlite_master
                     WHERE type = 'table' AND name = 'faber_integrals')",
                    [],
                    |row| row.get::<_, bool>(0),
                )
                .map(|ready| {
                    schema_ready = ready;
                    (false, false)
                })
        };
        match observed {
            Ok((_, true)) => {
                child.stop();
                println!("measurement_recovery_race genus={g} reason=root_committed_before_kill");
                return None;
            }
            Ok((true, false)) => {
                if !child.stop() {
                    println!(
                        "measurement_recovery_race genus={g} reason=child_completed_during_kill"
                    );
                    return None;
                }
                break;
            }
            Ok((false, false)) => {}
            Err(error) if observer_busy(&error) => {}
            Err(error) => panic!("observe recovery database: {error}"),
        }
        thread::sleep(Duration::from_millis(1));
    }
    drop(observer);

    // Reopen after the killed process has exited; a racing root commit makes
    // the attempt invalid even if the preceding observation lacked the root.
    let observer = recovery_observer(path);
    let (descendants, root) = recovery_state(&observer, g).expect("reopen interrupted database");
    if root {
        println!("measurement_recovery_race genus={g} reason=root_committed_during_kill");
        return None;
    }
    assert!(
        descendants,
        "committed nonzero genus >= 2 C-descendants vanished after interruption"
    );
    let rows: i64 = observer
        .query_row("SELECT count(*) FROM faber_integrals", [], |row| row.get(0))
        .expect("count committed rows after interruption");
    let expensive_descendant_rows: i64 = observer
        .query_row(
            "SELECT count(*) FROM faber_integrals WHERE k >= ?1 AND k < ?2
             AND substr(k,5,1) >= x'02' AND substr(k,5,1) <= ?3
             AND v != x'0100000101'",
            params![
                RECOVERY_DESC_LOWER,
                RECOVERY_DESC_UPPER,
                [u8::try_from(g).expect("one-byte recovery genus")].as_slice()
            ],
            |row| row.get(0),
        )
        .expect("count committed nonzero genus >= 2 C-descendants after interruption");
    Some(InterruptedCache {
        rows: rows.try_into().expect("nonnegative committed row count"),
        expensive_descendant_rows: expensive_descendant_rows
            .try_into()
            .expect("nonnegative committed C-descendant row count"),
        elapsed: start.elapsed(),
    })
}

fn finish_recovery_phase(path: &Path, phase: &str, g: u32) {
    println!(
        "measurement_recovery_start phase={phase} genus={g} requested_threads=24 storage_read_counts=StorageReport.read_keys/read_hits"
    );
    let start = Instant::now();
    let mut child = RecoveryChild::spawn(path, phase, g);
    let deadline = Duration::from_secs(900);
    let status = loop {
        if let Some(status) = child.child.try_wait().expect("poll recovery completion") {
            break status;
        }
        assert!(
            start.elapsed() < deadline,
            "{phase} child exceeded fifteen-minute recovery deadline"
        );
        thread::sleep(Duration::from_millis(20));
    };
    println!(
        "measurement_recovery_process phase={phase} genus={g} elapsed_s={:.9} includes_attach_math_flush_close=true",
        start.elapsed().as_secs_f64()
    );
    assert!(status.success(), "{phase} recovery child failed: {status}");
    assert!(
        recovery_state(&recovery_observer(path), g)
            .expect("verify complete root after recovery child")
            .1,
        "{phase} child exited without committing the complete root"
    );
}

// Declared before any child/observer in each attempt, so their guards have
// already released all handles before this attempts exact-path cleanup.
struct RecoveryFiles(DiskCache);

impl Drop for RecoveryFiles {
    fn drop(&mut self) {
        if let Err(error) = self.0.remove_files() {
            eprintln!(
                "recovery temporary SQLite cleanup failed path={}: {error}",
                self.0.path.display()
            );
        }
    }
}

#[test]
#[ignore = "manual committed-descendant interruption, restart and complete-root disk-hit measurement"]
fn measure_lambda1_recovery() {
    if cfg!(debug_assertions) {
        panic!("run cache_measurements with cargo test --release");
    }
    let g = setting(RECOVERY_GENUS, 15);
    assert!(
        matches!(g, 12 | 15),
        "{RECOVERY_GENUS} must be 12 or 15, with a known exact reference"
    );
    // Ignore caller-provided cache paths: termination and cleanup are limited
    // to the fresh files and child handles owned by this measurement.
    for attempt in 1..=3 {
        let cache = RecoveryFiles(DiskCache::fresh());
        println!(
            "measurement_recovery_attempt attempt={attempt} genus={g} requested_threads=24 cache_state=fresh path={}",
            cache.0.path.display()
        );
        let Some(interrupted) = interrupt_before_root(&cache.0.path, g) else {
            cache
                .0
                .remove_files()
                .expect("clean raced recovery attempt");
            continue;
        };
        println!(
            "measurement_recovery_interrupted attempt={attempt} genus={g} committed_rows={} expensive_desc_rows={} root_present=false process_until_reopen_s={:.9}",
            interrupted.rows,
            interrupted.expensive_descendant_rows,
            interrupted.elapsed.as_secs_f64()
        );
        finish_recovery_phase(&cache.0.path, "resume", g);
        finish_recovery_phase(&cache.0.path, "complete_root_hit", g);
        cache
            .0
            .remove_files()
            .expect("clean completed recovery measurement");
        return;
    }
    panic!("all three fresh attempts committed the root before a verified interruption");
}

#[test]
#[ignore = "manual different-genus sequence on one calling thread in each fresh process"]
fn measure_lambda1_sequence() {
    isolated("measure_lambda1_sequence", None, || {
        let first = setting("SC_FABER_MEASURE_START_GENUS", 10);
        let last = setting("SC_FABER_MEASURE_END_GENUS", 12);
        assert!(
            first >= 2 && first < last,
            "sequence requires 2 <= start < end"
        );
        let mut total = Duration::ZERO;
        for g in first..=last {
            let start = Instant::now();
            let value = lambda1_top(g);
            let elapsed = start.elapsed();
            total += elapsed;
            check_known_top(g, &value);
            println!(
                "measurement_query workload=lambda1_sequence genus={g} elapsed_s={:.9} cumulative_query_s={:.9} value={value}",
                elapsed.as_secs_f64(),
                total.as_secs_f64()
            );
        }
        println!(
            "measurement_sequence start_genus={first} end_genus={last} total_query_s={:.9}",
            total.as_secs_f64()
        );
    });
}

struct SmallCase {
    name: &'static str,
    g: u32,
    n: u32,
    kappa: &'static [(u32, u32)],
    psi: &'static [(u32, u32)],
    lambda: &'static [(u32, u32)],
    ch: &'static [(u32, u32)],
}

impl SmallCase {
    fn run(&self) -> Q {
        klp(self.g, self.n, self.kappa, self.psi, self.lambda, self.ch)
            .unwrap_or_else(|error| panic!("{}: {error}", self.name))
    }
}

#[test]
#[ignore = "manual serial formula, descendant and mixed workloads; repeats are warm"]
fn measure_small_repeated_workloads() {
    isolated("measure_small_repeated_workloads", Some(1), || {
        let cases = [
            SmallCase {
                name: "lambda_g_formula",
                g: 6,
                n: 1,
                kappa: &[],
                psi: &[(1, 10)],
                lambda: &[(6, 1)],
                ch: &[],
            },
            SmallCase {
                name: "one_point_psi",
                g: 4,
                n: 1,
                kappa: &[],
                psi: &[(1, 10)],
                lambda: &[],
                ch: &[],
            },
            SmallCase {
                name: "serial_lambda_psi",
                g: 3,
                n: 1,
                kappa: &[],
                psi: &[(1, 2)],
                lambda: &[(1, 5)],
                ch: &[],
            },
            SmallCase {
                name: "serial_pure_psi",
                g: 3,
                n: 3,
                kappa: &[],
                psi: &[(1, 2), (2, 3), (3, 4)],
                lambda: &[],
                ch: &[],
            },
            SmallCase {
                name: "mixed_kappa_psi_lambda_ch",
                g: 3,
                n: 2,
                kappa: &[(1, 1)],
                psi: &[(1, 2)],
                lambda: &[(1, 2)],
                ch: &[(1, 3)],
            },
            SmallCase {
                name: "lambda_g_squared_zero",
                g: 3,
                n: 0,
                kappa: &[],
                psi: &[],
                lambda: &[(3, 2)],
                ch: &[],
            },
        ];
        let repetitions = setting("SC_FABER_MEASURE_SMALL_REPEATS", 1000);
        let warm_route = if matches!(source_mode(), "owned" | "pull" | "persist") {
            "public_root_l1"
        } else {
            "original_public_path"
        };
        for case in cases {
            // Cases share one calling thread. This is the first call of this
            // case, not a claim that every mathematical memo table is empty.
            let start = Instant::now();
            let first = case.run();
            let first_elapsed = start.elapsed();
            let start = Instant::now();
            for _ in 0..repetitions {
                black_box(case.run());
            }
            let warm_elapsed = start.elapsed();
            assert_eq!(case.run(), first, "{} changed after repetition", case.name);
            println!(
                "measurement_small workload={} genus={} markings={} first_call_s={:.9} warm_route={warm_route} warm_repetitions={repetitions} warm_total_s={:.9} warm_mean_ns={:.3} value={first}",
                case.name,
                case.g,
                case.n,
                first_elapsed.as_secs_f64(),
                warm_elapsed.as_secs_f64(),
                warm_elapsed.as_nanos() as f64 / f64::from(repetitions)
            );
        }
    });
}
