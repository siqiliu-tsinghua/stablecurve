//! Manual fresh-process measurements of the unchanged general GRR kernel.
//!
//! `cargo test --release -p sc-faber --lib hodge::measure_pressure::raw_grr -- --ignored --exact --nocapture --test-threads=1`
//! SC_FABER_PRESSURE_GENUS accepts 10..=12 or 16; SC_FABER_PRESSURE_MIB
//! defaults to 600. Each run needs an absent SC_FABER_PRESSURE_DB, or reserves
//! its own name.
//! SC_FABER_PRESSURE_INSPECT=1 enables the existing test-only memory sampler;
//! those diagnostic runs must be kept separate from ordinary timing runs.
//! SC_FABER_PRESSURE_STORAGE=off is a no-persistence control, requiring at
//! least 2048 MiB; its engine has no pressure spill path. Disk is the default.
//! The owned child has a 30-second watchdog (SC_FABER_PRESSURE_TIMEOUT_SECS
//! accepts at most 60 for g<=12, or 300 for the explicit g=16 measurement).
//! A killed child's SQLite files are retained.
//! SC_FABER_PRESSURE_RECOVER=1 interrupts after a committed nonzero C-descendant
//! and restarts against the same database, before the integral root is saved.

use super::{KlpRoute, closed, unreduced_klp};
use crate::cache_source::{IntegralKey, encode_key_into};
use crate::socle::{CacheSession, set_cache_memory_limit};
use rusqlite::{Connection, OpenFlags, params};
use sc_num::{Q, Zero};
use std::path::{Path, PathBuf};
use std::process::{Child, Command};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

const NAME: &str = "hodge::measure_pressure::raw_grr";

fn setting(name: &str, default: u32) -> u32 {
    std::env::var(name).map_or(default, |value| {
        value.parse().ok().filter(|&v| v > 0).expect(name)
    })
}

struct OwnedChild(Child);

impl Drop for OwnedChild {
    fn drop(&mut self) {
        if !matches!(self.0.try_wait(), Ok(Some(_))) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
}

fn interrupt_after_commit(child: &mut OwnedChild, path: &Path, g: u32, timeout: u32) {
    let started = Instant::now();
    let mut observer: Option<Connection> = None;
    let mut root = Vec::new();
    encode_key_into(
        &IntegralKey::Integral {
            g,
            n: 0,
            kappa: vec![],
            psi: vec![],
            lambda: (g - 2..=g).map(|index| (index, 1)).collect(),
            ch: vec![],
        },
        &mut root,
    );
    let mut zero = Vec::new();
    sc_num::encode_rational_into(&Q::ZERO, &mut zero);
    let max_genus = [u8::try_from(g).unwrap()];
    loop {
        assert!(
            started.elapsed() < Duration::from_secs(u64::from(timeout)),
            "no committed C-descendant before the recovery deadline"
        );
        assert!(
            child.0.try_wait().unwrap().is_none(),
            "child completed before interruption; recovery was not tested"
        );
        if observer.is_none() && path.metadata().is_ok_and(|metadata| metadata.len() > 0) {
            let connection = Connection::open_with_flags(
                path,
                OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
            )
            .unwrap();
            connection.busy_timeout(Duration::from_millis(1)).unwrap();
            observer = Some(connection);
        }
        if let Some(connection) = &observer {
            let state = (|| -> rusqlite::Result<_> {
                if !connection.query_row(
                    "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE name='faber_integrals')",
                    [],
                    |row| row.get::<_, bool>(0),
                )? {
                    return Ok((false, false));
                }
                // The bounded g<=16 workload has one-byte C/Euler ranks. The
                // fifth key byte is genus; ignore pure psi, g<=1 and Q zero.
                connection.query_row(
                    "SELECT EXISTS(SELECT 1 FROM faber_integrals WHERE k>=x'010401' AND k<x'0105' AND substr(k,5,1)>=x'02' AND substr(k,5,1)<=?1 AND v!=?2), EXISTS(SELECT 1 FROM faber_integrals WHERE k=?3)",
                    params![max_genus.as_slice(), zero.as_slice(), root.as_slice()],
                    |row| Ok((row.get::<_, bool>(0)?, row.get::<_, bool>(1)?)),
                )
            })();
            match state {
                Ok((_, true)) => {
                    panic!("root committed before interruption; recovery was not tested")
                }
                Ok((true, false)) => break,
                Ok((false, false)) => {}
                Err(error)
                    if matches!(
                        error.sqlite_error_code(),
                        Some(
                            rusqlite::ErrorCode::DatabaseBusy | rusqlite::ErrorCode::DatabaseLocked
                        )
                    ) => {}
                Err(error) => panic!("observe recovery database: {error}"),
            }
        }
        std::thread::sleep(Duration::from_millis(1));
    }
    drop(observer);
    child.0.kill().expect("terminate the owned raw GRR child");
    assert!(
        !child.0.wait().unwrap().success(),
        "child finished before termination"
    );
    let connection = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY).unwrap();
    let (rows, root_present): (i64, bool) = connection.query_row(
        "SELECT (SELECT count(*) FROM faber_integrals WHERE k>=x'010401' AND k<x'0105' AND substr(k,5,1)>=x'02' AND substr(k,5,1)<=?1 AND v!=?2), EXISTS(SELECT 1 FROM faber_integrals WHERE k=?3)",
        params![max_genus.as_slice(), zero.as_slice(), root.as_slice()],
        |row| Ok((row.get(0)?, row.get(1)?)),
    ).unwrap();
    assert!(
        rows > 0 && !root_present,
        "interruption must preserve committed descendants without the root"
    );
    println!(
        "raw_grr_interrupted genus={g} elapsed_s={:.6} committed_nonzero_c_rows={rows} root_present={root_present}",
        started.elapsed().as_secs_f64()
    );
}

#[test]
#[ignore = "manual small raw-GRR pressure measurement; fresh cache and owned-child timeout"]
fn raw_grr() {
    if cfg!(debug_assertions) {
        panic!("use cargo test --release");
    }
    let g = setting("SC_FABER_PRESSURE_GENUS", 12);
    assert!(
        (10..=12).contains(&g) || g == 16,
        "development samples require 10 <= g <= 12, or the explicit g=16 case"
    );
    let mib = setting("SC_FABER_PRESSURE_MIB", 600);
    let storage = std::env::var("SC_FABER_PRESSURE_STORAGE").unwrap_or_else(|_| "disk".into());
    assert!(
        storage == "disk" || storage == "off",
        "storage must be disk or off"
    );
    assert!(
        storage != "off" || mib >= 2048,
        "the off control is not a pressure test"
    );
    let timeout = setting("SC_FABER_PRESSURE_TIMEOUT_SECS", 30);
    assert!(
        timeout <= if g == 16 { 300 } else { 60 },
        "development samples allow at most 60 seconds, or 300 for g=16"
    );
    let path = std::env::var_os("SC_FABER_PRESSURE_DB").map(PathBuf::from);
    let child_mode = std::env::var("SC_FABER_PRESSURE_CHILD").as_deref() == Ok("1");
    let restart = std::env::var("SC_FABER_PRESSURE_RESTART").as_deref() == Ok("1");
    assert!(
        !restart || child_mode,
        "only the owned recovery child may reuse a database"
    );
    if let Some(path) = &path {
        assert_eq!(
            path.exists(),
            restart,
            "require a fresh database, except for the explicit recovery restart"
        );
    }
    if !child_mode {
        let path = path.unwrap_or_else(|| {
            std::env::temp_dir().join(format!(
                "sc-faber-raw-grr-{}-{}.sqlite",
                std::process::id(),
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ))
        });
        let mut command = Command::new(std::env::current_exe().unwrap());
        command
            .args([
                "--ignored",
                "--exact",
                NAME,
                "--nocapture",
                "--test-threads=1",
            ])
            .env("SC_FABER_PRESSURE_CHILD", "1")
            .env("SC_FABER_PRESSURE_DB", &path)
            .env("SC_FABER_THREADS", "24")
            .env("SC_FABER_CACHE_MEASURE_REPORT", "1")
            .env_remove("SC_FABER_CACHE_MEASURE_PULL")
            .env_remove("SC_FABER_CACHE_MEASURE_SOURCE_ONLY");
        if std::env::var("SC_FABER_PRESSURE_INSPECT").as_deref() == Ok("1") {
            command.env("SC_FABER_MEASURE_MEMORY", "1");
        } else {
            command.env_remove("SC_FABER_MEASURE_MEMORY");
        }
        println!(
            "raw_grr_start genus={g} memory_mib={mib} storage={storage} path={}",
            path.display()
        );
        let mut started = Instant::now();
        let mut child = OwnedChild(command.spawn().expect("start raw GRR child"));
        if std::env::var("SC_FABER_PRESSURE_RECOVER").as_deref() == Ok("1") {
            assert_eq!(storage, "disk", "recovery requires SQLite");
            interrupt_after_commit(&mut child, &path, g, timeout);
            command.env("SC_FABER_PRESSURE_RESTART", "1");
            child = OwnedChild(command.spawn().expect("restart raw GRR child"));
            started = Instant::now();
        }
        loop {
            if let Some(status) = child.0.try_wait().expect("poll raw GRR child") {
                println!(
                    "raw_grr_process elapsed_s={:.6} status={status}",
                    started.elapsed().as_secs_f64()
                );
                assert!(
                    status.success(),
                    "raw GRR child failed; database retained at {}",
                    path.display()
                );
                return;
            }
            assert!(
                started.elapsed() < Duration::from_secs(u64::from(timeout)),
                "raw GRR child exceeded {timeout}s; terminating it and retaining {}",
                path.display()
            );
            std::thread::sleep(Duration::from_millis(10));
        }
    }
    set_cache_memory_limit(Some(u64::from(mib) * 1024 * 1024)).unwrap();
    let path = path.expect("parent supplied a SQLite path");
    let lambda = [(g - 2, 1), (g - 1, 1), (g, 1)];
    let started = Instant::now();
    let session = (storage == "disk").then(|| CacheSession::open(path.to_str().unwrap()).unwrap());
    let math = Instant::now();
    let calculate = || {
        let raw = || Ok(unreduced_klp(g, 0, &[], &[], &lambda, &[], KlpRoute::Pure));
        match crate::cache_source::current_registry() {
            Some(registry) => super::root_source::with_registry_policy(
                registry,
                g,
                0,
                &[],
                &[],
                &lambda,
                &[],
                false,
                raw,
            ),
            None => raw(),
        }
    };
    let result = match &session {
        Some(session) => session.run(calculate),
        None => CacheSession::run_uncached(calculate),
    };
    let elapsed = math.elapsed();
    let expected = closed::value(g, 0, &[], &[], &lambda, &[]).unwrap();
    println!(
        "raw_grr_query genus={g} memory_mib={mib} storage={storage} math_s={:.6} result={result:?} bypasses_closed=true",
        elapsed.as_secs_f64()
    );
    if let Some(session) = &session {
        session.flush().expect("flush raw GRR values");
        session.close().expect("close raw GRR cache");
    }
    assert_eq!(result.unwrap(), expected);
    println!(
        "raw_grr_durable elapsed_s={:.6}",
        started.elapsed().as_secs_f64()
    );
}
