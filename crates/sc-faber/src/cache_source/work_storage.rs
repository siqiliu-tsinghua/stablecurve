//! Pressure-only plan workspace IO, owned by the ordinary SQLite service.
//!
//! A logical immutable block is streamed through small rows in a separate
//! disposable SQLite file beside the persistent cache (or in the explicitly
//! selected workspace directory). SQLite's implicit TEMP directory is unused.
//! File backing follows that filesystem; a RAM filesystem still consumes RAM.
//! Neither producer nor backend needs a block-sized encoded allocation.

use super::budget::{Pressure, Reservation, ReservationError, ReservationKind, process_budget};
use super::work::{WorkBlock, WorkSource};
use rusqlite::{CachedStatement, Connection, Rows, params};
use std::ffi::OsStr;
use std::fs::{self, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

const CHUNK_BYTES: usize = 256 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct BatchPolicy {
    bytes: u64,
    age: Duration,
}

fn workspace_batch_policy(value: Option<&str>) -> Result<BatchPolicy, String> {
    match value {
        None | Some("block") => Ok(BatchPolicy {
            bytes: 0,
            age: Duration::ZERO,
        }),
        Some("256k") => Ok(BatchPolicy {
            bytes: 256 * 1024,
            age: Duration::from_millis(10),
        }),
        Some("1m") => Ok(BatchPolicy {
            bytes: 1024 * 1024,
            age: Duration::from_millis(25),
        }),
        _ => Err("SC_FABER_CACHE_MEASURE_WORK_BATCH must be block, 256k or 1m".into()),
    }
}

#[derive(Clone, Copy, Default)]
pub(super) struct WorkReport {
    pub(super) commits: usize,
    pub(super) written_blocks: usize,
    pub(super) written_bytes: u64,
    pub(super) loaded_blocks: usize,
    pub(super) loaded_bytes: u64,
    /// Backend block codec/SQL time plus deferred batch commit checks; excludes
    /// the following reclamation. Large-block probes exclude shared commits.
    pub(super) work_save_ns: u128,
    /// Detailed probes are enabled only by SC_FABER_CACHE_MEASURE_REPORT=1.
    pub(super) save_sql_ns: u128,
    pub(super) save_large_blocks: usize,
    pub(super) save_large_bytes: u64,
    pub(super) save_large_ns: u128,
    pub(super) cache_write_pages: Option<u64>,
    pub(super) cache_spill_pages: Option<u64>,
    pub(super) work_load_ns: u128,
    /// Query setup/reset plus per-chunk stepping, column access and blob copy.
    /// This is a subset of work_load_ns, separate from block decoding.
    pub(super) work_read_sql_ns: u128,
    /// Includes budget refresh, owner release and any post-eviction force sample.
    pub(super) work_reclaim_ns: u128,
    pub(super) max_loaded_block_resident_bytes: usize,
    pub(super) max_block_encoded_bytes: u64,
}

pub(super) struct WorkManager {
    cache_path: PathBuf,
    sources: Vec<WorkEntry>,
    next: usize,
    database: Option<WorkDatabase>,
    initialization_error: Option<String>,
    scratch: Option<Reservation>,
    report: WorkReport,
    diagnostic_timing: Option<bool>,
    cache_status_failed: bool,
    batch_policy: Option<BatchPolicy>,
    transaction_open: bool,
    transaction_bytes: u64,
    transaction_first: Option<Instant>,
    #[cfg(test)]
    transaction_age_override: Option<Duration>,
}

struct WorkEntry {
    source: Arc<dyn WorkSource>,
    reclaim_remaining: usize,
}

impl WorkEntry {
    fn admitted_load(&mut self) {
        // A production workspace has one coordinator. A validated readback
        // token means its whole next group has been admitted and is waiting
        // for these blocks, so its earlier finite admission-drain request is
        // satisfied. Preserve the explicit diagnostic all-future-block drain.
        let posted = self.source.take_reclaim_request();
        self.reclaim_remaining = if self.reclaim_remaining == usize::MAX || posted == usize::MAX {
            usize::MAX
        } else {
            0
        };
    }

    fn reclaim(&mut self, pressure_goal: usize) -> usize {
        self.reclaim_remaining = self
            .reclaim_remaining
            .max(self.source.take_reclaim_request());
        let goal = self.reclaim_remaining.max(pressure_goal);
        if goal == 0 {
            return 0;
        }
        let reclaimed = self.source.reclaim_stored(goal);
        if self.reclaim_remaining != usize::MAX {
            self.reclaim_remaining = self.reclaim_remaining.saturating_sub(reclaimed);
        }
        if reclaimed == 0 && !self.source.reclaim_pending() {
            self.reclaim_remaining = 0;
        }
        reclaimed
    }

    fn needs_save(&self, pressure_goal: usize) -> bool {
        pressure_goal != 0 || self.reclaim_remaining != 0
    }
}

fn workspace_pressure_goal(force: bool) -> usize {
    let budget = process_budget();
    let decision = if force {
        budget.force_check(0)
    } else {
        budget.refresh_if_due()
    };
    if matches!(decision.pressure, Pressure::Warning | Pressure::Hard)
        && !decision.stale
        && decision.observations_complete
    {
        usize::try_from(decision.reclaim_bytes).unwrap_or(usize::MAX)
    } else {
        0
    }
}

impl WorkManager {
    /// No filesystem access or environment reads until the first backend step.
    pub(super) fn new(path: &str) -> Self {
        Self {
            cache_path: PathBuf::from(path),
            sources: Vec::new(),
            next: 0,
            database: None,
            initialization_error: None,
            scratch: None,
            report: WorkReport::default(),
            diagnostic_timing: None,
            cache_status_failed: false,
            batch_policy: None,
            transaction_open: false,
            transaction_bytes: 0,
            transaction_first: None,
            #[cfg(test)]
            transaction_age_override: None,
        }
    }

    pub(super) fn register(&mut self, source: Arc<dyn WorkSource>) {
        source.bind_worker(std::thread::current());
        self.sources.push(WorkEntry {
            source,
            reclaim_remaining: 0,
        });
    }

    pub(super) fn report(&self) -> WorkReport {
        self.report
    }

    fn initialize(&mut self) -> Result<(), String> {
        if let Some(error) = &self.initialization_error {
            return Err(error.clone());
        }
        if self.database.is_some() {
            return Ok(());
        }
        // Read this optional override once, on the backend, including failures.
        let override_directory = std::env::var_os("SC_FABER_WORKSPACE_DIR");
        let initialized = workspace_directory(&self.cache_path, override_directory.as_deref())
            .and_then(|directory| WorkDatabase::open(&directory));
        match initialized {
            Ok(database) => {
                self.database = Some(database);
                Ok(())
            }
            Err(error) => {
                self.initialization_error = Some(error.clone());
                Err(error)
            }
        }
    }

    /// Service at most one immutable block. All source locks are released by
    /// `pull`/`next_load` before encoding, allocation, or SQLite access begins.
    pub(super) fn step(&mut self, _connection: &mut Connection) -> Result<bool, String> {
        if let Some(error) = &self.initialization_error {
            return Err(error.clone());
        }
        if self.database.is_none() {
            self.sources.retain(|entry| !entry.source.is_retired());
        }
        if self.sources.is_empty() {
            return Ok(false);
        }
        if !self.reserve_scratch()? {
            return Ok(false);
        }
        self.initialize()?;
        // Keep the standalone connection separate while dispatch mutates the
        // source rotation. On unwind its local owner still closes and cleans up.
        let mut database = self.database.take().unwrap();
        let result = self.step_on(database.connection.as_mut().unwrap());
        if self.sources.is_empty() {
            // Retiring the last plan drops its disposable file directly rather
            // than traversing all chunk rows just to delete the file afterwards.
            self.record_cache_status(database.connection.as_ref().unwrap());
            drop(database);
            // No live plan can request these disposable rows again. Do not
            // commit an unfinished batch merely to delete its file afterwards.
            self.transaction_open = false;
            self.transaction_bytes = 0;
            self.set_transaction_first(None);
            self.scratch.take();
        } else {
            self.database = Some(database);
        }
        result
    }

    fn reserve_scratch(&mut self) -> Result<bool, String> {
        if self.scratch.is_some() {
            return Ok(true);
        }
        let budget = process_budget();
        let bytes = 4 * 1024 * 1024;
        let mut credit = budget.try_reserve(bytes, ReservationKind::Progress);
        if matches!(
            &credit,
            Err(ReservationError::Stale(_) | ReservationError::Pressure { .. })
        ) {
            budget.force_check(0);
            credit = budget.try_reserve(bytes, ReservationKind::Progress);
        }
        match credit {
            Ok(credit) => {
                credit.start().map_err(|error| error.to_string())?;
                // One reusable envelope covers this backend's bounded stream
                // buffers and 2 MiB working page cache for its whole lifetime.
                self.scratch = Some(credit);
                Ok(true)
            }
            // Do not block storage waiting for its own progress. The controller
            // has bounded waits and can shrink or reject the pending group.
            Err(
                ReservationError::Busy
                | ReservationError::Stale(_)
                | ReservationError::Pressure { .. },
            ) => Ok(false),
            Err(error) => Err(error.to_string()),
        }
    }

    fn step_on(&mut self, connection: &mut Connection) -> Result<bool, String> {
        self.step_on_with_pressure(connection, &mut workspace_pressure_goal)
    }

    fn step_on_with_pressure(
        &mut self,
        connection: &mut Connection,
        pressure: &mut impl FnMut(bool) -> usize,
    ) -> Result<bool, String> {
        if let Some(error) = &self.initialization_error {
            return Err(error.clone());
        }
        let result = self
            .step_on_active(connection, pressure)
            .and_then(|progress| {
                // Read requests retain priority, and idle turns also service the
                // deadline. The last retired owner instead discards the whole DB.
                let started = self.transaction_open.then(Instant::now);
                let committed = !self.sources.is_empty() && self.commit_due(connection)?;
                if committed {
                    self.record_cache_status(connection);
                }
                if let Some(started) = started {
                    self.report.work_save_ns = self
                        .report
                        .work_save_ns
                        .saturating_add(started.elapsed().as_nanos());
                }
                Ok(progress || committed)
            });
        if let Err(error) = &result {
            // OFF mode cannot safely roll back an unsuccessful batch. Earlier
            // blocks may already have been acknowledged and evicted, so every
            // owner must stop; never retry SQL on this connection after failure.
            self.fail(error);
            self.initialization_error = Some(error.clone());
            #[cfg(test)]
            {
                self.transaction_age_override = None;
            }
        }
        result
    }

    fn step_on_active(
        &mut self,
        connection: &mut Connection,
        pressure: &mut impl FnMut(bool) -> usize,
    ) -> Result<bool, String> {
        let index = self.next % self.sources.len();
        let source = Arc::clone(&self.sources[index].source);
        self.next = (index + 1) % self.sources.len();
        let owner = sql_id(source.workspace_id())?;
        if source.is_retired() {
            if self.sources.len() > 1 {
                connection
                    .execute("DELETE FROM faber_work_blocks WHERE owner=?1", [owner])
                    .map_err(sql_error)?;
            }
            self.sources.remove(index);
            self.next = index;
            return Ok(true);
        }
        // A read unlocks a waiting calculation; dirty scalar persistence still
        // gets its normal bounded turn in the outer storage-service loop.
        if let Some(load) = source.next_load() {
            let started = Instant::now();
            self.report.max_loaded_block_resident_bytes = self
                .report
                .max_loaded_block_resident_bytes
                .max(load.resident_bytes);
            let credit = load.reservation.as_ref();
            if let Some(credit) = credit {
                if let Err(error) = credit.start() {
                    if source.is_retired() {
                        return Ok(true);
                    }
                    return Err(error.to_string());
                }
                if credit.reserved_bytes() < load.resident_bytes as u64 {
                    return Err("plan readback exceeds its coordinator reservation".into());
                }
                self.sources[index].admitted_load();
            } else if load.resident_bytes > 1024 * 1024 {
                // Legacy diagnostic sources can use the small scratch window.
                // Production loads must carry the coordinator's whole-group
                // credit instead of silently allocating another large block.
                return Err("plan readback requires a coordinator reservation".into());
            }
            let block = sql_block(load.block)?;
            let mut read_sql_ns = 0u128;
            let restored = (|| -> Result<u64, String> {
                let sql_started = Instant::now();
                let mut statement = connection
                    .prepare_cached(
                        "SELECT chunk,data FROM faber_work_blocks
                     WHERE owner=?1 AND block=?2 ORDER BY chunk",
                    )
                    .map_err(sql_error)?;
                let mut rows = statement.query(params![owner, block]).map_err(sql_error)?;
                read_sql_ns = sql_started.elapsed().as_nanos();
                let mut reader = ChunkReader::new(&mut rows);
                let result = source.restore(load.block, &mut reader);
                let consumed = reader.consumed;
                read_sql_ns = read_sql_ns.saturating_add(reader.sql_ns);
                let release_started = Instant::now();
                drop(reader);
                drop(rows);
                drop(statement);
                read_sql_ns = read_sql_ns.saturating_add(release_started.elapsed().as_nanos());
                result.map(|()| consumed)
            })();
            self.report.work_load_ns = self
                .report
                .work_load_ns
                .saturating_add(started.elapsed().as_nanos());
            self.report.work_read_sql_ns = self.report.work_read_sql_ns.saturating_add(read_sql_ns);
            let consumed = restored?;
            self.report.loaded_blocks += 1;
            self.report.loaded_bytes = self.report.loaded_bytes.saturating_add(consumed);
            self.report.max_block_encoded_bytes = self.report.max_block_encoded_bytes.max(consumed);
            self.reclaim(index, pressure);
            return Ok(true);
        }
        // Backing is an on-demand pressure tool, not a request to write the
        // complete remaining Plan. Prefer already saved, unpinned owners, and
        // recheck the actual native state after their memory has been dropped.
        let (reclaimed, needs_save) = self.reclaim(index, pressure);
        if !needs_save {
            return Ok(reclaimed != 0);
        }
        if let Some(block) = source.pull() {
            let started = Instant::now();
            let saved = self.save(connection, owner, &block);
            let elapsed = started.elapsed().as_nanos();
            self.report.work_save_ns = self.report.work_save_ns.saturating_add(elapsed);
            let written = saved?;
            if self.diagnostic_timing() && written > 2 * 1024 * 1024 {
                self.report.save_large_blocks += 1;
                self.report.save_large_bytes = self.report.save_large_bytes.saturating_add(written);
                self.report.save_large_ns = self.report.save_large_ns.saturating_add(elapsed);
            }
            self.report.written_blocks += 1;
            self.report.written_bytes = self.report.written_bytes.saturating_add(written);
            self.report.max_block_encoded_bytes = self.report.max_block_encoded_bytes.max(written);
            source.stored(block.block)?;
            // Release this backend pin before optional pressure eviction.
            // Without pressure the saved block remains owned by its source.
            drop(block);
            self.reclaim(index, pressure);
            return Ok(true);
        }
        Ok(reclaimed != 0)
    }

    fn reclaim(&mut self, index: usize, pressure: &mut impl FnMut(bool) -> usize) -> (usize, bool) {
        let started = Instant::now();
        let mut pressure_goal = pressure(false);
        let reclaimed = self.sources[index].reclaim(pressure_goal);
        if reclaimed != 0 {
            // Resample after a bounded eviction turn; an old high-pressure
            // observation must not drive every following turn to evict again.
            pressure_goal = pressure(true);
        }
        let needs_save = self.sources[index].needs_save(pressure_goal);
        self.report.work_reclaim_ns = self
            .report
            .work_reclaim_ns
            .saturating_add(started.elapsed().as_nanos());
        (reclaimed, needs_save)
    }

    fn save(
        &mut self,
        connection: &mut Connection,
        owner: i64,
        block: &WorkBlock,
    ) -> Result<u64, String> {
        let block_id = sql_block(block.block)?;
        let expected = block.payload.encoded_len()?;
        let detailed = self.diagnostic_timing();
        let policy = self.batch_policy()?;
        if !self.transaction_open {
            let started = detailed.then(Instant::now);
            let begun = connection
                .execute_batch("BEGIN DEFERRED")
                .map_err(sql_error);
            self.add_save_sql_time(started);
            begun?;
            self.transaction_open = true;
            // Batch timing is opt-in; the original per-block path does not
            // acquire an additional clock reading for this policy.
            self.set_transaction_first((policy.bytes != 0).then(Instant::now));
        }
        let started = detailed.then(Instant::now);
        let deleted = connection
            .execute(
                "DELETE FROM faber_work_blocks WHERE owner=?1 AND block=?2",
                params![owner, block_id],
            )
            .map_err(sql_error);
        self.add_save_sql_time(started);
        deleted?;
        {
            let statement = connection
                .prepare_cached(
                    "INSERT INTO faber_work_blocks(owner,block,chunk,data)
                     VALUES(?1,?2,?3,?4)",
                )
                .map_err(sql_error)?;
            let mut writer = ChunkWriter {
                statement,
                owner,
                block: block_id,
                chunk: 0,
                bytes: 0,
                buffer: Vec::with_capacity(CHUNK_BYTES),
                detailed_timing: detailed,
                sql_ns: 0,
            };
            let encoded = block
                .payload
                .encode(&mut writer)
                .and_then(|()| writer.flush().map_err(|error| error.to_string()));
            self.report.save_sql_ns = self.report.save_sql_ns.saturating_add(writer.sql_ns);
            encoded?;
            if writer.bytes != expected {
                return Err("plan workspace encoder returned an inconsistent length".into());
            }
        }
        self.transaction_bytes = self
            .transaction_bytes
            .checked_add(expected)
            .ok_or("plan workspace transaction byte count overflow")?;
        if policy.bytes == 0 {
            self.commit_work(connection)?;
        }
        // For opt-in batches, complete rows are already readable through this
        // same connection before COMMIT. The caller may acknowledge and release
        // this block immediately, without a pending payload or ACK collection.
        self.record_cache_status(connection);
        Ok(expected)
    }

    fn batch_policy(&mut self) -> Result<BatchPolicy, String> {
        if let Some(policy) = self.batch_policy {
            return Ok(policy);
        }
        let requested = match std::env::var("SC_FABER_CACHE_MEASURE_WORK_BATCH") {
            Ok(value) => Some(value),
            Err(std::env::VarError::NotPresent) => None,
            Err(error) => return Err(format!("SC_FABER_CACHE_MEASURE_WORK_BATCH: {error}")),
        };
        let policy = workspace_batch_policy(requested.as_deref())?;
        self.batch_policy = Some(policy);
        Ok(policy)
    }

    fn commit_due(&mut self, connection: &Connection) -> Result<bool, String> {
        if !self.transaction_open {
            return Ok(false);
        }
        let policy = self
            .batch_policy
            .expect("open workspace transaction has a policy");
        if self.transaction_bytes >= policy.bytes
            || self.transaction_age().is_some_and(|age| age >= policy.age)
        {
            self.commit_work(connection)?;
            return Ok(true);
        }
        Ok(false)
    }

    fn commit_work(&mut self, connection: &Connection) -> Result<(), String> {
        let started = self.diagnostic_timing().then(Instant::now);
        let committed = connection.execute_batch("COMMIT").map_err(sql_error);
        self.add_save_sql_time(started);
        committed?;
        self.transaction_open = false;
        self.transaction_bytes = 0;
        self.set_transaction_first(None);
        self.report.commits += 1;
        Ok(())
    }

    fn set_transaction_first(&mut self, first: Option<Instant>) {
        self.transaction_first = first;
        #[cfg(test)]
        {
            self.transaction_age_override = None;
        }
    }

    fn transaction_age(&self) -> Option<Duration> {
        let first = self.transaction_first?;
        #[cfg(test)]
        if let Some(age) = self.transaction_age_override {
            return Some(age);
        }
        Some(first.elapsed())
    }

    fn diagnostic_timing(&mut self) -> bool {
        *self.diagnostic_timing.get_or_insert_with(|| {
            std::env::var("SC_FABER_CACHE_MEASURE_REPORT").as_deref() == Ok("1")
        })
    }

    fn add_save_sql_time(&mut self, started: Option<Instant>) {
        if let Some(started) = started {
            self.report.save_sql_ns = self
                .report
                .save_sql_ns
                .saturating_add(started.elapsed().as_nanos());
        }
    }

    fn record_cache_status(&mut self, connection: &Connection) {
        if !self.diagnostic_timing() || self.cache_status_failed {
            return;
        }
        match workspace_cache_writes(connection) {
            Some((written, spilled)) => {
                self.report.cache_write_pages = Some(
                    self.report
                        .cache_write_pages
                        .unwrap_or(0)
                        .saturating_add(written),
                );
                self.report.cache_spill_pages = Some(
                    self.report
                        .cache_spill_pages
                        .unwrap_or(0)
                        .saturating_add(spilled),
                );
            }
            None => {
                self.cache_status_failed = true;
                self.report.cache_write_pages = None;
                self.report.cache_spill_pages = None;
            }
        }
    }

    pub(super) fn fail(&self, error: &str) {
        for entry in &self.sources {
            entry.source.fail(error.to_owned());
        }
    }
}

fn workspace_directory(
    cache: &Path,
    override_directory: Option<&OsStr>,
) -> Result<PathBuf, String> {
    if let Some(directory) = override_directory {
        if directory.is_empty() {
            return Err("SC_FABER_WORKSPACE_DIR must name a filesystem directory".into());
        }
        return Ok(PathBuf::from(directory));
    }
    if cache == Path::new(":memory:") || cache.as_os_str().is_empty() {
        return Err(
            "memory-backed cache has no disk workspace directory; set SC_FABER_WORKSPACE_DIR explicitly before pressure backing is needed".into(),
        );
    }
    Ok(cache
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf())
}

fn reserve_workspace_file(directory: &Path) -> Result<PathBuf, String> {
    let process = std::process::id();
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| format!("plan workspace clock: {error}"))?
        .as_nanos();
    for attempt in 0..100 {
        let path = directory.join(format!(
            "sc-faber-work-{process}-{timestamp}-{attempt}.sqlite"
        ));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => {
                drop(file);
                return Ok(path);
            }
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(format!("create plan workspace {}: {error}", path.display())),
        }
    }
    Err("cannot reserve a unique plan workspace file".into())
}

struct WorkDatabase {
    /// This path is owned only because create_new succeeded for this file.
    path: PathBuf,
    connection: Option<Connection>,
}

fn workspace_locking_mode(value: Option<&str>) -> Result<&'static str, String> {
    match value {
        None => Ok("exclusive"),
        Some(value) if value.eq_ignore_ascii_case("normal") => Ok("normal"),
        Some(value) if value.eq_ignore_ascii_case("exclusive") => Ok("exclusive"),
        _ => Err("SC_FABER_CACHE_MEASURE_WORK_LOCKING_MODE must be normal or exclusive".into()),
    }
}

fn workspace_page_bytes(value: Option<&str>) -> Result<u32, String> {
    match value {
        None | Some("65536") => Ok(65536),
        Some("4096") => Ok(4096),
        Some("16384") => Ok(16384),
        _ => Err("SC_FABER_CACHE_MEASURE_WORK_PAGE_BYTES must be 4096, 16384 or 65536".into()),
    }
}

impl WorkDatabase {
    fn open(directory: &Path) -> Result<Self, String> {
        let requested = match std::env::var("SC_FABER_CACHE_MEASURE_WORK_LOCKING_MODE") {
            Ok(value) => Some(value),
            Err(std::env::VarError::NotPresent) => None,
            Err(error) => return Err(format!("SC_FABER_CACHE_MEASURE_WORK_LOCKING_MODE: {error}")),
        };
        let mode = workspace_locking_mode(requested.as_deref())?;
        let requested_pages = match std::env::var("SC_FABER_CACHE_MEASURE_WORK_PAGE_BYTES") {
            Ok(value) => Some(value),
            Err(std::env::VarError::NotPresent) => None,
            Err(error) => return Err(format!("SC_FABER_CACHE_MEASURE_WORK_PAGE_BYTES: {error}")),
        };
        let page_bytes = workspace_page_bytes(requested_pages.as_deref())?;
        Self::open_with_settings(directory, mode, page_bytes)
    }

    fn open_with_settings(directory: &Path, mode: &str, page_bytes: u32) -> Result<Self, String> {
        let mut database = Self {
            path: reserve_workspace_file(directory)?,
            connection: None,
        };
        let connection = Connection::open(&database.path).map_err(sql_error)?;
        database.connection = Some(connection);
        // SQLite adopts page_size when an empty database is first created.
        // Set it before any schema or journal access; cache_size below remains
        // a fixed KiB target after the selected page size has taken effect.
        // https://www.sqlite.org/pragma.html#pragma_page_size
        database
            .connection
            .as_ref()
            .unwrap()
            .pragma_update(None, "page_size", i64::from(page_bytes))
            .map_err(sql_error)?;
        // Only this newly created, disposable file has one private connection.
        // The EXCLUSIVE candidate retains SQLite's native file locks between
        // transactions, avoiding repeated lock/unlock OS calls. The persistent
        // integral database is a different connection and retains WAL/FULL.
        // https://www.sqlite.org/pragma.html#pragma_locking_mode
        database
            .connection
            .as_ref()
            .unwrap()
            .pragma_update(None, "locking_mode", mode)
            .map_err(sql_error)?;
        database
            .connection
            .as_ref()
            .unwrap()
            .execute_batch(
                "PRAGMA journal_mode=OFF;
                 PRAGMA synchronous=OFF;
                 PRAGMA cache_size=-2048;
                 CREATE TABLE faber_work_blocks (
                     owner INTEGER NOT NULL, block INTEGER NOT NULL,
                     chunk INTEGER NOT NULL, data BLOB NOT NULL,
                     PRIMARY KEY(owner,block,chunk)
                 ) WITHOUT ROWID;",
            )
            .map_err(sql_error)?;
        Ok(database)
    }
}

impl Drop for WorkDatabase {
    fn drop(&mut self) {
        if let Some(connection) = self.connection.take() {
            if let Err((connection, error)) = connection.close() {
                eprintln!(
                    "sc-faber: plan workspace close failed ({error}); retaining {}",
                    self.path.display()
                );
                drop(connection);
                return;
            }
        }
        match fs::remove_file(&self.path) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => eprintln!(
                "sc-faber: plan workspace cleanup failed ({error}); retaining {}",
                self.path.display()
            ),
        }
    }
}

fn sql_id(value: u64) -> Result<i64, String> {
    i64::try_from(value).map_err(|_| "plan workspace ID exceeds SQLite integer range".into())
}

fn sql_block(value: usize) -> Result<i64, String> {
    i64::try_from(value).map_err(|_| "plan block ID exceeds SQLite integer range".into())
}

fn sql_error(error: rusqlite::Error) -> String {
    format!("plan workspace SQLite: {error}")
}

fn workspace_cache_writes(connection: &Connection) -> Option<(u64, u64)> {
    fn read(connection: &Connection, operation: i32) -> Option<u64> {
        let mut current = 0;
        let mut highwater = 0;
        // SAFETY: this private connection is live and used only by its owner
        // thread, after SQL statements/transactions and mathematical source
        // locks have been released. The raw handle never escapes this call.
        // sqlite3_db_status writes two initialized integers; reset=1 changes
        // diagnostic counters only, not database pages or settings.
        let result = unsafe {
            rusqlite::ffi::sqlite3_db_status(
                connection.handle(),
                operation,
                &mut current,
                &mut highwater,
                1,
            )
        };
        if result != rusqlite::ffi::SQLITE_OK {
            return None;
        }
        u64::try_from(current).ok()
    }
    Some((
        read(connection, rusqlite::ffi::SQLITE_DBSTATUS_CACHE_WRITE)?,
        read(connection, rusqlite::ffi::SQLITE_DBSTATUS_CACHE_SPILL)?,
    ))
}

fn io_error(error: impl std::fmt::Display) -> io::Error {
    io::Error::other(error.to_string())
}

struct ChunkWriter<'connection> {
    statement: CachedStatement<'connection>,
    owner: i64,
    block: i64,
    chunk: i64,
    bytes: u64,
    buffer: Vec<u8>,
    detailed_timing: bool,
    sql_ns: u128,
}

impl ChunkWriter<'_> {
    fn write_chunk(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let started = self.detailed_timing.then(Instant::now);
        let inserted = self
            .statement
            .execute(params![self.owner, self.block, self.chunk, &self.buffer])
            .map_err(io_error);
        if let Some(started) = started {
            self.sql_ns = self.sql_ns.saturating_add(started.elapsed().as_nanos());
        }
        inserted?;
        self.chunk = self
            .chunk
            .checked_add(1)
            .ok_or_else(|| io_error("too many plan chunks"))?;
        self.buffer.clear();
        Ok(())
    }
}

impl Write for ChunkWriter<'_> {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        let mut remaining = input;
        while !remaining.is_empty() {
            let count = remaining.len().min(CHUNK_BYTES - self.buffer.len());
            self.buffer.extend_from_slice(&remaining[..count]);
            remaining = &remaining[count..];
            self.bytes = self
                .bytes
                .checked_add(count as u64)
                .ok_or_else(|| io_error("encoded plan length exceeds u64"))?;
            if self.buffer.len() == CHUNK_BYTES {
                self.write_chunk()?;
            }
        }
        Ok(input.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.write_chunk()
    }
}

struct ChunkReader<'rows, 'statement> {
    rows: &'rows mut Rows<'statement>,
    chunk: i64,
    buffer: Vec<u8>,
    position: usize,
    finished: bool,
    consumed: u64,
    sql_ns: u128,
}

impl<'rows, 'statement> ChunkReader<'rows, 'statement> {
    fn new(rows: &'rows mut Rows<'statement>) -> Self {
        Self {
            rows,
            chunk: 0,
            buffer: Vec::new(),
            position: 0,
            finished: false,
            consumed: 0,
            sql_ns: 0,
        }
    }

    fn next_chunk(&mut self) -> io::Result<()> {
        let started = Instant::now();
        let result = (|| {
            self.buffer.clear();
            self.position = 0;
            let Some(row) = self.rows.next().map_err(io_error)? else {
                self.finished = true;
                return Ok(());
            };
            let index: i64 = row.get(0).map_err(io_error)?;
            if index != self.chunk {
                return Err(io_error("missing or unordered plan workspace chunk"));
            }
            let value = row.get_ref(1).map_err(io_error)?;
            let chunk = value.as_blob().map_err(io_error)?;
            if chunk.is_empty() || chunk.len() > CHUNK_BYTES {
                return Err(io_error("invalid plan workspace chunk length"));
            }
            // Reserve exactly: automatic doubling could retain more than one
            // chunk when an earlier, shorter SQL row was read first.
            self.buffer
                .try_reserve_exact(chunk.len())
                .map_err(io_error)?;
            if self.buffer.capacity() > CHUNK_BYTES {
                return Err(io_error("plan workspace chunk buffer exceeds its bound"));
            }
            self.buffer.extend_from_slice(chunk);
            self.chunk = self
                .chunk
                .checked_add(1)
                .ok_or_else(|| io_error("too many plan chunks"))?;
            Ok(())
        })();
        self.sql_ns = self.sql_ns.saturating_add(started.elapsed().as_nanos());
        result
    }
}

impl Read for ChunkReader<'_, '_> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if output.is_empty() || self.finished {
            return Ok(0);
        }
        if self.position == self.buffer.len() {
            self.next_chunk()?;
            if self.finished {
                return Ok(0);
            }
        }
        let count = output.len().min(self.buffer.len() - self.position);
        output[..count].copy_from_slice(&self.buffer[self.position..self.position + count]);
        self.position += count;
        self.consumed = self.consumed.saturating_add(count as u64);
        Ok(count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache_source::work::{WorkLoad, WorkPayload};
    use std::sync::Mutex;

    fn chunk_database(chunks: &[(i64, Vec<u8>)]) -> Connection {
        let connection = Connection::open_in_memory().unwrap();
        connection
            .execute_batch(
                "CREATE TABLE faber_work_blocks (
                owner INTEGER NOT NULL, block INTEGER NOT NULL,
                chunk INTEGER NOT NULL, data BLOB NOT NULL,
                PRIMARY KEY(owner,block,chunk)
             ) WITHOUT ROWID;",
            )
            .unwrap();
        for (chunk, bytes) in chunks {
            connection
                .execute(
                    "INSERT INTO faber_work_blocks VALUES(1,2,?1,?2)",
                    params![chunk, bytes],
                )
                .unwrap();
        }
        connection
    }

    fn read_chunks(connection: &Connection) -> io::Result<Vec<u8>> {
        let mut statement = connection
            .prepare_cached(
                "SELECT chunk,data FROM faber_work_blocks WHERE owner=1 AND block=2 ORDER BY chunk",
            )
            .unwrap();
        let mut rows = statement.query([]).unwrap();
        let mut reader = ChunkReader::new(&mut rows);
        let mut decoded = Vec::new();
        reader.read_to_end(&mut decoded)?;
        assert_eq!(reader.consumed, decoded.len() as u64);
        assert!(reader.buffer.capacity() <= CHUNK_BYTES);
        let sql_ns = reader.sql_ns;
        assert_eq!(reader.read(&mut [0; 1])?, 0);
        assert_eq!(
            reader.sql_ns, sql_ns,
            "finished cursor must not query again"
        );
        Ok(decoded)
    }

    #[test]
    fn ordered_chunk_cursor_crosses_rows_and_rejects_missing_or_invalid_chunks() {
        // SQL insertion order is irrelevant: the whole-block cursor orders
        // chunk identities and retains only one bounded owned chunk buffer.
        let connection = chunk_database(&[
            (2, vec![9; 17]),
            (0, vec![3; 37]),
            (1, vec![7; CHUNK_BYTES]),
        ]);
        let decoded = read_chunks(&connection).unwrap();
        let expected = [vec![3; 37], vec![7; CHUNK_BYTES], vec![9; 17]].concat();
        assert_eq!(decoded, expected);
        for chunks in [
            vec![(1, vec![1])],
            vec![(0, vec![1]), (2, vec![2])],
            vec![(0, Vec::new())],
            vec![(0, vec![0; CHUNK_BYTES + 1])],
        ] {
            assert!(read_chunks(&chunk_database(&chunks)).is_err());
        }
        let empty = chunk_database(&[]);
        let mut statement = empty
            .prepare("SELECT chunk,data FROM faber_work_blocks ORDER BY chunk")
            .unwrap();
        let mut rows = statement.query([]).unwrap();
        let mut reader = ChunkReader::new(&mut rows);
        assert_eq!(
            reader.read_exact(&mut [0; 1]).unwrap_err().kind(),
            io::ErrorKind::UnexpectedEof
        );
    }

    #[test]
    fn cursor_exposes_trailing_rows_and_resets_when_readback_is_abandoned() {
        let mut connection = chunk_database(&[(0, vec![1, 2, 3]), (1, vec![4])]);
        {
            let mut statement = connection.prepare_cached(
                "SELECT chunk,data FROM faber_work_blocks WHERE owner=1 AND block=2 ORDER BY chunk",
            ).unwrap();
            let mut rows = statement.query([]).unwrap();
            let mut reader = ChunkReader::new(&mut rows);
            let mut payload = [0; 3];
            reader.read_exact(&mut payload).unwrap();
            assert_eq!(payload, [1, 2, 3]);
            let mut trailing = [0];
            assert_eq!(reader.read(&mut trailing).unwrap(), 1);
            assert_eq!(trailing, [4]);
            // Leave the SELECT unfinished, as a retired workspace can do.
        }
        let transaction = connection.transaction().unwrap();
        transaction
            .execute("DELETE FROM faber_work_blocks", [])
            .unwrap();
        transaction.commit().unwrap();
        assert!(read_chunks(&connection).unwrap().is_empty());
    }

    #[test]
    fn private_workspace_locking_candidates_preserve_separate_database_settings() {
        assert_eq!(workspace_locking_mode(None).unwrap(), "exclusive");
        assert_eq!(workspace_locking_mode(Some("NORMAL")).unwrap(), "normal");
        assert_eq!(
            workspace_locking_mode(Some("exclusive")).unwrap(),
            "exclusive"
        );
        assert!(workspace_locking_mode(Some("off")).is_err());
        let main_path = reserve_workspace_file(&std::env::temp_dir()).unwrap();
        let main = Connection::open(&main_path).unwrap();
        main.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL;")
            .unwrap();
        let main_locking = pragma_text(&main, "locking_mode");
        let main_journal = pragma_text(&main, "journal_mode");
        for mode in ["normal", "exclusive"] {
            let database =
                WorkDatabase::open_with_settings(&std::env::temp_dir(), mode, 4096).unwrap();
            let path = database.path.clone();
            let connection = database.connection.as_ref().unwrap();
            assert_eq!(pragma_text(connection, "locking_mode"), mode);
            assert_eq!(pragma_text(connection, "journal_mode"), "off");
            assert_eq!(pragma_integer(connection, "synchronous"), 0);
            connection
                .execute("INSERT INTO faber_work_blocks VALUES(1,2,0,x'1234')", [])
                .unwrap();
            for _ in 0..2 {
                let data: Vec<u8> = connection
                    .query_row(
                        "SELECT data FROM faber_work_blocks WHERE owner=1 AND block=2 AND chunk=0",
                        [],
                        |row| row.get(0),
                    )
                    .unwrap();
                assert_eq!(data, [0x12, 0x34]);
            }
            assert_eq!(pragma_text(&main, "locking_mode"), main_locking);
            assert_eq!(pragma_text(&main, "journal_mode"), main_journal);
            assert_eq!(pragma_integer(&main, "synchronous"), 2);
            drop(database);
            assert!(!path.exists());
        }
        assert_eq!(main_journal, "wal");
        assert_eq!(main_locking, "normal");
        main.close().map_err(|(_, error)| error).unwrap();
        remove_test_file(&main_path);
    }

    #[test]
    fn private_workspace_page_candidates_are_applied_before_schema_creation() {
        assert_eq!(workspace_page_bytes(None).unwrap(), 65536);
        for (requested, expected) in [("4096", 4096), ("16384", 16384), ("65536", 65536)] {
            assert_eq!(workspace_page_bytes(Some(requested)).unwrap(), expected);
            let database =
                WorkDatabase::open_with_settings(&std::env::temp_dir(), "exclusive", expected)
                    .unwrap();
            let path = database.path.clone();
            let connection = database.connection.as_ref().unwrap();
            assert_eq!(pragma_integer(connection, "page_size"), i64::from(expected));
            assert_eq!(pragma_integer(connection, "cache_size"), -2048);
            assert_eq!(pragma_text(connection, "journal_mode"), "off");
            assert_eq!(pragma_integer(connection, "synchronous"), 0);
            connection
                .execute("INSERT INTO faber_work_blocks VALUES(1,2,0,x'1234')", [])
                .unwrap();
            let data: Vec<u8> = connection
                .query_row(
                    "SELECT data FROM faber_work_blocks WHERE owner=1 AND block=2 AND chunk=0",
                    [],
                    |row| row.get(0),
                )
                .unwrap();
            assert_eq!(data, [0x12, 0x34]);
            drop(database);
            assert!(!path.exists());
        }
        for invalid in ["", "0", "512", "8192", "65535", "65537", "-4096", "invalid"] {
            assert!(workspace_page_bytes(Some(invalid)).is_err());
        }
    }

    struct Payload(Vec<u8>);

    impl WorkPayload for Payload {
        fn resident_bytes(&self) -> usize {
            self.0.capacity()
        }

        fn encoded_len(&self) -> Result<u64, String> {
            u64::try_from(self.0.len()).map_err(|_| "test payload length overflow".into())
        }

        fn encode(&self, output: &mut dyn Write) -> Result<(), String> {
            // Writes cross both encoder and SQL chunk boundaries.
            for bytes in self.0.chunks(37_001) {
                output.write_all(bytes).map_err(|error| error.to_string())?;
            }
            Ok(())
        }
    }

    #[derive(Default)]
    struct SourceState {
        block: usize,
        claimed: bool,
        stored: bool,
        load: bool,
        load_credit: Option<Reservation>,
        retired: bool,
        restored: Vec<u8>,
        error: Option<String>,
        reclaim_request: usize,
        reclaim_pending: bool,
        reclaimable_bytes: usize,
        reclaim_calls: usize,
    }

    struct Source {
        payload: Arc<Payload>,
        state: Mutex<SourceState>,
    }

    impl WorkSource for Source {
        fn workspace_id(&self) -> u64 {
            1
        }

        fn pull(&self) -> Option<WorkBlock> {
            let mut state = self.state.lock().unwrap();
            if state.claimed || state.stored {
                return None;
            }
            state.claimed = true;
            Some(WorkBlock {
                block: state.block,
                payload: self.payload.clone(),
            })
        }

        fn next_load(&self) -> Option<WorkLoad> {
            let mut state = self.state.lock().unwrap();
            if !state.load {
                return None;
            }
            state.load = false;
            Some(WorkLoad {
                block: state.block,
                resident_bytes: self.payload.resident_bytes(),
                reservation: state.load_credit.take(),
            })
        }

        fn stored(&self, block: usize) -> Result<(), String> {
            let mut state = self.state.lock().unwrap();
            assert_eq!(block, state.block);
            state.stored = true;
            state.reclaimable_bytes = state.reclaimable_bytes.max(self.payload.resident_bytes());
            Ok(())
        }

        fn restore(&self, block: usize, input: &mut dyn Read) -> Result<(), String> {
            assert_eq!(block, self.state.lock().unwrap().block);
            let mut restored = Vec::new();
            input
                .read_to_end(&mut restored)
                .map_err(|error| error.to_string())?;
            self.state.lock().unwrap().restored = restored;
            Ok(())
        }

        fn fail(&self, error: String) {
            self.state.lock().unwrap().error = Some(error);
        }

        fn take_reclaim_request(&self) -> usize {
            std::mem::take(&mut self.state.lock().unwrap().reclaim_request)
        }

        fn reclaim_pending(&self) -> bool {
            let state = self.state.lock().unwrap();
            state.reclaim_pending || state.reclaimable_bytes != 0 || !state.stored
        }

        fn reclaim_stored(&self, byte_limit: usize) -> usize {
            let mut state = self.state.lock().unwrap();
            state.reclaim_calls += 1;
            let reclaimed = byte_limit.min(state.reclaimable_bytes);
            state.reclaimable_bytes -= reclaimed;
            reclaimed
        }

        fn is_retired(&self) -> bool {
            self.state.lock().unwrap().retired
        }
    }

    fn pragma_text(connection: &Connection, name: &str) -> String {
        connection
            .pragma_query_value(None, name, |row| row.get(0))
            .unwrap()
    }

    fn readback_credit(bytes: u64) -> Reservation {
        use crate::cache_source::budget::{Budget, BudgetConfig, MemorySnapshot};
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
            job_process_limit_bytes: None,
            job_memory_limit_bytes: None,
            job_headroom_bytes: None,
            job_headroom_peak_based: false,
            environment_limits_checked: true,
        };
        let budget = Budget::for_test(
            BudgetConfig {
                max_bytes: Some(10_000),
                warning_percent: 85,
                reclaim_percent: 70,
                system_headroom_bytes: 0,
                progress_reserve_bytes: 32,
                sample_interval: Duration::from_secs(60),
            },
            snapshot,
        )
        .unwrap();
        budget.try_reserve(bytes, ReservationKind::Compute).unwrap()
    }

    fn pending_readback(
        credit: Option<Reservation>,
        debt: usize,
        posted: usize,
    ) -> (WorkManager, Arc<Source>, Connection) {
        let data = vec![3, 5, 7, 9];
        let connection = chunk_database(&[(0, data.clone())]);
        let source = Arc::new(Source {
            payload: Arc::new(Payload(data)),
            state: Mutex::new(SourceState {
                block: 2,
                stored: true,
                load: true,
                load_credit: credit,
                reclaim_request: posted,
                reclaimable_bytes: 4,
                ..Default::default()
            }),
        });
        let mut manager = WorkManager::new(":memory:");
        manager.diagnostic_timing = Some(false);
        manager.register(source.clone());
        manager.sources[0].reclaim_remaining = debt;
        (manager, source, connection)
    }

    #[test]
    fn admitted_readback_cancels_only_old_finite_reclaim_debt() {
        for (debt, posted) in [(0, 23), (17, 0), (17, 23)] {
            let (mut manager, source, mut connection) =
                pending_readback(Some(readback_credit(4)), debt, posted);
            assert!(
                manager
                    .step_on_with_pressure(&mut connection, &mut |_| 0)
                    .unwrap()
            );
            assert_eq!(manager.sources[0].reclaim_remaining, 0);
            let state = source.state.lock().unwrap();
            assert_eq!(state.reclaim_request, 0);
            assert_eq!(state.reclaimable_bytes, 4);
            assert_eq!(state.reclaim_calls, 0);
            assert_eq!(state.restored, [3, 5, 7, 9]);
            assert_eq!(manager.report().loaded_blocks, 1);
        }
    }

    #[test]
    fn admitted_readback_preserves_max_drain_and_actual_memory_pressure() {
        for (debt, posted) in [(usize::MAX, 0), (17, usize::MAX), (usize::MAX, 23)] {
            let (mut manager, source, mut connection) =
                pending_readback(Some(readback_credit(4)), debt, posted);
            assert!(
                manager
                    .step_on_with_pressure(&mut connection, &mut |_| 0)
                    .unwrap()
            );
            assert_eq!(manager.sources[0].reclaim_remaining, usize::MAX);
            let state = source.state.lock().unwrap();
            assert_eq!(state.reclaimable_bytes, 0);
            assert_eq!(state.reclaim_calls, 1);
        }
        let (mut manager, source, mut connection) =
            pending_readback(Some(readback_credit(4)), 17, 23);
        assert!(
            manager
                .step_on_with_pressure(&mut connection, &mut |_| 4)
                .unwrap()
        );
        assert_eq!(manager.sources[0].reclaim_remaining, 0);
        let state = source.state.lock().unwrap();
        assert_eq!(state.reclaimable_bytes, 0);
        assert_eq!(state.reclaim_calls, 1);
    }

    #[test]
    fn absent_or_insufficient_readback_credit_does_not_cancel_reclaim() {
        let (mut manager, source, mut connection) = pending_readback(None, 17, 23);
        assert!(
            manager
                .step_on_with_pressure(&mut connection, &mut |_| 0)
                .unwrap()
        );
        assert_eq!(manager.sources[0].reclaim_remaining, 19);
        assert_eq!(source.state.lock().unwrap().reclaimable_bytes, 0);
        let (mut manager, source, mut connection) =
            pending_readback(Some(readback_credit(3)), 17, 23);
        let error = manager
            .step_on_with_pressure(&mut connection, &mut |_| 0)
            .unwrap_err();
        assert!(error.contains("exceeds its coordinator reservation"));
        assert_eq!(manager.sources[0].reclaim_remaining, 17);
        let state = source.state.lock().unwrap();
        assert_eq!(state.reclaim_request, 23);
        assert!(state.restored.is_empty());
        assert_eq!(manager.report().loaded_blocks, 0);
    }

    #[test]
    fn a_finished_readback_credit_does_not_cancel_reclaim() {
        let credit = readback_credit(4);
        credit.cancel_unused().unwrap();
        let (mut manager, source, mut connection) = pending_readback(Some(credit), 17, 23);
        let error = manager
            .step_on_with_pressure(&mut connection, &mut |_| 0)
            .unwrap_err();
        assert_eq!(error, ReservationError::Finished.to_string());
        assert_eq!(manager.sources[0].reclaim_remaining, 17);
        let state = source.state.lock().unwrap();
        assert_eq!(state.reclaim_request, 23);
        assert!(state.restored.is_empty());
        assert_eq!(manager.report().loaded_blocks, 0);
    }

    #[test]
    fn explicit_drain_survives_pins_and_future_saves_then_stops_evicting() {
        let source = Arc::new(Source {
            payload: Arc::new(Payload(Vec::new())),
            state: Mutex::new(SourceState {
                reclaim_request: usize::MAX,
                reclaim_pending: true,
                ..Default::default()
            }),
        });
        let mut entry = WorkEntry {
            source: source.clone(),
            reclaim_remaining: 0,
        };
        assert_eq!(entry.reclaim(0), 0);
        assert_eq!(entry.reclaim_remaining, usize::MAX);
        {
            let mut state = source.state.lock().unwrap();
            state.stored = true;
            state.reclaim_pending = false;
            state.reclaimable_bytes = 40;
        }
        assert_eq!(entry.reclaim(0), 40);
        assert_eq!(entry.reclaim(0), 0);
        assert_eq!(entry.reclaim_remaining, 0);
        {
            let mut state = source.state.lock().unwrap();
            state.reclaimable_bytes = 40;
            state.reclaim_calls = 0;
        }
        assert_eq!(entry.reclaim(0), 0);
        assert_eq!(source.state.lock().unwrap().reclaim_calls, 0);
        assert_eq!(source.state.lock().unwrap().reclaimable_bytes, 40);
        assert_eq!(entry.reclaim(7), 7);
        assert_eq!(source.state.lock().unwrap().reclaimable_bytes, 33);
    }

    #[test]
    fn finite_pressure_demand_stops_saving_and_leaves_future_blocks_resident() {
        let mut connection = chunk_database(&[]);
        let source = Arc::new(Source {
            payload: Arc::new(Payload(vec![1, 2, 3, 4])),
            state: Mutex::new(SourceState::default()),
        });
        let mut manager = WorkManager::new("unused-policy-test.sqlite");
        manager.register(source.clone());
        let mut normal = |_| 0;
        // A registered workspace alone does not request any new SQL writes.
        assert!(
            !manager
                .step_on_with_pressure(&mut connection, &mut normal)
                .unwrap()
        );
        assert!(!source.state.lock().unwrap().claimed);
        assert_eq!(manager.report().written_blocks, 0);

        source.state.lock().unwrap().reclaim_request = source.payload.resident_bytes();
        assert!(
            manager
                .step_on_with_pressure(&mut connection, &mut normal)
                .unwrap()
        );
        assert_eq!(manager.report().written_blocks, 1);
        assert_eq!(manager.sources[0].reclaim_remaining, 0);
        // A later sealed block is a new append-only identity, not an implicit
        // continuation of the earlier, now satisfied finite reclaim request.
        {
            let mut state = source.state.lock().unwrap();
            state.block = 1;
            state.claimed = false;
            state.stored = false;
        }
        for _ in 0..3 {
            assert!(
                !manager
                    .step_on_with_pressure(&mut connection, &mut normal)
                    .unwrap()
            );
        }
        assert!(!source.state.lock().unwrap().claimed);
        assert_eq!(manager.report().written_blocks, 1);
        let rows: i64 = connection
            .query_row(
                "SELECT count(*) FROM faber_work_blocks WHERE block=1",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(rows, 0);
    }

    #[test]
    fn native_pressure_relief_reclaims_saved_memory_before_saving_new_blocks() {
        let mut connection = chunk_database(&[]);
        let source = Arc::new(Source {
            payload: Arc::new(Payload(vec![1, 2, 3, 4])),
            state: Mutex::new(SourceState {
                reclaimable_bytes: 40,
                ..Default::default()
            }),
        });
        let mut manager = WorkManager::new("unused-policy-test.sqlite");
        manager.register(source.clone());
        let mut samples = Vec::new();
        let mut pressure = |force| {
            samples.push(force);
            if force { 0 } else { 40 }
        };
        assert!(
            manager
                .step_on_with_pressure(&mut connection, &mut pressure)
                .unwrap()
        );
        assert_eq!(samples, [false, true]);
        assert_eq!(source.state.lock().unwrap().reclaimable_bytes, 0);
        assert!(!source.state.lock().unwrap().claimed);
        assert_eq!(manager.report().written_blocks, 0);
    }

    fn pragma_integer(connection: &Connection, name: &str) -> i64 {
        connection
            .pragma_query_value(None, name, |row| row.get(0))
            .unwrap()
    }

    fn remove_test_file(path: &Path) {
        match fs::remove_file(path) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => panic!("remove test-owned file {}: {error}", path.display()),
        }
    }

    #[test]
    fn batch_candidates_are_explicit_and_the_default_remains_per_block() {
        assert_eq!(workspace_batch_policy(None).unwrap().bytes, 0);
        assert_eq!(workspace_batch_policy(Some("block")).unwrap().bytes, 0);
        let small = workspace_batch_policy(Some("256k")).unwrap();
        assert_eq!(small.bytes, 256 * 1024);
        assert_eq!(small.age, Duration::from_millis(10));
        let large = workspace_batch_policy(Some("1m")).unwrap();
        assert_eq!(large.bytes, 1024 * 1024);
        assert_eq!(large.age, Duration::from_millis(25));
        for invalid in ["", "0", "256", "all", "1M"] {
            assert!(workspace_batch_policy(Some(invalid)).is_err());
        }
    }

    #[test]
    fn an_acknowledged_uncommitted_block_reads_back_before_an_idle_commit() {
        let main_path = reserve_workspace_file(&std::env::temp_dir()).unwrap();
        let main = Connection::open(&main_path).unwrap();
        main.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL;")
            .unwrap();
        let mut database =
            WorkDatabase::open_with_settings(&std::env::temp_dir(), "exclusive", 65536).unwrap();
        let work_path = database.path.clone();
        let connection = database.connection.as_mut().unwrap();
        let source = Arc::new(Source {
            payload: Arc::new(Payload(vec![19; 7001])),
            state: Mutex::new(SourceState::default()),
        });
        let mut manager = WorkManager::new(main_path.to_str().unwrap());
        manager.diagnostic_timing = Some(false);
        // A controlled long age prevents scheduler pauses from turning this
        // pre-COMMIT visibility test into an ordinary committed read.
        manager.batch_policy = Some(BatchPolicy {
            bytes: 1024 * 1024,
            age: Duration::MAX,
        });
        manager.register(source.clone());
        assert!(
            manager
                .step_on_with_pressure(connection, &mut |_| usize::MAX)
                .unwrap()
        );
        assert!(source.state.lock().unwrap().stored);
        assert!(!connection.is_autocommit());
        assert_eq!(manager.report().commits, 0);
        assert_eq!(source.state.lock().unwrap().reclaimable_bytes, 0);
        source.state.lock().unwrap().load = true;
        assert!(
            manager
                .step_on_with_pressure(connection, &mut |_| 0)
                .unwrap()
        );
        assert_eq!(source.state.lock().unwrap().restored, source.payload.0);
        assert_eq!(manager.report().commits, 0);
        assert!(!connection.is_autocommit());
        let policy = workspace_batch_policy(Some("1m")).unwrap();
        assert!(policy.age > Duration::ZERO);
        manager.batch_policy = Some(policy);
        manager.transaction_age_override = Some(policy.age - Duration::from_nanos(1));
        assert!(!manager.commit_due(connection).unwrap());
        assert!(!connection.is_autocommit());
        manager.transaction_age_override = Some(policy.age);
        // No reads, writes or pressure requests: an idle turn still commits.
        assert!(
            manager
                .step_on_with_pressure(connection, &mut |_| 0)
                .unwrap()
        );
        assert!(connection.is_autocommit());
        assert_eq!(manager.report().commits, 1);
        assert!(manager.transaction_age_override.is_none());
        assert_eq!(pragma_text(&main, "journal_mode"), "wal");
        assert_eq!(pragma_integer(&main, "synchronous"), 2);
        assert_eq!(pragma_text(connection, "journal_mode"), "off");
        assert_eq!(pragma_integer(connection, "synchronous"), 0);
        drop(database);
        assert!(!work_path.exists());
        main.close().map_err(|(_, error)| error).unwrap();
        remove_test_file(&main_path);
    }

    #[test]
    fn byte_threshold_spans_steps_and_last_owner_discards_a_pending_batch() {
        let mut database =
            WorkDatabase::open_with_settings(&std::env::temp_dir(), "exclusive", 65536).unwrap();
        let path = database.path.clone();
        let connection = database.connection.as_mut().unwrap();
        let source = Arc::new(Source {
            payload: Arc::new(Payload(vec![7; CHUNK_BYTES / 2])),
            state: Mutex::new(SourceState::default()),
        });
        let mut manager = WorkManager::new(":memory:");
        manager.diagnostic_timing = Some(false);
        manager.batch_policy = Some(BatchPolicy {
            bytes: CHUNK_BYTES as u64,
            age: Duration::MAX,
        });
        manager.register(source.clone());
        for block in 0..3 {
            {
                let mut state = source.state.lock().unwrap();
                state.block = block;
                state.claimed = false;
                state.stored = false;
            }
            assert!(
                manager
                    .step_on_with_pressure(connection, &mut |_| 1)
                    .unwrap()
            );
            assert_eq!(manager.report().commits, block.div_ceil(2));
            assert_eq!(connection.is_autocommit(), block == 1);
        }
        let rows: i64 = connection
            .query_row("SELECT COUNT(*) FROM faber_work_blocks", [], |row| {
                row.get(0)
            })
            .unwrap();
        assert_eq!(rows, 3);
        source.state.lock().unwrap().retired = true;
        assert!(
            manager
                .step_on_with_pressure(connection, &mut |_| 0)
                .unwrap()
        );
        assert!(manager.sources.is_empty());
        assert_eq!(manager.report().commits, 1);
        // Retirement does not force COMMIT for a file which is being removed.
        assert!(!connection.is_autocommit());
        drop(database);
        assert!(!path.exists());
    }

    struct FailureWaiter {
        state: Mutex<Option<String>>,
        changed: std::sync::Condvar,
    }

    impl WorkSource for FailureWaiter {
        fn workspace_id(&self) -> u64 {
            2
        }
        fn pull(&self) -> Option<WorkBlock> {
            None
        }
        fn next_load(&self) -> Option<WorkLoad> {
            None
        }
        fn stored(&self, _: usize) -> Result<(), String> {
            unreachable!()
        }
        fn restore(&self, _: usize, _: &mut dyn Read) -> Result<(), String> {
            unreachable!()
        }
        fn fail(&self, error: String) {
            *self.state.lock().unwrap() = Some(error);
            self.changed.notify_all();
        }
        fn reclaim_stored(&self, _: usize) -> usize {
            0
        }
        fn is_retired(&self) -> bool {
            false
        }
    }

    #[test]
    fn a_failed_batch_wakes_every_owner_and_cannot_be_reused() {
        let mut database =
            WorkDatabase::open_with_settings(&std::env::temp_dir(), "exclusive", 65536).unwrap();
        let connection = database.connection.as_mut().unwrap();
        connection
            .execute_batch(
                "CREATE TRIGGER deny_second_block BEFORE INSERT ON faber_work_blocks
             WHEN NEW.block=1 BEGIN SELECT RAISE(ABORT,'injected workspace failure'); END;",
            )
            .unwrap();
        let source = Arc::new(Source {
            payload: Arc::new(Payload(vec![17; 100])),
            state: Mutex::new(SourceState::default()),
        });
        let waiting = Arc::new(FailureWaiter {
            state: Mutex::new(None),
            changed: std::sync::Condvar::new(),
        });
        let mut manager = WorkManager::new(":memory:");
        manager.diagnostic_timing = Some(false);
        manager.batch_policy = Some(BatchPolicy {
            bytes: 1024 * 1024,
            age: Duration::MAX,
        });
        manager.register(source.clone());
        manager.register(waiting.clone());
        assert!(
            manager
                .step_on_with_pressure(connection, &mut |_| 1)
                .unwrap()
        );
        assert!(source.state.lock().unwrap().stored);
        assert_eq!(manager.report().commits, 0);
        {
            let mut state = source.state.lock().unwrap();
            state.block = 1;
            state.claimed = false;
            state.stored = false;
        }
        let waiter = std::thread::spawn(move || {
            let state = waiting.state.lock().unwrap();
            let (state, _) = waiting
                .changed
                .wait_timeout_while(state, Duration::from_secs(5), |error| error.is_none())
                .unwrap();
            assert!(
                state
                    .as_ref()
                    .unwrap()
                    .contains("injected workspace failure")
            );
        });
        // The first turn visits the other owner; the next attempts block 1.
        assert!(
            !manager
                .step_on_with_pressure(connection, &mut |_| 1)
                .unwrap()
        );
        let error = manager
            .step_on_with_pressure(connection, &mut |_| 1)
            .unwrap_err();
        assert!(error.contains("injected workspace failure"));
        assert_eq!(
            source.state.lock().unwrap().error.as_deref(),
            Some(error.as_str())
        );
        waiter.join().unwrap();
        assert_eq!(manager.report().commits, 0);
        assert_eq!(
            manager
                .step_on_with_pressure(connection, &mut |_| 0)
                .unwrap_err(),
            error
        );
        // No rollback-and-continue: close and discard this private database.
        drop(database);
    }

    #[test]
    fn a_deferred_commit_failure_invalidates_previously_acknowledged_blocks() {
        let mut database =
            WorkDatabase::open_with_settings(&std::env::temp_dir(), "exclusive", 65536).unwrap();
        let connection = database.connection.as_mut().unwrap();
        connection
            .execute_batch(
                "PRAGMA foreign_keys=ON;
             CREATE TABLE test_parent (id INTEGER PRIMARY KEY);
             CREATE TABLE test_child (
                 parent INTEGER REFERENCES test_parent(id) DEFERRABLE INITIALLY DEFERRED
             );",
            )
            .unwrap();
        let source = Arc::new(Source {
            payload: Arc::new(Payload(vec![31; 100])),
            state: Mutex::new(SourceState::default()),
        });
        let waiting = Arc::new(FailureWaiter {
            state: Mutex::new(None),
            changed: std::sync::Condvar::new(),
        });
        let mut manager = WorkManager::new(":memory:");
        manager.diagnostic_timing = Some(false);
        manager.batch_policy = Some(BatchPolicy {
            bytes: 1024 * 1024,
            age: Duration::MAX,
        });
        manager.register(source.clone());
        manager.register(waiting.clone());
        assert!(
            manager
                .step_on_with_pressure(connection, &mut |_| 1)
                .unwrap()
        );
        assert!(source.state.lock().unwrap().stored);
        connection
            .execute("INSERT INTO test_child VALUES(404)", [])
            .unwrap();
        let policy = workspace_batch_policy(Some("1m")).unwrap();
        assert!(policy.age > Duration::ZERO);
        manager.batch_policy = Some(policy);
        manager.transaction_age_override = Some(policy.age);
        let error = manager
            .step_on_with_pressure(connection, &mut |_| 0)
            .unwrap_err();
        assert!(error.contains("FOREIGN KEY"));
        assert!(manager.transaction_age_override.is_none());
        assert_eq!(
            source.state.lock().unwrap().error.as_deref(),
            Some(error.as_str())
        );
        assert_eq!(
            waiting.state.lock().unwrap().as_deref(),
            Some(error.as_str())
        );
        assert_eq!(manager.report().commits, 0);
        assert_eq!(
            manager
                .step_on_with_pressure(connection, &mut |_| 0)
                .unwrap_err(),
            error
        );
        // A failed COMMIT cannot roll the batch back and resume the plans.
        drop(database);
    }

    #[test]
    fn large_save_probes_are_enabled_only_by_the_cached_diagnostic_flag() {
        for enabled in [false, true] {
            let source = Arc::new(Source {
                payload: Arc::new(Payload(vec![7; CHUNK_BYTES * 8 + 913])),
                state: Mutex::new(SourceState::default()),
            });
            let mut manager = WorkManager::new(":memory:");
            manager.diagnostic_timing = Some(enabled);
            manager.register(source.clone());
            let mut connection = chunk_database(&[]);
            assert!(
                manager
                    .step_on_with_pressure(&mut connection, &mut |_| 1)
                    .unwrap()
            );
            let report = manager.report();
            assert_eq!(report.written_blocks, 1);
            if enabled {
                assert!(report.save_sql_ns > 0);
                assert!(report.save_sql_ns <= report.work_save_ns);
                assert_eq!(report.save_large_blocks, 1);
                assert_eq!(report.save_large_bytes, source.payload.0.len() as u64);
                assert_eq!(report.save_large_ns, report.work_save_ns);
                assert!(report.cache_write_pages.is_some());
                assert!(report.cache_spill_pages.is_some());
            } else {
                assert_eq!(report.save_sql_ns, 0);
                assert_eq!(report.save_large_blocks, 0);
                assert_eq!(report.save_large_bytes, 0);
                assert_eq!(report.save_large_ns, 0);
                assert!(report.cache_write_pages.is_none());
                assert!(report.cache_spill_pages.is_none());
            }
        }
    }

    #[test]
    fn streaming_multiple_rows_uses_an_independent_database_and_cleans_only_its_file() {
        let directory = std::env::temp_dir();
        let main_path = reserve_workspace_file(&directory).unwrap();
        let sentinel = reserve_workspace_file(&directory).unwrap();
        fs::write(&sentinel, b"unrelated existing file").unwrap();
        let mut main = Connection::open(&main_path).unwrap();
        main.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL;")
            .unwrap();
        let main_journal = pragma_text(&main, "journal_mode");
        let main_synchronous = pragma_integer(&main, "synchronous");
        let source = Arc::new(Source {
            payload: Arc::new(Payload(
                (0..CHUNK_BYTES * 2 + 913)
                    .map(|index| (index % 251) as u8)
                    .collect(),
            )),
            state: Mutex::new(SourceState {
                reclaim_request: usize::MAX,
                ..Default::default()
            }),
        });
        let mut manager = WorkManager::new(main_path.to_str().unwrap());
        manager.diagnostic_timing = Some(true);
        assert!(manager.database.is_none());
        // Supply the same pure directory resolution explicitly, so this unit
        // test never changes or depends on the process-wide override variable.
        let selected = workspace_directory(&main_path, None).unwrap();
        assert_eq!(selected, main_path.parent().unwrap());
        manager.database = Some(WorkDatabase::open(&selected).unwrap());
        manager.register(source.clone());
        assert!(manager.step(&mut main).unwrap());
        let work = manager.database.as_ref().unwrap();
        let work_path = work.path.clone();
        assert_ne!(work_path, main_path);
        let connection = work.connection.as_ref().unwrap();
        let chunks: i64 = connection
            .query_row("SELECT COUNT(*) FROM faber_work_blocks", [], |row| {
                row.get(0)
            })
            .unwrap();
        assert_eq!(chunks, 3);
        let largest: i64 = connection
            .query_row(
                "SELECT MAX(length(data)) FROM faber_work_blocks",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(largest, CHUNK_BYTES as i64);
        assert_eq!(pragma_text(connection, "journal_mode"), "off");
        assert_eq!(pragma_integer(connection, "synchronous"), 0);
        assert_eq!(pragma_integer(connection, "cache_size"), -2048);
        source.state.lock().unwrap().load = true;
        assert!(manager.step(&mut main).unwrap());
        assert_eq!(source.state.lock().unwrap().restored, source.payload.0);
        let report = manager.report();
        assert_eq!(report.written_blocks, 1);
        assert_eq!(report.commits, 1);
        assert_eq!(report.loaded_blocks, 1);
        assert_eq!(
            report.max_block_encoded_bytes,
            source.payload.0.len() as u64
        );
        assert_eq!(
            report.max_loaded_block_resident_bytes,
            source.payload.0.len()
        );
        assert!(report.work_save_ns > 0);
        assert!(report.save_sql_ns > 0);
        assert!(report.save_sql_ns <= report.work_save_ns);
        assert_eq!(report.save_large_blocks, 0);
        assert_eq!(report.save_large_bytes, 0);
        assert_eq!(report.save_large_ns, 0);
        assert!(report.cache_write_pages.unwrap() > 0);
        assert!(report.cache_spill_pages.is_some());
        assert!(report.work_load_ns > 0);
        assert!(report.work_read_sql_ns > 0);
        assert!(report.work_read_sql_ns <= report.work_load_ns);
        assert!(report.work_reclaim_ns > 0);
        assert_eq!(pragma_text(&main, "journal_mode"), main_journal);
        assert_eq!(pragma_integer(&main, "synchronous"), main_synchronous);
        assert_eq!(main_journal, "wal");
        assert_eq!(main_synchronous, 2);
        let main_work_tables: i64 = main
            .query_row(
                "SELECT COUNT(*) FROM sqlite_master WHERE name='faber_work_blocks'",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(main_work_tables, 0);
        source.state.lock().unwrap().retired = true;
        assert!(manager.step(&mut main).unwrap());
        assert!(manager.database.is_none());
        assert!(manager.scratch.is_none());
        assert!(!work_path.exists());
        drop(manager);
        assert!(!work_path.exists());
        assert!(main_path.exists());
        assert_eq!(fs::read(&sentinel).unwrap(), b"unrelated existing file");
        main.close().map_err(|(_, error)| error).unwrap();
        remove_test_file(&main_path);
        remove_test_file(&sentinel);
    }

    #[test]
    fn memory_cache_requires_an_explicit_directory_and_constructor_is_lazy() {
        let mut manager = WorkManager::new(":memory:");
        assert!(manager.database.is_none());
        assert!(manager.initialization_error.is_none());
        let mut main = Connection::open_in_memory().unwrap();
        assert!(!manager.step(&mut main).unwrap());
        assert!(manager.database.is_none());
        assert!(workspace_directory(Path::new(":memory:"), None).is_err());
        let explicit = std::env::temp_dir();
        assert_eq!(
            workspace_directory(Path::new(":memory:"), Some(explicit.as_os_str())).unwrap(),
            explicit
        );
        assert!(workspace_directory(Path::new("cache.sqlite"), Some(OsStr::new(""))).is_err());
        assert_eq!(
            workspace_directory(Path::new("cache.sqlite"), None).unwrap(),
            Path::new(".")
        );
    }

    #[test]
    fn failed_file_reservation_leaves_existing_paths_untouched() {
        let sentinel = reserve_workspace_file(&std::env::temp_dir()).unwrap();
        fs::write(&sentinel, b"must remain").unwrap();
        assert!(WorkDatabase::open(&sentinel).is_err());
        assert_eq!(fs::read(&sentinel).unwrap(), b"must remain");
        remove_test_file(&sentinel);
    }
}
