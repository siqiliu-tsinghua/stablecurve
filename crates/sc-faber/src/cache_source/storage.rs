//! One SQLite owner, pulling completed values from their mathematical sources.
//!
//! Ordinary mathematical returns do not wait for this service. Explicit flush
//! captures finite per-source completion targets when its request is handled;
//! it is not an atomic snapshot of concurrently running queries. Only successful
//! transactions advance durable cursors. Drop stops without a durability claim.

use super::budget::{Pressure, process_budget};
use super::codec::encode_key_into;
use super::compat::LegacyReader;
use super::work::WorkSource;
use super::work_storage::WorkManager;
use super::{IntegralKey, LookupBatch, Registry, ResultSource};
use rusqlite::{Connection, OptionalExtension, params};
#[cfg(test)]
use sc_num::Zero;
use sc_num::{Q, decode_rational, encode_rational_into};
use std::any::Any;
use std::fmt;
use std::ops::Range;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock, mpsc};
use std::thread::{self, JoinHandle, Thread};
use std::time::{Duration, Instant};

const COPY_WINDOW_BYTES: usize = 64 * 1024;
const TRANSACTION_BYTES: usize = 8 * 1024 * 1024;
const TRANSACTION_ROWS: usize = 16_384;
const TRANSACTION_AGE: Duration = Duration::from_millis(100);
const IDLE_PAUSE: Duration = Duration::from_millis(1);
const BUSY_PAUSE: Duration = Duration::from_micros(50);
const HISTORY_KEY_LIMIT: usize = 256;
const HISTORY_KEY_BYTES: usize = 16 * 1024;

fn needs_descendant_history(connection: &Connection, legacy_rows: bool) -> Result<bool, String> {
    if legacy_rows {
        return Ok(true);
    }
    // Inspect only a small opening snapshot. The CASE avoids materializing a
    // huge/unknown key merely to decide that recovery must remain enabled.
    let mut statement = connection
        .prepare(
            "SELECT CASE WHEN length(k)<=?1 THEN k ELSE NULL END
         FROM faber_integrals ORDER BY k LIMIT ?2",
        )
        .map_err(sql_error)?;
    let mut rows = statement
        .query(params![
            HISTORY_KEY_BYTES as i64,
            (HISTORY_KEY_LIMIT + 1) as i64
        ])
        .map_err(sql_error)?;
    let mut count = 0;
    while let Some(row) = rows.next().map_err(sql_error)? {
        if count == HISTORY_KEY_LIMIT {
            return Ok(true);
        }
        count += 1;
        let rusqlite::types::ValueRef::Blob(key) = row.get_ref(0).map_err(sql_error)? else {
            return Ok(true);
        };
        if super::codec::descendant_needs_dag(key) != Some(false) {
            return Ok(true);
        }
    }
    Ok(false)
}

#[derive(Clone, Default)]
pub(crate) struct StorageReport {
    pub(crate) rows: usize,
    /// Committed key and rational bytes, excluding row descriptors and capacity.
    pub(crate) bytes: usize,
    pub(crate) commits: usize,
    pub(crate) pulls: usize,
    pub(crate) sources: usize,
    pub(crate) busy_pulls: usize,
    pub(crate) empty_pulls: usize,
    /// Live encoded bytes, row descriptors and the reusable key scratch buffer.
    pub(crate) max_buffer_bytes: usize,
    pub(crate) max_buffer_capacity_bytes: usize,
    /// Backend wall time in bounded source encoding windows, including identity
    /// preparation. This includes the pull now that original Q values are borrowed.
    pub(crate) encode_ns: u128,
    /// Backend wall time in transactions, including SQLite commit.
    pub(crate) commit_ns: u128,
    pub(crate) insert_ns: u128,
    pub(crate) transaction_commit_ns: u128,
    pub(crate) sort_ns: u128,
    pub(crate) read_keys: usize,
    pub(crate) read_hits: usize,
    /// Maximum native private/footprint observation; sampling can miss peaks.
    pub(crate) peak_sampled_private_bytes: u64,
    pub(crate) peak_native_commit_bytes: Option<u64>,
    pub(crate) last_available_bytes: Option<u64>,
    pub(crate) last_memory_limit_bytes: Option<u64>,
    pub(crate) last_warning_bytes: Option<u64>,
    pub(crate) last_cgroup_limit_bytes: Option<u64>,
    pub(crate) last_reclaim_target_bytes: u64,
    pub(crate) last_progress_at_risk: bool,
    pub(crate) last_active_reserved_bytes: u64,
    pub(crate) last_settled_reserved_bytes: u64,
    pub(crate) last_sample_generation: u64,
    pub(crate) reclaimed_bytes: usize,
    /// GNU allocator page-return attempts by the pressure worker, not credit.
    pub(crate) allocator_trim_calls: usize,
    pub(crate) allocator_trim_ns: u128,
    pub(crate) work_written_blocks: usize,
    pub(crate) work_written_bytes: u64,
    pub(crate) work_commits: usize,
    pub(crate) work_loaded_blocks: usize,
    pub(crate) work_loaded_bytes: u64,
    /// Backend block timings are separate from persistent integral SQL time.
    pub(crate) work_save_ns: u128,
    /// Optional detailed workspace probes; populated only in report mode.
    pub(crate) work_save_sql_ns: u128,
    pub(crate) work_save_large_blocks: usize,
    pub(crate) work_save_large_bytes: u64,
    pub(crate) work_save_large_ns: u128,
    pub(crate) work_cache_write_pages: Option<u64>,
    pub(crate) work_cache_spill_pages: Option<u64>,
    pub(crate) work_load_ns: u128,
    /// Workspace SELECT setup/reset, row steps, column access and blob copy.
    pub(crate) work_read_sql_ns: u128,
    pub(crate) work_reclaim_ns: u128,
    pub(crate) work_max_loaded_block_resident_bytes: usize,
    pub(crate) work_max_block_encoded_bytes: u64,
}

impl fmt::Debug for StorageReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StorageReport")
            .field("rows", &self.rows)
            .field("bytes", &self.bytes)
            .field("commits", &self.commits)
            .field("pulls", &self.pulls)
            .field("sources", &self.sources)
            .field("busy_pulls", &self.busy_pulls)
            .field("empty_pulls", &self.empty_pulls)
            .field("max_buffer_bytes", &self.max_buffer_bytes)
            .field("max_buffer_capacity_bytes", &self.max_buffer_capacity_bytes)
            .field("encode_ns", &self.encode_ns)
            .field("commit_ns", &self.commit_ns)
            .field("insert_ns", &self.insert_ns)
            .field("transaction_commit_ns", &self.transaction_commit_ns)
            .field("sort_ns", &self.sort_ns)
            .field("read_keys", &self.read_keys)
            .field("read_hits", &self.read_hits)
            .field(
                "peak_sampled_private_bytes",
                &self.peak_sampled_private_bytes,
            )
            .field("peak_native_commit_bytes", &self.peak_native_commit_bytes)
            .field("last_available_bytes", &self.last_available_bytes)
            .field("last_memory_limit_bytes", &self.last_memory_limit_bytes)
            .field("last_warning_bytes", &self.last_warning_bytes)
            .field("last_cgroup_limit_bytes", &self.last_cgroup_limit_bytes)
            .field("last_reclaim_target_bytes", &self.last_reclaim_target_bytes)
            .field("last_progress_at_risk", &self.last_progress_at_risk)
            .field(
                "last_active_reserved_bytes",
                &self.last_active_reserved_bytes,
            )
            .field(
                "last_settled_reserved_bytes",
                &self.last_settled_reserved_bytes,
            )
            .field("last_sample_generation", &self.last_sample_generation)
            .field("reclaimed_bytes", &self.reclaimed_bytes)
            .field("allocator_trim_calls", &self.allocator_trim_calls)
            .field("allocator_trim_ns", &self.allocator_trim_ns)
            .field("work_written_blocks", &self.work_written_blocks)
            .field("work_written_bytes", &self.work_written_bytes)
            .field("work_commits", &self.work_commits)
            .field("work_loaded_blocks", &self.work_loaded_blocks)
            .field("work_loaded_bytes", &self.work_loaded_bytes)
            .field("work_save_ns", &self.work_save_ns)
            .field("work_save_sql_ns", &self.work_save_sql_ns)
            .field("work_save_large_blocks", &self.work_save_large_blocks)
            .field("work_save_large_bytes", &self.work_save_large_bytes)
            .field("work_save_large_ns", &self.work_save_large_ns)
            .field("work_cache_write_pages", &self.work_cache_write_pages)
            .field("work_cache_spill_pages", &self.work_cache_spill_pages)
            .field("work_load_ns", &self.work_load_ns)
            .field("work_read_sql_ns", &self.work_read_sql_ns)
            .field("work_reclaim_ns", &self.work_reclaim_ns)
            .field(
                "work_max_loaded_block_resident_bytes",
                &self.work_max_loaded_block_resident_bytes,
            )
            .field(
                "work_max_block_encoded_bytes",
                &self.work_max_block_encoded_bytes,
            )
            .finish()
    }
}

type ReportReply = mpsc::Sender<Result<StorageReport, String>>;
type ValueReply = mpsc::Sender<Result<Option<Q>, String>>;
type BatchReply = mpsc::Sender<Result<Vec<Option<Q>>, String>>;

enum Request {
    Flush(ReportReply),
    Close(ReportReply),
    Get(IntegralKey, ValueReply),
    GetBatch(Arc<dyn LookupBatch>, BatchReply),
    DescendantHistory(mpsc::Sender<Result<bool, String>>),
    Workspace(Arc<dyn WorkSource>, mpsc::Sender<Result<(), String>>),
    #[cfg(test)]
    LegacyGet(Vec<u8>, ValueReply),
    Stop,
}

impl Request {
    fn fail(self, error: &str) {
        match self {
            Self::Flush(reply) | Self::Close(reply) => {
                let _ = reply.send(Err(error.to_owned()));
            }
            Self::Get(_, reply) => {
                let _ = reply.send(Err(error.to_owned()));
            }
            Self::GetBatch(_, reply) => {
                let _ = reply.send(Err(error.to_owned()));
            }
            Self::DescendantHistory(reply) => {
                let _ = reply.send(Err(error.to_owned()));
            }
            Self::Workspace(source, reply) => {
                source.fail(error.to_owned());
                let _ = reply.send(Err(error.to_owned()));
            }
            #[cfg(test)]
            Self::LegacyGet(_, reply) => {
                let _ = reply.send(Err(error.to_owned()));
            }
            Self::Stop => {}
        }
    }
}

#[derive(Default)]
struct Control {
    stop: AtomicBool,
    failed: AtomicBool,
    error: Mutex<Option<String>>,
    worker: OnceLock<Thread>,
}

impl Control {
    fn wake(&self) {
        if let Some(worker) = self.worker.get() {
            worker.unpark();
        }
    }

    fn error(&self) -> Option<String> {
        if !self.failed.load(Ordering::Acquire) {
            return None;
        }
        self.error
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clone()
    }

    fn fail(&self, error: String) {
        let mut stored = self.error.lock().unwrap_or_else(|error| error.into_inner());
        if stored.is_none() {
            *stored = Some(error);
        }
        drop(stored);
        self.failed.store(true, Ordering::Release);
    }
}

#[derive(Clone)]
pub(crate) struct StorageHandle {
    requests: mpsc::Sender<Request>,
    control: Arc<Control>,
}

impl StorageHandle {
    pub(crate) fn error(&self) -> Option<String> {
        self.control.error()
    }

    pub(crate) fn flush(&self) -> Result<StorageReport, String> {
        let (reply, receive) = mpsc::channel();
        self.send(Request::Flush(reply))?;
        self.receive(receive)
    }

    /// Read a committed SQLite value. Resident mathematical L1 lookups must
    /// not send requests; this handle remains bound to its original service.
    pub(crate) fn get(&self, key: IntegralKey) -> Result<Option<Q>, String> {
        let (reply, receive) = mpsc::channel();
        self.send(Request::Get(key, reply))?;
        self.receive(receive)
    }

    pub(crate) fn lookup_batch(
        &self,
        batch: Arc<dyn LookupBatch>,
    ) -> Result<Vec<Option<Q>>, String> {
        let (reply, receive) = mpsc::channel();
        self.send(Request::GetBatch(batch, reply))?;
        self.receive(receive)
    }

    /// Query-entry refresh only when generic history exists but its bounded
    /// opening snapshot contained no possible DAG descendants. L1 hits never
    /// need this request; own commits do not invalidate the opening token.
    pub(crate) fn descendant_history(&self) -> Result<bool, String> {
        let (reply, receive) = mpsc::channel();
        self.send(Request::DescendantHistory(reply))?;
        self.receive(receive)
    }

    pub(crate) fn register_workspace(&self, source: Arc<dyn WorkSource>) -> Result<(), String> {
        let (reply, receive) = mpsc::channel();
        self.send(Request::Workspace(source, reply))?;
        // A successful send alone is insufficient: a failing worker may drop
        // a queued registration after its final drain. The explicit reply
        // proves that the service owns this source and will wake its waiters.
        match receive.recv_timeout(Duration::from_secs(30)) {
            Ok(result) => result,
            Err(error) => Err(self
                .error()
                .unwrap_or_else(|| format!("plan workspace registration: {error}"))),
        }
    }

    #[cfg(test)]
    pub(crate) fn legacy_get(&self, key: Vec<u8>) -> Result<Option<Q>, String> {
        let (reply, receive) = mpsc::channel();
        self.send(Request::LegacyGet(key, reply))?;
        self.receive(receive)
    }

    fn send(&self, request: Request) -> Result<(), String> {
        if let Some(error) = self.error() {
            return Err(error);
        }
        self.requests.send(request).map_err(|_| {
            self.error()
                .unwrap_or_else(|| "cache storage thread is unavailable".to_owned())
        })?;
        self.control.wake();
        Ok(())
    }

    fn receive<T>(&self, receive: mpsc::Receiver<Result<T, String>>) -> Result<T, String> {
        receive.recv().unwrap_or_else(|_| {
            Err(self
                .error()
                .unwrap_or_else(|| "cache storage request was interrupted".to_owned()))
        })
    }
}

pub(crate) struct Storage {
    handle: StorageHandle,
    worker: Option<JoinHandle<()>>,
}

impl Storage {
    pub(crate) fn open(path: &str, registry: Arc<Registry>) -> Result<Self, String> {
        let (requests, receive) = mpsc::channel();
        let (ready, opened) = mpsc::channel();
        let control = Arc::new(Control::default());
        let worker_control = Arc::clone(&control);
        let path = path.to_owned();
        let worker = thread::Builder::new()
            .name("sc-faber-storage".into())
            .spawn(move || worker_main(&path, registry, receive, ready, worker_control))
            .map_err(|error| format!("start cache storage thread: {error}"))?;
        match opened.recv() {
            Ok(Ok(())) => Ok(Self {
                handle: StorageHandle { requests, control },
                worker: Some(worker),
            }),
            response => {
                let _ = worker.join();
                Err(match response {
                    Ok(Err(error)) => error,
                    Err(_) => control
                        .error()
                        .unwrap_or_else(|| "cache storage stopped while opening".to_owned()),
                    Ok(Ok(())) => unreachable!(),
                })
            }
        }
    }

    pub(crate) fn error(&self) -> Option<String> {
        self.handle.error()
    }

    pub(crate) fn handle(&self) -> StorageHandle {
        self.handle.clone()
    }

    #[cfg(test)]
    pub(crate) fn flush(&self) -> Result<StorageReport, String> {
        let (reply, receive) = mpsc::channel();
        self.handle.send(Request::Flush(reply))?;
        self.handle.receive(receive)
    }

    pub(crate) fn close(mut self) -> Result<StorageReport, String> {
        let (reply, receive) = mpsc::channel();
        let result = self
            .handle
            .send(Request::Close(reply))
            .and_then(|()| self.handle.receive(receive));
        if let Some(worker) = self.worker.take() {
            if worker.join().is_err() {
                self.handle
                    .control
                    .fail("cache storage thread panicked while closing".to_owned());
            }
        }
        if let Some(error) = self.error() {
            Err(error)
        } else {
            result
        }
    }

    /// Read a committed SQLite value. Mathematical L1 misses may use this;
    /// ordinary resident lookups must not send storage requests.
    #[cfg(test)]
    pub(crate) fn get(&self, key: IntegralKey) -> Result<Option<Q>, String> {
        self.handle.get(key)
    }

    #[cfg(test)]
    pub(crate) fn legacy_get(&self, key: Vec<u8>) -> Result<Option<Q>, String> {
        self.handle.legacy_get(key)
    }
}

impl Drop for Storage {
    fn drop(&mut self) {
        self.handle.control.stop.store(true, Ordering::Release);
        let _ = self.handle.requests.send(Request::Stop);
        self.handle.control.wake();
        if let Some(worker) = self.worker.take() {
            // Drop may run during unwinding. Explicit close is the durability
            // boundary; Drop merely releases the worker and its connection.
            let _ = worker.join();
        }
    }
}

struct SourceCursor {
    source: Arc<dyn ResultSource>,
    read: usize,
    durable: usize,
    pending: bool,
}

struct Barrier {
    targets: Vec<usize>,
    reply: ReportReply,
    close: bool,
}

struct EncodedRow {
    key: Range<usize>,
    value: Range<usize>,
}

#[derive(Default)]
struct Buffer {
    bytes: Vec<u8>,
    rows: Vec<EncodedRow>,
    sources: Vec<usize>,
    first: Option<Instant>,
}

struct DescendantBounds {
    first: Option<Vec<u8>>,
    last: Option<Vec<u8>>,
}

impl DescendantBounds {
    fn read(connection: &Connection) -> Result<Self, String> {
        let endpoint = |order| {
            connection
                .query_row(order, [], |row| row.get::<_, Vec<u8>>(0))
                .optional()
                .map_err(sql_error)
        };
        Ok(Self {
            first: endpoint(
                "SELECT k FROM faber_integrals WHERE k>=x'0104' AND k<x'0105' ORDER BY k ASC LIMIT 1",
            )?,
            last: endpoint(
                "SELECT k FROM faber_integrals WHERE k>=x'0104' AND k<x'0105' ORDER BY k DESC LIMIT 1",
            )?,
        })
    }

    fn eligible(key: &IntegralKey, encoded: &[u8]) -> bool {
        let IntegralKey::Descendant { ch, .. } = key else {
            return false;
        };
        encoded.starts_with(&[1, 4]) && ch.iter().any(|&(m, p)| m != 0 && p != 0)
    }

    fn excludes(&self, key: &IntegralKey, encoded: &[u8]) -> bool {
        if !Self::eligible(key, encoded) {
            return false;
        }
        self.first
            .as_ref()
            .is_none_or(|first| encoded < first.as_slice())
            || self
                .last
                .as_ref()
                .is_none_or(|last| encoded > last.as_slice())
    }
}

struct Worker {
    connection: Connection,
    legacy: LegacyReader,
    registry: Arc<Registry>,
    control: Arc<Control>,
    sources: Vec<SourceCursor>,
    next_source: usize,
    barriers: Vec<Barrier>,
    closing: bool,
    buffer: Buffer,
    key: Vec<u8>,
    report: StorageReport,
    sort_writes: bool,
    transaction_rows: usize,
    transaction_age: Duration,
    workspaces: WorkManager,
    last_memory_poll: Instant,
    descendant_history_version: Option<i64>,
}

impl Worker {
    fn open(path: &str, registry: Arc<Registry>, control: Arc<Control>) -> Result<Self, String> {
        let connection = Connection::open(path).map_err(|error| format!("open cache: {error}"))?;
        connection
            .busy_timeout(Duration::from_millis(5000))
            .map_err(sql_error)?;
        connection
            .pragma_update(None, "journal_mode", "WAL")
            .map_err(sql_error)?;
        connection
            .pragma_update(None, "synchronous", "FULL")
            .map_err(sql_error)?;
        // Diagnostic controls are read once by the sole SQLite owner. They
        // never change the arithmetic path or the durability policy.
        if let Ok(value) = std::env::var("SC_FABER_CACHE_MEASURE_PAGE_CACHE_MIB") {
            let mib: i64 = value
                .parse()
                .ok()
                .filter(|mib| (1..=512).contains(mib))
                .ok_or("SC_FABER_CACHE_MEASURE_PAGE_CACHE_MIB must be in 1..=512")?;
            connection
                .pragma_update(None, "cache_size", -mib * 1024)
                .map_err(sql_error)?;
        }
        let sort_writes = std::env::var("SC_FABER_CACHE_MEASURE_SORT_WRITES").as_deref() == Ok("1");
        let transaction_rows =
            diagnostic_usize("SC_FABER_CACHE_MEASURE_TRANSACTION_ROWS", TRANSACTION_ROWS)?;
        let transaction_ms = diagnostic_usize(
            "SC_FABER_CACHE_MEASURE_TRANSACTION_MS",
            TRANSACTION_AGE.as_millis() as usize,
        )?;
        if transaction_rows == 0 || transaction_ms == 0 {
            return Err("diagnostic transaction row/time thresholds must be positive".into());
        }
        let wal_pages = diagnostic_usize("SC_FABER_CACHE_MEASURE_WAL_PAGES", 1000)?;
        let wal_pages = i64::try_from(wal_pages).map_err(|_| "WAL page threshold exceeds i64")?;
        connection
            .pragma_update(None, "wal_autocheckpoint", wal_pages)
            .map_err(sql_error)?;
        let transaction_age = Duration::from_millis(
            u64::try_from(transaction_ms).map_err(|_| "transaction milliseconds exceed u64")?,
        );
        connection
            .execute_batch(
                "CREATE TABLE IF NOT EXISTS st (
                    k BLOB PRIMARY KEY, v TEXT NOT NULL
                ) WITHOUT ROWID;
                CREATE TABLE IF NOT EXISTS faber_integrals (
                    k BLOB PRIMARY KEY, v BLOB NOT NULL
                ) WITHOUT ROWID;",
            )
            .map_err(sql_error)?;
        // Capture before opening the read snapshot: a commit between this
        // token and that snapshot must conservatively invalidate a false hint.
        let opening_version = connection
            .pragma_query_value(None, "data_version", |row| row.get::<_, i64>(0))
            .ok();
        let snapshot = connection.unchecked_transaction().map_err(sql_error)?;
        let typed_history: bool = snapshot
            .query_row("SELECT EXISTS(SELECT 1 FROM faber_integrals)", [], |row| {
                row.get(0)
            })
            .map_err(sql_error)?;
        let st_history: bool = snapshot
            .query_row("SELECT EXISTS(SELECT 1 FROM st)", [], |row| row.get(0))
            .map_err(sql_error)?;
        let has_old_table: bool = snapshot.query_row(
            "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='faber_values')", [], |row| row.get(0),
        ).map_err(sql_error)?;
        let has_old_rows = has_old_table
            && snapshot
                .query_row(
                    "SELECT EXISTS(SELECT 1 FROM faber_values LIMIT 1)",
                    [],
                    |row| row.get::<_, bool>(0),
                )
                .map_err(sql_error)?;
        let history = typed_history || st_history || has_old_rows;
        let descendant_history = opening_version.is_none()
            || needs_descendant_history(&snapshot, st_history || has_old_rows).unwrap_or(true);
        let legacy = LegacyReader::new(&snapshot)?;
        snapshot.commit().map_err(sql_error)?;
        if history {
            // This is an open-time property. Results written in the current
            // session do not turn a new empty cache into a per-node L2 probe.
            registry.set_history(true);
            registry.set_descendant_history(descendant_history);
        }
        Ok(Self {
            connection,
            legacy,
            registry,
            control,
            sources: Vec::new(),
            next_source: 0,
            barriers: Vec::new(),
            closing: false,
            buffer: Buffer::default(),
            key: Vec::new(),
            report: StorageReport::default(),
            sort_writes,
            transaction_rows,
            transaction_age,
            workspaces: WorkManager::new(path),
            last_memory_poll: Instant::now(),
            descendant_history_version: (history && !descendant_history)
                .then_some(opening_version)
                .flatten(),
        })
    }

    fn discover(&mut self) {
        self.sources.extend(
            self.registry
                .sources_since(self.sources.len())
                .into_iter()
                .map(|source| {
                    let floor = source.floor();
                    SourceCursor {
                        source,
                        read: floor,
                        durable: floor,
                        pending: false,
                    }
                }),
        );
        self.report.sources = self.sources.len();
    }

    fn capture(&mut self, reply: ReportReply, close: bool) {
        self.discover();
        let targets = self
            .sources
            .iter()
            .map(|source| source.source.ready())
            .collect();
        self.barriers.push(Barrier {
            targets,
            reply,
            close,
        });
        self.closing |= close;
    }

    fn request(&mut self, request: Request) -> Result<(), String> {
        if self.closing {
            request.fail("cache storage is closing");
            return Ok(());
        }
        match request {
            Request::Flush(reply) => self.capture(reply, false),
            Request::Close(reply) => self.capture(reply, true),
            Request::Get(key, reply) => {
                let result = self.get(&key);
                let error = result.as_ref().err().cloned();
                let _ = reply.send(result);
                if let Some(error) = error {
                    return Err(error);
                }
            }
            Request::GetBatch(batch, reply) => {
                let result = self.get_batch(batch.as_ref());
                let error = result.as_ref().err().cloned();
                let _ = reply.send(result);
                if let Some(error) = error {
                    return Err(error);
                }
            }
            Request::DescendantHistory(reply) => {
                let _ = reply.send(Ok(self.descendant_history()));
            }
            Request::Workspace(source, reply) => {
                self.workspaces.register(source);
                let _ = reply.send(Ok(()));
            }
            #[cfg(test)]
            Request::LegacyGet(key, reply) => {
                let result = self.legacy_get(&key);
                let error = result.as_ref().err().cloned();
                let _ = reply.send(result);
                if let Some(error) = error {
                    return Err(error);
                }
            }
            Request::Stop => self.control.stop.store(true, Ordering::Release),
        }
        Ok(())
    }

    fn get(&mut self, key: &IntegralKey) -> Result<Option<Q>, String> {
        self.get_with_bounds(key, None)
    }

    fn descendant_history(&mut self) -> bool {
        if self.registry.has_history() && !self.registry.has_descendant_history() {
            let current = self
                .connection
                .pragma_query_value(None, "data_version", |row| row.get::<_, i64>(0))
                .ok();
            if current.is_none() || current != self.descendant_history_version {
                self.registry.set_history(true);
                self.descendant_history_version = None;
            }
        }
        self.registry.has_descendant_history()
    }

    fn get_batch(&mut self, batch: &dyn LookupBatch) -> Result<Vec<Option<Q>>, String> {
        if !batch.descendants_only() || batch.is_empty() {
            return (0..batch.len())
                .map(|index| self.get(&batch.key(index)))
                .collect();
        }
        self.connection
            .execute_batch("BEGIN DEFERRED")
            .map_err(sql_error)?;
        let result: Result<Vec<Option<Q>>, String> = {
            // These constant-size hints belong only to this read snapshot.
            // Probe errors disable the optimization: exact reads retain their
            // original compatibility/decoding error priority.
            let bounds = match self
                .legacy
                .has_positive_descendant_aliases(&self.connection)
            {
                Ok(false) => DescendantBounds::read(&self.connection).ok(),
                Ok(true) | Err(_) => None,
            };
            (0..batch.len())
                .map(|index| {
                    let key = batch.key(index);
                    if !matches!(key, IntegralKey::Descendant { .. }) {
                        return Err(
                            "descendant-only lookup batch contains another key namespace".into(),
                        );
                    }
                    self.get_with_bounds(&key, bounds.as_ref())
                })
                .collect()
        };
        let ended = self
            .connection
            .execute_batch(if result.is_ok() { "COMMIT" } else { "ROLLBACK" })
            .map_err(sql_error);
        match result {
            Ok(values) => ended.map(|()| values),
            Err(error) => Err(error),
        }
    }

    fn get_with_bounds(
        &mut self,
        key: &IntegralKey,
        bounds: Option<&DescendantBounds>,
    ) -> Result<Option<Q>, String> {
        self.report.read_keys += 1;
        encode_key_into(key, &mut self.key);
        self.track_buffer();
        if bounds.is_some_and(|bounds| bounds.excludes(key, &self.key)) {
            return Ok(None);
        }
        let value: Option<Vec<u8>> = self
            .connection
            .prepare_cached("SELECT v FROM faber_integrals WHERE k=?1")
            .map_err(sql_error)?
            .query_row([self.key.as_slice()], |row| row.get(0))
            .optional()
            .map_err(sql_error)?;
        let result = match value {
            Some(bytes) => decode_rational(&bytes).map(Some),
            // The same snapshot proved all positive-C alias namespaces empty.
            // Even an in-range hole needs no unranked/legacy point queries.
            None if bounds.is_some() && DescendantBounds::eligible(key, &self.key) => Ok(None),
            None => self
                .legacy
                .get_after_current_miss(&mut self.connection, key),
        };
        if result.as_ref().is_ok_and(|value| value.is_some()) {
            self.report.read_hits += 1;
        }
        result
    }

    #[cfg(test)]
    fn legacy_get(&self, key: &[u8]) -> Result<Option<Q>, String> {
        self.legacy.legacy_get(&self.connection, key)
    }

    fn covered(&self, targets: &[usize], durable: bool) -> bool {
        self.sources
            .iter()
            .zip(targets)
            .all(|(source, &target)| (if durable { source.durable } else { source.read }) >= target)
    }

    /// Reply only after every target is durable. A barrier already covered by
    /// an earlier transaction does not force an empty or unrelated transaction.
    fn complete_barriers(&mut self) -> Result<bool, String> {
        if self.barriers.iter().any(|barrier| {
            !self.covered(&barrier.targets, true) && self.covered(&barrier.targets, false)
        }) {
            self.commit()?;
        }
        let mut index = 0;
        while index < self.barriers.len() {
            if self.covered(&self.barriers[index].targets, true) {
                let barrier = self.barriers.remove(index);
                self.record_memory(process_budget().force_check(0));
                let _ = barrier.reply.send(Ok(self.report.clone()));
                if barrier.close {
                    return Ok(true);
                }
            } else {
                index += 1;
            }
        }
        Ok(false)
    }

    fn commit_due(&self) -> bool {
        !self.buffer.rows.is_empty()
            && (self.buffer.bytes.len() >= TRANSACTION_BYTES
                || self.buffer.rows.len() >= self.transaction_rows
                || self
                    .buffer
                    .first
                    .is_some_and(|first| first.elapsed() >= self.transaction_age))
    }

    fn commit(&mut self) -> Result<(), String> {
        if self.buffer.rows.is_empty() {
            return Ok(());
        }
        let started = Instant::now();
        if self.sort_writes {
            let sort = Instant::now();
            let bytes = &self.buffer.bytes;
            self.buffer.rows.sort_unstable_by(|left, right| {
                bytes[left.key.clone()]
                    .cmp(&bytes[right.key.clone()])
                    .then_with(|| left.key.start.cmp(&right.key.start))
            });
            self.report.sort_ns += sort.elapsed().as_nanos();
        }
        let transaction = self.connection.transaction().map_err(sql_error)?;
        let insert_started = Instant::now();
        {
            let mut insert = transaction
                .prepare_cached("INSERT OR REPLACE INTO faber_integrals (k, v) VALUES (?1, ?2)")
                .map_err(sql_error)?;
            for row in &self.buffer.rows {
                insert
                    .execute(params![
                        &self.buffer.bytes[row.key.clone()],
                        &self.buffer.bytes[row.value.clone()]
                    ])
                    .map_err(sql_error)?;
            }
        }
        self.report.insert_ns += insert_started.elapsed().as_nanos();
        let commit_started = Instant::now();
        transaction.commit().map_err(sql_error)?;
        self.report.transaction_commit_ns += commit_started.elapsed().as_nanos();
        self.report.commit_ns += started.elapsed().as_nanos();
        self.report.commits += 1;
        self.report.rows += self.buffer.rows.len();
        self.report.bytes += self.buffer.bytes.len();
        for index in self.buffer.sources.drain(..) {
            let source = &mut self.sources[index];
            source.durable = source.read;
            source.source.durable(source.durable);
            source.pending = false;
        }
        self.buffer.bytes.clear();
        self.buffer.rows.clear();
        self.buffer.first = None;
        Ok(())
    }

    fn track_buffer(&mut self) {
        let row_bytes = self
            .buffer
            .rows
            .len()
            .saturating_mul(size_of::<EncodedRow>());
        let row_capacity = self
            .buffer
            .rows
            .capacity()
            .saturating_mul(size_of::<EncodedRow>());
        let sources_bytes = self.buffer.sources.len().saturating_mul(size_of::<usize>());
        let sources_capacity = self
            .buffer
            .sources
            .capacity()
            .saturating_mul(size_of::<usize>());
        self.report.max_buffer_bytes = self.report.max_buffer_bytes.max(
            self.buffer
                .bytes
                .len()
                .saturating_add(row_bytes)
                .saturating_add(sources_bytes)
                .saturating_add(self.key.len()),
        );
        self.report.max_buffer_capacity_bytes = self.report.max_buffer_capacity_bytes.max(
            self.buffer
                .bytes
                .capacity()
                .saturating_add(row_capacity)
                .saturating_add(sources_capacity)
                .saturating_add(self.key.capacity()),
        );
    }

    fn pull_one(&mut self) -> Result<(bool, bool), String> {
        let count = self.sources.len();
        let mut busy = false;
        // One work-IO turn gets one real completion window, rather than only
        // one chance which an empty historical source can consume. Visit each
        // currently discovered source at most once and retain the round-robin
        // cursor after success, so control requests and other live sources get
        // their next bounded turn.
        for _ in 0..count {
            let index = self.next_source % count;
            self.next_source = (index + 1) % count;
            let (progressed, source_busy) = self.pull_source(index)?;
            busy |= source_busy;
            if progressed {
                return Ok((true, busy));
            }
        }
        Ok((false, busy))
    }

    fn pull_source(&mut self, index: usize) -> Result<(bool, bool), String> {
        let cursor = self.sources[index].read;
        let announced = self.sources[index].source.ready();
        self.report.pulls += 1;
        let started = Instant::now();
        let row_start = self.buffer.rows.len();
        let buffer = &mut self.buffer;
        // Even ready == cursor may have bounded planning identity preparation
        // to perform. Empty and busy sources still receive one visit per pass.
        let next = self.sources[index].source.pull_encoded(
            cursor,
            COPY_WINDOW_BYTES,
            &mut |key, value| {
                buffer.first.get_or_insert(started);
                let key_start = buffer.bytes.len();
                buffer.bytes.extend_from_slice(key);
                let value_start = buffer.bytes.len();
                encode_rational_into(value, &mut buffer.bytes);
                buffer.rows.push(EncodedRow {
                    key: key_start..value_start,
                    value: value_start..buffer.bytes.len(),
                });
            },
        );
        self.report.encode_ns += started.elapsed().as_nanos();
        let rows = self.buffer.rows.len() - row_start;
        let Some(next) = next else {
            if rows != 0 {
                return Err("busy cache source encoded a partial completion window".to_owned());
            }
            self.report.busy_pulls += 1;
            return Ok((false, true));
        };
        if next.checked_sub(cursor) != Some(rows) {
            return Err("cache source returned an invalid completion cursor".to_owned());
        }
        if rows == 0 {
            if announced > cursor {
                return Err(
                    "cache source could not read its published completion target".to_owned(),
                );
            }
            self.report.empty_pulls += 1;
            return Ok((false, false));
        }
        // The bounded callback only appends encoded bytes. All mathematical
        // source locks have been released before sorting or SQLite access.
        self.sources[index].read = next;
        if !self.sources[index].pending {
            self.sources[index].pending = true;
            self.buffer.sources.push(index);
        }
        self.track_buffer();
        Ok((true, false))
    }

    fn run(&mut self, requests: &mpsc::Receiver<Request>) -> Result<(), String> {
        loop {
            if self.control.stop.load(Ordering::Acquire) {
                return Err("cache storage stopped before pending requests completed".to_owned());
            }
            self.discover();
            // At most one ordinary control request, then one bounded source
            // window. Neither a stream of reads nor pull work starves the other.
            let handled_request = match requests.try_recv() {
                Ok(request) => {
                    self.request(request)?;
                    true
                }
                Err(mpsc::TryRecvError::Empty) => false,
                Err(mpsc::TryRecvError::Disconnected) => return Ok(()),
            };
            if self.complete_barriers()? {
                return Ok(());
            }
            if self.commit_due() {
                self.commit()?;
            }
            self.memory_step()?;
            let workspace_progress = self.workspaces.step(&mut self.connection)?;
            let work = self.workspaces.report();
            self.report.work_written_blocks = work.written_blocks;
            self.report.work_written_bytes = work.written_bytes;
            self.report.work_commits = work.commits;
            self.report.work_loaded_blocks = work.loaded_blocks;
            self.report.work_loaded_bytes = work.loaded_bytes;
            let (progressed, busy) = self.pull_one()?;
            if progressed || workspace_progress || handled_request {
                continue;
            }
            // pull_one already completed a bounded pass over all sources;
            // another sources.len() empty turns would repeat that work N times.
            let pause = if busy { BUSY_PAUSE } else { IDLE_PAUSE };
            // Requests and pressure readback/reclaim publish an unpark token.
            // A token arriving before this call also skips the timer, avoiding
            // a poll interval for each sequentially requested block.
            thread::park_timeout(pause);
        }
    }

    fn fail_waiters(&mut self, error: &str) {
        self.workspaces.fail(error);
        for barrier in self.barriers.drain(..) {
            let _ = barrier.reply.send(Err(error.to_owned()));
        }
    }

    fn memory_step(&mut self) -> Result<(), String> {
        if self.last_memory_poll.elapsed() < Duration::from_millis(100) {
            return Ok(());
        }
        self.last_memory_poll = Instant::now();
        let budget = process_budget();
        let decision = budget.refresh_if_due();
        self.record_memory(decision);
        if !matches!(decision.pressure, Pressure::Warning | Pressure::Hard) {
            return Ok(());
        }
        // A successful transaction precedes every reclamation. No source or
        // mathematical lock is held while committing or obtaining the gate.
        self.commit()?;
        let registry = Arc::clone(&self.registry);
        self.reclaim_sources(&registry);
        // GNU malloc may retain pages of already-dropped keys or term blocks.
        // Run only here, off the computation threads and at most every 100 ms
        // under pressure. All source guards from reclaim_sources are gone.
        #[cfg(all(target_os = "linux", target_env = "gnu"))]
        if matches!(budget.check(0).pressure, Pressure::Warning | Pressure::Hard) {
            let started = Instant::now();
            super::budget::release_allocator_pages();
            self.report.allocator_trim_calls = self.report.allocator_trim_calls.saturating_add(1);
            self.report.allocator_trim_ns = self
                .report
                .allocator_trim_ns
                .saturating_add(started.elapsed().as_nanos());
            self.record_memory(budget.force_check(0));
        }
        Ok(())
    }

    fn reclaim_sources(&mut self, registry: &Arc<Registry>) {
        let budget = process_budget();
        // The original global quiescent path remains available. Otherwise
        // only sources with their own borrower/durability protocol opt in.
        let reclaim = registry.try_reclaim();
        let mut order: Vec<_> = (0..self.sources.len()).collect();
        order.sort_unstable_by_key(|&index| (self.sources[index].source.heat(), index));
        for index in order {
            let source = &self.sources[index];
            if source.durable < source.source.ready() {
                continue;
            }
            let reclaimed = if reclaim.is_some() {
                source.source.reclaim()
            } else {
                source.source.reclaim_during_query(registry)
            };
            if reclaimed != 0 {
                self.report.reclaimed_bytes = self.report.reclaimed_bytes.saturating_add(reclaimed);
                registry.set_history(true);
                let decision = budget.force_check(0);
                if !matches!(decision.pressure, Pressure::Warning | Pressure::Hard) {
                    break;
                }
            }
        }
    }

    fn record_memory(&mut self, decision: super::budget::BudgetDecision) {
        let work = self.workspaces.report();
        self.report.work_save_ns = work.work_save_ns;
        self.report.work_commits = work.commits;
        self.report.work_save_sql_ns = work.save_sql_ns;
        self.report.work_save_large_blocks = work.save_large_blocks;
        self.report.work_save_large_bytes = work.save_large_bytes;
        self.report.work_save_large_ns = work.save_large_ns;
        self.report.work_cache_write_pages = work.cache_write_pages;
        self.report.work_cache_spill_pages = work.cache_spill_pages;
        self.report.work_load_ns = work.work_load_ns;
        self.report.work_read_sql_ns = work.work_read_sql_ns;
        self.report.work_reclaim_ns = work.work_reclaim_ns;
        self.report.work_max_loaded_block_resident_bytes = work.max_loaded_block_resident_bytes;
        self.report.work_max_block_encoded_bytes = work.max_block_encoded_bytes;
        self.report.last_available_bytes = decision.memory.available_physical_bytes;
        self.report.last_memory_limit_bytes = decision.limits.map(|limits| limits.hard_bytes);
        self.report.last_warning_bytes = decision.limits.map(|limits| limits.warning_bytes);
        self.report.last_cgroup_limit_bytes = decision.memory.cgroup_limit_bytes;
        self.report.last_reclaim_target_bytes = decision.reclaim_bytes;
        self.report.last_progress_at_risk = decision.progress_at_risk;
        self.report.last_active_reserved_bytes = decision.active_reserved_bytes;
        self.report.last_settled_reserved_bytes = decision.settled_reserved_bytes;
        self.report.last_sample_generation = decision.sample_generation;
        self.report.peak_sampled_private_bytes = self
            .report
            .peak_sampled_private_bytes
            .max(decision.known_private_bytes.unwrap_or(0));
        if let Some(peak) = decision.memory.process_peak_commit_bytes {
            self.report.peak_native_commit_bytes =
                Some(self.report.peak_native_commit_bytes.unwrap_or(0).max(peak));
        }
    }
}

fn diagnostic_usize(name: &str, default: usize) -> Result<usize, String> {
    match std::env::var(name) {
        Ok(value) => value
            .parse()
            .map_err(|_| format!("{name} must be an unsigned integer")),
        Err(std::env::VarError::NotPresent) => Ok(default),
        Err(error) => Err(format!("cannot read {name}: {error}")),
    }
}

fn sql_error(error: rusqlite::Error) -> String {
    format!("cache SQLite error: {error}")
}

#[cfg(test)]
fn parse_legacy(value: &str) -> Result<Q, String> {
    let (numerator, denominator) = value
        .split_once('/')
        .ok_or_else(|| "legacy cache rational is missing its denominator".to_owned())?;
    if numerator.contains('/') || denominator.contains('/') {
        return Err("legacy cache rational has extra separators".to_owned());
    }
    let numerator = numerator
        .parse::<Q>()
        .map_err(|_| "invalid legacy cache numerator".to_owned())?;
    let denominator = denominator
        .parse::<Q>()
        .map_err(|_| "invalid legacy cache denominator".to_owned())?;
    if denominator == Q::ZERO {
        return Err("zero denominator in legacy cache rational".to_owned());
    }
    Ok(numerator / denominator)
}

fn panic_error(payload: Box<dyn Any + Send>) -> String {
    let message = payload
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| payload.downcast_ref::<&str>().copied())
        .unwrap_or("unknown panic");
    format!("cache storage worker panicked: {message}")
}

fn worker_main(
    path: &str,
    registry: Arc<Registry>,
    requests: mpsc::Receiver<Request>,
    ready: mpsc::Sender<Result<(), String>>,
    control: Arc<Control>,
) {
    let _ = control.worker.set(thread::current());
    let opened = catch_unwind(AssertUnwindSafe(|| {
        Worker::open(path, Arc::clone(&registry), Arc::clone(&control))
    }));
    let mut worker = match opened {
        Ok(Ok(worker)) => worker,
        result => {
            let error = match result {
                Ok(Err(error)) => error,
                Err(payload) => panic_error(payload),
                Ok(Ok(_)) => unreachable!(),
            };
            control.fail(error.clone());
            registry.fail(error.clone());
            let _ = ready.send(Err(error));
            return;
        }
    };
    if ready.send(Ok(())).is_err() {
        return;
    }
    let result = catch_unwind(AssertUnwindSafe(|| worker.run(&requests)));
    let error = match result {
        Ok(Ok(())) => return,
        Ok(Err(error)) => error,
        Err(payload) => panic_error(payload),
    };
    control.fail(error.clone());
    registry.fail(error.clone());
    worker.fail_waiters(&error);
    while let Ok(request) = requests.try_recv() {
        request.fail(&error);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache_source::{CompletedIntegral, Pulled};
    use std::path::PathBuf;
    use std::sync::RwLock;
    use std::sync::atomic::AtomicUsize;

    static NEXT_TEST_FILE: AtomicUsize = AtomicUsize::new(0);

    #[test]
    fn explicit_request_keeps_a_wake_token_when_sent_before_worker_parks() {
        let control = Arc::new(Control::default());
        let (requests, receive) = mpsc::channel();
        let handle = StorageHandle {
            requests,
            control: Arc::clone(&control),
        };
        let allowed = Arc::new(AtomicBool::new(false));
        let (bound, binding) = mpsc::channel();
        let (woke, awake) = mpsc::channel();
        let worker_control = Arc::clone(&control);
        let worker_allowed = Arc::clone(&allowed);
        let worker = thread::spawn(move || {
            worker_control.worker.set(thread::current()).unwrap();
            bound.send(()).unwrap();
            // Keep the explicit unpark token intact until the actual park.
            while !worker_allowed.load(Ordering::Acquire) {
                thread::yield_now();
            }
            thread::park_timeout(Duration::from_secs(5));
            assert!(matches!(receive.try_recv(), Ok(Request::Flush(_))));
            woke.send(()).unwrap();
        });
        binding.recv_timeout(Duration::from_secs(5)).unwrap();
        let (reply, _response) = mpsc::channel();
        handle.send(Request::Flush(reply)).unwrap();
        allowed.store(true, Ordering::Release);
        let notified = awake.recv_timeout(Duration::from_secs(1));
        // Cleanup remains bounded even when the wakeup assertion fails.
        worker.thread().unpark();
        worker.join().unwrap();
        assert!(
            notified.is_ok(),
            "explicit request did not wake the parked backend"
        );
    }

    struct Database(PathBuf);

    impl Database {
        fn new() -> Self {
            Self(std::env::temp_dir().join(format!(
                "sc-faber-storage-{}-{}.sqlite",
                std::process::id(),
                NEXT_TEST_FILE.fetch_add(1, Ordering::Relaxed),
            )))
        }

        fn path(&self) -> &str {
            self.0.to_str().unwrap()
        }
    }

    impl Drop for Database {
        fn drop(&mut self) {
            for suffix in ["", "-wal", "-shm"] {
                let mut path = self.0.as_os_str().to_os_string();
                path.push(suffix);
                let _ = std::fs::remove_file(PathBuf::from(path));
            }
        }
    }

    #[derive(Default)]
    struct Values(RwLock<Vec<CompletedIntegral>>);

    impl Values {
        fn append(&self, key: IntegralKey, value: Q) {
            self.0
                .write()
                .unwrap()
                .push(CompletedIntegral { key, value });
        }
    }

    impl ResultSource for Values {
        fn ready(&self) -> usize {
            self.0.read().unwrap().len()
        }

        fn pull(&self, cursor: usize, _byte_limit: usize) -> Option<Pulled> {
            let values = self.0.try_read().ok()?;
            let rows: Vec<_> = values.get(cursor).cloned().into_iter().collect();
            Some(Pulled {
                next: cursor + rows.len(),
                rows,
            })
        }
    }

    fn key(degrees: &[u32]) -> IntegralKey {
        IntegralKey::Descendant {
            g: 2,
            d: degrees.to_vec(),
            ch: Vec::new(),
        }
    }

    struct LookupKeys(Vec<IntegralKey>, bool);

    impl LookupBatch for LookupKeys {
        fn len(&self) -> usize {
            self.0.len()
        }
        fn key(&self, index: usize) -> IntegralKey {
            self.0[index].clone()
        }
        fn descendants_only(&self) -> bool {
            self.1
        }
    }

    fn c_key(weight: u32) -> IntegralKey {
        IntegralKey::Descendant {
            g: 2,
            d: vec![],
            ch: vec![(1, weight)],
        }
    }

    fn put_typed(connection: &Connection, key: &IntegralKey, value: &Q, old: bool) {
        let mut encoded = Vec::new();
        if old {
            super::super::codec::encode_unranked_key_into(key, &mut encoded);
        } else {
            encode_key_into(key, &mut encoded);
        }
        let mut bytes = Vec::new();
        encode_rational_into(value, &mut bytes);
        connection
            .execute(
                "INSERT OR REPLACE INTO faber_integrals VALUES(?1,?2)",
                params![encoded, bytes],
            )
            .unwrap();
    }

    fn empty_worker() -> Worker {
        Worker::open(
            ":memory:",
            Arc::new(Registry::new()),
            Arc::new(Control::default()),
        )
        .unwrap()
    }

    fn seeded_worker(seed: impl FnOnce(&Connection)) -> (Database, Worker) {
        let database = Database::new();
        {
            let writer = Worker::open(
                database.path(),
                Arc::new(Registry::new()),
                Arc::new(Control::default()),
            )
            .unwrap();
            seed(&writer.connection);
        }
        let worker = Worker::open(
            database.path(),
            Arc::new(Registry::new()),
            Arc::new(Control::default()),
        )
        .unwrap();
        (database, worker)
    }

    fn seed_formula_history(connection: &Connection) {
        let mut one_point = Q::from(1);
        for g in 1..=18 {
            one_point /= Q::from(24 * g);
            put_typed(
                connection,
                &IntegralKey::Descendant {
                    g,
                    d: vec![3 * g - 2],
                    ch: vec![],
                },
                &one_point,
                false,
            );
        }
        for d in [vec![0], vec![0, 0, 2]] {
            put_typed(
                connection,
                &IntegralKey::Descendant {
                    g: 1,
                    d,
                    ch: vec![(1, 1)],
                },
                &Q::from_signeds(1, 2),
                false,
            );
        }
        for weight in 4..=36 {
            put_typed(
                connection,
                &IntegralKey::Descendant {
                    g: 2,
                    d: vec![0; weight as usize - 3],
                    ch: vec![(1, weight)],
                },
                &Q::ZERO,
                false,
            );
        }
    }

    #[test]
    fn fifty_three_immediate_formula_rows_keep_generic_history_without_dag_recovery() {
        let (_database, mut worker) = seeded_worker(seed_formula_history);
        let count: i64 = worker
            .connection
            .query_row("SELECT COUNT(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 53);
        assert!(worker.registry.has_history());
        assert!(!worker.registry.has_descendant_history());
        assert!(!worker.descendant_history());
        assert_eq!(
            worker.get(&key(&[4])).unwrap(),
            Some(Q::from_signeds(1, 1152))
        );
        assert_eq!(
            worker
                .get(&IntegralKey::Descendant {
                    g: 1,
                    d: vec![0],
                    ch: vec![(1, 1)]
                })
                .unwrap(),
            Some(Q::from_signeds(1, 2))
        );
        let root = IntegralKey::Integral {
            g: 2,
            n: 0,
            kappa: vec![],
            psi: vec![],
            lambda: vec![(1, 3)],
            ch: vec![],
        };
        put_typed(&worker.connection, &root, &Q::from_signeds(1, 2880), false);
        assert_eq!(worker.get(&root).unwrap(), Some(Q::from_signeds(1, 2880)));
        assert!(!worker.descendant_history());
        drop(worker);
    }

    #[test]
    fn canonical_dag_history_stays_enabled_for_both_zero_and_nonzero_values() {
        for value in [Q::ZERO, Q::from_signeds(7, 19)] {
            for key in [
                key(&[2, 3]),
                IntegralKey::Descendant {
                    g: 3,
                    d: vec![2, 2, 2],
                    ch: vec![(1, 3)],
                },
            ] {
                let (_database, worker) =
                    seeded_worker(|connection| put_typed(connection, &key, &value, false));
                assert!(worker.registry.has_history());
                assert!(worker.registry.has_descendant_history());
                drop(worker);
            }
        }
    }

    #[test]
    fn unknown_old_oversized_or_numerous_history_never_disables_dag_recovery() {
        for mode in 0..6 {
            let (_database, worker) = seeded_worker(|connection| match mode {
                0 => put_typed(connection, &key(&[4]), &Q::ZERO, true),
                1 => {
                    connection
                        .execute("INSERT INTO faber_integrals VALUES(x'090807',x'00')", [])
                        .unwrap();
                }
                2 => {
                    connection
                        .execute("INSERT INTO st VALUES(x'00','0')", [])
                        .unwrap();
                }
                3 => {
                    connection.execute_batch("CREATE TABLE faber_values(k BLOB PRIMARY KEY,v BLOB NOT NULL); INSERT INTO faber_values VALUES(x'00',x'00')").unwrap();
                }
                4 => {
                    for weight in 4..=260 {
                        put_typed(
                            connection,
                            &IntegralKey::Descendant {
                                g: 2,
                                d: vec![0; weight as usize - 3],
                                ch: vec![(1, weight)],
                            },
                            &Q::ZERO,
                            false,
                        );
                    }
                }
                5 => put_typed(
                    connection,
                    &IntegralKey::Descendant {
                        g: 1,
                        d: vec![0; HISTORY_KEY_BYTES + 1],
                        ch: vec![(1, (HISTORY_KEY_BYTES + 1) as u32)],
                    },
                    &Q::ZERO,
                    false,
                ),
                _ => unreachable!(),
            });
            assert!(worker.registry.has_descendant_history(), "mode={mode}");
            drop(worker);
        }
    }

    #[test]
    fn external_history_commits_are_seen_before_open_and_on_query_entry() {
        let (_database, mut worker) = seeded_worker(seed_formula_history);
        let path: String = worker
            .connection
            .query_row("PRAGMA database_list", [], |row| row.get(2))
            .unwrap();
        assert!(!worker.descendant_history());
        let external = Connection::open(&path).unwrap();
        put_typed(&external, &key(&[2, 3]), &Q::ZERO, false);
        let (reply, receive) = mpsc::channel();
        worker.request(Request::DescendantHistory(reply)).unwrap();
        assert!(receive.recv().unwrap().unwrap());
        assert!(worker.registry.has_descendant_history());
        drop(external);
        drop(worker);
        // The same committed eligible zero is also recognized at a later open.
        let reopened = Worker::open(
            &path,
            Arc::new(Registry::new()),
            Arc::new(Control::default()),
        )
        .unwrap();
        assert!(reopened.registry.has_descendant_history());
        drop(reopened);
    }

    #[test]
    fn own_commits_keep_the_resident_mode_and_missing_probe_tokens_are_conservative() {
        let (_database, mut worker) = seeded_worker(seed_formula_history);
        let source = Arc::new(Values::default());
        source.append(key(&[2, 3]), Q::ZERO);
        worker.registry.register(source);
        worker.discover();
        assert!(worker.pull_one().unwrap().0);
        worker.commit().unwrap();
        assert!(!worker.descendant_history());
        assert_eq!(worker.get(&key(&[2, 3])).unwrap(), Some(Q::ZERO));
        worker.descendant_history_version = None;
        assert!(worker.descendant_history());
        drop(worker);
        let (_database, worker) = seeded_worker(|connection| {
            connection.execute_batch("DROP TABLE faber_integrals; CREATE TABLE faber_integrals(unrelated INTEGER); INSERT INTO faber_integrals VALUES(1)").unwrap();
        });
        assert!(worker.registry.has_descendant_history());
        drop(worker);
    }

    #[test]
    fn backend_can_reclaim_a_durable_cold_source_inside_an_active_registry() {
        let registry = Arc::new(Registry::new());
        let source = Arc::new(crate::cache_source::manual::ManualSource::default());
        source.put("test-active-reclaim", &[1], &Q::from(7));
        registry.register(source.clone());
        let mut worker =
            Worker::open(":memory:", registry.clone(), Arc::new(Control::default())).unwrap();
        worker.discover();
        assert!(worker.pull_one().unwrap().0);
        worker.commit().unwrap();
        assert!(!registry.has_history());
        crate::cache_source::with_registry(registry.clone(), || {
            assert!(registry.try_reclaim().is_none());
            worker.reclaim_sources(&registry);
            assert!(registry.has_history());
            assert!(source.get("test-active-reclaim", &[1]).is_none());
            assert_eq!(source.floor(), 1);
            assert_eq!(source.ready(), 1);
            assert!(worker.report.reclaimed_bytes > 0);
        });
    }

    #[test]
    fn ranked_batch_bounds_preserve_holes_zero_values_and_namespace_priority() {
        let mut worker = empty_worker();
        put_typed(&worker.connection, &c_key(3), &Q::from(1), false);
        put_typed(&worker.connection, &c_key(9), &Q::ZERO, false);
        let bounds = DescendantBounds::read(&worker.connection).unwrap();
        for (weight, excluded) in [(1, true), (3, false), (6, false), (9, false), (12, true)] {
            let mut encoded = Vec::new();
            encode_key_into(&c_key(weight), &mut encoded);
            assert_eq!(bounds.excludes(&c_key(weight), &encoded), excluded);
        }
        assert_eq!(
            worker
                .get_batch(&LookupKeys(
                    vec![c_key(1), c_key(6), c_key(9), c_key(12)],
                    true
                ))
                .unwrap(),
            vec![None, None, Some(Q::ZERO), None]
        );
        let pure = key(&[4]);
        put_typed(
            &worker.connection,
            &IntegralKey::Legacy {
                kind: "st".into(),
                nums: vec![4],
            },
            &Q::from_signeds(7, 19),
            false,
        );
        put_typed(
            &worker.connection,
            &IntegralKey::Integral {
                g: 2,
                n: 0,
                kappa: vec![],
                psi: vec![],
                lambda: vec![],
                ch: vec![(1, 12)],
            },
            &Q::from(99),
            false,
        );
        assert_eq!(
            worker
                .get_batch(&LookupKeys(vec![pure, c_key(12)], true))
                .unwrap(),
            vec![Some(Q::from_signeds(7, 19)), None]
        );
        put_typed(&worker.connection, &c_key(12), &Q::from(5), true);
        assert_eq!(
            worker
                .get_batch(&LookupKeys(vec![c_key(12)], true))
                .unwrap(),
            vec![Some(Q::from(5))]
        );
        // Current exact values outrank even malformed old typed aliases.
        let mut old = Vec::new();
        super::super::codec::encode_unranked_key_into(&c_key(3), &mut old);
        worker
            .connection
            .execute("INSERT INTO faber_integrals VALUES(?1,x'00')", [old])
            .unwrap();
        assert_eq!(
            worker.get_batch(&LookupKeys(vec![c_key(3)], true)).unwrap(),
            vec![Some(Q::from(1))]
        );
        let mut missing = Vec::new();
        encode_key_into(&c_key(6), &mut missing);
        worker
            .connection
            .execute("INSERT INTO faber_integrals VALUES(?1,x'00')", [missing])
            .unwrap();
        assert!(worker.get_batch(&LookupKeys(vec![c_key(6)], true)).is_err());
        assert!(worker.connection.is_autocommit());
    }

    #[test]
    fn stable_old_typed_ranges_exclude_only_current_unstable_auxiliaries() {
        for (g, n, stable) in [
            (0, 0, false),
            (0, 1, false),
            (0, 2, false),
            (1, 0, false),
            (0, 3, true),
            (1, 1, true),
            (2, 0, true),
            (128, 0, true),
            (0, 128, true),
        ] {
            let worker = empty_worker();
            let key = IntegralKey::Descendant {
                g,
                d: vec![0; n],
                ch: vec![(1, 9)],
            };
            put_typed(&worker.connection, &key, &Q::from(1), true);
            assert_eq!(
                worker
                    .legacy
                    .has_positive_descendant_aliases(&worker.connection)
                    .unwrap(),
                stable
            );
        }
    }

    #[test]
    fn legacy_grr_namespaces_disable_bounds_and_keep_exact_error_priority() {
        for prefix in [b"psi-grr|".as_slice(), b"psi-normalized-grr|".as_slice()] {
            let mut worker = empty_worker();
            worker
                .connection
                .execute_batch(
                    "CREATE TABLE faber_values(k BLOB PRIMARY KEY,v BLOB NOT NULL) WITHOUT ROWID",
                )
                .unwrap();
            worker.legacy = LegacyReader::new(&worker.connection).unwrap();
            assert!(
                !worker
                    .legacy
                    .has_positive_descendant_aliases(&worker.connection)
                    .unwrap()
            );
            let mut raw = prefix.to_vec();
            if prefix == b"psi-grr|" {
                raw.extend_from_slice(&9u32.to_be_bytes());
                raw.extend_from_slice(&2u32.to_be_bytes());
            }
            raw.extend_from_slice(&[2, 0, 1, 1, 9]);
            worker
                .connection
                .execute(
                    "INSERT INTO faber_values VALUES(?1,?2)",
                    params![&raw, [0u8, 1, 0, 0, 0, 1, 2].as_slice()],
                )
                .unwrap();
            put_typed(&worker.connection, &c_key(3), &Q::from(1), false);
            assert_eq!(
                worker.get_batch(&LookupKeys(vec![c_key(9)], true)).unwrap(),
                vec![Some(Q::from_signeds(1, 2))]
            );
            worker
                .connection
                .execute("UPDATE faber_values SET v=x'00'", [])
                .unwrap();
            assert_eq!(
                worker.get_batch(&LookupKeys(vec![c_key(3)], true)).unwrap(),
                vec![Some(Q::from(1))]
            );
            assert!(worker.get_batch(&LookupKeys(vec![c_key(9)], true)).is_err());
            assert!(worker.connection.is_autocommit());
        }
    }

    #[test]
    fn bounds_are_rebuilt_after_own_commits_and_external_writers() {
        let database = Database::new();
        let registry = Arc::new(Registry::new());
        let mut worker = Worker::open(
            database.path(),
            registry.clone(),
            Arc::new(Control::default()),
        )
        .unwrap();
        put_typed(&worker.connection, &c_key(3), &Q::from(1), false);
        assert_eq!(
            worker.get_batch(&LookupKeys(vec![c_key(9)], true)).unwrap(),
            vec![None]
        );
        let values = Arc::new(Values::default());
        values.append(c_key(9), Q::from(2));
        registry.register(values);
        worker.discover();
        assert!(worker.pull_one().unwrap().0);
        worker.commit().unwrap();
        assert_eq!(
            worker.get_batch(&LookupKeys(vec![c_key(9)], true)).unwrap(),
            vec![Some(Q::from(2))]
        );
        let other = Connection::open(database.path()).unwrap();
        worker.connection.execute_batch("BEGIN DEFERRED").unwrap();
        let bounds = DescendantBounds::read(&worker.connection).unwrap();
        put_typed(&other, &c_key(12), &Q::from(3), false);
        assert_eq!(
            worker.get_with_bounds(&c_key(12), Some(&bounds)).unwrap(),
            None
        );
        assert_eq!(worker.get(&c_key(12)).unwrap(), None);
        worker.connection.execute_batch("COMMIT").unwrap();
        assert_eq!(
            worker
                .get_batch(&LookupKeys(vec![c_key(12)], true))
                .unwrap(),
            vec![Some(Q::from(3))]
        );
    }

    #[test]
    fn unmarked_mixed_batches_keep_legacy_integral_imports_outside_read_snapshots() {
        let mut worker = empty_worker();
        worker
            .connection
            .execute_batch(
                "CREATE TABLE faber_values(k BLOB PRIMARY KEY,v BLOB NOT NULL) WITHOUT ROWID",
            )
            .unwrap();
        worker.legacy = LegacyReader::new(&worker.connection).unwrap();
        let integral = IntegralKey::Integral {
            g: 2,
            n: 0,
            kappa: vec![],
            psi: vec![],
            lambda: vec![(1, 3)],
            ch: vec![],
        };
        assert_eq!(
            worker
                .get_batch(&LookupKeys(vec![c_key(9), integral], false))
                .unwrap(),
            vec![None, None]
        );
        assert!(worker.connection.is_autocommit());
        // Hint probing may fail on an unrelated malformed historical schema;
        // an exact current hit retains priority and remains readable.
        worker
            .connection
            .execute_batch("DROP TABLE faber_values; CREATE TABLE faber_values(unrelated INTEGER)")
            .unwrap();
        put_typed(&worker.connection, &c_key(9), &Q::from(7), false);
        assert_eq!(
            worker.get_batch(&LookupKeys(vec![c_key(9)], true)).unwrap(),
            vec![Some(Q::from(7))]
        );
    }

    #[derive(Default)]
    struct PlanningSource {
        visits: AtomicUsize,
    }

    impl ResultSource for PlanningSource {
        fn ready(&self) -> usize {
            0
        }

        fn pull(&self, cursor: usize, byte_limit: usize) -> Option<Pulled> {
            assert_eq!(cursor, 0);
            assert_eq!(byte_limit, COPY_WINDOW_BYTES);
            // A real ParallelSource advances its identity directory here,
            // including before the first mathematical layer is published.
            self.visits.fetch_add(1, Ordering::Relaxed);
            Some(Pulled {
                next: cursor,
                rows: Vec::new(),
            })
        }
    }

    #[test]
    fn one_pull_turn_reaches_real_rows_past_147_empty_sources() {
        let registry = Arc::new(Registry::new());
        let empty: Vec<_> = (0..147)
            .map(|_| {
                let source = Arc::new(PlanningSource::default());
                registry.register(source.clone());
                source
            })
            .collect();
        let values = Arc::new(Values::default());
        for degree in 0..4 {
            values.append(key(&[degree]), Q::from(degree + 1));
        }
        registry.register(values);
        let mut worker = Worker::open(":memory:", registry, Arc::new(Control::default())).unwrap();
        worker.discover();
        // run() performs one workspace IO step before each of these calls.
        // Every such turn must reach a real window instead of requiring 148
        // workspace operations to get past the empty-source prefix.
        for turn in 1..=4 {
            assert_eq!(worker.pull_one().unwrap(), (true, false));
            assert_eq!(worker.buffer.rows.len(), turn);
            assert_eq!(worker.sources[147].read, turn);
            assert_eq!(worker.report.pulls, 148 * turn);
            for source in &empty {
                assert_eq!(source.visits.load(Ordering::Relaxed), turn);
            }
        }
        // A full empty pass terminates and visits each planning source once.
        assert_eq!(worker.pull_one().unwrap(), (false, false));
        assert_eq!(worker.report.pulls, 148 * 5);
        assert_eq!(worker.report.empty_pulls, 147 * 5 + 1);
        worker.commit().unwrap();
        assert_eq!(worker.report.rows, 4);
        assert_eq!(worker.report.commits, 1);
        assert_eq!(worker.sources[147].durable, 4);
        assert!(worker.sources[..147].iter().all(|source| source.read == 0));
    }

    #[test]
    fn a_pull_pass_skips_busy_sources_and_still_prepares_ready_zero_sources() {
        let registry = Arc::new(Registry::new());
        registry.register(Arc::new(Busy));
        let planning = Arc::new(PlanningSource::default());
        registry.register(planning.clone());
        let values = Arc::new(Values::default());
        values.append(key(&[4]), Q::from_signeds(1, 1152));
        registry.register(values);
        let mut worker = Worker::open(":memory:", registry, Arc::new(Control::default())).unwrap();
        worker.discover();
        assert_eq!(worker.pull_one().unwrap(), (true, true));
        assert_eq!(worker.report.pulls, 3);
        assert_eq!(worker.report.busy_pulls, 1);
        assert_eq!(planning.visits.load(Ordering::Relaxed), 1);
        assert_eq!(worker.pull_one().unwrap(), (false, true));
        assert_eq!(worker.report.pulls, 6);
        assert_eq!(worker.report.busy_pulls, 2);
        assert_eq!(planning.visits.load(Ordering::Relaxed), 2);
        assert_eq!(worker.sources[2].read, 1);
    }

    #[test]
    fn one_window_per_turn_preserves_round_robin_and_completion_barriers() {
        let registry = Arc::new(Registry::new());
        for offset in [0, 10] {
            let values = Arc::new(Values::default());
            values.append(key(&[offset]), Q::from(1));
            values.append(key(&[offset + 1]), Q::from(2));
            registry.register(values);
        }
        let mut worker = Worker::open(":memory:", registry, Arc::new(Control::default())).unwrap();
        let (reply, completed) = mpsc::channel();
        worker.capture(reply, false);
        for expected in [[1, 0], [1, 1], [2, 1], [2, 2]] {
            assert_eq!(worker.pull_one().unwrap(), (true, false));
            assert_eq!([worker.sources[0].read, worker.sources[1].read], expected);
            assert!(matches!(
                completed.try_recv(),
                Err(mpsc::TryRecvError::Empty)
            ));
            assert!(!worker.complete_barriers().unwrap());
            if expected != [2, 2] {
                assert!(matches!(
                    completed.try_recv(),
                    Err(mpsc::TryRecvError::Empty)
                ));
            }
        }
        let report = completed.try_recv().unwrap().unwrap();
        assert_eq!(report.rows, 4);
        assert_eq!(report.commits, 1);
        assert!(worker.sources.iter().all(|source| source.durable == 2));
    }

    #[test]
    fn flush_and_close_persist_exact_values_without_empty_transactions() {
        let database = Database::new();
        let registry = Arc::new(Registry::new());
        let values = Arc::new(Values::default());
        values.append(key(&[4]), Q::from_signeds(1, 1152));
        values.append(key(&[0, 5]), Q::ZERO);
        registry.register(values);
        let storage = Storage::open(database.path(), registry.clone()).unwrap();
        assert!(!registry.has_history());
        let first = storage.flush().unwrap();
        assert!(!registry.has_history());
        assert_eq!(first.rows, 2);
        assert!(first.commits > 0);
        assert!(first.bytes > 0);
        assert!(first.max_buffer_capacity_bytes >= first.max_buffer_bytes);
        assert_eq!(
            storage.get(key(&[4])).unwrap(),
            Some(Q::from_signeds(1, 1152))
        );
        assert_eq!(storage.get(key(&[0, 5])).unwrap(), Some(Q::ZERO));
        let second = storage.flush().unwrap();
        assert_eq!(second.rows, first.rows);
        assert_eq!(second.commits, first.commits);
        storage.close().unwrap();
        let next_registry = Arc::new(Registry::new());
        let reopened = Storage::open(database.path(), next_registry.clone()).unwrap();
        assert!(next_registry.has_history());
        assert_eq!(
            reopened.get(key(&[4])).unwrap(),
            Some(Q::from_signeds(1, 1152))
        );
        assert_eq!(reopened.get(key(&[0, 5])).unwrap(), Some(Q::ZERO));
        assert_eq!(reopened.get(key(&[9])).unwrap(), None);
        assert_eq!(reopened.close().unwrap().commits, 0);
    }

    #[test]
    fn flush_discovers_later_values_and_later_registered_sources() {
        let database = Database::new();
        let registry = Arc::new(Registry::new());
        let values = Arc::new(Values::default());
        registry.register(values.clone());
        let storage = Storage::open(database.path(), registry.clone()).unwrap();
        assert_eq!(storage.flush().unwrap().commits, 0);
        values.append(key(&[4]), Q::from_signeds(1, 1152));
        assert_eq!(storage.flush().unwrap().rows, 1);
        values.append(key(&[2, 3]), Q::from_signeds(-7, 19));
        let another = Arc::new(Values::default());
        another.append(key(&[1, 4]), Q::from_signeds(1, 384));
        registry.register(another);
        let report = storage.flush().unwrap();
        assert_eq!(report.rows, 3);
        assert_eq!(report.sources, 2);
        assert_eq!(
            storage.get(key(&[2, 3])).unwrap(),
            Some(Q::from_signeds(-7, 19))
        );
        assert_eq!(storage.close().unwrap().rows, 3);
    }

    #[test]
    fn legacy_text_is_read_exactly_and_kept_separate_from_binary_values() {
        let database = Database::new();
        let expected = Q::from_signeds(-7, 19);
        {
            let connection = Connection::open(database.path()).unwrap();
            connection
                .execute_batch(
                    "CREATE TABLE st (k BLOB PRIMARY KEY, v TEXT NOT NULL) WITHOUT ROWID;",
                )
                .unwrap();
            connection
                .execute(
                    "INSERT INTO st VALUES (?1, ?2)",
                    params![b"old-key".as_slice(), "-7/19"],
                )
                .unwrap();
        }
        let storage = Storage::open(database.path(), Arc::new(Registry::new())).unwrap();
        assert_eq!(
            storage.legacy_get(b"old-key".to_vec()).unwrap(),
            Some(expected)
        );
        assert_eq!(storage.legacy_get(b"missing".to_vec()).unwrap(), None);
        assert_eq!(storage.get(key(&[4])).unwrap(), None);
        storage.close().unwrap();
        assert_eq!(
            parse_legacy("-123456789012345678901234567890/3").unwrap(),
            "-123456789012345678901234567890/3".parse::<Q>().unwrap()
        );
        assert!(parse_legacy("1/0").is_err());
        assert!(parse_legacy("1/2/3").is_err());
    }

    #[test]
    fn cloned_handles_keep_their_original_service_after_another_database_opens() {
        let first_database = Database::new();
        let second_database = Database::new();
        for (database, value) in [(&first_database, "1/1"), (&second_database, "2/1")] {
            let connection = Connection::open(database.path()).unwrap();
            connection
                .execute_batch(
                    "CREATE TABLE st (k BLOB PRIMARY KEY, v TEXT NOT NULL) WITHOUT ROWID;",
                )
                .unwrap();
            connection
                .execute(
                    "INSERT INTO st VALUES (?1, ?2)",
                    params![b"same-key".as_slice(), value],
                )
                .unwrap();
        }
        let first = Storage::open(first_database.path(), Arc::new(Registry::new())).unwrap();
        let first_handle = first.handle();
        let second = Storage::open(second_database.path(), Arc::new(Registry::new())).unwrap();
        let second_handle = second.handle();
        assert_eq!(
            first_handle.legacy_get(b"same-key".to_vec()).unwrap(),
            Some(Q::from(1))
        );
        assert_eq!(
            second_handle.legacy_get(b"same-key".to_vec()).unwrap(),
            Some(Q::from(2))
        );
        first.close().unwrap();
        assert!(first_handle.legacy_get(b"same-key".to_vec()).is_err());
        assert_eq!(
            second_handle.legacy_get(b"same-key".to_vec()).unwrap(),
            Some(Q::from(2))
        );
        second.close().unwrap();
    }

    #[test]
    fn failed_transaction_reports_error_and_wakes_flush_and_close() {
        let database = Database::new();
        {
            let connection = Connection::open(database.path()).unwrap();
            connection.execute_batch(
                "CREATE TABLE faber_integrals (k BLOB PRIMARY KEY, v BLOB NOT NULL) WITHOUT ROWID;
                 CREATE TRIGGER deny_write BEFORE INSERT ON faber_integrals BEGIN
                    SELECT RAISE(ABORT, 'test write failure');
                 END;",
            ).unwrap();
        }
        let registry = Arc::new(Registry::new());
        let values = Arc::new(Values::default());
        values.append(key(&[4]), Q::from_signeds(1, 1152));
        registry.register(values);
        let storage = Storage::open(database.path(), registry).unwrap();
        assert!(storage.flush().is_err());
        assert!(storage.error().is_some());
        assert!(storage.get(key(&[4])).is_err());
        assert!(storage.close().is_err());
    }

    struct Fault(bool);

    impl ResultSource for Fault {
        fn ready(&self) -> usize {
            1
        }

        fn pull(&self, cursor: usize, _byte_limit: usize) -> Option<Pulled> {
            assert!(!self.0, "test source panic");
            Some(Pulled {
                next: cursor,
                rows: Vec::new(),
            })
        }
    }

    #[test]
    fn unreadable_published_targets_and_pull_panics_fail_waiters() {
        for panic in [false, true] {
            let database = Database::new();
            let registry = Arc::new(Registry::new());
            registry.register(Arc::new(Fault(panic)));
            let storage = Storage::open(database.path(), registry).unwrap();
            assert!(storage.flush().is_err());
            assert!(storage.close().is_err());
        }
    }

    #[test]
    fn open_failure_is_returned_after_the_worker_exits() {
        let directory = std::env::temp_dir();
        assert!(Storage::open(directory.to_str().unwrap(), Arc::new(Registry::new())).is_err());
    }

    struct Busy;

    impl ResultSource for Busy {
        fn ready(&self) -> usize {
            1
        }

        fn pull(&self, _cursor: usize, _byte_limit: usize) -> Option<Pulled> {
            None
        }
    }

    #[test]
    fn drop_during_unwind_stops_without_waiting_for_a_busy_completion_target() {
        let database = Database::new();
        let registry = Arc::new(Registry::new());
        registry.register(Arc::new(Busy));
        let storage = Storage::open(database.path(), registry).unwrap();
        let result = catch_unwind(AssertUnwindSafe(move || {
            let _storage = storage;
            panic!("test storage owner unwind");
        }));
        assert!(result.is_err());
    }
}
