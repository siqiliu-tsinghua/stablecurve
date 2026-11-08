//! Optional, temporary workspace backing for memory-pressure slow paths.
//!
//! A logical block may exceed SQLite's single-BLOB limit. The service streams
//! its bytes through bounded chunk rows without changing mathematical block or
//! equation boundaries. This is temporary working storage, not a durable
//! checkpoint or a persistent mathematical-value namespace.

use super::budget::Reservation;
use std::io::{Read, Write};
use std::sync::Arc;
use std::thread::Thread;

pub(crate) trait WorkPayload: Send + Sync {
    /// Vector capacities and live numeric storage estimate; not process RSS.
    fn resident_bytes(&self) -> usize;
    fn encoded_len(&self) -> Result<u64, String>;
    /// Called only by storage, after all source locks have been released.
    fn encode(&self, output: &mut dyn Write) -> Result<(), String>;
}

pub(crate) struct WorkBlock {
    pub(crate) block: usize,
    pub(crate) payload: Arc<dyn WorkPayload>,
}

#[derive(Clone)]
pub(crate) struct WorkLoad {
    pub(crate) block: usize,
    pub(crate) resident_bytes: usize,
    /// The coordinator owns this whole job group's admission and settlement.
    /// A backend borrows the same credit; it must not reserve it a second time
    /// or settle the group. None is retained for legacy internal/test callers.
    pub(crate) reservation: Option<Reservation>,
}

/// One append-only collection of sealed blocks from an active calculation.
///
/// The service must not hold SQL, budget, or source locks while calling another
/// source. `pull` pins existing data without copying mathematical values.
/// After the whole logical block is successfully stored, call `stored`, then
/// drop WorkBlock. The source may retain its owner until pressure requests
/// eviction; a successful save alone does not require immediate release. Failed
/// storage/readback must call `fail` on affected sources, waking blocked jobs.
pub(crate) trait WorkSource: Send + Sync {
    /// Bind the sole backend before acknowledging workspace registration.
    /// Only pressure preparation/control requests use this wakeup channel.
    fn bind_worker(&self, _worker: Thread) {}
    fn workspace_id(&self) -> u64;
    fn pull(&self) -> Option<WorkBlock>;
    fn next_load(&self) -> Option<WorkLoad>;
    fn stored(&self, block: usize) -> Result<(), String>;
    /// Decode on the backend without source locks; publish only a complete block.
    fn restore(&self, block: usize, input: &mut dyn Read) -> Result<(), String>;
    fn fail(&self, error: String);
    /// Pressure coordinators change metadata; storage alone releases owners.
    fn take_reclaim_request(&self) -> usize {
        0
    }
    /// O(1) signal that candidates/pins or an unsaved suffix still need a drain.
    fn reclaim_pending(&self) -> bool {
        false
    }
    /// Drop stored, unpinned resident blocks outside the source lock. Whole
    /// blocks may overshoot byte_limit; returned bytes are capacity estimates.
    fn reclaim_stored(&self, byte_limit: usize) -> usize;
    /// The service can delete all workspace chunk rows after retirement. Workspaces
    /// retire after their jobs have joined; in-flight backend pins may outlive it.
    fn is_retired(&self) -> bool;
}
