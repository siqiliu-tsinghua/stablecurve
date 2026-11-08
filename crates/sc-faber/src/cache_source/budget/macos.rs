//! macOS observations for the shared budget; selected by the parent's cfg.
//!
//! Process usage is Apple's physical-footprint ledger, including compressed
//! accounting; resident size is reported separately. Neither is private commit.
//! The available-growth estimate is free + max(purgeable, pageable external -
//! speculative). Speculative is already in free. File-backed pages can require
//! writeback or later faults: this is reclaimable growth, not immediate free RAM
//! or a hard allocation guarantee. The normal system/progress margins still apply.
//!
//! ABI/semantics are from Apple's XNU sources:
//! - bsd/sys/resource.h and libsyscall/wrappers/libproc/libproc.h
//! - osfmk/mach/{host_info,vm_statistics}.h and osfmk/kern/host.c
//! - libsyscall/mach/mach_init.c (host_page_size uses kernel page size).
//! No swap availability is presented as an enforceable system commit limit.

use super::{MemorySnapshot, PrivateMemoryKind};
use std::ffi::{c_char, c_int, c_void};

const RUSAGE_INFO_V2: c_int = 2;
const HOST_VM_INFO64: c_int = 4;
const KERN_SUCCESS: c_int = 0;

#[repr(C)]
#[derive(Default)]
struct RusageInfoV2 {
    _uuid: [u8; 16],
    _user_time: u64,
    _system_time: u64,
    _package_idle_wakeups: u64,
    _interrupt_wakeups: u64,
    _pageins: u64,
    _wired_size: u64,
    resident_size: u64,
    physical_footprint: u64,
    _process_start_abstime: u64,
    _process_exit_abstime: u64,
    _child_user_time: u64,
    _child_system_time: u64,
    _child_package_idle_wakeups: u64,
    _child_interrupt_wakeups: u64,
    _child_pageins: u64,
    _child_elapsed_abstime: u64,
    _disk_bytes_read: u64,
    _disk_bytes_written: u64,
}

/// Rev1 is sufficient for external/internal counts and is supported by older
/// kernels. Newer revisions append fields; our count advertises only this size.
#[repr(C, align(8))]
#[derive(Default)]
struct VmStatistics64Rev1 {
    free_count: u32,
    _active_count: u32,
    _inactive_count: u32,
    _wire_count: u32,
    _zero_fill_count: u64,
    _reactivations: u64,
    _pageins: u64,
    _pageouts: u64,
    _faults: u64,
    _cow_faults: u64,
    _lookups: u64,
    _hits: u64,
    _purges: u64,
    purgeable_count: u32,
    speculative_count: u32,
    _decompressions: u64,
    _compressions: u64,
    _swapins: u64,
    _swapouts: u64,
    _compressor_page_count: u32,
    _throttled_count: u32,
    external_page_count: u32,
    _internal_page_count: u32,
    _uncompressed_pages_in_compressor: u64,
}

// Mach counts are integer_t words, not bytes (38 words for this 152-byte ABI).
const HOST_VM_INFO64_REV1_COUNT: u32 =
    (std::mem::size_of::<VmStatistics64Rev1>() / std::mem::size_of::<c_int>()) as u32;

#[link(name = "System")]
unsafe extern "C" {
    // XNU builds the libproc wrappers into libsystem_kernel, reexported by
    // libSystem. rusage_info_t is void*, so rusage_info_t* addresses struct
    // storage as void**, not a variable containing a struct pointer.
    fn proc_pid_rusage(pid: c_int, flavor: c_int, buffer: *mut *mut c_void) -> c_int;
    fn sysctlbyname(
        name: *const c_char,
        old_value: *mut c_void,
        old_length: *mut usize,
        new_value: *mut c_void,
        new_length: usize,
    ) -> c_int;
    fn mach_host_self() -> u32;
    fn host_page_size(host: u32, page_size: *mut usize) -> c_int;
    fn host_statistics64(host: u32, flavor: c_int, info: *mut c_int, count: *mut u32) -> c_int;
    fn mach_port_deallocate(task: u32, name: u32) -> c_int;
    static mach_task_self_: u32;
}

struct HostPort(u32);

impl HostPort {
    fn acquire() -> Option<Self> {
        // SAFETY: parameterless Mach call returns a send right for this task.
        let port = unsafe { mach_host_self() };
        if port == 0 || port == u32::MAX {
            None
        } else {
            Some(Self(port))
        }
    }
}

impl Drop for HostPort {
    fn drop(&mut self) {
        // SAFETY: this object owns the send right from mach_host_self. libSystem
        // initializes the current task port; releasing this send right does not
        // destroy the host or release any other caller's reference.
        unsafe { mach_port_deallocate(mach_task_self_, self.0) };
    }
}

fn process_memory() -> Option<RusageInfoV2> {
    let process = c_int::try_from(std::process::id()).ok()?;
    let mut usage = RusageInfoV2::default();
    // SAFETY: the full V2 C layout is provided for flavor V2. The unusual
    // void** parameter is Apple's buffer typedef, not an extra pointer level.
    let result = unsafe {
        proc_pid_rusage(
            process,
            RUSAGE_INFO_V2,
            (&mut usage as *mut RusageInfoV2).cast::<*mut c_void>(),
        )
    };
    (result == 0).then_some(usage)
}

fn physical_memory() -> Option<u64> {
    let mut total = 0u64;
    let mut length = std::mem::size_of::<u64>();
    // SAFETY: the fixed C string is NUL terminated and old_value has exactly
    // the documented uint64_t storage for hw.memsize. This is a read-only call.
    let result = unsafe {
        sysctlbyname(
            c"hw.memsize".as_ptr(),
            (&mut total as *mut u64).cast::<c_void>(),
            &mut length,
            std::ptr::null_mut(),
            0,
        )
    };
    (result == 0 && length == std::mem::size_of::<u64>() && total != 0).then_some(total)
}

fn vm_statistics() -> Option<(VmStatistics64Rev1, u64)> {
    let host = HostPort::acquire()?;
    let mut page_size = 0usize;
    // SAFETY: the live host send right and vm_size_t storage match the Mach ABI.
    // Apple's wrapper returns vm_kernel_page_size, including under Rosetta.
    if unsafe { host_page_size(host.0, &mut page_size) } != KERN_SUCCESS
        || page_size == 0
        || !page_size.is_power_of_two()
    {
        return None;
    }
    let mut statistics = VmStatistics64Rev1::default();
    let mut count = HOST_VM_INFO64_REV1_COUNT;
    // SAFETY: aligned Rev1 storage is advertised in 32-bit integer_t units.
    // The kernel fills only the advertised compatible revision. No fields from
    // later revisions are inspected or assumed to have been returned.
    if unsafe {
        host_statistics64(
            host.0,
            HOST_VM_INFO64,
            (&mut statistics as *mut VmStatistics64Rev1).cast::<c_int>(),
            &mut count,
        )
    } != KERN_SUCCESS
        || count != HOST_VM_INFO64_REV1_COUNT
    {
        return None;
    }
    Some((statistics, u64::try_from(page_size).ok()?))
}

fn growth_estimate(
    free: u64,
    purgeable: u64,
    external: u64,
    speculative: u64,
    page_size: u64,
    physical_total: u64,
) -> Option<u64> {
    if page_size == 0 || speculative > free {
        return None;
    }
    // The external pageable total includes speculative file pages. Remove
    // them before adding any reclaimable estimate to free. Using max rather
    // than sum does not assume purgeable and external categories are disjoint.
    let reclaimable = purgeable.max(external.saturating_sub(speculative));
    let bytes = free.checked_add(reclaimable)?.checked_mul(page_size)?;
    (bytes <= physical_total).then_some(bytes)
}

pub(super) fn capture() -> MemorySnapshot {
    let mut snapshot = MemorySnapshot::unknown();
    let process_ok = if let Some(usage) = process_memory() {
        snapshot.process_private_bytes = Some(usage.physical_footprint);
        snapshot.process_private_kind = Some(PrivateMemoryKind::PhysicalFootprint);
        snapshot.process_resident_bytes = Some(usage.resident_size);
        true
    } else {
        false
    };
    snapshot.total_physical_bytes = physical_memory();
    if let Some((statistics, page_size)) = vm_statistics() {
        snapshot.available_physical_bytes = snapshot.total_physical_bytes.and_then(|total| {
            growth_estimate(
                u64::from(statistics.free_count),
                u64::from(statistics.purgeable_count),
                u64::from(statistics.external_page_count),
                u64::from(statistics.speculative_count),
                page_size,
                total,
            )
        });
    }
    // Only the documented footprint/system observations are covered, not
    // private task quotas or a promise about immediately reclaimable pages.
    snapshot.environment_limits_checked = process_ok
        && snapshot.total_physical_bytes.is_some()
        && snapshot.available_physical_bytes.is_some();
    snapshot
}

#[cfg(all(test, target_os = "macos"))]
mod tests {
    use super::*;

    #[test]
    fn ffi_layouts_match_the_versioned_apple_abi() {
        assert_eq!(std::mem::size_of::<RusageInfoV2>(), 160);
        assert_eq!(std::mem::size_of::<VmStatistics64Rev1>(), 152);
        assert_eq!(std::mem::align_of::<VmStatistics64Rev1>(), 8);
        assert_eq!(HOST_VM_INFO64_REV1_COUNT, 38);
    }

    #[test]
    fn growth_does_not_count_speculation_or_reclaimable_categories_twice() {
        assert_eq!(
            growth_estimate(100, 30, 80, 20, 16_384, 10_000_000),
            Some(160 * 16_384)
        );
        assert_eq!(
            growth_estimate(100, 90, 80, 20, 16_384, 10_000_000),
            Some(190 * 16_384)
        );
        assert_eq!(
            growth_estimate(100, 30, 10, 20, 4096, 10_000_000),
            Some(130 * 4096)
        );
        assert!(growth_estimate(10, 0, 0, 11, 4096, 10_000_000).is_none());
        assert!(growth_estimate(u64::MAX, 1, 0, 0, 4096, u64::MAX).is_none());
        assert!(growth_estimate(100, 0, 0, 0, 16_384, 100).is_none());
    }

    #[test]
    fn native_footprint_and_system_growth_have_explicit_consistent_units() {
        let process = process_memory().expect("proc_pid_rusage V2 unavailable");
        assert!(process.resident_size > 0);
        let (statistics, page_size) = vm_statistics().expect("host_statistics64 Rev1 unavailable");
        assert!(page_size.is_power_of_two());
        let total = physical_memory().expect("hw.memsize unavailable");
        assert!(u64::from(statistics.speculative_count) <= u64::from(statistics.free_count));
        let snapshot = capture();
        assert_eq!(
            snapshot.process_private_kind,
            Some(PrivateMemoryKind::PhysicalFootprint)
        );
        assert!(snapshot.process_private_bytes.is_some());
        assert!(snapshot.process_resident_bytes.is_some());
        assert_eq!(snapshot.total_physical_bytes, Some(total));
        assert!(
            snapshot
                .available_physical_bytes
                .is_some_and(|available| available <= total)
        );
        assert!(snapshot.environment_limits_checked);
        assert_eq!(snapshot.commit_headroom_bytes, None);
        assert!(!snapshot.commit_limit_enforced);
    }
}
