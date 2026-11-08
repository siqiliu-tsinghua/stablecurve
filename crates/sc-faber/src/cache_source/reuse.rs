//! Coarse L1 reuse observations at public-root and query-batch boundaries.
//!
//! These are the observed root requests, not all recursive descendant lookups.
//! One packed atomic keeps requests and hits coherent. Counts are halved only
//! on overflow; no clock, global epoch, or arithmetic-item instrumentation is
//! needed. A wholly reclaimed source starts a fresh observation history.

use std::sync::atomic::{AtomicU64, Ordering};

const COUNT_MAX: u64 = u32::MAX as u64;
const SCALE: u64 = 1_000_000;

#[derive(Default)]
pub(crate) struct ReuseStats(AtomicU64);

impl ReuseStats {
    /// Record one public lookup or one already counted query-root batch.
    /// The low word holds hits; the high word holds their shared denominator.
    pub(crate) fn record(&self, hits: usize, requests: usize) {
        if requests == 0 {
            return;
        }
        let mut hits = u64::try_from(hits.min(requests)).unwrap_or(u64::MAX);
        let mut requests = u64::try_from(requests).unwrap_or(u64::MAX);
        while requests > COUNT_MAX {
            hits = hits.div_ceil(2);
            requests = requests.div_ceil(2);
        }
        let mut previous = self.0.load(Ordering::Relaxed);
        loop {
            let mut total = (previous >> 32) + requests;
            let mut reused = (previous & COUNT_MAX) + hits;
            if total > COUNT_MAX {
                total = total.div_ceil(2);
                reused = reused.div_ceil(2);
            }
            let next = (total << 32) | reused;
            match self
                .0
                .compare_exchange_weak(previous, next, Ordering::Relaxed, Ordering::Relaxed)
            {
                Ok(_) => return,
                Err(current) => previous = current,
            }
        }
    }

    /// Root L1 hit proportion in millionths; no observations rank as zero.
    pub(crate) fn score(&self) -> u64 {
        let packed = self.0.load(Ordering::Relaxed);
        let requests = packed >> 32;
        ((packed & COUNT_MAX) * SCALE)
            .checked_div(requests)
            .unwrap_or(0)
    }

    /// The caller has excluded queries and reclaimed the entire source.
    pub(crate) fn reset(&self) {
        self.0.store(0, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hit_rate_uses_request_count_instead_of_total_hits() {
        let frequent_cold = ReuseStats::default();
        let rare_hot = ReuseStats::default();
        assert_eq!(frequent_cold.score(), 0);
        frequent_cold.record(10, 100);
        rare_hot.record(1, 1);
        assert_eq!(frequent_cold.score(), 100_000);
        assert_eq!(rare_hot.score(), 1_000_000);
        frequent_cold.record(0, 0);
        assert_eq!(frequent_cold.score(), 100_000);
        frequent_cold.record(15, 25);
        assert_eq!(frequent_cold.score(), 200_000);
        frequent_cold.reset();
        assert_eq!(frequent_cold.score(), 0);
    }

    #[test]
    fn overflow_decays_both_counts_without_wrapping_or_sticking() {
        let stats = ReuseStats(AtomicU64::new((COUNT_MAX << 32) | (COUNT_MAX / 2)));
        stats.record(1, 1);
        assert_eq!(stats.score(), 500_000);
        let packed = stats.0.load(Ordering::Relaxed);
        assert_eq!(packed >> 32, 1u64 << 31);
        assert_eq!(packed & COUNT_MAX, 1u64 << 30);
        stats.record(usize::MAX / 2, usize::MAX);
        assert!((499_999..=500_001).contains(&stats.score()));
        stats.reset();
        stats.record(usize::MAX, usize::MAX);
        assert_eq!(stats.score(), 1_000_000);
    }

    #[test]
    fn concurrent_root_batches_keep_a_coherent_denominator() {
        let stats = ReuseStats::default();
        std::thread::scope(|scope| {
            for _ in 0..4 {
                let stats = &stats;
                scope.spawn(move || {
                    for _ in 0..1000 {
                        stats.record(1, 2);
                    }
                });
            }
        });
        assert_eq!(stats.score(), 500_000);
        assert_eq!(stats.0.load(Ordering::Relaxed) >> 32, 8000);
    }
}
