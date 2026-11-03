//! 按曲线次数共享分拆、钩长和全部后代因子。查询只借用不可变表，
//! 延长阶数时在锁外建表，再发布新的 Arc；求和期间不持有全局锁。

use sc_num::{One, Q, Zero};
use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock, RwLock};

struct Row {
    part: Vec<u32>,
    hook: Q,
    factors: Vec<Q>,
}

struct Spectrum {
    degree: u32,
    max_tau: u32,
    rows: Vec<Row>,
}

pub(crate) struct Stationary {
    spectrum: Arc<Spectrum>,
    weights: Arc<Vec<Q>>,
}

impl Stationary {
    pub(crate) fn len(&self) -> usize {
        self.weights.len()
    }

    pub(crate) fn row(&self, i: usize) -> (&Q, &[Q]) {
        (&self.weights[i], &self.spectrum.rows[i].factors)
    }

    pub(crate) fn sum(
        &self,
        operations_per_row: usize,
        evaluate: impl Fn(&Q, &[Q]) -> Q + Sync,
    ) -> Q {
        let threads = if self.len().saturating_mul(operations_per_row) < 16_384 {
            1
        } else {
            threads()
        };
        self.sum_with_threads(threads, evaluate)
    }

    pub(crate) fn sum_with_threads(
        &self,
        threads: usize,
        evaluate: impl Fn(&Q, &[Q]) -> Q + Sync,
    ) -> Q {
        let values = sc_num::parallel::chunks(
            self.len(),
            128,
            threads,
            || (),
            |_, range| {
                range.fold(Q::ZERO, |sum, i| {
                    let (weight, factors) = self.row(i);
                    sum + evaluate(weight, factors)
                })
            },
        );
        values.into_iter().fold(Q::ZERO, |sum, value| sum + value)
    }
}

fn threads() -> usize {
    sc_num::parallel::thread_count("SC_GWALL_THREADS", 16)
}

pub(crate) fn prepare(g: u32, d: u32, max_tau: u32) -> Stationary {
    prepare_with_threads(g, d, max_tau, threads())
}

fn prepare_with_threads(g: u32, d: u32, max_tau: u32, threads: usize) -> Stationary {
    static SPECTRA: OnceLock<RwLock<BTreeMap<u32, Arc<Spectrum>>>> = OnceLock::new();
    static WEIGHTS: OnceLock<RwLock<BTreeMap<(u32, u32), Arc<Vec<Q>>>>> = OnceLock::new();
    let spectra = SPECTRA.get_or_init(|| RwLock::new(BTreeMap::new()));
    let previous = spectra.read().unwrap().get(&d).cloned();
    let spectrum = if let Some(previous) = previous.as_ref().filter(|p| p.max_tau >= max_tau) {
        Arc::clone(previous)
    } else {
        let parts = previous.as_ref().map_or_else(
            || crate::comb::partition(d),
            |p| p.rows.iter().map(|row| row.part.clone()).collect(),
        );
        let corrections: Vec<_> = (0..=max_tau)
            .map(|k| crate::qseries::z_over_2sinh_coeff(k + 2))
            .collect();
        let workers = if parts.len().saturating_mul(max_tau as usize + 1) < 4096 {
            1
        } else {
            threads
        };
        let rows = sc_num::parallel::chunks(
            parts.len(),
            64,
            workers,
            || (),
            |_, range| {
                range
                    .map(|i| {
                        let part = &parts[i];
                        let hook = previous.as_ref().map_or_else(
                            || crate::gwpoint::hook_product(part),
                            |p| p.rows[i].hook.clone(),
                        );
                        Row {
                            part: part.clone(),
                            hook,
                            factors: factors(part, &corrections),
                        }
                    })
                    .collect::<Vec<_>>()
            },
        )
        .into_iter()
        .flatten()
        .collect();
        let built = Arc::new(Spectrum {
            degree: d,
            max_tau,
            rows,
        });
        let mut cache = spectra.write().unwrap();
        // 并发查询可能已经发布了更长的表。
        Arc::clone(
            cache
                .entry(d)
                .and_modify(|row| {
                    if row.max_tau < max_tau {
                        *row = Arc::clone(&built);
                    }
                })
                .or_insert(built),
        )
    };
    let weights = WEIGHTS.get_or_init(|| RwLock::new(BTreeMap::new()));
    let cached = weights.read().unwrap().get(&(g, d)).cloned();
    let weights = cached.unwrap_or_else(|| {
        let exponent = 2 * i64::from(g) - 2;
        let values = sc_num::parallel::chunks(
            spectrum.rows.len(),
            128,
            if spectrum.rows.len() < 4096 {
                1
            } else {
                threads
            },
            || (),
            |_, range| {
                range
                    .map(|i| {
                        let value =
                            sc_num::pow_unsigned(&spectrum.rows[i].hook, exponent.unsigned_abs());
                        if exponent < 0 { Q::ONE / value } else { value }
                    })
                    .collect::<Vec<_>>()
            },
        )
        .into_iter()
        .flatten()
        .collect();
        let built = Arc::new(values);
        Arc::clone(
            weights
                .write()
                .unwrap()
                .entry((g, spectrum.degree))
                .or_insert(built),
        )
    });
    Stationary { spectrum, weights }
}

fn factors(part: &[u32], corrections: &[Q]) -> Vec<Q> {
    // 直接递推 a^l/l!，每个阶数只乘 a/l，避免反复求幂和查表。
    let bases: Vec<_> = part
        .iter()
        .rev()
        .enumerate()
        .map(|(i, &p)| {
            let b = Q::from_signeds(-(2 * i as i64 + 1), 2);
            (Q::from(p) + &b, b)
        })
        .collect();
    let mut powers = vec![(Q::ONE, Q::ONE); bases.len()];
    corrections
        .iter()
        .enumerate()
        .map(|(k, correction)| {
            let divisor = Q::from(k + 1);
            let mut value = correction.clone();
            for ((a, b), (pa, pb)) in bases.iter().zip(&mut powers) {
                *pa = &*pa * a / &divisor;
                *pb = &*pb * b / &divisor;
                value += &*pa - &*pb;
            }
            value
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_tables_match_scalar_formula_and_extend_safely() {
        for d in (0..=12).chain([20]) {
            let short = prepare_with_threads(2, d, 2, 1);
            let long = std::thread::scope(|s| {
                let a = s.spawn(|| prepare_with_threads(2, d, 7, 4));
                let b = s.spawn(|| prepare_with_threads(2, d, 5, 2));
                b.join().unwrap();
                a.join().unwrap()
            });
            for (i, part) in crate::comb::partition(d).iter().enumerate() {
                assert_eq!(short.row(i).1, &long.row(i).1[..short.row(i).1.len()]);
                for k in 0..=7 {
                    assert_eq!(
                        long.row(i).1[k],
                        crate::gwpoint::power(part, k as u32 + 1) / sc_num::factorial(k as u32 + 1)
                    );
                }
            }
            let evaluate = |w: &Q, row: &[Q]| w * &row[2] * &row[4];
            assert_eq!(
                long.sum_with_threads(1, evaluate),
                long.sum_with_threads(4, evaluate)
            );
        }
    }
}
