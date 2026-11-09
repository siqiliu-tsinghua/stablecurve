//! Bounded normal-form reduction using c(E)c(E*) = 1.
//!
//! The identity holds on every marked moduli space, so kappa, psi and ch
//! factors can remain outside this calculation. Sparse lambda powers avoid
//! allocating one entry for each factor of a high power.

use sc_num::{One, Q, Zero};
use std::collections::{BTreeMap, btree_map::Entry};

type Powers = Vec<(u32, u32)>;

// A complete squarefree expansion of lambda_1^(3g-3) can visit essentially the
// whole partition DAG before the existing shared GRR evaluator even starts.
// These bounds keep speculative simplification small. On exhaustion the caller
// evaluates the original monomial, never a partially simplified approximation.
const MAX_EXPANDED_STATES: usize = 4096;
const MAX_TERMS: usize = 512;

/// One square relation, shared by sparse normalization and ordered batches.
/// lambda_0 is the unit; callers omit its factor when j is zero.
pub(super) fn square_terms(g: u32, index: u32) -> impl ExactSizeIterator<Item = (u32, u32, i32)> {
    debug_assert!(index > 0 && index <= g);
    let twice = 2 * u64::from(index);
    let first = twice.saturating_sub(u64::from(g)) as u32;
    (first..index).map(move |j| {
        // The clipped range guarantees 2i-j <= g, hence fits u32.
        let larger = (twice - u64::from(j)) as u32;
        let coefficient = if (u64::from(index) + u64::from(j) + 1) % 2 == 0 {
            2
        } else {
            -2
        };
        (j, larger, coefficient)
    })
}

/// Reduce a normalized ascending sparse lambda monomial to squarefree terms.
///
/// `None` means no square was present, or complete reduction exceeded the
/// bounded work allowance. `Some([])` is an exact zero. Unit lambda_0 factors
/// and zero powers are ignored. Coefficients are collected before a node is
/// expanded; this is a per-call polynomial, not a persistent memo table.
pub(super) fn reduce(g: u32, lambda: &[(u32, u32)]) -> Option<Vec<(Powers, Q)>> {
    debug_assert!(lambda.windows(2).all(|pair| pair[0].0 < pair[1].0));
    if lambda
        .iter()
        .any(|&(index, power)| power != 0 && (index > g || (index == g && g > 0 && power >= 2)))
    {
        return Some(Vec::new());
    }
    let lambda: Powers = lambda
        .iter()
        .copied()
        .filter(|&(index, power)| index != 0 && power != 0)
        .collect();
    if !lambda.iter().any(|&(_, power)| power >= 2) {
        return None;
    }

    let energy = lambda.iter().try_fold(0u128, |sum, &(index, power)| {
        let index = u128::from(index);
        sum.checked_add(index.checked_mul(index)?.checked_mul(u128::from(power))?)
    })?;
    let mut pending = BTreeMap::from([((energy, lambda), Q::ONE)]);
    let mut output = BTreeMap::new();
    let mut expanded = 0;

    while let Some(((energy, powers), coefficient)) = pending.pop_first() {
        let Some(pair) = powers.iter().rposition(|&(_, power)| power >= 2) else {
            accumulate(&mut output, powers, coefficient);
            continue;
        };
        if expanded == MAX_EXPANDED_STATES {
            return None;
        }
        expanded += 1;

        let index = powers[pair].0;
        if index == g {
            continue;
        }
        let terms = square_terms(g, index);
        if terms.len() > MAX_TERMS {
            return None;
        }
        let mut remaining = powers;
        remaining[pair].1 -= 2;
        if remaining[pair].1 == 0 {
            remaining.remove(pair);
        }

        // lambda_i^2 = 2 sum_j (-1)^(i+j+1) lambda_j lambda_(2i-j).
        // Only one pair is replaced, even when its original power exceeds two.
        for (j, larger, sign) in terms {
            let mut term = remaining.clone();
            if j != 0 {
                increment(&mut term, j)?;
            }
            increment(&mut term, larger)?;
            // Every edge strictly raises sum(i^2 * power) by 2(i-j)^2.
            // Popping the smallest energy therefore collects all incoming
            // coefficients before expanding a monomial, including paths of
            // different lengths. The secondary sparse key makes ties stable.
            let difference = u128::from(index - j);
            let next_energy = energy.checked_add(2 * difference * difference)?;
            accumulate(
                &mut pending,
                (next_energy, term),
                &coefficient * Q::from(sign),
            );
            if pending.len().checked_add(output.len())? > MAX_TERMS {
                return None;
            }
        }
    }
    Some(output.into_iter().collect())
}

fn increment(powers: &mut Powers, index: u32) -> Option<()> {
    match powers.binary_search_by_key(&index, |&(index, _)| index) {
        Ok(position) => powers[position].1 = powers[position].1.checked_add(1)?,
        Err(position) => powers.insert(position, (index, 1)),
    }
    Some(())
}

fn accumulate<K: Ord>(terms: &mut BTreeMap<K, Q>, key: K, coefficient: Q) {
    if coefficient == Q::ZERO {
        return;
    }
    match terms.entry(key) {
        Entry::Occupied(mut entry) => {
            *entry.get_mut() += coefficient;
            if *entry.get() == Q::ZERO {
                entry.remove();
            }
        }
        Entry::Vacant(entry) => {
            entry.insert(coefficient);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signs_rank_and_pair_multiplicities_are_exact() {
        assert_eq!(
            reduce(4, &[(2, 2)]),
            Some(vec![
                (vec![(1, 1), (3, 1)], Q::from(2)),
                (vec![(4, 1)], Q::from(-2)),
            ])
        );
        assert_eq!(
            reduce(3, &[(2, 3)]),
            Some(vec![(vec![(1, 1), (2, 1), (3, 1)], Q::from(2))])
        );
        assert_eq!(reduce(3, &[(3, 2)]), Some(vec![]));
        assert_eq!(reduce(2, &[(1, 4)]), Some(vec![]));
        assert_eq!(reduce(2, &[(3, 1)]), Some(vec![]));
    }

    #[test]
    fn converging_reduction_paths_collect_before_expansion() {
        assert_eq!(
            reduce(3, &[(1, 6)]),
            Some(vec![(vec![(1, 1), (2, 1), (3, 1)], Q::from(16))])
        );
        assert_eq!(
            reduce(2, &[(0, 17), (1, 3), (2, 0)]),
            Some(vec![(vec![(1, 1), (2, 1)], Q::from(2))])
        );
        assert_eq!(reduce(4, &[(0, 5), (1, 1), (3, 1)]), None);
        assert_eq!(reduce(0, &[(0, 5)]), None);
    }

    #[test]
    fn squarefree_normal_forms_preserve_weight_and_order() {
        for g in 2..=6 {
            for powers in [
                vec![(1, 5), (2, 3)],
                vec![(1, 2), (g - 1, 2)],
                vec![(g - 1, 4)],
            ] {
                let mut normalized = BTreeMap::<u32, u32>::new();
                for (index, power) in powers {
                    *normalized.entry(index).or_default() += power;
                }
                let powers: Powers = normalized.into_iter().collect();
                let weight = powers
                    .iter()
                    .map(|&(i, p)| u64::from(i) * u64::from(p))
                    .sum::<u64>();
                let result = reduce(g, &powers).expect("small reduction remains bounded");
                for (term, coefficient) in result {
                    assert_ne!(coefficient, Q::ZERO);
                    assert!(term.windows(2).all(|pair| pair[0].0 < pair[1].0));
                    assert!(term.iter().all(|&(i, p)| i > 0 && i <= g && p == 1));
                    assert_eq!(term.iter().map(|&(i, _)| u64::from(i)).sum::<u64>(), weight);
                }
            }
        }
    }

    #[test]
    fn normal_forms_match_independent_ch_expansions_before_rank_truncation() {
        fn expand(powers: &[(u32, u32)]) -> Vec<super::super::ChMono> {
            let mut polynomial = vec![(vec![], Q::ONE)];
            for &(index, power) in powers {
                let factor = super::super::lambda_expansion(index);
                for _ in 0..power {
                    polynomial = super::super::ch_mul(&polynomial, &factor);
                }
            }
            polynomial
        }

        // Weight at most rank means no relation here discards a lambda above
        // the rank. Both sides can be compared as formal odd-ch polynomials.
        for powers in [
            vec![(1, 8)],
            vec![(2, 4)],
            vec![(1, 2), (3, 2)],
            vec![(2, 2), (4, 1)],
            vec![(4, 2)],
        ] {
            let mut simplified = Vec::new();
            for (term, coefficient) in reduce(8, &powers).unwrap() {
                simplified.extend(
                    expand(&term)
                        .into_iter()
                        .map(|(ch, value)| (ch, coefficient.clone() * value)),
                );
            }
            assert_eq!(
                super::super::normalize_ch(simplified),
                expand(&powers),
                "lambda powers {powers:?}"
            );
        }
    }

    #[test]
    fn large_lambda_one_normal_form_falls_back_without_grr_work() {
        assert_eq!(reduce(20, &[(1, 57)]), None);
        // The rank-zero rule precedes any bounded symbolic work.
        assert_eq!(reduce(20, &[(1, 57), (20, 2)]), Some(vec![]));
    }
}
