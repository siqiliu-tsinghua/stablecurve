//! Hodge evaluations that do not construct a descendant work graph.
//!
//! Faber--Pandharipande, *Hodge integrals and Gromov--Witten theory*,
//! Theorem 4 and Section 4.3: <https://arxiv.org/abs/math/9810173>.
//! The general lambda_g lambda_(g-1) pairing uses *Logarithmic Series and
//! Hodge Integrals*, Sections 0.4--0.7 and Theorem 1:
//! <https://webspace.science.uu.nl/~faber010/MMJ2000fp.pdf>.
//! The explicit Bernoulli descendant formula is Liu--Xu, Theorem 10:
//! <https://arxiv.org/abs/0705.3564>.

use sc_num::{One, Q, Zero};
use std::collections::BTreeMap;

// Inputs are normalized by klp. lambda_0 is the unit and is never expanded.
fn direct_factor(g: u32, lambda: &[(u32, u32)]) -> Option<u32> {
    if g < 2 {
        return None;
    }
    let factors = lambda.iter().copied().filter(|&(index, _)| index != 0);
    if (g == 2 && factors.clone().eq([(1, 1), (2, 1)]))
        || (g > 2 && factors.clone().eq([(g - 2, 1), (g - 1, 1), (g, 1)]))
    {
        Some(1)
    } else if factors.eq([(g - 1, 3)]) {
        Some(2)
    } else {
        None
    }
}

/// These fixed formulas are cheaper than requesting a historical SQLite root.
pub(super) fn inexpensive(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> bool {
    n == 0
        && kappa.is_empty()
        && psi.is_empty()
        && ch.is_empty()
        && direct_factor(g, lambda).is_some()
}

/// Try exact evaluations before the general GRR recursion.
///
/// The triple and cubic formulas apply on the unmarked space. For the general
/// lambda_g lambda_(g-1) family, Hodge classes pull back under forgetting marks:
/// first push the kappa/psi factor to (g,0), then use boundary annihilation and
/// Mumford GRR there. In particular, we do not discard rational-tail boundary
/// classes on pointed spaces. None delegates to the ordinary evaluator.
pub(super) fn value(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Option<Q> {
    if g < 2 {
        return None;
    }
    let weight = kappa
        .iter()
        .chain(lambda)
        .chain(ch)
        .map(|&(index, power)| u128::from(index) * u128::from(power))
        .sum::<u128>()
        + psi
            .iter()
            .map(|&(_, power)| u128::from(power))
            .sum::<u128>();
    if weight != 3 * u128::from(g) - 3 + u128::from(n) {
        return Some(Q::ZERO);
    }
    let hodge_weight = lambda
        .iter()
        .chain(ch)
        .map(|&(index, power)| u128::from(index) * u128::from(power))
        .sum::<u128>();
    if hodge_weight > 3 * u128::from(g) - 3 {
        return Some(Q::ZERO);
    }
    if inexpensive(g, n, kappa, psi, lambda, ch) {
        return top_triple(g).map(|value| Q::from(direct_factor(g, lambda).unwrap()) * value);
    }
    let top_power = lambda.iter().find(|&&(index, _)| index == g)?.1;
    lambda.iter().find(|&&(index, _)| index == g - 1)?;
    if top_power >= 2 {
        return Some(Q::ZERO);
    }
    // Factorial indices and the existing marked pushforward use u32.
    let twice_g = g.checked_mul(2)?;
    twice_g.checked_add(n)?;

    let marked = super::fund_poly(g, n, kappa, psi).ok()?;
    let pushed = crate::push::push(g, n, &marked).ok()?;
    if pushed.terms.is_empty() {
        return Some(Q::ZERO);
    }

    let remaining = lambda
        .iter()
        .filter(|&&(index, _)| index != 0)
        .flat_map(|&(index, power)| {
            let consumed = u32::from(index == g || index == g - 1);
            std::iter::repeat_n(index, power.saturating_sub(consumed) as usize)
        })
        .collect::<Vec<_>>();
    let lambda_polynomial = lambda_to_kappa(&remaining);
    let mut ch_coefficient = Q::ONE;
    let mut ch_kappa = Vec::new();
    for &(degree, power) in ch {
        let grr = crate::socle::bernoulli(u64::from(degree) + 1)
            / sc_num::factorial(degree.checked_add(1)?);
        ch_coefficient *= sc_num::pow_unsigned(&grr, u64::from(power));
        ch_kappa.extend(std::iter::repeat_n(degree, power as usize));
    }
    let mut total = Q::ZERO;
    for term in pushed.terms {
        let mut base = term.kappa.into_iter().next().unwrap_or_default();
        base.extend_from_slice(&ch_kappa);
        for (remaining, coefficient) in &lambda_polynomial {
            let mut indices = base.clone();
            indices.extend_from_slice(remaining);
            indices.sort_unstable();
            total += &term.coeff * coefficient * kappa_pairing(g, &indices)?;
        }
    }
    let bernoulli = crate::socle::bernoulli(u64::from(twice_g));
    let bernoulli = if bernoulli < Q::ZERO {
        -bernoulli
    } else {
        bernoulli
    };
    Some(
        total * ch_coefficient * bernoulli
            / (sc_num::pow_unsigned(&Q::from(2), u64::from(twice_g - 1))
                * sc_num::factorial(twice_g)),
    )
}

fn top_triple(g: u32) -> Option<Q> {
    let twice_g = g.checked_mul(2)?;
    // Consecutive even Bernoulli numbers have opposite signs.
    let numerator = -(crate::socle::bernoulli(u64::from(twice_g))
        * crate::socle::bernoulli(u64::from(twice_g - 2)));
    let denominator =
        Q::from(2) * Q::from(twice_g) * Q::from(twice_g - 2) * sc_num::factorial(twice_g - 2);
    Some(numerator / denominator)
}

// Boundary terms vanish after pairing with lambda_g lambda_(g-1) on unmarked Mbar_g.
// Keys are sorted kappa indices, with repeated factors retained.
fn lambda_to_kappa(parts: &[u32]) -> BTreeMap<Vec<u32>, Q> {
    let mut polynomial = BTreeMap::from([(Vec::new(), Q::ONE)]);
    for &index in parts {
        let factor: Vec<_> = super::lambda_expansion(index)
            .into_iter()
            .map(|(chern, mut coefficient)| {
                let mut kappa = Vec::new();
                for (degree, power) in chern {
                    let grr = crate::socle::bernoulli(u64::from(degree) + 1)
                        / sc_num::factorial(degree + 1);
                    coefficient *= sc_num::pow_unsigned(&grr, u64::from(power));
                    kappa.extend(std::iter::repeat_n(degree, power as usize));
                }
                kappa.sort_unstable();
                (kappa, coefficient)
            })
            .collect();
        let mut product = BTreeMap::<Vec<u32>, Q>::new();
        for (left, a) in &polynomial {
            for (right, b) in &factor {
                let mut kappa = left.clone();
                kappa.extend_from_slice(right);
                kappa.sort_unstable();
                *product.entry(kappa).or_insert(Q::ZERO) += a * b;
            }
        }
        product.retain(|_, coefficient| *coefficient != Q::ZERO);
        polynomial = product;
    }
    polynomial
}

// The common |B_(2g)|/[2^(2g-1)(2g)!] is applied by the caller. The multiset
// multiplicity counts labeled set partitions without enumerating them all.
fn kappa_pairing(g: u32, kappa: &[u32]) -> Option<Q> {
    let mut result = Q::ZERO;
    for (blocks, multiplicity) in crate::socle::setparts_with_auts(kappa) {
        let factorial_index = g
            .checked_mul(2)?
            .checked_sub(3)?
            .checked_add(u32::try_from(blocks.len()).ok()?)?;
        let mut term = multiplicity * sc_num::factorial(factorial_index);
        for block in &blocks {
            let sum = block
                .iter()
                .try_fold(0u32, |sum, &index| sum.checked_add(index))?;
            let last = sum.checked_mul(2)?.checked_add(1)?;
            term /= (1..=last)
                .step_by(2)
                .fold(Q::ONE, |value, odd| value * Q::from(odd));
        }
        if (kappa.len() - blocks.len()) % 2 == 0 {
            result += term;
        } else {
            result -= term;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lambda_powers(parts: &[u32]) -> Vec<(u32, u32)> {
        let mut lambda = BTreeMap::<u32, u32>::new();
        for &index in parts.iter().filter(|&&index| index != 0) {
            *lambda.entry(index).or_default() += 1;
        }
        lambda.into_iter().collect()
    }

    fn oracle(
        g: u32,
        n: u32,
        kappa: &[(u32, u32)],
        psi: &[(u32, u32)],
        lambda: &[(u32, u32)],
        ch: &[(u32, u32)],
    ) -> Q {
        super::super::klp_unreduced(g, n, kappa, psi, lambda, ch).unwrap()
    }

    #[test]
    fn triples_match_independent_small_genus_anchors_and_genus_twenty() {
        for (g, expected) in [
            (2, "1/5760"),
            (3, "1/1451520"),
            (4, "1/87091200"),
            (
                20,
                "40261619800200260590461963277216711/6793436173750188520428381832391044953538560000000000",
            ),
        ] {
            let lambda = lambda_powers(&(g - 2..=g).collect::<Vec<_>>());
            assert!(inexpensive(g, 0, &[], &[], &lambda, &[]));
            assert_eq!(
                value(g, 0, &[], &[], &lambda, &[]).unwrap().to_string(),
                expected
            );
        }
    }

    #[test]
    fn closed_families_match_the_unreduced_evaluator_through_genus_five() {
        for g in 2..=5 {
            for parts in [(g - 2..=g).collect::<Vec<_>>(), vec![g - 1; 3]] {
                let lambda = lambda_powers(&parts);
                assert_eq!(
                    value(g, 0, &[], &[], &lambda, &[]),
                    Some(oracle(g, 0, &[], &[], &lambda, &[])),
                    "g={g}, {parts:?}"
                );
            }
        }
    }

    #[test]
    fn every_squarefree_top_family_matches_the_unreduced_evaluator_through_genus_six() {
        let mut checked = 0;
        for g in 2..=6 {
            for mask in 0..1u32 << g {
                let parts: Vec<u32> = (1..=g)
                    .filter(|&index| mask & (1 << (index - 1)) != 0)
                    .collect();
                if parts.iter().sum::<u32>() != 3 * g - 3
                    || !parts.contains(&g)
                    || !parts.contains(&(g - 1))
                {
                    continue;
                }
                let lambda = lambda_powers(&parts);
                assert_eq!(
                    value(g, 0, &[], &[], &lambda, &[]),
                    Some(oracle(g, 0, &[], &[], &lambda, &[])),
                    "g={g}, {parts:?}"
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 7);
    }

    #[test]
    fn repeated_factors_preserve_multiset_pushforward_multiplicities() {
        for (g, parts) in [
            (4, vec![4, 3, 1, 1]),
            (5, vec![5, 4, 1, 1, 1]),
            (6, vec![6, 5, 2, 2]),
            (6, vec![6, 5, 1, 1, 1, 1]),
        ] {
            let lambda = lambda_powers(&parts);
            assert_eq!(
                value(g, 0, &[], &[], &lambda, &[]),
                Some(oracle(g, 0, &[], &[], &lambda, &[])),
                "g={g}, {parts:?}"
            );
        }
    }

    #[test]
    fn marked_and_mixed_pairings_match_the_unreduced_evaluator() {
        let cases: &[(
            u32,
            u32,
            &[(u32, u32)],
            &[(u32, u32)],
            &[(u32, u32)],
            &[(u32, u32)],
        )] = &[
            (3, 0, &[(1, 1)], &[], &[(2, 1), (3, 1)], &[]),
            (3, 1, &[], &[(1, 2)], &[(2, 1), (3, 1)], &[]),
            (4, 2, &[], &[(1, 2), (2, 2)], &[(3, 1), (4, 1)], &[]),
            (4, 1, &[], &[(1, 3)], &[(3, 1), (4, 1)], &[]),
            (3, 2, &[], &[(1, 3)], &[(2, 1), (3, 1)], &[]),
            (4, 2, &[(1, 1)], &[(1, 2), (2, 1)], &[(3, 1), (4, 1)], &[]),
            (4, 1, &[], &[(1, 2)], &[(3, 1), (4, 1)], &[(1, 1)]),
            (4, 0, &[(1, 1)], &[], &[(3, 1), (4, 1)], &[(1, 1)]),
            (5, 0, &[], &[], &[(4, 1), (5, 1)], &[(3, 1)]),
        ];
        for &(g, n, kappa, psi, lambda, ch) in cases {
            assert_eq!(
                value(g, n, kappa, psi, lambda, ch),
                Some(oracle(g, n, kappa, psi, lambda, ch)),
                "g={g}, n={n}, kappa={kappa:?}, psi={psi:?}, ch={ch:?}"
            );
        }
    }

    #[test]
    fn general_pairing_has_no_genus_twenty_cutoff() {
        // The one-kappa Faber formula is independent of lambda expansion.
        let g = 21;
        let lambda = [(g - 1, 1), (g, 1)];
        let actual = value(g, 0, &[(g - 2, 1)], &[], &lambda, &[]).unwrap();
        let expected = crate::socle::bernoulli(2 * u64::from(g)) * sc_num::factorial(g - 1)
            / (sc_num::pow_unsigned(&Q::from(2), u64::from(g)) * sc_num::factorial(2 * g));
        let expected = if expected < Q::ZERO {
            -expected
        } else {
            expected
        };
        assert_eq!(actual, expected);
    }
}
