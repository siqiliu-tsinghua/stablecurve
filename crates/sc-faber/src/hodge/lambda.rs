//! Ordered batches of pure-lambda top-degree integrals.
//!
//! A single integral benefits from bounded normal-form reduction. Enumerating
//! every monomial admits a cheaper triangular algorithm: replace one highest
//! repeated index and reuse integrals that have already been computed.

use sc_num::{Q, Zero};
use std::collections::BTreeMap;

fn degree(g: u32) -> Result<u32, String> {
    g.checked_sub(1)
        .filter(|_| g >= 2)
        .and_then(|value| value.checked_mul(3))
        .ok_or_else(|| "require stable unmarked genus and 3g-3 fitting u32".to_owned())
}

/// Count the monomials in lambda_1,...,lambda_g of weight 3g-3.
///
/// The count includes monomials whose integral vanishes. Genus must be at
/// least two. An unrepresentable degree, allocation or count returns an error.
pub fn lambda_monomial_count(g: u32) -> Result<u64, String> {
    let degree = usize::try_from(degree(g)?).map_err(|_| "degree exceeds usize")?;
    let length = degree.checked_add(1).ok_or("count allocation overflow")?;
    let mut counts = Vec::<u64>::new();
    counts
        .try_reserve_exact(length)
        .map_err(|error| format!("allocate monomial counts: {error}"))?;
    counts.resize(length, 0);
    counts[0] = 1;
    for part in 1..=usize::try_from(g).map_err(|_| "genus exceeds usize")? {
        for weight in part..=degree {
            counts[weight] = counts[weight]
                .checked_add(counts[weight - part])
                .ok_or("monomial count exceeds u64")?;
        }
    }
    Ok(counts[degree])
}

/// Compute every pure-lambda top-degree integral on unmarked genus-g curves.
///
/// `parts` is a descending list of lambda indices with sum 3g-3. The callback
/// receives every monomial, including zeros, in descending lexicographic order;
/// its borrowed arguments are valid only during that call. Returning an error
/// stops the batch immediately. The successful result is the number visited.
///
/// Mumford relations reuse earlier answers; only distinct-index monomials need
/// the general `klp` entry. Completed integrals use the current cache session's
/// ordinary root source, including values obtained solely from relations.
/// Bind a `CacheSession::run` scope to persist and resume backend work. Callback
/// code runs without result-source locks and may make other mathematical calls.
///
/// Answers for the current genus remain in memory until the batch returns, so
/// memory usage grows with the monomial count. The caller owns output, progress
/// reporting, cache flushing and closing. Genus must be at least two.
pub fn lambda_integrals(
    g: u32,
    mut visit: impl FnMut(&[u32], &Q) -> Result<(), String>,
) -> Result<u64, String> {
    let degree = degree(g)?;
    let calculate = || {
        crate::cache_source::catch_resource(|| {
            let mut previous = BTreeMap::new();
            let mut visited = 0u64;
            visit_remaining(degree, g, &mut Vec::new(), &mut |parts| {
                let powers = lambda_powers(parts);
                let value = if let Some(value) = relation_value(g, parts, &previous)? {
                    super::root_source::complete_lambda(g, &powers, value)?
                } else {
                    super::klp(g, 0, &[], &[], &powers, &[])?
                };
                // Publish before the callback: even an output failure must leave
                // completed mathematics available to the background consumer.
                visit(parts, &value)?;
                previous.insert(parts.to_vec(), value);
                visited = visited.checked_add(1).ok_or("visited count exceeds u64")?;
                Ok(())
            })?;
            Ok(visited)
        })
    };
    match crate::cache_source::current_registry() {
        Some(registry) => crate::cache_source::try_with_registry(registry, calculate)?,
        None => crate::cache_source::with_optional_registry(None, calculate),
    }
}

fn lambda_powers(parts: &[u32]) -> Vec<(u32, u32)> {
    let mut powers: Vec<(u32, u32)> = Vec::new();
    for &index in parts.iter().rev() {
        if let Some(last) = powers.last_mut().filter(|last| last.0 == index) {
            last.1 += 1;
        } else {
            powers.push((index, 1));
        }
    }
    powers
}

fn relation_value(
    g: u32,
    parts: &[u32],
    previous: &BTreeMap<Vec<u32>, Q>,
) -> Result<Option<Q>, String> {
    let Some(pair) = parts.windows(2).position(|pair| pair[0] == pair[1]) else {
        return Ok(None);
    };
    let mut value = Q::ZERO;
    for (j, larger, coefficient) in super::relations::square_terms(g, parts[pair]) {
        let mut rhs = Vec::with_capacity(parts.len());
        rhs.extend_from_slice(&parts[..pair]);
        rhs.extend_from_slice(&parts[pair + 2..]);
        if j != 0 {
            rhs.push(j);
        }
        rhs.push(larger);
        rhs.sort_unstable_by(|a, b| b.cmp(a));
        let earlier = previous.get(rhs.as_slice()).ok_or_else(|| {
            format!("missing earlier Mumford term at genus {g}: {parts:?} -> {rhs:?}")
        })?;
        value += earlier * Q::from(coefficient);
    }
    Ok(Some(value))
}

fn visit_remaining(
    remaining: u32,
    cap: u32,
    parts: &mut Vec<u32>,
    visit: &mut impl FnMut(&[u32]) -> Result<(), String>,
) -> Result<(), String> {
    if remaining == 0 {
        return visit(parts);
    }
    for next in (1..=remaining.min(cap)).rev() {
        parts.push(next);
        let result = visit_remaining(remaining - next, next, parts, visit);
        parts.pop();
        result?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache_source::{IntegralKey, Registry, with_registry};
    use std::sync::Arc;

    #[test]
    fn batches_match_the_original_backend_through_genus_six() {
        crate::cache_source::with_optional_registry(None, || {
            for g in 2..=6 {
                let mut last: Option<Vec<u32>> = None;
                let count = lambda_integrals(g, |parts, value| {
                    assert!(parts.windows(2).all(|pair| pair[0] >= pair[1]));
                    assert_eq!(parts.iter().sum::<u32>(), 3 * g - 3);
                    if let Some(last) = &last {
                        assert!(last.as_slice() > parts);
                    }
                    last = Some(parts.to_vec());
                    let expected =
                        super::super::klp_unreduced(g, 0, &[], &[], &lambda_powers(parts), &[])?;
                    assert_eq!(*value, expected, "g={g}, parts={parts:?}");
                    Ok(())
                })
                .unwrap();
                assert_eq!(count, lambda_monomial_count(g).unwrap());
            }
        });
    }

    #[test]
    fn counts_and_invalid_genera_preserve_the_enumeration_contract() {
        assert_eq!(lambda_monomial_count(2), Ok(2));
        assert_eq!(lambda_monomial_count(3), Ok(7));
        assert_eq!(lambda_monomial_count(20), Ok(516_054));
        for g in [0, 1, u32::MAX] {
            assert!(lambda_monomial_count(g).is_err());
            assert!(lambda_integrals(g, |_, _| panic!("invalid genus must not visit")).is_err());
        }
    }

    #[test]
    fn callback_error_stops_before_any_expensive_genus_twenty_query() {
        crate::cache_source::with_optional_registry(None, || {
            let mut visited = 0;
            assert_eq!(
                lambda_integrals(20, |parts, value| {
                    visited += 1;
                    assert_eq!(parts, [20, 20, 17]);
                    assert_eq!(*value, Q::ZERO);
                    Err("stop after the first zero".to_owned())
                }),
                Err("stop after the first zero".to_owned())
            );
            assert_eq!(visited, 1);
        });
    }

    #[test]
    fn callback_can_switch_owners_without_changing_batch_result_ownership() {
        let outer = Arc::new(Registry::new());
        let inner = Arc::new(Registry::new());
        with_registry(outer.clone(), || {
            assert_eq!(
                lambda_integrals(2, |_, _| {
                    with_registry(inner.clone(), || {
                        assert_eq!(
                            super::super::klp(1, 1, &[], &[], &[(1, 1)], &[]).unwrap(),
                            Q::from_signeds(1, 24)
                        );
                    });
                    Ok(())
                }),
                Ok(2)
            );
        });
        let rows = |registry: &Registry| {
            registry
                .sources()
                .into_iter()
                .flat_map(|source| source.pull(0, usize::MAX).unwrap().rows)
                .collect::<Vec<_>>()
        };
        let outer_rows = rows(&outer);
        assert_eq!(outer_rows.len(), 2);
        assert!(
            outer_rows
                .iter()
                .all(|row| matches!(row.key, IntegralKey::Integral { g: 2, n: 0, .. }))
        );
        let inner_rows = rows(&inner);
        assert_eq!(inner_rows.len(), 1);
        assert!(matches!(
            inner_rows[0].key,
            IntegralKey::Integral { g: 1, n: 1, .. }
        ));
    }

    #[test]
    fn relation_only_answers_persist_as_ordinary_integral_roots() {
        use crate::socle::CacheSession;
        use std::path::PathBuf;
        use std::sync::atomic::{AtomicUsize, Ordering};

        static NEXT_DATABASE: AtomicUsize = AtomicUsize::new(0);
        struct Database(PathBuf);
        impl Drop for Database {
            fn drop(&mut self) {
                for suffix in ["", "-wal", "-shm"] {
                    let mut path = self.0.as_os_str().to_os_string();
                    path.push(suffix);
                    let _ = std::fs::remove_file(PathBuf::from(path));
                }
            }
        }
        let database = Database(std::env::temp_dir().join(format!(
            "sc-faber-lambda-batch-{}-{}.sqlite",
            std::process::id(),
            NEXT_DATABASE.fetch_add(1, Ordering::Relaxed)
        )));
        let session = CacheSession::open(database.0.to_str().unwrap()).unwrap();
        session.run(|| lambda_integrals(2, |_, _| Ok(()))).unwrap();
        session.flush().unwrap();
        session.close().unwrap();
        let connection = rusqlite::Connection::open(&database.0).unwrap();
        let count: i64 = connection
            .query_row("SELECT COUNT(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 2);
        drop(connection);
        let reopened = CacheSession::open(database.0.to_str().unwrap()).unwrap();
        assert_eq!(
            reopened
                .run(|| super::super::klp(2, 0, &[], &[], &[(1, 3)], &[]))
                .unwrap(),
            Q::from_signeds(1, 2880)
        );
        reopened.close().unwrap();
    }
}
