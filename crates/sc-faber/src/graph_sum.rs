//! 独立稳定图的流式求值。数学预计算由调用方完成，工作线程只借用；
//! 每个线程保留自己的配对与后代缓存，外层并行期间禁用内层线程池。
use crate::classes::Pairing;
use crate::graph::StableGraph;
use crate::kppoly::KappaPsiPolynomial as KP;
use crate::taut::DecStratum;
use sc_num::{Q, Zero};

pub(crate) type Emit<'a> = dyn FnMut(&StableGraph, KP) -> Result<(), String> + 'a;
pub(crate) type Visit<'a> = dyn Fn(&StableGraph, &mut Emit<'_>) -> Result<(), String> + Sync + 'a;

pub(crate) fn collect(
    g: u32,
    n: u32,
    d: u32,
    visit: &Visit<'_>,
) -> Result<Vec<DecStratum>, String> {
    // 完整类保留原生成顺序；不为并行而提前物化装饰。
    let mut terms = Vec::new();
    for c in 0..=d {
        for graph in crate::generators::layer(g, n, c).iter() {
            visit(graph, &mut |graph, poly| {
                terms.push(DecStratum {
                    graph: graph.clone(),
                    poly,
                });
                Ok(())
            })?;
        }
    }
    Ok(terms)
}

pub(crate) fn integrate(
    g: u32,
    n: u32,
    d: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    visit: &Visit<'_>,
) -> Result<Q, String> {
    integrate_with_threads(
        g,
        n,
        d,
        kappa,
        psi,
        visit,
        sc_num::parallel::thread_count("SC_FABER_THREADS", 16),
    )
}

fn integrate_with_threads(
    g: u32,
    n: u32,
    d: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    visit: &Visit<'_>,
    threads: usize,
) -> Result<Q, String> {
    let layers: Vec<_> = (0..=d).map(|c| crate::generators::layer(g, n, c)).collect();
    let graphs: Vec<_> = layers.iter().flat_map(|layer| layer.iter()).collect();
    crate::cache_source::with_inherited_registry(|context| {
        if graphs.len() < 8 {
            // 少量图不建立外层并行区，让昂贵的单顶点后代查询自行并行。
            let mut pairing = Pairing::new(kappa, psi);
            let mut total = Q::ZERO;
            for graph in graphs {
                visit(graph, &mut |graph, poly| {
                    total += pairing.evaluate(graph, &poly)?;
                    Ok(())
                })?;
            }
            return Ok(total);
        }
        // 小块动态领取，避免图的圈秩、装饰数差异导致静态切分失衡。
        sc_num::parallel::chunks(
            graphs.len(),
            1,
            threads,
            || Pairing::new(kappa, psi),
            |pairing, range| {
                context.run(|| {
                    let mut total = Q::ZERO;
                    for i in range {
                        visit(graphs[i], &mut |graph, poly| {
                            total += pairing.evaluate(graph, &poly)?;
                            Ok(())
                        })?;
                    }
                    Ok(total)
                })
            },
        )
        .into_iter()
        .try_fold(Q::ZERO, |sum, part: Result<Q, String>| Ok(sum + part?))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn graph_pool_inherits_the_callers_owner_for_real_socle_results() {
        assert!(
            (0..=4)
                .map(|c| crate::generators::layer(2, 1, c).len())
                .sum::<usize>()
                >= 8
        );
        let registry = std::sync::Arc::new(crate::cache_source::Registry::new());
        let first_visits = std::sync::atomic::AtomicUsize::new(0);
        let barrier = std::sync::Barrier::new(4);
        let visit = |graph: &StableGraph, emit: &mut Emit<'_>| {
            if first_visits.fetch_add(1, std::sync::atomic::Ordering::Relaxed) < 4 {
                barrier.wait();
            }
            assert_eq!(
                crate::cache_source::current_registry().unwrap().id(),
                registry.id()
            );
            crate::generators::visit_monomials(
                graph,
                4 - graph.codim(),
                false,
                Some(&[]),
                &mut |poly| emit(graph, poly),
            )
        };
        let actual = crate::cache_source::with_registry(registry.clone(), || {
            integrate_with_threads(2, 1, 4, &[], &[], &visit, 4)
        })
        .unwrap();
        assert_ne!(actual, Q::ZERO);
        assert!(first_visits.load(std::sync::atomic::Ordering::Relaxed) >= 4);
        assert!(
            registry
                .sources()
                .iter()
                .any(|source| source.ready() > source.floor())
        );
    }

    #[test]
    fn graph_workers_match_serial_pairing_and_propagate_errors() {
        let visit = |graph: &StableGraph, emit: &mut Emit<'_>| {
            assert_eq!(sc_num::parallel::thread_count("SC_FABER_THREADS", 16), 1);
            crate::generators::visit_monomials(
                graph,
                4 - graph.codim(),
                false,
                Some(&[]),
                &mut |poly| emit(graph, poly),
            )
        };
        let one = integrate_with_threads(2, 1, 4, &[], &[], &visit, 1).unwrap();
        assert_ne!(one, Q::ZERO);
        assert_eq!(
            one,
            integrate_with_threads(2, 1, 4, &[], &[], &visit, 4).unwrap()
        );
        let fail = |_: &StableGraph, _: &mut Emit<'_>| Err("graph failure".to_string());
        assert_eq!(
            integrate_with_threads(2, 1, 4, &[], &[], &fail, 4),
            Err("graph failure".into())
        );
    }
}
