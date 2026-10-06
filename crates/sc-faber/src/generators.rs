//! 重言生成元枚举：R^d(M̄_{g,n}) 的稳定图 + κ-ψ 单项式生成元。
//!
//! 对应 admcycles 的 `tautgens(g, n, deg)`。枚举分两步：
//! 1. 余维数 ≤ d 的全部连通稳定图（自环扩张 + 顶点分裂，逐余维递推，
//!    规范形去重）；2. 对每个余维 c 的图，枚举总权重 d−c 的 κ-ψ 单项式
//!    （κ 指标 ≤ 3g_v−3+n_v 逐顶点，ψ 可落任意腿）。
//!
//! 注意：枚举的是生成元（spanning set），不是基——FZ 关系（M3）才给约化。

use crate::graph::{Leg, StableGraph};
use crate::kppoly::{KappaPsiPolynomial, KpTerm};
use crate::taut::DecStratum;
use sc_num::{One, Q};
use std::collections::HashSet;

/// 枚举 (g, n) 余维数恰为 codim 的全部连通稳定图（无装饰）。
pub fn graphs(g: u32, n: u32, codim: u32) -> Vec<StableGraph> {
    let base = match StableGraph::new(&[g], &[(1..=n).collect::<Vec<Leg>>()], &[]) {
        Ok(gr) => gr,
        Err(_) => return vec![],
    };
    let mut layer: Vec<StableGraph> = vec![base];
    let mut seen: HashSet<StableGraph> = layer.iter().cloned().collect();
    for _ in 0..codim {
        let mut next: Vec<StableGraph> = Vec::new();
        for gr in &layer {
            for h in one_step_expansions(gr) {
                if h.codim() == gr.codim() + 1 && seen.insert(h.clone()) {
                    next.push(h);
                }
            }
        }
        layer = next;
        // seen 需要跨层累积（codim 严格递增，此处天然不重）
    }
    layer
}

/// 余维数 ≤ codim 的全部连通稳定图。
pub fn graphs_upto(g: u32, n: u32, codim: u32) -> Vec<StableGraph> {
    (0..=codim).flat_map(|c| graphs(g, n, c)).collect()
}

/// 对图 Γ 做一步"余维 +1"扩张：在某顶点加自环（genus−1），或分裂某顶点。
fn one_step_expansions(gr: &StableGraph) -> Vec<StableGraph> {
    let mut out = Vec::new();
    let vs = gr.vertices();
    let mut next_leg = vs
        .iter()
        .map(|v| v.legs.iter().max().copied().unwrap_or(0))
        .max()
        .unwrap_or(gr.num_markings());
    let mut new_leg = || {
        next_leg += 1;
        next_leg
    };
    // 记录原始顶点表与边表（每次扩张从原图重建）
    for (vi, v) in vs.iter().enumerate() {
        // (a) 自环：顶点亏格 ≥ 1
        if v.genus >= 1 {
            let (l1, l2) = (new_leg(), new_leg());
            let mut genera: Vec<u32> = vs.iter().map(|x| x.genus).collect();
            let mut legs: Vec<Vec<Leg>> = vs.iter().map(|x| x.legs.clone()).collect();
            genera[vi] -= 1;
            legs[vi].push(l1);
            legs[vi].push(l2);
            let mut edges = gr.edges().to_vec();
            edges.push((l1, l2));
            if let Ok(h) = StableGraph::new(&genera, &legs, &edges) {
                out.push(h);
            }
        }
        // (b) 分裂顶点 vi 为两个顶点，中间加一条边
        // 顶点内容（亏格 g_v，腿集 L）→ (g1, L1) + (g2, L2)，g1+g2 = g_v；
        // L1 可为空（如 M̄_2 的 Δ_1：两个 g1 顶点无边外腿）。(L1, g1) 与
        // (L2, g2) 是同一个图，故只取含首条腿的 L1（无腿时只取 g1 ≤ g2）；
        // 两个新顶点须稳定（2g + 腿数 ≥ 3，腿数含新结点半边）。其余重复由
        // 规范形 seen 去重。
        let legs_v = v.legs.clone();
        let m = legs_v.len();
        for mask in 0u32..(1u32 << m) {
            if m > 0 && mask & 1 == 0 {
                continue;
            }
            let k1 = mask.count_ones();
            let k2 = m as u32 - k1;
            let l1: Vec<Leg> = legs_v
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(_, &l)| l)
                .collect();
            let l2: Vec<Leg> = legs_v
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) == 0)
                .map(|(_, &l)| l)
                .collect();
            for g1 in 0..=v.genus {
                let g2 = v.genus - g1;
                if (m == 0 && g1 > g2) || 2 * g1 + k1 < 2 || 2 * g2 + k2 < 2 {
                    continue;
                }
                let (e1, e2) = (new_leg(), new_leg());
                let mut genera: Vec<u32> = vs.iter().map(|x| x.genus).collect();
                let mut legs: Vec<Vec<Leg>> = vs.iter().map(|x| x.legs.clone()).collect();
                genera[vi] = g1;
                legs[vi] = l1.clone();
                legs[vi].push(e1);
                genera.insert(vi + 1, g2);
                legs.insert(vi + 1, {
                    let mut l = l2.clone();
                    l.push(e2);
                    l
                });
                let mut edges = gr.edges().to_vec();
                edges.push((e1, e2));
                if let Ok(h) = StableGraph::new(&genera, &legs, &edges) {
                    out.push(h);
                }
            }
        }
    }
    out
}

/// R^d(M̄_{g,n}) 的生成元列表（图 + 权重 d−codim(图) 的 κ-ψ 单项式）。
pub fn tautgens(g: u32, n: u32, d: u32) -> Vec<DecStratum> {
    let mut out = Vec::new();
    for gr in graphs_upto(g, n, d) {
        let w = d - gr.codim();
        for mono in monomials(&gr, w) {
            out.push(DecStratum {
                graph: gr.clone(),
                poly: mono,
            });
        }
    }
    out
}

/// 图 gr 上总权重 w 的全部 κ-ψ 单项式（规范形单项）。
fn monomials(gr: &StableGraph, w: u32) -> Vec<KappaPsiPolynomial> {
    // 逐顶点 κ 指标候选 + 腿候选，做权重恰好 w 的组合枚举
    let vs = gr.vertices();
    let nv = vs.len();
    let kappa_bounds: Vec<u32> = vs
        .iter()
        .map(|v| (3 * v.genus as i64 - 3 + v.legs.len() as i64).max(0) as u32)
        .collect();
    let mut all_legs: Vec<Leg> = vs.iter().flat_map(|v| v.legs.iter().copied()).collect();
    all_legs.sort_unstable();

    let mut out = Vec::new();
    let mut kappa_part: Vec<Vec<u32>> = vec![vec![]; nv];
    let mut legs_part: Vec<(Leg, u32)> = vec![];
    enum_kappa(&kappa_bounds, 0, w, &mut kappa_part, &mut |kappa, rem| {
        enum_psi(&all_legs, 0, rem, &mut legs_part, &mut |psi| {
            let term = KpTerm {
                kappa: kappa.to_vec(),
                psi: psi.to_vec(),
                coeff: Q::ONE,
            };
            let poly = KappaPsiPolynomial::new(vec![term]);
            if poly.weight() == Some(w) {
                out.push(poly);
            }
        });
    });
    out
}

/// 枚举逐顶点 κ 分划，总 κ 权重 ≤ budget（剩余给 ψ）。
fn enum_kappa(
    bounds: &[u32],
    v: usize,
    budget: u32,
    acc: &mut Vec<Vec<u32>>,
    k: &mut impl FnMut(&[Vec<u32>], u32),
) {
    if v == bounds.len() {
        k(acc, budget);
        return;
    }
    // 顶点 v 的分划：递归地选升序指标序列，权重 ≤ budget
    let mut part: Vec<u32> = vec![];
    enum_partition(bounds[v], budget, 1, &mut part, &mut |part, rem| {
        acc[v] = part.to_vec();
        enum_kappa(bounds, v + 1, rem, acc, k);
    });
    acc[v] = vec![];
}

/// 枚举 ≤ bound 的升序分划，总重 ≤ budget（回调带剩余额度）。
fn enum_partition(
    bound: u32,
    budget: u32,
    min_idx: u32,
    acc: &mut Vec<u32>,
    k: &mut impl FnMut(&[u32], u32),
) {
    k(acc, budget); // 空分划（可以随时停）
    let mut idx = min_idx;
    while idx <= bound && idx <= budget {
        acc.push(idx);
        enum_partition(bound, budget - idx, idx, acc, k);
        acc.pop();
        idx += 1;
    }
}

/// 枚举腿上的 ψ 幂分配，总幂恰好 rem。
fn enum_psi(
    legs: &[Leg],
    i: usize,
    rem: u32,
    acc: &mut Vec<(Leg, u32)>,
    k: &mut impl FnMut(&[(Leg, u32)]),
) {
    if rem == 0 {
        k(acc);
        return;
    }
    if i == legs.len() {
        return;
    }
    // 腿 i 分到 1..rem 幂，或 0（跳过）
    enum_psi(legs, i + 1, rem, acc, k);
    for p in 1..=rem {
        acc.push((legs[i], p));
        enum_psi(legs, i + 1, rem - p, acc, k);
        acc.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn count(g: u32, n: u32, d: u32) -> usize {
        tautgens(g, n, d).len()
    }

    #[test]
    fn graphs_g2_codim1() {
        // M̄_2 余维 1：Δ_0（g1 自环）、Δ_1（两个 g1 顶点）→ 2 个图
        let gs = graphs(2, 0, 1);
        assert_eq!(gs.len(), 2, "graphs: {gs:?}");
    }

    #[test]
    fn graphs_g1_1_codim1() {
        // M̄_{1,1} 余维 1：Δ_0（g0 顶点自环）
        assert_eq!(graphs(1, 1, 1).len(), 1);
    }

    #[test]
    fn graphs_g0_4_codim1() {
        // M̄_{0,4} 余维 1：3 个分离边界（Δ_{12},Δ_{13},Δ_{14}）
        assert_eq!(graphs(0, 4, 1).len(), 3);
    }

    #[test]
    fn tautgens_codim0_is_fundclass() {
        // d=0：只有基本类 1 个生成元
        assert_eq!(count(2, 3, 0), 1);
        assert_eq!(count(0, 4, 0), 1);
    }

    #[test]
    fn tautgens_span_all_graphs() {
        // 每个余维 ≤ d 的图至少有一个装饰生成元（κ 空、ψ 至少配满）
        for &(g, n, d) in &[(2u32, 0u32, 1u32), (1, 1, 1), (0, 4, 1), (2, 1, 2)] {
            let gens = tautgens(g, n, d);
            assert!(!gens.is_empty());
            for gr in graphs_upto(g, n, d) {
                assert!(
                    gens.iter().any(|t| t.graph == gr),
                    "图 {gr} 无装饰生成元（g={g},n={n},d={d}）"
                );
            }
        }
    }

    #[test]
    fn nodal_psi_generators_present() {
        // R^2(M̄_{1,1}) 中余维 1 的 Δ_0 图（g0 顶点 3 腿：标记 1 + 自环两半边）
        // 带权重 1 装饰：κ 界 0−3+3=0 → 无 κ；ψ_1、ψ_L1、ψ_L2 共 3 个
        let gs = graphs(1, 1, 1);
        assert_eq!(gs.len(), 1);
        let gens: Vec<_> = tautgens(1, 1, 2)
            .into_iter()
            .filter(|t| t.graph == gs[0])
            .collect();
        assert_eq!(gens.len(), 3);
    }

    #[test]
    fn codim_d_graph_gets_unit_decoration() {
        // R^1(M̄_{1,1})：Δ_0 图只带权重 0 装饰 [Δ_0]·1（ψ 装饰属 R^2）
        let gs = graphs(1, 1, 1);
        let gens: Vec<_> = tautgens(1, 1, 1)
            .into_iter()
            .filter(|t| t.graph == gs[0])
            .collect();
        assert_eq!(gens.len(), 1);
        assert_eq!(gens[0].poly.weight(), Some(0));
    }

    #[test]
    fn fundclass_kappa_generators() {
        // (2,0,d=2)：基本类上的权重 2 装饰 = κ_2、κ_1κ_1（无边外腿，无 ψ）
        let trivial: Vec<_> = tautgens(2, 0, 2)
            .into_iter()
            .filter(|t| t.graph.vertices().len() == 1 && t.graph.edges().is_empty())
            .collect();
        assert_eq!(trivial.len(), 2);
        // (2,0,d=1)：只有 κ_1（κ_2 权重已是 2）
        let trivial: Vec<_> = tautgens(2, 0, 1)
            .into_iter()
            .filter(|t| t.graph.vertices().len() == 1 && t.graph.edges().is_empty())
            .collect();
        assert_eq!(trivial.len(), 1);
    }

    #[test]
    fn weights_are_consistent() {
        // 所有生成元余维 = 图余维 + 多项式权重 = d
        for &(g, n, d) in &[(2u32, 0u32, 2u32), (1, 2, 2), (0, 5, 2)] {
            for t in tautgens(g, n, d) {
                assert_eq!(
                    t.graph.codim() + t.poly.weight().unwrap(),
                    d,
                    "g={g} n={n} d={d}"
                );
            }
        }
    }
}
