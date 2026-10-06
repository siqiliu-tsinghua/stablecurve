//! 平衡权重系统求和：DR（Pixton）、Chiodo、Witten r-旋三处共用的内核。
//!
//! 给定图（顶点 0..nv、有向边 e = (起点, 终点)）、各顶点给定值 given 与
//! 模数 m，求
//!
//!   Σ_{w: E → {0..m−1}，平衡} Π_e c_e(w_e)，
//!
//! 平衡 = 每个顶点 v：given[v] + Σ_{起点为 v} w_e − Σ_{终点为 v} w_e ≡ 0 (mod m)。
//!
//! 逐边枚举 m^{|E|} 个权重再检查平衡是指数级的。这里利用平衡条件的结构：
//! - 自环（起点 = 终点）对平衡无贡献，其因子独立求和 Σ_w c_e(w)；
//! - 取非自环边的一棵生成树，非树边自由枚举（m^{b₁'}，b₁' 为非自环边的
//!   圈秩，树形图只有一项）；
//! - 树边的权由平衡条件自叶向根逐条唯一确定（mod m）。
//!
//! 每条边的因子 c_e 由调用方预先按 w ∈ 0..m 制表（`table[e][w]`），昂贵的
//! 求值（Bernoulli 多项式、r-旋边因子等）因此每个 (边类型, w) 只做一次。

use sc_num::{One, Q, Zero, is_zero};

/// Σ_{平衡 w} Π_e table[e][w_e]。`edges[e] = (起点, 终点)`，顶点号 < `given.len()`，
/// 图须连通；`table[e].len() == m`。
pub(crate) fn balanced_flow_sum(
    edges: &[(usize, usize)],
    given: &[i64],
    m: i64,
    table: &[&[Q]],
) -> Q {
    let nv = given.len();
    debug_assert_eq!(edges.len(), table.len());
    // 自环：独立因子
    let mut loop_factor = Q::ONE;
    let mut nonloop: Vec<usize> = Vec::new();
    for (e, &(u, v)) in edges.iter().enumerate() {
        if u == v {
            let s = table[e].iter().fold(Q::ZERO, |acc, c| acc + c);
            if is_zero(&s) {
                return Q::ZERO;
            }
            loop_factor *= s;
        } else {
            nonloop.push(e);
        }
    }
    // 非自环边的 BFS 生成树：order 为访问序，parent_edge[y] 为 y 到父顶点的树边
    let mut parent_edge: Vec<Option<usize>> = vec![None; nv];
    let mut seen = vec![false; nv];
    let mut in_tree = vec![false; edges.len()];
    let mut order: Vec<usize> = vec![0];
    seen[0] = true;
    let mut k = 0;
    while k < order.len() {
        let x = order[k];
        k += 1;
        for &e in &nonloop {
            let (u, v) = edges[e];
            let y = if u == x {
                v
            } else if v == x {
                u
            } else {
                continue;
            };
            if !seen[y] {
                seen[y] = true;
                in_tree[e] = true;
                parent_edge[y] = Some(e);
                order.push(y);
            }
        }
    }
    debug_assert_eq!(order.len(), nv, "图须连通");
    let free: Vec<usize> = nonloop.iter().copied().filter(|&e| !in_tree[e]).collect();
    let base: Vec<i64> = given.iter().map(|&x| x.rem_euclid(m)).collect();

    let mut total = Q::ZERO;
    let mut fw = vec![0i64; free.len()];
    let mut state = vec![0i64; nv];
    'outer: loop {
        state.copy_from_slice(&base);
        let mut term = Q::ONE;
        let mut zero = false;
        for (j, &e) in free.iter().enumerate() {
            let c = &table[e][fw[j] as usize];
            if is_zero(c) {
                zero = true;
                break;
            }
            let (u, v) = edges[e];
            state[u] = (state[u] + fw[j]).rem_euclid(m);
            state[v] = (state[v] - fw[j]).rem_euclid(m);
            term *= c;
        }
        if !zero {
            // 自叶向根确定树边：处理 y 时其子树边与全部自由边都已计入
            for &y in order.iter().skip(1).rev() {
                let e = parent_edge[y].expect("非根顶点必有树边");
                let (u, v) = edges[e];
                // y 为起点（+w）或终点（−w），取 w 使 state[y] 归零
                let w = if u == y {
                    (-state[y]).rem_euclid(m)
                } else {
                    state[y]
                };
                let c = &table[e][w as usize];
                if is_zero(c) {
                    zero = true;
                    break;
                }
                state[u] = (state[u] + w).rem_euclid(m);
                state[v] = (state[v] - w).rem_euclid(m);
                term *= c;
            }
            // 根处的平衡即 Σ given ≡ 0 (mod m)
            if !zero && state[order[0]] == 0 {
                total += term;
            }
        }
        // 自由边权的里程表式递增
        let mut j = 0;
        while j < fw.len() {
            if fw[j] + 1 < m {
                fw[j] += 1;
                continue 'outer;
            }
            fw[j] = 0;
            j += 1;
        }
        break;
    }
    total * loop_factor
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 参考实现：逐边枚举全部 m^{|E|} 个权重再检查平衡。
    fn brute(edges: &[(usize, usize)], given: &[i64], m: i64, table: &[&[Q]]) -> Q {
        fn rec(
            i: usize,
            state: &mut [i64],
            edges: &[(usize, usize)],
            m: i64,
            table: &[&[Q]],
            acc: &Q,
        ) -> Q {
            if i == edges.len() {
                return if state.iter().all(|&s| s == 0) {
                    acc.clone()
                } else {
                    Q::ZERO
                };
            }
            let (u, v) = edges[i];
            let mut total = Q::ZERO;
            for w in 0..m {
                state[u] = (state[u] + w).rem_euclid(m);
                state[v] = (state[v] - w).rem_euclid(m);
                total += rec(
                    i + 1,
                    state,
                    edges,
                    m,
                    table,
                    &(acc * &table[i][w as usize]),
                );
                state[u] = (state[u] - w).rem_euclid(m);
                state[v] = (state[v] + w).rem_euclid(m);
            }
            total
        }
        let mut state: Vec<i64> = given.iter().map(|&x| x.rem_euclid(m)).collect();
        rec(0, &mut state, edges, m, table, &Q::ONE)
    }

    #[test]
    fn matches_brute_force() {
        // 树、环、重边、自环、多个圈；给定值含负数；表含零（模拟 spin 的奇偶过滤）
        let graphs: Vec<(Vec<(usize, usize)>, Vec<i64>)> = vec![
            (vec![], vec![0]),
            (vec![(0, 0)], vec![3]),
            (vec![(0, 1)], vec![2, -2]),
            (vec![(1, 0), (0, 1)], vec![1, -1]),
            (vec![(0, 1), (1, 2), (2, 0)], vec![4, -1, -3]),
            (vec![(0, 1), (0, 1), (1, 1), (1, 2)], vec![1, 2, -3]),
            (
                vec![(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (2, 2)],
                vec![5, -2, 0, -3],
            ),
            (vec![(0, 1), (0, 1)], vec![1, 0]), // Σ given ≢ 0 → 0
        ];
        for m in [1i64, 2, 3, 5] {
            for (edges, given) in &graphs {
                let table: Vec<Vec<Q>> = (0..edges.len())
                    .map(|e| {
                        (0..m)
                            .map(|w| {
                                if (w + e as i64) % 3 == 2 {
                                    Q::ZERO
                                } else {
                                    Q::from((w + 1) * (e as i64 + 2)) / Q::from(e as i64 + 1)
                                }
                            })
                            .collect()
                    })
                    .collect();
                let table: Vec<&[Q]> = table.iter().map(|t| t.as_slice()).collect();
                assert_eq!(
                    balanced_flow_sum(edges, given, m, &table),
                    brute(edges, given, m, &table),
                    "m={m} edges={edges:?} given={given:?}"
                );
            }
        }
    }
}
