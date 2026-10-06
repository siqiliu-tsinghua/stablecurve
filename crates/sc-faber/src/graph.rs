//! 稳定图：M̄_{g,n} 边界 strata 的组合骨架。
//!
//! 一个稳定图 =（顶点集，边集）：顶点对应稳定曲线的归一化分支（带亏格与
//! 腿多重集），边对应结点。腿号 `1..=n` 固定为标记点，`n+1..` 为结点半边。
//! 对应 admcycles 的 `StableGraph`。
//!
//! 构造时即规范形：顶点内腿升序、边以 (小, 大) 表出并排序、结点腿按确定性
//! 顺序重编号、顶点按键 (genus, legs) 排序。因此 [`PartialEq`] 即编码相等。
//! 注意这不是完整的图同构规范形（同构但初始编号不同的图可能编码不同），
//! 完整规范形在 M1 乘法化简时引入。

use std::collections::BTreeSet;
use std::fmt;

/// 半边（腿）编号：`1..=n` 为标记，`n+1..` 为结点腿。
pub type Leg = u32;

/// 顶点：归一化分支，`legs` 构造后保证升序且含该分支的全部标记与结点半边。
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Vertex {
    pub genus: u32,
    pub legs: Vec<Leg>,
}

/// 稳定图（连通）。
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StableGraph {
    vertices: Vec<Vertex>,
    edges: Vec<(Leg, Leg)>, // 每条边 (a, b) 保证 a < b；整体升序
    markings: u32,
}

impl StableGraph {
    /// 由（各顶点亏格, 各顶点腿表, 边表）构造并验证。
    ///
    /// 这是 admcycles `StableGraph(genera, legs, edges)` 的对应接口。
    /// 腿的初始编号只需自洽，构造后会被规范化；标记数 n 由
    /// 「不出现在任何边中的最长前缀 1,2,...」自动确定。
    pub fn new(genera: &[u32], legs: &[Vec<Leg>], edges: &[(Leg, Leg)]) -> Result<Self, String> {
        if genera.len() != legs.len() {
            return Err(format!(
                "顶点数不一致：{} 个亏格 vs {} 个腿表",
                genera.len(),
                legs.len()
            ));
        }
        for &(a, b) in edges {
            if a == b {
                return Err(format!("边 ({a},{b}) 自配对同一半边"));
            }
        }
        for vs in legs {
            for &l in vs {
                if l == 0 {
                    return Err("腿号从 1 开始".into());
                }
            }
        }
        // 标记集 = 出现在顶点但不出现在任何边中的腿；必须恰为 {1..=n}
        let edge_legs: BTreeSet<Leg> = edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        let mut mark_set: BTreeSet<Leg> = BTreeSet::new();
        for l in legs.iter().flat_map(|v| v.iter().copied()) {
            if !edge_legs.contains(&l) {
                mark_set.insert(l);
            }
        }
        let markings = mark_set.len() as u32;
        if (1..=markings).collect::<BTreeSet<Leg>>() != mark_set {
            return Err("标记编号不连续（混入未配对结点腿？）".into());
        }

        let mut g = StableGraph {
            vertices: genera
                .iter()
                .zip(legs)
                .map(|(&genus, ls)| {
                    let mut ls = ls.clone();
                    ls.sort_unstable();
                    Vertex { genus, legs: ls }
                })
                .collect(),
            edges: edges
                .iter()
                .map(|&(a, b)| if a < b { (a, b) } else { (b, a) })
                .collect(),
            markings,
        };
        g.canonicalize();
        g.canonicalize_isomorphic();
        g.validate()?;
        Ok(g)
    }

    /// 原样构造（仅验证，不规范化）——用于外部表（如 admcycles 导出）的
    /// 腿编号与其装饰配对：规范化会重编号结点腿，使按旧编号书写的 ψ 失配。
    pub fn new_raw(
        genera: &[u32],
        legs: &[Vec<Leg>],
        edges: &[(Leg, Leg)],
    ) -> Result<Self, String> {
        let mut g = StableGraph {
            vertices: genera
                .iter()
                .zip(legs)
                .map(|(&genus, ls)| {
                    let mut ls = ls.clone();
                    ls.sort_unstable();
                    Vertex { genus, legs: ls }
                })
                .collect(),
            edges: edges.to_vec(),
            markings: 0,
        };
        // 标记集推导（与 new 相同）
        let edge_legs: BTreeSet<Leg> = edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        let mut mark_set: BTreeSet<Leg> = BTreeSet::new();
        for l in legs.iter().flat_map(|v| v.iter().copied()) {
            if !edge_legs.contains(&l) {
                mark_set.insert(l);
            }
        }
        g.markings = mark_set.len() as u32;
        if (1..=g.markings).collect::<BTreeSet<Leg>>() != mark_set {
            return Err("标记编号不连续".into());
        }
        g.validate()?;
        Ok(g)
    }

    /// 标记数 n。
    pub fn num_markings(&self) -> u32 {
        self.markings
    }

    /// 总亏格：Σ 顶点亏格 + 圈秩 b₁（连通时 b₁ = E − V + 1）。
    pub fn genus(&self) -> u32 {
        let v = self.vertices.len() as u32;
        let e = self.edges.len() as u32;
        self.vertices.iter().map(|x| x.genus).sum::<u32>() + e + 1 - v
    }

    /// 该图对应 stratum 的余维数（在 ambient M̄_{g,n} 中）。等于边数（结点数）。
    pub fn codim(&self) -> u32 {
        self.edges.len() as u32
    }

    /// stratum 自身维数：Σ_v (3g_v − 3 + n_v)。
    pub fn stratum_dim(&self) -> i64 {
        self.vertices
            .iter()
            .map(|x| 3 * x.genus as i64 - 3 + x.legs.len() as i64)
            .sum()
    }

    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }

    pub fn edges(&self) -> &[(Leg, Leg)] {
        &self.edges
    }

    /// 图自同构数（标记点逐点固定）：保边、保亏格的（顶点置换 × 结点腿
    /// 双射）个数。stratum 类约定 [Γ] = ξ_*fund/|Aut|。
    pub fn automorphism_number(&self) -> u64 {
        let n = self.vertices.len();
        let mut good_perms: Vec<Vec<usize>> = Vec::new();
        let mut perm: Vec<usize> = (0..n).collect();
        let mut used = vec![false; n];
        permute(0, &mut perm, &mut used, n, &mut |p: &[usize]| {
            for (i, &j) in p.iter().enumerate() {
                let a = &self.vertices[i];
                let b = &self.vertices[j];
                let ma: Vec<Leg> = a
                    .legs
                    .iter()
                    .copied()
                    .filter(|l| *l <= self.markings)
                    .collect();
                let mb: Vec<Leg> = b
                    .legs
                    .iter()
                    .copied()
                    .filter(|l| *l <= self.markings)
                    .collect();
                if a.genus != b.genus || ma != mb {
                    return;
                }
            }
            good_perms.push(p.to_vec());
        });
        let leg_vertex = |l: Leg| -> usize {
            self.vertices
                .iter()
                .position(|v| v.legs.contains(&l))
                .unwrap()
        };
        let mut total: u64 = 0;
        for p in &good_perms {
            let m = self.edges.len();
            let mut used_img = vec![false; m];
            total += count_edge_maps(0, &mut used_img, &self.edges, &|l| leg_vertex(l), p);
        }
        total
    }

    /// 规范化：结点腿按「顶点规范序内首现顺序」重编号，再重排顶点与边。
    fn canonicalize(&mut self) {
        let markings = self.markings;
        let mut order: Vec<usize> = (0..self.vertices.len()).collect();
        order.sort_by(|&i, &j| {
            let key = |k: usize| {
                let v = &self.vertices[k];
                (
                    v.genus,
                    v.legs
                        .iter()
                        .copied()
                        .filter(|&l| l <= markings)
                        .collect::<Vec<_>>(),
                    v.legs.len(),
                )
            };
            key(i).cmp(&key(j))
        });
        let mut next = markings + 1;
        let mut mapping = std::collections::BTreeMap::new();
        for &vi in &order {
            let legs: Vec<Leg> = self.vertices[vi].legs.clone();
            for l in legs {
                if l > markings && !mapping.contains_key(&l) {
                    mapping.insert(l, next);
                    next += 1;
                }
            }
        }
        let map = |l: Leg| -> Leg { *mapping.get(&l).unwrap_or(&l) };
        for v in &mut self.vertices {
            for l in &mut v.legs {
                *l = map(*l);
            }
            v.legs.sort_unstable();
        }
        for (a, b) in &mut self.edges {
            let (na, nb) = (map(*a), map(*b));
            *a = na.min(nb);
            *b = na.max(nb);
        }
        self.edges.sort_unstable();
        self.vertices.sort_unstable();
    }

    /// 同构不变规范形（结构化证书）：顶点按色 (genus, 标记集) 排序，
    /// 同色组内枚举置换，取使边对多重集（(i,i)=自环、(i,j)=分离边，按位序）
    /// 最小的安排；随后按获奖安排确定性重建（结点半边按边序成对编号）。
    /// 修复 M5 定位的同构变体重复计数 bug。
    fn canonicalize_isomorphic(&mut self) {
        let n = self.vertices.len();
        if n <= 1 {
            return;
        }
        // 顶点色与初始排序
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let key = |k: usize| {
                (
                    self.vertices[k].genus,
                    self.vertices[k]
                        .legs
                        .iter()
                        .copied()
                        .filter(|l| *l <= self.markings)
                        .collect::<Vec<Leg>>(),
                )
            };
            key(a).cmp(&key(b))
        });
        // 同色分组（在 order 上的连续区段）
        let color = |k: usize| {
            (
                self.vertices[k].genus,
                self.vertices[k]
                    .legs
                    .iter()
                    .copied()
                    .filter(|l| *l <= self.markings)
                    .collect::<Vec<Leg>>(),
            )
        };
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut i = 0;
        while i < order.len() {
            let mut g = vec![order[i]];
            let mut j = i + 1;
            while j < order.len() && color(order[j]) == color(order[i]) {
                g.push(order[j]);
                j += 1;
            }
            groups.push(g);
            i = j;
        }
        // 原顶点 → 其边表（位对）
        let vertex_of = |l: Leg| -> usize {
            self.vertices
                .iter()
                .position(|v| v.legs.contains(&l))
                .unwrap()
        };
        let orig_pairs: Vec<(usize, usize)> = self
            .edges
            .iter()
            .map(|&(a, b)| {
                let (va, vb) = (vertex_of(a), vertex_of(b));
                (va.min(vb), va.max(vb))
            })
            .collect();
        // 组内枚举：pos_of[原顶点] = 最终位序
        let mut best: Option<(Vec<(usize, usize)>, Vec<usize>)> = None; // (边对, 安排)
        let mut pos = vec![0usize; n];
        rec_arrange(&groups, 0, &mut pos, &mut |arr: &[usize]| {
            // arr: 长度 n 的数组，arr[p] = 位于位序 p 的原顶点
            let mut pairs: Vec<(usize, usize)> = orig_pairs
                .iter()
                .map(|&(a, b)| {
                    let pa = arr.iter().position(|&x| x == a).unwrap();
                    let pb = arr.iter().position(|&x| x == b).unwrap();
                    (pa.min(pb), pa.max(pb))
                })
                .collect();
            pairs.sort_unstable();
            if best.as_ref().map(|(b, _)| pairs < *b).unwrap_or(true) {
                best = Some((pairs, arr.to_vec()));
            }
        });
        let Some((pairs, arr)) = best else { return };
        // 确定性重建：按排序后的最终边对，逐边分配两个新号
        let genera: Vec<u32> = arr.iter().map(|&v| self.vertices[v].genus).collect();
        let mut new_legs: Vec<Vec<Leg>> = arr
            .iter()
            .map(|&v| {
                self.vertices[v]
                    .legs
                    .iter()
                    .copied()
                    .filter(|l| *l <= self.markings)
                    .collect()
            })
            .collect();
        let mut next = self.markings + 1;
        let mut edges_new: Vec<(Leg, Leg)> = Vec::new();
        for &(i, j) in &pairs {
            let a = next;
            let b = next + 1;
            next += 2;
            new_legs[i].push(a);
            new_legs[j].push(b);
            edges_new.push((a, b));
        }
        for l in &mut new_legs {
            l.sort_unstable();
        }
        edges_new.sort_unstable();
        self.vertices = (0..n)
            .map(|v| Vertex {
                genus: genera[v],
                legs: new_legs[v].clone(),
            })
            .collect();
        self.vertices.sort_unstable();
        self.edges = edges_new;
    }

    /// 结构验证：标记恰出现一次、结点腿恰配对一次、逐点稳定、连通。
    fn validate(&self) -> Result<(), String> {
        let mut seen: BTreeSet<Leg> = BTreeSet::new();
        for (vi, v) in self.vertices.iter().enumerate() {
            for w in v.legs.windows(2) {
                if w[0] == w[1] {
                    return Err(format!("顶点 {vi} 腿 {} 重复", w[0]));
                }
            }
            for &l in &v.legs {
                if !seen.insert(l) {
                    return Err(format!("腿 {l} 出现在多个顶点"));
                }
            }
            if 2 * v.genus as i64 - 2 + v.legs.len() as i64 <= 0 {
                return Err(format!(
                    "顶点 {vi} 不稳定：g={}，{} 条腿",
                    v.genus,
                    v.legs.len()
                ));
            }
        }
        for m in 1..=self.markings {
            if !seen.contains(&m) {
                return Err(format!("标记 {m} 缺失"));
            }
        }
        let mut used: BTreeSet<Leg> = BTreeSet::new();
        for &(a, b) in &self.edges {
            if a <= self.markings || b <= self.markings {
                return Err(format!("边 ({a},{b}) 引用了标记腿"));
            }
            if !seen.contains(&a) || !seen.contains(&b) {
                return Err(format!("边 ({a},{b}) 引用不存在的腿"));
            }
            if !used.insert(a) || !used.insert(b) {
                return Err(format!("边半边 {a}/{b} 被重复使用"));
            }
        }
        for &l in &seen {
            if l > self.markings && !used.contains(&l) {
                return Err(format!("结点腿 {l} 未配对"));
            }
        }
        let mut parent: Vec<usize> = (0..self.vertices.len()).collect();
        fn find(p: &mut Vec<usize>, i: usize) -> usize {
            if p[i] != i {
                let r = find(p, p[i]);
                p[i] = r;
            }
            p[i]
        }
        let leg2vtx: std::collections::BTreeMap<Leg, usize> = self
            .vertices
            .iter()
            .enumerate()
            .flat_map(|(i, v)| v.legs.iter().map(move |&l| (l, i)))
            .collect();
        for &(a, b) in &self.edges {
            let (i, j) = (leg2vtx[&a], leg2vtx[&b]);
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            parent[ri] = rj;
        }
        let roots: BTreeSet<usize> = (0..self.vertices.len())
            .map(|i| find(&mut parent, i))
            .collect();
        if roots.len() != 1 {
            return Err(format!("图不连通：{} 个连通分量", roots.len()));
        }
        Ok(())
    }
}

/// 枚举分组排列的笛卡尔积：groups[k] 的元素在位序上占据第 k 段。
fn rec_arrange(
    groups: &[Vec<usize>],
    gi: usize,
    pos: &mut Vec<usize>,
    f: &mut impl FnMut(&[usize]),
) {
    if gi == groups.len() {
        f(pos);
        return;
    }
    let g = &groups[gi];
    let mut perm = g.clone();
    permute(
        0,
        &mut perm,
        &mut vec![false; g.len()],
        g.len(),
        &mut |p: &[usize]| {
            let base: Vec<usize> = p.iter().map(|&x| g[x]).collect();
            let saved: Vec<usize> = pos.to_vec();
            let start = groups[..gi].iter().map(|gr| gr.len()).sum::<usize>();
            for (k, &v) in base.iter().enumerate() {
                pos[start + k] = v;
            }
            rec_arrange(groups, gi + 1, pos, f);
            for (k, v) in saved.iter().enumerate() {
                pos[k] = *v;
            }
        },
    );
}

/// 顶点全排列枚举（回调式）。
fn permute(
    k: usize,
    perm: &mut Vec<usize>,
    used: &mut Vec<bool>,
    n: usize,
    f: &mut impl FnMut(&[usize]),
) {
    if k == n {
        f(perm);
        return;
    }
    for i in 0..n {
        if !used[i] {
            used[i] = true;
            perm[k] = i;
            permute(k + 1, perm, used, n, f);
            used[i] = false;
        }
    }
}

/// 在顶点映射 p 下，把源边逐条映射到未用目标边（两种定向），计数保边双射。
fn count_edge_maps(
    i: usize,
    used_img: &mut Vec<bool>,
    edges: &[(Leg, Leg)],
    leg_vertex: &dyn Fn(Leg) -> usize,
    p: &[usize],
) -> u64 {
    if i == edges.len() {
        return 1;
    }
    let (a, b) = edges[i];
    let mut sum = 0;
    for (j, &(c, d)) in edges.iter().enumerate() {
        if used_img[j] {
            continue;
        }
        // 定向 1：a→c, b→d；定向 2：a→d, b→c
        if p[leg_vertex(a)] == leg_vertex(c) && p[leg_vertex(b)] == leg_vertex(d) {
            used_img[j] = true;
            sum += count_edge_maps(i + 1, used_img, edges, leg_vertex, p);
            used_img[j] = false;
        }
        if p[leg_vertex(a)] == leg_vertex(d) && p[leg_vertex(b)] == leg_vertex(c) {
            used_img[j] = true;
            sum += count_edge_maps(i + 1, used_img, edges, leg_vertex, p);
            used_img[j] = false;
        }
    }
    sum
}

impl fmt::Display for StableGraph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let vs: Vec<String> = self
            .vertices
            .iter()
            .map(|v| format!("{}{}", v.genus, v.legs.len()))
            .collect();
        write!(
            f,
            "[{}],{},[{}]",
            vs.join(","),
            (1..=self.markings)
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(""),
            self.edges
                .iter()
                .map(|(a, b)| format!("{a}{b}"))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trivial_graph() {
        // M̄_{2,3} 基本类：单顶点
        let g = StableGraph::new(&[2], &[vec![1, 2, 3]], &[]).unwrap();
        assert_eq!(g.genus(), 2);
        assert_eq!(g.num_markings(), 3);
        assert_eq!(g.codim(), 0);
        assert_eq!(g.stratum_dim(), 3 * 2 - 3 + 3);
    }

    #[test]
    fn irreducible_divisor() {
        // M̄_{1,2} 的 Δ_0：g=0 顶点带标记 1,2 与自环（0-3+4=1 维，余维 1）
        let g = StableGraph::new(&[0], &[vec![1, 2, 3, 4]], &[(3, 4)]).unwrap();
        assert_eq!(g.genus(), 1); // Σg=0，b1=1
        assert_eq!(g.codim(), 1);
        assert_eq!(g.num_markings(), 2);
        assert_eq!(g.stratum_dim(), 1);
    }

    #[test]
    fn irreducible_divisor_g2() {
        // M̄_{2,0} 的 Δ_0：g=1 顶点 + 自环
        let g = StableGraph::new(&[1], &[vec![1, 2]], &[(1, 2)]).unwrap();
        assert_eq!(g.genus(), 2);
        assert_eq!(g.num_markings(), 0);
        assert_eq!(g.codim(), 1);
        assert_eq!(g.stratum_dim(), 2);
    }

    #[test]
    fn separating_divisor() {
        // M̄_{2,2} 的分离边界：两个 g=1 顶点各带一个标记（Σg=2，b1=0）
        let g = StableGraph::new(&[1, 1], &[vec![1, 3], vec![2, 4]], &[(3, 4)]).unwrap();
        assert_eq!(g.genus(), 2);
        assert_eq!(g.codim(), 1);
        assert_eq!(g.stratum_dim(), 4); // (3-3+2)+(3-3+2)
    }

    #[test]
    fn canonicalization_is_input_order_independent() {
        // 同一图三种写法：顶点顺序、边方向、结点腿编号不同 → 规范形相等
        let a = StableGraph::new(&[1, 1], &[vec![1, 3], vec![2, 4]], &[(3, 4)]).unwrap();
        let b = StableGraph::new(&[1, 1], &[vec![3, 1], vec![4, 2]], &[(4, 3)]).unwrap();
        let c = StableGraph::new(&[1, 1], &[vec![4, 2], vec![1, 7]], &[(7, 4)]).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, c);
        // 引用标记腿的边必须被拒（标记集因此不连续，先命中该检查）
        let e = StableGraph::new(&[1, 1], &[vec![1, 3], vec![2, 4]], &[(1, 3)]).unwrap_err();
        assert!(e.contains("不连续"));
    }

    #[test]
    fn rejects_unstable() {
        // g=0 顶点只有 2 条腿（标记 1 + 结点半边）
        let e = StableGraph::new(&[0, 1], &[vec![1, 3], vec![2, 4]], &[(3, 4)]).unwrap_err();
        assert!(e.contains("不稳定"));
    }

    #[test]
    fn rejects_disconnected() {
        let e = StableGraph::new(&[1, 1], &[vec![1], vec![2]], &[]).unwrap_err();
        assert!(e.contains("连通"));
    }

    #[test]
    fn rejects_duplicate_marking() {
        let e = StableGraph::new(&[1, 1], &[vec![1, 3], vec![1, 4]], &[(3, 4)]).unwrap_err();
        assert!(e.contains("多个顶点"));
    }

    #[test]
    fn rejects_dangling_leg() {
        // 结点腿 5 出现在顶点但不在任何边 → 标记集 {1,2,5} 不连续
        let e = StableGraph::new(&[1, 0], &[vec![1, 3], vec![2, 4, 5]], &[(3, 4)]).unwrap_err();
        assert!(e.contains("不连续"));
    }
}
