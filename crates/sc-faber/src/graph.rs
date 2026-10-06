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
        let Some(markings) = infer_markings(legs, edges) else {
            return Err("标记编号不连续（混入未配对结点腿？）".into());
        };

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
        // 同构规范化从颜色与边对整体重建，与初始编号无关
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
        let Some(markings) = infer_markings(legs, edges) else {
            return Err("标记编号不连续".into());
        };
        g.markings = markings;
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
    ///
    /// 顶点置换 p 保持边重数矩阵 m 时，结点腿双射的个数与 p 无关：u≠v 之间
    /// 的 m_uv 条边整体置换（定向由端点决定），u 处 m_uu 个自环置换且各自可
    /// 翻转，共 Π_{u<v} m_uv! · Π_u m_uu!·2^{m_uu}。顶点置换只在细化颜色类内
    /// 回溯（自同构保持细化颜色），逐顶点检查与已指派顶点的边重数相容。
    pub fn automorphism_number(&self) -> u64 {
        let n = self.vertices.len();
        let owner = self.leg_owner();
        let mut mult = vec![vec![0u32; n]; n];
        for &(a, b) in &self.edges {
            let (u, v) = (owner[a as usize], owner[b as usize]);
            mult[u][v] += 1;
            if u != v {
                mult[v][u] += 1;
            }
        }
        let fact = |m: u32| (1..=m as u64).product::<u64>();
        let mut edge_maps: u64 = 1;
        for (u, row) in mult.iter().enumerate() {
            edge_maps *= fact(row[u]) << row[u];
            for &m in &row[u + 1..] {
                edge_maps *= fact(m);
            }
        }
        let color = self.refined_colors(&mult);
        fn count(
            v: usize,
            img: &mut Vec<usize>,
            used: &mut Vec<bool>,
            color: &[usize],
            mult: &[Vec<u32>],
        ) -> u64 {
            let n = color.len();
            if v == n {
                return 1;
            }
            let mut total = 0;
            for w in 0..n {
                if used[w]
                    || color[w] != color[v]
                    || mult[w][w] != mult[v][v]
                    || (0..v).any(|x| mult[v][x] != mult[w][img[x]])
                {
                    continue;
                }
                used[w] = true;
                img[v] = w;
                total += count(v + 1, img, used, color, mult);
                used[w] = false;
            }
            total
        }
        count(0, &mut vec![0; n], &mut vec![false; n], &color, &mult) * edge_maps
    }

    /// 顶点颜色的 1 维 Weisfeiler–Leman 细化：初始色 (亏格, 标记集)，逐轮
    /// 加入（自环数, 邻居颜色及边重数的多重集），直到颜色数不再增加。
    /// 颜色号按签名排序给出，只依赖图结构（同构不变）。
    fn refined_colors(&self, mult: &[Vec<u32>]) -> Vec<usize> {
        let n = self.vertices.len();
        let rank = |sigs: &[Vec<u64>]| -> Vec<usize> {
            let mut distinct: Vec<&Vec<u64>> = sigs.iter().collect();
            distinct.sort();
            distinct.dedup();
            sigs.iter()
                .map(|s| distinct.binary_search(&s).unwrap())
                .collect()
        };
        let init: Vec<Vec<u64>> = self
            .vertices
            .iter()
            .map(|v| {
                let mut s = vec![v.genus as u64];
                s.extend(
                    v.legs
                        .iter()
                        .filter(|&&l| l <= self.markings)
                        .map(|&l| l as u64),
                );
                s
            })
            .collect();
        let mut color = rank(&init);
        let mut classes = color.iter().max().map_or(0, |&c| c + 1);
        loop {
            let sigs: Vec<Vec<u64>> = (0..n)
                .map(|v| {
                    let mut nb: Vec<(u64, u64)> = (0..n)
                        .filter(|&w| w != v && mult[v][w] > 0)
                        .map(|w| (color[w] as u64, mult[v][w] as u64))
                        .collect();
                    nb.sort_unstable();
                    let mut s = vec![color[v] as u64, mult[v][v] as u64];
                    s.extend(nb.into_iter().flat_map(|(c, m)| [c, m]));
                    s
                })
                .collect();
            let next = rank(&sigs);
            let k = next.iter().max().map_or(0, |&c| c + 1);
            color = next;
            if k == classes {
                return color;
            }
            classes = k;
        }
    }

    /// 腿 → 所在顶点下标（以腿号直接索引；不在任何顶点上为 `usize::MAX`）。
    fn leg_owner(&self) -> Vec<usize> {
        let max_leg = self
            .vertices
            .iter()
            .flat_map(|v| v.legs.iter().copied())
            .chain(self.edges.iter().flat_map(|&(a, b)| [a, b]))
            .max()
            .unwrap_or(0) as usize;
        let mut owner = vec![usize::MAX; max_leg + 1];
        for (i, v) in self.vertices.iter().enumerate() {
            for &l in &v.legs {
                owner[l as usize] = i;
            }
        }
        owner
    }

    /// 规范化：结点腿按「顶点规范序内首现顺序」重编号，再重排顶点与边。
    /// （改写前 new 的第一步；现仅测试参考实现使用。）
    #[cfg(test)]
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
    ///
    /// 结果只取决于最小边对表（同组顶点同色，故各位序的亏格与标记与安排
    /// 无关），与初始编号及枚举顺序无关。颜色与腿 → 顶点表只算一次，每个
    /// 安排用逆位置数组 O(E log E) 得到边对。
    ///
    /// 单顶点也走此路径：自环按 (n+1,n+2), (n+3,n+4), … 编号。原先单顶点
    /// 只做腿重编号，多个自环的配对随输入编号而变（同构图编码不同）。
    fn canonicalize_isomorphic(&mut self) {
        let n = self.vertices.len();
        if n == 0 {
            return;
        }
        let mk = self.markings;
        let colors: Vec<(u32, Vec<Leg>)> = self
            .vertices
            .iter()
            .map(|v| {
                (
                    v.genus,
                    v.legs.iter().copied().filter(|&l| l <= mk).collect(),
                )
            })
            .collect();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| colors[a].cmp(&colors[b]));
        // group_end[p] = 位序 p 所在同色组的结束位序（不含）
        let mut group_end = vec![0usize; n];
        let mut i = 0;
        while i < n {
            let mut j = i + 1;
            while j < n && colors[order[j]] == colors[order[i]] {
                j += 1;
            }
            group_end[i..j].fill(j);
            i = j;
        }
        let owner = self.leg_owner();
        let mut orig_pairs: Vec<(usize, usize)> = Vec::with_capacity(self.edges.len());
        for &(a, b) in &self.edges {
            let (u, v) = (owner[a as usize], owner[b as usize]);
            if u == usize::MAX || v == usize::MAX {
                return; // 边引用不存在的腿：留给 validate 报错
            }
            orig_pairs.push((u, v));
        }
        // arr[p] = 位于位序 p 的原顶点；pos[v] = 原顶点 v 的位序
        let mut arr = order.clone();
        let mut pos = vec![0usize; n];
        for (p, &v) in arr.iter().enumerate() {
            pos[v] = p;
        }
        fn search(
            k: usize,
            group_end: &[usize],
            arr: &mut Vec<usize>,
            pos: &mut Vec<usize>,
            orig_pairs: &[(usize, usize)],
            cur: &mut Vec<(usize, usize)>,
            best: &mut Vec<(usize, usize)>,
            found: &mut bool,
        ) {
            if k == arr.len() {
                cur.clear();
                cur.extend(orig_pairs.iter().map(|&(a, b)| {
                    let (pa, pb) = (pos[a], pos[b]);
                    (pa.min(pb), pa.max(pb))
                }));
                cur.sort_unstable();
                if !*found || *cur < *best {
                    best.clone_from(cur);
                    *found = true;
                }
                return;
            }
            // 位序 k 取组内尚未安排的任一顶点（交换法枚举组内全排列）
            for i in k..group_end[k] {
                arr.swap(k, i);
                pos[arr[k]] = k;
                pos[arr[i]] = i;
                search(k + 1, group_end, arr, pos, orig_pairs, cur, best, found);
                arr.swap(k, i);
                pos[arr[k]] = k;
                pos[arr[i]] = i;
            }
        }
        let mut cur = Vec::with_capacity(orig_pairs.len());
        let mut pairs = Vec::with_capacity(orig_pairs.len());
        let mut found = false;
        search(
            0,
            &group_end,
            &mut arr,
            &mut pos,
            &orig_pairs,
            &mut cur,
            &mut pairs,
            &mut found,
        );
        // 确定性重建：按排序后的最终边对，逐边分配两个新号
        let mut new_legs: Vec<Vec<Leg>> = order.iter().map(|&v| colors[v].1.clone()).collect();
        let mut next = mk + 1;
        let mut edges_new: Vec<(Leg, Leg)> = Vec::with_capacity(pairs.len());
        for &(i, j) in &pairs {
            new_legs[i].push(next);
            new_legs[j].push(next + 1);
            edges_new.push((next, next + 1));
            next += 2;
        }
        let mut vertices: Vec<Vertex> = order
            .iter()
            .zip(new_legs)
            .map(|(&v, mut legs)| {
                legs.sort_unstable();
                Vertex {
                    genus: colors[v].0,
                    legs,
                }
            })
            .collect();
        vertices.sort_unstable();
        self.vertices = vertices;
        self.edges = edges_new;
    }

    /// 参考实现（改写前的原版，测试对照用）。同构不变规范形（结构化证书）：顶点按色 (genus, 标记集) 排序，
    /// 同色组内枚举置换，取使边对多重集（(i,i)=自环、(i,j)=分离边，按位序）
    /// 最小的安排；随后按获奖安排确定性重建（结点半边按边序成对编号）。
    /// 修复 M5 定位的同构变体重复计数 bug。
    #[cfg(test)]
    fn canonicalize_isomorphic_reference(&mut self) {
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
        // 腿号很小：以腿号直接索引（owner = 腿所在顶点，usize::MAX = 未出现）
        let max_leg = self
            .vertices
            .iter()
            .flat_map(|v| v.legs.iter().copied())
            .chain(self.edges.iter().flat_map(|&(a, b)| [a, b]))
            .chain(std::iter::once(self.markings))
            .max()
            .unwrap_or(0) as usize;
        let mut owner = vec![usize::MAX; max_leg + 1];
        for (vi, v) in self.vertices.iter().enumerate() {
            for w in v.legs.windows(2) {
                if w[0] == w[1] {
                    return Err(format!("顶点 {vi} 腿 {} 重复", w[0]));
                }
            }
            for &l in &v.legs {
                if owner[l as usize] != usize::MAX {
                    return Err(format!("腿 {l} 出现在多个顶点"));
                }
                owner[l as usize] = vi;
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
            if owner[m as usize] == usize::MAX {
                return Err(format!("标记 {m} 缺失"));
            }
        }
        let mut used = vec![false; max_leg + 1];
        for &(a, b) in &self.edges {
            if a <= self.markings || b <= self.markings {
                return Err(format!("边 ({a},{b}) 引用了标记腿"));
            }
            if owner[a as usize] == usize::MAX || owner[b as usize] == usize::MAX {
                return Err(format!("边 ({a},{b}) 引用不存在的腿"));
            }
            for l in [a, b] {
                if used[l as usize] {
                    return Err(format!("边半边 {a}/{b} 被重复使用"));
                }
                used[l as usize] = true;
            }
        }
        for (l, &o) in owner.iter().enumerate() {
            if o != usize::MAX && l > self.markings as usize && !used[l] {
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
        for &(a, b) in &self.edges {
            let (i, j) = (owner[a as usize], owner[b as usize]);
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            parent[ri] = rj;
        }
        let mut roots: Vec<usize> = (0..self.vertices.len())
            .map(|i| find(&mut parent, i))
            .collect();
        roots.sort_unstable();
        roots.dedup();
        if roots.len() != 1 {
            return Err(format!("图不连通：{} 个连通分量", roots.len()));
        }
        Ok(())
    }
}

/// 标记数推导：标记集 = 出现在顶点但不出现在任何边中的腿，须恰为 {1..=n}；
/// 否则返回 None。
fn infer_markings(legs: &[Vec<Leg>], edges: &[(Leg, Leg)]) -> Option<u32> {
    let max_leg = legs
        .iter()
        .flatten()
        .copied()
        .chain(edges.iter().flat_map(|&(a, b)| [a, b]))
        .max()
        .unwrap_or(0) as usize;
    let mut in_edge = vec![false; max_leg + 1];
    for &(a, b) in edges {
        in_edge[a as usize] = true;
        in_edge[b as usize] = true;
    }
    let mut is_mark = vec![false; max_leg + 1];
    let mut count = 0usize;
    for &l in legs.iter().flatten() {
        if !in_edge[l as usize] && !is_mark[l as usize] {
            is_mark[l as usize] = true;
            count += 1;
        }
    }
    // 标记互异且共 count 个：恰为 {1..=count} ⟺ 每个标记都落在 1..=count
    is_mark
        .iter()
        .enumerate()
        .all(|(l, &m)| !m || (1..=count).contains(&l))
        .then_some(count as u32)
}

/// 枚举分组排列的笛卡尔积：groups[k] 的元素在位序上占据第 k 段。
#[cfg(test)]
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
#[cfg(test)]
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
#[cfg(test)]
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

    /// 改写前的自同构数：全部 n! 顶点置换按色过滤，再回溯数保边双射。
    fn automorphism_number_reference(g: &StableGraph) -> u64 {
        let n = g.vertices.len();
        let mark = |v: &Vertex| -> Vec<Leg> {
            v.legs
                .iter()
                .copied()
                .filter(|&l| l <= g.markings)
                .collect()
        };
        let mut good_perms: Vec<Vec<usize>> = Vec::new();
        let mut perm: Vec<usize> = (0..n).collect();
        permute(
            0,
            &mut perm,
            &mut vec![false; n],
            n,
            &mut |p: &[usize]| {
                let ok = p.iter().enumerate().all(|(i, &j)| {
                    g.vertices[i].genus == g.vertices[j].genus
                        && mark(&g.vertices[i]) == mark(&g.vertices[j])
                });
                if ok {
                    good_perms.push(p.to_vec());
                }
            },
        );
        let leg_vertex =
            |l: Leg| -> usize { g.vertices.iter().position(|v| v.legs.contains(&l)).unwrap() };
        good_perms
            .iter()
            .map(|p| {
                let mut used = vec![false; g.edges.len()];
                count_edge_maps(0, &mut used, &g.edges, &leg_vertex, p)
            })
            .sum()
    }

    /// 打乱顶点顺序、把结点腿改成稀疏的新编号、随机翻转边方向。
    fn relabel(g: &StableGraph, seed: &mut u64) -> (Vec<u32>, Vec<Vec<Leg>>, Vec<(Leg, Leg)>) {
        let mut rnd = |k: usize| -> usize {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) as usize) % k.max(1)
        };
        let n = g.vertices.len();
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            order.swap(i, rnd(i + 1));
        }
        let node_legs: Vec<Leg> = g.edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        let mut fresh: Vec<Leg> = (0..node_legs.len() as u32)
            .map(|i| g.markings + 1 + 3 * i)
            .collect();
        for i in (1..fresh.len()).rev() {
            fresh.swap(i, rnd(i + 1));
        }
        let map = |l: Leg| -> Leg {
            node_legs
                .iter()
                .position(|&x| x == l)
                .map_or(l, |i| fresh[i])
        };
        let genera = order.iter().map(|&v| g.vertices[v].genus).collect();
        let legs = order
            .iter()
            .map(|&v| {
                let mut ls: Vec<Leg> = g.vertices[v].legs.iter().map(|&l| map(l)).collect();
                ls.reverse();
                ls
            })
            .collect();
        let edges = g
            .edges
            .iter()
            .map(|&(a, b)| {
                if rnd(2) == 0 {
                    (map(a), map(b))
                } else {
                    (map(b), map(a))
                }
            })
            .collect();
        (genera, legs, edges)
    }

    /// 回归：单顶点多自环的编码与腿的配对方式无关（原先 (3,5),(4,6) 与
    /// (3,4),(5,6) 被当作两个图）。
    #[test]
    fn single_vertex_loops_canonical() {
        let a = StableGraph::new(&[0], &[vec![1, 2, 3, 4, 5, 6]], &[(3, 5), (4, 6)]).unwrap();
        let b = StableGraph::new(&[0], &[vec![1, 2, 3, 4, 5, 6]], &[(3, 4), (5, 6)]).unwrap();
        let c = StableGraph::new(&[0], &[vec![9, 1, 7, 2, 11, 8]], &[(11, 7), (9, 8)]).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, c);
        assert_eq!(a.edges(), &[(3, 4), (5, 6)]);
        assert_eq!(a.automorphism_number(), 8); // 2!·2²
    }

    /// 改写后的同构规范化与自同构数，逐图对照改写前的参考实现。
    #[test]
    fn canonical_form_and_automorphisms_match_reference() {
        let mut seed = 0x9e3779b97f4a7c15u64;
        let mut checked = 0;
        for &(g, n, c) in &[
            (0u32, 6u32, 3u32),
            (1, 3, 3),
            (2, 2, 4),
            (3, 0, 4),
            (2, 3, 3),
            (4, 0, 4),
        ] {
            for gr in crate::generators::graphs_upto(g, n, c) {
                assert_eq!(
                    gr.automorphism_number(),
                    automorphism_number_reference(&gr),
                    "{gr}"
                );
                for _ in 0..3 {
                    let (genera, legs, edges) = relabel(&gr, &mut seed);
                    // 新路径：规范形与输入编号无关
                    assert_eq!(StableGraph::new(&genera, &legs, &edges).unwrap(), gr);
                    // 参考路径：原版两步规范化得到同一编码（原版单顶点多自环的
                    // 编码随输入而变，见 single_vertex_loops_canonical）
                    if gr.vertices.len() == 1 {
                        continue;
                    }
                    let mut raw = StableGraph::new_raw(&genera, &legs, &edges).unwrap();
                    raw.edges = raw
                        .edges
                        .iter()
                        .map(|&(a, b)| (a.min(b), a.max(b)))
                        .collect();
                    raw.canonicalize();
                    raw.canonicalize_isomorphic_reference();
                    assert_eq!(raw, gr);
                }
                checked += 1;
            }
        }
        assert!(checked > 500, "只检查了 {checked} 个图");
    }
}
