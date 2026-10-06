//! 第一档显式类：Hain 除子（后续：Θ 类、Hyperell/Biell、Chiodo、Witten r-旋）。
//!
//! Hain 除子：Abel–Jacobi 像 θ 除子在紧型上的延拓（admcycles `Hain_divisor`）。
//! 公式：−¼ Σ_{l=0}^{g} Σ_{I⊆{1..n}} (Σ_{i∈I} A_i)² · δ_{l,I}，
//! 其中 δ_{l,I} 为广义边界除子（含 l=0/|I|≤1 与 l=g/|Iᶜ|≤1 的 ψ 退化情形）。
//!
//! 口径：admcycles 的 decstratum 求值不除 |Aut|，其 δ 系数含 1/aut；
//! 我们 evaluate 统一除 aut，故 δ 图项系数为 1。

use crate::graph::StableGraph;
use crate::kppoly::{KappaPsiPolynomial as KP, KpTerm};
use crate::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero, is_zero};

/// δ_{l,I}（广义边界除子，返回 TautClass 项表；空 = 0）。
fn delt(g: u32, a: &[i64], l: u32, subset: &[u32]) -> Vec<DecStratum> {
    let n = a.len() as u32;
    let full: Vec<u32> = (1..=n).collect();
    // ψ 退化情形
    if l == 0 && subset.len() == 1 {
        let i = subset[0];
        return vec![fund_psi(g, n, i, -1)];
    }
    if l == g && subset.len() == (n - 1) as usize {
        let comp: Vec<u32> = full
            .iter()
            .copied()
            .filter(|x| !subset.contains(x))
            .collect();
        return vec![fund_psi(g, n, comp[0], -1)];
    }
    if (l == 0 && subset.is_empty()) || (l == g && subset.len() == n as usize) {
        return vec![];
    }
    // 一般分离边界图（顶点 [l: I∪{n+1}], [g−l: Iᶜ∪{n+2}]，边 (n+1,n+2)）
    let comp: Vec<u32> = full
        .iter()
        .copied()
        .filter(|x| !subset.contains(x))
        .collect();
    let mut legs1 = subset.to_vec();
    legs1.push(n + 1);
    let mut legs2 = comp;
    legs2.push(n + 2);
    let graph = match StableGraph::new(&[l, g - l], &[legs1, legs2], &[(n + 1, n + 2)]) {
        Ok(gr) => gr,
        Err(_) => return vec![],
    };
    vec![DecStratum {
        graph,
        poly: KP::one(),
    }]
}

fn fund_psi(g: u32, n: u32, leg: u32, coeff: i64) -> DecStratum {
    let legs: Vec<u32> = (1..=n).collect();
    DecStratum {
        graph: StableGraph::new(&[g], &[legs], &[]).unwrap(),
        poly: KP::new(vec![KpTerm {
            kappa: vec![vec![]],
            psi: vec![(leg, 1)],
            coeff: Q::from(coeff),
        }]),
    }
}

/// Hain 除子 D_{g,A}（ΣA = 0）。
pub fn hain_divisor(g: u32, a: &[i64]) -> Result<TautClass, String> {
    let n = a.len() as u32;
    if a.iter().sum::<i64>() != 0 {
        return Err(format!("Hain 要求 ΣA = 0，得到 {a:?}"));
    }
    if 2 * g as i64 - 2 + n as i64 <= 0 {
        return Err(format!("M̄_({g},{n}) 不稳定"));
    }
    let mut terms: Vec<DecStratum> = Vec::new();
    let idx: Vec<u32> = (1..=n).collect();
    for l in 0..=g {
        for mask in 0u64..(1u64 << n) {
            let subset: Vec<u32> = idx
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(_, &x)| x)
                .collect();
            let s: i64 = subset.iter().map(|&i| a[(i - 1) as usize]).sum();
            for d in delt(g, a, l, &subset) {
                terms.push(DecStratum {
                    graph: d.graph,
                    poly: d.poly.scale(&(Q::from_signeds(-1, 4) * Q::from(s * s))),
                });
            }
        }
    }
    Ok(TautClass::new(terms))
}

/// 通用：∫_{M̄_{g,n}} cl · X，X 为基本类上的 κ-ψ 单项式
/// （κ_j 限制到覆盖顶点之和；ψ_i 落其顶点）。
/// ψ 的腿号须在标记 1..=n 内——n+1.. 是结点半边的编号，越界会被
/// 误当作结点 ψ。
pub fn int_against(
    cl: &TautClass,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<sc_num::Q, String> {
    use sc_num::Zero;
    if let Some(t) = cl.terms.first() {
        crate::kppoly::check_marking_legs(psi, t.graph.num_markings())?;
    }
    let mut total = sc_num::Q::ZERO;
    for t in &cl.terms {
        let graph = &t.graph;
        let nv = graph.vertices().len();
        for term in &t.poly.terms {
            let factors: Vec<u32> = kappa
                .iter()
                .flat_map(|&(j, p)| vec![j; p as usize])
                .collect();
            let ncomb = nv.pow(factors.len() as u32);
            for comb in 0..ncomb {
                // 原类项自身的 κ 逐顶点保留，X 的 κ 因子再分配到顶点
                let mut kverts: Vec<Vec<u32>> = term.kappa.clone();
                kverts.resize(nv, Vec::new());
                let mut c = comb;
                for &j in &factors {
                    kverts[c % nv].push(j);
                    c /= nv;
                }
                let mut xpsi = term.psi.clone();
                for &(l, p) in psi {
                    if let Some(e) = xpsi.iter_mut().find(|e| e.0 == l) {
                        e.1 += p;
                    } else {
                        xpsi.push((l, p));
                    }
                }
                xpsi.retain(|&(_, p)| p > 0);
                xpsi.sort_unstable();
                let poly = KP::new(vec![KpTerm {
                    kappa: kverts,
                    psi: xpsi,
                    coeff: term.coeff.clone(),
                }]);
                let sub = TautClass::new(vec![DecStratum {
                    graph: graph.clone(),
                    poly,
                }]);
                total += sub.evaluate()?;
            }
        }
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Θ 类初始条件独立可验的近亲：Hain(1,(1,−1)) = ψ₁+ψ₂+δ/2 …admcycles 差分见
    /// tools；先验已知：g=1 时 D = (1/2)(Σ A_i² ψ_i + …) 紧型上 = DR。
    #[test]
    fn hain_basic() {
        let d = hain_divisor(1, &[1, -1]).unwrap();
        // 展开系数非空且余维 1
        assert!(!d.is_zero());
        assert_eq!(d.codim(), Some(1));
    }
}

#[cfg(test)]
mod differential {
    use super::*;

    /// admcycles Hain·X 实测锚点（2026-10-05）。
    #[test]
    fn hain_vs_admcycles() {
        // g=1: ∫D·ψ₁ψ₂ = 0
        let d = hain_divisor(1, &[1, -1]).unwrap();
        assert_eq!(int_against(&d, &[], &[(1, 1), (2, 1)]).unwrap(), Q::ZERO);
        // g=2: ∫D·ψ₁²ψ₂ψ₃² = 9/160
        let d = hain_divisor(2, &[1, 1, -2]).unwrap();
        assert_eq!(
            int_against(&d, &[], &[(1, 2), (2, 1), (3, 2)])
                .unwrap()
                .to_string(),
            "9/160"
        );
        // g=2: ∫D·κ₁ψ₁ψ₂² = 33/160
        let d = hain_divisor(2, &[3, -3]).unwrap();
        assert_eq!(
            int_against(&d, &[(1, 1)], &[(1, 1), (2, 2)])
                .unwrap()
                .to_string(),
            "33/160"
        );
    }
}

// ===================== Chiodo DR（固定 r）与 Θ 类 =====================

use crate::generators::tautgens;
use crate::socle::bernoulli;

/// B_n(x)（Bernoulli 多项式求值）。
fn bern_poly(n: u32, x: &Q) -> Q {
    // B_n(x) = Σ_{k=0}^{n} C(n,k) B_k x^{n−k}
    let mut acc = Q::ZERO;
    for k in 0..=n {
        let mut binom = 1i64;
        for i in 0..k {
            // C(n,k) 逐步（n 小）
            binom = binom * (n - i) as i64 / (i + 1) as i64;
        }
        let mut term = bernoulli(k as u64) * Q::from(binom);
        for _ in k..n {
            term *= x;
        }
        acc += term;
    }
    acc
}

/// A-多项式：A 指数向量 → 系数。
type APoly = std::collections::BTreeMap<Vec<u32>, Q>;

fn apoly_eval(p: &APoly, a: &[Q]) -> Q {
    let mut acc = Q::ZERO;
    for (exp, c) in p {
        let mut t = c.clone();
        for (i, &e) in exp.iter().enumerate() {
            for _ in 0..e {
                t *= &a[i];
            }
        }
        acc += t;
    }
    acc
}

/// edge_power_series：h = (1−exp(Σ_i A_i(X^{i+1}−(−Y)^{i+1})))/(X+Y)
/// 的 (X^{e1}Y^{e2})-系数（A-多项式）。用 h = −Σ_{k≥1}(X+Y)^{k−1}S^k/k!，
/// S = Σ_i A_i Σ_{j=0}^{i} X^{i−j}(−Y)^j。
fn edge_coeff(e1: u32, e2: u32) -> APoly {
    let d = e1 + e2;
    // S 的单项式表：(xdeg, ydeg, A_i, sign)
    let mut s_mons: Vec<(u32, u32, u32, i64)> = Vec::new();
    for i in 0..=d {
        for j in 0..=i {
            let sign = if j % 2 == 0 { 1 } else { -1 };
            s_mons.push((i - j, j, i, sign));
        }
    }
    // acc: HashMap<(xdeg,ydeg), APoly>
    let mut acc: std::collections::BTreeMap<(u32, u32), APoly> = std::collections::BTreeMap::new();
    // 幂 (X+Y)^{k−1} S^k，截断到度 d；k = 1..=d+1
    for k in 1..=(d + 1) {
        // 先算 S^k（截断单项式度 ≤ d），再乘 (X+Y)^{k−1}
        // S^k：k 个 S 单项式的组合，X度+Y度 ≤ d
        let mut sk: Vec<(u32, u32, Vec<u32>, i64)> = vec![(0, 0, vec![0u32; d as usize + 1], 1)];
        for _ in 0..k {
            let mut next = Vec::new();
            for (x, y, mut exp, sg) in sk {
                for &(mx, my, mi, ms) in &s_mons {
                    let (nx, ny) = (x + mx, y + my);
                    if nx + ny > d {
                        continue;
                    }
                    exp[mi as usize] += 1;
                    next.push((nx, ny, exp.clone(), sg * ms));
                    exp[mi as usize] -= 1;
                }
            }
            // 合并同键
            let mut merged: std::collections::BTreeMap<(u32, u32, Vec<u32>), i64> =
                std::collections::BTreeMap::new();
            for (x, y, exp, sg) in next {
                *merged.entry((x, y, exp)).or_insert(0) += sg;
            }
            sk = merged
                .into_iter()
                .filter(|(_, v)| *v != 0)
                .map(|((x, y, e), v)| (x, y, e, v))
                .collect();
        }
        // 乘 (X+Y)^{k−1} / k!，带负号
        let mut fact = Q::ONE;
        for i in 2..=k {
            fact *= Q::from(i as i64);
        }
        let pref = Q::from(-1) / fact;
        for (x, y, exp, sg) in sk {
            for j in 0..k {
                // (X+Y)^{k−1} 的 X^j Y^{k−1−j} 项系数为 C(k−1, j)
                let mut binom = 1i64;
                for t in 0..j {
                    binom = binom * ((k - 1 - t) as i64) / (t + 1) as i64;
                }
                let (nx, ny) = (x + j, y + (k - 1 - j));
                if nx == e1 && ny == e2 {
                    // 额外乘 C(k−1,j)（在 pref 已含 −1/k!）
                    let mut key = vec![0u32; d as usize + 1];
                    for (i, &c) in exp.iter().enumerate() {
                        key[i] = c;
                    }
                    *acc.entry((nx, ny))
                        .or_default()
                        .entry(key)
                        .or_insert(Q::ZERO) += pref.clone() * Q::from(sg * binom);
                }
            }
        }
    }
    acc.get(&(e1, e2)).cloned().unwrap_or_default()
}

/// econtrib(e1, e2, w1, r)。
fn econtrib(e1: u32, e2: u32, w1: i64, r: i64) -> Q {
    let d = (e1 + e2) as usize;
    let mut a: Vec<Q> = Vec::with_capacity(d + 1);
    for m in 0..=d as u32 {
        let b = bern_poly(m + 2, &(Q::from(w1) / Q::from(r)));
        let mut v = b * Q::from(if m % 2 == 0 { 1 } else { -1 });
        v = v / Q::from((m + 1) as i64) / Q::from((m + 2) as i64);
        a.push(v);
    }
    apoly_eval(&edge_coeff(e1, e2), &a)
}

/// Chiodo 形式的 DR（固定 r，无插值；无 2^{−d}）。返回类（我们的 aut 口径）。
/// 单生成元的 Chiodo 系数（与 admcycles DR_coeff_new 同口径：含 1/|Aut|）。
/// 经 dr_chiodo 全类计算后按图匹配；None = 零系数或图不在类中。
pub fn dr_coeff_one(
    generator: &DecStratum,
    g: u32,
    a: &[i64],
    d: u32,
    k: i64,
    r: i64,
) -> Option<Q> {
    let cls = dr_chiodo(g, a, d, k, r).ok()?;
    let want = &generator.poly.terms[0];
    let mut acc: Option<Q> = None;
    for t in &cls.terms {
        if t.graph != generator.graph {
            continue;
        }
        // 图相同时多项式多单项式——按 (κ, ψ) 键精确匹配
        for term in &t.poly.terms {
            let kap_eq =
                term.kappa == want.kappa
                    || (term.kappa.len() <= want.kappa.len()
                        && term.kappa.iter().enumerate().all(|(i, k)| {
                            want.kappa.get(i).map(|w| w == k).unwrap_or(k.is_empty())
                        }));
            if kap_eq && term.psi == want.psi {
                acc = Some(match acc {
                    Some(x) => x + &term.coeff,
                    None => term.coeff.clone(),
                });
            }
        }
    }
    acc.map(|c| c / Q::from(generator.graph.automorphism_number() as i64))
}

pub fn dr_chiodo(g: u32, a: &[i64], d: u32, k: i64, r: i64) -> Result<TautClass, String> {
    let n = a.len() as u32;
    if r <= 0 {
        return Err("r 必须为正".into());
    }
    let asum: i64 = a.iter().sum();
    if (asum - k * (2 * g as i64 - 2 + n as i64)) % r != 0 {
        return Err(format!(
            "ΣA={asum} 与 k(2g−2+n)={} 需模 {r} 同余",
            k * (2 * g as i64 - 2 + n as i64)
        ));
    }
    let mut terms: Vec<DecStratum> = Vec::new();
    for generator in tautgens(g, n, d) {
        let (coeff, _) = chiodo_coeff_at(&generator, g, a, d, k, r, false)?;
        if !is_zero(&coeff) {
            terms.push(DecStratum {
                graph: generator.graph.clone(),
                poly: KP::new(vec![KpTerm {
                    kappa: generator.poly.terms[0].kappa.clone(),
                    psi: generator.poly.terms[0].psi.clone(),
                    coeff,
                }]),
            });
        }
    }
    Ok(TautClass::new(terms))
}

/// 单生成元的 Chiodo 系数在固定 m 处的值 = scalar(0,m)·total·m^{2d−h0}。
/// spin=true 时仅数奇权重（每条边流 w 均为奇；配合偶 m 使用——与
/// admcycles mpoly_special 的 odd 过滤一致，m 偶时 −w 与 w 同奇偶）。
fn chiodo_coeff_at(
    gs: &DecStratum,
    _g: u32,
    a: &[i64],
    d: u32,
    k: i64,
    m: i64,
    spin: bool,
) -> Result<(Q, i64), String> {
    let n = a.len() as u32;
    let graph = &gs.graph;
    let nv = graph.vertices().len();
    let t = &gs.poly.terms[0];
    let mut scalar = Q::ONE;
    let mut given = vec![0i64; nv];
    // κ：B 项（逐顶点逐 κ 指标）
    for (vi, v) in graph.vertices().iter().enumerate() {
        given[vi] = -k * (2 * v.genus as i64 - 2 + v.legs.len() as i64);
        let kap = &t.kappa[vi];
        // 统计各指标出现次数
        let mut counts: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
        for &idx in kap {
            *counts.entry(idx).or_insert(0) += 1;
        }
        for (&idx, &ex) in &counts {
            let mm = idx - 1; // 代码中 m 为 0 基
            let b = bern_poly(mm + 2, &(Q::from(k) / Q::from(m)));
            let mut vfac = b * Q::from(if (mm + 1) % 2 == 0 { 1 } else { -1 });
            vfac = vfac / Q::from((mm + 1) as i64) / Q::from((mm + 2) as i64);
            let mut p = Q::ONE;
            for _ in 0..ex {
                p *= &vfac;
            }
            let mut exf = Q::ONE;
            for i in 2..=ex {
                exf *= Q::from(i as i64);
            }
            scalar = scalar * p / exf;
        }
    }
    // 标记：A 权重无条件计入 given（所有标记，与装饰无关——差分定位的
    // 真根因：此前误绑在 ψ 装饰循环里，未装饰标记的 A 贡献被丢弃）
    for (li, &av) in a.iter().enumerate() {
        let l = li as u32 + 1;
        let vi = graph
            .vertices()
            .iter()
            .position(|v| v.legs.contains(&l))
            .ok_or(format!("标记 {l} 不在图中"))?;
        given[vi] += av;
    }
    // 标记 ψ 装饰：分划和因子（仅对有装饰的标记）
    for &(l, pw) in &t.psi {
        if l <= n {
            // Σ_{par ⊢ pw} Π_m [c_m^{ex}/ex!]，c_m = (−1)^m B_{m+2}(a_l/r)/((m+1)(m+2))
            let cm = |mm: u32| -> Q {
                let b = bern_poly(mm + 2, &(Q::from(a[(l - 1) as usize]) / Q::from(m)));
                let v = b * Q::from(if mm % 2 == 0 { 1 } else { -1 });
                v / Q::from((mm + 1) as i64) / Q::from((mm + 2) as i64)
            };
            // 递归枚举分划的指数向量（c_0^{e0} c_1^{e1} ...，Σ(m+1)e_m = pw）
            fn partsum(cm: &dyn Fn(u32) -> Q, rem: u32, mm: u32) -> Q {
                if rem == 0 {
                    return Q::ONE;
                }
                if mm > rem {
                    return Q::ZERO;
                }
                let mut acc = Q::ZERO;
                let mut e = 0u32;
                let mut cfac = Q::ONE;
                loop {
                    let used = e * (mm + 1);
                    if used > rem {
                        break;
                    }
                    if e > 0 {
                        let sub = partsum(cm, rem - used, mm + 1);
                        if sub != Q::ZERO {
                            let mut exf = Q::ONE;
                            for i in 2..=e {
                                exf *= Q::from(i as i64);
                            }
                            acc += cfac.clone() / exf * sub;
                        }
                    } else {
                        acc += partsum(cm, rem, mm + 1);
                    }
                    e += 1;
                    cfac *= cm(mm);
                }
                acc
            }
            scalar *= partsum(&cm, pw, 0);
        }
    }
    // 结点 exp 与权重和
    let mut edges_v: Vec<(usize, usize)> = Vec::new();
    let mut exps_fine: Vec<(u32, u32)> = Vec::new();
    let leg_vertex = |l: u32| -> usize {
        graph
            .vertices()
            .iter()
            .position(|v| v.legs.contains(&l))
            .unwrap()
    };
    for &(ea, eb) in graph.edges() {
        edges_v.push((leg_vertex(ea), leg_vertex(eb)));
        exps_fine.push((t.psi_power_at(ea), t.psi_power_at(eb)));
    }

    let s: i64 = given.iter().sum();
    given[0] -= s;
    let ne = edges_v.len();
    let h0 = (graph.edges().len() as i64) - (nv as i64) + 1;
    let mut total = Q::ZERO;
    let mut w = vec![0i64; ne];
    fn recw(
        i: usize,
        w: &mut Vec<i64>,
        state: &mut Vec<i64>,
        edges: &[(usize, usize)],
        fine: &[(u32, u32)],
        m: i64,
        acc: &Q,
        total: &mut Q,
        spin: bool,
    ) {
        if i == edges.len() {
            if state.iter().all(|&x| x == 0) {
                if spin && w.iter().any(|&x| x % 2 == 0) {
                    return;
                }
                let mut term = acc.clone();
                for (j, &(e1, e2)) in fine.iter().enumerate() {
                    // 我们枚举的 w 满足 Bw = −given，参考的 v 满足 Bv = given，
                    // 故 w = −v；参考代入 −v，我们应代入 w（codex 仲裁 2026-10-05）
                    let w1 = w[j].rem_euclid(m);
                    term *= econtrib(e1, e2, w1, m);
                }
                *total += term;
            }
            return;
        }
        let (u, v) = edges[i];
        for x in 0..m {
            w[i] = x;
            state[u] = (state[u] + x) % m;
            state[v] = ((state[v] - x) % m + m) % m;
            recw(i + 1, w, state, edges, fine, m, acc, total, spin);
            state[u] = ((state[u] - x) % m + m) % m;
            state[v] = (state[v] + x) % m;
        }
    }
    let mut state: Vec<i64> = given.iter().map(|&x| x.rem_euclid(m)).collect();
    recw(
        0,
        &mut w,
        &mut state,
        &edges_v,
        &exps_fine,
        m,
        &Q::ONE,
        &mut total,
        spin,
    );
    let mut mpow = Q::ONE;
    let e = 2 * d as i64 - h0;
    if e >= 0 {
        for _ in 0..e {
            mpow *= Q::from(m);
        }
    } else {
        for _ in 0..(-e) {
            mpow /= Q::from(m);
        }
    }
    Ok((scalar * total * mpow, h0))
}

/// Pixton 多项式插值模式（= admcycles DR_cycle(chiodo_coeff=True) 主路径）：
/// 在 m = start..（2d+2 个点，chiodo 保险多点）逐点求 chiodo_coeff_at，
/// 拉格朗日插值 P(r) 后取 P(0)。spin=true（CSS21 spin DR cycle）：m 取偶数、
/// 仅奇权重、因子 1/2^{g−h0}，要求 A 全奇且 k 奇。
pub fn dr_chiodo_poly(g: u32, a: &[i64], d: u32, k: i64, spin: bool) -> Result<TautClass, String> {
    let n = a.len() as u32;
    let asum: i64 = a.iter().sum();
    let chi = 2 * g as i64 - 2 + n as i64;
    if chi <= 0 {
        return Err(format!("M̄_({g},{n}) 不稳定"));
    }
    // r → 0 的多项式插值要求等式精确成立（模 r 同余对一切大 r 成立 ⇔ 相等）
    if asum != k * chi {
        return Err(format!("ΣA={asum} 须等于 k(2g−2+n)={}", k * chi));
    }
    if spin {
        if a.iter().any(|&x| x.rem_euclid(2) != 1) {
            return Err("spin DR 的分歧数据须全为奇数".into());
        }
        if k.rem_euclid(2) != 1 {
            return Err("spin DR 的 k 须为奇数".into());
        }
    }
    let abs_sum: i64 = a.iter().map(|x| x.abs()).sum();
    let m0 = (abs_sum + 1) / 2 + g as i64 * k.abs();
    let start = m0 + 1;
    let npts = 2 * d as usize + 2;
    let mut terms: Vec<DecStratum> = Vec::new();
    for generator in tautgens(g, n, d) {
        // 逐点求值
        let mut pts: Vec<i64> = Vec::with_capacity(npts);
        let mut vals: Vec<Q> = Vec::with_capacity(npts);
        for i in 0..npts {
            let m = if spin {
                2 * (start + i as i64)
            } else {
                start + i as i64
            };
            let (c, h0) = chiodo_coeff_at(&generator, g, a, d, k, m, spin)?;
            // spin 因子 1/2^{g−h0}
            let mut c = c;
            if spin {
                for _ in 0..(g as i64 - h0) {
                    c /= Q::from(2);
                }
            }
            pts.push(m);
            vals.push(c);
        }
        // 拉格朗日 P(0) = Σ_i v_i Π_{j≠i} (−m_j)/(m_i−m_j)
        let mut p0 = Q::ZERO;
        for i in 0..npts {
            let mut wt = Q::ONE;
            for j in 0..npts {
                if i != j {
                    wt = wt * Q::from(-pts[j]) / Q::from(pts[i] - pts[j]);
                }
            }
            p0 += vals[i].clone() * wt;
        }
        if !is_zero(&p0) {
            terms.push(DecStratum {
                graph: generator.graph.clone(),
                poly: KP::new(vec![KpTerm {
                    kappa: generator.poly.terms[0].kappa.clone(),
                    psi: generator.poly.terms[0].psi.clone(),
                    coeff: p0,
                }]),
            });
        }
    }
    Ok(TautClass::new(terms))
}

/// Norbury Θ_{g,n} = 2^{g−1+n} · 2^{2g−2d−1} · (−1)^d · DR_chiodo(g, (1,..,1), d=2g−2+n, k=−1, r=2)。
/// 注意 2g−2d−1 = −(2g+2n−3) 恒为负，故 2 的幂为分母。
pub fn theta_class(g: u32, n: u32) -> Result<TautClass, String> {
    if 2 * g as i64 - 2 + n as i64 <= 0 {
        return Err(format!("M̄_({g},{n}) 不稳定"));
    }
    let d = 2 * g + n - 2;
    let a = vec![1i64; n as usize];
    let cls = dr_chiodo(g, &a, d, -1, 2)?;
    // 2^{g−1+n} / 2^{2g+2n−3} · (−1)^d
    let mut factor = Q::ONE;
    for _ in 0..(g + n - 1) {
        factor *= Q::from(2i64);
    }
    for _ in 0..(2 * g + 2 * n - 3) {
        factor /= Q::from(2i64);
    }
    if d % 2 == 1 {
        factor = -factor;
    }
    Ok(cls.scale(&factor))
}

#[cfg(test)]
mod theta_tests {
    use super::*;

    /// 逐项锚点（admcycles 实测 2026-10-05）：DR_chiodo(1,(1,),1,−1,2) =
    /// −11/6 κ₁ − 1/6 ψ₁ − 1/24[Δ₀]（其 aut 口径含 1/2，我们口径 Δ₀ 系数 −1/12）；
    /// Θ_{1,1} = 上述取反。注意 Θ_{1,1}=3ψ₁ 只在模 FZ 关系下成立，逐项不等。
    #[test]
    fn theta11_terms_vs_admcycles() {
        let d = dr_chiodo(1, &[1], 1, -1, 2).unwrap();
        assert_eq!(d.terms.len(), 2);
        for t in &d.terms {
            let g = &t.graph;
            let c = &t.poly.terms[0].coeff;
            if g.vertices().len() == 1 && g.edges().is_empty() {
                // fundclass：κ₁ −11/6 + ψ₁ −1/6
                assert_eq!(t.poly.terms.len(), 2);
                assert!(
                    t.poly
                        .terms
                        .iter()
                        .any(|x| x.psi == vec![(1, 1)] && x.coeff.to_string() == "-1/6")
                );
                assert!(
                    t.poly
                        .terms
                        .iter()
                        .any(|x| x.kappa == vec![vec![1]] && x.coeff.to_string() == "-11/6")
                );
                let _ = c;
            } else {
                // Δ₀（aut=2）：我们口径 −1/12
                assert_eq!(c.to_string(), "-1/12");
            }
        }
        let th = theta_class(1, 1).unwrap();
        assert_eq!(th.terms.len(), 2);
        let _ = d;
    }

    /// 积分锚点（admcycles 实测；g=1 用环关系 Θ=3ψ 的 ∫=1/8 交叉验证）。
    #[test]
    fn theta_integrals_vs_admcycles() {
        let t11 = theta_class(1, 1).unwrap();
        assert_eq!(int_against(&t11, &[], &[]).unwrap().to_string(), "1/8");
        let t21 = theta_class(2, 1).unwrap();
        assert_eq!(
            int_against(&t21, &[], &[(1, 1)]).unwrap().to_string(),
            "3/128"
        );
        let t20 = theta_class(2, 0).unwrap();
        assert_eq!(
            int_against(&t20, &[(1, 1)], &[]).unwrap().to_string(),
            "3/128"
        );
        let t30 = theta_class(3, 0).unwrap();
        assert_eq!(
            int_against(&t30, &[(1, 2)], &[]).unwrap().to_string(),
            "111/1024"
        );
        let t31 = theta_class(3, 1).unwrap();
        assert_eq!(
            int_against(&t31, &[], &[(1, 2)]).unwrap().to_string(),
            "15/1024"
        );
    }
}

// ===================== 超椭圆配对（Z/2 Hurwitz 特化） =====================

/// ∫_{\bar H_{hyperell(g,n,m)}} π^*X，X 为 M̄_{g,N} 基本类上的 κ-ψ 单项式。
///
/// \bar H = Z/2 双重覆盖（源亏格 g，2g+2 个分支点，m 对共轭点；前 n 个分支
/// 点与全部共轭点标记）。传输公式（经 oracle Hpullpush 三锚点验证）：
///   pair = deltadeg(=1/2) · 2^{Σκ幂次} · (1/2)^{Σψ_fixed} · ∫_{M̄_{0,2g+2+2m}} tr(X)
/// tr：κ_j → κ_j（指标保持，商顶点，因子 2/每个幂次）；ψ（固定/分支标记）→ 对应
/// 分支腿（因子 1/2/幂）；ψ（共轭标记）→ 对应腿（无因子）。经 oracle
/// Hpullpush 五锚点验证：g=2,n=6: ψ₁³→1/16、ψ₁ψ₂ψ₃→3/8、κ₃→1、κ₁ψ₁²→3/2、
/// κ₁ψ₁ψ₂→3、κ₂κ₁→18。
/// ∫_{ar H} π^*X（Z/2 平凡覆盖配对，共通内核）。
///
/// 商空间腿：branch 个分支腿（前 n 个是 X 的固定标记）+ m 个轨道腿（每对共轭
/// 合并到一腿）。传输：κ_j 指标保持 ×2^{Σ幂}；ψ_固定 → 分支腿 ×(1/2)^{幂}；
/// ψ_共轭 → 轨道腿（两腿求和，无因子）。deltadeg = 2^{2·g_q+m−1}。
fn z2_pair(
    g_quot: u32,
    branch: u32, // 分支点数（=商上标记的固定腿数上限）
    _n: u32,     // X 的固定标记数（≤ branch；腿号已是全空间编号，此处不需要）
    m: u32,      // 共轭对数
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<sc_num::Q, String> {
    let full = branch + 2 * m; // 全空间腿数（分支 1..branch，共轭 branch+1..branch+2m）
    let q_legs = branch + m; // 商空间腿数
    // 腿映射 → 商腿幂
    let mut qpsi = vec![0u32; q_legs as usize];
    let mut psi_fixed_p = 0u32;
    for &(l, p) in psi {
        if l >= 1 && l <= branch {
            qpsi[(l - 1) as usize] += p; // 分支腿（全部带 1/2 因子）
            psi_fixed_p += p;
        } else if l > branch && l <= full {
            // 共轭腿：第 ((l-branch-1)/2)+1 对 → 轨道腿 branch + pair_idx
            let pair_idx = (l - branch - 1) / 2;
            qpsi[(branch + pair_idx) as usize] += p;
        } else {
            return Err(format!("标记腿 {l} 超出 1..={full}"));
        }
    }
    // 传输系数：deltadeg · 2^{Σκ幂} · (1/2)^{Σψ_固定幂}
    let kexp: u32 = kappa.iter().map(|&(_, p)| p).sum();
    let exp: i64 = 2 * g_quot as i64 + m as i64 - 1;
    let mut factor = Q::ONE;
    if exp >= 0 {
        for _ in 0..exp {
            factor *= Q::from(2);
        }
    } else {
        for _ in 0..(-exp) {
            factor /= Q::from(2);
        }
    }
    for _ in 0..kexp {
        factor *= Q::from(2);
    }
    for _ in 0..psi_fixed_p {
        factor /= Q::from(2);
    }
    // 商空间求值：κ^κ ψ^{qpsi}（度不符 = 配对对象度数错 → 0）
    let v = match crate::socle::socle_formula(g_quot, &qpsi, &kappa_parts(kappa)) {
        Ok(v) => v,
        Err(e) if e.contains("维度不符") => Q::ZERO,
        Err(e) => return Err(format!("z2_pair socle: {e}")),
    };
    Ok(factor * v)
}

/// 超椭圆配对：ar H = Z/2 双覆盖（源亏格 g，2g+2 分支，m 共轭对；商 = P¹）。
pub fn hyperell_pair(
    g: u32,
    n: u32,
    m: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<sc_num::Q, String> {
    if n > 2 * g + 2 {
        return Err(format!("固定标记 {n} > 分支点数 {}", 2 * g + 2));
    }
    z2_pair(0, 2 * g + 2, n, m, kappa, psi)
}

/// 双椭圆配对：源亏格 g 的 C 双覆盖椭圆曲线（2g−2 分支，m 共轭对；商亏格 1）。
pub fn biell_pair(
    g: u32,
    n: u32,
    m: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<sc_num::Q, String> {
    if g < 2 {
        return Err("双椭圆需 g ≥ 2".into());
    }
    if n > 2 * g - 2 {
        return Err(format!("固定标记 {n} > 分支点数 {}", 2 * g - 2));
    }
    z2_pair(1, 2 * g - 2, n, m, kappa, psi)
}

fn kappa_parts(kappa: &[(u32, u32)]) -> Vec<u32> {
    let mut p = Vec::new();
    for &(i, c) in kappa {
        for _ in 0..c {
            p.push(i);
        }
    }
    p.sort_unstable();
    p
}

#[cfg(test)]
mod hyperell_tests {
    use super::*;

    /// oracle Hpullpush 锚点。超椭圆 g=2 m=0（六）+ g=1 m=1 共轭对（二）。
    #[test]
    fn hyperell_pair_anchors() {
        assert_eq!(
            hyperell_pair(2, 6, 0, &[], &[(1, 3)]).unwrap().to_string(),
            "1/16"
        );
        assert_eq!(
            hyperell_pair(2, 6, 0, &[], &[(1, 1), (2, 1), (3, 1)])
                .unwrap()
                .to_string(),
            "3/8"
        );
        assert_eq!(
            hyperell_pair(2, 6, 0, &[(3, 1)], &[]).unwrap().to_string(),
            "1"
        );
        assert_eq!(
            hyperell_pair(2, 6, 0, &[(1, 1)], &[(1, 2)])
                .unwrap()
                .to_string(),
            "3/2"
        );
        assert_eq!(
            hyperell_pair(2, 6, 0, &[(1, 1)], &[(1, 1), (2, 1)])
                .unwrap()
                .to_string(),
            "3"
        );
        assert_eq!(
            hyperell_pair(2, 6, 0, &[(2, 1), (1, 1)], &[])
                .unwrap()
                .to_string(),
            "18"
        );
        // g=1 m=1 共轭对（X 标记：1 固定、5/6 共轭）
        assert_eq!(
            hyperell_pair(1, 1, 1, &[], &[(5, 1), (6, 1)])
                .unwrap()
                .to_string(),
            "1"
        );
        assert_eq!(
            hyperell_pair(1, 1, 1, &[], &[(5, 2)]).unwrap().to_string(),
            "1"
        );
        // 度不符 → 0
        assert_eq!(
            hyperell_pair(1, 1, 1, &[], &[(5, 1)]).unwrap().to_string(),
            "0"
        );
    }

    /// spin DR 锚点：admcycles doctest DR_cycle(1,(2a+1,−2a+1),spin=True) 取 a=0。
    /// oracle 展示（generator 系数含 1/aut）：主图 −1/4κ₁+1/4ψ₁+1/4ψ₂；
    /// loop 图 1/24（aut=2）；分离图 −1/4。我们的类系数 = aut × oracle。
    #[test]
    fn spin_dr_anchor() {
        let cls = dr_chiodo_poly(1, &[1, 1], 1, 1, true).unwrap();
        let mut got: Vec<(
            usize,
            usize,
            Vec<u32>,
            Vec<Vec<u32>>,
            Vec<(u32, u32)>,
            String,
        )> = Vec::new();
        for t in &cls.terms {
            let g = &t.graph;
            for term in &t.poly.terms {
                got.push((
                    g.vertices().len(),
                    g.edges().len(),
                    g.vertices().iter().map(|v| v.genus).collect(),
                    term.kappa.clone(),
                    term.psi.clone(),
                    term.coeff.to_string(),
                ));
            }
        }
        let coeff_of = |pred: &dyn Fn(
            &(
                usize,
                usize,
                Vec<u32>,
                Vec<Vec<u32>>,
                Vec<(u32, u32)>,
                String,
            ),
        ) -> bool|
         -> String {
            got.iter()
                .filter(|t| pred(t))
                .map(|(_, _, _, _, _, c)| c.clone())
                .next()
                .unwrap_or_else(|| "<缺失>".into())
        };
        // 主图 [1][[1,2]]：κ₁ → −1/4；ψ₁ → 1/4；ψ₂ → 1/4
        assert_eq!(
            coeff_of(&|t| t.0 == 1 && t.1 == 0 && t.2 == vec![1] && t.3 == vec![vec![1]]),
            "-1/4",
            "{got:?}"
        );
        assert_eq!(
            coeff_of(&|t| t.0 == 1
                && t.1 == 0
                && t.3 == vec![Vec::<u32>::new()]
                && t.4 == vec![(1, 1)]),
            "1/4",
            "{got:?}"
        );
        assert_eq!(
            coeff_of(&|t| t.0 == 1
                && t.1 == 0
                && t.3 == vec![Vec::<u32>::new()]
                && t.4 == vec![(2, 1)]),
            "1/4",
            "{got:?}"
        );
        // loop 图 [0][[4,5,1,2]] 自环：oracle 1/24 × aut 2 = 1/12
        assert_eq!(
            coeff_of(&|t| t.0 == 1
                && t.1 == 1
                && t.2 == vec![0]
                && t.3 == vec![Vec::<u32>::new()]
                && t.4.is_empty()),
            "1/12",
            "{got:?}"
        );
        // 分离图 [0,1] 一边无装饰：oracle −1/4 × aut 1
        assert_eq!(
            coeff_of(&|t| t.0 == 2
                && t.1 == 1
                && t.2 == vec![0, 1]
                && t.3 == vec![Vec::<u32>::new(), Vec::<u32>::new()]
                && t.4.is_empty()),
            "-1/4",
            "{got:?}"
        );
    }

    /// drpoly/drspin 积分差分（tools/spin_oracle.py + g=1 deg-1 批，15 例全对）。
    #[test]
    fn drpoly_drspin_int_oracle() {
        // Chiodo 插值主路径（= DR_cycle(chiodo_coeff=True)）
        let dp = |a: &[i64], k: i64, kap: &[(u32, u32)], psi: &[(u32, u32)]| {
            let cls = dr_chiodo_poly(1, a, 1, k, false).unwrap();
            int_against(&cls, kap, psi).unwrap().to_string()
        };
        assert_eq!(dp(&[2, 3, -5], 0, &[], &[(1, 1), (2, 1)]), "3/2");
        assert_eq!(dp(&[2, 3, -5], 0, &[(2, 1)], &[]), "37/24");
        assert_eq!(dp(&[2, 3, -5], 0, &[], &[(3, 2)]), "1/2");
        assert_eq!(dp(&[2, 3, -5], 0, &[(1, 1)], &[(2, 1)]), "31/8");
        // spin DR g=1（X deg 1）
        let ds = |a: &[i64], kap: &[(u32, u32)], psi: &[(u32, u32)]| {
            let cls = dr_chiodo_poly(1, a, 1, 1, true).unwrap();
            int_against(&cls, kap, psi).unwrap().to_string()
        };
        assert_eq!(ds(&[1, 1], &[(1, 1)], &[]), "1/24");
        assert_eq!(ds(&[1, 1], &[], &[(1, 1)]), "1/24");
        assert_eq!(ds(&[1, 1], &[], &[(2, 1)]), "1/24");
        assert_eq!(ds(&[3, -1], &[(1, 1)], &[]), "5/24");
        assert_eq!(ds(&[3, -1], &[], &[(1, 1)]), "1/8");
        assert_eq!(ds(&[3, -1], &[], &[(2, 1)]), "1/8");
    }

    /// 16 例 oracle 全量对拍（tools/hyb_oracle.py，admcycles Hyperell/Biell
    /// 类 evaluate）。合成公式：1/(branch−n)! × Σ z2_pair(ρ^*X 单顶点展开)。
    #[test]
    fn z2_int_oracle_grid() {
        // Hyperell(1,1,1)，deg X = 2（psi-only + κ 混合 + κ-only）
        let h =
            |k: &[(u32, u32)], p: &[(u32, u32)]| hyperell_int(1, 1, 1, k, p).unwrap().to_string();
        assert_eq!(h(&[], &[(1, 2)]), "1/24");
        assert_eq!(h(&[], &[(1, 1), (2, 1)]), "1/6");
        assert_eq!(h(&[], &[(2, 2)]), "1/6");
        assert_eq!(h(&[], &[(2, 1), (3, 1)]), "1/6");
        assert_eq!(h(&[(2, 1)], &[]), "5/24");
        assert_eq!(h(&[(1, 1)], &[(1, 1)]), "1/4");
        assert_eq!(h(&[(1, 1)], &[(2, 1)]), "1/2");
        // Hyperell(2,0,0)（= [M̄₂] 数值），deg X = 3：κ-only
        assert_eq!(
            hyperell_int(2, 0, 0, &[(3, 1)], &[]).unwrap().to_string(),
            "1/1152"
        );
        assert_eq!(
            hyperell_int(2, 0, 0, &[(1, 3)], &[]).unwrap().to_string(),
            "43/2880"
        );
        assert_eq!(
            hyperell_int(2, 0, 0, &[(1, 1), (2, 1)], &[])
                .unwrap()
                .to_string(),
            "1/240"
        );
        // Biell(2,0,0)（含 admcycles 双重计数 /2 特判），deg X = 2
        assert_eq!(
            biell_int(2, 0, 0, &[(2, 1)], &[]).unwrap().to_string(),
            "1/32"
        );
        assert_eq!(
            biell_int(2, 0, 0, &[(1, 2)], &[]).unwrap().to_string(),
            "5/48"
        );
        // Biell(3,0,0)，deg X = 4
        assert_eq!(
            biell_int(3, 0, 0, &[(4, 1)], &[]).unwrap().to_string(),
            "7/1152"
        );
        assert_eq!(
            biell_int(3, 0, 0, &[(2, 2)], &[]).unwrap().to_string(),
            "109/1152"
        );
        // 放大差分（batch 7/8 摘要）：固定标记 ψ 深度、m=1 混合、g=3 高次 κ
        assert_eq!(
            hyperell_int(2, 1, 0, &[], &[(1, 3)]).unwrap().to_string(),
            "1/1920"
        );
        assert_eq!(
            hyperell_int(3, 0, 0, &[(5, 1)], &[]).unwrap().to_string(),
            "1/46080"
        );
        assert_eq!(
            hyperell_int(2, 1, 1, &[], &[(1, 2), (2, 1), (3, 1)])
                .unwrap()
                .to_string(),
            "1/80"
        );
        assert_eq!(
            biell_int(2, 1, 0, &[], &[(1, 2)]).unwrap().to_string(),
            "1/48"
        );
        assert_eq!(
            biell_int(2, 1, 0, &[(2, 1)], &[]).unwrap().to_string(),
            "7/48"
        );
        assert_eq!(
            biell_int(2, 1, 0, &[(1, 1)], &[(1, 1)])
                .unwrap()
                .to_string(),
            "7/48"
        );
    }
}

/// X 紧凑腿号（固定 1..n，共轭 n+1..n+2m）→ 全空间（分支 1..branch，共轭 branch+1..）。
fn remap_conj(psi: &[(u32, u32)], n: u32, branch: u32) -> Vec<(u32, u32)> {
    psi.iter()
        .map(|&(l, p)| if l <= n { (l, p) } else { (l - n + branch, p) })
        .collect()
}

/// weak compositions of `total` into `parts` 份（有序）。
fn comps(total: u32, parts: u32) -> Vec<Vec<u32>> {
    if parts == 1 {
        return vec![vec![total]];
    }
    let mut out = Vec::new();
    for e in 0..=total {
        for mut rest in comps(total - e, parts - 1) {
            rest.insert(0, e);
            out.push(rest);
        }
    }
    out
}

fn fact_u32(f: u32) -> Q {
    let mut v = Q::ONE;
    for i in 2..=f {
        v *= Q::from(i as i64);
    }
    v
}

/// Z/2 配对积分合成（Hyperell/Biell 共通；16 例 oracle 对拍全中）：
///   ∫ Locus(g,n,m)·X = 1/(branch−n)! · Σ_terms coeff · z2_pair(term)
/// ρ^*X 展开（忘却拉回忘掉未标记分支腿 n+1..branch）：ψ_l 稳定；
/// κ_j^c → (κ_j − Σ_a ψ_a^j)^c 的多项式展开（admcycles
/// forgetful_pullback_list 的单顶点部分；其有理尾巴图项在 H̄ 上为空层
/// ——σ 不动点计数矛盾——配对为 0，实测 16 例均无需）。Biell(2,0,0) 需
/// 额外 ×1/2（admcycles 双重计数特判：无标记 g=2 双椭圆曲线有两种覆盖结构）。
fn z2_int(
    g_quot: u32,
    branch: u32,
    n: u32,
    m: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    extra_half: bool,
) -> Result<sc_num::Q, String> {
    // 合并同指标的 κ 计数
    let mut kap: Vec<(u32, u32)> = Vec::new();
    for &(j, c) in kappa {
        if let Some(e) = kap.iter_mut().find(|(jj, _)| *jj == j) {
            e.1 += c;
        } else {
            kap.push((j, c));
        }
    }
    let forgotten: Vec<u32> = (n + 1..=branch).collect();
    // (coeff, kappa 项, 追加 ψ 项) 逐 κ 指标展开
    let mut terms: Vec<(Q, Vec<(u32, u32)>, Vec<(u32, u32)>)> =
        vec![(Q::ONE, vec![], remap_conj(psi, n, branch))];
    for &(j, c) in &kap {
        let mut next = Vec::new();
        for (coeff, kterm, pterm) in &terms {
            for comp in comps(c, forgotten.len() as u32 + 1) {
                let moved = c - comp[0];
                let mut cf = coeff * fact_u32(c);
                for &e in &comp {
                    cf /= fact_u32(e);
                }
                if moved % 2 == 1 {
                    cf *= Q::from(-1);
                }
                let mut k2 = kterm.clone();
                if comp[0] > 0 {
                    k2.push((j, comp[0]));
                }
                let mut p2 = pterm.clone();
                for (ai, &e) in comp[1..].iter().enumerate() {
                    if e > 0 {
                        p2.push((forgotten[ai], j * e));
                    }
                }
                next.push((cf, k2, p2));
            }
        }
        terms = next;
    }
    let mut acc = Q::ZERO;
    for (coeff, kterm, pterm) in &terms {
        let v = z2_pair(g_quot, branch, n, m, kterm, pterm)?;
        acc += coeff * v;
    }
    let mut res = acc / fact_u32(branch - n);
    if extra_half {
        res /= Q::from(2);
    }
    Ok(res)
}

/// ∫_{M̄_{g,n+2m}} Hyperell(g,n,m) · X（X 为 κ-ψ 单项式，紧凑腿号）。
pub fn hyperell_int(
    g: u32,
    n: u32,
    m: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<sc_num::Q, String> {
    if n > 2 * g + 2 {
        return Err(format!("固定标记 {n} > 分支点数 {}", 2 * g + 2));
    }
    z2_int(0, 2 * g + 2, n, m, kappa, psi, false)
}

/// ∫_{M̄_{g,n+2m}} Biell(g,n,m) · X（X 为 κ-ψ 单项式，紧凑腿号）。
pub fn biell_int(
    g: u32,
    n: u32,
    m: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<sc_num::Q, String> {
    if g < 2 {
        return Err("双椭圆需 g ≥ 2".into());
    }
    if n > 2 * g - 2 {
        return Err(format!("固定标记 {n} > 分支点数 {}", 2 * g - 2));
    }
    z2_int(1, 2 * g - 2, n, m, kappa, psi, g == 2 && n == 0 && m == 0)
}
