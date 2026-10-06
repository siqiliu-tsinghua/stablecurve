//! GRRcomp：广义 Hodge Chern 类（PRvZ20 Theorem 1；admcycles/GRRcomp.py 移植）。
//!
//! D = l·K̃ + Σ d_i σ_i + Σ a_{h,S} C_{h,S} 上的 ch(R^•π_* O(D))。
//! 链图构造用 new_raw（保持源码腿号 n+1..n+2L+2，ψ 装饰按其书写——Θ 疑案
//! 的老坑：StableGraph::new 规范化重编号结点腿致装饰失配）。

use crate::graph::StableGraph;
use crate::kppoly::{KappaPsiPolynomial as KP, KpTerm};
use crate::socle::bernoulli;
use crate::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero};

/// (h, S, a_{h,S})：泛曲线边界除子 C_{h,S} 的系数。
pub type ATriple = (u32, Vec<u32>, i64);

fn fact(n: u32) -> Q {
    let mut v = Q::ONE;
    for i in 2..=n {
        v *= Q::from(i as i64);
    }
    v
}

fn binom(a: u32, b: u32) -> Q {
    if b > a {
        return Q::ZERO;
    }
    let mut v = Q::ONE;
    for i in 0..b {
        v = v * Q::from((a - i) as i64) / Q::from((b - i) as i64);
    }
    v
}

fn q_powi(base: i64, e: u32) -> Q {
    let mut v = Q::ONE;
    for _ in 0..e {
        v *= Q::from(base);
    }
    v
}

/// B_n/n!（divided_bernoulli）。
pub fn divided_bernoulli(n: u32) -> Q {
    bernoulli(n as u64) / fact(n)
}

/// Bernoulli 多项式 B_m(x) 数值。
fn bern_poly_val(m: u32, x: i64) -> Q {
    let mut v = Q::ZERO;
    for k in 0..=m {
        let t = binom(m, k) * bernoulli(k as u64);
        let mut xp = Q::ONE;
        for _ in 0..(m - k) {
            xp *= Q::from(x);
        }
        v += t * xp;
    }
    v
}

/// B_n(l)/n!（divided_bernoulli_polynomial）。
pub fn divided_bernoulli_poly(l: i64, n: u32) -> Q {
    bern_poly_val(n, l) / fact(n)
}

/// 有序链与正分配：a 的下标升序子列 i_1 < … < i_L（相邻两项 h 不减、S 包含，
/// 见 [`chain_step`]），连同 total 分成 L 个正整数的全部有序分配。
///
/// 参考实现枚举 total 分进 |a| 份的全部弱组合（sage IntegerVectors），再取
/// 非零分量为活跃项、检查链条件；两者一一对应，但 |a| 大时绝大多数组合不成链
/// （|a| = 27、total = 6 时组合数约 10⁶，成链的只有千余），故直接枚举链。
fn chains_with_parts(a: &[ATriple], total: u32) -> Vec<(Vec<usize>, Vec<u32>)> {
    fn extend(
        a: &[ATriple],
        total: u32,
        chain: &mut Vec<usize>,
        out: &mut Vec<(Vec<usize>, Vec<u32>)>,
    ) {
        let l = chain.len() as u32;
        if l > total {
            return;
        }
        // total 分成 l 个正整数的全部有序分配
        let mut parts = vec![0u32; l as usize];
        fn split(
            rem: u32,
            i: usize,
            parts: &mut Vec<u32>,
            chain: &[usize],
            out: &mut Vec<(Vec<usize>, Vec<u32>)>,
        ) {
            if i + 1 == parts.len() {
                parts[i] = rem;
                out.push((chain.to_vec(), parts.clone()));
                return;
            }
            let later = (parts.len() - i - 1) as u32;
            for p in 1..=rem - later {
                parts[i] = p;
                split(rem - p, i + 1, parts, chain, out);
            }
        }
        split(total, 0, &mut parts, chain, out);
        let last = *chain.last().unwrap();
        for j in last + 1..a.len() {
            if chain_step(&a[last], &a[j]) {
                chain.push(j);
                extend(a, total, chain, out);
                chain.pop();
            }
        }
    }
    let mut out = Vec::new();
    if total == 0 {
        return out;
    }
    for i in 0..a.len() {
        extend(a, total, &mut vec![i], &mut out);
    }
    out
}

fn multinomial(v: &[u32]) -> Q {
    let total: u32 = v.iter().sum();
    fact(total) / v.iter().map(|&x| fact(x)).fold(Q::ONE, |a, b| a * b)
}

/// 升序 powerset。
fn powerset_sorted(s: &[u32]) -> Vec<Vec<u32>> {
    let mut out = vec![Vec::new()];
    for &x in s {
        let mut add = Vec::new();
        for v in &out {
            let mut w = v.clone();
            w.push(x);
            add.push(w);
        }
        out.extend(add);
    }
    out
}

// ===================== kppoly 构造（admcycles kappacl/psicl 语义） =====================

fn kp_term(numvert: usize, kappa: Vec<Vec<u32>>, psi: Vec<(u32, u32)>, coeff: Q) -> KP {
    let _ = numvert;
    KP::new(vec![KpTerm { kappa, psi, coeff }])
}

fn kp_one(numvert: usize) -> KP {
    kp_term(numvert, vec![Vec::new(); numvert], Vec::new(), Q::ONE)
}

/// kappacl(vertex, index, numvert, g, n)：index==0 → (2g−2+n)·1。
fn kp_kappa(vertex: usize, index: u32, numvert: usize, g: u32, n: u32) -> KP {
    if index == 0 {
        return kp_term(
            numvert,
            vec![Vec::new(); numvert],
            Vec::new(),
            Q::from(2 * g as i64 - 2 + n as i64),
        );
    }
    let mut kappa = vec![Vec::new(); numvert];
    kappa[vertex] = vec![index];
    kp_term(numvert, kappa, Vec::new(), Q::ONE)
}

fn psi_pow(nv: usize, leg: u32, p: u32) -> KP {
    if p == 0 {
        return kp_one(nv);
    }
    let base = kp_term(nv, vec![Vec::new(); nv], vec![(leg, 1)], Q::ONE);
    let mut r = base.clone();
    for _ in 1..p {
        r = r.mul(&base);
    }
    r
}

fn kp_mul(a: &KP, b: &KP) -> KP {
    a.mul(b)
}

// ===================== 度数过滤 =====================

/// 保留度数恰为 t 的项（图边数 + κψ 度 = t）。
fn degree_part(cl: &TautClass, t: u32) -> TautClass {
    let mut terms = Vec::new();
    for ds in &cl.terms {
        let edges = ds.graph.edges().len() as u32;
        let mut keep = Vec::new();
        for term in &ds.poly.terms {
            let wd: u32 = term.kappa.iter().flatten().sum::<u32>()
                + term.psi.iter().map(|&(_, p)| p).sum::<u32>();
            if edges + wd == t {
                keep.push(term.clone());
            }
        }
        if !keep.is_empty() {
            terms.push(DecStratum {
                graph: ds.graph.clone(),
                poly: KP::new(keep),
            });
        }
    }
    TautClass::new(terms)
}

// ===================== todd_inv_sigma =====================

/// todd_inv_sigma(t)：奇 t 的 σ 图（自环 ×½ + 分离 ×1），
/// decsum = Σ_i (−ψ_{n+1})^i ψ_{n+2}^{t−1−i}。
fn todd_inv_sigma(t: u32, g: u32, n: u32) -> TautClass {
    if t % 2 == 0 {
        return TautClass::zero();
    }
    let decsum = |numvert: usize| -> KP {
        let mut acc = KP::zero();
        for i in 0..t {
            let mut term = kp_one(numvert);
            if i > 0 {
                let neg = psi_pow(numvert, n + 1, i).scale(&Q::from(-1));
                term = kp_mul(&term, &neg);
            }
            if t - 1 - i > 0 {
                term = kp_mul(&term, &psi_pow(numvert, n + 2, t - 1 - i));
            }
            acc = acc.add(&term);
        }
        acc
    };
    let mut parts: Vec<TautClass> = Vec::new();
    if g >= 1 {
        let legs: Vec<u32> = (1..=n + 2).collect();
        let graph = StableGraph::new(&[g - 1], &[legs], &[(n + 1, n + 2)]).unwrap();
        parts.push(
            TautClass::new(vec![DecStratum {
                graph,
                poly: decsum(1),
            }])
            .scale(&(Q::ONE / Q::from(2))),
        );
    }
    for h in 0..=g {
        for s in powerset_sorted(&(2..=n).collect::<Vec<u32>>()) {
            if (h == 0 && s.is_empty()) || (h == g && s.len() as i64 > n as i64 - 3) {
                continue;
            }
            let mut sl: Vec<u32> = vec![1];
            sl.extend(s.iter().copied());
            let comp: Vec<u32> = (1..=n).filter(|x| !sl.contains(x)).collect();
            let mut l1 = sl.clone();
            l1.push(n + 1);
            let mut l2 = comp;
            l2.push(n + 2);
            let graph = StableGraph::new(&[h, g - h], &[l1, l2], &[(n + 1, n + 2)]).unwrap();
            parts.push(TautClass::new(vec![DecStratum {
                graph,
                poly: decsum(2),
            }]));
        }
    }
    degree_part(&TautClass::sum(parts).scale(&divided_bernoulli(t + 1)), t)
}

// ===================== 链图构造 =====================

/// ψ 二项式和：Σ_j C(tv−1, j) ψ_{左}^j ψ_{右}^{tv−1−j}（度 tv−1）。
fn psi_binom(nv: usize, left: u32, right: u32, tv: u32) -> KP {
    let mut acc = KP::zero();
    for j in 0..tv {
        let mut term = kp_one(nv);
        if j > 0 {
            term = kp_mul(&term, &psi_pow(nv, left, j));
        }
        if tv - 1 - j > 0 {
            term = kp_mul(&term, &psi_pow(nv, right, tv - 1 - j));
        }
        acc = acc.add(&term.scale(&binom(tv - 1, j)));
    }
    acc
}

/// 公共头：Π_{i<upto} (−1)^{tv_i−1} · psi_binom(i)。
/// a_i^{tv_i} 不在此处——todd_inv_exp 外层统一乘（graphs 的 ψ 不含 a 因子）。
fn chain_psi_head(nv: usize, n: u32, tvn: &[u32], upto: usize) -> KP {
    let mut poly = kp_one(nv);
    for (i, &tv) in tvn.iter().enumerate() {
        if i >= upto {
            break;
        }
        let factor = psi_binom(nv, n + 2 * i as u32 + 1, n + 2 * i as u32 + 2, tv)
            .scale(&q_powi(-1, tv - 1));
        poly = kp_mul(&poly, &factor);
    }
    poly
}

/// 尾段 Σ_j (−ψ_{左})^j ψ_{右}^{t1−1−j}。
fn psi_tail(nv: usize, left: u32, right: u32, t1: u32) -> KP {
    let mut acc = KP::zero();
    for j in 0..t1 {
        let mut term = kp_one(nv);
        if j > 0 {
            term = kp_mul(&term, &psi_pow(nv, left, j).scale(&Q::from(-1)));
        }
        if t1 - 1 - j > 0 {
            term = kp_mul(&term, &psi_pow(nv, right, t1 - 1 - j));
        }
        acc = acc.add(&term);
    }
    acc
}

/// todd_inv_exp_graphs(a, t1, tvect, g, n)：σ 展开的链图组。
/// 锚点：([[1,[1],−1]], 1, [1], 3, 2) → 5 图（−ψ₃ / ½ / 1 / 1 / 1）。
fn todd_inv_exp_graphs(a: &[ATriple], t1: u32, tvect: &[u32], g: u32, n: u32) -> TautClass {
    let mut atn: Vec<(u32, Vec<u32>, i64)> = Vec::new();
    let mut tvn: Vec<u32> = Vec::new();
    for (ai, &tvi) in a.iter().zip(tvect) {
        if tvi != 0 {
            atn.push(ai.clone());
            tvn.push(tvi);
        }
    }
    let l = atn.len();
    if l == 0 {
        return TautClass::zero();
    }
    let ordered = (0..l.saturating_sub(1))
        .all(|i| atn[i].0 <= atn[i + 1].0 && atn[i].1.iter().all(|x| atn[i + 1].1.contains(x)));
    if !ordered {
        return TautClass::zero();
    }
    let mut parts: Vec<TautClass> = Vec::new();
    let push = |parts: &mut Vec<TautClass>,
                genera: Vec<u32>,
                labels: Vec<Vec<u32>>,
                edges: Vec<(u32, u32)>,
                poly: KP,
                scale: Q| {
        if let Ok(graph) = StableGraph::new_raw(&genera, &labels, &edges) {
            parts.push(TautClass::new(vec![DecStratum { graph, poly }]).scale(&scale));
        }
    };
    // 头段 labels/generas（L 顶点链 + 末顶点占位）
    let head_genera = |extra_last: u32| -> Vec<u32> {
        let mut v: Vec<u32> = vec![atn[0].0];
        for i in 1..l {
            v.push(atn[i].0 - atn[i - 1].0);
        }
        v.push(extra_last);
        v
    };
    let head_labels = |nextra: &[u32]| -> Vec<Vec<u32>> {
        let mut out: Vec<Vec<u32>> = Vec::new();
        let mut first = atn[0].1.clone();
        first.push(n + 1);
        out.push(first);
        for i in 1..l {
            let mut li: Vec<u32> = atn[i]
                .1
                .iter()
                .copied()
                .filter(|x| !atn[i - 1].1.contains(x))
                .collect();
            li.push(n + 2 * i as u32);
            li.push(n + 2 * i as u32 + 1);
            out.push(li);
        }
        let scomp: Vec<u32> = (1..=n).filter(|x| !atn[l - 1].1.contains(x)).collect();
        let mut last = scomp;
        for &x in nextra {
            last.push(x);
        }
        out.push(last);
        out
    };
    let head_edges = |count: u32| -> Vec<(u32, u32)> {
        (0..count).map(|i| (n + 2 * i + 1, n + 2 * i + 2)).collect()
    };

    let scomp: Vec<u32> = (1..=n).filter(|x| !atn[l - 1].1.contains(x)).collect();
    for h in atn[l - 1].0..=g {
        for s in powerset_sorted(&scomp) {
            let strictly = !(h == atn[l - 1].0 && s.is_empty())
                && !(h == g && (n as i64 - atn[l - 1].1.len() as i64 - s.len() as i64) < 2);
            if strictly {
                let mid_extra = vec![n + 2 * l as u32, n + 2 * l as u32 + 1];
                let fin_extra = vec![n + 2 * l as u32 + 2];
                let mut genera = head_genera(h - atn[l - 1].0);
                genera.push(g - h);
                // 拆 labels：头（末顶点替换为中继）+ 末顶点
                let mut labels = head_labels(&mid_extra);
                labels.pop();
                let mut mid: Vec<u32> = s.clone();
                mid.extend(mid_extra);
                labels.push(mid);
                let mut fin: Vec<u32> = scomp.iter().copied().filter(|x| !s.contains(x)).collect();
                fin.extend(fin_extra);
                labels.push(fin);
                let mut edges = head_edges(l as u32);
                edges.push((n + 2 * l as u32 + 1, n + 2 * l as u32 + 2));
                let nv = genera.len();
                let poly = kp_mul(
                    &chain_psi_head(nv, n, &tvn, l),
                    &psi_tail(nv, n + 2 * l as u32 + 1, n + 2 * l as u32 + 2, t1),
                );
                push(&mut parts, genera, labels, edges, poly, Q::ONE);
            }
            if h == atn[l - 1].0 && s.is_empty() {
                let genera = head_genera(g - atn[l - 1].0);
                let labels = head_labels(&[n + 2 * l as u32]);
                let edges = head_edges(l as u32);
                let nv = genera.len();
                let lp = n + 2 * l as u32 - 1;
                let rp = n + 2 * l as u32;
                let mut poly = chain_psi_head(nv, n, &tvn, l - 1);
                let fin = psi_pow(nv, lp, tvn[l - 1]).scale(&q_powi(-1, tvn[l - 1]));
                poly = kp_mul(&poly, &fin);
                poly = kp_mul(&poly, &psi_tail(nv, lp, rp, t1));
                push(&mut parts, genera, labels, edges, poly, Q::ONE);
            }
        }
    }
    if g > atn[l - 1].0 {
        let genera = head_genera(g - atn[l - 1].0 - 1);
        let extra = vec![n + 2 * l as u32, n + 2 * l as u32 + 1, n + 2 * l as u32 + 2];
        let labels = head_labels(&extra);
        let mut edges = head_edges(l as u32);
        edges.push((n + 2 * l as u32 + 1, n + 2 * l as u32 + 2));
        let nv = genera.len();
        let poly = kp_mul(
            &chain_psi_head(nv, n, &tvn, l),
            &psi_tail(nv, n + 2 * l as u32 + 1, n + 2 * l as u32 + 2, t1),
        );
        push(&mut parts, genera, labels, edges, poly, Q::ONE / Q::from(2));
    }
    TautClass::sum(parts)
}

// ===================== e^C / K / S 合成 =====================

/// 单顶点图上的纯 ψ^p 类（M̄_{g,n}）。
fn psi_only(g: u32, n: u32, leg: u32, p: u32, coeff: Q) -> TautClass {
    let legs: Vec<u32> = (1..=n).collect();
    let graph = StableGraph::new(&[g], &[legs], &[]).unwrap();
    TautClass::new(vec![DecStratum {
        graph,
        poly: KP::new(vec![KpTerm {
            kappa: vec![Vec::new()],
            psi: if p > 0 { vec![(leg, p)] } else { Vec::new() },
            coeff,
        }]),
    }])
}

/// 链图布局（exp/Kin/Sin 共用）：L+1 顶点、末顶点 extra 结点半边。
fn chain_parts(
    n: u32,
    g: u32,
    atn: &[(u32, Vec<u32>, i64)],
    last_extra: &[u32],
) -> (Vec<u32>, Vec<Vec<u32>>, Vec<(u32, u32)>) {
    let l = atn.len();
    let mut genera: Vec<u32> = vec![atn[0].0];
    for i in 1..l {
        genera.push(atn[i].0 - atn[i - 1].0);
    }
    genera.push(g - atn[l - 1].0);
    let scomp: Vec<u32> = (1..=n).filter(|x| !atn[l - 1].1.contains(x)).collect();
    let mut labels: Vec<Vec<u32>> = Vec::new();
    let mut first = atn[0].1.clone();
    first.push(n + 1);
    labels.push(first);
    for i in 1..l {
        let mut li: Vec<u32> = atn[i]
            .1
            .iter()
            .copied()
            .filter(|x| !atn[i - 1].1.contains(x))
            .collect();
        li.push(n + 2 * i as u32);
        li.push(n + 2 * i as u32 + 1);
        labels.push(li);
    }
    let mut last = scomp;
    for &x in last_extra {
        last.push(x);
    }
    labels.push(last);
    let edges: Vec<(u32, u32)> = (0..l as u32)
        .map(|i| (n + 2 * i + 1, n + 2 * i + 2))
        .collect();
    (genera, labels, edges)
}

/// 链的相邻关系：h 不减且 S 包含。
fn chain_step(x: &ATriple, y: &ATriple) -> bool {
    x.0 <= y.0 && x.1.iter().all(|s| y.1.contains(s))
}

/// generalized_hodge_chern_exp(t, a)：e^C 项（末 tvn > 1；最后段指数 tv−1−j−1）。
fn chern_exp(t: u32, a: &[ATriple], g: u32, n: u32) -> TautClass {
    let mut parts: Vec<TautClass> = Vec::new();
    for (idx, tvn) in chains_with_parts(a, t + 1) {
        let atn: Vec<ATriple> = idx.iter().map(|&i| a[i].clone()).collect();
        let l = atn.len();
        if tvn[l - 1] <= 1 {
            continue;
        }
        let (genera, labels, edges) = chain_parts(n, g, &atn, &[n + 2 * l as u32]);
        let nv = genera.len();
        let mut poly = chain_psi_head(nv, n, &tvn, l - 1);
        // a^{tv}(−1)^{tv−1} Σ_{j<tv−1} C(tv−1,j) ψ_{n+2L−1}^j ψ_{n+2L}^{tv−1−j−1}
        let lp = n + 2 * l as u32 - 1;
        let rp = n + 2 * l as u32;
        let tv = tvn[l - 1];
        let mut fin = KP::zero();
        for j in 0..tv - 1 {
            let mut term = kp_one(nv);
            if j > 0 {
                term = kp_mul(&term, &psi_pow(nv, lp, j));
            }
            if tv - 1 - j - 1 > 0 {
                term = kp_mul(&term, &psi_pow(nv, rp, tv - 1 - j - 1));
            }
            fin = fin.add(&term.scale(&binom(tv - 1, j)));
        }
        let fin = fin.scale(&(q_powi(atn[l - 1].2, tv) * q_powi(-1, tv - 1)));
        poly = kp_mul(&poly, &fin);
        if let Ok(graph) = StableGraph::new_raw(&genera, &labels, &edges) {
            parts.push(TautClass::new(vec![DecStratum { graph, poly }]).scale(&multinomial(&tvn)));
        }
    }
    TautClass::sum(parts).scale(&(Q::ONE / fact(t + 1)))
}

/// generalized_hodge_chern_Kin(alpha, beta, a)。
fn chern_kin(alpha: u32, beta: u32, a: &[ATriple], g: u32, n: u32) -> TautClass {
    if beta == 0 {
        return chern_exp(alpha - 1, a, g, n).scale(&fact(alpha));
    }
    if alpha == 0 {
        let legs: Vec<u32> = (1..=n).collect();
        let graph = StableGraph::new(&[g], &[legs], &[]).unwrap();
        return TautClass::new(vec![DecStratum {
            graph,
            poly: kp_kappa(0, beta - 1, 1, g, n),
        }]);
    }
    let mut parts: Vec<TautClass> = Vec::new();
    for (idx, tvn) in chains_with_parts(a, alpha) {
        let atn: Vec<ATriple> = idx.iter().map(|&i| a[i].clone()).collect();
        let l = atn.len();
        let scomp: Vec<u32> = (1..=n).filter(|x| !atn[l - 1].1.contains(x)).collect();
        let (genera, labels, edges) = chain_parts(n, g, &atn, &[n + 2 * l as u32]);
        let nv = genera.len();
        let mut poly = chain_psi_head(nv, n, &tvn, l);
        for (&tv, ai) in tvn.iter().zip(&atn) {
            poly = poly.scale(&q_powi(ai.2, tv));
        }
        let kap = kp_kappa(l, beta - 1, nv, g - atn[l - 1].0, scomp.len() as u32 + 1);
        poly = kp_mul(&poly, &kap);
        if let Ok(graph) = StableGraph::new_raw(&genera, &labels, &edges) {
            parts.push(TautClass::new(vec![DecStratum { graph, poly }]).scale(&multinomial(&tvn)));
        }
    }
    TautClass::sum(parts)
}

/// todd_inv_exp(a, t)：todd⁻¹ 的 e^D 展开（奇 t1 图组 + σ 项）。
pub fn todd_inv_exp(a: &[ATriple], t: u32, g: u32, n: u32) -> TautClass {
    let mut parts: Vec<TautClass> = Vec::new();
    let mut t1 = 1;
    while t1 < t {
        let coeff = divided_bernoulli(t1 + 1) / fact(t - t1);
        // 不成链的组合在 todd_inv_exp_graphs 中为零，只需枚举链；零分量对
        // Π a_i^{tv_i} 与多项式系数都贡献 1
        for (idx, tv) in chains_with_parts(a, t - t1) {
            let mut tvect = vec![0u32; a.len()];
            let mut c = coeff.clone();
            for (&i, &p) in idx.iter().zip(&tv) {
                tvect[i] = p;
                c *= q_powi(a[i].2, p);
            }
            let graphs = todd_inv_exp_graphs(a, t1, &tvect, g, n);
            parts.push(graphs.scale(&(c * multinomial(&tv))));
        }
        t1 += 2;
    }
    parts.push(todd_inv_sigma(t, g, n));
    TautClass::sum(parts)
}

/// generalized_hodge_chern_Sin(alpha, beta, p, a)。
fn chern_sin(alpha: u32, beta: u32, p: u32, a: &[ATriple], g: u32, n: u32) -> TautClass {
    if beta == 0 {
        return chern_exp(alpha - 1, a, g, n).scale(&fact(alpha));
    }
    if alpha == 0 {
        return psi_only(g, n, p, beta - 1, q_powi(-1, beta - 1));
    }
    let mut parts: Vec<TautClass> = Vec::new();
    for (idx, tvn) in chains_with_parts(a, alpha) {
        let atn: Vec<ATriple> = idx.iter().map(|&i| a[i].clone()).collect();
        let l = atn.len();
        let scomp: Vec<u32> = (1..=n).filter(|x| !atn[l - 1].1.contains(x)).collect();
        if !scomp.contains(&p) {
            continue;
        }
        let (genera, labels, edges) = chain_parts(n, g, &atn, &[n + 2 * l as u32]);
        let nv = genera.len();
        let mut poly = chain_psi_head(nv, n, &tvn, l);
        for (&tv, ai) in tvn.iter().zip(&atn) {
            poly = poly.scale(&q_powi(ai.2, tv));
        }
        if beta - 1 > 0 {
            poly = kp_mul(
                &poly,
                &psi_pow(nv, p, beta - 1).scale(&q_powi(-1, beta - 1)),
            );
        }
        if let Ok(graph) = StableGraph::new_raw(&genera, &labels, &edges) {
            parts.push(TautClass::new(vec![DecStratum { graph, poly }]).scale(&multinomial(&tvn)));
        }
    }
    TautClass::sum(parts)
}

/// generalized_hodge_chern_single(t, l, d, a)：ch_t = K + S + todd⁻¹e^D，
/// 度 t 截断。锚点：(1, 0, [1,−1], [(1,[],−1)], 1, 2) =
/// 1/12κ₁−13/12ψ₁−1/12ψ₂ / loop 1/24 / 分离 −11/12。
pub fn generalized_hodge_chern_single(
    t: u32,
    l: i64,
    d: &[i64],
    a: &[(u32, Vec<u32>, i64)],
    g: u32,
    n: u32,
) -> Result<TautClass, String> {
    Ok(merge_isomorphic(&hodge_chern_unmerged(t, l, d, a, g, n)?))
}

/// 同构项合并前的 ch_t（系数为 oracle 生成元口径，未乘 |Aut|）。
fn hodge_chern_unmerged(
    t: u32,
    l: i64,
    d: &[i64],
    a: &[(u32, Vec<u32>, i64)],
    g: u32,
    n: u32,
) -> Result<TautClass, String> {
    if n == 0 {
        return Err("generalized_hodge_chern 需 n ≥ 1".into());
    }
    if d.len() != n as usize {
        return Err("d 长度须为 n".into());
    }
    let mut asort: Vec<ATriple> = a.to_vec();
    asort.sort_by(|x, y| x.1.cmp(&y.1).then(x.2.cmp(&y.2)));
    let mut parts: Vec<TautClass> = Vec::new();
    for alpha in 0..=t + 1 {
        let beta = t + 1 - alpha;
        parts.push(
            chern_kin(alpha, beta, &asort, g, n)
                .scale(&(divided_bernoulli_poly(l, beta) / fact(alpha))),
        );
    }
    for alpha in 0..=t {
        let beta = t + 1 - alpha;
        for p in 1..=n {
            for b in 0..=beta {
                let c = q_powi(-1, b)
                    * q_powi(d[(p - 1) as usize], beta - b)
                    * divided_bernoulli_poly(l, b)
                    / (fact(alpha) * fact(beta - b));
                parts.push(chern_sin(alpha, beta, p, &asort, g, n).scale(&c));
            }
        }
    }
    parts.push(todd_inv_exp(&asort, t, g, n));
    Ok(degree_part(&TautClass::sum(parts), t))
}

/// twist_div_to_zero(l, d)：使 twist 除子 D = l·K̃ + Σd_iσ_i 拉回到零的
/// a 表：a_{h,S} = −l(2h−2+|S|+1) − Σ_{p∈S} d_p。
/// 锚点：twist_div_to_zero(3,[3,4],2,2) = [[0,[1,2],−10],[1,[1],−9],[1,[1,2],−16]]。
pub fn twist_div_to_zero(l: i64, d: &[i64], g: u32, n: u32) -> Vec<ATriple> {
    let mut out = Vec::new();
    for h in 0..=g {
        for s in powerset_sorted(&(2..=n).collect::<Vec<u32>>()) {
            if (h == 0 && s.is_empty()) || (h == g && s.len() as i64 > n as i64 - 3) {
                continue;
            }
            let mut sl: Vec<u32> = vec![1];
            sl.extend(s.iter().copied());
            let ah = -l * (2 * h as i64 - 2 + sl.len() as i64 + 1)
                - sl.iter().map(|&p| d[(p - 1) as usize]).sum::<i64>();
            if ah != 0 {
                out.push((h, sl, ah));
            }
        }
    }
    out
}

/// 同构合并：按装饰单项式的完整同构证书（[`StableGraph::decorated_certificate`]）
/// 聚合——链图 new_raw 与分离图 new 构造的同构项合并为一项（oracle 的类展示
/// 即合并后形式），代表项取首次出现者。
///
/// 原先的键是"顶点色（亏格、标记、κ）多重集 + 边两端色与 ψ 多重集"，不是
/// 完整的同构不变量（如四个同色顶点连成 4-圈与"三角形 + 悬挂点"键相同），
/// 理论上会把不同构的项合并。
fn merge_isomorphic(cl: &TautClass) -> TautClass {
    use std::collections::BTreeMap;
    // 每项生成 (证书, 系数和, 代表项)
    let mut acc: BTreeMap<Vec<u64>, (DecStratum, Q)> = BTreeMap::new();
    for ds in &cl.terms {
        for term in &ds.poly.terms {
            let key = ds.graph.decorated_certificate(&term.kappa, &term.psi);
            let e = acc.entry(key).or_insert_with(|| {
                (
                    DecStratum {
                        graph: ds.graph.clone(),
                        poly: KP::new(vec![term.clone()]),
                    },
                    Q::ZERO,
                )
            });
            e.1 = e.1.clone() + term.coeff.clone();
        }
    }
    // 口径：push 系数 = aut × oracle 展示系数（evaluate 统一除 aut——
    // Θ/Witten/drspin 四链路同约定；grr 的 doctest 系数为 oracle 展示值）
    let mut terms: Vec<DecStratum> = Vec::new();
    for (_, (mut ds, c)) in acc {
        if c != Q::ZERO {
            let aut = Q::from(ds.graph.automorphism_number() as i64);
            for t in ds.poly.terms.iter_mut() {
                t.coeff = c.clone() * aut.clone();
            }
            ds.poly = KP::new(ds.poly.terms.clone());
            terms.push(ds);
        }
    }
    TautClass::new(terms)
}

/// generalized_hodge_chern(l, d, a, dmax, g, n)：Σ_{t=1..dmax} ch_t。
pub fn generalized_hodge_chern(
    l: i64,
    d: &[i64],
    a: &[(u32, Vec<u32>, i64)],
    dmax: u32,
    g: u32,
    n: u32,
) -> Result<TautClass, String> {
    let mut parts: Vec<TautClass> = Vec::new();
    for t in 1..=dmax {
        parts.push(generalized_hodge_chern_single(t, l, d, a, g, n)?);
    }
    Ok(TautClass::sum(parts))
}

#[cfg(test)]
mod merge_tests {
    use super::*;

    /// 回归：不同构但颜色多重集相同的两个图不被合并（旧键会把 4-圈与
    /// "三角形 + 悬挂点" 合为一项）；同构的两种写法合并、系数相加。
    #[test]
    fn merge_isomorphic_uses_complete_invariant() {
        let one = |graph: StableGraph, c: i64| DecStratum {
            graph,
            poly: KP::new(vec![KpTerm {
                kappa: vec![vec![]; 4],
                psi: vec![],
                coeff: Q::from(c),
            }]),
        };
        let cycle = |legs: [[u32; 2]; 4], edges: [(u32, u32); 4]| {
            StableGraph::new_raw(&[1, 1, 1, 1], &legs.map(|l| l.to_vec()), &edges).unwrap()
        };
        let c1 = cycle(
            [[1, 2], [3, 4], [5, 6], [7, 8]],
            [(2, 3), (4, 5), (6, 7), (8, 1)],
        );
        let c2 = cycle(
            [[5, 2], [3, 8], [1, 6], [7, 4]],
            [(2, 3), (8, 1), (6, 7), (4, 5)],
        );
        let tadpole = StableGraph::new_raw(
            &[1, 1, 1, 1],
            &[vec![1, 2, 3], vec![4, 5], vec![6, 7], vec![8]],
            &[(1, 4), (5, 6), (7, 2), (3, 8)],
        )
        .unwrap();
        let merged = merge_isomorphic(&TautClass::new(vec![
            one(c1.clone(), 1),
            one(tadpole.clone(), 1),
            one(c2, 2),
        ]));
        assert_eq!(merged.terms.len(), 2);
        // 系数口径：合并后的系数和 × 代表图的 |Aut|（4-圈为二面体群，8 阶）
        assert_eq!(c1.automorphism_number(), 8);
        for t in &merged.terms {
            let expect = if t.graph == tadpole {
                tadpole.automorphism_number() as i64
            } else {
                3 * 8
            };
            assert_eq!(t.poly.terms[0].coeff, Q::from(expect), "{}", t.graph);
        }
    }
}

#[cfg(test)]
mod merge_consistency {
    use super::*;
    use crate::classes::int_against;

    /// 合并同构项不改变类：未合并的项各自换算到 ξ_*/|Aut| 口径（乘其图的
    /// |Aut|）后，与合并后的类对一组 X 的积分逐个相等。旧颜色键在此例把链图
    /// 中间顶点与末端椭圆尾上的 ψ²（系数 ±1/180）误并为一项而相消。
    #[test]
    fn merge_preserves_integrals() {
        let (t, l, d, g, n) = (4u32, 1i64, vec![1i64, 1, -2], 3u32, 3u32);
        let a = twist_div_to_zero(l, &d, g, n);
        let raw = hodge_chern_unmerged(t, l, &d, &a, g, n).unwrap();
        let scaled = TautClass::new(
            raw.terms
                .iter()
                .map(|ds| DecStratum {
                    graph: ds.graph.clone(),
                    poly: ds
                        .poly
                        .scale(&Q::from(ds.graph.automorphism_number() as i64)),
                })
                .collect(),
        );
        let merged = merge_isomorphic(&raw);
        // X：余维 3g−3+n−t = 5 的 κ-ψ 单项式
        let mut xs: Vec<(Vec<(u32, u32)>, Vec<(u32, u32)>)> = Vec::new();
        for p1 in 0..=5u32 {
            for p2 in 0..=5 - p1 {
                let p3 = 5 - p1 - p2;
                let psi: Vec<(u32, u32)> = [(1, p1), (2, p2), (3, p3)]
                    .into_iter()
                    .filter(|&(_, p)| p > 0)
                    .collect();
                xs.push((vec![], psi));
            }
        }
        for kappa in [
            vec![(1, 5)],
            vec![(1, 3), (2, 1)],
            vec![(5, 1)],
            vec![(2, 2), (1, 1)],
        ] {
            xs.push((kappa, vec![]));
        }
        xs.push((vec![(1, 2)], vec![(1, 2), (3, 1)]));
        xs.push((vec![(3, 1)], vec![(2, 2)]));
        let mut nonzero = 0;
        for (kappa, psi) in &xs {
            let lhs = int_against(&scaled, kappa, psi).unwrap();
            let rhs = int_against(&merged, kappa, psi).unwrap();
            assert_eq!(lhs, rhs, "κ={kappa:?} ψ={psi:?}");
            if lhs != Q::ZERO {
                nonzero += 1;
            }
        }
        assert!(nonzero > 5, "非零积分过少：{nonzero}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// chern_single(1, 0, [1,−1], [(1,[],−1)], 1, 2)：oracle doctest 三图。
    #[test]
    fn chern_single_anchor() {
        let cl = generalized_hodge_chern_single(1, 0, &[1, -1], &[(1, vec![], -1)], 1, 2).unwrap();
        let mut rows: Vec<(Vec<u32>, Vec<Vec<u32>>, String)> = Vec::new();
        for t in &cl.terms {
            let gs: Vec<u32> = t.graph.vertices().iter().map(|v| v.genus).collect();
            let ls: Vec<Vec<u32>> = t
                .graph
                .vertices()
                .iter()
                .map(|v| {
                    let mut l: Vec<u32> = v.legs.clone().into_iter().filter(|x| *x <= 2).collect();
                    l.sort_unstable();
                    l
                })
                .collect();
            for term in &t.poly.terms {
                let kappa: String = term
                    .kappa
                    .iter()
                    .enumerate()
                    .flat_map(|(vi, k)| k.iter().map(move |&j| (vi, j)))
                    .map(|(vi, j)| format!("k{}_v{}", j, vi))
                    .collect::<Vec<_>>()
                    .join(",");
                let psi: String = term
                    .psi
                    .iter()
                    .map(|&(l, p)| format!("p{}^{}", l, p))
                    .collect::<Vec<_>>()
                    .join(",");
                rows.push((
                    gs.clone(),
                    ls.clone(),
                    format!("{}*({};{})", term.coeff, kappa, psi),
                ));
            }
        }
        let get = |pred: &dyn Fn(&(Vec<u32>, Vec<Vec<u32>>, String)) -> bool| -> String {
            rows.iter()
                .filter(|r| pred(r))
                .map(|(_, _, s)| s.clone())
                .next()
                .unwrap_or_default()
        };
        // 主图 [1]：1/12κ₁、−13/12ψ₁、−1/12ψ₂
        assert_eq!(
            get(&|r| r.0 == vec![1] && r.1 == vec![vec![1, 2]] && r.2.contains("k1_v0")),
            "1/12*(k1_v0;)",
            "{rows:?}"
        );
        assert_eq!(
            get(&|r| r.0 == vec![1] && r.1 == vec![vec![1, 2]] && r.2.contains("p1^1")),
            "-13/12*(;p1^1)",
            "{rows:?}"
        );
        assert_eq!(
            get(&|r| r.0 == vec![1] && r.1 == vec![vec![1, 2]] && r.2.contains("p2^1")),
            "-1/12*(;p2^1)",
            "{rows:?}"
        );
        // loop 图 [0]（自环，aut=2）：类系数 = 2×1/24；分离图 {0,1}：−11/12
        assert_eq!(
            get(&|r| r.0 == vec![0] && r.2 == "1/12*(;)"),
            "1/12*(;)",
            "{rows:?}"
        );
        assert_eq!(
            get(&|r| r.0 == vec![0, 1] && r.2 == "-11/12*(;)"),
            "-11/12*(;)",
            "{rows:?}"
        );
    }

    /// λ₁ = generalized_hodge_chern_single(1, 1, [0], [], 2, 1)
    /// （oracle generalized_lambda(1,1,[0],[],2,1) doctest）：
    /// 1/12κ₁−1/12ψ₁ / loop 2×1/24（aut 口径）/ 分离 1/12。
    #[test]
    fn lambda1_anchor() {
        let cl = generalized_hodge_chern_single(1, 1, &[0], &[], 2, 1).unwrap();
        let rows: Vec<(Vec<u32>, String)> = cl
            .terms
            .iter()
            .flat_map(|t| {
                let gs: Vec<u32> = t.graph.vertices().iter().map(|v| v.genus).collect();
                t.poly.terms.iter().map(move |term| {
                    let kappa: Vec<String> = term
                        .kappa
                        .iter()
                        .flatten()
                        .map(|&j| format!("k{}", j))
                        .collect();
                    let psi: Vec<String> = term
                        .psi
                        .iter()
                        .map(|&(l, p)| format!("p{}^{}", l, p))
                        .collect();
                    (
                        gs.clone(),
                        format!("{}*({};{})", term.coeff, kappa.join(","), psi.join(",")),
                    )
                })
            })
            .collect();
        assert!(
            rows.contains(&(vec![2u32], "1/12*(k1;)".to_string())),
            "{rows:?}"
        );
        assert!(
            rows.contains(&(vec![2u32], "-1/12*(;p1^1)".to_string())),
            "{rows:?}"
        );
        assert!(
            rows.contains(&(vec![1u32], "1/12*(;)".to_string())),
            "{rows:?}"
        );
        assert!(
            rows.contains(&(vec![1u32, 1u32], "1/12*(;)".to_string())),
            "{rows:?}"
        );
    }

    /// GRRcomp 积分差分摘要（tools/grr_diff*.py，18 例全对）强例。
    #[test]
    fn grr_int_oracle() {
        use crate::classes::int_against;
        let ev = |t: u32,
                  l: i64,
                  d: &[i64],
                  a: &[(u32, Vec<u32>, i64)],
                  g: u32,
                  n: u32,
                  kap: &[(u32, u32)],
                  psi: &[(u32, u32)]| {
            let cl = generalized_hodge_chern_single(t, l, d, a, g, n).unwrap();
            int_against(&cl, kap, psi).unwrap().to_string()
        };
        assert_eq!(
            ev(1, 0, &[1, -1], &[(1, vec![], -1)], 1, 2, &[(1, 1)], &[]),
            "-1/12"
        );
        assert_eq!(ev(1, 1, &[1, 0], &[], 1, 2, &[], &[(1, 1)]), "1/24");
        assert_eq!(
            ev(
                1,
                0,
                &[2, -1],
                &[(1, vec![1], -1), (0, vec![2], 1)],
                1,
                2,
                &[],
                &[(2, 1)]
            ),
            "-1/12"
        );
        assert_eq!(ev(1, 1, &[0], &[], 2, 1, &[(1, 1)], &[(1, 2)]), "5/576");
        assert_eq!(ev(1, 1, &[0], &[], 2, 1, &[(3, 1)], &[]), "1/480");
        assert_eq!(
            ev(
                1,
                0,
                &[1, 1, -2],
                &[(1, vec![1], -1)],
                1,
                3,
                &[],
                &[(1, 1), (2, 1)]
            ),
            "-1/6"
        );
        // twist_div_to_zero doctest
        let a = twist_div_to_zero(3, &[3, 4], 2, 2);
        assert_eq!(
            a,
            vec![(0, vec![1, 2], -10), (1, vec![1], -9), (1, vec![1, 2], -16),]
        );
    }

    #[test]
    fn basics() {
        assert_eq!(divided_bernoulli(12).to_string(), "-691/1307674368000");
        assert_eq!(divided_bernoulli_poly(0, 4).to_string(), "-1/720");
        assert_eq!(divided_bernoulli_poly(2, 4).to_string(), "119/720");
        assert_eq!(divided_bernoulli_poly(1, 2).to_string(), "1/12");
        assert_eq!(divided_bernoulli_poly(-1, 3).to_string(), "-1/2");
        assert_eq!(divided_bernoulli_poly(3, 1).to_string(), "5/2");
    }

    /// 回归：直接枚举链 == 枚举 total 分进 |a| 份的全部弱组合后按链条件过滤。
    #[test]
    fn chains_match_composition_filter() {
        fn comps(t: u32, k: usize) -> Vec<Vec<u32>> {
            if k == 0 {
                return if t == 0 { vec![vec![]] } else { vec![] };
            }
            let mut out = Vec::new();
            for e in 0..=t {
                for mut rest in comps(t - e, k - 1) {
                    rest.insert(0, e);
                    out.push(rest);
                }
            }
            out
        }
        for (l, d, g, n) in [
            (1i64, vec![1i64, -1], 1u32, 2u32),
            (1, vec![2, 1, -3], 2, 3),
            (3, vec![3, 4], 2, 2),
            (1, vec![1, 1, 1, -3], 2, 4),
        ] {
            let mut a = twist_div_to_zero(l, &d, g, n);
            a.sort_by(|x, y| x.1.cmp(&y.1).then(x.2.cmp(&y.2)));
            for total in 0..=4 {
                let mut brute: Vec<(Vec<usize>, Vec<u32>)> = comps(total, a.len())
                    .into_iter()
                    .filter_map(|tv| {
                        let idx: Vec<usize> = (0..a.len()).filter(|&i| tv[i] > 0).collect();
                        let ok = !idx.is_empty()
                            && idx.windows(2).all(|w| chain_step(&a[w[0]], &a[w[1]]));
                        ok.then(|| (idx.clone(), idx.iter().map(|&i| tv[i]).collect()))
                    })
                    .collect();
                let mut fast = chains_with_parts(&a, total);
                brute.sort();
                fast.sort();
                assert_eq!(fast, brute, "a={a:?} total={total}");
            }
        }
    }

    /// todd_inv_exp([],1,1,2)：loop 1/24 + 分离 1/12（oracle doctest）。
    #[test]
    fn todd_sigma_anchor() {
        let cl = todd_inv_exp(&[], 1, 1, 2);
        // 收集 (nv, ne, genera, coeff)
        let mut got: Vec<(usize, usize, Vec<u32>, String)> = Vec::new();
        for t in &cl.terms {
            for term in &t.poly.terms {
                got.push((
                    t.graph.vertices().len(),
                    t.graph.edges().len(),
                    t.graph.vertices().iter().map(|v| v.genus).collect(),
                    term.coeff.to_string(),
                ));
            }
        }
        assert_eq!(
            got.iter()
                .find(|(nv, ne, gs, _)| *nv == 1 && *ne == 1 && gs == &vec![0u32])
                .map(|(_, _, _, c)| c.clone()),
            Some("1/24".into()),
            "{got:?}"
        );
        assert_eq!(
            got.iter()
                .find(|(nv, ne, gs, _)| *nv == 2 && *ne == 1 && gs == &vec![0u32, 1u32])
                .map(|(_, _, _, c)| c.clone()),
            Some("1/12".into()),
            "{got:?}"
        );
    }

    /// todd_inv_exp_graphs([[1,[1],−1]],1,[1],3,2) → 5 图（oracle 实测）。
    #[test]
    fn todd_graphs_anchor() {
        let cl = todd_inv_exp_graphs(&[(1, vec![1], -1)], 1, &[1], 3, 2);
        assert_eq!(cl.terms.len(), 5, "应有 5 个图");
        // 图 1：[1,2] [[1,3],[2,4]] (3,4) −ψ₃
        // 图 2：[1,1] [[1,3],[2,4,5,6]] (3,4),(5,6) ½
        // 图 3：[1,0,2] [[1,3],[2,4,5],[6]] 1
        // 图 4：[1,1,1] [[1,3],[2,4,5],[6]] 1
        // 图 5：[1,1,1] [[1,3],[4,5],[2,6]] 1
        let desc: Vec<(Vec<u32>, Vec<Vec<u32>>, Vec<(u32, u32)>, String)> = cl
            .terms
            .iter()
            .map(|t| {
                (
                    t.graph.vertices().iter().map(|v| v.genus).collect(),
                    t.graph
                        .vertices()
                        .iter()
                        .map(|v| {
                            let mut l = v.legs.clone();
                            l.sort_unstable();
                            l
                        })
                        .collect(),
                    {
                        let mut e = t.graph.edges().to_vec();
                        e.sort_unstable();
                        e
                    },
                    t.poly
                        .terms
                        .iter()
                        .map(|x| x.coeff.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                )
            })
            .collect();
        let has = |gs: Vec<u32>, ls: Vec<Vec<u32>>, es: Vec<(u32, u32)>, c: &str| {
            desc.iter()
                .any(|(g, l, e, cc)| *g == gs && *l == ls && *e == es && cc == c)
        };
        assert!(
            has(vec![1, 2], vec![vec![1, 3], vec![2, 4]], vec![(3, 4)], "-1"),
            "{desc:?}"
        );
        assert!(
            has(
                vec![1, 1],
                vec![vec![1, 3], vec![2, 4, 5, 6]],
                vec![(3, 4), (5, 6)],
                "1/2"
            ),
            "{desc:?}"
        );
        assert!(
            has(
                vec![1, 0, 2],
                vec![vec![1, 3], vec![2, 4, 5], vec![6]],
                vec![(3, 4), (5, 6)],
                "1"
            ),
            "{desc:?}"
        );
        assert!(
            has(
                vec![1, 1, 1],
                vec![vec![1, 3], vec![2, 4, 5], vec![6]],
                vec![(3, 4), (5, 6)],
                "1"
            ),
            "{desc:?}"
        );
        assert!(
            has(
                vec![1, 1, 1],
                vec![vec![1, 3], vec![4, 5], vec![2, 6]],
                vec![(3, 4), (5, 6)],
                "1"
            ),
            "{desc:?}"
        );
    }
}
