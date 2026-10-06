//! Witten r-spin 极限类（PPZ16 附录；admcycles/witten.py 移植）。
//!
//! P_d(X,a)（PPZ16 §4.5 递推）实现为 a-多项式（系数 ∈ QQ[X]）；符号求和用
//! Bernoulli 闭式 Σ_{b=1..N} b^k = [B_{k+1}(N+1) − B_{k+1}(1)]/(k+1)
//! （负 N 解析延拓——锚点 PM(1,6,−2) = −53/6 证实）。边因子 rspin_edge_factor：
//! S = Σ c_i X^i（c_i = leg 因子在 X=m 处的值），除以 X+1（综合除法）取
//! −T[d1]。权重和 mod (m−1)；插值点步长 step = lcm(C_i 分母)；P(0) 拉格朗日；
//! 整体 ×(−1)^{d−g}。

use crate::generators::tautgens;
use crate::kppoly::{KappaPsiPolynomial as KP, KpTerm};
use crate::socle::bernoulli;
use crate::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero};

type PX = Vec<Q>; // X 多项式（升幂，可含尾部零）

fn px_add(p: &PX, q: &PX) -> PX {
    let mut r = p.clone();
    r.resize(r.len().max(q.len()), Q::ZERO);
    for (i, x) in q.iter().enumerate() {
        r[i] = r[i].clone() + x.clone();
    }
    r
}

fn px_mul(p: &PX, q: &PX) -> PX {
    if p.is_empty() || q.is_empty() {
        return Vec::new();
    }
    let mut r = vec![Q::ZERO; p.len() + q.len() - 1];
    for (i, a) in p.iter().enumerate() {
        if *a == Q::ZERO {
            continue;
        }
        for (j, b) in q.iter().enumerate() {
            r[i + j] = r[i + j].clone() + a.clone() * b.clone();
        }
    }
    r
}

fn px_eval(p: &PX, x: &Q) -> Q {
    let mut v = Q::ZERO;
    for c in p.iter().rev() {
        v = v * x + c;
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

/// Bernoulli 多项式 B_n(x) 升幂系数：B_n(x) = Σ_k C(n,k) B_k x^{n−k}。
/// 约定 B_1 = −1/2（与 sage/admcycles 一致；测试锚定）。
fn bern_poly_coeffs(n: u32) -> PX {
    let mut c = vec![Q::ZERO; n as usize + 1];
    for k in 0..=n {
        c[(n - k) as usize] = binom(n, k) * bernoulli(k as u64);
    }
    c
}

/// Σ_{b=1..N} b^k 的 N-多项式（升幂，度 k+1）。
fn sum_pow_coeffs(k: u32) -> PX {
    let b = bern_poly_coeffs(k + 1);
    // B_{k+1}(N+1)：N 幂系数
    let mut shifted = vec![Q::ZERO; b.len()];
    for (j, c) in b.iter().enumerate() {
        if *c == Q::ZERO {
            continue;
        }
        for (t, s) in shifted.iter_mut().enumerate().take(j + 1) {
            *s += c.clone() * binom(j as u32, t as u32);
        }
    }
    // B_{k+1}(1) = Σ b_j
    let b1: Q = b.iter().fold(Q::ZERO, |acc, c| acc + c.clone());
    let mut out = shifted;
    out[0] = out[0].clone() - b1;
    for s in out.iter_mut() {
        *s = s.clone() / Q::from((k + 1) as i64);
    }
    out
}

/// P_d(X,a)：a-多项式（外层升幂），系数 ∈ QQ[X]。
/// [`pm_poly_calc`] 的记忆化：P_d 只依赖 d，而每层递推都是整套多项式运算；
/// leg_factor/edge_factor 在权重枚举中被反复调用，不缓存时每次都从 P_0 重建。
fn pm_poly(d: u32) -> std::rc::Rc<Vec<PX>> {
    thread_local! {
        static MEMO: std::cell::RefCell<std::collections::HashMap<u32, std::rc::Rc<Vec<PX>>>> =
            std::cell::RefCell::new(std::collections::HashMap::new());
    }
    if let Some(p) = MEMO.with(|m| m.borrow().get(&d).cloned()) {
        return p;
    }
    let p = std::rc::Rc::new(pm_poly_calc(d));
    MEMO.with(|m| m.borrow_mut().insert(d, p.clone()));
    p
}

fn pm_poly_calc(d: u32) -> Vec<PX> {
    if d == 0 {
        return vec![vec![Q::ONE]];
    }
    let prev = pm_poly(d - 1);
    // P̃(b) = P_{d−1}(X, b−1)：b^k 系数 = Σ_{i≥k} prev[i]·C(i,k)·(−1)^{i−k}
    let mut pb: Vec<PX> = vec![Vec::new(); prev.len()];
    for (i, coef) in prev.iter().enumerate() {
        for (k, pbk) in pb.iter_mut().enumerate().take(i + 1) {
            let sgn = Q::from(if (i - k) % 2 == 0 { 1 } else { -1 });
            let term: PX = coef
                .iter()
                .map(|c| c.clone() * binom(i as u32, k as u32) * sgn.clone())
                .collect();
            *pbk = px_add(pbk, &term);
        }
    }
    // q(b) = (2d−1)X·P̃(b) − 2b·P̃(b)：q_k = L·pb_k − 2·pb_{k−1}
    let lin: PX = vec![Q::ZERO, Q::from(2 * d as i64 - 1)];
    let mut q: Vec<PX> = vec![Vec::new(); pb.len() + 1];
    for k in 0..=pb.len() {
        if let Some(coef) = pb.get(k) {
            let t = px_mul(&lin, coef);
            q[k] = px_add(&q[k], &t);
        }
        if k > 0 {
            let t2: PX = pb[k - 1].iter().map(|c| c.clone() * Q::from(-2)).collect();
            q[k] = px_add(&q[k], &t2);
        }
    }
    // T1 = ½ Σ_k q_k·S_k(a)
    let mut t1: Vec<PX> = Vec::new();
    for (k, coef) in q.iter().enumerate() {
        if coef.is_empty() {
            continue;
        }
        let sk = sum_pow_coeffs(k as u32);
        for (ai, sc) in sk.iter().enumerate() {
            let term: PX = coef
                .iter()
                .map(|c| c.clone() * sc.clone() * Q::from(1) / Q::from(2))
                .collect();
            while t1.len() <= ai {
                t1.push(Vec::new());
            }
            t1[ai] = px_add(&t1[ai], &term);
        }
    }
    // 三次式 C3(b) = (X−1−b)(2dX−b)(2dX−X−2b) 的 b^t 系数（QQ[X]）：
    // 令 u = −b：(u+A)(u+B)(2u+C)，A = X−1，B = 2dX，C = (2d−1)X，
    // 逐括号乘（每个括号 = c_lin·u + c_const，系数 ∈ QQ[X]）
    let mut c3: Vec<PX> = vec![vec![Q::ONE]];
    for (clin, cconst) in [
        (vec![Q::ONE], vec![Q::from(-1), Q::ONE]), // u + X − 1
        (vec![Q::ONE], vec![Q::ZERO, Q::from(2 * d as i64)]), // u + 2dX
        (vec![Q::from(2)], vec![Q::ZERO, Q::from(2 * d as i64 - 1)]), // 2u + (2d−1)X
    ] {
        let mut nxt: Vec<PX> = vec![Vec::new(); c3.len() + 1];
        for (t, coef) in c3.iter().enumerate() {
            if !coef.is_empty() {
                nxt[t] = px_add(&nxt[t], &px_mul(coef, &cconst));
                nxt[t + 1] = px_add(&nxt[t + 1], &px_mul(coef, &clin));
            }
        }
        c3 = nxt;
    }
    // r(b) = C3(b)·P̃(b)：u = −b，故 b^t 系数 = (−1)^t·c3_u[t]；
    // r_k = Σ_t c3_b[t]·pb[k−t]
    let mut r: Vec<PX> = vec![Vec::new(); pb.len() + 4];
    for (t, cf) in c3.iter().enumerate() {
        let sgn = Q::from(if t % 2 == 0 { 1 } else { -1 });
        let cb: PX = cf.iter().map(|c| c.clone() * sgn.clone()).collect();
        for (k, coef) in pb.iter().enumerate() {
            if cf.is_empty() || coef.is_empty() {
                continue;
            }
            let term = px_mul(&cb, coef);
            r[k + t] = px_add(&r[k + t], &term);
        }
    }
    // U = Σ_k r_k·S_k(X−2)：S_k(N) 多项式代 N = X−2
    let mut u: PX = Vec::new();
    for (k, coef) in r.iter().enumerate() {
        if coef.is_empty() {
            continue;
        }
        let sk = sum_pow_coeffs(k as u32);
        for (ai, sc) in sk.iter().enumerate() {
            let mut shift: PX = vec![Q::ONE];
            for _ in 0..ai {
                shift = px_mul(&shift, &vec![Q::from(-2), Q::ONE]);
            }
            let term = px_mul(
                &coef.iter().map(|c| c.clone() * sc.clone()).collect::<PX>(),
                &shift,
            );
            u = px_add(&u, &term);
        }
    }
    // T2 = U / (4dX(X−1)) = U / (4dX² − 4dX)：长除断言整除
    let den: PX = vec![Q::ZERO, Q::from(-(4 * d as i64)), Q::from(4 * d as i64)];
    let t2 = px_div_exact(&u, &den).expect("P_d 第二项整除失败");
    // P_d：a^{j>0} 系数 = T1[j]；a⁰ = T1[0] − T2
    let mut out = t1;
    if out.is_empty() {
        out.push(Vec::new());
    }
    out[0] = px_add(
        &out[0],
        &t2.iter()
            .map(|c| Q::ZERO.clone() - c.clone())
            .collect::<PX>(),
    );
    out
}

/// 多项式长除（断言整除，返回商）。
fn px_div_exact(num: &PX, den: &PX) -> Option<PX> {
    let dn = den.len() as i64 - 1;
    let lead = den.last()?;
    let mut rem = num.clone();
    while !rem.iter().all(|c| *c == Q::ZERO) {
        let rn = rem.len() as i64 - 1;
        if rn < dn {
            return None;
        }
        let fac = rem.last().unwrap().clone() / lead.clone();
        let shift = (rn - dn) as usize;
        for (i, dc) in den.iter().enumerate() {
            rem[shift + i] = rem[shift + i].clone() - dc.clone() * fac.clone();
        }
        while rem.len() > 1 && *rem.last().unwrap() == Q::ZERO {
            rem.pop();
        }
    }
    // 商的重构：简单法——重做一遍记录因子
    let mut rem2 = num.clone();
    let mut quo = vec![Q::ZERO; num.len().saturating_sub(den.len()) + 1];
    while rem2.len() > 1 || rem2[0] != Q::ZERO {
        let rn = rem2.len() as i64 - 1;
        if rn < dn {
            if rem2.iter().all(|c| *c == Q::ZERO) {
                break;
            }
            return None;
        }
        let fac = rem2.last().unwrap().clone() / lead.clone();
        let shift = (rn - dn) as usize;
        quo[shift] = fac.clone();
        for (i, dc) in den.iter().enumerate() {
            rem2[shift + i] = rem2[shift + i].clone() - dc.clone() * fac.clone();
        }
        while rem2.len() > 1 && *rem2.last().unwrap() == Q::ZERO {
            rem2.pop();
        }
    }
    Some(quo)
}

/// P_d(X,·) 在数值 (x, a) 处的值。
pub fn pm_eval(d: u32, x: i64, a: i64) -> Q {
    let tbl = pm_poly(d);
    let xv = Q::from(x);
    let av = Q::from(a);
    let mut acc = Q::ZERO;
    let mut ap = Q::ONE;
    for coef in tbl.iter() {
        if !coef.is_empty() {
            acc += px_eval(coef, &xv) * ap.clone();
        }
        ap *= av.clone();
    }
    acc
}

/// rspin_leg_factor(d, a, x)：a < 0 时 a → X+a（数值 x+a）。
pub fn leg_factor(d: u32, a: i64, x: i64) -> Q {
    let aa = if a < 0 { x + a } else { a };
    pm_eval(d, x, aa)
}

/// rspin_edge_factor(w1, m, d1, d2)。
pub fn edge_factor(w1: i64, m: i64, d1: u32, d2: u32) -> Q {
    let d = d1 + d2 + 1;
    let mm = m - 1;
    let mut c = vec![Q::ZERO; d as usize + 1];
    for i in 0..=d {
        let a1 = (w1 + i as i64).rem_euclid(mm);
        let a2 = (m - 2 - i as i64 - w1).rem_euclid(mm);
        c[i as usize] = leg_factor(i, a1, m) * leg_factor(d - i, a2, m);
    }
    // T = S/(X+1)：T_{j−1} = c_j − T_j（从高位降），取 −T[d1]
    let mut t = Q::ZERO;
    for j in ((d1 as usize + 1)..=d as usize).rev() {
        t = c[j].clone() - t;
    }
    Q::ZERO - t
}

/// 边因子表 [edge_factor(x, m, d1, d2)]_{x=0..m−2}（按 (d1, d2, m) 跨生成元缓存）。
fn edge_table(d1: u32, d2: u32, m: i64) -> std::rc::Rc<Vec<Q>> {
    type Key = (u32, u32, i64);
    thread_local! {
        static MEMO: std::cell::RefCell<std::collections::HashMap<Key, std::rc::Rc<Vec<Q>>>> =
            std::cell::RefCell::new(std::collections::HashMap::new());
    }
    if let Some(t) = MEMO.with(|c| c.borrow().get(&(d1, d2, m)).cloned()) {
        return t;
    }
    let t: std::rc::Rc<Vec<Q>> =
        std::rc::Rc::new((0..m - 1).map(|x| edge_factor(x, m, d1, d2)).collect());
    MEMO.with(|c| c.borrow_mut().insert((d1, d2, m), t.clone()));
    t
}

fn num_lcm(a: i64, b: i64) -> i64 {
    let (a, b) = (a.abs(), b.abs());
    let g = num_gcd(a, b);
    if g == 0 { 0 } else { a / g * b }
}

fn num_gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { num_gcd(b, a % b) }
}

/// Witten r-spin 极限类：A_i = C_i·r + D_i（(C_num, C_den, D)，C ∈ Q、D ∈ Z；
/// 整数 A_i 记 (0,1,A_i)）。返回 r^{g−1+C}·W_{g,n}^r(A) 的 r→∞ 极限（r=0 值）。
///
/// 移植自 admcycles Wittenrspin/rspin_coeff（PPZ16 附录）：给定权
/// w_v = 1 − Σj·e_{κ_j,v} + Σ_{i∈v}(A_i − e_{ψ,i})（无 genus 项），
/// 边枚举 x ∈ 0..m−2、w_{v1} += x、w_{v2} −= x+2+e₁+e₂，接受 ⟺ 全部
/// ≡ 0 mod (m−1)；v(m) = total·(m−1)^{−h₀}·scalar(m)；插值点步长
/// step = lcm(C 分母)，P(0) 拉格朗日，整体 ×(−1)^{d−g}（d = g−1+ΣC）。
pub fn witten_class(g: u32, a: &[(i64, i64, i64)]) -> Result<TautClass, String> {
    let n = a.len() as u32;
    let dsum: i64 = a.iter().map(|&(_, _, d)| d).sum();
    if dsum != 2 * g as i64 - 2 {
        return Err(format!("ΣD={dsum} ≠ 2g−2={}", 2 * g as i64 - 2));
    }
    if a.iter().any(|&(cn, cd, d)| cn == cd && d == -1) {
        return Ok(TautClass::default());
    }
    let mut step: i64 = 1;
    for &(_, cd, _) in a {
        if cd <= 0 {
            return Err("C 分母须为正".into());
        }
        step = num_lcm(step, cd);
    }
    let cnum_sum: i64 = a.iter().map(|&(cn, _, _)| cn).sum();
    if cnum_sum % step != 0 {
        return Err("ΣC 须为整数".into());
    }
    let c_sum = cnum_sum / step;
    let d_deg = (g as i64 - 1 + c_sum).max(0) as u32;
    let abs_d: i64 = a.iter().map(|&(_, _, d)| d.abs()).sum();
    let m0 = ((abs_d + 1) / 2 + g as i64 + 1) * step;

    let mut terms: Vec<DecStratum> = Vec::new();
    for gs in tautgens(g, n, d_deg) {
        let graph = &gs.graph;
        let nv = graph.vertices().len();
        let t = &gs.poly.terms[0];
        let ne = graph.edges().len();
        let h0 = ne as i64 - nv as i64 + 1;
        // 边表（顶点对 + 半边 ψ 幂）
        let leg_vertex = |l: u32| -> usize {
            graph
                .vertices()
                .iter()
                .position(|v| v.legs.contains(&l))
                .unwrap()
        };
        let mut edges_v: Vec<(usize, usize)> = Vec::new();
        let mut exps: Vec<(u32, u32)> = Vec::new();
        for &(ea, eb) in graph.edges() {
            edges_v.push((leg_vertex(ea), leg_vertex(eb)));
            exps.push((t.psi_power_at(ea), t.psi_power_at(eb)));
        }
        let exp_sum: u32 = exps.iter().map(|&(e1, e2)| e1 + e2).sum();
        let deg = 2 * exp_sum as i64 + 2 * ne as i64;
        let npts = deg as usize + 1;
        let mut pts: Vec<i64> = Vec::with_capacity(npts);
        let mut vals: Vec<Q> = Vec::with_capacity(npts);
        // scalar 在 X=0 处的符号因子值（腿/κ 因子不参与插值——参考
        // interp(X)·scalar(X) 结构：插值只对 total·(m−1)^{−h0}）
        let mut scalar0 = Q::ONE;
        for vi in 0..nv {
            for &idx in &t.kappa[vi] {
                scalar0 = scalar0 * leg_factor(idx, 0, 0) * Q::from(-1);
            }
        }
        for (li, &(_, _, dd)) in a.iter().enumerate().take(n as usize) {
            let l = li as u32 + 1;
            let e = t.psi_power_at(l);
            scalar0 *= leg_factor(e, dd, 0);
        }
        for i in 0..npts {
            let m = m0 + step * (i as i64 + 1);
            // given 权重与 scalar（每 m 重算）
            // 给定权初值 = 1 − g_v（Pixton M[i,0][0] = genus；g=1 主图
            // 锚点反推验证接受条件成立）
            let mut given: Vec<i64> = graph
                .vertices()
                .iter()
                .map(|v| 1 - v.genus as i64)
                .collect();
            for (gv, kv) in given.iter_mut().zip(&t.kappa) {
                // kappa 列表 = 指标多重集，每个出现即幂 1
                *gv -= kv.iter().map(|&idx| idx as i64).sum::<i64>();
            }
            for (li, &(cn, cd, dd)) in a.iter().enumerate().take(n as usize) {
                let l = li as u32 + 1;
                let vi = graph
                    .vertices()
                    .iter()
                    .position(|v| v.legs.contains(&l))
                    .ok_or(format!("标记 {l} 不在图中"))?;
                let e = t.psi_power_at(l) as i64;
                let ai_m = cn * m / cd + dd; // step | m 保证整除
                given[vi] += ai_m - e;
            }
            // 边权枚举 mod (m−1)：起点 += x、终点 −= x+2+e₁+e₂。常数部分
            // 并入给定值，剩下的 ±x 即平衡流；边因子按 x ∈ 0..m−1 制表
            // （同 (e₁, e₂) 的边共用一张表）
            let mm = m - 1;
            for (&(_, v), &(e1, e2)) in edges_v.iter().zip(&exps) {
                given[v] -= 2 + e1 as i64 + e2 as i64;
            }
            let tables: Vec<std::rc::Rc<Vec<Q>>> =
                exps.iter().map(|&(e1, e2)| edge_table(e1, e2, m)).collect();
            let table: Vec<&[Q]> = tables.iter().map(|t| t.as_slice()).collect();
            let total = crate::flow::balanced_flow_sum(&edges_v, &given, mm, &table);
            // total·(m−1)^{−h0}：参考 mvalues.append(total * ZZ(m-1)**(-h0))
            let mut mpow = Q::ONE;
            for _ in 0..h0 {
                mpow /= Q::from(mm);
            }
            for _ in 0..(-h0) {
                mpow *= Q::from(mm);
            }
            pts.push(m);
            vals.push(total * mpow);
        }
        // 拉格朗日 P(0)
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
        let sgn = Q::from(if ((d_deg as i64 - g as i64) % 2 + 2) % 2 == 0 {
            1
        } else {
            -1
        });
        let coeff = p0 * sgn * scalar0;
        if coeff != Q::ZERO {
            // coeff = mpoly(0)·sgn 已是类系数口径（= aut × oracle 生成元
            // 系数）；evaluate 统一除 aut，勿在此再乘
            terms.push(DecStratum {
                graph: graph.clone(),
                poly: KP::new(vec![KpTerm {
                    kappa: t.kappa.clone(),
                    psi: t.psi.clone(),
                    coeff,
                }]),
            });
        }
    }
    Ok(TautClass::new(terms))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bernoulli_convention() {
        assert_eq!(crate::socle::bernoulli(1).to_string(), "-1/2");
        assert_eq!(crate::socle::bernoulli(2).to_string(), "1/6");
    }

    #[test]
    fn pm_anchors() {
        assert_eq!(pm_eval(0, 5, 1).to_string(), "1");
        assert_eq!(pm_eval(1, 5, 3).to_string(), "3/8");
        assert_eq!(pm_eval(1, 7, 0).to_string(), "-65/24");
        assert_eq!(pm_eval(1, 6, -2).to_string(), "-53/6");
        assert_eq!(pm_eval(2, 5, 4).to_string(), "441/128");
        assert_eq!(pm_eval(2, 9, 0).to_string(), "39865/1152");
        assert_eq!(pm_eval(3, 7, 2).to_string(), "5977855/82944");
    }

    #[test]
    fn leg_anchors() {
        assert_eq!(leg_factor(1, 0, 9).to_string(), "-119/24");
        assert_eq!(leg_factor(2, 0, 9).to_string(), "39865/1152");
        assert_eq!(leg_factor(1, 3, 9).to_string(), "61/24");
        assert_eq!(leg_factor(2, -1, 9).to_string(), "39865/1152");
        assert_eq!(leg_factor(3, 2, 11).to_string(), "15181047/5120");
    }

    /// Witten 极限类积分锚点（tools/witten_anchors.py oracle 实测）。
    #[test]
    fn witten_int_anchors() {
        use crate::classes::int_against;
        // A = (2, r−2)：g=1，d = 1
        let w1 = witten_class(1, &[(0, 1, 2), (1, 1, -2)]).unwrap();
        assert_eq!(
            int_against(&w1, &[], &[(1, 1)]).unwrap().to_string(),
            "-1/8"
        );
        assert_eq!(
            int_against(&w1, &[(1, 1)], &[]).unwrap().to_string(),
            "-1/4"
        );
        assert_eq!(
            int_against(&w1, &[], &[(2, 1)]).unwrap().to_string(),
            "-1/8"
        );
        // A = (3, r/2−1, r/2−2)：g=1，d = 1
        let w2 = witten_class(1, &[(0, 1, 3), (1, 2, -1), (1, 2, -2)]).unwrap();
        assert_eq!(
            int_against(&w2, &[], &[(1, 1), (2, 1)])
                .unwrap()
                .to_string(),
            "-1/2"
        );
        assert_eq!(
            int_against(&w2, &[(1, 1)], &[(3, 1)]).unwrap().to_string(),
            "-5/4"
        );
        // g=2 例（oracle 亦为 0）
        let w3 = witten_class(2, &[(0, 1, 2)]).unwrap();
        assert_eq!(int_against(&w3, &[(2, 1)], &[]).unwrap().to_string(), "0");
        assert_eq!(int_against(&w3, &[(1, 2)], &[]).unwrap().to_string(), "0");
        let w4 = witten_class(2, &[(0, 1, 1), (0, 1, 1)]).unwrap();
        assert_eq!(int_against(&w4, &[(2, 1)], &[]).unwrap().to_string(), "0");
        // g=2 A=(2,) 差分（tools/witten_diff.py，oracle 5 强例全中）
        let w5 = witten_class(2, &[(0, 1, 2)]).unwrap();
        assert_eq!(
            int_against(&w5, &[(3, 1)], &[]).unwrap().to_string(),
            "11/1920"
        );
        assert_eq!(
            int_against(&w5, &[(1, 1), (2, 1)], &[])
                .unwrap()
                .to_string(),
            "73/1920"
        );
        assert_eq!(
            int_against(&w5, &[(2, 1)], &[(1, 1)]).unwrap().to_string(),
            "17/1920"
        );
        assert_eq!(
            int_against(&w5, &[(1, 1)], &[(1, 2)]).unwrap().to_string(),
            "3/640"
        );
        assert_eq!(
            int_against(&w5, &[], &[(1, 3)]).unwrap().to_string(),
            "1/1920"
        );
    }

    #[test]
    fn edge_anchors() {
        assert_eq!(edge_factor(0, 7, 0, 0).to_string(), "5/24");
        assert_eq!(edge_factor(1, 7, 0, 0).to_string(), "-31/24");
        assert_eq!(edge_factor(2, 9, 1, 0).to_string(), "10439/1152");
        assert_eq!(edge_factor(0, 11, 1, 1).to_string(), "2176461/2560");
        assert_eq!(edge_factor(3, 11, 2, 0).to_string(), "17161833/5120");
    }
}
