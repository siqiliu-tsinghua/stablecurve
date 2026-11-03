//! κ、λ、ch 与 ψ 类的积分。klp 先用带符号集合分划把 κ 化为后代，
//! 再走 Mumford GRR 的直接后代递归。int_class 保留 Faber 旧路线供对照。
//!
//! int_class 的旧路线（KaLa5）：
//! - λ_i 按 Mumford 生成函数展开为奇 ch：
//!   `λ_n = [t^n] exp( Σ_i ch_{2i−1}·(2i−2)!·t^{2i−1} )`（偶 ch 在 M_g 消没）；
//! - g=1 带一个 ch₁：恒等式 1/24·∫_{M̄_{0,n+2}}（类不变、双标记结构）；
//! - g≥2 带 ch_m：Mumford 公式
//!   `ch_m·α = B_{m+1}/(m+1)! · [ κ_m·α' + ½Σ_j(−1)^j(Δ₀ 项) + ½Σ_h(可约项) ]`，
//!   其中 ψκ 部分先 push 到 (g,0)；Δ₀ 项落到 (g−1,2) 带 ψ₁^jψ₂^{m−1−j}；
//!   可约项按 (h, g−h) 分裂，各类在 D/F 两侧二项式分配、结点 ψ 带 (−1)^j，
//!   D 侧恰好顶维（3h−2）的分配保留，递归为两个 (·,1) 积分之积。
//!
//! 与原版的实现差异：类以**语义指标**传递（κ_j、ch_c），省去原版的线性
//! 下标平移（i+2/i+1/−2n 等在两侧同义）。

use crate::kppoly::{KappaPsiPolynomial as KP, KpTerm};
use crate::push::push;
use crate::socle::bernoulli;
use crate::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero, is_zero};

mod psi;

pub(crate) fn psi_integral(g: u32, powers: &[u32]) -> Q {
    psi::pure(g, powers)
}

pub(crate) fn psi_sum(g: u32, roots: &[(Vec<u32>, Q)]) -> Q {
    psi::pure_sum(g, roots)
}

fn fact(n: u64) -> Q {
    let mut r = Q::ONE;
    for i in 2..=n {
        r *= Q::from(i as i64);
    }
    r
}

fn binom(e: u32, k: u32) -> Q {
    if k > e {
        return Q::ZERO;
    }
    let mut r = Q::ONE;
    for i in 0..k.min(e - k) {
        r = r * Q::from(e - i) / Q::from(i + 1);
    }
    r
}

/// ch-单项式：`Vec<(ch指标, 幂)>` + 系数。
type ChMono = (Vec<(u32, u32)>, Q);

/// λ_n 展开为奇 ch 多项式（Faber `lamb`）：
/// `λ_n = [t^n] exp(Σ_i ch_{2i−1}(2i−2)! t^{2i−1})`——逐生成元乘 Σ_k f^k/k!。
pub fn lambda_expansion(n: u32) -> Vec<ChMono> {
    let mut e: Vec<(u32, Vec<ChMono>)> = vec![(0, vec![(vec![], Q::ONE)])];
    for deg in (1..=n).step_by(2) {
        let fac = fact((deg - 1) as u64); // (2i−2)!
        let snap = e.clone();
        for (d0, poly) in &snap {
            let mut k = 1u32;
            loop {
                let nd = d0 + k * deg;
                if nd > n {
                    break;
                }
                // fac^k / k!
                let mut cf = Q::ONE;
                for _ in 0..k {
                    cf *= &fac;
                }
                cf /= fact(k as u64);
                for (mono, c) in poly {
                    let mut m2 = mono.clone();
                    m2.push((deg, k));
                    m2.sort_unstable();
                    if let Some(slot) = e.iter_mut().find(|(d, _)| *d == nd) {
                        slot.1.push((m2, c * &cf));
                    } else {
                        e.push((nd, vec![(m2, c * &cf)]));
                    }
                }
                k += 1;
            }
        }
    }
    e.sort_by_key(|(d, _)| *d);
    e.into_iter()
        .find(|(d, _)| *d == n)
        .map(|(_, v)| normalize_ch(v))
        .unwrap_or_else(|| vec![(vec![], Q::ZERO)])
}

fn normalize_ch(mut v: Vec<ChMono>) -> Vec<ChMono> {
    for (mono, _) in &mut v {
        *mono = normalized_powers(mono);
    }
    v.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    v.dedup_by(|term, previous| {
        if term.0 == previous.0 {
            previous.1 += &term.1;
            true
        } else {
            false
        }
    });
    v.retain(|(_, c)| !is_zero(c));
    v
}

fn normalized_powers(powers: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let mut out = powers.to_vec();
    out.retain(|&(_, p)| p != 0);
    out.sort_unstable_by_key(|v| v.0);
    out.dedup_by(|v, prev| {
        if v.0 == prev.0 {
            prev.1 = prev.1.checked_add(v.1).expect("指数之和超出 u32");
            true
        } else {
            false
        }
    });
    out
}

/// ch-多项式相乘。
fn ch_mul(a: &[ChMono], b: &[ChMono]) -> Vec<ChMono> {
    let mut out = Vec::new();
    for (ma, ca) in a {
        for (mb, cb) in b {
            let mut m = ma.clone();
            m.extend_from_slice(mb);
            m.sort_unstable();
            m.sort_by_key(|x| x.0);
            // 合并同指标
            let mut merged: Vec<(u32, u32)> = Vec::new();
            for (c, p) in m {
                if merged.last().map(|e| e.0 == c).unwrap_or(false) {
                    merged.last_mut().unwrap().1 += p;
                } else {
                    merged.push((c, p));
                }
            }
            out.push((merged, ca * cb));
        }
    }
    normalize_ch(out)
}

// int_class 的进程内记忆化。
thread_local! {
    static IC_MEMO: std::cell::RefCell<
        std::collections::HashMap<
            (u32, u32, Vec<(u32, u32)>, Vec<(u32, u32)>, Vec<(u32, u32)>),
            Q,
        >,
    > = std::cell::RefCell::new(std::collections::HashMap::new());
}

/// ∫_{M̄_{g,n}} κ^κ · ψ^ψ · ch^ch（全部语义指标；单项式）。
pub fn int_class(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<Q, String> {
    use crate::socle::{cache_get, cache_put};
    // 键规范化：各分量排序，避免同子问题因构造路径不同而 miss
    let kappa = normalized_powers(kappa);
    let psi = normalized_powers(psi);
    let ch = normalized_powers(ch);
    // 权重守卫（Faber klmgn 的 sx 检查）：非顶权 → 积分为 0。
    // 递归内部调用自顶权输入而来均为顶权，不受影响。
    let w: u32 = kappa.iter().map(|&(i, p)| i * p).sum::<u32>()
        + psi.iter().map(|&(_, p)| p).sum::<u32>()
        + ch.iter().map(|&(c, p)| c * p).sum::<u32>();
    if w as i64 != 3 * g as i64 - 3 + n as i64 {
        return Ok(Q::ZERO);
    }
    let key = (g, n, kappa.clone(), psi.clone(), ch.clone());
    if let Some(v) = IC_MEMO.with(|m| m.borrow().get(&key).cloned()) {
        return Ok(v);
    }
    // 键编码必须含全部三个分量的长度（否则 ψ/ch 分界歧义 → 键碰撞 → 脏值）
    let flat: Vec<u32> = [
        g,
        n,
        (kappa.len() * 2) as u32,
        (psi.len() * 2) as u32,
        (ch.len() * 2) as u32,
    ]
    .into_iter()
    .chain(kappa.iter().flat_map(|&(a, b)| [a, b]))
    .chain(psi.iter().flat_map(|&(a, b)| [a, b]))
    .chain(ch.iter().flat_map(|&(a, b)| [a, b]))
    .collect();
    if let Some(v) = cache_get("ic", &flat) {
        IC_MEMO.with(|m| m.borrow_mut().insert(key.clone(), v.clone()));
        return Ok(v);
    }
    let r = int_class_inner(g, n, &kappa, &psi, &ch);
    if std::env::var_os("SC_FABER_TRACE").is_some() {
        eprintln!(
            "TRACE int_class(g={g},n={n},k={kappa:?},p={psi:?},ch={ch:?}) = {:?}",
            r.as_ref()
                .map(|q| q.to_string())
                .unwrap_or_else(|e| e.clone())
        );
    }
    if let Ok(v) = &r {
        IC_MEMO.with(|m| m.borrow_mut().insert(key, v.clone()));
        cache_put("ic", &flat, v);
    }
    r
}

fn int_class_inner(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<Q, String> {
    if ch.is_empty() {
        // 基本类上的 κ-ψ 单项式 → socle
        let cls = fund_monomial(g, n, kappa, psi)?;
        return cls.evaluate();
    }
    for &(c, _) in ch {
        if c % 2 == 0 {
            return Err(format!("ch_{c} 无定义（Hodge 丛只用奇数次 ch）"));
        }
        // ch_c(E) = 0（c ≥ 2g）。Mumford 递归把 ch 分到低亏格分支时会出现这种
        // 项（如 λ_g² 型单项式），原先报错；admcycles 的 Hodge 数据给出 0。
        if c + 1 > 2 * g {
            return Ok(Q::ZERO);
        }
    }
    if g == 1 {
        // Faber g=1 分支：恰一个 ch₁ → 1/24·(0,n+2)；否则 0（含 ch₁² 的
        // 可疑情形——见 FABER_REDESIGN 对拍记录）
        if ch.len() == 1 && ch[0] == (1, 1) {
            let inner = int_class(0, n + 2, kappa, psi, &[])?;
            return Ok(Q::from_signeds(1, 24) * inner);
        }
        return Ok(Q::ZERO);
    }
    // g ≥ 2
    if n > 0 {
        // ψκ 部分推到 (g,0)，逐 κ-单项式带上全部 ch 递归
        let poly = fund_poly(g, n, kappa, psi)?;
        let pushed = push(g, n, &poly)?;
        let mut acc = Q::ZERO;
        for t in &pushed.terms {
            let part = t.kappa.first().cloned().unwrap_or_default();
            // 分划 → (指标, 幂) 表
            let mut k2: Vec<(u32, u32)> = Vec::new();
            for &i in &part {
                if let Some(e) = k2.iter_mut().find(|e| e.0 == i) {
                    e.1 += 1;
                } else {
                    k2.push((i, 1));
                }
            }
            acc += &t.coeff * &int_class(g, 0, &k2, &[], ch)?;
        }
        return Ok(acc);
    }
    // Mumford 核心（n = 0）
    let m: u32 = ch.iter().map(|&(c, _)| c).max().unwrap();
    let facm = bernoulli((m + 1) as u64) / fact((m + 1) as u64);
    // 消耗一个 ch_m
    let mut rest_ch: Vec<(u32, u32)> = Vec::new();
    let mut consumed = false;
    for &(c, p) in ch {
        if c == m && !consumed {
            consumed = true;
            if p > 1 {
                rest_ch.push((c, p - 1));
            }
        } else {
            rest_ch.push((c, p));
        }
    }
    // κ_m 项
    let mut k1 = kappa.to_vec();
    if let Some(e) = k1.iter_mut().find(|e| e.0 == m) {
        e.1 += 1;
    } else {
        k1.push((m, 1));
    }
    k1.sort_unstable();
    let mut ans = int_class(g, 0, &k1, &[], &rest_ch)?;
    // Δ₀ 项：½ Σ_j (−1)^j · (g−1,2)
    for j in 0..m {
        let mut psi2: Vec<(u32, u32)> = Vec::new();
        if j > 0 {
            psi2.push((1, j));
        }
        if m - 1 - j > 0 {
            psi2.push((2, m - 1 - j));
        }
        let v = int_class(g - 1, 2, kappa, &psi2, &rest_ch)?;
        let sign = if j % 2 == 0 { 1 } else { -1 };
        ans += Q::from_signeds(sign, 2) * v;
    }
    // 可约项：½ Σ_h Σ_分配
    ans += Q::from_signeds(1, 2) * reducible_sum(g, m, kappa, &rest_ch)?;
    Ok(facm * ans)
}

/// Mumford 可约项内层：Σ_{h=1..g−1} Σ_{D/F 分配} (−1)^j·∏C(e,k)·
/// int(h,1,D)·int(g−h,1,F)，D 侧权恰为 3h−2。
fn reducible_sum(
    g: u32,
    m: u32,
    kappa: &[(u32, u32)],
    rest_ch: &[(u32, u32)],
) -> Result<Q, String> {
    let mut total = Q::ZERO;
    for h in 1..g {
        let target = 3 * h - 2;
        for j in 0..m {
            // 待分配的类表：(kind, index, exponent)：kind 0=κ, 1=ch
            let mut items: Vec<(u8, u32, u32)> = Vec::new();
            for &(i, p) in kappa {
                items.push((0, i, p));
            }
            for &(c, p) in rest_ch {
                items.push((1, c, p));
            }
            let mut d_side: Vec<(u8, u32, u32)> = Vec::new();
            let mut f_side: Vec<(u8, u32, u32)> = Vec::new();
            let mut coeff = Q::from(if j % 2 == 0 { 1 } else { -1 });
            total += dist_rec(
                &items,
                0,
                &mut d_side,
                &mut f_side,
                &mut coeff,
                h,
                g,
                m,
                j,
                target,
            )?;
        }
    }
    Ok(total)
}

fn dist_rec(
    items: &[(u8, u32, u32)],
    idx: usize,
    d: &mut Vec<(u8, u32, u32)>,
    f: &mut Vec<(u8, u32, u32)>,
    coeff: &mut Q,
    h: u32,
    g: u32,
    m: u32,
    j: u32,
    target: u32,
) -> Result<Q, String> {
    if idx == items.len() {
        // 结点 ψ 计入 D 侧权重
        let dw: u32 = j + d.iter().map(|&(_, i, p)| i * p).sum::<u32>();
        if dw != target {
            return Ok(Q::ZERO);
        }
        let d_kappa: Vec<(u32, u32)> = d
            .iter()
            .filter(|&&(k, _, _)| k == 0)
            .map(|&(_, i, p)| (i, p))
            .collect();
        let d_ch: Vec<(u32, u32)> = d
            .iter()
            .filter(|&&(k, _, _)| k == 1)
            .map(|&(_, i, p)| (i, p))
            .collect();
        let f_kappa: Vec<(u32, u32)> = f
            .iter()
            .filter(|&&(k, _, _)| k == 0)
            .map(|&(_, i, p)| (i, p))
            .collect();
        let f_ch: Vec<(u32, u32)> = f
            .iter()
            .filter(|&&(k, _, _)| k == 1)
            .map(|&(_, i, p)| (i, p))
            .collect();
        let mut d_psi = vec![(1u32, j)];
        let mut f_psi = vec![(1u32, m - 1 - j)];
        d_psi.retain(|&(_, p)| p > 0);
        f_psi.retain(|&(_, p)| p > 0);
        let vd = int_class(h, 1, &d_kappa, &d_psi, &d_ch)?;
        let vf = int_class(g - h, 1, &f_kappa, &f_psi, &f_ch)?;
        return Ok(coeff.clone() * vd * vf);
    }
    let (kind, index, e) = items[idx];
    let mut sum = Q::ZERO;
    for k in 0..=e {
        let c0 = coeff.clone();
        *coeff *= binom(e, k);
        if k > 0 {
            d.push((kind, index, k));
        }
        if e - k > 0 {
            f.push((kind, index, e - k));
        }
        sum += dist_rec(items, idx + 1, d, f, coeff, h, g, m, j, target)?;
        if k > 0 {
            d.pop();
        }
        if e - k > 0 {
            f.pop();
        }
        *coeff = c0;
    }
    Ok(sum)
}

/// 基本类上的 κ-ψ 单项式 → KappaPsiPolynomial。
pub fn fund_poly(_g: u32, _n: u32, kappa: &[(u32, u32)], psi: &[(u32, u32)]) -> Result<KP, String> {
    let mut part: Vec<u32> = Vec::new();
    for &(i, p) in kappa {
        if i == 0 {
            return Err("κ_0 无定义".into());
        }
        for _ in 0..p {
            part.push(i);
        }
    }
    part.sort_unstable();
    let mut kterm = KpTerm {
        kappa: vec![part],
        psi: psi.to_vec(),
        coeff: Q::ONE,
    };
    kterm.psi.retain(|&(_, p)| p > 0);
    kterm.psi.sort_unstable();
    Ok(KP::new(vec![kterm]))
}

/// 基本类上的 κ-ψ 单项式 → TautClass。
fn fund_monomial(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<TautClass, String> {
    let poly = fund_poly(g, n, kappa, psi)?;
    let legs: Vec<u32> = (1..=n).collect();
    let graph = crate::graph::StableGraph::new(&[g], &[legs], &[])?;
    Ok(TautClass::new(vec![DecStratum { graph, poly }]))
}

/// 完整 klp：∫_{M̄_{g,n}} κ^κ·ψ^ψ·λ^λ·ch^ch。
///
/// 纯 λ-ψ 的较大查询自动并行（默认最多 24 个可用线程）。
/// `SC_FABER_THREADS=1` 使用串行递归；正整数可指定并行线程数上限。
/// 并行任务共享已完成的子积分，同一调用线程的后续查询继续复用结果。
pub fn klp(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<Q, String> {
    // 输入校验：ψ 腿 ⊆ 1..=n（κ 指标 >0 由 fund_poly 保证）；同腿 ψ 幂合并
    crate::kppoly::check_marking_legs(psi, n)?;
    let psi = &crate::kppoly::merge_psi(psi)[..];
    let kappa = &normalized_powers(kappa)[..];
    let ch = &normalized_powers(ch)[..];
    let lambda = &normalized_powers(lambda)[..];
    if 2 * g as i64 - 2 + n as i64 <= 0 {
        return Err(format!("M̄_({g},{n}) 不稳定"));
    }
    // Hodge 丛秩为 g：λ_l = 0（l > g）；g = 0 时 ch_c = 0（c ≥ 1）
    if lambda.iter().any(|&(l, p)| p > 0 && l > g) {
        return Ok(Q::ZERO);
    }
    if g == 0 && ch.iter().any(|&(c, p)| p > 0 && c > 0) {
        return Ok(Q::ZERO);
    }
    if kappa.is_empty() && ch.is_empty() {
        return Ok(psi::integrate(g, n, psi, lambda));
    }
    if kappa.iter().any(|&(i, _)| i == 0) {
        return Err("κ_0 无定义".into());
    }
    let weight = kappa
        .iter()
        .chain(lambda)
        .chain(ch)
        .map(|&(i, p)| u64::from(i) * u64::from(p))
        .sum::<u64>()
        + psi.iter().map(|&(_, p)| u64::from(p)).sum::<u64>();
    if weight as i64 != 3 * i64::from(g) - 3 + i64::from(n) {
        return Ok(Q::ZERO);
    }
    for &(c, _) in ch {
        if c % 2 == 0 {
            return Err(format!("ch_{c} 无定义（Hodge 丛只用奇数次 ch）"));
        }
        if u64::from(c) >= 2 * u64::from(g) {
            return Ok(Q::ZERO);
        }
    }
    Ok(psi::integrate_mixed(g, n, kappa, psi, lambda, ch))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_powers_and_duplicate_factors() {
        let expected = Q::from_signeds(1, 24);
        assert_eq!(klp(1, 1, &[], &[(1, 1)], &[], &[(1, 0)]).unwrap(), expected);
        assert_eq!(
            int_class(1, 1, &[], &[(1, 1)], &[(1, 0)]).unwrap(),
            expected
        );
        assert_eq!(
            klp(2, 0, &[], &[], &[(1, 1), (1, 2)], &[]).unwrap(),
            klp(2, 0, &[], &[], &[(1, 3)], &[]).unwrap()
        );
        assert_eq!(binom(62, 28), Q::from(349615716557887465i64));
    }

    #[test]
    fn mixed_direct_recursion_matches_legacy_grr() {
        for g in 1..=3 {
            for n in 1..=2 {
                let dimension = 3 * g - 3 + n;
                for k in 0..=3.min(dimension) {
                    for l in 0..=g.min(dimension - k) {
                        for c in 0..=2.min(dimension - k - l) {
                            let kappa = if k == 0 { vec![] } else { vec![(1, k)] };
                            let lambda = if l == 0 { vec![] } else { vec![(l, 1)] };
                            let ch = if c == 0 { vec![] } else { vec![(1, c)] };
                            let psi = vec![(1, dimension - k - l - c)];
                            let expansion = ch_mul(&lambda_expansion(l), &[(ch.clone(), Q::ONE)]);
                            let mut expected = Q::ZERO;
                            for (mono, coefficient) in expansion {
                                expected +=
                                    coefficient * int_class(g, n, &kappa, &psi, &mono).unwrap();
                            }
                            assert_eq!(
                                klp(g, n, &kappa, &psi, &lambda, &ch).unwrap(),
                                expected,
                                "g={g} n={n} k={k} l={l} ch1={c}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// Instructions（Faber 2004）的全部数值例——M4 验收。
    #[test]
    fn instructions_examples() {
        // > klp(2,2,ps1^2*la1^3) = 1/1440
        assert_eq!(
            klp(2, 2, &[], &[(1, 2)], &[(1, 3)], &[]).unwrap(),
            Q::from_signeds(1, 1440)
        );
        // > klp(2,2,ps2*ka2*la1^2) = 7/720
        assert_eq!(
            klp(2, 2, &[(2, 1)], &[(2, 1)], &[(1, 2)], &[]).unwrap(),
            Q::from_signeds(7, 720)
        );
        // > klp(3,0,ka1*la3*la2) = 1/120960（= 1/42·2/8/720）
        assert_eq!(
            klp(3, 0, &[(1, 1)], &[], &[(3, 1), (2, 1)], &[]).unwrap(),
            Q::from_signeds(1, 120960)
        );
        // > klp(3,1,ps1^2*la3*la2) = 1/120960
        assert_eq!(
            klp(3, 1, &[], &[(1, 2)], &[(3, 1), (2, 1)], &[]).unwrap(),
            Q::from_signeds(1, 120960)
        );
        // > klp(4,2,ps1^2*ps2^2*la4*la3) = 1/276480
        assert_eq!(
            klp(4, 2, &[], &[(1, 2), (2, 2)], &[(4, 1), (3, 1)], &[]).unwrap(),
            Q::from_signeds(1, 276480)
        );
        // > klp(4,1,ps1^3*la4*la3) = 1/3225600
        assert_eq!(
            klp(4, 1, &[], &[(1, 3)], &[(4, 1), (3, 1)], &[]).unwrap(),
            Q::from_signeds(1, 3225600)
        );
        // λ_g^k 系列之轻量部分
        assert_eq!(
            klp(3, 0, &[], &[], &[(1, 6)], &[]).unwrap(),
            Q::from_signeds(1, 90720)
        );
        assert_eq!(
            klp(4, 0, &[], &[], &[(1, 9)], &[]).unwrap(),
            Q::from_signeds(1, 113400)
        );
        // ch 直接输入：klmgn(2,2,[1,1,8$3]) = 1/1440（ψ₁²ch₁³）；[1,1,9] = −1/17280（ψ₁²ch₃）
        assert_eq!(
            klp(2, 2, &[], &[(1, 2)], &[], &[(1, 3)]).unwrap(),
            Q::from_signeds(1, 1440)
        );
        assert_eq!(
            klp(2, 2, &[], &[(1, 2)], &[], &[(3, 1)]).unwrap(),
            Q::from_signeds(-1, 17280)
        );
        // 度数不符 → 0（klmgn(2,2,[1,1,1]) 的演示例）
        assert_eq!(int_class(2, 2, &[], &[(1, 3)], &[]).unwrap(), Q::ZERO);
    }

    #[test]
    fn lambda_expansions() {
        // λ₁ = ch₁
        assert_eq!(lambda_expansion(1), vec![(vec![(1, 1)], Q::ONE)]);
        // λ₂ = ch₁²/2（偶 ch 消没）
        assert_eq!(
            lambda_expansion(2),
            vec![(vec![(1, 2)], Q::from_signeds(1, 2))]
        );
        // λ₃ = 2ch₃ + ch₁³/6
        let e3 = lambda_expansion(3);
        assert_eq!(e3.len(), 2);
        assert!(e3.contains(&(vec![(3, 1)], Q::from(2))));
        assert!(e3.contains(&(vec![(1, 3)], Q::from_signeds(1, 6))));
    }

    /// 较重的 λ 例子（g=5..8；Faber 本人在 2004 服务器上 g=7 需 35s）。
    /// socle 分裂项与 push 的 κ 展开按多重集枚举后，g=8 也只需不到一秒
    /// （此前分别为指数级，g=7/8 约 2.5s/29s，故曾标 ignore）。
    #[test]
    fn instructions_heavy() {
        assert_eq!(
            klp(5, 0, &[], &[], &[(1, 12)], &[]).unwrap(),
            Q::from_signeds(31, 680400)
        );
        assert_eq!(
            klp(6, 0, &[], &[], &[(1, 15)], &[]).unwrap(),
            Q::from_signeds(431, 481140)
        );
        assert_eq!(
            klp(7, 0, &[], &[], &[(1, 18)], &[]).unwrap(),
            Q::from_signeds(32017001, 638512875)
        );
        // g=8：与 DLYZ Hodge 势 H_8 的提取值、admcycles 的 Schmitt–Sun 数据一致
        assert_eq!(
            klp(8, 0, &[], &[], &[(1, 21)], &[]).unwrap(),
            Q::from_signeds(38176003213i64, 5746615875)
        );
    }

    /// 回归：递归把 ch_c 分到亏格 h < (c+1)/2 的分支时该项为 0（ch_c(E) = 0，
    /// c ≥ 2h），原先报"ch_7 在 g=3 无定义"。λ_g² = 0（Mumford）型单项式即
    /// 触发此路径；admcycles 数据：∫_{M̄_7} λ_4λ_7² = ∫_{M̄_8} λ_5λ_8² = 0。
    #[test]
    fn high_ch_on_low_genus_branch_vanishes() {
        assert_eq!(
            klp(7, 0, &[], &[], &[(4, 1), (7, 2)], &[]).unwrap(),
            Q::ZERO
        );
        assert_eq!(
            klp(8, 0, &[], &[], &[(5, 1), (8, 2)], &[]).unwrap(),
            Q::ZERO
        );
        // 直接输入：ch_9 在 M̄_4 上为 0；偶数次 ch 仍报错
        assert_eq!(klp(4, 0, &[], &[], &[], &[(9, 1)]).unwrap(), Q::ZERO);
        // （偶数次须用顶权输入才走到检查：非顶权在入口即记 0）
        assert!(klp(2, 1, &[], &[(1, 2)], &[], &[(2, 1)]).is_err());
    }

    /// 回归：非顶权含 λ 输入必须为 0（曾返回垃圾值，admcycles 差分发现）。
    #[test]
    fn nontop_weight_is_zero() {
        // (3,1) 顶权 7；λ₂λ₃ψ₁ 权 6、λ₁λ₂ψ₁κ₂ 权 6 → 均 0
        assert_eq!(
            klp(3, 1, &[], &[(1, 1)], &[(2, 1), (3, 1)], &[]).unwrap(),
            Q::ZERO
        );
        assert_eq!(
            klp(3, 1, &[(2, 1)], &[(1, 1)], &[(1, 1), (2, 1)], &[]).unwrap(),
            Q::ZERO
        );
        // 顶权对照（不被守卫误伤）
        assert_eq!(
            klp(3, 1, &[], &[(1, 2)], &[(3, 1), (2, 1)], &[]).unwrap(),
            Q::from_signeds(1, 120960)
        );
    }

    /// 输入规范化回归：同腿 ψ 合并、越界腿报错、秩外 λ 与 g=0 的 Hodge 类为 0。
    #[test]
    fn input_normalization() {
        // ψ₁²·ψ₁² = ψ₁⁴（曾因只取第一项而返回 0）
        assert_eq!(
            klp(2, 1, &[], &[(1, 2), (1, 2)], &[], &[]).unwrap(),
            Q::from_signeds(1, 1152)
        );
        assert!(klp(2, 1, &[], &[(2, 4)], &[], &[]).is_err());
        assert!(klp(1, 0, &[], &[], &[], &[]).is_err()); // 不稳定
        // λ_3 在 M̄_{2,n} 上为 0（Hodge 丛秩 2）
        assert_eq!(klp(2, 2, &[], &[(1, 2)], &[(3, 1)], &[]).unwrap(), Q::ZERO);
        // g = 0：Hodge 丛为零
        assert_eq!(klp(0, 4, &[], &[], &[(1, 1)], &[]).unwrap(), Q::ZERO);
        assert_eq!(klp(0, 4, &[], &[], &[], &[(1, 1)]).unwrap(), Q::ZERO);
    }

    #[test]
    fn g1_lambda() {
        // ∫λ₁ M̄_{1,1} = 1/24；∫λ₁ψ₁ M̄_{1,2} = 1/24；∫λ₁ψ₁ψ₂ M̄_{1,3} = 1/12
        assert_eq!(
            klp(1, 1, &[], &[], &[(1, 1)], &[]).unwrap(),
            Q::from_signeds(1, 24)
        );
        assert_eq!(
            klp(1, 2, &[], &[(1, 1)], &[(1, 1)], &[]).unwrap(),
            Q::from_signeds(1, 24)
        );
        assert_eq!(
            klp(1, 3, &[], &[(1, 1), (2, 1)], &[(1, 1)], &[]).unwrap(),
            Q::from_signeds(1, 12)
        );
    }
}
