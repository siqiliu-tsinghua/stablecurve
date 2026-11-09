//! κ、λ、ch 与 ψ 类的积分。klp 先尝试闭式与 Mumford 关系化简，再用带符号
//! 集合分划把 κ 化为后代，走 GRR 的直接后代递归。int_class 保留旧路线供对照。
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

use crate::cache_source::legacy_source::LegacySource;
use crate::cache_source::{IntegralKey, Registry};
use crate::kppoly::{KappaPsiPolynomial as KP, KpTerm};
use crate::push::push;
use crate::socle::bernoulli;
use crate::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero, is_zero};
use std::sync::{Arc, Weak};

mod closed;
mod lambda;
mod psi;
mod relations;
pub(crate) mod root_source;

pub use lambda::{lambda_integrals, lambda_monomial_count};

pub(crate) fn clear_cache_context(id: u64) {
    psi::clear_cache_context(id);
    IC_SOURCES.with(|sources| {
        if let Ok(mut sources) = sources.try_borrow_mut() {
            sources.remove(&id);
        }
    });
}

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

type IcKey = (u32, u32, Vec<(u32, u32)>, Vec<(u32, u32)>, Vec<(u32, u32)>);

// int_class 的原内存模式与按缓存归属隔离的可读取数学表。
thread_local! {
    static IC_MEMO: std::cell::RefCell<std::collections::HashMap<IcKey, Q>> =
        std::cell::RefCell::new(std::collections::HashMap::new());
    static IC_SOURCES: std::cell::RefCell<std::collections::BTreeMap<u64, Weak<LegacySource<IcKey>>>> =
        const { std::cell::RefCell::new(std::collections::BTreeMap::new()) };
}

fn ic_source(registry: &Registry) -> Arc<LegacySource<IcKey>> {
    IC_SOURCES.with(|cell| {
        let mut sources = cell.borrow_mut();
        if let Some(source) = sources.get(&registry.id()).and_then(Weak::upgrade) {
            return source;
        }
        sources.retain(|_, source| source.strong_count() != 0);
        let source = Arc::new(LegacySource::new("ic", ic_nums));
        registry.register(source.clone());
        sources.insert(registry.id(), Arc::downgrade(&source));
        source
    })
}

// 三个分量的长度均是扁平 u32 元素数，保持旧 st 表 ic 键的边界语义。
// 正常计算只在确需 L2 查找时编码；后台 pull 自行调用同一个函数。
fn ic_nums(key: &IcKey) -> Vec<u32> {
    let (g, n, kappa, psi, ch) = key;
    [
        *g,
        *n,
        (kappa.len() * 2) as u32,
        (psi.len() * 2) as u32,
        (ch.len() * 2) as u32,
    ]
    .into_iter()
    .chain(
        kappa
            .iter()
            .chain(psi)
            .chain(ch)
            .flat_map(|&(index, power)| [index, power]),
    )
    .collect()
}

/// ∫_{M̄_{g,n}} κ^κ · ψ^ψ · ch^ch（全部语义指标；单项式）。
pub fn int_class(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<Q, String> {
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
        return crate::cache_source::with_current_registry(|registry| {
            let Some(registry) = registry else {
                return Ok(Q::ZERO);
            };
            // The historical non-top zero return precedes input rejection.
            // Keep that return order, but publish only actual stable-space
            // integrals: malformed auxiliary zeros are not mathematical roots.
            if 2 * i64::from(g) - 2 + i64::from(n) <= 0
                || kappa.iter().any(|&(index, _)| index == 0)
                || psi.iter().any(|&(leg, _)| leg == 0 || leg > n)
                || ch.iter().any(|&(index, _)| index == 0 || index % 2 == 0)
            {
                return Ok(Q::ZERO);
            }
            root_source::with_registry_policy(
                registry.clone(),
                g,
                n,
                &kappa,
                &psi,
                &[],
                &ch,
                false,
                || Ok(Q::ZERO),
            )
        });
    }
    let key = (g, n, kappa.clone(), psi.clone(), ch.clone());
    crate::cache_source::with_current_registry(|registry| {
        let source = if let Some(registry) = registry {
            registry.check_error()?;
            let source = ic_source(registry);
            if let Some(value) = source.get(&key) {
                return Ok(value);
            }
            if registry.has_history() {
                if let Some(value) = registry.lookup(IntegralKey::Legacy {
                    kind: "ic".to_owned(),
                    nums: ic_nums(&key),
                })? {
                    source.insert(key, &value);
                    return Ok(value);
                }
            }
            Some(source)
        } else {
            if let Some(value) = IC_MEMO.with(|memo| memo.borrow().get(&key).cloned()) {
                return Ok(value);
            }
            None
        };
        // No mathematical source lock or TLS memo borrow spans recursion.
        let result = int_class_inner(g, n, &kappa, &psi, &ch);
        if std::env::var_os("SC_FABER_TRACE").is_some() {
            eprintln!(
                "TRACE int_class(g={g},n={n},k={kappa:?},p={psi:?},ch={ch:?}) = {:?}",
                result
                    .as_ref()
                    .map(|q| q.to_string())
                    .unwrap_or_else(|error| error.clone())
            );
        }
        if let Ok(value) = &result {
            if let Some(source) = source {
                source.insert(key, value);
            } else {
                IC_MEMO.with(|memo| memo.borrow_mut().insert(key, value.clone()));
            }
        }
        result
    })
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
/// 先尝试 Hodge 闭式与 Mumford 关系；未消去的项合并后交给通用 GRR。
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
    // Errors and the established zero/pure/mixed route are decided before any
    // root cache can return a value for an otherwise rejected input.
    let route = klp_route(g, n, kappa, psi, lambda, ch)?;
    let calculate = || calculate_klp(g, n, kappa, psi, lambda, ch, route);
    // 固定本次调用的结果归属。未启用缓存时不建立根结果目录或后台读取口。
    if let Some(registry) = crate::cache_source::current_registry() {
        let read_history = !(inexpensive_klp(g, n, psi, lambda, ch, route)
            || closed::inexpensive(g, n, kappa, psi, lambda, ch));
        return root_source::with_registry_policy(
            registry,
            g,
            n,
            kappa,
            psi,
            lambda,
            ch,
            read_history,
            calculate,
        );
    }
    crate::cache_source::with_optional_registry(None, calculate)
}

fn unreduced_klp(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
    route: KlpRoute,
) -> Q {
    match route {
        KlpRoute::Zero => Q::ZERO,
        KlpRoute::Pure => psi::integrate(g, n, psi, lambda),
        KlpRoute::Mixed => psi::integrate_mixed(g, n, kappa, psi, lambda, ch),
    }
}

fn calculate_klp(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
    route: KlpRoute,
) -> Result<Q, String> {
    if matches!(route, KlpRoute::Zero) {
        return Ok(Q::ZERO);
    }
    if !lambda.iter().any(|&(index, _)| index > 0) {
        return Ok(unreduced_klp(g, n, kappa, psi, lambda, ch, route));
    }
    // Preserve the original degree/pullback guards before spending work on a
    // lambda normal form. Input errors have already been checked by klp_route.
    let hodge_weight = lambda
        .iter()
        .chain(ch)
        .map(|&(index, power)| u128::from(index) * u128::from(power))
        .sum::<u128>();
    let weight = hodge_weight
        + kappa
            .iter()
            .map(|&(index, power)| u128::from(index) * u128::from(power))
            .sum::<u128>()
        + psi
            .iter()
            .map(|&(_, power)| u128::from(power))
            .sum::<u128>();
    if weight as i128 != 3 * i128::from(g) - 3 + i128::from(n)
        || hodge_weight
            > if g > 1 {
                3 * u128::from(g) - 3
            } else {
                u128::from(g)
            }
        || (g > 0
            && lambda
                .iter()
                .any(|&(index, power)| index == g && power >= 2))
    {
        return Ok(Q::ZERO);
    }
    if let Some(value) = closed::value(g, n, kappa, psi, lambda, ch) {
        return Ok(value);
    }
    // lambda_1=ch_1 and lambda_2=ch_1^2/2 already give one GRR monomial.
    // Expanding them into higher lambda indices usually creates extra ch
    // roots. Only do so when the lambda weight forces every squarefree leaf
    // to contain both top classes, hence all leaves have a closed pairing.
    // The largest squarefree weight missing g or g-1 is g(g-1)/2 + 1.
    if g >= 2
        && lambda.iter().all(|&(index, _)| index <= 2)
        && lambda
            .iter()
            .map(|&(index, power)| u128::from(index) * u128::from(power))
            .sum::<u128>()
            <= u128::from(g) * u128::from(g - 1) / 2 + 1
    {
        return Ok(unreduced_klp(g, n, kappa, psi, lambda, ch, route));
    }
    let Some(reduced) = relations::reduce(g, lambda) else {
        return Ok(unreduced_klp(g, n, kappa, psi, lambda, ch, route));
    };
    let mut result = Q::ZERO;
    let mut closed_terms = Vec::new();
    let mut remaining = Vec::new();
    for (term, coefficient) in reduced {
        let has_top_pair = g >= 2
            && term.iter().any(|&(index, _)| index == g)
            && term.iter().any(|&(index, _)| index == g - 1);
        if has_top_pair || inexpensive_klp(g, n, psi, &term, ch, route) {
            closed_terms.push((term, coefficient));
        } else {
            remaining.push((term, coefficient));
        }
    }
    if !remaining.is_empty() {
        // Expand and combine all remaining leaves before planning. Running a
        // separate descendant graph per leaf loses the shared-GRR advantage.
        match psi::integrate_reduced(g, n, kappa, psi, lambda, &remaining, ch) {
            psi::ReducedIntegral::Reduced(value) => result += value,
            // The backend reused the original ch expansion when the reduced
            // representation would grow the graph. Its value already covers
            // the entire monomial, so do not evaluate the closed leaves.
            psi::ReducedIntegral::Original(value) => return Ok(value),
        }
    }
    for (term, coefficient) in closed_terms {
        // Publish individually evaluated leaves only after accepting the
        // reduction. No source lock spans these ordinary public calls.
        result += coefficient * klp(g, n, kappa, psi, &term, ch)?;
    }
    Ok(result)
}

// Independent pre-reduction route for mathematical oracle tests. It bypasses
// optimized public root values as well as both new simplification modules.
#[cfg(test)]
pub(crate) fn klp_unreduced(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<Q, String> {
    crate::kppoly::check_marking_legs(psi, n)?;
    let psi = crate::kppoly::merge_psi(psi);
    let kappa = normalized_powers(kappa);
    let lambda = normalized_powers(lambda);
    let ch = normalized_powers(ch);
    let route = klp_route(g, n, &kappa, &psi, &lambda, &ch)?;
    crate::cache_source::with_optional_registry(None, || {
        Ok(unreduced_klp(g, n, &kappa, &psi, &lambda, &ch, route))
    })
}

#[derive(Clone, Copy)]
enum KlpRoute {
    Zero,
    Pure,
    Mixed,
}

// Select a cache access policy, without evaluating a second mathematical
// formula. These cases already have a short exact formula or a zero guard in
// the original evaluator, so a cold L1 should not wait for an SQLite request.
fn inexpensive_klp(
    g: u32,
    n: u32,
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
    route: KlpRoute,
) -> bool {
    match route {
        KlpRoute::Zero => true,
        KlpRoute::Mixed => false,
        KlpRoute::Pure => {
            let weight: u64 = lambda
                .iter()
                .map(|&(index, power)| u64::from(index) * u64::from(power))
                .sum();
            let psi_weight: u64 = psi.iter().map(|&(_, power)| u64::from(power)).sum();
            let factors = lambda.iter().filter(|&&(index, _)| index > 0);
            let top_power = lambda
                .iter()
                .filter(|&&(index, _)| index == g && g > 0)
                .map(|&(_, power)| u64::from(power))
                .sum::<u64>();
            weight + psi_weight != (3 * i64::from(g) - 3 + i64::from(n)) as u64
                || top_power >= 2
                || weight
                    > if g > 1 {
                        3 * u64::from(g) - 3
                    } else {
                        u64::from(g)
                    }
                || (ch.is_empty()
                    && ((factors.clone().count() == 0 && (g <= 1 || n == 1))
                        || (g > 0 && factors.copied().eq([(g, 1)]))))
        }
    }
}

fn klp_route(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<KlpRoute, String> {
    if 2 * i64::from(g) - 2 + i64::from(n) <= 0 {
        return Err(format!("M̄_({g},{n}) 不稳定"));
    }
    // Hodge 丛秩为 g：λ_l = 0（l > g）；g = 0 时 ch_c = 0（c ≥ 1）
    if lambda.iter().any(|&(l, p)| p > 0 && l > g) {
        return Ok(KlpRoute::Zero);
    }
    if g == 0 && ch.iter().any(|&(c, p)| p > 0 && c > 0) {
        return Ok(KlpRoute::Zero);
    }
    if kappa.is_empty() && ch.is_empty() {
        return Ok(KlpRoute::Pure);
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
        return Ok(KlpRoute::Zero);
    }
    for &(c, _) in ch {
        if c % 2 == 0 {
            return Err(format!("ch_{c} 无定义（Hodge 丛只用奇数次 ch）"));
        }
        if u64::from(c) >= 2 * u64::from(g) {
            return Ok(KlpRoute::Zero);
        }
    }
    Ok(KlpRoute::Mixed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache_source::{ResultSource, with_optional_registry, with_registry};

    #[test]
    fn reduced_pure_hodge_monomials_match_the_original_backend() {
        fn visit(g: u32, remaining: u32, cap: u32, parts: &mut Vec<(u32, u32)>) {
            if remaining == 0 {
                assert_eq!(
                    klp(g, 0, &[], &[], parts, &[]),
                    klp_unreduced(g, 0, &[], &[], parts, &[]),
                    "g={g}, lambda={parts:?}"
                );
                return;
            }
            for index in (1..=remaining.min(cap)).rev() {
                if let Some(last) = parts.last_mut().filter(|last| last.0 == index) {
                    last.1 += 1;
                    visit(g, remaining - index, index, parts);
                    parts.last_mut().unwrap().1 -= 1;
                } else {
                    parts.push((index, 1));
                    visit(g, remaining - index, index, parts);
                    parts.pop();
                }
            }
        }
        for g in 2..=6 {
            visit(g, 3 * g - 3, g, &mut Vec::new());
        }
    }

    #[test]
    fn reduced_mixed_monomials_match_the_original_backend() {
        for g in 2..=4 {
            for n in 1..=2 {
                for lambda in [
                    vec![(1, 2)],
                    vec![(1, 3)],
                    vec![(2, 2)],
                    vec![(1, 2), (2, 2)],
                    vec![(g - 1, 3)],
                ] {
                    for kappa in [vec![], vec![(1, 1)], vec![(1, 2)]] {
                        for ch in [vec![], vec![(1, 1)], vec![(3, 1)]] {
                            let weight: u32 = lambda
                                .iter()
                                .chain(&kappa)
                                .chain(&ch)
                                .map(|&(index, power)| index * power)
                                .sum();
                            let Some(power) = (3 * g - 3 + n).checked_sub(weight) else {
                                continue;
                            };
                            let psi = [(1, power)];
                            assert_eq!(
                                klp(g, n, &kappa, &psi, &lambda, &ch),
                                klp_unreduced(g, n, &kappa, &psi, &lambda, &ch),
                                "g={g}, n={n}, k={kappa:?}, p={psi:?}, l={lambda:?}, ch={ch:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    #[ignore = "manual small cold comparison of algebraic reduction and the original backend"]
    fn benchmark_hodge_reduction() {
        for (g, lambda) in [
            (5, vec![(1, 12)]),
            (8, vec![(1, 21)]),
            (10, vec![(1, 27)]),
            (9, vec![(2, 12)]),
            (8, vec![(3, 7)]),
            (9, vec![(4, 6)]),
        ] {
            let normalization = std::time::Instant::now();
            let normal_form = relations::reduce(g, &lambda);
            eprintln!(
                "hodge_normal_form g={g} lambda={lambda:?} elapsed_s={:.6} leaves={:?}",
                normalization.elapsed().as_secs_f64(),
                normal_form.as_ref().map(Vec::len)
            );
            let run = |reduced: bool| {
                let lambda = lambda.clone();
                std::thread::spawn(move || {
                    with_optional_registry(None, || {
                        let started = std::time::Instant::now();
                        let value = if reduced {
                            klp(g, 0, &[], &[], &lambda, &[])
                        } else {
                            klp_unreduced(g, 0, &[], &[], &lambda, &[])
                        }
                        .unwrap();
                        (value, started.elapsed())
                    })
                })
                .join()
                .unwrap()
            };
            let (expected, original) = run(false);
            let (actual, reduced) = run(true);
            assert_eq!(actual, expected);
            eprintln!(
                "hodge_reduction g={g} lambda={lambda:?} original_s={:.6} reduced_s={:.6}",
                original.as_secs_f64(),
                reduced.as_secs_f64()
            );
        }
    }

    #[test]
    fn int_class_non_top_legal_roots_persist_zero_without_duplicate_completions() {
        use crate::cache_source::storage::Storage;
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
            "sc-faber-ic-nontop-{}-{}.sqlite",
            std::process::id(),
            NEXT_DATABASE.fetch_add(1, Ordering::Relaxed)
        )));
        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            assert_eq!(
                int_class(2, 2, &[(1, 2)], &[(1, 1)], &[(1, 1)]).unwrap(),
                Q::ZERO
            );
            assert_eq!(
                int_class(
                    2,
                    2,
                    &[(1, 1), (1, 1)],
                    &[(1, 0), (1, 1)],
                    &[(1, 0), (1, 1)]
                )
                .unwrap(),
                Q::ZERO
            );
            assert_eq!(int_class(1, 2, &[(1, 1)], &[], &[]).unwrap(), Q::ZERO);
            // KLP and the fundamental-monomial API share this same root key.
            assert_eq!(klp(1, 2, &[(1, 1)], &[], &[], &[]).unwrap(), Q::ZERO);
            assert_eq!(
                crate::int::int_fund_monomial(1, 2, &[(1, 1)], &[]).unwrap(),
                Q::ZERO
            );
            assert_eq!(int_class(0, 5, &[], &[], &[(1, 1)]).unwrap(), Q::ZERO);
        });
        let expected = [
            IntegralKey::Integral {
                g: 2,
                n: 2,
                kappa: vec![(1, 2)],
                psi: vec![(1, 1)],
                lambda: vec![],
                ch: vec![(1, 1)],
            },
            IntegralKey::Integral {
                g: 1,
                n: 2,
                kappa: vec![(1, 1)],
                psi: vec![],
                lambda: vec![],
                ch: vec![],
            },
            IntegralKey::Integral {
                g: 0,
                n: 5,
                kappa: vec![],
                psi: vec![],
                lambda: vec![],
                ch: vec![(1, 1)],
            },
        ];
        let sources = registry.sources();
        assert_eq!(sources.len(), 1);
        assert_eq!(sources[0].ready(), expected.len());
        let rows = sources[0].pull(0, usize::MAX).unwrap().rows;
        for (row, key) in rows.iter().zip(&expected) {
            assert_eq!(&row.key, key);
            assert_eq!(row.value, Q::ZERO);
        }
        let storage = Storage::open(database.0.to_str().unwrap(), registry).unwrap();
        assert_eq!(storage.close().unwrap().rows, expected.len());
        let reopened =
            Storage::open(database.0.to_str().unwrap(), Arc::new(Registry::new())).unwrap();
        for key in expected {
            assert_eq!(reopened.get(key).unwrap(), Some(Q::ZERO));
        }
        reopened.close().unwrap();
    }

    #[test]
    fn int_class_non_top_zero_precedence_does_not_publish_invalid_roots() {
        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            // These historically return zero before the top-degree path's
            // parameter checks. They must keep that behavior without caching
            // invalid input as a stable mathematical integral.
            assert_eq!(int_class(1, 0, &[], &[], &[(1, 1)]).unwrap(), Q::ZERO);
            assert_eq!(int_class(2, 1, &[(0, 1)], &[], &[]).unwrap(), Q::ZERO);
            assert_eq!(int_class(2, 1, &[], &[(2, 1)], &[]).unwrap(), Q::ZERO);
            assert_eq!(int_class(2, 1, &[], &[(0, 1)], &[]).unwrap(), Q::ZERO);
            assert_eq!(int_class(2, 1, &[], &[], &[(2, 1)]).unwrap(), Q::ZERO);
            assert_eq!(int_class(2, 1, &[], &[], &[(0, 1)]).unwrap(), Q::ZERO);
            assert!(registry.sources().is_empty());
            // Reaching top degree still follows the original error path.
            assert!(
                int_class(2, 1, &[(0, 1)], &[(1, 4)], &[])
                    .unwrap_err()
                    .contains("κ_0")
            );
            assert!(
                int_class(2, 1, &[], &[(1, 2)], &[(2, 1)])
                    .unwrap_err()
                    .contains("ch_2")
            );
            assert!(registry.sources().iter().all(|source| source.ready() == 0));
        });
    }

    #[test]
    fn ic_cache_context_owns_the_only_ic_memo_and_preserves_raw_key_boundaries() {
        let expected = with_optional_registry(None, || {
            IC_MEMO.with(|memo| memo.borrow_mut().clear());
            let value = int_class(2, 2, &[], &[(1, 2)], &[(1, 3)]).unwrap();
            IC_MEMO.with(|memo| memo.borrow_mut().clear());
            value
        });
        assert_eq!(expected, Q::from_signeds(1, 1440));
        let first = Arc::new(Registry::new());
        let key = (2, 2, vec![], vec![(1, 2)], vec![(1, 3)]);
        with_registry(first.clone(), || {
            assert_eq!(
                int_class(2, 2, &[], &[(1, 2)], &[(1, 3)]).unwrap(),
                expected
            );
            let source = ic_source(&first);
            assert_eq!(source.get(&key), Some(expected.clone()));
            let ready = source.ready();
            assert_eq!(
                int_class(2, 2, &[], &[(1, 1), (1, 1)], &[(1, 1), (1, 2)]).unwrap(),
                expected,
            );
            assert_eq!(source.ready(), ready);
            assert!(IC_MEMO.with(|memo| memo.borrow().is_empty()));
            let root = source
                .pull(0, usize::MAX)
                .unwrap()
                .rows
                .into_iter()
                .find(|row| {
                    row.key
                        == IntegralKey::Legacy {
                            kind: "ic".to_owned(),
                            nums: ic_nums(&key),
                        }
                })
                .unwrap();
            assert_eq!(root.value, expected);
            assert_eq!(ic_nums(&key), [2, 2, 0, 2, 2, 1, 2, 1, 3]);
        });
        let second = Arc::new(Registry::new());
        with_registry(second.clone(), || {
            assert_eq!(
                int_class(2, 2, &[], &[(1, 2)], &[(1, 3)]).unwrap(),
                expected
            );
            assert!(IC_MEMO.with(|memo| memo.borrow().is_empty()));
        });
        assert!(!Arc::ptr_eq(&ic_source(&first), &ic_source(&second)));
        first.fail("test IC storage failure".to_owned());
        with_registry(first, || {
            assert_eq!(
                int_class(2, 2, &[], &[(1, 2)], &[(1, 3)]).unwrap_err(),
                "test IC storage failure",
            );
        });
    }

    #[test]
    fn ic_history_hit_is_loaded_into_its_mathematical_source() {
        use crate::cache_source::storage::Storage;

        let registry = Arc::new(Registry::new());
        let persisted = Arc::new(LegacySource::new("ic", ic_nums));
        let key = (2, 2, vec![], vec![(1, 2)], vec![(1, 3)]);
        let expected = Q::from_signeds(1, 1440);
        persisted.insert(key.clone(), &expected);
        registry.register(persisted);
        let storage = Storage::open(":memory:", registry.clone()).unwrap();
        registry.bind_storage(storage.handle());
        storage.flush().unwrap();
        // Simulate an already populated backend with no IC source in this
        // calling context; production sets this flag only when opening a DB.
        registry.set_history(true);
        with_registry(registry.clone(), || {
            assert_eq!(
                int_class(2, 2, &[], &[(1, 2)], &[(1, 3)]).unwrap(),
                expected
            );
            let source = ic_source(&registry);
            assert_eq!(source.get(&key), Some(expected));
            // A recomputation would also complete recursive IC states.
            assert_eq!(source.ready(), 1);
        });
        storage.close().unwrap();
    }

    #[test]
    fn root_cache_cannot_bypass_klp_input_errors() {
        let registry = Arc::new(Registry::new());
        // Other internal APIs can publish a root before klp sees it. Even a
        // cached zero must not convert a rejected public input into Ok(zero).
        for (g, n, kappa, psi, ch) in [
            (2, 1, vec![(0, 1)], vec![(1, 4)], vec![]),
            (2, 1, vec![], vec![(1, 2)], vec![(2, 1)]),
            (1, 0, vec![], vec![], vec![]),
            (2, 1, vec![], vec![(2, 4)], vec![]),
        ] {
            root_source::with_registry(registry.clone(), g, n, &kappa, &psi, &[], &ch, || {
                Ok(Q::ZERO)
            })
            .unwrap();
            with_registry(registry.clone(), || {
                assert!(klp(g, n, &kappa, &psi, &[], &ch).is_err());
            });
        }
    }

    #[test]
    fn klp_route_preserves_existing_zero_and_error_precedence() {
        for registry in [None, Some(Arc::new(Registry::new()))] {
            with_optional_registry(registry, || {
                // Rank and g=0 Hodge vanishing precede kappa/ch validation.
                assert_eq!(
                    klp(2, 1, &[(0, 1)], &[], &[(3, 1)], &[(2, 1)]).unwrap(),
                    Q::ZERO
                );
                assert_eq!(klp(0, 4, &[(0, 1)], &[], &[], &[(2, 1)]).unwrap(), Q::ZERO);
                // Kappa zero is checked before weight; even ch after weight.
                assert!(klp(2, 1, &[(0, 1)], &[], &[], &[]).is_err());
                assert_eq!(klp(2, 1, &[], &[], &[], &[(2, 1)]).unwrap(), Q::ZERO);
                assert!(klp(2, 1, &[], &[(1, 2)], &[], &[(2, 1)]).is_err());
                // The first vanishing odd ch wins before a later even ch.
                assert_eq!(
                    klp(1, 7, &[], &[], &[], &[(3, 1), (4, 1)]).unwrap(),
                    Q::ZERO
                );
            });
        }
    }

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
