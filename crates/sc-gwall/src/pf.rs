//! local CY3 配分函数层（M4）：branch → Ztil → Z → Y + BPS 反演。
//!
//! branch(g,d,X) = ∫ c_X(Rπf*(L⊕L))（deg L = g−1）——chernnum 的数学：
//!   c_X(E⊕E) = Σ_{λ⊢X} [Π_k 2(−1)^{λ_k}(λ_k−1)! / permfact(λ)]
//!              · Π_k ch_{λ_k}(Rπf*L)
//! （与原版手写 branch2 = 2ch₂ + 2ch₁² 精确吻合，故运行时分拆和生成，
//! 不抄 branch2..10 的手写展开）。
//! Ztil = gw(g,d,[]) + Σ branch(2k)t^{2k}；Z = exp(−ΣN₀ t^{2j})·Ztil（去
//! 常映射）；Y = t^{2d(g−1)}·Z|_{t=2arcsin(√−y/2)}。
//! Z 的 t^{2k} 系数 = reduced 不变量 N^h_d(g)，h = k + (g−1)(1−d)
//! （vdim = 2(g+h−1) − 2(g−1)d 反解；g=1/2 各例自洽验证）。
//! BPS 反演（BP 论文 Prop general inversion）：
//!   n_d^h(g) = Σ_{h'≤h} Σ_{k|d} μ(k)·k^{2(g+h')−3}·α_{h+g,h'+g}·N^{h'}_{d/k}(g)
//!   α_{G,G'} = [r^{G−G'}] A(r)^{2G'−2}，A(r) = arcsin(√r/2)/(√r/2)
//! （A 同时是 Y 变换的核——GWall 的 y 变量设计正为对齐此反演）。

use crate::comb::partition;
use crate::gwl::{Ctx, Ins, gw_ctx};
use crate::gwpoint::n0;
use crate::qseries::QSeries;
use sc_num::{One, Q, Zero};

/// permfact(λ) = Π_m (m 在 λ 中的重数)!（分拆自同构阶）。
fn permfact(lambda: &[u32]) -> Q {
    let d: u32 = lambda.iter().sum();
    let mut counts = vec![0u32; d as usize + 1];
    for &x in lambda {
        counts[x as usize] += 1;
    }
    let mut v = Q::ONE;
    for &c in &counts {
        for i in 2..=c {
            v *= Q::from(i as i64);
        }
    }
    v
}

/// branch(g,d,X) = ∫ c_X(Rπf*(L⊕L))，deg L = g−1（chernnum 公式化）。
pub fn branch(ctx: &mut Ctx, g: u32, d: u32, x: u32) -> Q {
    if x == 0 {
        return gw_ctx(ctx, g, d, &[]);
    }
    let mut answer = Q::ZERO;
    for lambda in partition(x) {
        let mut coeff = Q::ONE / permfact(&lambda);
        for &lk in &lambda {
            // 2·(−1)^λ·(λ−1)!（合并 (−O−O) = −2O 与 (−1)^{λ−1}）；
            // (λ−1)! = Π_{i=2..λ−1} i（λ=1 为 1）
            let mut c = Q::from(2) * Q::from(if lk % 2 == 0 { 1 } else { -1 });
            for i in 2..lk {
                c *= Q::from(i as i64);
            }
            coeff *= c;
        }
        let mut ins: Vec<Ins> = lambda.iter().map(|&lk| Ins::ChL(lk)).collect();
        ins.push(Ins::LDeg(g as i64 - 1));
        answer += coeff * gw_ctx(ctx, g, d, &ins);
    }
    answer
}

/// Ztil(g,d,N) = gw(g,d,[]) + Σ_{k=1..N/2} branch(2k)·t^{2k}（t 截到 N）。
pub fn ztil(ctx: &mut Ctx, g: u32, d: u32, n: u32) -> Result<QSeries, String> {
    if n > 10 {
        return Err("branch 上限 10（原版 'No can do more than branch number 10'）".into());
    }
    let mut c = vec![Q::ZERO; n as usize + 1];
    c[0] = gw_ctx(ctx, g, d, &[]);
    let mut k = 1u32;
    while 2 * k <= n {
        c[2 * k as usize] = branch(ctx, g, d, 2 * k);
        k += 1;
    }
    Ok(QSeries { c })
}

/// Z(g,d,N) = exp(−Σ_{j≥1} N₀(g,2j)·t^{2j})·Ztil，截断 t^{N+1}。
/// 带 (g,d) 级缓存（截断不足时以更大的 n 重算并覆盖）。
pub fn z(ctx: &mut Ctx, g: u32, d: u32, n: u32) -> Result<QSeries, String> {
    if n > 10 && d > 1 {
        // branch 上限 10 仅约束 d≥2（度 1 走 Pandharipande 闭式无限制）
        return Err("branch 上限 10".into());
    }
    if let Some((c, n0_)) = ctx.zmem.get(&(g, d)) {
        if *n0_ >= n {
            return Ok(QSeries {
                c: c.iter().take(n as usize + 1).cloned().collect(),
            });
        }
    }
    // 度 1 闭式（Pandharipande，度 1 锚已全系数验证）：
    // Z_1 = (4sin²(t/2)/t²)^{g−1}，4sin²(t/2)/t² = Σ (−1)^k·2·t^{2k}/(2k+2)!
    // g = 0 时为负幂（级数求逆）。结果必须保持完整长度 n+1：QSeries::mul
    // 按较短因子截断，短级数会把连通化 L_d 的高次系数截掉。
    if d == 1 {
        let mut s = vec![Q::ZERO; n as usize + 1];
        for k in 0..=(n / 2) as usize {
            let mut f = Q::ONE;
            for i in 2..=(2 * k + 2) {
                f *= Q::from(i as i64);
            }
            s[2 * k] = Q::from(if k % 2 == 0 { 2 } else { -2 }) / f;
        }
        let z1 = QSeries { c: s }.pow_signed(g as i64 - 1)?;
        ctx.zmem.insert((g, d), (z1.c.clone(), n));
        return Ok(z1);
    }
    let mut n0s = vec![Q::ZERO; n as usize + 2];
    let mut j = 1usize;
    while 2 * j < n as usize + 2 {
        n0s[2 * j] = Q::ZERO - n0(g, 2 * j as u32);
        j += 1;
    }
    let e = QSeries { c: n0s }.exp().unwrap();
    let full = e.mul(&ztil(ctx, g, d, n)?).trunc(n + 1);
    ctx.zmem.insert((g, d), (full.c.clone(), n));
    Ok(full)
}

/// A(r) = arcsin(√r/2)/(√r/2)，截断阶 n（偶级数：只含 r^m）。
fn arcsin_over(n: u32) -> QSeries {
    // arcsin(w)/w = Σ_m a_m w^{2m}，a_m = C(2m,m)/(4^m(2m+1))
    let mut c = vec![Q::ZERO; n as usize + 1];
    for (m, cm) in c.iter_mut().enumerate() {
        // C(2m, m)
        let mut bin = Q::ONE;
        for i in 0..m {
            bin = bin * Q::from((2 * m - i) as i64) / Q::from((m - i) as i64);
        }
        let mut p4 = Q::ONE;
        for _ in 0..m {
            p4 *= Q::from(4);
        }
        *cm = bin / p4 / Q::from((2 * m + 1) as i64);
    }
    QSeries { c }
}

/// Y(g,d,N)：y = −4sin²(t/2) 变量下的重整化配分函数（截断 d(g−1)+N）。
/// u = t² = (−y)·A(−y/4)²，Y = u^{d(g−1)}·Z̃(u)|_{u←(−y)A(−y/4)²}。
/// g = 0 时 u 的幂为负（Laurent 级数），不支持。
pub fn y_series(ctx: &mut Ctx, g: u32, d: u32, n: u32) -> Result<QSeries, String> {
    if g == 0 {
        return Err("y_series 需要 g ≥ 1（g = 0 时为 Laurent 级数）".into());
    }
    let order = d * (g - 1) + n;
    // Z 的偶系数 → u 级数
    let zt = z(ctx, g, d, n)?;
    let mut zu = vec![Q::ZERO; (n / 2) as usize + 1];
    for (i, c) in zt.c.iter().enumerate() {
        if i % 2 == 0 && i / 2 < zu.len() {
            zu[i / 2] = c.clone();
        }
    }
    let mut su = QSeries { c: zu };
    // 乘 u^{d(g−1)}（左移）
    let shift = (d * (g - 1)) as usize;
    if shift > 0 {
        let mut c = vec![Q::ZERO; su.c.len() + shift];
        for (i, v) in su.c.iter().enumerate() {
            c[i + shift] = v.clone();
        }
        su = QSeries { c };
    }
    // u ← h(y) = (−y)·A(−y/4)²：A 在 r 变量截 order，代 r = −y/4
    let a = arcsin_over(order);
    let a2 = a.powi(2);
    // A(−y/4)：r ← −y/4（奇次为 0，m 次乘 (−1/4)^m）
    let mut af = vec![Q::ZERO; order as usize + 1];
    for (m, c) in a2.c.iter().enumerate() {
        if m <= order as usize {
            let mut s = Q::ONE;
            for _ in 0..m {
                s = s * Q::from(-1) / Q::from(4);
            }
            af[m] = c.clone() * s;
        }
    }
    let a2f = QSeries { c: af };
    // h(y) = (−y)·A2(y)
    let mut h = vec![Q::ZERO; order as usize + 1];
    for (hi, a) in h.iter_mut().skip(1).zip(&a2f.c) {
        *hi = -a.clone();
    }
    su.compose(&QSeries { c: h })
}

/// 连通化：L_d = [q^d]log(1+Σ Z_e q^e)，递推
/// L_d = Z_d − (1/d)Σ_{j=1..d−1} j·L_j·Z_{d−j}（指数公式：连通分量）。
/// 论文表的 N^h_d(g) 从 L_d 读（M4 悬案根因之二：初版直读 disconnected
/// 的 Z 系数——n_2^1(2) 因误差相消侥幸通过，g≥3 暴露）。
pub fn l_series(ctx: &mut Ctx, g: u32, d: u32, n: u32) -> Result<QSeries, String> {
    // 逐度递推（t 截断 n）
    let mut l: Vec<QSeries> = Vec::with_capacity(d as usize + 1);
    for dd in 0..=d {
        if dd == 0 {
            l.push(QSeries::constant(Q::ONE, n));
            continue;
        }
        let zd = z(ctx, g, dd, n)?;
        // L_d = Z_d − (1/d)Σ_{j=1..d−1} j·L_j·Z_{d−j}（Z_d 不除 d）
        let mut s = QSeries::constant(Q::ZERO, n);
        for j in 1..dd {
            let lj = l[j as usize].clone();
            let zdj = z(ctx, g, dd - j, n)?;
            s = s.add(&lj.mul(&zdj).scale(&Q::from(j as i64)));
        }
        let cf = Q::from(-1) / Q::from(dd as i64);
        l.push(zd.add(&s.scale(&cf)));
    }
    Ok(l[d as usize].clone())
}

/// N^h_d(g)：论文表的（连通、reduced）不变量 = L_d 的 t^{2(h−(d−1)(g−1))} 系数。
pub fn n_reduced(ctx: &mut Ctx, g: u32, d: u32, h: u32) -> Result<Q, String> {
    let kk = h as i64 - (d as i64 - 1) * (g as i64 - 1);
    if kk < 0 {
        return Ok(Q::ZERO);
    }
    let idx = 2 * kk as u32;
    if idx > 10 && d > 1 {
        return Err(format!(
            "N^{h}_{d}(g={g}) 需要 L 截断 t^{idx}，超出 branch 上限 10（原版同限）"
        ));
    }
    let lt = l_series(ctx, g, d, idx)?;
    Ok(lt.c.get(idx as usize).cloned().unwrap_or(Q::ZERO))
}

/// Möbius μ(k)。
fn moebius(mut k: u32) -> i64 {
    if k == 1 {
        return 1;
    }
    let mut res: i64 = 1;
    let mut p = 2u32;
    while p * p <= k {
        if k % p == 0 {
            k /= p;
            if k % p == 0 {
                return 0;
            }
            res = -res;
        }
        p += 1;
    }
    if k > 1 {
        res = -res;
    }
    res
}

/// α_{G,G'} = [r^{G−G'}] A(r)^{2G'−2}，A(r) = arcsin(√r/2)/(√r/2)。
/// 注意 A(r) = B(r/4)（B = arcsin(√r)/√r）——M4 悬案根因之一：
/// 初版误用 B，α 偏差 4^{G−G'}（G=G' 时不可见，恰使 h=0 行全对的假象）。
pub fn alpha(big: u32, small: u32) -> Q {
    let e = 2 * small as i64 - 2; // A 的幂（small = 0 即 g = 0 时为 −2）
    let n = big.saturating_sub(small);
    let b = arcsin_over(n);
    // A 的系数 = B_m / 4^m
    let mut af = vec![Q::ZERO; n as usize + 1];
    for (m, c) in b.c.iter().enumerate() {
        let mut p4 = Q::ONE;
        for _ in 0..m {
            p4 *= Q::from(4);
        }
        af[m] = c.clone() / p4;
    }
    // A 的常数项为 1，负幂恒可逆
    let ap = QSeries { c: af }.pow_signed(e).expect("A(r) 常数项为 1");
    ap.c.get(n as usize).cloned().unwrap_or(Q::ZERO)
}

/// local BPS：n_d^h(g)（论文 eq inverted GV）。
/// 所需 Z 截断逐项恰取（n_reduced 内部控制，≤ branch 上限 10）。
pub fn bps_local(ctx: &mut Ctx, g: u32, d: u32, h: u32) -> Result<Q, String> {
    // 先检查截断上限，免得算完低阶项才报错：下面每个用到的 d/k > 1 在
    // h' = h 处需要的截断最高（与 n_reduced 的判据一致）
    for k in (1..=d).filter(|k| d % *k == 0 && moebius(*k) != 0 && d / k > 1) {
        let dd = d / k;
        let kk = h as i64 - (dd as i64 - 1) * (g as i64 - 1);
        if 2 * kk > 10 {
            return Err(format!(
                "n^{h}_{d}(g={g}) 需要 N^{h}_{dd} 的 L 截断 t^{}，超出 branch 上限 10（原版同限）",
                2 * kk
            ));
        }
    }
    let mut answer = Q::ZERO;
    for k in 1..=d {
        if d % k != 0 {
            continue;
        }
        let mu = moebius(k);
        if mu == 0 {
            continue;
        }
        let dd = d / k;
        for hh in 0..=h {
            let ncoef = n_reduced(ctx, g, dd, hh)?;
            // m^{2(g+h')−3}，指数可为负（g=1, h'=0 → m⁻¹）
            let e = 2 * (g + hh) as i64 - 3;
            let mut kp = Q::ONE;
            if e >= 0 {
                for _ in 0..e {
                    kp *= Q::from(k as i64);
                }
            } else {
                for _ in 0..(-e) {
                    kp /= Q::from(k as i64);
                }
            }
            answer += Q::from(mu) * kp * alpha(h + g, hh + g) * ncoef;
        }
    }
    Ok(answer)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// branch2 对照原版手写展开：2ch₂ + 2ch₁²（ch = ChL，deg = g−1）。
    #[test]
    fn branch2_formula() {
        let mut ctx = Ctx::new();
        for (g, d) in [(2u32, 1u32), (2, 2), (3, 1)] {
            let b = branch(&mut ctx, g, d, 2);
            let hand = gw_ctx(&mut ctx, g, d, &[Ins::ChL(2), Ins::LDeg(g as i64 - 1)]) * Q::from(2)
                + gw_ctx(
                    &mut ctx,
                    g,
                    d,
                    &[Ins::ChL(1), Ins::ChL(1), Ins::LDeg(g as i64 - 1)],
                ) * Q::from(2);
            assert_eq!(b.to_string(), hand.to_string(), "g={g} d={d}");
        }
    }

    /// BP 论文 local BPS 表 n_3/n_4/n_5(g)（tex 1870-1915 行逐格核对；
    /// d=5 g=2 的 h≥4 格含论文未知参数 a_{5,2}/*，跳过）。
    #[test]
    fn bps_d345_table() {
        let mut ctx = Ctx::new();
        #[rustfmt::skip]
        let table: &[(u32, u32, u32, i64)] = &[
            // (d, g, h, want)——每行对应论文表的一行 (d, g)
            (3, 2, 0, -3), (3, 2, 1, 2), (3, 2, 2, 73), (3, 2, 3, 50),
            (3, 3, 0, -27), (3, 3, 1, 36), (3, 3, 2, -18), (3, 3, 3, 4),
            (3, 3, 4, 2641), (3, 3, 5, 9604),
            (3, 4, 0, -243), (3, 4, 1, 486), (3, 4, 2, -405), (3, 4, 3, 180),
            (3, 4, 4, -45), (3, 4, 5, 6), (3, 4, 7, 836310),
            // (3,4,6) 见 bps_paper_erratum（论文勘误，真值 93913）
            (3, 5, 0, -2187), (3, 5, 1, 5832), (3, 5, 2, -6804), (3, 5, 3, 4536),
            (3, 5, 4, -1890), (3, 5, 5, 504), (3, 5, 6, -84), (3, 5, 7, 8),
            (3, 6, 0, -19683), (3, 6, 1, 65610), (3, 6, 2, -98415), (3, 6, 3, 87480),
            (3, 6, 4, -51030), (3, 6, 5, 20412), (3, 6, 6, -5670), (3, 6, 7, 1080),
            (3, 7, 0, -177147), (3, 7, 1, 708588), (3, 7, 2, -1299078), (3, 7, 3, 1443420),
            (3, 7, 4, -1082565), (3, 7, 5, 577368), (3, 7, 6, -224532), (3, 7, 7, 64152),
            (4, 2, 0, 0), (4, 2, 1, -60), (4, 2, 2, 30),
            (4, 3, 0, 0), (4, 3, 1, 0), (4, 3, 2, -4032),
            (4, 4, 0, 0), (4, 4, 1, 0), (4, 4, 2, 0), (4, 4, 3, -261120),
            (4, 5, 0, 0), (4, 5, 1, 0), (4, 5, 2, 0), (4, 5, 3, 0), (4, 5, 4, -16760832),
            // 行尾格 (4,3,3)/(4,4,4)/(4,5,5) 见 bps_paper_erratum（论文勘误）
            (5, 2, 0, -5), (5, 2, 1, 10), (5, 2, 2, -7), (5, 2, 3, 2),
            (5, 3, 0, -125), (5, 3, 1, 500), (5, 3, 2, -850), (5, 3, 3, 800),
            (5, 3, 4, -455), (5, 3, 5, 160), (5, 3, 6, -34), (5, 3, 7, 4),
            (5, 4, 0, -3125), (5, 4, 1, 18750), (5, 4, 2, -50625), (5, 4, 3, 81250),
            (5, 4, 4, -86250), (5, 4, 5, 63750), (5, 4, 6, -33625), (5, 4, 7, 12750),
            (5, 5, 0, -78125), (5, 5, 1, 625000), (5, 5, 2, -2312500), (5, 5, 3, 5250000),
            (5, 5, 4, -8181250), (5, 5, 5, 9275000), (5, 5, 6, -7910000), (5, 5, 7, 5175000),
        ];
        for &(d, g, h, want) in table {
            let got = bps_local(&mut ctx, g, d, h).unwrap_or_else(|e| panic!("{e}"));
            assert_eq!(got.to_string(), want.to_string(), "n_{d}^{h}(g={g})");
        }
    }

    /// 论文勘误（codex 仲裁 + 本方数值验收）：4 格以我们的全量 GW 反演
    /// 为准——论文 local 表的印值有误：
    /// 1. n_3^6(4)：étale 界引理（tex 491-502）h ≤ (d_min−1)(g−1) = 6 时
    ///    local = étale；独立复算 étale 生成函数给 93913（表印 692352
    ///    与论文自身引理矛盾）；
    /// 2. n_4 三个行尾格：tex 622 行 GRR 展开的 ε 指数印刷错误——
    ///    正确 GV 权重 μ(m)m^{2d_min(g−1)+1}，论文印 m^{d_min(g−1)+2}；
    ///    正确指数复算得 −64、−394752、−67141632（= 我们的值；表印值
    ///    是错误指数的产物）。d=2 时 m=1 两指数相同——n₂ 表全对的原因。
    #[test]
    fn bps_paper_erratum() {
        let mut ctx = Ctx::new();
        for &(d, g, h, want) in &[
            (3u32, 4u32, 6u32, 93913i64),
            (4, 3, 3, -64),
            (4, 4, 4, -394752),
            (4, 5, 5, -67141632),
        ] {
            let got = bps_local(&mut ctx, g, d, h).unwrap();
            assert_eq!(got.to_string(), want.to_string(), "n_{d}^{h}(g={g})");
        }
    }

    /// A(r) 首系数：arcsin(w)/w = 1 + w²/6 + 3w⁴/40 + ...
    #[test]
    fn arcsin_series() {
        let a = arcsin_over(3);
        assert_eq!(a.c[0].to_string(), "1");
        assert_eq!(a.c[1].to_string(), "1/6");
        assert_eq!(a.c[2].to_string(), "3/40");
        assert_eq!(a.c[3].to_string(), "5/112");
    }

    /// 度 1 的 Y 有精确闭式：u^{g−1}·Z_1 = (4sin²(t/2))^{g−1} = (−y)^{g−1}
    /// （y = −4sin²(t/2)）。同时检验 y 代换核 h(y) 与度 1 闭式。
    #[test]
    fn y_series_degree1_closed_form() {
        let mut ctx = Ctx::new();
        let y2 = y_series(&mut ctx, 2, 1, 4).unwrap();
        let s2: Vec<String> = y2.c.iter().map(|v| v.to_string()).collect();
        assert_eq!(s2, ["0", "-1", "0", "0"]);
        let y3 = y_series(&mut ctx, 3, 1, 4).unwrap();
        let s3: Vec<String> = y3.c.iter().map(|v| v.to_string()).collect();
        assert_eq!(s3, ["0", "0", "1", "0", "0"]);
        assert!(y_series(&mut ctx, 0, 1, 4).is_err());
    }

    /// BP 论文 §1（tex 355-367 行）引述的完整结果：
    /// g = 0（Faber–Pandharipande）：n_d^h(0) = 1 当且仅当 d = 1、h = 0；
    /// g = 1（Pandharipande）：n_d^h(1) = 1 当且仅当 h = 0（d ≥ 1）。
    /// g = 0 走负幂分支（α 中 A^{−2}、度 1 闭式求逆）；g = 1 的 d ≥ 3 检验
    /// 度 1 闭式保持完整截断长度（否则连通化 L_d 丢高次系数）。
    #[test]
    fn bps_genus0_and_genus1() {
        let mut ctx = Ctx::new();
        for d in 1..=4u32 {
            for h in 0..=3u32 {
                // g = 0 时所需截断 t^{2(h+d−1)}，受 branch ≤ 10 限制
                if h + d <= 6 {
                    let want0 = if d == 1 && h == 0 { 1 } else { 0 };
                    let got0 = bps_local(&mut ctx, 0, d, h).unwrap();
                    assert_eq!(got0.to_string(), want0.to_string(), "n_{d}^{h}(g=0)");
                } else {
                    assert!(bps_local(&mut ctx, 0, d, h).is_err());
                }
                let want1 = if h == 0 { 1 } else { 0 };
                let got1 = bps_local(&mut ctx, 1, d, h).unwrap();
                assert_eq!(got1.to_string(), want1.to_string(), "n_{d}^{h}(g=1)");
            }
        }
    }

    /// 度 1 全系锚（Pandharipande 闭式，BP 论文 eq）：
    /// Σ_h N^h_1(g)·t^{2(g+h−1)} = (4sin²(t/2))^{g−1}，
    /// 4sin²(t/2) = Σ 2(−1)^{k+1}t^{2k}/(2k)!。
    #[test]
    fn degree1_anchor() {
        let mut ctx = Ctx::new();
        for g in [2u32, 3, 4] {
            // (4sin²(t/2))^{g−1} 的系数
            let n = 8usize;
            let mut s2 = vec![Q::ZERO; n + 1];
            for k in 1..=(n / 2) {
                let mut f = Q::ONE;
                for i in 2..=(2 * k) {
                    f *= Q::from(i as i64);
                }
                s2[2 * k] = Q::from(if k % 2 == 1 { 2 } else { -2 }) / f;
            }
            let mut want = QSeries { c: vec![Q::ONE] };
            for _ in 0..g - 1 {
                let mut nc = vec![Q::ZERO; n + 1];
                for (i, a) in want.c.iter().enumerate() {
                    for (j, b) in s2.iter().enumerate() {
                        if i + j <= n {
                            nc[i + j] = nc[i + j].clone() + a.clone() * b.clone();
                        }
                    }
                }
                want = QSeries { c: nc };
            }
            for h in 0..=((n / 2).saturating_sub(g as usize - 1)).min(3) as u32 {
                let got = n_reduced(&mut ctx, g, 1, h).unwrap();
                let w = want
                    .c
                    .get(2 * (g as usize - 1 + h as usize))
                    .cloned()
                    .unwrap_or(Q::ZERO);
                assert_eq!(got.to_string(), w.to_string(), "N^{h}_1(g={g})");
            }
        }
    }

    /// BPS 对拍：BP 论文 local 表 n_2^h(g)（g=2..7 行）。
    #[test]
    fn bps_d2_table() {
        let mut ctx = Ctx::new();
        // BP 论文 local 表 n_2^h(g) 全部数值格（tex 1858-1864 行）。
        // 注意 local 表（第二块）每行比 étale 表（第一块）多一个非零格
        // （h = (d−1)(g−1) 处），勿混抄。
        #[rustfmt::skip]
        let table: &[(u32, u32, i64)] = &[
            // (g, h, want)——每行对应论文表的一行 g
            (2, 0, -2), (2, 1, 8), (2, 2, 0),
            (3, 0, -8), (3, 1, 4), (3, 2, 31), (3, 3, 8),
            (4, 0, -32), (4, 1, 24), (4, 2, -6), (4, 3, 128), (4, 4, 96),
            (5, 0, -128), (5, 1, 128), (5, 2, -48), (5, 3, 8), (5, 4, 511), (5, 5, 768),
            (6, 0, -512), (6, 1, 640), (6, 2, -320), (6, 3, 80), (6, 4, -10), (6, 5, 2048),
            (6, 6, 5120),
            (7, 0, -2048), (7, 1, 3072), (7, 2, -1920), (7, 3, 640), (7, 4, -120), (7, 5, 12),
            (7, 6, 8191), (7, 7, 30720),
        ];
        for &(g, h, want) in table {
            let got = bps_local(&mut ctx, g, 2, h).unwrap_or_else(|e| panic!("{e}"));
            assert_eq!(got.to_string(), want.to_string(), "n_2^{h}(g={g})");
        }
    }
}
