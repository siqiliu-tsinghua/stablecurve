//! Virasoro 条件——原版 `libgw/virasoro.cc` 的直译。
//!
//! 对亏格 g ≥ 1、末点 ψ < 3g−1 的 Pⁿ 绝对不变量：取 L_k 条件的
//! (A)+(B)+(C)+(D) 项（(E) 项因 k ≥ 1 不出现），组装成线性方程，
//! 逐步加方程直至唯一可解。未知不变量编码为变量名 `(xdim+1)·d + Hpow`。
//! A/D 项经 [`crate::trr::trr_variable`] 把算不出的子不变量留作变量。
//!
//! `VIRASORO_TRY`/`VIRASORO_CHECK` 取原版 config.h 的值（5/0）。

use crate::ctx::Ctx;
use crate::inv::{Inv, Pdata};
use crate::linsys::{Linsys, LinsysT};
use crate::num::{One, Q, Zero, q};
use crate::split2::idiv;
use crate::symfunc::symfunc;

use crate::trr::trr_variable;

const VIRASORO_TRY: i32 = 5;

fn vdbg(msg: &str) {
    if crate::env_flag("SC_GROWI_VIRASORO") {
        eprintln!("[vir] {msg}");
    }
}

const VIRASORO_CHECK: i32 = 0;

/// A 项：新点插在最前，经 trr_variable 处理（可产生变量）。
fn virasoro_a(ctx: &mut Ctx, sys: &mut Linsys, iv: &Inv, k: i32, d: i32) -> Option<Q> {
    let mut iw_pts = vec![Pdata::default()];
    iw_pts.extend(iv.pts.iter().copied());
    iw_pts[0].mult = 0;
    let mut iw = Inv::new(iv.dim, iv.deg, iv.genus, d as i16, iw_pts);

    let mut r = Q::ONE;
    let mut result = Q::ZERO;
    for i in 0..=(k + 1).min(iv.dim as i32) {
        iw.pts[0] = Pdata::new(i as i16, (k + 1 - i) as i16, 0);
        let mut term = -&r * &symfunc(&(q(3 - iv.dim as i64) / Q::from(2)), k, i);
        r *= Q::from((iv.dim as i32 + 2 - iv.deg as i32) as i64);
        if crate::num::is_zero(&term) {
            continue;
        }
        vdbg(&format!("A need term {term}*{iw}"));
        term = trr_variable(ctx, sys, &iw, &term)?;
        vdbg(&format!("A got {term}"));
        result += term;
    }
    Some(result)
}

/// B 项：逐点抬 H 降 ψ（同类点去重计 cnt）。
fn virasoro_b(ctx: &mut Ctx, iv: &Inv, k: i32, d: i32) -> Option<Q> {
    let mut iw = iv.clone();
    iw.cdeg = d as i16;
    let n = iw.pts.len();

    let mut result = Q::ZERO;
    let mut cnt: i64 = 1;
    for j in 0..n {
        if j < n - 1 && iw.pts[j] == iw.pts[j + 1] {
            cnt += 1;
            continue;
        }
        let a = iw.pts[j].h;
        let m = iw.pts[j].psi;
        let mut r = Q::ONE;
        for i in 0..=(k + 1).min(iv.dim as i32 - a as i32) {
            if i > 0 {
                r *= Q::from((iv.dim as i32 + 2 - iv.deg as i32) as i64);
            }
            iw.pts[j].h = a + i as i16;
            let newpsi = k + m as i32 - i;
            if newpsi >= 0 {
                iw.pts[j].psi = newpsi as i16;
                let sf = symfunc(
                    &(q(1 - iv.dim as i64) / Q::from(2) + Q::from(a as i64 + m as i64)),
                    k,
                    i,
                );
                let mut term = &r * &Q::from(cnt) * &sf;
                if crate::num::is_zero(&term) {
                    continue;
                }
                vdbg(&format!("B need term {term}*{iw}"));
                term *= ctx.compute(&iw)?;
                result += term;
            }
        }
        iw.pts[j].h = a;
        iw.pts[j].psi = m;
        cnt = 1;
    }
    Some(result)
}

/// C 项： genus−1、前插两点（利用两点对称去重）。
fn virasoro_c(ctx: &mut Ctx, iv: &Inv, k: i32, d: i32) -> Option<Q> {
    let mut iw_pts = vec![Pdata::default(); 2];
    iw_pts.extend(iv.pts.iter().copied());
    iw_pts[0].mult = 0;
    iw_pts[1].mult = 0;
    let mut iw = Inv::new(iv.dim, iv.deg, iv.genus - 1, d as i16, iw_pts);

    let mut r = Q::ONE;
    let mut result = Q::ZERO;
    for i in 0..=(k + 1).min(iv.dim as i32) {
        for m in idiv(i - k - 1, 2)..=-1 {
            for a in i..=iv.dim as i32 {
                iw.pts[0] = Pdata::new(a as i16, (-m - 1) as i16, 0);
                iw.pts[1] = Pdata::new((iv.dim as i32 + i - a) as i16, (m + k - i) as i16, 0);
                // 两点对称：只取 p0 ≤ p1
                if iw.pts[0] > iw.pts[1] {
                    continue;
                }
                let mut term = &r
                    * &symfunc(
                        &(q(iv.dim as i64 + 1) / Q::from(2) - Q::from(a as i64 - m as i64)),
                        k,
                        i,
                    );
                if m % 2 != 0 {
                    term = -term;
                }
                if iw.pts[0] == iw.pts[1] {
                    term /= Q::from(2);
                }
                if crate::num::is_zero(&term) {
                    continue;
                }
                vdbg(&format!("B need term {term}*{iw}"));
                term *= ctx.compute(&iw)?;
                result += term;
            }
        }
        r *= Q::from((iv.dim as i32 + 2 - iv.deg as i32) as i64);
    }
    Some(result)
}

/// D 项：亏格分裂 h ⊔ g−h，经 trr_variable（首类分裂保留变量）。
fn virasoro_d(ctx: &mut Ctx, sys: &mut Linsys, iv: &Inv, k: i32, d: i32) -> Option<Q> {
    let (dim, deg) = (iv.dim as i32, iv.deg as i32);
    let dm = dim + 2 - deg;

    // ar[0] 骨架：[新点(兼对角), iv 各点]；ar[1]：[对角点]
    let mut v_pts = vec![Pdata::default()];
    v_pts.extend(iv.pts.iter().copied());
    v_pts[0].mult = 0;
    let v_hdr = Inv::new(iv.dim, iv.deg, iv.genus, iv.cdeg, v_pts);
    let w_hdr = Inv::new(iv.dim, iv.deg, 0, 0, [Pdata::default()]);

    let mut result = Q::ZERO;
    let mut coeff = Q::ONE / Q::from(deg as i64);
    let mut ddw: i32 = dim - 3;
    let mut skip = true;
    let mut sp = crate::split2::PairSplit::new(&v_hdr, &w_hdr, 1, 1);

    loop {
        // 亏格分裂（利用 h ↔ g−h 对称只取一半）
        for g in 0..=iv.genus as i32 / 2 {
            let mut r = Q::ONE;
            for i in 0..=(k + 1).min(dim) {
                for m in (i - k)..=-1 {
                    let dw = ddw - (dim - 3) * g - m - k + i;
                    let d1 = idiv(dim - deg - dw + i, dm).max(0);
                    let d2 = idiv(dim - 1 - dw, dm).min(d);

                    for wd in d1..=d2 {
                        let hw = wd * dm + dw + 1;
                        let hv = dim + i - hw;

                        let (mut av, mut aw) = sp.materialize(
                            &v_hdr,
                            &w_hdr,
                            (d - wd) as i16,
                            wd as i16,
                            hv as i16,
                            hw as i16,
                            0,
                        );
                        av.genus = iv.genus - g as i16;
                        aw.genus = g as i16;
                        av.pts[0].psi = (-m - 1) as i16;
                        aw.pts[0].psi = (m + k - i) as i16;

                        let mut fact = &r
                            * &symfunc(
                                &(q(dim as i64 + 1) / Q::from(2) - Q::from(hv as i64 - m as i64)),
                                k,
                                i,
                            );
                        if 2 * g == iv.genus as i32 {
                            fact /= Q::from(2);
                        }
                        if m % 2 != 0 {
                            fact = -fact;
                        }
                        if crate::env_flag("SC_GROWI_VIRASORO") {
                            eprintln!(
                                "[vir] D g={g} i={i} m={m} wd={wd} hv={hv} hw={hw} fact={fact}"
                            );
                        }
                        if crate::num::is_zero(&fact) {
                            continue;
                        }
                        vdbg(&format!("D need term {fact}*{coeff}*{av}*{aw}"));

                        // 原版先算 w 侧
                        let mut tmpw = ctx.compute(&aw)?;
                        if crate::num::is_zero(&tmpw) {
                            continue;
                        }
                        tmpw *= &(&fact * &coeff);

                        if skip {
                            // 全亏格 + 全部标记点：走 trr_variable（留变量）
                            let tmpv = trr_variable(ctx, sys, &av, &tmpw)?;
                            result += tmpv;
                            continue;
                        }

                        let tmpv = ctx.compute(&av)?;
                        if crate::num::is_zero(&tmpv) {
                            continue;
                        }
                        result += &tmpv * &tmpw;
                    }
                }
                r *= Q::from(dm as i64);
            }
            skip = false;
        }
        if !sp.next(&mut ddw, &mut coeff) {
            break;
        }
    }

    Some(result)
}

/// 组装一条 L_k 方程并提交。各项算不出 → 不定（原版 UNDETERMINED）。
fn virasoro_k(ctx: &mut Ctx, sys: &mut Linsys, iv: &Inv, k: i32, d: i32) -> LinsysT {
    let va = virasoro_a(ctx, sys, iv, k, d);
    let vb = virasoro_b(ctx, iv, k, d);
    let vc = virasoro_c(ctx, iv, k, d);
    let vd = virasoro_d(ctx, sys, iv, k, d);
    if crate::env_flag("SC_GROWI_VIRASORO") {
        eprintln!("[vir] k={k} d={d}: A={va:?} B={vb:?} C={vc:?} D={vd:?}");
    }
    let x = match (va, vb, vc, vd) {
        (Some(a), Some(b), Some(c), Some(dd)) => a + b + c + dd,
        _ => return LinsysT::Undetermined,
    };
    sys.add_const_term(x);
    if crate::env_flag("SC_GROWI_VIRASORO") {
        eprintln!("[vir] row: {}", sys.debug_current_row());
    }
    sys.send_row()
}

/// Virasoro 入口。原版 `virasoro`。
pub fn virasoro(ctx: &mut Ctx, iv: &mut Inv) -> Option<Q> {
    // 必须有标记点；相对/原始类/非正亏格/非 Pⁱ/ψ 过多均不支持
    if iv.pts.is_empty() {
        return None;
    }
    if iv.pts[iv.pts.len() - 1].mult != 0 {
        return None;
    }
    if iv.deg > 1 && iv.reduce() != 0 {
        return None;
    }
    if iv.genus <= 0 {
        return None;
    }
    if iv.deg != 1 {
        return None;
    }
    if iv.pts[iv.pts.len() - 1].psi as i32 >= 3 * iv.genus as i32 - 1 {
        return None;
    }

    let (dim, deg, g) = (iv.dim as i32, iv.deg as i32, iv.genus as i32);
    let dm = dim + 2 - deg;

    // 变量数上界（度范围 × 每度 H 幂个数）
    let mx = (dim + 1) * ((3 * g + 2 * dim - deg) / dm);
    let mut sys = Linsys::new(mx as usize);
    sys.add_var((dim + 1) * iv.cdeg as i32 + iv.pts[iv.pts.len() - 1].h as i32);

    // 最小可用的 k 与对应度 d
    let last = iv.pts[iv.pts.len() - 1];
    let mut k = last.h as i32 + last.psi as i32 - 1;
    let mut d = iv.cdeg as i32;
    d += 1;
    k += dm;
    d -= k / dm;
    k %= dm;
    if k == 0 {
        d += 1;
        k += dm;
    }
    if d < 0 {
        k -= dm * d;
        d = 0;
    }

    // 去掉最后一点后的不变量（L_k 的 t 幂）
    let mut iw = iv.clone();
    let n = iw.pts.len();
    iw.pts = iw.pts[..n - 1].to_vec().into();

    let mut t = LinsysT::Undetermined;
    let mut j = 0;
    for _ in 1..=(mx + VIRASORO_TRY) {
        t = virasoro_k(ctx, &mut sys, &iw, k, d);
        if t == LinsysT::Contradiction {
            eprintln!("*virasoro* got a contradiction for {iv}");
            return None;
        }
        if t == LinsysT::Solvable {
            j += 1;
        }
        if j > VIRASORO_CHECK {
            break;
        }
        d += 1;
        k += dm;
    }

    if t != LinsysT::Solvable {
        return None;
    }

    // 把解出的全部不变量入库
    for i in 0..sys.numvars() {
        let name = sys.varname(i);
        let deg2 = name / (dim + 1);
        let h2 = name % (dim + 1);
        let mut w = iv.clone();
        let ln = w.pts.len();
        w.pts[ln - 1].psi += w.pts[ln - 1].h - h2 as i16 + (dm * (deg2 - w.cdeg as i32)) as i16;
        w.cdeg = deg2 as i16;
        w.pts[ln - 1].h = h2 as i16;
        let val = sys.value_by_number(i)?;
        w.canonicalize();
        ctx.cache.put(w, val);
    }

    sys.value_by_name((dim + 1) * iv.cdeg as i32 + iv.pts[iv.pts.len() - 1].h as i32)
}
