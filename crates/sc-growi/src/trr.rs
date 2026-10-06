//! 任意亏格拓扑递归——原版 `libgw/trr.cc` 的直译（M2 部分）。
//!
//! 两个引擎：
//! - [`negative_genus`]：负亏格（链）不变量的递推——ψ 项减分裂项，
//!   见 `trr.h` 头注的链图说明；
//! - [`trr_reduce`]：亏格 g 不变量在 ψ ≥ 3g−1 的点上展开——类分裂 ×
//!   ψ/节点数 k × 度分裂，对角类在 v 侧**下标 1**（链不变量的第二个
//!   特殊点），w 侧的亏格由 k 决定。
//!
//! `trr_variable`（把未知不变量留作线性方程组变量的版本）供 Virasoro
//! 使用，M3 随 linsys 一起移植。

use crate::ctx::Ctx;
use crate::inv::{Inv, Pdata};
use crate::num::{One, Q, Zero};

use crate::split2::{PairSplit, idiv};

/// 负亏格（链）不变量递推。原版 `negative_genus`。
///
/// 链不变量 <p0, p1, ...>_{-k} = <p0ψ, p1, ...>_{1-k} − Σ 分裂项，
/// 其中分裂把 [p0 | p1, 其余类]，v 侧亏格 0、w 侧亏格 −(k−1)。
fn negative_genus(ctx: &mut Ctx, x: &Inv) -> Option<Q> {
    let mut iv = x.clone();

    // ψ 项：亏格 +1、首点 ψ+1
    iv.genus += 1;
    iv.pts[0].psi += 1;
    let mut result = ctx.compute(&iv)?;
    iv.pts[0].psi -= 1;

    // 分裂项：iv = [p0, p0, 其余...]（亏格 0），iw = [_, p1]（亏格 g+1）
    let iw = Inv::new(
        x.dim,
        x.deg,
        iv.genus,
        x.cdeg,
        [Pdata::default(), iv.pts[1]],
    );
    iv.pts[1] = iv.pts[0];
    iv.genus = 0;
    result -= crate::split2::split2_sum(ctx, &iv, &iw, 2, 2, 0)?;

    Some(result)
}

/// 亏格 g 的拓扑递归。原版 `trr_reduce`。
/// 前提：最后一点的 ψ 幂 ≥ 3g−1（规范形下即任意点足够多的 ψ）。
fn trr_reduce(ctx: &mut Ctx, iv: &Inv) -> Option<Q> {
    let n = iv.pts.len();
    let g = iv.genus as i32;
    let (dim, deg, d) = (iv.dim as i32, iv.deg as i32, iv.cdeg as i32);
    let dm = dim + 2 - deg; // 注意：原版此处不做 dm≤0 分支（E>n+2 时除零）

    // v 侧骨架：[末点(ψ+1−3g), 对角, x0..x_{n-2}]——注意 trr_reduce 用的是
    // **最后**一个（ψ 重的）标记点；trr_variable 才用首点。
    let mut v_pts = vec![Pdata::default(); n + 1];
    v_pts[0] = iv.pts[n - 1];
    v_pts[0].psi += 1 - 3 * iv.genus;
    for i in 0..n - 1 {
        v_pts[i + 2] = iv.pts[i];
    }
    let v_hdr = Inv::new(iv.dim, iv.deg, 0, iv.cdeg, v_pts);
    let w_hdr = Inv::new(iv.dim, iv.deg, iv.genus, 0, [Pdata::default()]);

    let mut result = Q::ZERO;
    let mut coeff = Q::ONE / Q::from(iv.deg as i64); // 对角类因子 1/E
    let mut dw: i32 = (dim - 3) * (1 - g);

    let mut sp = PairSplit::new(&v_hdr, &w_hdr, 2, 1);

    loop {
        // ψ/节点数 k 从 3g−2 降到 0
        for k in (0..=3 * g - 2).rev() {
            let v_genus = k + 2 - 3 * g;
            let d1 = idiv(dim - deg - dw + k, dm).max(0);
            let d2 = idiv(dim - 1 - dw + k, dm).min(d);

            for wd in d1..=d2 {
                // 对角类：v 侧在**下标 1**，w 侧在下标 0
                let hw = wd * dm + dw - k + 1;
                let hv = dim - hw;

                let (av, mut aw) = sp.materialize(
                    &v_hdr,
                    &w_hdr,
                    (d - wd) as i16,
                    wd as i16,
                    hv as i16,
                    hw as i16,
                    1,
                );
                let av = Inv {
                    genus: v_genus as i16,
                    ..av
                };
                aw.pts[0].psi = k as i16;

                if crate::env_flag("SC_GROWI_VIRASORO") {
                    eprintln!("[trr-red] need term {coeff}*{av}*{aw}");
                }
                let tmpv = ctx.compute(&av)?;
                if tmpv == Q::ZERO {
                    continue;
                }
                let tmpw = ctx.compute(&aw)?;
                if tmpw == Q::ZERO {
                    continue;
                }
                result += &(&coeff * &tmpv) * &tmpw;
            }
        }
        if !sp.next(&mut dw, &mut coeff) {
            break;
        }
    }

    Some(result)
}

/// trr_reduce 的"变量版"——原版 `trr_variable`（供 Virasoro 使用）。
///
/// 与 [`trr_reduce`] 的差异：在**第一**个标记点上展开；能被 simplify
/// 直接算出的当常数；ψ 不足的当线性变量（名字 `(xdim+1)·d+Hpow`）；
/// 2 点特例里 w 侧只查 simplify，算不出的也当变量。返回常数部分。
pub(crate) fn trr_variable(
    ctx: &mut Ctx,
    sys: &mut crate::linsys::Linsys,
    iv: &Inv,
    cof: &Q,
) -> Option<Q> {
    // 已知？小心不要排序原不变量（simplify 排序的是拷贝）
    let mut probe = iv.clone();
    if let Some(r) = crate::simplify::simplify(ctx, &mut probe) {
        return Some(cof * r);
    }

    // ψ 不足 → 留作变量
    if (iv.pts[0].psi as i32) < 3 * (iv.genus as i32) - 1 {
        let name = iv.cdeg as i32 * (iv.dim as i32 + 1) + iv.pts[0].h as i32;
        if !sys.add_var_term(name, -cof) {
            eprintln!("*trr* too many variables!");
            return None;
        }
        return Some(Q::ZERO);
    }

    // 分裂结构（与 trr_reduce 相同的骨架）
    let n = iv.pts.len();
    let g = iv.genus as i32;
    let (dim, deg, d) = (iv.dim as i32, iv.deg as i32, iv.cdeg as i32);
    let dm = dim + 2 - deg;

    let mut v_pts = vec![Pdata::default(); n + 1];
    v_pts[0] = iv.pts[0];
    v_pts[0].psi += 1 - 3 * iv.genus;
    for i in (1..n).rev() {
        v_pts[i + 1] = iv.pts[i];
    }
    let v_hdr = Inv::new(iv.dim, iv.deg, 0, iv.cdeg, v_pts);
    let w_hdr = Inv::new(iv.dim, iv.deg, iv.genus, 0, [Pdata::default()]);

    let mut result = Q::ZERO;
    let mut coeff = Q::ONE / Q::from(deg as i64);
    let mut dw: i32 = (dim - 3) * (1 - g);
    let mut sp = PairSplit::new(&v_hdr, &w_hdr, 2, 1);

    loop {
        for k in (0..=3 * g - 2).rev() {
            let v_genus = k + 2 - 3 * g;
            let d1 = idiv(dim - deg - dw + k, dm).max(0);
            let d2 = idiv(dim - 1 - dw + k, dm).min(d);

            for wd in d1..=d2 {
                let hw = wd * dm + dw - k + 1;
                let hv = dim - hw;

                let (av, mut aw) = sp.materialize(
                    &v_hdr,
                    &w_hdr,
                    (d - wd) as i16,
                    wd as i16,
                    hv as i16,
                    hw as i16,
                    1,
                );
                let av = Inv {
                    genus: v_genus as i16,
                    ..av
                };
                aw.pts[0].psi = k as i16;

                if crate::env_flag("SC_GROWI_VIRASORO") {
                    eprintln!("[trr-var] need term {coeff}*{av}*{aw}");
                }
                let tmpv = ctx.compute(&av)?;
                if tmpv == Q::ZERO {
                    continue;
                }

                if av.pts.len() != 2 {
                    let tmpw = ctx.compute(&aw)?;
                    if tmpw == Q::ZERO {
                        continue;
                    }
                    result += &(&coeff * &tmpv) * &tmpw;
                } else {
                    // 2 点特例：w 侧只用 simplify（勿排序原骨架）
                    let mut probe2 = aw.clone();
                    let tmpw = crate::simplify::simplify(ctx, &mut probe2);
                    match tmpw {
                        Some(v) if v != Q::ZERO => {
                            result += &(&coeff * &tmpv) * &v;
                        }
                        Some(_) => {} // 零：跳过
                        None => {
                            // 留作变量（名字取自未排序的 aw）
                            let name = aw.cdeg as i32 * (dim + 1) + aw.pts[0].h as i32;
                            if !sys.add_var_term(name, -&(cof * &coeff * &tmpv)) {
                                eprintln!("*trr* too many variables!");
                                return None;
                            }
                        }
                    }
                }
            }
        }
        if !sp.next(&mut dw, &mut coeff) {
            break;
        }
    }

    Some(cof * result)
}

/// 任意亏格 TRR 入口。原版 `trr_genusg`。
pub fn trr_genusg(ctx: &mut Ctx, iv: &mut Inv) -> Option<Q> {
    // 相对不变量不支持
    if !iv.pts.is_empty() && iv.pts[iv.pts.len() - 1].mult != 0 {
        return None;
    }
    // 原始类暂不支持
    if iv.deg > 1 && iv.reduce() != 0 {
        return None;
    }

    // 负亏格递推
    if iv.genus < 0 {
        if let Some(r) = negative_genus(ctx, iv) {
            return store(ctx, iv, r);
        }
    }

    // 某点 ψ ≥ 3g−1 → 一般拓扑递归
    if iv.genus > 0
        && !iv.pts.is_empty()
        && iv.pts[iv.pts.len() - 1].psi as i32 >= 3 * iv.genus as i32 - 1
    {
        if let Some(r) = trr_reduce(ctx, iv) {
            return store(ctx, iv, r);
        }
    }

    None
}

fn store(ctx: &mut Ctx, iv: &Inv, r: Q) -> Option<Q> {
    ctx.cache.put(iv.clone(), r.clone());
    Some(r)
}
