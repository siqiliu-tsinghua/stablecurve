//! 相对分裂入口——原版 `libgw/relsplit.cc` 的直译。
//!
//! 两条路：
//! - 相对亏格 0 不变量（最后一点带重数）：[`rs_genus0_rel`]，
//!   枚举全部分量分裂（relinv）+ 主 H 项 + 主 ψ 项
//! - 绝对亏格 0 不变量（超曲面）：[`rs_genus0_abs`]——Gathmann 论文的
//!   技巧：把最后一点重数设成"超限一档"（cdeg·E+1），化为相对问题再算

use crate::ctx::Ctx;
use crate::inv::{Inv, Pdata};
use crate::num::{Q, Zero};

use crate::relinv::RelInv;

/// 相对亏格 0 不变量。最后一点必须是相对点。`skip0`：跳过开头的修正项数。
fn rs_genus0_rel(ctx: &mut Ctx, iv: &mut Inv, mut skip0: i32) -> Option<Q> {
    let mut result = Q::ZERO;
    let j = iv.pts.len() - 1;
    let ivdeg = iv.deg as i64;
    let ivdim = iv.dim;
    let mut xd = 0i32; // 是否切到环境射影空间（记录原 xdeg）

    // 修正项：重数降一档后的相对分裂
    iv.pts[j].mult -= 1;
    let mut rv = RelInv::new(iv, j as i32);

    loop {
        skip0 -= 1;
        if crate::env_flag("SC_GROWI_REL") && ctx.depth() <= 8 {
            eprintln!(
                "[rel cand depth={} skip0={skip0} mult={}] cof={}",
                ctx.depth(),
                rv.mult_cond(),
                rv.cof()
            );
            for c in 0..=rv.cpts() as i32 {
                eprintln!("    cpt{c} = {}", rv.cpt(c));
            }
        }
        if skip0 < 0 && rv.mult_cond() {
            // 内部分量（绝对不变量）
            let mut iw = rv.cpt(0);
            iw.delmult();
            let mut term = ctx.compute(&iw)?;
            if term != Q::ZERO {
                // 外部分量从后往前（最容易为零，可提前剪枝）
                for i in (1..=rv.cpts() as i32).rev() {
                    let iw = rv.cpt(i);
                    term = &term * &ctx.compute(&iw)?;
                    if term == Q::ZERO {
                        break;
                    }
                }
                if term != Q::ZERO {
                    let contrib = &(rv.cof() * &term);
                    if crate::env_flag("SC_GROWI_REL") && ctx.depth() <= 8 {
                        eprintln!("    CONTRIB -{contrib}");
                    }
                    result -= contrib;
                }
            }
        }
        if !rv.next() {
            break;
        }
    }

    // 重数降到 0 且前一点也无重数 → 换到环境 P^(n+1)（E=1）
    if iv.pts[j].mult == 0 && (j == 0 || iv.pts[j - 1].mult == 0) {
        xd = iv.deg as i32;
        iv.deg = 1;
        iv.dim += 1;
    }

    // 主 H 项
    if iv.pts[j].h >= 0 && (iv.pts[j].h < ivdim || (iv.pts[j].h <= ivdim && iv.pts[j].mult == 0)) {
        iv.pts[j].h += 1;
        let mh = ctx.compute(iv)?;
        if crate::env_flag("SC_GROWI_REL") && ctx.depth() <= 2 {
            eprintln!("[mainH depth={}] {ivdeg} * {mh}", ctx.depth());
        }
        result += Q::from(ivdeg) * mh;
        iv.pts[j].h -= 1;
    }

    // 主 ψ 项
    if iv.pts[j].mult != 0 {
        iv.pts[j].psi += 1;
        let mp = ctx.compute(iv)?;
        if crate::env_flag("SC_GROWI_REL") && ctx.depth() <= 2 {
            eprintln!("[mainPsi depth={}] {} * {mp}", ctx.depth(), iv.pts[j].mult);
        }
        result += Q::from(iv.pts[j].mult as i64) * mp;
        iv.pts[j].psi -= 1;
    }

    // 复原空间
    if xd != 0 {
        iv.deg = xd as i16;
        iv.dim -= 1;
    }
    iv.pts[j].mult += 1;

    Some(result)
}

/// 绝对亏格 0 不变量（超曲面）。可以无标记点。
fn rs_genus0_abs(ctx: &mut Ctx, iv: &mut Inv) -> Option<Q> {
    let mut f = false;

    // 无标记点：用除子公理加一个 H 点
    if iv.pts.is_empty() {
        f = true;
        iv.pts = vec![Pdata::new(1, 0, 0)].into();
    }

    // 末点重数设为超限一档
    let n = iv.pts.len();
    iv.pts[n - 1].mult = (iv.cdeg as i32 * iv.deg as i32 + 1) as i16;

    let mut result = rs_genus0_rel(ctx, iv, 1)?;

    // 复原
    iv.pts[n - 1].mult = 0;
    if f {
        iv.pts = Vec::new().into();
        result /= Q::from(iv.cdeg as i64); // 补偿除子公理
    }

    Some(result)
}

/// 相对分裂入口。原版 `relsplit`。成功入库返回；失败返回 `None`。
pub fn relsplit(ctx: &mut Ctx, iv: &mut Inv) -> Option<Q> {
    // 相对亏格 0
    if iv.genus == 0 && !iv.pts.is_empty() && iv.pts[iv.pts.len() - 1].mult != 0 {
        let r = rs_genus0_rel(ctx, iv, 0)?;
        return store(ctx, iv, r);
    }

    // 绝对亏格 0
    if iv.genus == 0 && (iv.pts.is_empty() || iv.pts[iv.pts.len() - 1].mult == 0) {
        let r = rs_genus0_abs(ctx, iv)?;
        return store(ctx, iv, r);
    }

    None
}

fn store(ctx: &mut Ctx, iv: &Inv, r: Q) -> Option<Q> {
    ctx.cache.put(iv.clone(), r.clone());
    Some(r)
}
