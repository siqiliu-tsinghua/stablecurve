//! 二分量分裂——原版 `libgw/split2.cc` 的直译。
//!
//! 核心 [`PairSplit`] 是原版 `next_pair`（4 指针双缓冲交换）的安全重写：
//! 枚举器自持双缓冲与两侧完整点集，`next()` 返回后调用方用 `v()/w()`
//! 读当前分裂并直接在 w 侧填写次数与对角类。`dw`（w 侧虚维数）与
//! `coeff`（组合系数）由调用方持有，`next()` 做相对更新（原版语义）。
//! 枚举顺序：初始状态（全部可动类在 v 侧）是第一个分裂；耗尽时复原。
//!
//! 其上是亏格 0 的 WDVV 重建阶梯：两点单/双 descendant → Kontsevich–Manin
//! 首重建定理 → g=0 TRR。g=1 TRR（环路项）在 M2 随 trr.rs 加入。

use crate::ctx::Ctx;
use crate::inv::{Inv, Pdata};
use crate::num::{One, Q, Zero};

/// 真正向下取整的整除（原版 `idiv`；m 可为负，n 为正）。
pub(crate) fn idiv(m: i32, n: i32) -> i32 {
    if m >= 0 { m / n } else { -((n - m - 1) / n) }
}

/// 类在两个分量间的分裂枚举器（原版 `next_pair` 的双缓冲直译）。
///
/// 缓冲按 v 初始 + w 初始的总长分配（原版靠调用方给 inv 超额分配），
/// 当前长度由 `vn`/`wn` 显式跟踪——w 侧在枚举中会增长。
pub(crate) struct PairSplit {
    v: Vec<Pdata>,
    w: Vec<Pdata>,
    v2: Vec<Pdata>,
    w2: Vec<Pdata>,
    vn: usize,
    wn: usize,
    fv: usize,
    fw: usize,
    xdim: i16,
}

impl PairSplit {
    pub(crate) fn new(v: &Inv, w: &Inv, fv: usize, fw: usize) -> Self {
        let total = v.pts.len() + w.pts.len();
        let pad = |src: &[Pdata]| -> Vec<Pdata> {
            let mut b = src.to_vec();
            b.resize(total, Pdata::default());
            b
        };
        PairSplit {
            v: pad(&v.pts),
            w: pad(&w.pts),
            v2: pad(&v.pts),
            w2: pad(&w.pts),
            vn: v.pts.len(),
            wn: w.pts.len(),
            fv,
            fw,
            xdim: v.dim,
        }
    }

    /// 走到下一个分裂。返回 false 表示耗尽（缓冲已复原）。
    pub(crate) fn next(&mut self, dw: &mut i32, coeff: &mut Q) -> bool {
        if self.vn == self.fv {
            return false;
        }

        // 1) 把 w 侧开头小于 v 首类的整块类挪回 v 侧（成块移动不产生二项式
        //    系数，但改变虚维数）
        let mut i = 0usize;
        while self.fw + i < self.wn && self.w[self.fw + i] < self.v[self.fv] {
            *dw += self.w[self.fw + i].vdim(self.xdim) / 2;
            i += 1;
        }

        // 2) 从 v 侧挪一个类到 w 侧
        let moving = self.v[self.fv];
        *dw -= moving.vdim(self.xdim) / 2;

        // 3) 二项式系数：×（v 侧同类总数）、÷（w 侧既有同类数 +1）
        let mut n = 1usize;
        while self.fv + n < self.vn && self.v[self.fv + n] == moving {
            n += 1;
        }
        *coeff *= Q::from(n as i64);
        let mut m = 1usize;
        while m < self.wn - self.fw - i + 1 && self.w[self.fw + i + m - 1] == moving {
            m += 1;
        }
        *coeff /= Q::from(m as i64);

        // 4) 组装新分裂到双缓冲。固定前缀（含对角槽位 0）保持不变。
        //    注意：moving 落在 w 可动区前端（位置 fw）——原版 iw_dst 不随
        //    回挪的 i 前进，只有 iw_src 前进。
        let w_tail = self.wn - self.fw - i;
        let v_tail = self.vn - self.fv - 1;

        let (v2, w2) = (&mut self.v2, &mut self.w2);
        w2[self.fw] = moving;
        w2[self.fw + 1..self.fw + 1 + w_tail]
            .copy_from_slice(&self.w[self.fw + i..self.fw + i + w_tail]);
        v2[self.fv..self.fv + i].copy_from_slice(&self.w[self.fw..self.fw + i]);
        v2[self.fv + i..self.fv + i + v_tail]
            .copy_from_slice(&self.v[self.fv + 1..self.fv + 1 + v_tail]);

        let new_v = self.vn + i - 1;
        let new_w = self.wn - i + 1;

        // 5) 交换缓冲
        std::mem::swap(&mut self.v, &mut self.v2);
        std::mem::swap(&mut self.w, &mut self.w2);
        self.vn = new_v;
        self.wn = new_w;

        true
    }

    /// 把当前分裂物化为一对 Inv（填对角类与两侧次数）。
    /// `hdr_v`/`hdr_w` 提供头部（维数、亏格等）与固定前缀模板。
    /// 把当前分裂物化为一对 Inv。`vdiag` 是 v 侧对角类槽位
    /// （split2_sum 约定 0；trr 路径 fv=2、对角在 1）。w 侧恒为 0。
    pub(crate) fn materialize(
        &self,
        hdr_v: &Inv,
        hdr_w: &Inv,
        dv: i16,
        dw_deg: i16,
        hv: i16,
        hw: i16,
        vdiag: usize,
    ) -> (Inv, Inv) {
        let mut a = hdr_v.clone();
        a.cdeg = dv;
        a.pts = self.v[..self.vn].to_vec().into();
        a.pts[vdiag] = Pdata::new(hv, 0, 0);
        let mut b = hdr_w.clone();
        b.cdeg = dw_deg;
        b.pts = self.w[..self.wn].to_vec().into();
        b.pts[0] = Pdata::new(hw, 0, 0);
        (a, b)
    }
}

/// 对所有（类分裂 × 度分裂）求和，对应原版 `split2_sum`。
///
/// `iv`/`iw` 给出两侧骨架（含固定前缀类、亏格）；`fv`/`fw` 为固定前缀长
/// 度（槽位 0 恒为对角类）；`skip` 跳过前若干个度候选项（跨类分裂计数）。
pub(crate) fn split2_sum(
    ctx: &mut Ctx,
    iv: &Inv,
    iw: &Inv,
    fv: usize,
    fw: usize,
    mut skip: i32,
) -> Option<Q> {
    let mut result = Q::ZERO;
    let mut coeff = Q::ONE / Q::from(iv.deg as i64); // 对角类因子 1/E
    let d = iv.cdeg;

    // w 侧头部（原版：xdim/xdeg 由 iv 补齐，cdeg 由度循环填，cgenus 用 iw 的）
    let hdr_w = Inv {
        dim: iv.dim,
        deg: iv.deg,
        genus: iw.genus,
        cdeg: 0,
        pts: iw.pts.to_vec().into(),
    };
    let hdr_v = Inv {
        genus: iv.genus,
        pts: iv.pts.to_vec().into(),
        ..iv.clone()
    };

    // w 侧虚维数初值
    let mut dw: i32 = if iw.genus >= 0 {
        (iv.dim as i32 - 3) * (1 - iw.genus as i32)
    } else {
        iv.dim as i32 - 3 + iw.genus as i32
    };
    for i in 1..fw {
        dw -= iw.pts[i].vdim(iv.dim) / 2;
    }

    // 是否有 descendant 或负亏格（读目标不变量的最后一点 + 两侧固定前缀）
    let mut desc = iv.pts[iv.pts.len() - 1].psi != 0 || iv.genus < 0 || iw.genus < 0;
    if !desc && iv.pts[..fv].iter().any(|p| p.psi != 0) {
        desc = true;
    }
    if !desc && iw.pts[..fw].iter().any(|p| p.psi != 0) {
        desc = true;
    }

    let mut sp = PairSplit::new(iv, iw, fv, fw);

    loop {
        // 度可行界（三分支直译）
        let (d1, d2) = {
            let dm_full = iv.dim as i32 + 2 - iv.deg as i32;
            if dm_full > 0 {
                (
                    idiv(iv.dim as i32 - iv.deg as i32 - dw, dm_full).max(0),
                    idiv(iv.dim as i32 - 1 - dw, dm_full).min(d as i32),
                )
            } else if dm_full < 0 {
                let dm = iv.deg as i32 - iv.dim as i32 - 2;
                (
                    idiv(dw - 2 * iv.dim as i32 + iv.deg as i32 - 2, dm).max(0),
                    idiv(1 + dw, dm).min(d as i32),
                )
            } else {
                (0, d as i32)
            }
        };

        for wd in d1..=d2 {
            skip -= 1;
            if skip >= 0 {
                continue;
            }

            // 对角类：H^dim 分裂为两侧
            let hw = wd * (iv.dim as i32 + 2 - iv.deg as i32) + dw + 1;
            let hv = iv.dim as i32 - hw;

            // 正次数分量上出现基本类且无 descendant → 0
            if hv == 0 && (d - wd as i16) != 0 && !desc {
                continue;
            }
            if hw == 0 && wd != 0 && !desc {
                continue;
            }

            let (av, aw) = sp.materialize(
                &hdr_v,
                &hdr_w,
                d - wd as i16,
                wd as i16,
                hv as i16,
                hw as i16,
                0,
            );

            if crate::env_flag("SC_GROWI_TERMS") {
                eprintln!("[term] coeff={coeff} dw={dw} wd={wd} hv={hv} hw={hw} v={av} w={aw}");
            }

            // 先算小侧（亏格低优先，其次度低），为零即跳过另一侧
            let v_first =
                iv.genus == 0 || iv.genus < iw.genus || (iv.genus == iw.genus && av.cdeg < aw.cdeg);
            let (tmpv, tmpw) = if v_first {
                let tv = ctx.compute(&av)?;
                if tv == Q::ZERO {
                    continue;
                }
                let tw = ctx.compute(&aw)?;
                if tw == Q::ZERO {
                    continue;
                }
                (tv, tw)
            } else {
                let tw = ctx.compute(&aw)?;
                if tw == Q::ZERO {
                    continue;
                }
                let tv = ctx.compute(&av)?;
                if tv == Q::ZERO {
                    continue;
                }
                (tv, tw)
            };
            result += &(&coeff * &tmpv) * &tmpw;
        }

        if !sp.next(&mut dw, &mut coeff) {
            break;
        }
    }

    Some(result)
}

/// 两点、至多一个 descendant 的不变量（WDVV）。原版 `twopoint_onedesc`。
fn twopoint_onedesc(ctx: &mut Ctx, x: &Inv) -> Option<Q> {
    if x.pts[0].h == 0 {
        return None; // a=0：本应先用基本类公理
    }

    // 第二个（近乎平凡的）分裂项
    let iv = Inv::new(
        x.dim,
        x.deg,
        0,
        x.cdeg,
        [
            Pdata::new(1, 0, 0),
            Pdata::new(x.pts[1].h, x.pts[1].psi + 1, 0),
            Pdata::new(x.pts[0].h - 1, 0, 0),
        ],
    );
    let mut result = ctx.compute(&iv)?;

    // 第一个分裂项
    let iv2 = Inv::new(
        x.dim,
        x.deg,
        0,
        x.cdeg,
        [
            Pdata::new(1, 0, 0),
            Pdata::new(x.pts[1].h, x.pts[1].psi + 1, 0),
            Pdata::new(0, 0, 0),
        ],
    );
    let iw2 = Inv::new(
        x.dim,
        x.deg,
        0,
        0,
        [
            Pdata::default(),
            Pdata::new(1, 0, 0),
            Pdata::new(x.pts[0].h - 1, 0, 0),
        ],
    );
    result -= split2_sum(ctx, &iv2, &iw2, 3, 3, 1)?;

    Some(result)
}

/// 两点、双 descendant。原版 `twopoint_twodesc`。
fn twopoint_twodesc(ctx: &mut Ctx, x: &Inv) -> Option<Q> {
    let mut result = Q::ZERO;
    let mut iv = x.clone();

    // 两个除子公理修正项
    if iv.pts[0].h < iv.dim {
        iv.pts[0].h += 1;
        iv.pts[0].psi -= 1;
        let t = ctx.compute(&iv)?;
        iv.pts[0].h -= 1;
        iv.pts[0].psi += 1;
        result = -t;
    }
    if iv.pts[1].h < iv.dim {
        iv.pts[1].h += 1;
        iv.pts[1].psi -= 1;
        let t = ctx.compute(&iv)?;
        iv.pts[1].h -= 1;
        iv.pts[1].psi += 1;
        result -= t;
    }

    // 主分裂项：iv = [x0, x1, 1]，iw = [_, x0, 1]
    let mut iv_pts = iv.pts.to_vec();
    iv_pts.push(Pdata::new(0, 0, 0));
    let iv2 = Inv {
        pts: iv_pts.into(),
        ..x.clone()
    };
    let iw2 = Inv::new(
        x.dim,
        x.deg,
        0,
        0,
        [Pdata::default(), x.pts[0], Pdata::new(1, 0, 0)],
    );
    result += split2_sum(ctx, &iv2, &iw2, 3, 3, 0)?;

    Some(result / Q::from(x.cdeg as i64 * x.deg as i64))
}

/// Kontsevich–Manin 首重建定理。原版 `first_reconstruction`。
fn first_reconstruction(ctx: &mut Ctx, x: &Inv) -> Option<Q> {
    if x.pts[0].h == 0 {
        return None;
    }

    let n = x.pts.len();
    let mut iv_pts = vec![Pdata::default(); n];
    for i in 1..n - 2 {
        iv_pts[i + 2] = x.pts[i];
    }
    iv_pts[1] = x.pts[n - 2];
    iv_pts[2] = x.pts[n - 1];

    // WDVV 左边：iv = [_, x_{n-2}, x_{n-1}, x_1.., 0]，iw = [_, 1, x_0-1]
    let iv_l = Inv::new(x.dim, x.deg, 0, x.cdeg, iv_pts.clone());
    let iw_l = Inv::new(
        x.dim,
        x.deg,
        0,
        0,
        [
            Pdata::default(),
            Pdata::new(1, 0, 0),
            Pdata::new(x.pts[0].h - 1, 0, 0),
        ],
    );
    let mut result = -split2_sum(ctx, &iv_l, &iw_l, 3, 3, 1)?;

    // WDVV 右边：iv = [_, 1, x_{n-1}, x_1.., 0]，iw = [_, x_{n-2}, x_0-1]
    let mut iv_pts_r = iv_pts;
    iv_pts_r[1] = Pdata::new(1, 0, 0);
    let iv_r = Inv::new(x.dim, x.deg, 0, x.cdeg, iv_pts_r);
    let iw_r = Inv::new(
        x.dim,
        x.deg,
        0,
        0,
        [
            Pdata::default(),
            x.pts[n - 2],
            Pdata::new(x.pts[0].h - 1, 0, 0),
        ],
    );
    result += split2_sum(ctx, &iv_r, &iw_r, 3, 3, 0)?;

    Some(result)
}

/// 亏格 0 拓扑递归关系。原版 `trr_genus0`。
/// 前提：≥3 个标记点，最后一点是 descendant。
fn trr_genus0(ctx: &mut Ctx, x: &Inv) -> Option<Q> {
    let n = x.pts.len();
    let mut iv_pts = vec![Pdata::default(); n];
    for i in 0..n - 1 {
        iv_pts[i + 1] = x.pts[i];
    }
    let iv = Inv::new(x.dim, x.deg, 0, x.cdeg, iv_pts);
    let iw = Inv::new(
        x.dim,
        x.deg,
        0,
        0,
        [
            Pdata::default(),
            Pdata::new(x.pts[n - 1].h, x.pts[n - 1].psi - 1, 0),
        ],
    );
    split2_sum(ctx, &iv, &iw, 3, 2, 0)
}

/// 二分量分裂入口。原版 `split2`。
pub fn split2(ctx: &mut Ctx, iv: &mut Inv) -> Option<Q> {
    // 相对不变量、负亏格不支持
    if !iv.pts.is_empty() && iv.pts[iv.pts.len() - 1].mult != 0 {
        return None;
    }
    if iv.genus < 0 {
        return None;
    }
    // 归约后仍有原始类 → 放弃（原版 reduce() 非 0 即有原始类）
    if iv.deg > 1 && iv.reduce() != 0 {
        return None;
    }

    let npts = iv.pts.len();
    let last_psi = if npts > 0 { iv.pts[npts - 1].psi } else { 0 };
    let first_psi = if npts > 1 { iv.pts[0].psi } else { 0 };

    // 两点、至多一个 descendant
    if iv.genus == 0 && npts == 2 && first_psi == 0 {
        if let Some(r) = twopoint_onedesc(ctx, iv) {
            return store(ctx, iv, r);
        }
    }

    // 两点、双 descendant
    if iv.genus == 0 && npts == 2 && first_psi != 0 {
        if let Some(r) = twopoint_twodesc(ctx, iv) {
            return store(ctx, iv, r);
        }
    }

    // n 点、无 descendant
    if iv.genus == 0 && npts > 2 && last_psi == 0 {
        if let Some(r) = first_reconstruction(ctx, iv) {
            return store(ctx, iv, r);
        }
    }

    // n 点、有 descendant：g=0 TRR
    if iv.genus == 0 && npts > 2 && last_psi != 0 {
        if let Some(r) = trr_genus0(ctx, iv) {
            return store(ctx, iv, r);
        }
    }

    // g=1 TRR（Getzler）。E>1 时环路项有原始类 complications，原版即不支持。
    // 原版从最后一点向前遍历 descendant 点，逐个尝试。
    if iv.genus == 1 && iv.deg == 1 {
        for i in (0..npts).rev() {
            if iv.pts[i].psi == 0 {
                break;
            }
            if let Some(r) = trr_genus1(ctx, iv, i) {
                return store(ctx, iv, r);
            }
        }
    }

    None
}

/// 亏格 1 拓扑递归（Getzler）。原版 `trr_genus1`。
/// `pos`：descendant 点的下标。第一项（分裂）+ 环路项（/24，亏格 0）。
fn trr_genus1(ctx: &mut Ctx, x: &Inv, pos: usize) -> Option<Q> {
    let n = x.pts.len();

    // 第一项：iv = [diag, x[pos]ψ−1, x[0..pos), x[pos+1..)]（亏格 0），
    // iw = [diag]（亏格 1）
    let mut iv_pts = vec![Pdata::default(); n + 1];
    for i in 0..pos {
        iv_pts[i + 2] = x.pts[i];
    }
    for i in pos + 1..n {
        iv_pts[i + 1] = x.pts[i];
    }
    iv_pts[1] = x.pts[pos];
    iv_pts[1].psi -= 1;
    let iv = Inv::new(x.dim, x.deg, 0, x.cdeg, iv_pts);
    let iw = Inv::new(x.dim, x.deg, 1, 0, [Pdata::default()]);

    let mut result = split2_sum(ctx, &iv, &iw, 2, 1, 0)?;

    // 环路项：亏格 0，追加对角点对 H^i × H^(dim−i)
    let mut loop_pts = x.pts.to_vec();
    loop_pts[pos].psi -= 1;
    loop_pts.push(Pdata::default());
    loop_pts.push(Pdata::default());
    let iv2 = Inv::new(x.dim, x.deg, 0, x.cdeg, loop_pts);
    for i in 0..=x.dim as usize {
        let mut iv2 = iv2.clone();
        iv2.pts[n].h = i as i16;
        iv2.pts[n + 1].h = x.dim - i as i16;
        result += ctx.compute(&iv2)? / Q::from(24);
    }

    Some(result)
}

fn store(ctx: &mut Ctx, iv: &Inv, r: Q) -> Option<Q> {
    ctx.cache.put(iv.clone(), r.clone());
    Some(r)
}

/// 阶乘（relinv 的组合因子用）。
pub(crate) fn factorial(n: i64) -> Q {
    let mut r = Q::ONE;
    for i in 2..=n {
        r *= Q::from(i);
    }
    r
}
