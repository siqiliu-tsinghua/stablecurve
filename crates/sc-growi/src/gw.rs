//! GW() 入口——原版 `libgw/gwinv.cc` 的直译。
//!
//! 职责：输入校验 + 多线性展开（多项式插入类按单项式笛卡尔积展开）+
//! 逐组合构造不变量并调 `compute`。原版用"里程计"循环增量维护系数与
//! 虚维数；本实现每个组合直接重算（边界层组合数很小，语义等价、清晰）。

use crate::ctx::Ctx;
use crate::inv::{Inv, Pdata};
use crate::num::{One, Q, Zero};
use crate::poly::{Mono, Poly, Var};

/// 一个插入类：多项式 + 相对重数 + 重复次数（`H^2:9` 的 9）。
pub struct GwClass {
    pub poly: Poly,
    pub mult: i16,
    pub count: u32,
}

impl GwClass {
    pub fn new(poly: Poly) -> Self {
        GwClass {
            poly,
            mult: 0,
            count: 1,
        }
    }

    pub fn with_mult(poly: Poly, mult: i16) -> Self {
        GwClass {
            poly,
            mult,
            count: 1,
        }
    }
}

/// GW 查询的设置项（N/E/G/D）。
#[derive(Clone, Debug)]
pub struct GwSettings {
    pub dim: i16,
    pub deg: i16,
    pub genus: i16,
    pub cdeg: i16,
}

impl Default for GwSettings {
    /// 原版缺省：N=0, E=1, G=0, D=0
    fn default() -> Self {
        GwSettings {
            dim: 0,
            deg: 1,
            genus: 0,
            cdeg: 0,
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum GwError {
    Invalid(String),
    BadDimension(i32),
    Unknown(Inv),
}

impl std::fmt::Display for GwError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GwError::Invalid(m) => write!(f, "invalid input: {m}"),
            GwError::BadDimension(vd) => write!(f, "invariant has dimension {vd}"),
            GwError::Unknown(iv) => write!(f, "don't know how to compute {iv}"),
        }
    }
}

impl std::error::Error for GwError {}

/// set_point 的结果。原版用整数：0=零贡献、-1=非法、d>0=系数。
enum SetPoint {
    Zero,
    Invalid,
    Coeff(i32),
}

/// 把一个单项式写成标记点的类。原版 `set_point` 直译。
fn set_point(pt: &mut Pdata, mono: &Mono, rel: bool, xdim: i16, xdeg: i16) -> SetPoint {
    let mut cof = 1i32;

    pt.h = 0;
    pt.psi = 0;
    for &(v, e) in mono.pairs() {
        match v {
            Var::H => {
                if pt.h < 0 {
                    return SetPoint::Zero; // 原始类 × H
                }
                pt.h += e as i16;
            }
            Var::Psi => {
                pt.psi += e as i16;
            }
            Var::Prim(k) => {
                if xdeg == 1 {
                    return SetPoint::Invalid; // P^n 无原始类
                }
                if pt.h < 0 || e > 2 {
                    return SetPoint::Zero;
                }
                if e == 1 {
                    if pt.h != 0 {
                        return SetPoint::Zero;
                    }
                    pt.h = -(k as i16);
                } else {
                    // 同一原始类的平方 = (1/E)·H^top
                    pt.h += xdim;
                    cof = xdeg as i32;
                }
            }
        }
    }

    // H 幂上界
    if rel && pt.mult == 0 {
        if pt.h > xdim + 1 {
            return SetPoint::Zero;
        }
    } else if pt.h > xdim {
        return SetPoint::Zero;
    }

    SetPoint::Coeff(cof)
}

/// GW 主函数。等价原版 `GW(gwlist)`。
pub fn gw(ctx: &mut Ctx, s: &GwSettings, classes: &[GwClass]) -> Result<Q, GwError> {
    // —— 校验（原版逐条对应）——
    if s.dim < 0 {
        return Err(GwError::Invalid(format!("dimension N={} invalid", s.dim)));
    }
    if s.deg < 1 {
        return Err(GwError::Invalid(format!("degree E={} invalid", s.deg)));
    }
    if s.genus < 0 {
        return Err(GwError::Invalid(format!("genus G={} invalid", s.genus)));
    }
    if s.cdeg < 0 {
        return Err(GwError::Invalid(format!(
            "curve degree D={} invalid",
            s.cdeg
        )));
    }

    let mut pts: Vec<Pdata> = Vec::new();
    let mut polys: Vec<&Poly> = Vec::new();
    let mut tm: i64 = 0;
    let mut rel = false;
    for c in classes {
        if c.count == 0 {
            continue;
        }
        // 原版逐点先查零类：任一插入类为零 ⇒ 不变量为 0（多线性），
        // 不能把该点删掉（那会变成另一个不变量）
        if c.poly.is_zero() {
            return Ok(Q::ZERO);
        }
        if c.mult < 0 {
            return Err(GwError::Invalid("negative multiplicity".into()));
        }
        for _ in 0..c.count {
            pts.push(Pdata::new(0, 0, c.mult));
            polys.push(&c.poly);
        }
        tm += c.mult as i64 * c.count as i64;
        if c.mult != 0 {
            rel = true;
        }
    }

    if s.cdeg != 0 && s.dim == 0 && !rel {
        return Err(GwError::Invalid(format!(
            "degree D={} invalid for N=0",
            s.cdeg
        )));
    }
    let max_mult = s.deg as i64 * s.cdeg as i64;
    if tm > max_mult {
        return Err(GwError::Invalid(format!(
            "total multiplicity {tm} exceeds maximum {}",
            max_mult
        )));
    }

    let npts = pts.len();

    // —— 多线性展开：遍历各多项式单项式的全组合 ——
    let mut result = Q::ZERO;
    let mut idx = vec![0usize; npts];

    'combos: loop {
        // 当前组合：各点取 polys[i] 的第 idx[i] 个单项式
        let mut coeff = Q::ONE;
        let mut zero_combo = false;

        for i in 0..npts {
            let (c, m) = &polys[i].terms()[idx[i]];
            coeff *= c;
            match set_point(&mut pts[i], m, rel, s.dim, s.deg) {
                SetPoint::Invalid => {
                    return Err(GwError::Invalid(format!(
                        "polynomial {} contains invalid variables",
                        polys[i]
                    )));
                }
                SetPoint::Zero => {
                    zero_combo = true; // 本组合零贡献，但仍要推进里程计
                }
                SetPoint::Coeff(k) => {
                    if k > 1 {
                        coeff /= Q::from(k as i64);
                    }
                }
            }
        }

        if !zero_combo {
            let iv = Inv::new(s.dim, s.deg, s.genus, s.cdeg, pts.iter().copied());
            // 负亏格：第二点不得有 ψ 类（本入口 genus≥0，防御性保留）
            if iv.genus < 0 && npts >= 2 && iv.pts[1].psi != 0 {
                return Err(GwError::Invalid(format!(
                    "psi class on second point of {iv}"
                )));
            }
            let vd = iv.virtual_dim();
            if vd != 0 {
                return Err(GwError::BadDimension(vd));
            }
            match ctx.compute(&iv) {
                Some(v) => result += &coeff * &v,
                None => return Err(GwError::Unknown(iv)),
            }
        }

        // 里程计推进（零类时只有空组合，走一次即止）
        if npts == 0 {
            break 'combos;
        }
        for i in 0..npts {
            idx[i] += 1;
            if idx[i] < polys[i].terms().len() {
                break;
            }
            idx[i] = 0;
            if i == npts - 1 {
                break 'combos;
            }
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::q;

    /// lines in P^2 thru H^2:2 → 1
    #[test]
    fn lines_p2() {
        let mut ctx = Ctx::memory();
        let s = GwSettings {
            dim: 2,
            deg: 1,
            genus: 0,
            cdeg: 1,
        };
        let cls = vec![GwClass {
            poly: Poly::from(Var::H).pow(2),
            mult: 0,
            count: 2,
        }];
        assert_eq!(gw(&mut ctx, &s, &cls), Ok(q(1)));
    }

    /// 零类 ⇒ 不变量为 0（原版 gwinv.cc 语义）；曾被当成"删掉这个点"。
    #[test]
    fn zero_class_gives_zero() {
        let mut ctx = Ctx::memory();
        for e in [
            "lines in P^2 thru H^2:2, H-H",
            "lines in P^2 thru H^2:2, 0",
            "lines in P^2 thru 0",
        ] {
            let qy = crate::query::parse(e).unwrap();
            assert_eq!(gw(&mut ctx, &qy.settings, &qy.classes), Ok(Q::ZERO), "{e}");
        }
    }
}
