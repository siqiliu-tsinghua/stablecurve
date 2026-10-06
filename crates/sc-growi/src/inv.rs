//! 不变量数据与规范形——原版 `libgw/inv.{h,cc}` 的移植。
//!
//! 一个不变量 = (流形维数 `dim`、流形次数 `deg`、曲线亏格 `genus`、曲线次数
//! `cdeg`、标记点列表)。每个标记点 [`Pdata`] 记 `(H 幂, psi 幂, 相对重数)`；
//! `h < 0` 表示第 `|h|` 号原始上同调类（只出现在 deg>1 的超曲面里，同一负数
//! 出现多次即同一个原始类）。
//!
//! 规范形：点按 (mult, psi, h) 字典序升序——**mult 是最高位**（原版
//! `inv.cc:16` 注释强调这一点），相对点因此排最后、descendant 排倒数第二
//! 区域，算法各处"最后一个点"的判断依赖此性质。负亏格（链不变量）时前两个
//! 点固定不参与排序（原版 `inv::sort`）。
//!
//! 虚维数（原版 `inv.h vdim` 与 `gwinv.cc`）：注意相对不变量实为环境
//! 射影空间 P^(n+1) 的不变量，维数公式因此分三种情形。

use std::fmt;

/// 标记点数据。字段顺序即规范排序优先级：mult > psi > h。
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Default)]
pub struct Pdata {
    pub mult: i16,
    pub psi: i16,
    pub h: i16,
}

impl Pdata {
    pub fn new(h: i16, psi: i16, mult: i16) -> Self {
        Pdata { mult, psi, h }
    }

    /// 该点对虚维数的贡献（取"实"维数以容纳非代数的原始类；`n` 为流形维数）。
    /// 原版 `pdata::vdim`。
    pub fn vdim(&self, n: i16) -> i32 {
        if self.h < 0 {
            n as i32 + 2 * (self.psi as i32 + self.mult as i32 - 1)
        } else {
            2 * (self.h as i32 + self.psi as i32 + self.mult as i32 - 1)
        }
    }

    /// 同上但忽略重数。原版 `pdata::vdim_nomult`。
    pub fn vdim_nomult(&self, n: i16) -> i32 {
        if self.h < 0 {
            n as i32 + 2 * (self.psi as i32 - 1)
        } else {
            2 * (self.h as i32 + self.psi as i32 - 1)
        }
    }
}

impl fmt::Display for Pdata {
    /// 原版 `operator<<(pdata)` 的直译（psi 显示为 'c'）。
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.mult != 0 {
            write!(f, "[")?;
        }
        if self.h > 0 {
            write!(f, "H")?;
            if self.h > 1 {
                write!(f, "{}", self.h)?;
            }
        }
        if self.h < 0 {
            write!(f, "pc({})", -self.h)?;
        }
        if self.psi != 0 {
            write!(f, "c")?;
            if self.psi > 1 {
                write!(f, "{}", self.psi)?;
            }
        }
        if self.h == 0 && self.psi == 0 {
            write!(f, "1")?;
        }
        if self.mult != 0 {
            write!(f, ",{}]", self.mult)?;
        }
        Ok(())
    }
}

/// 一个（绝对或相对的）Gromov-Witten 不变量。
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Inv {
    pub dim: i16,   // 流形维数（原版 xdim）
    pub deg: i16,   // 流形次数（原版 xdeg）
    pub genus: i16, // 曲线亏格；负数 = 链不变量 g=-k（原版 cgenus）
    pub cdeg: i16,  // 曲线次数（原版 cdeg）
    pub pts: Box<[Pdata]>,
}

impl Inv {
    pub fn new<I: IntoIterator<Item = Pdata>>(
        dim: i16,
        deg: i16,
        genus: i16,
        cdeg: i16,
        pts: I,
    ) -> Self {
        Inv {
            dim,
            deg,
            genus,
            cdeg,
            pts: pts.into_iter().collect(),
        }
    }

    /// 是否相对不变量（存在重数为正的点）。
    pub fn is_relative(&self) -> bool {
        self.pts.iter().any(|p| p.mult > 0)
    }

    /// 排成规范形。原版 `inv::sort`：负亏格时前两个点固定。
    pub fn canonicalize(&mut self) {
        if self.genus >= 0 {
            self.pts.sort();
        } else {
            let (head, tail) = self.pts.split_at_mut(2.min(self.pts.len()));
            let _ = head;
            tail.sort();
        }
    }

    pub fn into_canonical(mut self) -> Self {
        self.canonicalize();
        self
    }

    /// 把原始类重命名为最小可用编号（-1, -2, ...）。返回出现的最小原始类
    /// （无原始类则为 0）。原版 `inv::reduce` 的直译。
    pub fn reduce(&mut self) -> i16 {
        let mut j: i16 = 0;
        for i in 0..self.pts.len() {
            let k = self.pts[i].h;
            if k < j {
                j -= 1;
                if k != j {
                    for p in &mut self.pts[i..] {
                        if p.h == k {
                            p.h = j;
                        } else if p.h == j {
                            p.h = k;
                        }
                    }
                }
            }
        }
        j
    }

    /// 删除第 `i` 个标记点。原版 `inv::del`。
    pub fn del(&mut self, i: usize) {
        self.pts.copy_within(i + 1.., i);
        let n = self.pts.len() - 1;
        self.pts = self.pts[..n].to_vec().into_boxed_slice();
    }

    /// 清除所有点的相对重数。原版 `inv::delmult`。
    pub fn delmult(&mut self) {
        for p in self.pts.iter_mut() {
            p.mult = 0;
        }
    }

    /// 虚维数。原版 `gwinv.cc`（与 `inv.cc` check）的公式：点贡献取负，
    /// 再按相对/负亏格/普通三种情形加基项。规范形下结果应为 0。
    pub fn virtual_dim(&self) -> i32 {
        let mut vd: i32 = -self.pts.iter().map(|p| p.vdim(self.dim)).sum::<i32>();
        let (dim, deg, genus, cdeg) = (
            self.dim as i32,
            self.deg as i32,
            self.genus as i32,
            self.cdeg as i32,
        );
        if self.is_relative() {
            vd += 2 * ((2 + dim) * cdeg + (dim - 2) * (1 - genus));
        } else if genus < 0 {
            vd += 2 * ((2 + dim - deg) * cdeg + dim - 3 + genus);
        } else {
            vd += 2 * ((2 + dim - deg) * cdeg + (dim - 3) * (1 - genus));
        }
        vd
    }

    /// 合法性检查（原版 INV_DEBUG 下 `inv::check` 的直译，收集全部违规）。
    pub fn check(&self) -> Result<(), Vec<String>> {
        let mut errs = Vec::new();
        let bad = |errs: &mut Vec<String>, m: String| errs.push(m);

        if self.dim < 0 {
            bad(&mut errs, format!("dimension N={} invalid", self.dim));
        }
        if self.deg < 1 {
            bad(&mut errs, format!("degree E={} invalid", self.deg));
        }
        if self.cdeg < 0 || (self.cdeg != 0 && self.dim == 0) {
            bad(&mut errs, format!("curve degree D={} invalid", self.cdeg));
        }
        let rel = self.is_relative();
        let mut tm: i32 = 0;
        for (i, p) in self.pts.iter().enumerate() {
            if p.psi < 0 {
                bad(&mut errs, format!("point {i}: psi power {} invalid", p.psi));
            }
            if p.mult < 0 {
                bad(
                    &mut errs,
                    format!("point {i}: multiplicity {} invalid", p.mult),
                );
            }
            if p.h < 0 && self.deg == 1 {
                bad(&mut errs, format!("point {i}: primitive class in P^n"));
            }
            if p.h > self.dim + 1 || (p.h > self.dim && (!rel || p.mult != 0)) {
                bad(&mut errs, format!("point {i}: H power {} invalid", p.h));
            }
            tm += p.mult as i32;
        }
        if rel && self.genus < 0 {
            bad(&mut errs, "relative negative-genus invariant".into());
        }
        let vd = self.virtual_dim();
        if vd != 0 {
            bad(&mut errs, format!("virtual dimension is {vd}"));
        }
        if self.genus < 0 {
            if self.pts.len() < 2 {
                bad(
                    &mut errs,
                    "negative-genus invariant needs >= 2 points".into(),
                );
            } else if self.pts[1].psi != 0 {
                bad(&mut errs, "psi class on second point".into());
            }
        }
        if tm > self.deg as i32 * self.cdeg as i32 {
            bad(
                &mut errs,
                format!("total multiplicity {tm} exceeds {}", self.deg * self.cdeg),
            );
        }
        if errs.is_empty() { Ok(()) } else { Err(errs) }
    }
}

impl fmt::Display for Inv {
    /// 原版 `operator<<(inv)` 的直译（含游程压缩显示 `H:2`）。
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "<n={},e={},g={},d={}",
            self.dim, self.deg, self.genus, self.cdeg
        )?;
        let mut i = 0;
        while i < self.pts.len() {
            write!(f, ",{}", self.pts[i])?;
            let mut n = 1;
            while i + n < self.pts.len() && self.pts[i + n] == self.pts[i] {
                n += 1;
            }
            if n > 1 {
                write!(f, ":{n}")?;
            }
            i += n;
        }
        write!(f, ">")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// P^2 里过两点的直线：<n=2,e=1,g=0,d=1,H:2>
    fn lines_p2() -> Inv {
        Inv::new(2, 1, 0, 1, [Pdata::new(2, 0, 0); 2])
    }

    #[test]
    fn display_matches_original_format() {
        assert_eq!(lines_p2().to_string(), "<n=2,e=1,g=0,d=1,H2:2>");
        let iv = Inv::new(3, 1, 1, 1, [Pdata::new(0, 0, 0), Pdata::new(1, 2, 0)]);
        assert_eq!(iv.to_string(), "<n=3,e=1,g=1,d=1,1,Hc2>");
        let iv = Inv::new(2, 3, 0, 4, [Pdata::new(-1, 0, 0), Pdata::new(0, 0, 2)]);
        assert_eq!(iv.to_string(), "<n=2,e=3,g=0,d=4,pc(1),[1,2]>");
    }

    #[test]
    fn canonical_order_mult_first() {
        // mult 为最高位：相对点排最后；同为 mult=0 时 psi 升序、再 h 升序
        let mut iv = Inv::new(
            3,
            1,
            0,
            1,
            [
                Pdata::new(0, 0, 3),
                Pdata::new(0, 1, 0),
                Pdata::new(2, 0, 0),
                Pdata::new(1, 0, 0),
            ],
        );
        iv.canonicalize();
        let hs: Vec<i16> = iv.pts.iter().map(|p| p.h).collect();
        let psis: Vec<i16> = iv.pts.iter().map(|p| p.psi).collect();
        let mults: Vec<i16> = iv.pts.iter().map(|p| p.mult).collect();
        assert_eq!(mults, vec![0, 0, 0, 3]);
        assert_eq!(psis, vec![0, 0, 1, 0]);
        assert_eq!(hs, vec![1, 2, 0, 0]);
    }

    #[test]
    fn negative_genus_keeps_first_two() {
        let fixed0 = Pdata::new(0, 2, 0);
        let fixed1 = Pdata::new(0, 0, 0);
        let mut iv = Inv::new(
            3,
            1,
            -2,
            2,
            [fixed0, fixed1, Pdata::new(3, 0, 0), Pdata::new(1, 0, 0)],
        );
        iv.canonicalize();
        assert_eq!(iv.pts[0], fixed0);
        assert_eq!(iv.pts[1], fixed1);
        assert!(iv.pts[2] <= iv.pts[3]);
    }

    #[test]
    fn reduce_renames_primitives() {
        // 单个 -3 → -1
        let mut iv = Inv::new(3, 2, 0, 1, [Pdata::new(-3, 0, 0), Pdata::new(2, 0, 0)]);
        assert_eq!(iv.reduce(), -1);
        assert_eq!(iv.pts[0].h, -1);
        // 乱序的两个原始类 -5,-2 → -1,-2
        let mut iv = Inv::new(
            3,
            2,
            0,
            1,
            [
                Pdata::new(-5, 0, 0),
                Pdata::new(-2, 0, 0),
                Pdata::new(0, 1, 0),
            ],
        );
        assert_eq!(iv.reduce(), -2);
        let mut hs: Vec<i16> = iv.pts.iter().map(|p| p.h).collect();
        hs.sort();
        assert_eq!(hs, vec![-2, -1, 0]);
        // 无原始类
        let mut iv = lines_p2();
        assert_eq!(iv.reduce(), 0);
    }

    #[test]
    fn virtual_dim_known_cases() {
        // P^2 过两点的直线：0
        assert_eq!(lines_p2().virtual_dim(), 0);
        // P^3 过两点的直线（两个 H^3 点条件）：0
        let iv = Inv::new(3, 1, 0, 1, [Pdata::new(3, 0, 0); 2]);
        assert_eq!(iv.virtual_dim(), 0);
        // 椭圆亏格、P^3、过两点的直线（gwcheck: -1/12 那个不变量）：0
        assert_eq!(
            Inv::new(3, 1, 1, 1, [Pdata::new(3, 0, 0); 2]).virtual_dim(),
            0
        );
        // 单个 <psi> 点、亏格 1：vd = 0 + 2*4 = 8 ≠ 0 → 非法
        assert_eq!(Inv::new(3, 1, 1, 1, [Pdata::new(0, 1, 0)]).virtual_dim(), 8);
    }

    #[test]
    fn check_valid_and_invalid() {
        assert!(lines_p2().check().is_ok());
        // 维数不为零
        let bad = Inv::new(2, 1, 0, 1, [Pdata::new(2, 0, 0)]);
        assert!(bad.check().is_err());
        // 原始类出现在 P^n（E=1）
        let bad = Inv::new(
            3,
            1,
            0,
            1,
            [
                Pdata::new(-1, 1, 0),
                Pdata::new(3, 0, 0),
                Pdata::new(3, 0, 0),
            ],
        );
        assert!(bad.check().is_err());
    }

    #[test]
    fn delete_point() {
        let mut iv = lines_p2();
        iv.del(0);
        assert_eq!(iv.pts.len(), 1);
        assert_eq!(iv.pts[0].h, 2);
    }
}
