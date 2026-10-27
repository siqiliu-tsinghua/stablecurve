//! 变量、单项式与多项式——平坦 `Vec` 实现。
//!
//! 对应原版 `math/{var,monomial,poly}.cc`。原版用"逐节点堆分配的排序链表 +
//! 引用计数 + 对象内游标"（见 REDESIGN §1.3 第 1–3 条）；本实现的项存在排序
//! `Vec` 里，值语义、`iter()` 遍历。本程序里多项式只出现在用户边界（插入类
//! 是 H、psi 与原始类的多项式，实际只有几个单项式），无需更花哨的结构。
//!
//! 语义对齐注意：
//! - 原版多项式可以 UNDEFINED（如除法"除不尽"）；这里用 `Option` 显式表达。
//! - 排序序按原版 `cmp(monomial)` 逐位移植（`monomial.cc:171`）：变量反序、
//!   指数正序、前缀在前。即规范形中变量按 `Var` **降序**（原始类最前、psi 次
//!   之、H 最后）排列。

use std::cmp::Ordering;
use std::fmt;

use crate::num::{One, Q, Zero, q};

/// 变量。对应原版 `var('H')`、`var('c')`（psi）与 `var('P', k)`（第 k 号
/// 原始上同调类，只用于 deg>1 的超曲面）。
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum Var {
    H,
    Psi,
    Prim(u16),
}

impl fmt::Display for Var {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Var::H => write!(f, "H"),
            // 原版内部记 psi 为 'c'，显示沿用，便于与原版日志对照
            Var::Psi => write!(f, "c"),
            Var::Prim(k) => write!(f, "P{k}"),
        }
    }
}

/// 单项式：`(变量, 指数)` 对的规范序列——变量按 [`Var`] 降序、指数升序，
/// 同一变量不重复出现、指数为零的项已剔除。
#[derive(Clone, Default, PartialEq, Eq, Hash, Debug)]
pub struct Mono {
    pairs: Vec<(Var, u8)>,
}

impl Mono {
    /// 常数 1。
    pub fn one() -> Self {
        Mono { pairs: Vec::new() }
    }

    /// 单个变量的一次方。
    pub fn from_var(v: Var) -> Self {
        Mono {
            pairs: vec![(v, 1)],
        }
    }

    pub fn is_one(&self) -> bool {
        self.pairs.is_empty()
    }

    pub fn pairs(&self) -> &[(Var, u8)] {
        &self.pairs
    }

    /// 由任意 (变量, 指数) 集合构造规范形：合并同变量指数、剔除零指数、
    /// 按规范序排列。
    ///
    /// 合并后指数超过 255（u8）时 panic，不截断也不丢项；调用方须先检查
    /// 次数（查询解析器按 `MAX_EXP` 预检）。乘法 [`Mono::mul_mono`] 同此。
    pub fn from_pairs<I: IntoIterator<Item = (Var, u8)>>(pairs: I) -> Self {
        let mut acc: Vec<(Var, u32)> = Vec::new();
        for (v, e) in pairs {
            if e == 0 {
                continue;
            }
            match acc.iter_mut().find(|(w, _)| *w == v) {
                Some((_, fe)) => *fe += e as u32,
                None => acc.push((v, e as u32)),
            }
        }
        if let Some(&(v, e)) = acc.iter().find(|&&(_, e)| e > u8::MAX as u32) {
            panic!("单项式指数溢出：{v}^{e} 超出 u8（≤ 255）");
        }
        // 规范序 = 原版 cmp：变量降序
        acc.sort_by_key(|a| std::cmp::Reverse(a.0));
        Mono {
            pairs: acc.into_iter().map(|(v, e)| (v, e as u8)).collect(),
        }
    }

    /// 原版 `cmp(monomial)` 的直译。
    pub fn cmp_mono(&self, other: &Self) -> Ordering {
        for i in 0..self.pairs.len().max(other.pairs.len()) {
            match (self.pairs.get(i), other.pairs.get(i)) {
                (None, None) => return Ordering::Equal,
                (None, Some(_)) => return Ordering::Less, // 前缀在前
                (Some(_), None) => return Ordering::Greater,
                (Some((vx, ex)), Some((vy, ey))) => {
                    // 变量反序（大变量靠前）
                    match vy.cmp(vx) {
                        Ordering::Equal => {}
                        ne => return ne,
                    }
                    // 指数正序
                    match ex.cmp(ey) {
                        Ordering::Equal => {}
                        ne => return ne,
                    }
                }
            }
        }
        Ordering::Equal
    }

    /// 乘法 = 合并指数。
    pub fn mul_mono(&self, other: &Self) -> Self {
        Mono::from_pairs(
            self.pairs
                .iter()
                .chain(other.pairs.iter())
                .map(|&(v, e)| (v, e as u32 as u8)),
        )
    }

    /// 除法：指数逐项相减；任何一项除不尽则返回 `None`（对应原版
    /// monomial 除法的 UNDEFINED）。
    pub fn div_mono(&self, other: &Self) -> Option<Self> {
        let mut pairs = self.pairs.clone();
        for &(v, e) in &other.pairs {
            let slot = pairs.iter_mut().find(|(w, _)| *w == v)?;
            if slot.1 < e {
                return None;
            }
            slot.1 -= e;
        }
        pairs.retain(|&(_, e)| e > 0);
        Some(Mono { pairs })
    }

    /// 幂：指数相乘；溢出 u8 则 `None`。
    pub fn pow_mono(&self, n: u32) -> Option<Self> {
        let mut pairs = Vec::with_capacity(self.pairs.len());
        for &(v, e) in &self.pairs {
            let e = (e as u32).checked_mul(n)?;
            if e > u8::MAX as u32 {
                return None;
            }
            pairs.push((v, e as u8));
        }
        Some(Mono { pairs })
    }
}

impl PartialOrd for Mono {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Mono {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_mono(other)
    }
}

impl From<Var> for Mono {
    fn from(v: Var) -> Self {
        Mono::from_var(v)
    }
}

impl fmt::Display for Mono {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.pairs.is_empty() {
            return write!(f, "1");
        }
        for &(v, e) in &self.pairs {
            write!(f, "{v}")?;
            if e > 1 {
                write!(f, "{e}")?;
            }
        }
        Ok(())
    }
}

/// 多项式：`(系数, 单项式)` 的规范序列——按 [`Mono`] 的序升序排列，
/// 系数恒非零（空序列即多项式 0）。对应原版 `poly` 类。
#[derive(Clone, Default, PartialEq, Eq, Debug)]
pub struct Poly {
    terms: Vec<(Q, Mono)>,
}

impl Poly {
    pub fn zero() -> Self {
        Poly { terms: Vec::new() }
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn terms(&self) -> &[(Q, Mono)] {
        &self.terms
    }

    /// 规范化：按单项式序排序、合并同类项、剔零系数。
    fn from_terms(terms: Vec<(Q, Mono)>) -> Self {
        let mut terms: Vec<(Q, Mono)> = terms.into_iter().filter(|(c, _)| *c != Q::ZERO).collect();
        terms.sort_by(|a, b| a.1.cmp(&b.1));
        let mut out: Vec<(Q, Mono)> = Vec::with_capacity(terms.len());
        for (c, m) in terms {
            match out.last_mut() {
                Some((pc, pm)) if *pm == m => {
                    *pc = &*pc + &c;
                    if *pc == Q::ZERO {
                        out.pop();
                    }
                }
                _ => out.push((c, m)),
            }
        }
        Poly { terms: out }
    }

    /// 加法（两排序序列归并，`paddsub` 的直译）。
    pub fn add(&self, other: &Self) -> Self {
        Poly::from_terms(
            self.terms
                .iter()
                .chain(other.terms.iter())
                .map(|(c, m)| (c.clone(), m.clone()))
                .collect(),
        )
    }

    pub fn sub(&self, other: &Self) -> Self {
        Poly::from_terms(
            self.terms
                .iter()
                .chain(other.terms.iter())
                .map(|(c, m)| (c.clone(), m.clone()))
                .enumerate()
                .map(|(i, (c, m))| {
                    // 第二个多项式的系数取负
                    if i >= self.terms.len() {
                        (-c, m)
                    } else {
                        (c, m)
                    }
                })
                .collect(),
        )
    }

    pub fn neg(&self) -> Self {
        Poly {
            terms: self.terms.iter().map(|(c, m)| (-c, m.clone())).collect(),
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut terms = Vec::with_capacity(self.terms.len() * other.terms.len());
        for (c1, m1) in &self.terms {
            for (c2, m2) in &other.terms {
                terms.push((c1 * c2, m1.mul_mono(m2)));
            }
        }
        Poly::from_terms(terms)
    }

    /// 乘以标量。
    pub fn scale(&self, s: &Q) -> Self {
        if *s == Q::ZERO {
            return Poly::zero();
        }
        Poly {
            terms: self.terms.iter().map(|(c, m)| (c * s, m.clone())).collect(),
        }
    }

    /// 除以单项式：所有项除尽才有效（对应原版 `pdiv` 的 UNDEFINED 语义）。
    pub fn div_mono(&self, m: &Mono) -> Option<Self> {
        let mut terms = Vec::with_capacity(self.terms.len());
        for (c, t) in &self.terms {
            terms.push((c.clone(), t.div_mono(m)?));
        }
        Some(Poly::from_terms(terms))
    }

    /// 非负整数幂（0 次幂为 1）。
    pub fn pow(&self, n: u32) -> Self {
        let mut result = Poly::from(q(1));
        let mut base = self.clone();
        let mut n = n;
        while n > 0 {
            if n & 1 == 1 {
                result = result.mul(&base);
            }
            n >>= 1;
            if n > 0 {
                base = base.mul(&base);
            }
        }
        result
    }
}

impl From<Var> for Poly {
    fn from(v: Var) -> Self {
        Poly {
            terms: vec![(q(1), Mono::from_var(v))],
        }
    }
}

impl From<Mono> for Poly {
    fn from(m: Mono) -> Self {
        Poly {
            terms: vec![(q(1), m)],
        }
    }
}

impl From<Q> for Poly {
    fn from(c: Q) -> Self {
        Poly::from_terms(vec![(c, Mono::one())])
    }
}

impl From<i64> for Poly {
    fn from(n: i64) -> Self {
        Poly::from(q(n))
    }
}

impl std::ops::Add for &Poly {
    type Output = Poly;
    fn add(self, rhs: &Poly) -> Poly {
        Poly::add(self, rhs)
    }
}

impl std::ops::Sub for &Poly {
    type Output = Poly;
    fn sub(self, rhs: &Poly) -> Poly {
        Poly::sub(self, rhs)
    }
}

impl std::ops::Mul for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Poly) -> Poly {
        Poly::mul(self, rhs)
    }
}

impl std::ops::Neg for &Poly {
    type Output = Poly;
    fn neg(self) -> Poly {
        Poly::neg(self)
    }
}

impl fmt::Display for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        for (i, (c, m)) in self.terms.iter().enumerate() {
            let neg = *c < Q::ZERO;
            if i == 0 {
                if neg {
                    write!(f, "-")?;
                }
            } else if neg {
                write!(f, "-")?;
            } else {
                write!(f, "+")?;
            }
            let ac = if neg { -c } else { c.clone() };
            if m.is_one() {
                write!(f, "{ac}")?;
            } else if ac == Q::ONE {
                write!(f, "{m}")?;
            } else {
                write!(f, "{ac}{m}")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::q_frac;

    fn h() -> Poly {
        Poly::from(Var::H)
    }
    fn psi() -> Poly {
        Poly::from(Var::Psi)
    }

    #[test]
    fn mono_order_matches_original() {
        // 原版 cmp：变量反序（P < c < H），指数正序，前缀在前
        let hh = Mono::from_var(Var::H);
        let h2 = Mono::from_pairs([(Var::H, 2)]);
        let c = Mono::from_var(Var::Psi);
        let p1 = Mono::from_var(Var::Prim(1));
        assert!(p1 < c);
        assert!(c < hh);
        assert!(hh < h2);
        // 相同变量前缀 + 一个更多变量：多变量者靠前（H*H=H2 前缀情形）
        let hc = Mono::from_pairs([(Var::H, 1), (Var::Psi, 1)]);
        assert!(hc < hh); // 同前缀 (H,1)，hc 是 h 的"后缀更长"……见 cmp：i=1 时 x 无 → x 更小
    }

    #[test]
    fn mono_merge_and_div() {
        let m = Mono::from_pairs([(Var::H, 1), (Var::H, 2), (Var::Psi, 0)]);
        assert_eq!(m.to_string(), "H3");
        let d = m.div_mono(&Mono::from_pairs([(Var::H, 2)])).unwrap();
        assert_eq!(d.to_string(), "H");
        assert!(m.div_mono(&Mono::from_pairs([(Var::H, 4)])).is_none());
    }

    /// 指数溢出 u8 须报错，不得静默丢掉变量（H^200·H^100 曾变成 1）。
    #[test]
    #[should_panic(expected = "单项式指数溢出")]
    fn mono_exponent_overflow_panics() {
        let _ = Mono::from_pairs([(Var::H, 200), (Var::H, 100)]);
    }

    #[test]
    fn poly_basics() {
        let half_h = h().scale(&q_frac(1, 2));
        let p = &half_h + &psi();
        assert_eq!(p.to_string(), "c+1/2H"); // c < H：变量反序
        assert_eq!((&p - &p).to_string(), "0");
        let q_ = &h() + &Poly::from(1);
        assert_eq!(q_.to_string(), "1+H"); // 常数项（空单项式）排最前
        assert_eq!((&q_ * &q_).to_string(), "1+2H+H2");
    }

    #[test]
    fn binomial_expansion() {
        // (5H + psi)^5：6 项，系数 C(5,k) 5^(5-k)
        let base = &h().scale(&q(5)) + &psi();
        let p5 = base.pow(5);
        assert_eq!(p5.terms().len(), 6);
        let terms: std::collections::BTreeMap<String, String> = p5
            .terms()
            .iter()
            .map(|(c, m)| (m.to_string(), c.to_string()))
            .collect();
        assert_eq!(terms["H5"], "3125");
        assert_eq!(terms["cH4"], "3125"); // 单项式显示：变量降序，c 在 H 前
        assert_eq!(terms["c5"], "1");
    }

    #[test]
    fn poly_div() {
        let p = Poly::from_terms(vec![
            (q(1), Mono::from_pairs([(Var::H, 2)])),
            (q(1), Mono::from_pairs([(Var::H, 1), (Var::Psi, 1)])),
        ]);
        let d = p.div_mono(&Mono::from_var(Var::H)).unwrap();
        assert_eq!(d.to_string(), "c+H");
        // 除不尽 → None（原版 UNDEFINED）
        let bad = &p + &psi();
        assert!(bad.div_mono(&Mono::from_var(Var::H)).is_none());
    }

    #[test]
    fn zero_and_one() {
        assert!(Poly::zero().is_zero());
        assert_eq!(Poly::zero().pow(0).to_string(), "1");
        assert_eq!(Poly::zero().pow(3).to_string(), "0");
        assert_eq!((&Poly::zero() + &Poly::from(5)).to_string(), "5");
    }
}
