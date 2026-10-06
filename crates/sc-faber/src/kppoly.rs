//! κ-ψ 多项式：strata 上的示性类代数。
//!
//! 装饰在固定稳定图上：κ 类逐顶点（`kappa[v]` 为顶点 v 的 κ 指标分划，
//! 升序），ψ 类逐腿（`(腿号, 幂)`，标记腿与结点半边皆可，幂为正——零幂
//! 在规范形中删除）。对应 admcycles 的 `KappaPsiPolynomial`（kppoly）。
//!
//! 规范形：键 (kappa, psi) 升序、无重复、系数非零。
//! 权重（= 余维数贡献）：κ_i 贡献 i，ψ_j 贡献 1。
//!
//! 乘法 [`KappaPsiPolynomial::mul`] 仅对同一图上的多项式有意义（顶点数
//! 一致、腿号一致），跨图乘法需要图粘合，在 M3 实现。

use sc_num::{One, Q, is_zero};
use std::cmp::Ordering;

/// 单项式：系数 + (逐顶点 κ 分划, 逐腿 ψ 幂)。
#[derive(Clone, Debug)]
pub struct KpTerm {
    pub kappa: Vec<Vec<u32>>, // kappa[v]：顶点 v 的 κ 指标，升序
    pub psi: Vec<(u32, u32)>, // (腿号, 幂)，按腿号升序，幂 ≥ 1
    pub coeff: Q,
}

impl KpTerm {
    /// 该项的权重：Σ i·k_i + Σ ψ 幂。
    pub fn weight(&self) -> u32 {
        self.kappa.iter().flatten().sum::<u32>() + self.psi.iter().map(|&(_, p)| p).sum::<u32>()
    }

    fn key(&self) -> (&[Vec<u32>], &[(u32, u32)]) {
        (&self.kappa, &self.psi)
    }

    /// ψ 幂（腿号 l 处；未出现为 0）。
    pub fn psi_power_at(&self, l: u32) -> u32 {
        self.psi
            .iter()
            .find(|&&(lg, _)| lg == l)
            .map_or(0, |&(_, p)| p)
    }
}

// 项的比较只看键（κ, ψ）：规范形保证同键项已合并，系数不参与序。
impl PartialEq for KpTerm {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}
impl Eq for KpTerm {}
impl PartialOrd for KpTerm {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for KpTerm {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key().cmp(&other.key())
    }
}

/// ψ 幂表的规范形：按腿号升序，同一条腿的幂相加（ψ_l^a·ψ_l^b = ψ_l^{a+b}），
/// 丢弃零幂。
pub fn merge_psi(psi: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let mut v: Vec<(u32, u32)> = psi.iter().copied().filter(|&(_, p)| p > 0).collect();
    v.sort_unstable();
    let mut out: Vec<(u32, u32)> = Vec::with_capacity(v.len());
    for (l, p) in v {
        match out.last_mut() {
            Some(last) if last.0 == l => last.1 += p,
            _ => out.push((l, p)),
        }
    }
    out
}

/// 用户输入的 ψ 腿号须为标记 1..=n（n+1.. 是结点半边，不对外开放）。
pub fn check_marking_legs(psi: &[(u32, u32)], n: u32) -> Result<(), String> {
    for &(l, _) in psi {
        if l == 0 || l > n {
            return Err(format!("ψ 的标记号 {l} 超出 1..={n}"));
        }
    }
    Ok(())
}

/// κ-ψ 多项式（规范形）。
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct KappaPsiPolynomial {
    pub terms: Vec<KpTerm>,
}

impl KappaPsiPolynomial {
    /// 由任意项表构造规范形：合并同键项、丢弃零系数与零 ψ 幂、按键排序。
    pub fn new(terms: Vec<KpTerm>) -> Self {
        let nv = terms.iter().map(|t| t.kappa.len()).max().unwrap_or(0);
        let mut ts: Vec<KpTerm> = terms
            .into_iter()
            .map(|mut t| {
                t.kappa.resize(nv, vec![]);
                t.psi = merge_psi(&t.psi);
                t
            })
            .collect();
        ts.sort();
        let mut merged: Vec<KpTerm> = Vec::with_capacity(ts.len());
        for t in ts {
            if merged.last().map(|m| m == &t).unwrap_or(false) {
                merged.last_mut().unwrap().coeff += &t.coeff;
            } else {
                merged.push(t);
            }
        }
        merged.retain(|t| !is_zero(&t.coeff));
        KappaPsiPolynomial { terms: merged }
    }

    /// 常值 1（规范形：单顶点 kappa=[[]]，与基本类上多项式一致）。
    pub fn one() -> Self {
        Self::new(vec![KpTerm {
            kappa: vec![vec![]],
            psi: vec![],
            coeff: Q::ONE,
        }])
    }

    pub fn zero() -> Self {
        KappaPsiPolynomial { terms: vec![] }
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// 顶点 v 上的 κ_i^p。
    pub fn kappa_power(v: usize, i: u32, p: u32) -> Self {
        let mut kappa = vec![vec![]; v];
        kappa.push(vec![i; p as usize]);
        Self::new(vec![KpTerm {
            kappa,
            psi: vec![],
            coeff: Q::ONE,
        }])
    }

    /// 腿 l 上的 ψ^p（l 可为标记或结点半边；规范形带单顶点壳）。
    pub fn psi_power(l: u32, p: u32) -> Self {
        Self::new(vec![KpTerm {
            kappa: vec![vec![]],
            psi: vec![(l, p)],
            coeff: Q::ONE,
        }])
    }

    /// 是否齐次；若是返回公共权重。
    pub fn weight(&self) -> Option<u32> {
        let mut it = self.terms.iter().map(KpTerm::weight);
        let first = it.next()?;
        it.all(|w| w == first).then_some(first)
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut ts = self.terms.clone();
        ts.extend(other.terms.iter().cloned());
        Self::new(ts)
    }

    pub fn scale(&self, c: &Q) -> Self {
        let ts = self
            .terms
            .iter()
            .map(|t| KpTerm {
                kappa: t.kappa.clone(),
                psi: t.psi.clone(),
                coeff: &t.coeff * c,
            })
            .collect();
        Self::new(ts)
    }

    /// 同图乘法：κ 逐顶点归并，ψ 逐腿相加，系数相乘。
    pub fn mul(&self, other: &Self) -> Self {
        let mut ts = Vec::with_capacity(self.terms.len() * other.terms.len().max(1));
        for a in &self.terms {
            for b in &other.terms {
                let nv = a.kappa.len().max(b.kappa.len());
                let mut kappa = vec![vec![]; nv];
                for (i, ka) in a.kappa.iter().enumerate() {
                    kappa[i].extend_from_slice(ka);
                }
                for (i, kb) in b.kappa.iter().enumerate() {
                    kappa[i].extend_from_slice(kb);
                    kappa[i].sort_unstable();
                }
                let mut psi = a.psi.clone();
                psi.extend_from_slice(&b.psi);
                psi.sort_unstable();
                let mut merged: Vec<(u32, u32)> = Vec::with_capacity(psi.len());
                for (l, p) in psi {
                    if merged.last().map(|m| m.0 == l).unwrap_or(false) {
                        merged.last_mut().unwrap().1 += p;
                    } else {
                        merged.push((l, p));
                    }
                }
                merged.retain(|&(_, p)| p > 0);
                ts.push(KpTerm {
                    kappa,
                    psi: merged,
                    coeff: &a.coeff * &b.coeff,
                });
            }
        }
        Self::new(ts)
    }
}

impl std::ops::Add for &KappaPsiPolynomial {
    type Output = KappaPsiPolynomial;
    fn add(self, rhs: &KappaPsiPolynomial) -> KappaPsiPolynomial {
        self.add(rhs)
    }
}

impl std::ops::Mul for &KappaPsiPolynomial {
    type Output = KappaPsiPolynomial;
    fn mul(self, rhs: &KappaPsiPolynomial) -> KappaPsiPolynomial {
        self.mul(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::q;

    #[test]
    fn one_and_zero() {
        let one = KappaPsiPolynomial::one();
        assert!(!one.is_zero());
        assert_eq!(one.weight(), Some(0));
        assert!(KappaPsiPolynomial::zero().is_zero());
        assert_eq!(&one * &one, one);
    }

    #[test]
    fn normalization_merges_and_sorts() {
        let p = KappaPsiPolynomial::new(vec![
            KpTerm {
                kappa: vec![vec![1]],
                psi: vec![],
                coeff: q(1),
            },
            KpTerm {
                kappa: vec![vec![1]],
                psi: vec![],
                coeff: q(1),
            },
            KpTerm {
                kappa: vec![vec![2]],
                psi: vec![],
                coeff: q(3),
            },
        ]);
        assert_eq!(p.terms.len(), 2);
        assert_eq!(p.terms[0].kappa[0], vec![1]);
        assert_eq!(p.terms[0].coeff, q(2));
        assert_eq!(p.terms[1].kappa[0], vec![2]);
        let z = KappaPsiPolynomial::new(vec![
            KpTerm {
                kappa: vec![vec![1]],
                psi: vec![],
                coeff: q(5),
            },
            KpTerm {
                kappa: vec![vec![1]],
                psi: vec![],
                coeff: q(-5),
            },
        ]);
        assert!(z.is_zero());
        // 零幂 ψ 被删：ψ_l^0 不出现
        let p = KappaPsiPolynomial::new(vec![KpTerm {
            kappa: vec![vec![]],
            psi: vec![(1, 0), (2, 3)],
            coeff: q(1),
        }]);
        assert_eq!(p.terms[0].psi, vec![(2, 3)]);
    }

    #[test]
    fn kappa_per_vertex() {
        // 顶点 0 的 κ1 · 顶点 1 的 κ2：分顶点存放；同顶点归并为分划
        let a = KappaPsiPolynomial::kappa_power(0, 1, 1);
        let b = KappaPsiPolynomial::kappa_power(1, 2, 1);
        let p = &a * &b;
        assert_eq!(p.terms.len(), 1);
        assert_eq!(p.terms[0].kappa, vec![vec![1], vec![2]]);
        assert_eq!(p.weight(), Some(3));
        let p = &a * &a;
        assert_eq!(p.terms[0].kappa, vec![vec![1, 1]]);
    }

    #[test]
    fn psi_at_nodal_legs() {
        // ψ 可落在结点半边（腿号 > n），同腿幂相加，不同腿共存
        let p1 = KappaPsiPolynomial::psi_power(7, 1);
        let p2 = KappaPsiPolynomial::psi_power(7, 2);
        let p = &p1 * &p2;
        assert_eq!(p.terms[0].psi, vec![(7, 3)]);
        assert_eq!(p.weight(), Some(3));
        let p3 = KappaPsiPolynomial::psi_power(2, 1);
        let p = &p1 * &p3;
        assert_eq!(p.terms[0].psi, vec![(2, 1), (7, 1)]);
    }

    #[test]
    fn homogeneity_and_distributivity() {
        let k1 = KappaPsiPolynomial::kappa_power(0, 1, 1);
        let k2 = KappaPsiPolynomial::kappa_power(0, 2, 1);
        let p1 = KappaPsiPolynomial::psi_power(1, 1);
        assert_eq!((&k1 + &p1).weight(), Some(1));
        assert_eq!((&k1 + &k2).weight(), None);
        let lhs = &k1 + &k2;
        let rhs = &(&p1 * &k1) + &(&p1 * &k2);
        assert_eq!(&(&lhs * &p1), &rhs);
    }
}
