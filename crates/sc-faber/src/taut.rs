//! 重言类：Σ (稳定图, κ-ψ 多项式)，及基本代数运算。
//!
//! 对应 admcycles 的 `TautologicalClass`（decstratum 列表）。M0 提供
//! 构造器（基本类 / ψ / κ）与线性代数；乘法与前推在 M1/M2 加入。

use crate::graph::StableGraph;
use crate::kppoly::KappaPsiPolynomial;
use sc_num::{Q, Zero};

/// 装饰 stratum：稳定图 + 其上的 κ-ψ 多项式。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecStratum {
    pub graph: StableGraph,
    pub poly: KappaPsiPolynomial,
}

/// 重言类（规范形：图按编码排序、同图项已合并、空多项式项已删）。
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct TautClass {
    pub terms: Vec<DecStratum>,
}

impl TautClass {
    /// 由项表构造规范形。
    pub fn new(terms: Vec<DecStratum>) -> Self {
        let mut ts = terms;
        ts.sort_by(|a, b| a.graph.cmp(&b.graph));
        let mut merged: Vec<DecStratum> = Vec::with_capacity(ts.len());
        for t in ts {
            if let Some(last) = merged.last_mut() {
                if last.graph == t.graph {
                    last.poly = last.poly.add(&t.poly);
                    continue;
                }
            }
            merged.push(t);
        }
        merged.retain(|t| !t.poly.is_zero());
        TautClass { terms: merged }
    }

    pub fn zero() -> Self {
        TautClass { terms: vec![] }
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// M̄_{g,n} 的基本类（单顶点平凡图）。
    pub fn fundclass(g: u32, n: u32) -> Result<Self, String> {
        if 2 * g as i64 - 2 + n as i64 <= 0 {
            return Err(format!("M̄_({g},{n}) 不稳定"));
        }
        let legs: Vec<u32> = (1..=n).collect();
        let graph = StableGraph::new(&[g], &[legs], &[])?;
        Ok(TautClass::new(vec![DecStratum {
            graph,
            poly: KappaPsiPolynomial::one(),
        }]))
    }

    /// ψ_j ∈ H²(M̄_{g,n})。
    pub fn psi_class(g: u32, n: u32, j: u32) -> Result<Self, String> {
        if j == 0 || j > n {
            return Err(format!("标记号 {j} 超出 1..={n}"));
        }
        Ok(Self::fundclass(g, n)?.term_mul(&KappaPsiPolynomial::psi_power(j, 1)))
    }

    /// κ_i ∈ H^{2i}(M̄_{g,n})。
    pub fn kappa_class(g: u32, n: u32, i: u32) -> Result<Self, String> {
        if i == 0 || i as i64 > 3 * g as i64 - 3 + n as i64 {
            return Err(format!("κ_{i} 超出 M̄_({g},{n}) 的上同调范围"));
        }
        Ok(Self::fundclass(g, n)?.term_mul(&KappaPsiPolynomial::kappa_power(0, i, 1)))
    }

    /// 乘一个 κ-ψ 多项式（同图；公开供求值组合用）。
    pub fn term_mul_pub(&self, p: &KappaPsiPolynomial) -> Self {
        self.term_mul(p)
    }

    fn term_mul(&self, p: &KappaPsiPolynomial) -> Self {
        TautClass::new(
            self.terms
                .iter()
                .map(|t| DecStratum {
                    graph: t.graph.clone(),
                    poly: t.poly.mul(p),
                })
                .collect(),
        )
    }

    /// 求值 ∫_{M̄_{g,n}}（要求余维 = 3g−3+n）。逐顶点 socle 公式
    /// （admcycles `decstratum.evaluate` 同构），再除以图自同构数
    /// （stratum 类约定 [Γ] = ξ_*fund/|Aut|——admcycles 的裸 decstratum
    /// 不除，使用其 irrbdiv 等构造时需注意此差异）。
    pub fn evaluate(&self) -> Result<Q, String> {
        use crate::socle::socle_formula;
        let mut acc = Q::ZERO;
        for t in &self.terms {
            let aut_inv = Q::from_signeds(1, t.graph.automorphism_number() as i64);
            for term in &t.poly.terms {
                let mut v = term.coeff.clone() * &aut_inv;
                for (vi, vtx) in t.graph.vertices().iter().enumerate() {
                    let psilist: Vec<u32> =
                        vtx.legs.iter().map(|&l| term.psi_power_at(l)).collect();
                    let kappalist = term.kappa.get(vi).cloned().unwrap_or_default();
                    let dv: i64 = 3 * vtx.genus as i64 - 3 + vtx.legs.len() as i64;
                    let deg: i64 = psilist.iter().map(|&x| x as i64).sum::<i64>()
                        + kappalist.iter().map(|&x| x as i64).sum::<i64>();
                    if deg != dv {
                        v = Q::ZERO;
                        break;
                    }
                    v *= socle_formula(vtx.genus, &psilist, &kappalist)?;
                }
                acc += v;
            }
        }
        Ok(acc)
    }

    /// 余维数：图余维 + 多项式权重（各项一致时）；零类返回 None。
    pub fn codim(&self) -> Option<u32> {
        let mut it = self
            .terms
            .iter()
            .map(|t| Some(t.graph.codim() + t.poly.weight()?));
        let first = it.next()??;
        it.all(|w| w == Some(first)).then_some(first)
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
            .map(|t| DecStratum {
                graph: t.graph.clone(),
                poly: t.poly.scale(c),
            })
            .collect();
        Self::new(ts)
    }
}

impl std::ops::Add for &TautClass {
    type Output = TautClass;
    fn add(self, rhs: &TautClass) -> TautClass {
        self.add(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::q;

    #[test]
    fn fundclass_basics() {
        let f = TautClass::fundclass(2, 3).unwrap();
        assert_eq!(f.terms.len(), 1);
        assert_eq!(f.terms[0].graph.genus(), 2);
        assert_eq!(f.codim(), Some(0));
        // 不稳定模空间拒绝
        assert!(TautClass::fundclass(0, 2).is_err());
        assert!(TautClass::fundclass(1, 0).is_err());
        // M̄_{0,3}、M̄_{1,1} 合法
        assert!(TautClass::fundclass(0, 3).is_ok());
        assert!(TautClass::fundclass(1, 1).is_ok());
    }

    #[test]
    fn psi_and_kappa_codim() {
        let g = TautClass::fundclass(2, 3).unwrap();
        assert_eq!(TautClass::psi_class(2, 3, 1).unwrap().codim(), Some(1));
        assert_eq!(TautClass::kappa_class(2, 3, 2).unwrap().codim(), Some(2));
        // 越界拒绝：M̄_{2,3} 维数 6 → κ_7 非法、ψ_4 非法
        assert!(TautClass::kappa_class(2, 3, 7).is_err());
        assert!(TautClass::psi_class(2, 3, 4).is_err());
        assert!(TautClass::kappa_class(2, 3, 6).is_ok());
        let _ = g;
    }

    #[test]
    fn linear_structure() {
        // ψ1 + ψ2 合并为同图两项；2·ψ1 - ψ1 = ψ1；零类
        let p1 = TautClass::psi_class(1, 2, 1).unwrap();
        let p2 = TautClass::psi_class(1, 2, 2).unwrap();
        let sum = &p1 + &p2;
        assert_eq!(sum.terms.len(), 1);
        assert_eq!(sum.terms[0].poly.terms.len(), 2);
        let scaled = p1.scale(&q(2));
        let diff = scaled.add(&p1.scale(&q(-1)));
        assert_eq!(diff, p1);
        assert!(TautClass::zero().is_zero());
        assert_eq!(TautClass::zero().codim(), None);
    }

    #[test]
    fn evaluate_strata_and_fundclass() {
        use crate::generators::graphs;
        // ∫_{M̄_{1,1}} Δ_0 = 1/2（经典值）
        let g0 = &graphs(1, 1, 1)[0];
        let d0 = TautClass::new(vec![DecStratum {
            graph: g0.clone(),
            poly: KappaPsiPolynomial::one(),
        }]);
        assert_eq!(d0.evaluate().unwrap(), Q::from_signeds(1, 2));
        // ∫κ_1³ M̄_2 = 43/2880；∫ψ_1⁴ M̄_{2,1} = 1/1152（M3 真值表）
        let k = TautClass::kappa_class(2, 0, 1).unwrap();
        let k3 = k
            .term_mul_pub(&KappaPsiPolynomial::kappa_power(0, 1, 1))
            .term_mul_pub(&KappaPsiPolynomial::kappa_power(0, 1, 1));
        assert_eq!(k3.evaluate().unwrap(), Q::from_signeds(43, 2880));
        let p = TautClass::psi_class(2, 1, 1).unwrap();
        let p4 = p
            .term_mul_pub(&KappaPsiPolynomial::psi_power(1, 1))
            .term_mul_pub(&KappaPsiPolynomial::psi_power(1, 1))
            .term_mul_pub(&KappaPsiPolynomial::psi_power(1, 1));
        assert_eq!(p4.evaluate().unwrap(), Q::from_signeds(1, 1152));
        // 非顶维类 → 0
        assert_eq!(
            TautClass::psi_class(2, 1, 1).unwrap().evaluate().unwrap(),
            Q::ZERO
        );
    }

    #[test]
    fn mixed_homogeneity_detected() {
        // κ1 + ψ1 齐次（余维 1）；κ1 + κ2 不齐次 → codim None
        let a = TautClass::kappa_class(1, 2, 1).unwrap();
        let b = TautClass::psi_class(1, 2, 1).unwrap();
        assert_eq!((&a + &b).codim(), Some(1));
        let c = TautClass::kappa_class(1, 2, 2).unwrap();
        assert_eq!((&a + &c).codim(), None);
    }
}
