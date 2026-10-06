//! 基本类上 κ-ψ 多项式的交数求值。
//!
//! 路线：κ_m = π_*ψ_{n+1}^{m+1}（投影公式逐因子展开）把 κ-单项式化为
//! 更多标记上的纯 ψ 单项式，再走 [`crate::tau`]。
//!
//! 覆盖（M2）：g=0 全部（Kontsevich）；g=1 可被 dilaton 约化的部分；
//! 其余（g≥1 的零幂/高幂列表，Faber 由 Zaal TAU 表覆盖）返回 Err 等 M3
//! 的 FZ 归约。非顶权单项式直接贡献 0。

use crate::kppoly::{KappaPsiPolynomial, KpTerm};
use crate::tau::Tau;
use sc_num::{One, Q, Zero};

/// ∫_{M̄_{g,n}} p，p 为 (g,n) 基本类上的 κ-ψ 多项式（单顶点、腿 ⊆ 1..=n）。
pub fn int_fund(g: u32, n: u32, p: &KappaPsiPolynomial) -> Result<Q, String> {
    let top = 3 * g as i64 - 3 + n as i64;
    let mut tau = Tau::new();
    let mut acc = Q::ZERO;
    for t in &p.terms {
        if t.kappa.len() > 1 {
            return Err("int_fund 仅支持基本类上的多项式（单顶点）".into());
        }
        if t.weight() as i64 != top {
            continue; // 非顶权 → 0
        }
        for &(l, _) in &t.psi {
            if l == 0 || l > n {
                return Err(format!("腿 {l} 超出标记范围 1..={n}"));
            }
        }
        let part = t.kappa.first().cloned().unwrap_or_default();
        let mut d = vec![0u32; n as usize];
        for &(l, pw) in &t.psi {
            d[(l - 1) as usize] = pw;
        }
        for &i in &part {
            d.push(i + 1); // κ_i → 新标记上的 ψ^{i+1}
        }
        acc += &t.coeff * tau.eval(&d)?;
    }
    Ok(acc)
}

/// 便捷构造+求值：∫_{M̄_{g,n}} Πκ_i^{k_i} Πψ_j^{p_j}。
pub fn int_fund_monomial(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<Q, String> {
    let mut term = KpTerm {
        kappa: vec![
            kappa
                .iter()
                .flat_map(|&(i, p)| vec![i; p as usize])
                .collect(),
        ],
        psi: psi.to_vec(),
        coeff: Q::ONE,
    };
    term.kappa[0].sort_unstable();
    let p = KappaPsiPolynomial::new(vec![term]);
    int_fund(g, n, &p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_kappa_psi_values() {
        // ∫_{0,4}κ_1 = ∫_{0,5}ψ_5² = 1
        assert_eq!(int_fund_monomial(0, 4, &[(1, 1)], &[]).unwrap(), Q::ONE);
        // ∫_{0,5}κ_1² = ∫_{0,7}ψ_6²ψ_7² = 6
        assert_eq!(int_fund_monomial(0, 5, &[(1, 2)], &[]).unwrap(), Q::from(6));
        // ∫_{0,5}κ_1ψ_1 = 3
        assert_eq!(
            int_fund_monomial(0, 5, &[(1, 1)], &[(1, 1)]).unwrap(),
            Q::from(3)
        );
        // ∫_{0,4}ψ_1 = 1（Kontsevich）
        assert_eq!(int_fund_monomial(0, 4, &[], &[(1, 1)]).unwrap(), Q::ONE);
    }

    #[test]
    fn g1_reachable_part() {
        // ∫_{1,1}ψ = 1/24；∫_{1,2}ψ_1ψ_2 = 1/24（dilaton 链）
        assert_eq!(
            int_fund_monomial(1, 1, &[], &[(1, 1)]).unwrap(),
            Q::from_signeds(1, 24)
        );
        assert_eq!(
            int_fund_monomial(1, 2, &[], &[(1, 1), (2, 1)]).unwrap(),
            Q::from_signeds(1, 24)
        );
        // ∫_{1,1}κ_1 = ∫_{1,2}ψ_2² = 1/24（admcycles 验证；string 全亏格生效）
        assert_eq!(
            int_fund_monomial(1, 1, &[(1, 1)], &[]).unwrap(),
            Q::from_signeds(1, 24)
        );
        // ∫_{1,2}κ_1：权重 1 ≠ 顶 2 → 0（非 Err）
        assert_eq!(int_fund_monomial(1, 2, &[(1, 1)], &[]).unwrap(), Q::ZERO);
    }

    #[test]
    fn weight_filter_blocks_wrong_genus_chain() {
        // 非顶权 → 0（不是走错亏格的 dilaton 链）：∫_{2,2}ψ_1ψ_2 权 2 ≠ 顶 4
        assert_eq!(
            int_fund_monomial(2, 2, &[], &[(1, 1), (2, 1)]).unwrap(),
            Q::ZERO
        );
        // ∫_{2,1}ψ_1：权 1 ≠ 顶 4 → 0
        assert_eq!(int_fund_monomial(2, 1, &[], &[(1, 1)]).unwrap(), Q::ZERO);
        // 顶权 g=2 但全幂 ≥1：dilaton 走到 [4] → FZ 需求
        assert!(int_fund_monomial(2, 1, &[], &[(1, 4)]).is_err());
    }

    #[test]
    fn push_tau_consistency_via_api() {
        // ∫_{0,n}p == ∫_{0,n-1}pushn(p)（与 push.rs 的性质测试同一恒等式，走公开 API）
        use crate::kppoly::KappaPsiPolynomial as KP;
        use crate::push::pushn;
        let p = &(&KP::kappa_power(0, 1, 1) * &KP::psi_power(1, 1)) + &KP::psi_power(2, 2);
        // (0,5)：κ_1ψ_1 权 2 顶 ✓ + ψ_2² 权 2 顶 ✓
        let lhs = int_fund(0, 5, &p).unwrap();
        let q = pushn(0, 5, &p).unwrap();
        let rhs = int_fund(0, 4, &q).unwrap();
        assert_eq!(lhs, rhs);
        assert_eq!(lhs, Q::from(4)); // κ_1ψ_1 → 3，ψ_2² → 2!/2! = 1
    }
}
