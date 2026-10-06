//! 基本类上 κ-ψ 多项式的交数求值。
//!
//! 路线：逐项走 [`crate::socle::socle_formula`]——κ-单项式经集合分划求和
//! 化为纯 ψ 交数（κ 的乘积须带分划修正，例如
//! π_*(ψ_{n+1}²ψ_{n+2}²) = κ₁² + κ₂），g=0 用 Kontsevich 多项系数，
//! g≥1 用纯 ψ 顶度递归。覆盖全部 (g,n)；非顶权单项式直接贡献 0。
//!
//! 早期版本（M2）把每个 κ_i 独立换成新标记上的 ψ^{i+1} 再走 [`crate::tau`]，
//! 对两个以上 κ 因子的单项式是错的（如给出 ∫_{0,5}κ₁² = 6，应为 5）。

use crate::kppoly::{KappaPsiPolynomial, KpTerm};
use crate::socle::socle_formula;
use sc_num::{One, Q, Zero};

/// ∫_{M̄_{g,n}} p，p 为 (g,n) 基本类上的 κ-ψ 多项式（单顶点、腿 ⊆ 1..=n）。
pub fn int_fund(g: u32, n: u32, p: &KappaPsiPolynomial) -> Result<Q, String> {
    if 2 * g + n < 3 {
        return Err(format!("M̄_{{{g},{n}}} 不稳定"));
    }
    let top = 3 * g as i64 - 3 + n as i64;
    let mut acc = Q::ZERO;
    for t in &p.terms {
        if t.kappa.len() > 1 {
            return Err("int_fund 仅支持基本类上的多项式（单顶点）".into());
        }
        for &(l, _) in &t.psi {
            if l == 0 || l > n {
                return Err(format!("腿 {l} 超出标记范围 1..={n}"));
            }
        }
        if t.weight() as i64 != top {
            continue; // 非顶权 → 0
        }
        let psilist: Vec<u32> = (1..=n).map(|l| t.psi_power_at(l)).collect();
        let kappalist = t.kappa.first().cloned().unwrap_or_default();
        acc += &t.coeff * socle_formula(g, &psilist, &kappalist)?;
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
        // ∫_{0,5}κ_1² = ⟨τ₂²τ₀⁵⟩ − ⟨τ₃τ₀⁵⟩ = 6 − 1 = 5（= 2!·V_{0,5}/(2π²)²，V_{0,5} = 10π⁴）
        assert_eq!(int_fund_monomial(0, 5, &[(1, 2)], &[]).unwrap(), Q::from(5));
        // ∫_{0,6}κ_1³ = 61（V_{0,6} = 244π⁶/3）；∫_{0,6}κ_1κ_2 = 10 − 1 = 9
        assert_eq!(
            int_fund_monomial(0, 6, &[(1, 3)], &[]).unwrap(),
            Q::from(61)
        );
        assert_eq!(
            int_fund_monomial(0, 6, &[(1, 1), (2, 1)], &[]).unwrap(),
            Q::from(9)
        );
        // ∫_{0,5}κ_1ψ_1 = 3
        assert_eq!(
            int_fund_monomial(0, 5, &[(1, 1)], &[(1, 1)]).unwrap(),
            Q::from(3)
        );
        // ∫_{0,4}ψ_1 = 1（Kontsevich）
        assert_eq!(int_fund_monomial(0, 4, &[], &[(1, 1)]).unwrap(), Q::ONE);
    }

    #[test]
    fn g1_values() {
        // ∫_{1,1}ψ = 1/24；∫_{1,2}ψ_1ψ_2 = 1/24（dilaton 链）
        assert_eq!(
            int_fund_monomial(1, 1, &[], &[(1, 1)]).unwrap(),
            Q::from_signeds(1, 24)
        );
        assert_eq!(
            int_fund_monomial(1, 2, &[], &[(1, 1), (2, 1)]).unwrap(),
            Q::from_signeds(1, 24)
        );
        // ∫_{1,1}κ_1 = 1/24（admcycles 验证）
        assert_eq!(
            int_fund_monomial(1, 1, &[(1, 1)], &[]).unwrap(),
            Q::from_signeds(1, 24)
        );
        // ∫_{1,2}κ_1：权重 1 ≠ 顶 2 → 0（非 Err）
        assert_eq!(int_fund_monomial(1, 2, &[(1, 1)], &[]).unwrap(), Q::ZERO);
    }

    #[test]
    fn higher_genus_matches_golden() {
        // 与 tests/golden_m3.rs 中 admcycles 对拍值一致
        assert_eq!(
            int_fund_monomial(2, 1, &[], &[(1, 4)]).unwrap(),
            Q::from_signeds(1, 1152)
        );
        assert_eq!(
            int_fund_monomial(2, 0, &[(1, 3)], &[]).unwrap(),
            Q::from_signeds(43, 2880)
        );
        assert_eq!(
            int_fund_monomial(2, 0, &[(1, 1), (2, 1)], &[]).unwrap(),
            Q::from_signeds(1, 240)
        );
        assert_eq!(
            int_fund_monomial(3, 0, &[(1, 6)], &[]).unwrap(),
            Q::from_signeds(176557, 107520)
        );
    }

    #[test]
    fn weight_filter_and_errors() {
        // 非顶权 → 0：∫_{2,2}ψ_1ψ_2 权 2 ≠ 顶 4；∫_{2,1}ψ_1 权 1 ≠ 顶 4
        assert_eq!(
            int_fund_monomial(2, 2, &[], &[(1, 1), (2, 1)]).unwrap(),
            Q::ZERO
        );
        assert_eq!(int_fund_monomial(2, 1, &[], &[(1, 1)]).unwrap(), Q::ZERO);
        // 不稳定空间、越界的腿 → Err
        assert!(int_fund_monomial(1, 0, &[], &[]).is_err());
        assert!(int_fund_monomial(0, 4, &[], &[(5, 1)]).is_err());
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
