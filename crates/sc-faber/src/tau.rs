//! 纯 ψ 类交数，复用 Hodge 引擎的 string/dilaton 与 DVV 递归。

use sc_num::{Q, Zero};

/// τ 求值入口；同一线程内共享后代积分结果。
#[derive(Default)]
pub struct Tau;

impl Tau {
    pub fn new() -> Self {
        Self
    }

    /// 幂列表允许含 0；由 Σd = 3g−3+n 反解亏格。
    /// 非整数、负亏格或不稳定输入返回 0。
    pub fn eval(&mut self, d: &[u32]) -> Result<Q, String> {
        let triple_g = d.iter().map(|&x| i64::from(x)).sum::<i64>() - d.len() as i64 + 3;
        if triple_g < 0 || triple_g % 3 != 0 || triple_g / 3 > i64::from(u32::MAX) {
            return Ok(Q::ZERO);
        }
        Ok(crate::hodge::psi_integral((triple_g / 3) as u32, d))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::One;

    fn t(d: &[u32]) -> Q {
        Tau::new().eval(d).unwrap()
    }

    #[test]
    fn kontsevich_basics() {
        // ∫_{0,n} ψ_1^{n-3} = 1 对一切 n ≥ 3
        for n in 3..=8u32 {
            let mut d = vec![0u32; (n - 1) as usize];
            d.push(n - 3);
            assert_eq!(t(&d), Q::ONE, "n={n}");
        }
        // Kontsevich 离散值抽查：(n−3)!/∏d_i!
        assert_eq!(t(&[2, 0, 0, 0, 0]), Q::ONE); // 2!/2!
        assert_eq!(t(&[2, 1, 0, 0, 0, 0]), Q::from(3)); // n=6：3!/(2!·1!)
        assert_eq!(t(&[1, 1, 1, 0, 0, 0]), Q::from(6)); // n=6：3!/1
        assert_eq!(t(&[2, 2, 0, 0, 0, 0, 0]), Q::from(6)); // n=7：4!/(2!·2!)
        // ∫_{0,4} ψ_1ψ_2：Σ=2 > dim=1，过度数 → 0
        assert_eq!(t(&[1, 1, 0, 0]), Q::ZERO);
        // ∫_{0,4} ψ_1 = 1
        assert_eq!(t(&[1, 0, 0, 0]), Q::ONE);
    }

    #[test]
    fn dilaton_g1() {
        assert_eq!(t(&[1]), Q::from_signeds(1, 24)); // ∫_{1,1}ψ = 1/24
        assert_eq!(t(&[1, 1]), Q::from_signeds(1, 24)); // ∫_{1,2}ψ_1ψ_2
        assert_eq!(t(&[1, 1, 1]), Q::from_signeds(1, 12)); // 系数 2g+n−3 = 2
        // 一般 [1,...,1]（n 个）= (n−1)!/24
        for n in 1..=5u32 {
            let d = vec![1u32; n as usize];
            let expect = Q::from_signeds(1, 24) * Q::from((1..n).product::<u32>());
            assert_eq!(t(&d), expect, "n={n}");
        }
    }

    #[test]
    fn string_all_genera_admcycles_verified() {
        // 以下真值全部来自 admcycles（WSL oracle，2026-10）
        assert_eq!(t(&[2, 0]), Q::from_signeds(1, 24)); // ∫_{1,2}ψ_1²（= ∫κ_1）
        assert_eq!(t(&[1, 0, 2]), Q::from_signeds(1, 12)); // [0,1,2]
        assert_eq!(t(&[1, 0, 3, 0]), Q::from_signeds(1, 8)); // [0,1,3] on (1,4)
        assert_eq!(t(&[3, 0, 0]), Q::from_signeds(1, 24)); // [0,3] on (1,3)
        assert_eq!(t(&[2, 2, 0, 0]), Q::from_signeds(1, 6)); // [0,0,2,2] on (1,4)
        assert_eq!(t(&[1, 1, 0, 2]), Q::from_signeds(1, 4)); // [0,1,1,2] = 1/4
        assert_eq!(t(&[0, 0, 0, 2, 3]), Q::from_signeds(7, 24)); // [0,2,3] on (1,5)
    }

    #[test]
    fn higher_genus_descendants() {
        assert_eq!(t(&[4]), Q::from_signeds(1, 1152));
        assert_eq!(t(&[5, 0]), Q::from_signeds(1, 1152));
        assert_eq!(t(&[1, 0, 5]), Q::from_signeds(1, 288));
        assert_eq!(t(&[0; 6]), Q::ZERO);
    }

    #[test]
    fn degree_mismatch_is_zero() {
        assert_eq!(t(&[2, 0, 0, 0]), Q::ZERO); // Σd=2 ≠ dim M̄_{0,4}=1
        assert_eq!(t(&[3, 0, 0, 0]), Q::ZERO); // Σd=3 ≠ dim M̄_{1,4}=4 且非 g 整数
        assert_eq!(t(&[3, 0, 0]), Q::from_signeds(1, 24)); // (1,3) 合法度数，string 可解
        assert_eq!(t(&[1, 1, 0, 0, 0]), Q::from(2)); // Σ=2 = dim M̄_{0,5}：Kontsevich 2!/(1!1!)
        assert_eq!(t(&[]), Q::ZERO);
    }
}
