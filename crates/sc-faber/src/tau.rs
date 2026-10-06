//! 纯 ψ 类交数 ∫_{M̄_{g,n}} Πψ_i^{d_i}——Faber `tau` 过程的移植（修正版）。
//!
//! 递归（输入幂列表 x，g 由 Σd = 3g−3+n 反解；度数不符 → 0）：
//! - g=0：Kontsevich 公式 (n−3)!/∏d_i!（覆盖一切合法输入）；
//! - x[1]=1（g≥1）：dilaton 方程 (2g−2+n−1)·tau(x 去掉首元素)，
//!   特例 (g,n)=(1,1) → 1/24；
//! - 其余（g≥2 且所有 d_i≥2——g=1 由鸽笼原理不出现此情形）：Faber 用
//!   Zaal 的 TAU 表 → 本实现走 FZ（M3），当前返回 Err。
//!
//! string 方程（消去零幂标记）对**一切亏格**有效（admcycles 双点验证：
//! [0,1,2]_{(1,3)} = 1/12、[0,1,5]_{(2,3)} = 1/288 = string 预测精确一致），
//! 与 Faber 原版的无守卫实现一致。开发中曾两度怀疑其为 g=0 专属（第一轮
//! 记错 Mumford 关系推出 ∫ψ²=1/2，第二轮手推漏项），均被 admcycles 裁决
//! 推翻——详见 FABER_REDESIGN 的对拍记录。

use sc_num::{One, Q, Zero};
use std::collections::HashMap;

/// τ 求值器（记忆化）。
pub struct Tau {
    memo: HashMap<(u32, Vec<u32>), Q>,
}

fn factorial(n: u32) -> Q {
    let mut r = Q::ONE;
    for i in 2..=n {
        r *= Q::from(i);
    }
    r
}

impl Tau {
    pub fn new() -> Self {
        Tau {
            memo: HashMap::new(),
        }
    }

    /// 求 ∫_{M̄_{g,n}} Πψ_i^{d_i}；`d` 为各标记的幂（可含 0）。
    /// 度数不符（Σd ≠ 3g−3+n）返回 0。
    pub fn eval(&mut self, d: &[u32]) -> Result<Q, String> {
        let mut x = d.to_vec();
        x.sort_unstable();
        let n = x.len() as u32;
        let sum: u32 = x.iter().sum();
        let s = sum as i64 - n as i64;
        let g = if s % 3 == 0 {
            (s / 3 + 1) as u32
        } else {
            return Ok(Q::ZERO); // 度数不符（g 非整数）→ 0
        };
        self.tau_sorted(g, n, &x)
    }

    fn tau_sorted(&mut self, g: u32, n: u32, x: &[u32]) -> Result<Q, String> {
        let key = (g, x.to_vec());
        if let Some(v) = self.memo.get(&key) {
            return Ok(v.clone());
        }
        let v = self.compute(g, n, x)?;
        self.memo.insert(key, v.clone());
        Ok(v)
    }

    fn compute(&mut self, g: u32, n: u32, x: &[u32]) -> Result<Q, String> {
        if g == 0 {
            // Kontsevich：(n−3)!/∏d_i!
            if n < 3 {
                return Ok(Q::ZERO);
            }
            let mut den = Q::ONE;
            for &di in x {
                den *= factorial(di);
            }
            return Ok(factorial(n - 3) / den);
        }
        if x.is_empty() {
            return Ok(Q::ZERO); // 空积式（防御）
        }
        if x[0] == 0 {
            // string（一切亏格）：Σ_{j: x_j>0} tau(删一个零、x_j−1)
            let mut acc = Q::ZERO;
            for j in 0..x.len() {
                if x[j] == 0 {
                    continue;
                }
                let mut y: Vec<u32> = Vec::with_capacity(x.len() - 1);
                let mut zero_skipped = false;
                for (i, &xi) in x.iter().enumerate() {
                    if i == j {
                        y.push(xi - 1);
                    } else if xi == 0 && !zero_skipped {
                        zero_skipped = true;
                    } else {
                        y.push(xi);
                    }
                }
                y.sort_unstable();
                acc += self.tau_sorted(g, n - 1, &y)?;
            }
            return Ok(acc);
        }
        if x[0] == 1 {
            // dilaton：(2g−2+n−1)·tau(x 去掉首元素)
            if g == 1 && n == 1 {
                return Ok(Q::from_signeds(1, 24));
            }
            let rest = &x[1..];
            let c = Q::from(2 * g + n - 3);
            let inner = self.tau_sorted(g, n - 1, rest)?;
            return Ok(c * inner);
        }
        // g≥1 且所有 d_i≥2：Faber 走 TAU 表 → 待 FZ（M3）
        Err(format!(
            "∫ψ^{x:?}（g={g}, n={n}）：所有幂 ≥2，需要 FZ 归约（M3）"
        ))
    }
}

impl Default for Tau {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn needs_fz_errors() {
        // g=1 已被 string+dilaton 完全覆盖（鸽笼：顶度列表 Σd=n 必含 0 或 1）；
        // 待 FZ 的是 g≥2 的"全幂 ≥2"列表：
        assert!(Tau::new().eval(&[4]).is_err()); // [4] on (2,1)：真值 1/1152
        assert!(Tau::new().eval(&[5, 0]).is_err()); // [0,5] on (2,2)：真值 1/1152（string→[4]）
        assert!(Tau::new().eval(&[1, 0, 5]).is_err()); // [0,1,5] on (2,3)：真值 1/288
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
