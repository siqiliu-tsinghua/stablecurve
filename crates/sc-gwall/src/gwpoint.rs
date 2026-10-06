//! GWpoint / power / N_0：纯点后代 local curve 不变量的闭式层（M1）。
//!
//! `power(part, l)`：Faber 闭式——分拆 λ 的幂和 + z/(2sinh(z/2)) 的
//! Bernoulli 修正（GWall `power` 逐式照抄；λ 升序，负下标按倒数展开）。
//! `gwpoint(g, d, tau)`：⟨Π τ_{k_i}(pt)⟩_{g,d} =
//!   Σ_{λ ⊢ d} (d!/dim λ)^{2g−2} · Π_i power(λ, k_i+1)/(k_i+1)!
//! （Maple combinat `Chi(λ, parts[1])` 只在 parts[1] = [1^d]（升序首项）
//! 处被调用 ⇒ d!/χ = Π 钩长；此假设由 M4 的论文 BPS 表对拍全局背书）。
//! `n0(g, b)`：常映射贡献的 Faber 闭式——与 sc-faber klp 的 λ³ 积分
//! 互证：N_0(g, 2g−4)/(3g−3) = (−1)^{g+1}·∫_{M̄_{g−1}} λ_{g−2}³
//! （h = 2,3,4 三点精确咬合，见测试）。

use crate::comb::partition;
use crate::qseries::{bernoulli, z_over_2sinh_coeff};
use sc_num::{One, Q, Zero};

fn fact(n: u32) -> Q {
    let mut v = Q::ONE;
    for i in 2..=n {
        v *= Q::from(i as i64);
    }
    v
}

/// q^e（整数指数，e ≥ 0）。
fn q_pow(base: &Q, e: u32) -> Q {
    let mut v = Q::ONE;
    for _ in 0..e {
        v *= base.clone();
    }
    v
}

/// 幂和 + Bernoulli 修正（GWall `power`）：
/// Σ_{i=1..len} (λ_{len−i} − i + 1/2)^l − (−i + 1/2)^l + l!·[z^{l+1}] z/(2sinh(z/2))。
pub fn power(part: &[u32], l: u32) -> Q {
    let len = part.len();
    let mut pow = Q::ZERO;
    for i in 1..=len {
        let a = Q::from(part[len - i] as i64) - Q::from(2 * i as i64 - 1) / Q::from(2);
        let b = Q::from(-(2 * i as i64 - 1)) / Q::from(2);
        pow = pow + q_pow(&a, l) - q_pow(&b, l);
    }
    pow + fact(l) * z_over_2sinh_coeff(l + 1)
}

/// Π 钩长 = d!/dim λ（Young 图钩子公式；0 基：
/// h(i,j) = (λ_i−1−j) + (λ'_j−1−i) + 1，λ'_j = #{i : λ_i ≥ j+1}）。
fn hook_product(part: &[u32]) -> Q {
    let mut rows: Vec<u32> = part.to_vec();
    rows.sort_unstable();
    rows.reverse(); // 降序行
    let col_len = |j: usize| -> usize { rows.iter().filter(|&&r| r > j as u32).count() };
    let mut h = Q::ONE;
    for (i, &ri) in rows.iter().enumerate() {
        for j in 0..ri as usize {
            let hook = (ri as i64 - 1 - j as i64) + (col_len(j) as i64 - 1 - i as i64) + 1;
            h *= Q::from(hook);
        }
    }
    h
}

/// gwpoint 的分拆展开 ⟨Π τ_k(pt)⟩_{g,d} = Σ_{λ ⊢ d} w(λ)·Π_k f_λ(k) 中的
/// (w(λ), λ)，w(λ) = (Π 钩长)^{2g−2}；f_λ 见 [`tau_factor`]。
pub(crate) fn partition_weights(g: u32, d: u32) -> Vec<(Q, Vec<u32>)> {
    let e = 2 * g as i64 - 2;
    partition(d)
        .into_iter()
        .map(|part| {
            let base = hook_product(&part);
            let w = if e >= 0 {
                q_pow(&base, e as u32)
            } else {
                Q::ONE / q_pow(&base, (-e) as u32)
            };
            (w, part)
        })
        .collect()
}

/// 单个点插入 τ_k(pt) 在分拆 λ 处的因子 f_λ(k) = power(λ, k+1)/(k+1)!。
pub(crate) fn tau_factor(part: &[u32], k: u32) -> Q {
    power(part, k + 1) / fact(k + 1)
}

/// ⟨Π τ_{k_i}(pt)⟩_{g,d}（不连通；g 为 target 曲线亏格）。GWall `GWpoint`。
pub fn gwpoint(g: u32, d: u32, tau: &[u32]) -> Q {
    let mut answer = Q::ZERO;
    for (w, part) in partition_weights(g, d) {
        let mut prod = w;
        for &k in tau {
            prod *= tau_factor(&part, k);
        }
        answer += prod;
    }
    answer
}

/// 常映射贡献 N_0(g, b)（b 偶）。GWall `N_0`。
pub fn n0(g: u32, b: u32) -> Q {
    let sgn = Q::from(if (b / 2 + 1) % 2 == 0 { 1 } else { -1 });
    let prod = bernoulli(b) * bernoulli(b + 2);
    let absp = if prod >= Q::ZERO {
        prod
    } else {
        Q::ZERO - prod
    };
    sgn * Q::from(3 * g as i64 - 3) * absp / Q::from((b + 2) as i64) / Q::from(b as i64) / fact(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_hand() {
        // l=1：每项差恰为分量值 → Σλ_i − 1/24（[z²] 系数）
        assert_eq!(power(&[1], 1).to_string(), "23/24");
        assert_eq!(power(&[2], 1).to_string(), "47/24");
        assert_eq!(power(&[1, 1], 1).to_string(), "47/24");
        assert_eq!(power(&[1, 2], 1).to_string(), "71/24");
        // l=2：[z³] = 0；(1/2)²−(1/2)² = 0
        assert_eq!(power(&[1], 2).to_string(), "0");
        assert_eq!(power(&[2], 2).to_string(), "2"); // (3/2)²−(1/2)²
        // l=3：1/4 + 3!·7/5760 = 1/4 + 7/960 = 247/960
        assert_eq!(power(&[1], 3).to_string(), "247/960");
    }

    #[test]
    fn hooks_basics() {
        // [1^d]：全部钩长 1..d → Π = d!；dim [1^d] = 1
        assert_eq!(hook_product(&[1, 1, 1]).to_string(), "6");
        // [2,1]（d=3）：钩长 3,1,1,1 → 3；dim [2,1] = 2，3!/2 = 3 ✓
        assert_eq!(hook_product(&[2, 1]).to_string(), "3");
        // [2,2]（d=4）：钩 3,2,2,1 → 12；dim = 2，4!/2 = 12 ✓
        assert_eq!(hook_product(&[2, 2]).to_string(), "12");
    }

    #[test]
    fn gwpoint_d1_degenerate() {
        // d=1：唯一分拆 [1]、钩积 1 → gwpoint = Π power([1],k+1)/(k+1)!
        assert_eq!(gwpoint(0, 1, &[0, 0]).to_string(), "529/576"); // (23/24)²
        assert_eq!(gwpoint(5, 1, &[0]).to_string(), "23/24"); // 指数 (1)^8 = 1
        // τ_1（l=2）power = 0 → 含 τ_1 的度 1 不变量为 0
        assert_eq!(gwpoint(0, 1, &[1, 1, 1]).to_string(), "0");
    }

    #[test]
    fn n0_lambda3_cross() {
        // N_0(g,2g−4)/(3g−3) = ±∫_{M̄_{g−1}}λ_{g−2}³（右值来自 sc-faber klp
        // 独立积分：h=2,3,4 = 1/2880、1/725760、1/43545600）
        let rel = |g: u32| -> Q { n0(g, 2 * g - 4) / Q::from(3 * g as i64 - 3) };
        assert_eq!(rel(3).to_string(), "1/2880");
        assert_eq!(rel(4).to_string(), "-1/725760");
        assert_eq!(rel(5).to_string(), "1/43545600");
        // N_0 绝对值抽查（g=2,b=2 手算：(−1)²·3·|B₄B₂|/(4·2·2) = 3/(30·6·16)·? ）
        // |B4·B2| = 1/180；3·(1/180)/(4·2·2) = 3/2880 = 1/960
        assert_eq!(n0(2, 2).to_string(), "1/960");
    }
}
