//! 遗忘映射前推 π_*：Faber `pushn/push` 的移植。
//!
//! `pushn(g,n,p)` 把 (g,n) 基本类上的 κ-ψ 多项式 p（单顶点、腿=标记
//! 1..n）沿遗忘第 n 个标记的 π: M̄_{g,n}→M̄_{g,n−1} 前推，满足
//! ∫_{g,n} p = ∫_{g,n−1} pushn(p)（**g=0 严格成立**，M2 验收即此）。
//!
//! 步骤（KaLa5 `pushn`）：
//! 1. 替换 κ_j → κ_j + ψ_n^j（拉回关系 κ_j = π^*κ'_j + ψ_n^j）；
//! 2. 拆出无 ψ_n 的部分 q0 与有 ψ_n 的部分 q1；
//! 3. q1：ψ_n^k → κ_{k−1}（k≥2），ψ_n¹ → 标量 2g−3+n（dilaton）；
//! 4. q0：string 方程——对每个整除 ψ_j（j<n）的单项式加 a/ψ_j。
//!
//! ⚠️ 第 4 步的 string 方程仅在 g=0 无修正成立；g≥1 有亏格分裂修正项。
//! Faber 原版无条件应用（其 ch 递归在 g≥2 使用 push）——数值上其验证例
//! 全对，但逐步语义待 M4 搭好 admcycles 后对拍定论。本实现如实移植，
//! 调用方负责 g=0 场景的正确性边界。

use crate::kppoly::{KappaPsiPolynomial, KpTerm};
use sc_num::Q;

/// 校验：单顶点、腿 ⊆ 1..=n。
fn check_fund(g: u32, n: u32, p: &KappaPsiPolynomial) -> Result<(), String> {
    for t in &p.terms {
        if t.kappa.len() > 1 {
            return Err("pushn 仅支持基本类上的多项式（单顶点）".into());
        }
        for &(l, _) in &t.psi {
            if l == 0 || l > n {
                return Err(format!("腿 {l} 超出标记范围 1..={n}"));
            }
        }
    }
    let _ = g;
    Ok(())
}

/// 把单项式的 κ 分划的每个因子独立地展开为 κ 或 ψ_n^幂。
/// 相同指标按重数合并：κ_i^e → Σ_k C(e,k)·κ_i^{e−k}·ψ_n^{ik}，
/// 项数 ∏(e+1)，而非 2^{分划长度}。
fn expand_kappa_psi_n(term: &KpTerm, n: u32, out: &mut Vec<KpTerm>) {
    let part = term.kappa.first().cloned().unwrap_or_default();
    let mut groups: Vec<(u32, u32)> = Vec::new();
    for &i in &part {
        match groups.last_mut() {
            Some(e) if e.0 == i => e.1 += 1,
            _ => groups.push((i, 1)),
        }
    }
    // 递归：对每个 (指标, 重数) 决定有几个因子转成 ψ_n^指标
    fn rec(
        term_coeff: &Q,
        psi: &[(u32, u32)],
        groups: &[(u32, u32)],
        i: usize,
        keep: &mut Vec<u32>,
        psi_n_pow: u32,
        n: u32,
        out: &mut Vec<KpTerm>,
    ) {
        if i == groups.len() {
            let kappa = vec![keep.clone()];
            let mut psi = psi.to_vec();
            if psi_n_pow > 0 {
                // 同腿幂必须合并（否则 psi_power_at 误读分支）
                if let Some(e) = psi.iter_mut().find(|e| e.0 == n) {
                    e.1 += psi_n_pow;
                } else {
                    psi.push((n, psi_n_pow));
                    psi.sort_unstable();
                }
            }
            out.push(KpTerm {
                kappa,
                psi,
                coeff: term_coeff.clone(),
            });
            return;
        }
        let (idx, e) = groups[i];
        let mut binom = Q::from(1);
        for k in 0..=e {
            // 留 e−k 个 κ_idx，转 k 个为 ψ_n^idx
            let kept = keep.len();
            keep.extend(std::iter::repeat_n(idx, (e - k) as usize));
            let c = term_coeff * &binom;
            rec(&c, psi, groups, i + 1, keep, psi_n_pow + idx * k, n, out);
            keep.truncate(kept);
            binom = binom * Q::from(e - k) / Q::from(k + 1);
        }
    }
    rec(
        &term.coeff,
        &term.psi,
        &groups,
        0,
        &mut Vec::new(),
        0,
        n,
        out,
    );
}

/// 遗忘第 n 个标记的前推（见模块注释）。
pub fn pushn(g: u32, n: u32, p: &KappaPsiPolynomial) -> Result<KappaPsiPolynomial, String> {
    check_fund(g, n, p)?;
    if p.is_zero() {
        return Ok(KappaPsiPolynomial::zero());
    }
    // 1. κ_j → κ_j + ψ_n^j
    let mut expanded: Vec<KpTerm> = Vec::new();
    for t in &p.terms {
        expand_kappa_psi_n(t, n, &mut expanded);
    }
    // 2/3. 有 ψ_n 的部分：ψ_n^k → κ_{k−1}（k≥2）或标量（k=1）
    let mut result: Vec<KpTerm> = Vec::new();
    let mut q0: Vec<KpTerm> = Vec::new();
    for mut t in expanded {
        let p_n = t.psi_power_at(n);
        if p_n == 0 {
            t.psi.retain(|&(l, _)| l != n);
            q0.push(t);
        } else {
            t.psi.retain(|&(l, _)| l != n);
            if p_n >= 2 {
                t.kappa[0].push(p_n - 1);
                t.kappa[0].sort_unstable();
            } else {
                t.coeff *= Q::from(2 * g + n - 3); // dilaton
            }
            result.push(t);
        }
    }
    // 4. q0 的 string 方程（g=0 语义；见模块注释）。
    //    与 Faber 一致：q0 单项本身保留在结果中（下游过权自然积分为 0），
    //    另加每个可整除 ψ_j 的 string 项。
    for t in q0 {
        result.push(t.clone());
        for j in 1..n {
            let pj = t.psi_power_at(j);
            if pj == 0 {
                continue;
            }
            let mut a = t.clone();
            for e in a.psi.iter_mut() {
                if e.0 == j {
                    e.1 -= 1;
                }
            }
            a.psi.retain(|&(_, p)| p > 0);
            result.push(a);
        }
    }
    Ok(KappaPsiPolynomial::new(result))
}

/// 逐步遗忘到 (g,0)：Faber `push`。
pub fn push(g: u32, n: u32, p: &KappaPsiPolynomial) -> Result<KappaPsiPolynomial, String> {
    let mut q = p.clone();
    for m in (1..=n).rev() {
        q = pushn(g, m, &q)?;
    }
    Ok(q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kppoly::KappaPsiPolynomial as KP;
    use sc_num::{One, Zero};

    fn kpw(i: u32, p: u32) -> KP {
        KP::kappa_power(0, i, p)
    }
    fn psw(l: u32, p: u32) -> KP {
        KP::psi_power(l, p)
    }

    /// 参考求值（g=0）：走 socle（κ 乘积带集合分划修正；过权项记 0）。
    fn int0(n: u32, p: &KP) -> Q {
        crate::int::int_fund(0, n, p).unwrap()
    }

    #[test]
    fn pushn_psi_dilaton_scalar() {
        // g=0, n=4：pushn(ψ_4) = 2·0+4−3 = 1（dilaton 标量）
        let q = pushn(0, 4, &psw(4, 1)).unwrap();
        assert_eq!(q, KP::one());
        // g=1, n=2：pushn(ψ_2ψ_1) = 1·ψ_1（∫ 应 = 1/24 = ∫_{1,2}ψ_1ψ_2）
        let q = pushn(1, 2, &(&psw(1, 1) * &psw(2, 1))).unwrap();
        assert_eq!(q, psw(1, 1));
    }

    #[test]
    fn pushn_kappa_to_psi_conversion() {
        // κ_2 的 ψ 分支：ψ_4³ → κ_2（g=0,n=4）
        // 直接构造 ψ_4³ 输入验证转换：pushn(0,4, ψ_4³) = κ_2
        let q = pushn(0, 4, &psw(4, 3)).unwrap();
        assert_eq!(q, kpw(2, 1));
        // ψ_4² → κ_1
        let q = pushn(0, 4, &psw(4, 2)).unwrap();
        assert_eq!(q, kpw(1, 1));
    }

    #[test]
    fn pushn_string_g0() {
        // g=0：pushn(0,4, ψ_1) = 1 + ψ_1（string 项 1；q0 的 ψ_1 按 Faber
        // 语义保留，其下方 (0,3) 过权积分为 0，不影响 ∫ 恒等式）
        let q = pushn(0, 4, &psw(1, 1)).unwrap();
        assert_eq!(q, &KP::one() + &psw(1, 1));
        // pushn(0,5, ψ_1ψ_2) = ψ_1ψ_2 + ψ_1 + ψ_2（q0 保留项 + 两个 string 项；
        // 保留项在 (0,4) 过权、积分贡献 0——与 Faber 的 q 完全一致）
        let q = pushn(0, 5, &(&psw(1, 1) * &psw(2, 1))).unwrap();
        let expect = &(&psw(1, 1) * &psw(2, 1)) + &psw(1, 1);
        let expect = &expect + &psw(2, 1);
        assert_eq!(q, expect);
    }

    #[test]
    fn push_consistency_g0_random() {
        // 性质测试：∫_{0,n} p == ∫_{0,n-1} pushn(p)（g=0，随机多项式）
        let mut seed = 0x5eedu64;
        let mut rnd = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        for case in 0..40 {
            let n = 4 + rnd() % 3; // n ∈ 4..6（pushn 目标 ≥ M̄_{0,3} 稳定）
            let w = n - 3; // 顶权
            let mut terms = Vec::new();
            for _ in 0..3 {
                // 随机 ψ-单项式（总幂 w），带随机小系数
                let mut d = vec![0u32; n as usize];
                for _ in 0..w {
                    d[(rnd() % n) as usize] += 1;
                }
                let c = Q::from_signeds((rnd() % 7) as i64 - 3, 1 + (rnd() % 3) as i64);
                let psi: Vec<(u32, u32)> = d
                    .iter()
                    .enumerate()
                    .filter(|&(_, &p)| p > 0)
                    .map(|(i, &p)| (i as u32 + 1, p))
                    .collect();
                terms.push(KpTerm {
                    kappa: vec![vec![]],
                    psi,
                    coeff: c,
                });
            }
            let p = KP::new(terms);
            let lhs = int0(n, &p);
            let q = pushn(0, n, &p).unwrap();
            let rhs = int0(n - 1, &q);
            assert_eq!(lhs, rhs, "case={case} n={n} p={p:?}");
        }
    }

    #[test]
    fn kappa_to_psi_expansion_values() {
        // ∫_{0,4} κ_1 = ∫_{0,5} ψ_5² = 2!/2! = 1
        assert_eq!(int0(4, &kpw(1, 1)), Q::ONE);
        // ∫_{0,5} κ_1² = ∫_{0,7} ψ_6²ψ_7² − ∫_{0,6} ψ_6³ = 6 − 1 = 5
        // （π_*(ψ_6²ψ_7²) = κ_1² + κ_2，κ 乘积不能逐因子换成 ψ）
        assert_eq!(int0(5, &(&kpw(1, 1) * &kpw(1, 1))), Q::from(5));
        // ∫_{0,5} κ_2 = ∫_{0,6} ψ_6³ = 3!/3! = 1
        assert_eq!(int0(5, &kpw(2, 1)), Q::ONE);
        // ∫_{0,5} κ_1 ψ_1 = ∫_{0,6} ψ_6²ψ_1 = 3!/2! = 3
        assert_eq!(int0(5, &(&kpw(1, 1) * &psw(1, 1))), Q::from(3));
        // 超度 → 0：κ_1 在 M̄_{0,3}（点）上
        assert_eq!(int0(3, &kpw(1, 1)), Q::ZERO);
    }

    #[test]
    fn pushn_repeated_kappa_indices() {
        // 回归：κ_i^e 按重数展开（二项式权重）后，∫ 恒等式须保持
        // WP 锚点：∫_{0,5}κ₁² = 5，∫_{0,6}κ₁³ = 61，∫_{0,7}κ₁⁴ = 1379
        assert_eq!(int0(5, &kpw(1, 2)), Q::from(5));
        assert_eq!(int0(6, &kpw(1, 3)), Q::from(61));
        assert_eq!(int0(7, &kpw(1, 4)), Q::from(1379));
        // (0,7) 顶权 4；(0,8) 顶权 5
        let cases: Vec<(u32, KP)> = vec![
            (7, kpw(1, 4)),
            (7, &kpw(1, 2) * &kpw(2, 1)),
            (7, &kpw(1, 3) * &psw(1, 1)),
            (7, &(&kpw(1, 2) * &psw(1, 1)) * &psw(7, 1)),
            (7, &kpw(2, 2) + &(&kpw(1, 2) * &psw(7, 2))),
            (8, &(&kpw(1, 3) * &kpw(2, 1)) + &(&kpw(1, 2) * &psw(8, 3))),
        ];
        for (n, p) in &cases {
            let q = pushn(0, *n, p).unwrap();
            assert_eq!(int0(*n, p), int0(n - 1, &q), "n={n} p={p:?}");
        }
    }
}
