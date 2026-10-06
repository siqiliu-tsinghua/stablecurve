//! GWloop / fastGWloop：环路后代（H¹ 外代数）的消去层（M2）。
//!
//! 数学：a₁..a_g, b₁..b_g 为 H¹(C_g) 的辛基（a_i∪b_j = δ_ij，a∪a = b∪b = 0）。
//! 环路后代反对易 → 不变量 = 置换和 + 反对称符号：
//!   ⟨Πτ(pt)·Π τ_{x_i}(a_{α_i})·τ_{y_i}(b_{β_i})⟩ = Σ_σ sign·Π binom(x+y_σ,x)·cup·⟨纯点⟩
//! 两实现互为内建差分（原版即双实现）：
//!
//! - [`gwloop`]：全 S_m 置换和（慢参考版，GWall `GWloop`）
//! - [`fast_gwloop`]：按环路指标分组 + 组内重复剪枝（GWall `fastGWloop`）
//!
//! 符号统一用 (−1)^{m(m−1)/2}·parity(σ)（原版 evils 的魔数依赖 Maple
//! permute 枚举顺序且作者自注仅 m<9 有效——弃用，见 comb::parity）。

use crate::comb::{CartProd, parity, permute};
use crate::gwpoint::gwpoint;
use sc_num::{One, Q, Zero};

/// 环路插入 τ_k(a_idx) / τ_k(b_idx)（idx ∈ 1..=g）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LoopIns {
    pub tau: u32,
    pub idx: u32,
}

/// 整数二项式 C(n, k)。
pub fn binom(n: u32, k: u32) -> Q {
    if k > n {
        return Q::ZERO;
    }
    let mut v = Q::ONE;
    for i in 0..k {
        v = v * Q::from((n - i) as i64) / Q::from((k - i) as i64);
    }
    v
}

fn q_pow_neg1(e: i64) -> Q {
    Q::from(if e % 2 == 0 { 1 } else { -1 })
}

/// 慢参考版：Σ_{σ∈S_m} (−1)^{m(m−1)/2}·parity(σ)·Π_i [binom(x_i+y_{σi},x_i)·cup]·
/// gwpoint([x_i+y_{σi}−1] ++ tau_point)。
/// 环路指标越界（> g）或 a/b 未配对 → 不变量为 0（外代数消没）。
pub fn gwloop(g: u32, d: u32, tau_point: &[u32], ta: &[LoopIns], tb: &[LoopIns]) -> Q {
    let m = ta.len();
    if m != tb.len() {
        return Q::ZERO;
    }
    if ta.iter().any(|t| t.idx == 0 || t.idx > g) || tb.iter().any(|t| t.idx == 0 || t.idx > g) {
        return Q::ZERO;
    }
    if m == 0 {
        return gwpoint(g, d, tau_point);
    }
    let mut answer = Q::ZERO;
    for sigma in permute(m as u32) {
        let mut prod = Q::ONE;
        let mut new_tau: Vec<u32> = Vec::with_capacity(m + tau_point.len());
        for (i1, si) in sigma.iter().enumerate() {
            let j = (*si - 1) as usize;
            let x = ta[i1].tau;
            let y = tb[j].tau;
            let cup = if ta[i1].idx == tb[j].idx { 1 } else { 0 };
            prod = prod * binom(x + y, x) * Q::from(cup);
            if x + y == 0 {
                // 后代数 −1：负指标 → 不变量为零（GW 层的 member(−1) 检查下沉）
                prod = Q::ZERO;
                break;
            }
            new_tau.push(x + y - 1);
        }
        if prod == Q::ZERO {
            continue;
        }
        new_tau.extend_from_slice(tau_point);
        let sig = q_pow_neg1(m as i64 * (m as i64 - 1) / 2) * Q::from(parity(&sigma));
        answer += prod * sig * gwpoint(g, d, &new_tau);
    }
    answer
}

/// 快版：按环路指标分组，组内枚举置换（组间独立）；组内后代重复 → 0。
pub fn fast_gwloop(g: u32, d: u32, tau_point: &[u32], ta: &[LoopIns], tb: &[LoopIns]) -> Q {
    let m = ta.len();
    if m != tb.len() {
        return Q::ZERO;
    }
    if ta.iter().any(|t| t.idx == 0 || t.idx > g) || tb.iter().any(|t| t.idx == 0 || t.idx > g) {
        return Q::ZERO;
    }
    if m == 0 {
        return gwpoint(g, d, tau_point);
    }
    // 分组（头插语义不影响结果：GWpoint 对 τ 乘性对称）
    let mut alist: Vec<Vec<u32>> = vec![Vec::new(); g as usize + 1];
    let mut blist: Vec<Vec<u32>> = vec![Vec::new(); g as usize + 1];
    let mut reordersign: i64 = 0; // 累计指数（模 2）
    for i in 0..m {
        let (ca, cb) = (ta[i].idx as usize, tb[i].idx as usize);
        let mut expo = 0i64;
        expo += alist[ca..].iter().map(|v| v.len() as i64).sum::<i64>();
        expo += blist[cb..].iter().map(|v| v.len() as i64).sum::<i64>();
        reordersign += expo;
        alist[ca].push(ta[i].tau);
        blist[cb].push(tb[i].tau);
    }
    // a/b 分组大小不对称 → cup 全零 → 不变量为零
    for i in 1..=g as usize {
        if alist[i].len() != blist[i].len() {
            return Q::ZERO;
        }
    }
    // 组内重复后代 → 外代数消没
    for i in 1..=g as usize {
        let mut sa = alist[i].clone();
        sa.sort_unstable();
        sa.dedup();
        if sa.len() != alist[i].len() {
            return Q::ZERO;
        }
        let mut sb = blist[i].clone();
        sb.sort_unstable();
        sb.dedup();
        if sb.len() != blist[i].len() {
            return Q::ZERO;
        }
    }
    // 各组置换的笛卡尔积
    let perms: Vec<Vec<Vec<u32>>> = (1..=g as usize)
        .map(|i| permute(alist[i].len() as u32))
        .collect();
    let mut answer = Q::ZERO;
    for sig in CartProd::new(&perms) {
        let mut evil_par: i64 = 1; // Π 组奇偶
        let mut new_tau: Vec<u32> = tau_point.to_vec();
        let mut bin = Q::ONE;
        for i in 1..=g as usize {
            let p = &perms[i - 1][sig[i - 1]];
            if parity(p) < 0 {
                evil_par = -evil_par;
            }
            for k in 0..alist[i].len() {
                let x = alist[i][k];
                let y = blist[i][p[k] as usize - 1];
                bin *= binom(x + y, x);
                if x + y == 0 {
                    bin = Q::ZERO;
                    break;
                }
                new_tau.push(x + y - 1);
            }
        }
        if bin != Q::ZERO {
            answer += bin * Q::from(evil_par) * gwpoint(g, d, &new_tau);
        }
    }
    answer * q_pow_neg1(m as i64 * (m as i64 - 1) / 2 + reordersign)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn li(tau: u32, idx: u32) -> LoopIns {
        LoopIns { tau, idx }
    }

    /// 双实现互验：GWloop == fastGWloop（系统枚举）。
    #[test]
    fn gwloop_vs_fast() {
        let mut count = 0;
        for g in 1..=3u32 {
            for d in 1..=2u32 {
                for ta in [
                    vec![],
                    vec![li(0, 1)],
                    vec![li(1, 1)],
                    vec![li(0, 1), li(0, 1)],
                    vec![li(0, 1), li(1, 1)],
                    vec![li(0, 1), li(0, 2)],
                    vec![li(2, 1), li(0, 2)],
                ] {
                    for tb in [
                        vec![li(0, 1)],
                        vec![li(1, 1)],
                        vec![li(0, 2)],
                        vec![li(0, 1), li(0, 1)],
                        vec![li(1, 1), li(0, 1)],
                        vec![li(0, 2), li(1, 1)],
                        vec![li(0, 2), li(0, 2)],
                    ] {
                        if ta.len() != tb.len() {
                            continue;
                        }
                        for tp in [vec![], vec![0u32], vec![0u32, 1u32]] {
                            let a = gwloop(g, d, &tp, &ta, &tb);
                            let b = fast_gwloop(g, d, &tp, &ta, &tb);
                            assert_eq!(
                                a.to_string(),
                                b.to_string(),
                                "g={g} d={d} tp={tp:?} ta={ta:?} tb={tb:?}"
                            );
                            count += 1;
                        }
                    }
                }
            }
        }
        assert!(count > 200, "应覆盖 {count} 组");
    }

    /// 外代数消没：同组重复后代 / 未配对环路 → 双侧皆 0。
    #[test]
    fn antisym_vanishing() {
        // 同 (idx=1) 组重复 tau=0：fast 直接剪枝，GWloop 置换和消没
        let ta = vec![li(0, 1), li(0, 1)];
        let tb = vec![li(0, 1), li(1, 1)];
        assert_eq!(gwloop(2, 2, &[0], &ta, &tb).to_string(), "0");
        assert_eq!(fast_gwloop(2, 2, &[0], &ta, &tb).to_string(), "0");
        // a、b 环路指标集不同 → cup 全零
        let ta2 = vec![li(0, 1)];
        let tb2 = vec![li(0, 2)];
        assert_eq!(gwloop(2, 1, &[], &ta2, &tb2).to_string(), "0");
        assert_eq!(fast_gwloop(2, 1, &[], &ta2, &tb2).to_string(), "0");
        // 指标越界（> g）
        assert_eq!(gwloop(1, 1, &[], &[li(0, 2)], &[li(0, 2)]).to_string(), "0");
    }

    /// m=1 退化：binom(x+y,x)·gwpoint([x+y−1]++pt)，符号恒正。
    #[test]
    fn m1_degenerate() {
        let ta = vec![li(1, 2)];
        let tb = vec![li(0, 2)];
        let direct = binom(1, 1) * gwpoint(3, 1, &[1 - 1, 0]);
        assert_eq!(gwloop(3, 1, &[0], &ta, &tb).to_string(), direct.to_string());
        assert_eq!(
            fast_gwloop(3, 1, &[0], &ta, &tb).to_string(),
            direct.to_string()
        );
    }
}
