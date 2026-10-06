//! GWloop / fastGWloop：环路后代（H¹ 外代数）的消去层（M2）。
//!
//! 数学：a₁..a_g, b₁..b_g 为 H¹(C_g) 的辛基（a_i∪b_j = δ_ij，a∪a = b∪b = 0）。
//! 环路后代反对易 → 不变量 = 置换和 + 反对称符号：
//!   ⟨Πτ(pt)·Π τ_{x_i}(a_{α_i})·τ_{y_i}(b_{β_i})⟩ = Σ_σ sign·Π binom(x+y_σ,x)·cup·⟨纯点⟩
//! 两实现互为内建差分（原版即双实现）：
//!
//! - [`gwloop`]：全 S_m 置换和（慢参考版，GWall `GWloop`）
//! - [`fast_gwloop`]：行列式化——对每个分拆 λ，置换和是一个 m×m 行列式
//!   （GWall `fastGWloop` 按环路指标分组枚举组内置换，仍是阶乘级）
//!
//! 符号统一用 (−1)^{m(m−1)/2}·parity(σ)（原版 evils 的魔数依赖 Maple
//! permute 枚举顺序且作者自注仅 m<9 有效——弃用，见 comb::parity）。

use crate::comb::{parity, permute};
use crate::gwpoint::{gwpoint, partition_weights, tau_factor};
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

/// 快版：行列式化。gwpoint 对每个分拆 λ 关于 τ 是乘性的，
///   ⟨Π τ_k(pt)⟩ = Σ_λ w(λ)·Π_k f_λ(k)，
/// 故置换和 Σ_σ sgn σ·Π_i binom·cup·⟨…⟩ = Σ_λ w(λ)·Π_{τ_point} f_λ·det M_λ，
///   M_λ[i][j] = C(x_i+y_j, x_i)·[α_i = β_j]·f_λ(x_i+y_j−1)（x_i+y_j = 0 时为 0）。
/// 行（a 侧）、列（b 侧）各按环路指标稳定排序后 M_λ 为块对角，det 等于两个
/// 重排的符号乘各块行列式之积。大小为 1 的块直接取元素（常见情形，与按组
/// 枚举置换同价），同指标多次出现的块做 O(k³) 消元而非 k! 项置换和。
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
    let mut rows: Vec<usize> = (0..m).collect();
    rows.sort_by_key(|&i| ta[i].idx);
    let mut cols: Vec<usize> = (0..m).collect();
    cols.sort_by_key(|&j| tb[j].idx);
    // a/b 分组大小不对称 → cup 矩阵奇异 → 不变量为零
    if rows
        .iter()
        .map(|&i| ta[i].idx)
        .ne(cols.iter().map(|&j| tb[j].idx))
    {
        return Q::ZERO;
    }
    // 同组块（排序后位置区间）；组内后代重复 → 两行（列）相同 → 外代数消没
    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let mut s = 0;
    while s < m {
        let mut e = s + 1;
        while e < m && ta[rows[e]].idx == ta[rows[s]].idx {
            e += 1;
        }
        let distinct = |taus: Vec<u32>| {
            let mut t = taus;
            t.sort_unstable();
            t.windows(2).all(|w| w[0] != w[1])
        };
        if !distinct(rows[s..e].iter().map(|&i| ta[i].tau).collect())
            || !distinct(cols[s..e].iter().map(|&j| tb[j].tau).collect())
        {
            return Q::ZERO;
        }
        blocks.push((s, e));
        s = e;
    }
    let sign = perm_sign(&rows) * perm_sign(&cols);
    let entry = |part: &[u32], i: usize, j: usize| -> Q {
        let (x, y) = (ta[i].tau, tb[j].tau);
        if x + y == 0 {
            Q::ZERO // 后代数 −1：负指标 → 零
        } else {
            binom(x + y, x) * tau_factor(part, x + y - 1)
        }
    };
    let mut answer = Q::ZERO;
    'part: for (w, part) in partition_weights(g, d) {
        let mut prod = w;
        for &k in tau_point {
            prod *= tau_factor(&part, k);
        }
        for &(s, e) in &blocks {
            let dv = if e - s == 1 {
                entry(&part, rows[s], cols[s])
            } else {
                det(rows[s..e]
                    .iter()
                    .map(|&i| cols[s..e].iter().map(|&j| entry(&part, i, j)).collect())
                    .collect())
            };
            if dv == Q::ZERO {
                continue 'part;
            }
            prod *= dv;
        }
        answer += prod;
    }
    answer * Q::from(sign) * q_pow_neg1(m as i64 * (m as i64 - 1) / 2)
}

/// 排列 p（p[k] = 位置 k 上的原下标）的符号。
fn perm_sign(p: &[usize]) -> i64 {
    let mut inv = 0usize;
    for i in 0..p.len() {
        for j in i + 1..p.len() {
            if p[i] > p[j] {
                inv += 1;
            }
        }
    }
    if inv % 2 == 0 { 1 } else { -1 }
}

/// 有理方阵的行列式（高斯消元，选非零主元）。
fn det(mut a: Vec<Vec<Q>>) -> Q {
    let n = a.len();
    let mut det = Q::ONE;
    for c in 0..n {
        let Some(p) = (c..n).find(|&r| a[r][c] != Q::ZERO) else {
            return Q::ZERO;
        };
        if p != c {
            a.swap(p, c);
            det = -det;
        }
        let (top, rest) = a.split_at_mut(c + 1);
        let pivot = &top[c];
        det *= &pivot[c];
        for row in rest.iter_mut() {
            if row[c] == Q::ZERO {
                continue;
            }
            let f = &row[c] / &pivot[c];
            for k in c..n {
                row[k] -= &f * &pivot[k];
            }
        }
    }
    det
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

    /// 行列式版对照慢参考版：m = 3、4，指标与后代次数混合（含同组重复后代、
    /// a/b 分组不对称、x+y = 0 的负后代）。
    #[test]
    fn determinant_vs_permutation_sum() {
        let mut seed = 0x2545f4914f6cdd1du64;
        let mut rnd = |k: u32| -> u32 {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed % k as u64) as u32
        };
        let mut nonzero = 0;
        for _ in 0..300 {
            let g = 2 + rnd(3);
            let d = 1 + rnd(3);
            let m = 3 + rnd(2) as usize;
            let ta: Vec<LoopIns> = (0..m).map(|_| li(rnd(6), 1 + rnd(g))).collect();
            // b 侧指标多数取 a 侧指标的随机排列（分组对上才可能非零），少数随机
            let mut idx: Vec<u32> = ta.iter().map(|t| t.idx).collect();
            for i in (1..m).rev() {
                idx.swap(i, rnd(i as u32 + 1) as usize);
            }
            if rnd(4) == 0 {
                idx[0] = 1 + rnd(g);
            }
            let tb: Vec<LoopIns> = idx.iter().map(|&i| li(rnd(6), i)).collect();
            let tp: Vec<u32> = (0..rnd(3)).map(|_| rnd(3)).collect();
            let a = gwloop(g, d, &tp, &ta, &tb);
            let b = fast_gwloop(g, d, &tp, &ta, &tb);
            assert_eq!(a, b, "g={g} d={d} tp={tp:?} ta={ta:?} tb={tb:?}");
            if a != Q::ZERO {
                nonzero += 1;
            }
        }
        assert!(nonzero > 40, "非零样本过少：{nonzero}");
    }
}
