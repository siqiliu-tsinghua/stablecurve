//! 组合件：Maple `combinat` 库中被 GWall 实际用到的三个函数的自实现。
//!
//! `permute`（全排列）、`cartprod`（笛卡尔积，逐组取一的枚举器）、
//! `partition`（整数分拆，**升序**——GWall 的 `power()` 依赖此约定做负下标）。
//! 另附置换奇偶 `parity`（Maple `parity(convert(σ,'disjcyc'))` 的等价实现，
//! 逆序数计数——原版 `evils` 的魔数公式依赖 Maple 枚举顺序，弃用，见 REDESIGN）。

/// {1..m} 的全部排列（顺序不保证与 Maple 一致——所有调用方按 σ 的内容
/// 计算，不依赖枚举顺序；唯一例外是被弃用的 evils 魔数）。
pub fn permute(m: u32) -> Vec<Vec<u32>> {
    let mut out = Vec::new();
    let mut cur: Vec<u32> = (1..=m).collect();
    perm_rec(&mut cur, 0, &mut out);
    out
}

fn perm_rec(cur: &mut Vec<u32>, i: usize, out: &mut Vec<Vec<u32>>) {
    if i == cur.len() {
        out.push(cur.clone());
        return;
    }
    for j in i..cur.len() {
        cur.swap(i, j);
        perm_rec(cur, i + 1, out);
        cur.swap(i, j);
    }
}

/// 置换奇偶性（1 = 偶，−1 = 奇）：逆序数计数，等价于
/// Maple 的 parity(convert(σ, 'disjcyc'))。
pub fn parity(sigma: &[u32]) -> i64 {
    let mut inv = 0u64;
    for i in 0..sigma.len() {
        for j in (i + 1)..sigma.len() {
            if sigma[i] > sigma[j] {
                inv += 1;
            }
        }
    }
    if inv % 2 == 0 { 1 } else { -1 }
}

/// 笛卡尔积枚举器（Maple `cartprod`）：每组一个取值的所有组合，
/// 按"后组先动"的字典序（与 Maple 一致；调用方求和，顺序不敏感）。
pub struct CartProd {
    sizes: Vec<usize>,
    cur: Vec<usize>,
    finished: bool,
}

impl CartProd {
    /// groups：每组一个取值集合（任意元素类型，只关心长度）。
    pub fn new<T>(groups: &[Vec<T>]) -> Self {
        let sizes: Vec<usize> = groups.iter().map(|v| v.len()).collect();
        let finished = sizes.contains(&0);
        let cur = vec![0usize; sizes.len()];
        CartProd {
            sizes,
            cur,
            finished,
        }
    }
}

impl Iterator for CartProd {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Vec<usize>> {
        if self.finished {
            return None;
        }
        let out = self.cur.clone();
        // 最后一位 +1 进位
        let mut i = self.cur.len();
        loop {
            if i == 0 {
                self.finished = true;
                break;
            }
            i -= 1;
            self.cur[i] += 1;
            if self.cur[i] < self.sizes[i] {
                break;
            }
            self.cur[i] = 0;
        }
        Some(out)
    }
}

/// 整数 d 的全部分拆，**每项升序**（与 Maple partition 一致），
/// 分拆之间按 GWall 所需即可（调用方只做无序求和）。
pub fn partition(d: u32) -> Vec<Vec<u32>> {
    let mut out = Vec::new();
    part_rec(d, 1, &mut Vec::new(), &mut out);
    out
}

fn part_rec(rem: u32, min: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
    if rem == 0 {
        out.push(cur.clone());
        return;
    }
    // 升序：下一个分量 ≥ min 且 ≤ rem
    for x in min..=rem {
        cur.push(x);
        part_rec(rem - x, x, cur, out);
        cur.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permute_basic() {
        assert_eq!(permute(0), vec![Vec::<u32>::new()]);
        assert_eq!(permute(1), vec![vec![1]]);
        let p3 = permute(3);
        assert_eq!(p3.len(), 6);
        // 恰为 S_3 全体（互异）
        let mut s = p3.clone();
        s.sort();
        s.dedup();
        assert_eq!(s.len(), 6);
        assert!(p3.contains(&vec![1, 2, 3]));
        assert!(p3.contains(&vec![3, 2, 1]));
    }

    #[test]
    fn parity_basics() {
        assert_eq!(parity(&[1, 2, 3]), 1);
        assert_eq!(parity(&[2, 1, 3]), -1);
        assert_eq!(parity(&[2, 3, 1]), 1); // 3-循环 = 偶
        assert_eq!(parity(&[3, 2, 1]), -1); // 3 逆序
        // 与 disjcyc 公式一致：parity = (−1)^{m−#cycles}
        let m = 4u32;
        let p = vec![2, 1, 4, 3u32]; // 两个 2-循环 → 偶
        let mut seen = [false; 4];
        let mut cycles = 0;
        for i in 0..4u32 {
            if !seen[i as usize] {
                cycles += 1;
                let mut j = i;
                while !seen[j as usize] {
                    seen[j as usize] = true;
                    j = p[j as usize] - 1;
                }
            }
        }
        assert_eq!(
            parity(&p),
            if (m as i64 - cycles) % 2 == 0 { 1 } else { -1 }
        );
    }

    #[test]
    fn cartprod_basic() {
        let groups = vec![vec![10, 11], vec![20, 21, 22]];
        let got: Vec<Vec<usize>> = CartProd::new(&groups).collect();
        assert_eq!(got.len(), 6);
        assert!(got.contains(&vec![0, 0]));
        assert!(got.contains(&vec![1, 2]));
    }

    #[test]
    fn partition_ascending() {
        let ps = partition(5);
        assert_eq!(ps.len(), 7);
        for p in &ps {
            // 升序验证
            assert!(p.windows(2).all(|w| w[0] <= w[1]), "{p:?}");
            assert_eq!(p.iter().sum::<u32>(), 5);
        }
        assert!(ps.contains(&vec![1, 1, 1, 1, 1]));
        assert!(ps.contains(&vec![5]));
        assert!(ps.contains(&vec![2, 3]));
        assert_eq!(partition(0), vec![Vec::<u32>::new()]);
    }
}
