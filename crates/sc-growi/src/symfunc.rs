//! Virasoro 条件的对称函数系数；逐次乘线性因子，保存完整系数行。
//!
//! `symfunc(x,k,i)` = 数 x,x+1,…,x+k 的第 (k+1−i) 个初等对称函数，即
//! `∏(s+x+j) = Σ symfunc(x,k,i)·s^i`。递推：
//! `symfunc(x,k,i) = Σ_j (x+j)·symfunc(x+j+1, k−1−j, i−j)`。
//! 原版存 symfunc.dbm；这里直接以有理数为键复用系数表，省去字符串转换。

use std::cell::RefCell;
use std::collections::BTreeMap;

use crate::num::{One, Q, Zero};

thread_local! {
    static CACHE: RefCell<BTreeMap<Q, Vec<Vec<Q>>>> = const { RefCell::new(BTreeMap::new()) };
}

pub fn symfunc(x: &Q, k: i32, i: i32) -> Q {
    // 初等情形
    if i < 0 || k < -1 || i > k + 1 {
        return Q::ZERO;
    }
    if i == k + 1 {
        return Q::ONE;
    }

    CACHE.with(|cell| {
        let mut cache = cell.borrow_mut();
        let rows = cache.entry(x.clone()).or_insert_with(|| vec![vec![Q::ONE]]);
        // 第 r 行是 ∏_{j<r}(s+x+j)；一次延长得到该行所有系数。
        while rows.len() <= k as usize + 1 {
            let r = rows.len() - 1;
            let shift = x + Q::from(r);
            let previous = &rows[r];
            let mut next = vec![Q::ZERO; previous.len() + 1];
            for (j, c) in previous.iter().enumerate() {
                next[j] += c * &shift;
                next[j + 1] += c;
            }
            rows.push(next);
        }
        rows[k as usize + 1][i as usize].clone()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::{q, q_frac};

    #[test]
    fn base_cases() {
        assert_eq!(symfunc(&q(1), -1, 0), Q::ONE);
        assert_eq!(symfunc(&q(1), 0, 0), q(1));
        assert_eq!(symfunc(&q(1), 0, 1), Q::ONE);
        assert_eq!(symfunc(&q(1), 5, 7), Q::ZERO); // i > k+1
    }

    #[test]
    fn product_identity() {
        // ∏_{j=0..k}(s+x+j) = Σ_i symfunc(x,k,i) s^i，取 x=1/2, k=2 验证 s=1
        let x = q_frac(1, 2);
        let k = 2i32;
        let mut lhs = Q::ONE;
        for j in 0..=k {
            lhs *= q_frac((3 + 2 * j) as i64, 2); // s=1 + x + j
        }
        let mut rhs = Q::ZERO;
        for i in 0..=k + 1 {
            rhs += symfunc(&x, k, i);
        }
        assert_eq!(lhs, rhs);
    }
}
