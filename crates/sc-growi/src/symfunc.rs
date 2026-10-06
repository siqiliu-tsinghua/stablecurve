//! Virasoro 条件的对称函数系数——原版 `libgw/symfunc.{h,cc}` 的直译。
//!
//! `symfunc(x,k,i)` = 数 x,x+1,…,x+k 的第 (k+1−i) 个初等对称函数，即
//! `∏(s+x+j) = Σ symfunc(x,k,i)·s^i`。递推：
//! `symfunc(x,k,i) = Σ_j (x+j)·symfunc(x+j+1, k−1−j, i−j)`。
//! 原版存 symfunc.dbm；这里用进程内缓存（同查询跨 Virasoro 调用复用）。

use std::cell::RefCell;
use std::collections::HashMap;

use crate::num::{One, Q, Zero};

thread_local! {
    static CACHE: RefCell<HashMap<(String, i32, i32), Q>> = RefCell::new(HashMap::new());
}

pub fn symfunc(x: &Q, k: i32, i: i32) -> Q {
    // 初等情形
    if i < 0 || k < -1 || i > k + 1 {
        return Q::ZERO;
    }
    if i == k + 1 {
        return Q::ONE;
    }

    let key = (x.to_string(), k, i);
    let hit = CACHE.with(|c| c.borrow().get(&key).cloned());
    if let Some(v) = hit {
        return v;
    }

    // 递推计算
    let mut result = Q::ZERO;
    for j in 0..=i {
        let xj = x + Q::from(j);
        let xj1 = xj.clone() + Q::ONE;
        result += &xj * &symfunc(&xj1, k - 1 - j, i - j);
    }

    CACHE.with(|c| c.borrow_mut().insert(key, result.clone()));
    result
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
