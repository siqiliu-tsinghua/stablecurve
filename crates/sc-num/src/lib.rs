//! 有理数后端。
//!
//! 默认使用 malachite（纯 Rust）；若日后 profile 显示大数运算是瓶颈，可在
//! `rug` feature 下切换到 GMP 绑定（REDESIGN §2.1）。全部代码只通过 [`Q`]
//! 别名和本文件的辅助函数接触有理数，切换后端时不动算法层。
//!
//! 与原版 `math/rational.cc` 的语义差异：原版用空指针表示 UNDEFINED 并让
//! 运算传染；本 crate 用 `Option`/`Result` 显式传播（见各算法模块）。

pub use malachite_q::Rational as Q;

// ZERO/ONE 常量经由 traits 提供，重导出给本 crate 各模块使用
pub use malachite_base::num::basic::traits::{One, Zero};

/// 从整数构造。
#[inline]
pub fn q(n: i64) -> Q {
    Q::from(n)
}

/// 是否为零（未定义视为非零，与原版 is_zero 一致）。
#[inline]
pub fn is_zero(x: &Q) -> bool {
    *x == Q::ZERO
}

/// 从分数 `n/d` 构造（自动规范化符号与既约性）。
#[inline]
pub fn q_frac(n: i64, d: i64) -> Q {
    Q::from_signeds(n, d)
}

/// 精确阶乘表，按所需最大指标逐步延长。
pub fn factorial(n: u32) -> Q {
    thread_local! {
        static TABLE: std::cell::RefCell<Vec<Q>> = std::cell::RefCell::new(vec![Q::ONE]);
    }
    TABLE.with(|cell| {
        let mut table = cell.borrow_mut();
        while table.len() <= n as usize {
            let next = table.last().unwrap() * Q::from(table.len());
            table.push(next);
        }
        table[n as usize].clone()
    })
}

/// 精确二项式；中间乘积同样使用任意精度。
pub fn binomial(n: u32, k: u32) -> Q {
    if k > n {
        return Q::ZERO;
    }
    (0..k.min(n - k)).fold(Q::ONE, |v, j| v * Q::from(n - j) / Q::from(j + 1))
}

/// Bernoulli 数，B₁ = −1/2。每行的二项系数用相邻比递推。
pub fn bernoulli(n: u64) -> Q {
    thread_local! {
        static TABLE: std::cell::RefCell<Vec<Q>> = std::cell::RefCell::new(vec![Q::ONE]);
    }
    let n = usize::try_from(n).expect("Bernoulli 指标超出 usize");
    TABLE.with(|cell| {
        let mut b = cell.borrow_mut();
        while b.len() <= n {
            let m = b.len();
            let mut coefficient = Q::ONE;
            let mut sum = Q::ZERO;
            for (k, bk) in b.iter().enumerate() {
                sum += &coefficient * bk;
                coefficient = coefficient * Q::from(m + 1 - k) / Q::from(k + 1);
            }
            b.push(-sum / Q::from(m + 1));
        }
        b[n].clone()
    })
}

/// 平方求幂，不缩窄指数。
pub fn pow_unsigned(base: &Q, mut exponent: u64) -> Q {
    let mut base = base.clone();
    let mut value = Q::ONE;
    while exponent != 0 {
        if exponent & 1 != 0 {
            value *= &base;
        }
        exponent >>= 1;
        if exponent != 0 {
            base = &base * &base;
        }
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        let x = q_frac(1, 2);
        let y = q_frac(1, 3);
        assert_eq!(&x + &y, q_frac(5, 6));
        assert_eq!(&x * &y, q_frac(1, 6));
        assert_eq!(&x / &y, q_frac(3, 2));
        assert_eq!(-&x, q_frac(-1, 2));
    }

    #[test]
    fn normalization() {
        // 符号归一到分母、自动约分
        assert_eq!(q_frac(2, -4), q_frac(-1, 2));
        assert_eq!(q_frac(-6, -3), q(2));
    }

    #[test]
    fn big_values() {
        // 大数必须精确（gwcheck 里最大的 Kontsevich 数有 39 位）
        let mut x = q(1);
        for _ in 0..128 {
            x = &x * &q(2);
        }
        assert_eq!(x.to_string(), "340282366920938463463374607431768211456");
    }
}
