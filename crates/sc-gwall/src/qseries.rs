//! 截断有理系数幂级数（Z/Y 配分函数组装用的最小级数机器）。
//!
//! `QSeries { c: Vec<Q> }`：x⁰..x^N 的系数，长度即截断阶 + 1。
//! 提供 add / scale / mul（截断）/ exp（ODE 递推）/ 代换 x ← f（内层无常数项）。
//! Bernoulli 数（Maple `bernoulli`；B₁ = −1/2 约定，与 sc-faber socle 一致）。

use sc_num::{One, Q, Zero};

/// Bernoulli 数 B_n（B₁ = −1/2）：递推 Σ_{k} C(n+1,k) B_k = 0。
pub fn bernoulli(n: u32) -> Q {
    let mut b: Vec<Q> = vec![Q::ONE]; // B_0 = 1
    for m in 1..=n as usize {
        // B_m 由 Σ_{k=0..m} C(m+1,k) B_k = 0 解出
        let mut acc = Q::ZERO;
        for (k, bk) in b.iter().enumerate() {
            acc += binom_q(m + 1, k) * bk;
        }
        b.push(Q::ZERO - acc / binom_q(m + 1, m));
    }
    b[n as usize].clone()
}

fn binom_q(a: usize, b: usize) -> Q {
    if b > a {
        return Q::ZERO;
    }
    let mut v = Q::ONE;
    for i in 0..b {
        v = v * Q::from((a - i) as i64) / Q::from((b - i) as i64);
    }
    v
}

fn fact(n: u32) -> Q {
    let mut v = Q::ONE;
    for i in 2..=n {
        v *= Q::from(i as i64);
    }
    v
}

#[derive(Clone, Debug, PartialEq)]
pub struct QSeries {
    pub c: Vec<Q>, // c[i] = x^i 系数；末尾零系数可有可无
}

impl QSeries {
    /// 截断到阶 N（丢弃高次）。
    pub fn trunc(&self, n: u32) -> Self {
        QSeries {
            c: self.c.iter().take(n as usize + 1).cloned().collect(),
        }
    }

    pub fn constant(v: Q, n: u32) -> Self {
        let mut c = vec![Q::ZERO; n as usize + 1];
        c[0] = v;
        QSeries { c }
    }

    /// 变量 x 本身。
    pub fn var(n: u32) -> Self {
        let mut c = vec![Q::ZERO; n as usize + 1];
        if n >= 1 {
            c[1] = Q::ONE;
        }
        QSeries { c }
    }

    pub fn add(&self, o: &Self) -> Self {
        let n = self.c.len().max(o.c.len());
        let c = (0..n)
            .map(|i| {
                let a = self.c.get(i).cloned().unwrap_or(Q::ZERO);
                let b = o.c.get(i).cloned().unwrap_or(Q::ZERO);
                a + b
            })
            .collect();
        QSeries { c }
    }

    pub fn scale(&self, k: &Q) -> Self {
        QSeries {
            c: self.c.iter().map(|x| x.clone() * k.clone()).collect(),
        }
    }

    pub fn mul(&self, o: &Self) -> Self {
        // 截断级数语义：结果截到两者较短长度
        let n = self.c.len().min(o.c.len());
        let mut c = vec![Q::ZERO; n];
        for (i, a) in self.c.iter().enumerate() {
            if *a == Q::ZERO || i >= n {
                continue;
            }
            for (j, b) in o.c.iter().enumerate() {
                if i + j >= n {
                    break;
                }
                c[i + j] = c[i + j].clone() + a.clone() * b.clone();
            }
        }
        QSeries { c }
    }

    /// exp(this)（截断到原阶数）：S' = S'·f 的 ODE 递推
    /// E = Σ e_n x^n，e_0 = exp(f_0)（要求 f_0 = 0，否则只支持 f_0 使
    /// exp(f_0) 为有理数——GWall 的 Z 中 N₀(0) = 0，恒满足 f_0 = 0）。
    pub fn exp(&self) -> Result<Self, String> {
        let f0 = self.c.first().cloned().unwrap_or(Q::ZERO);
        if f0 != Q::ZERO {
            return Err("exp 要求常数项为 0".into());
        }
        let n = self.c.len();
        let mut e = vec![Q::ONE]; // e_0 = 1
        for m in 1..n {
            // m·e_m = Σ_{k=1..m} k·f_k·e_{m−k}
            let mut acc = Q::ZERO;
            for k in 1..=m {
                let fk = self.c.get(k).cloned().unwrap_or(Q::ZERO);
                let emk = e.get(m - k).cloned().unwrap_or(Q::ZERO);
                acc += Q::from(k as i64) * fk * emk;
            }
            e.push(acc / Q::from(m as i64));
        }
        Ok(QSeries { c: e })
    }

    /// 乘法逆元 1/this（截断到自身阶数）：要求常数项非零。
    /// 递推 b_0 = 1/a_0，b_m = −(Σ_{k=1..m} a_k·b_{m−k})/a_0。
    pub fn inverse(&self) -> Result<Self, String> {
        let a0 = self.c.first().cloned().unwrap_or(Q::ZERO);
        if a0 == Q::ZERO {
            return Err("inverse 要求常数项非零".into());
        }
        let n = self.c.len();
        let mut b: Vec<Q> = Vec::with_capacity(n);
        b.push(Q::ONE / &a0);
        for m in 1..n {
            let mut acc = Q::ZERO;
            for k in 1..=m {
                acc += &self.c[k] * &b[m - k];
            }
            b.push(-acc / &a0);
        }
        Ok(QSeries { c: b })
    }

    /// 整数幂（可为负；负幂要求常数项非零）。
    pub fn pow_signed(&self, e: i64) -> Result<Self, String> {
        if e >= 0 {
            Ok(self.powi(e as u32))
        } else {
            Ok(self.inverse()?.powi((-e) as u32))
        }
    }

    /// 代换 x ← f(y)（f 常数项须为 0；结果截断到自身阶数）。
    pub fn compose(&self, f: &Self) -> Result<Self, String> {
        let f0 = f.c.first().cloned().unwrap_or(Q::ZERO);
        if f0 != Q::ZERO {
            return Err("compose 要求内层级数无常数项".into());
        }
        let n = self.c.len();
        if n == 0 {
            return Ok(QSeries { c: Vec::new() });
        }
        // Horner：S = c_0 + x(c_1 + x(c_2 + ...))，中间保持目标截断长度
        let mut acc = QSeries::constant(self.c.last().cloned().unwrap_or(Q::ZERO), n as u32 - 1);
        for i in (0..n - 1).rev() {
            acc = QSeries::constant(self.c[i].clone(), n as u32 - 1).add(&acc.mul(f));
        }
        Ok(acc)
    }

    /// 幂 series^k（k ≥ 0 整数幂，截断长度保持与自身一致）。
    pub fn powi(&self, k: u32) -> Self {
        if self.c.is_empty() {
            return QSeries { c: Vec::new() };
        }
        let mut acc = QSeries::constant(Q::ONE, self.c.len() as u32 - 1);
        for _ in 0..k {
            acc = acc.mul(self);
        }
        acc
    }

    /// 数值求值（对给定有理 x，截断内的多项式值）——测试用。
    pub fn eval(&self, x: &Q) -> Q {
        let mut v = Q::ZERO;
        for c in self.c.iter().rev() {
            v = v * x + c;
        }
        v
    }
}

/// z/(2·sinh(z/2)) 的 x^l 系数（GWall `power()` 的 coeftayl）。
/// 恒等式：t·e^{t/2}/(e^t−1) = t/(2·sinh(t/2)) 是 Bernoulli 多项式的
/// 母函数取 x = 1/2，故系数 = B_l(1/2)/l!（奇 l 自动为 0——中心对称）。
pub fn z_over_2sinh_coeff(l: u32) -> Q {
    // B_l(1/2) = Σ_k C(l,k) B_k 2^{k−l}
    let mut acc = Q::ZERO;
    for k in 0..=l {
        let mut pow2 = Q::ONE;
        for _ in 0..(l - k) {
            pow2 *= Q::from(2);
        }
        acc += binom_q(l as usize, k as usize) * bernoulli(k) / pow2;
    }
    acc / fact(l)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(v: i64) -> Q {
        Q::from(v)
    }

    #[test]
    fn bernoulli_vals() {
        // 与 Maple/sage 一致（B₁ = −1/2）
        assert_eq!(bernoulli(0).to_string(), "1");
        assert_eq!(bernoulli(1).to_string(), "-1/2");
        assert_eq!(bernoulli(2).to_string(), "1/6");
        assert_eq!(bernoulli(4).to_string(), "-1/30");
        assert_eq!(bernoulli(6).to_string(), "1/42");
        assert_eq!(bernoulli(10).to_string(), "5/66");
        assert_eq!(bernoulli(12).to_string(), "-691/2730");
    }

    #[test]
    fn sinh_coeff_basics() {
        // z/(2 sinh(z/2)) = 1 − z²/24 + 7z⁴/5760 − ...
        assert_eq!(z_over_2sinh_coeff(0).to_string(), "1");
        assert_eq!(z_over_2sinh_coeff(1).to_string(), "0");
        assert_eq!(z_over_2sinh_coeff(2).to_string(), "-1/24");
        assert_eq!(z_over_2sinh_coeff(4).to_string(), "7/5760");
    }

    #[test]
    fn series_exp() {
        // exp(x) 截 4：1 + x + x²/2 + x³/6 + x⁴/24
        let x = QSeries::var(4);
        let e = x.exp().unwrap();
        assert_eq!(e.c[0].to_string(), "1");
        assert_eq!(e.c[1].to_string(), "1");
        assert_eq!(e.c[2].to_string(), "1/2");
        assert_eq!(e.c[3].to_string(), "1/6");
        assert_eq!(e.c[4].to_string(), "1/24");
        // exp(2x) = exp(x)²
        let e2 = x.scale(&q(2)).exp().unwrap();
        let sq = e.mul(&e);
        assert_eq!(e2.trunc(4), sq.trunc(4));
        // exp(x)·exp(−x) = 1
        let em = x.scale(&q(-1)).exp().unwrap();
        let prod = e.mul(&em).trunc(4);
        assert_eq!(prod.c[0].to_string(), "1");
        for i in 1..5 {
            assert_eq!(prod.c[i].to_string(), "0", "i={i}");
        }
    }

    #[test]
    fn series_compose() {
        // exp(x) 代 x ← 2y：= exp(2y)
        let x = QSeries::var(4);
        let e = x.exp().unwrap();
        let mut f = vec![Q::ZERO; 5];
        f[1] = q(2);
        let e2 = e.compose(&QSeries { c: f }).unwrap();
        assert_eq!(e2.c[1].to_string(), "2");
        assert_eq!(e2.c[2].to_string(), "2");
        assert_eq!(e2.c[3].to_string(), "4/3");
        // 数值：exp(2) 的截断值（级数层验证）
        assert_eq!(e2.eval(&q(1)).to_string(), e.eval(&q(2)).to_string());
    }

    #[test]
    fn inverse_and_negative_powers() {
        // 1/(1−x) = 1 + x + x² + …
        let mut c = vec![Q::ZERO; 5];
        c[0] = q(1);
        c[1] = q(-1);
        let inv = QSeries { c }.inverse().unwrap();
        for i in 0..5 {
            assert_eq!(inv.c[i].to_string(), "1", "i={i}");
        }
        // exp(x)^{-2} = exp(−2x)
        let x = QSeries::var(4);
        let lhs = x.exp().unwrap().pow_signed(-2).unwrap();
        let rhs = x.scale(&q(-2)).exp().unwrap();
        assert_eq!(lhs, rhs);
        // 常数项为零不可逆
        assert!(x.inverse().is_err());
    }

    #[test]
    fn powi_basics() {
        let x = QSeries::var(3);
        let x2 = x.powi(2);
        assert_eq!(x2.c[2].to_string(), "1");
        assert_eq!(x2.c.len(), 4);
        let x0 = x.powi(0);
        assert_eq!(x0.c[0].to_string(), "1");
    }
}
