//! Q 上的线性方程组——原版 `math/linsys.{h,cc}` 的直译。
//!
//! 增量式高斯消元：方程逐行加入，随时报告 矛盾 / 不定 / 唯一可解。
//! 布局与原版一致：`(maxsize+1)²` 稠密矩阵，第 `maxsize` 列为右端项，
//! 第 `rows` 行是正在构造的当前方程。变量以整数命名（Virasoro 用
//! `(xdim+1)·d + Hpow` 编码不变量），消元过程中会交换列（即重命名）。

use std::collections::HashMap;

use crate::num::{One, Q, Zero};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LinsysT {
    Contradiction,
    Undetermined,
    Solvable,
}

pub struct Linsys {
    maxsize: usize,
    varnames: Vec<i32>,
    entries: Vec<Q>, // (maxsize+1)*(maxsize+1)；列 maxsize 为 RHS
    rows: usize,
    cols: usize,
    /// 变量名 → 当前列号（加速 add_var 的线性查找；原版为 O(n) 扫描，
    /// 语义一致）
    by_name: HashMap<i32, usize>,
}

use std::fmt::Write as _;

impl Linsys {
    pub fn new(size: usize) -> Self {
        let n = size + 1;
        Linsys {
            maxsize: size,
            varnames: vec![0; size],
            entries: vec![Q::ZERO; n * n],
            rows: 0,
            cols: 0,
            by_name: HashMap::new(),
        }
    }

    #[inline]
    fn e(&self, r: usize, c: usize) -> &Q {
        &self.entries[r * (self.maxsize + 1) + c]
    }

    #[inline]
    fn e_mut(&mut self, r: usize, c: usize) -> &mut Q {
        let n = self.maxsize + 1;
        &mut self.entries[r * n + c]
    }

    /// 引入（或找到）变量。返回内部列号；超容量返回 `None`。
    pub fn add_var(&mut self, v: i32) -> Option<usize> {
        if let Some(&c) = self.by_name.get(&v) {
            return Some(c);
        }
        let i = self.cols;
        if i >= self.maxsize {
            return None;
        }
        self.cols += 1;
        self.varnames[i] = v;
        self.by_name.insert(v, i);
        for j in 0..=self.rows {
            *self.e_mut(j, i) = Q::ZERO;
        }
        Some(i)
    }

    /// 当前方程加一项 `coeff·x[v]`（同变量累加）。
    pub fn add_var_term(&mut self, v: i32, coeff: Q) -> bool {
        match self.add_var(v) {
            Some(i) => {
                *self.e_mut(self.rows, i) += coeff;
                true
            }
            None => false,
        }
    }

    /// 当前方程加常数项（右端）。
    pub fn add_const_term(&mut self, c: Q) {
        *self.e_mut(self.rows, self.maxsize) += c;
    }

    /// 提交当前方程并消元。
    pub fn send_row(&mut self) -> LinsysT {
        let rows = self.rows;
        let cols = self.cols;
        let maxsize = self.maxsize;

        // 用已有方程消去新方程的前 rows 列
        for i in 0..rows {
            let c = self.e(rows, i).clone() / self.e(i, i);
            *self.e_mut(rows, i) = Q::ZERO;
            for j in rows..cols {
                let sub = &c * self.e(i, j);
                *self.e_mut(rows, j) -= sub;
            }
            let sub = &c * self.e(i, maxsize);
            *self.e_mut(rows, maxsize) -= sub;
        }

        // 找新方程中第一个非零系数
        let mut found = None;
        for i in rows..cols {
            if !crate::num::is_zero(self.e(rows, i)) {
                found = Some(i);
                break;
            }
        }

        if crate::env_flag("SC_GROWI_VIRASORO") {
            eprintln!(
                "[linsys] send_row: rows={} cols={} found={found:?} rhs={:?}",
                rows,
                cols,
                self.e(rows, maxsize)
            );
        }
        let Some(i) = found else {
            if crate::num::is_zero(self.e(rows, maxsize)) {
                return if cols > rows {
                    LinsysT::Undetermined
                } else {
                    LinsysT::Solvable
                };
            }
            return LinsysT::Contradiction;
        };

        // 交换变量列 i 与 rows（重命名）
        if i != rows {
            self.varnames.swap(i, rows);
            self.by_name.insert(self.varnames[i], i);
            self.by_name.insert(self.varnames[rows], rows);
            for j in 0..=rows {
                self.entries
                    .swap(j * (maxsize + 1) + i, j * (maxsize + 1) + rows);
            }
        }

        // 消去新主元上方的项
        for i in 0..rows {
            let c = self.e(i, rows).clone() / self.e(rows, rows);
            *self.e_mut(i, rows) = Q::ZERO;
            for j in rows + 1..cols {
                let sub = &c * self.e(rows, j);
                *self.e_mut(i, j) -= sub;
            }
            let sub = &c * self.e(rows, maxsize);
            *self.e_mut(i, maxsize) -= sub;
        }

        // 开新行
        self.rows += 1;
        for i in 0..self.cols {
            *self.e_mut(self.rows, i) = Q::ZERO;
        }
        *self.e_mut(self.rows, maxsize) = Q::ZERO;

        if cols > self.rows {
            LinsysT::Undetermined
        } else {
            LinsysT::Solvable
        }
    }

    /// 变量个数。
    pub fn numvars(&self) -> usize {
        self.cols
    }

    /// 第 i 个（消元后顺序的）变量的名字。
    pub fn varname(&self, i: usize) -> i32 {
        self.varnames[i]
    }

    /// 变量的值（系统唯一可解时）。未定 → `None`。
    pub fn value_by_number(&self, i: usize) -> Option<Q> {
        if self.cols > self.rows {
            return None;
        }
        Some(self.e(i, self.maxsize).clone() / self.e(i, i))
    }

    /// 按变量名取值。
    pub fn value_by_name(&self, v: i32) -> Option<Q> {
        if self.cols > self.rows {
            return None;
        }
        self.by_name
            .get(&v)
            .map(|&j| self.e(j, self.maxsize).clone() / self.e(j, j))
    }

    /// 调试：当前（未提交）方程行。
    pub fn debug_current_row(&self) -> String {
        let mut out = String::new();
        for j in 0..self.cols {
            out.push_str(&format!(
                "{:?}x[{}] ",
                self.e(self.rows, j),
                self.varnames[j]
            ));
        }
        let _ = write!(out, "= {:?}", self.e(self.rows, self.maxsize));
        out
    }

    /// 行列式（方阵时）。
    pub fn det(&self) -> Option<Q> {
        if self.cols != self.rows {
            return None;
        }
        let mut result = Q::ONE;
        for i in (0..self.rows).rev() {
            result *= self.e(i, i);
        }
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::q;

    #[test]
    fn solves_unique_system() {
        let mut s = Linsys::new(3);
        // x1 + x2 = 3 ; x1 - x2 = 1 → x1=2, x2=1
        assert!(s.add_var_term(1, q(1)));
        assert!(s.add_var_term(2, q(1)));
        s.add_const_term(q(3));
        assert_eq!(s.send_row(), LinsysT::Undetermined);
        assert!(s.add_var_term(1, q(1)));
        assert!(s.add_var_term(2, q(-1)));
        s.add_const_term(q(1));
        assert_eq!(s.send_row(), LinsysT::Solvable);
        assert_eq!(s.value_by_name(1), Some(q(2)));
        assert_eq!(s.value_by_name(2), Some(q(1)));
    }

    #[test]
    fn detects_contradiction() {
        let mut s = Linsys::new(2);
        assert!(s.add_var_term(1, q(1)));
        s.add_const_term(q(1));
        // 单变量首方程：cols==rows==1，立即可解
        assert_eq!(s.send_row(), LinsysT::Solvable);
        assert!(s.add_var_term(1, q(1)));
        s.add_const_term(q(2));
        // 化简后 0=1 → 矛盾
        assert_eq!(s.send_row(), LinsysT::Contradiction);
    }

    #[test]
    fn dependent_equation() {
        let mut s = Linsys::new(2);
        assert!(s.add_var_term(1, q(1)));
        s.add_const_term(q(1));
        assert_eq!(s.send_row(), LinsysT::Solvable);
        assert!(s.add_var_term(1, q(2)));
        s.add_const_term(q(2));
        // 2x=2 与 x=1 相依赖：化简后 0=0，系统保持可解
        assert_eq!(s.send_row(), LinsysT::Solvable);
        assert_eq!(s.value_by_name(1), Some(q(1)));
    }
}
