//! GW 递归族：GW（L_k 算子消 τ(1) 后代）→ GWhodge（D_{2k−1} 消 ch 插入）
//! → GWomega（H_k 消 ω 类插入）→ gw 前端（五种插入类型分流）。
//!
//! 递归内部 τ 指数用 i64（−1 出现于递归中间态，命中即不变量为零——
//! GWall 的 `member(-1,tlist)` 检查）。记忆化对应原版 `option remember`。
//! 环路插入可任意排序：gw 前端把交错的 a/b 环路换到规范序（a 在前）并计入符号。
//!
//! 记忆化键规范化（[`Normal`]）：偶类插入（τ(1)、τ(pt)、ch(O)、ch(L)）可交换，
//! 排序；奇类环路插入反交换，a、b 各自排序并记下置换符号，两个相同的奇类插入
//! 使不变量为零。原版只排序 τ₁（'reorder tau_1 for maximum efficiency'），同一
//! 不变量以不同插入顺序反复出现时各算一遍（z(3,2,8) 的未命中 410 万 → 15 万）。

use crate::gwloop::{LoopIns, fast_gwloop, gwloop as gwloop_slow};
use crate::qseries::bernoulli;
use sc_num::{One, Q, Zero};
use std::collections::HashMap;

type Key = (
    u8, // 0=GW, 1=GWhodge, 2=GWomega
    u32,
    u32,
    Vec<i64>, // tau_omega（仅 GWomega）
    Vec<i64>, // tau_hodge（GWhodge/GWomega）
    Vec<i64>, // tau_1
    Vec<(i64, i64)>,
    Vec<(i64, i64)>,
    Vec<i64>, // tau_p
    i64,      // deg（仅 GWomega）
);

pub struct Ctx {
    mem: HashMap<Key, Q>,
    pub fast: bool,
    /// (g, d) → (已算 Z 级数, 其截断阶)——pf 层跨调用复用（扩大截断时重算）
    pub zmem: std::collections::HashMap<(u32, u32), (Vec<Q>, u32)>,
}

impl Ctx {
    pub fn new() -> Self {
        Ctx {
            mem: HashMap::new(),
            fast: true,
            zmem: std::collections::HashMap::new(),
        }
    }
}

impl Default for Ctx {
    fn default() -> Self {
        Self::new()
    }
}

fn fact(n: u32) -> Q {
    let mut v = Q::ONE;
    for i in 2..=n {
        v *= Q::from(i as i64);
    }
    v
}

fn binom(n: u32, k: u32) -> Q {
    if k > n {
        return Q::ZERO;
    }
    let mut v = Q::ONE;
    for i in 0..k {
        v = v * Q::from((n - i) as i64) / Q::from((k - i) as i64);
    }
    v
}

fn sgn(e: i64) -> Q {
    Q::from(if e.rem_euclid(2) == 0 { 1 } else { -1 })
}

fn harm(a: u32, b: u32) -> Q {
    // Σ_{r=a..b} 1/r（a ≥ 1）
    let mut v = Q::ZERO;
    for r in a..=b {
        v += Q::ONE / Q::from(r as i64);
    }
    v
}

/// 递归入口参数的规范形（见模块注释）：值 = sign · ⟨规范序插入⟩。
struct Normal {
    t1: Vec<i64>,
    ta: Vec<(i64, i64)>,
    tb: Vec<(i64, i64)>,
    tp: Vec<i64>,
    sign: Q,
}

impl Normal {
    /// 不变量为零时返回 None：两个相同的奇类插入（反交换），或 a、b 环路个数
    /// 不等（递归只成对增删环路，闭式层个数不等即为零；GWomega 的 H 分支
    /// 按同一下标同时改 a、b，个数不等时还会越界）。
    fn new(t1: &[i64], ta: &[(i64, i64)], tb: &[(i64, i64)], tp: &[i64]) -> Option<Self> {
        if ta.len() != tb.len() {
            return None;
        }
        let (ta, sa) = sort_odd(ta)?;
        let (tb, sb) = sort_odd(tb)?;
        Some(Normal {
            t1: sorted(t1),
            ta,
            tb,
            tp: sorted(tp),
            sign: Q::from(sa * sb),
        })
    }
}

fn sorted(v: &[i64]) -> Vec<i64> {
    let mut v = v.to_vec();
    v.sort_unstable();
    v
}

/// 奇类插入按 (τ, 指标) 排序，返回置换符号；有两个相同的插入 → None。
fn sort_odd(v: &[(i64, i64)]) -> Option<(Vec<(i64, i64)>, i64)> {
    let mut w = v.to_vec();
    let mut sign = 1i64;
    // 插入排序（列表很短），每次相邻对换翻转符号
    for i in 1..w.len() {
        let mut j = i;
        while j > 0 && w[j - 1] > w[j] {
            w.swap(j - 1, j);
            sign = -sign;
            j -= 1;
        }
    }
    if w.windows(2).any(|p| p[0] == p[1]) {
        return None;
    }
    Some((w, sign))
}

/// tlist 检查：任一 τ 指数为 −1 → 不变量为零。
fn has_neg1(t1: &[i64], ta: &[(i64, i64)], tb: &[(i64, i64)], tp: &[i64]) -> bool {
    t1.contains(&-1) || tp.contains(&-1) || ta.iter().chain(tb).any(|&(t, _)| t == -1)
}

impl Ctx {
    fn base(&self, g: u32, d: u32, ta: &[(i64, i64)], tb: &[(i64, i64)], tp: &[i64]) -> Option<Q> {
        // 闭式层：i64 → u32（此处无 −1 已由上层检查）
        let ta_u: Vec<LoopIns> = ta
            .iter()
            .map(|&(t, i)| LoopIns {
                tau: t as u32,
                idx: i as u32,
            })
            .collect();
        let tb_u: Vec<LoopIns> = tb
            .iter()
            .map(|&(t, i)| LoopIns {
                tau: t as u32,
                idx: i as u32,
            })
            .collect();
        let tp_u: Vec<u32> = tp.iter().map(|&x| x as u32).collect();
        if ta_u.iter().any(|t| t.tau == u32::MAX) || tp_u.contains(&u32::MAX) {
            return Some(Q::ZERO);
        }
        Some(if self.fast {
            fast_gwloop(g, d, &tp_u, &ta_u, &tb_u)
        } else {
            gwloop_slow(g, d, &tp_u, &ta_u, &tb_u)
        })
    }

    /// GW：L_k 递归（k = τ₁[0]−1 ∈ {−1, 0, >0}）。
    pub fn gw_rec(
        &mut self,
        g: u32,
        d: u32,
        t1: &[i64],
        ta: &[(i64, i64)],
        tb: &[(i64, i64)],
        tp: &[i64],
    ) -> Q {
        let Some(nf) = Normal::new(t1, ta, tb, tp) else {
            return Q::ZERO;
        };
        let key = (
            0u8,
            g,
            d,
            vec![],
            vec![],
            nf.t1.clone(),
            nf.ta.clone(),
            nf.tb.clone(),
            nf.tp.clone(),
            0,
        );
        let val = match self.mem.get(&key) {
            Some(v) => v.clone(),
            None => {
                let v = self.gw_body(g, d, &nf.t1, &nf.ta, &nf.tb, &nf.tp);
                self.mem.insert(key, v.clone());
                v
            }
        };
        val * nf.sign
    }

    fn gw_body(
        &mut self,
        g: u32,
        d: u32,
        t1: &[i64],
        ta: &[(i64, i64)],
        tb: &[(i64, i64)],
        tp: &[i64],
    ) -> Q {
        if has_neg1(t1, ta, tb, tp) {
            return Q::ZERO;
        }
        let n1 = t1.len();
        if n1 == 0 {
            return self.base(g, d, ta, tb, tp).unwrap_or(Q::ZERO);
        }
        let k = t1[0] - 1;
        let tnew: Vec<i64> = t1[1..].to_vec();
        let mut answer = Q::ZERO;
        let rec = |s: &mut Self, t1: &[i64], ta: &[(i64, i64)], tb: &[(i64, i64)], tp: &[i64]| {
            s.gw_rec(g, d, t1, ta, tb, tp)
        };
        if k == -1 {
            // L_{−1}（string）：导数项——各 τ 指数减一
            for j in 0..tnew.len() {
                let mut v = tnew.clone();
                v[j] -= 1;
                answer += rec(self, &v, ta, tb, tp);
            }
            for j in 0..ta.len() {
                let mut va = ta.to_vec();
                va[j].0 -= 1;
                answer += rec(self, &tnew, &va, tb, tp);
            }
            for j in 0..tb.len() {
                let mut vb = tb.to_vec();
                vb[j].0 -= 1;
                answer += rec(self, &tnew, ta, &vb, tp);
            }
            for j in 0..tp.len() {
                let mut vp = tp.to_vec();
                vp[j] -= 1;
                answer += rec(self, &tnew, ta, tb, &vp);
            }
            // 常数项：τ₁ 与 τ_p 各一个零
            let (mut z1, mut zp) = (0usize, 0usize);
            let (mut i1, mut ip) = (0usize, 0usize);
            for (j, &v) in tnew.iter().enumerate() {
                if v == 0 {
                    z1 += 1;
                    i1 = j;
                }
            }
            for (j, &v) in tp.iter().enumerate() {
                if v == 0 {
                    zp += 1;
                    ip = j;
                }
            }
            if z1 > 0 && zp > 0 {
                let mut v1 = tnew.clone();
                v1.remove(i1);
                let mut vp = tp.to_vec();
                vp.remove(ip);
                answer += Q::from(z1 as i64 * zp as i64) * rec(self, &v1, ta, tb, &vp);
            }
            // 奇常数项：τ_a=τ_b=0 且 idx 匹配的对
            for ia in 0..ta.len() {
                if ta[ia].0 != 0 {
                    continue;
                }
                for ib in 0..tb.len() {
                    if tb[ib].0 == 0 && tb[ib].1 == ta[ia].1 {
                        let mut va = ta.to_vec();
                        va.remove(ia);
                        let mut vb = tb.to_vec();
                        vb.remove(ib);
                        // 源码 (−1)^(ia−1+na+ib)（Maple 1 基 ia,ib）：
                        // 0 基下 = ia + na + ib + 1（代数恒等，勿再"简化"）
                        let e = ia as i64 + ta.len() as i64 + ib as i64 + 1;
                        answer += sgn(e) * rec(self, &tnew, &va, &vb, tp);
                    }
                }
            }
        } else if k == 0 {
            // L_0（dilaton）
            let mut vp = tp.to_vec();
            vp.insert(0, 0);
            answer = Q::from(2 * g as i64 - 2) * rec(self, &tnew, ta, tb, &vp);
            let mult: i64 = tnew.iter().sum::<i64>()
                + ta.iter().map(|&x| x.0 + 1).sum::<i64>()
                + tb.iter().map(|&x| x.0).sum::<i64>()
                + tp.iter().map(|&x| x + 1).sum::<i64>();
            answer += Q::from(mult) * rec(self, &tnew, ta, tb, tp);
            for i in 0..tnew.len() {
                let mut v1 = tnew.clone();
                let l = v1.remove(i);
                let mut vp = tp.to_vec();
                vp.push(l - 1);
                answer += Q::from(2 - 2 * g as i64) * rec(self, &v1, ta, tb, &vp);
            }
            // 常数项：τ₁ 中零的个数 > 1
            let zeros = tnew.iter().filter(|&&x| x == 0).count();
            if zeros > 1 {
                let mut v1 = tnew.clone();
                let p1 = v1.iter().position(|&x| x == 0).unwrap();
                v1.remove(p1);
                let p2 = v1.iter().position(|&x| x == 0).unwrap();
                v1.remove(p2);
                answer +=
                    Q::from(2 - 2 * g as i64) * binom(zeros as u32, 2) * rec(self, &v1, ta, tb, tp);
            }
        } else {
            // L_k, k>0
            let k = k as u32;
            let mut vp = tp.to_vec();
            vp.insert(0, k as i64);
            answer = Q::from(2 * g as i64 - 2) * harm(1, k + 1) * rec(self, &tnew, ta, tb, &vp);
            for i in 0..tnew.len() {
                let l = tnew[i];
                let mut v = tnew.clone();
                v[i] = l + k as i64;
                answer += binom(k + l as u32, k + 1) * rec(self, &v, ta, tb, tp);
            }
            for i in 0..tp.len() {
                let l = tp[i];
                let mut vp = tp.to_vec();
                vp[i] = l + k as i64;
                answer += binom(k + l as u32 + 1, k + 1) * rec(self, &tnew, ta, tb, &vp);
            }
            for i in 0..ta.len() {
                let l = ta[i].0;
                let mut va = ta.to_vec();
                va[i].0 = l + k as i64;
                answer += binom(k + l as u32, k + 1) * rec(self, &tnew, &va, tb, tp);
            }
            for i in 0..tb.len() {
                let l = tb[i].0;
                let mut vb = tb.to_vec();
                vb[i].0 = l + k as i64;
                answer += binom(k + l as u32 + 1, k + 1) * rec(self, &tnew, ta, &vb, tp);
            }
            for i in 0..tnew.len() {
                let l = tnew[i];
                let mut v1 = tnew.clone();
                v1.remove(i);
                let mut vp = tp.to_vec();
                vp.push(k as i64 + l - 1);
                let mult = if l == 0 {
                    Q::ONE / Q::from((k + 1) as i64)
                } else {
                    binom(k + l as u32, k + 1) * harm(l as u32, k + l as u32)
                };
                answer += Q::from(2 - 2 * g as i64) * mult * rec(self, &v1, ta, tb, &vp);
            }
            for l in 0..k.saturating_sub(1) {
                let mut vp = tp.to_vec();
                vp.push(l as i64);
                vp.push((k - l - 2) as i64);
                answer += Q::from(1 - g as i64) / Q::from((k + 1) as i64) / binom(k, l + 1)
                    * rec(self, &tnew, ta, tb, &vp);
            }
        }
        answer
    }

    /// GWhodge：D_{2k−1} 递归（非奇 hodge → 0）。
    pub fn gwhodge_rec(
        &mut self,
        g: u32,
        d: u32,
        th: &[i64],
        t1: &[i64],
        ta: &[(i64, i64)],
        tb: &[(i64, i64)],
        tp: &[i64],
    ) -> Q {
        if th.iter().any(|&h| h <= 0 || h % 2 == 0) {
            return Q::ZERO;
        }
        let Some(nf) = Normal::new(t1, ta, tb, tp) else {
            return Q::ZERO;
        };
        let th = sorted(th);
        let key = (
            1u8,
            g,
            d,
            vec![],
            th.clone(),
            nf.t1.clone(),
            nf.ta.clone(),
            nf.tb.clone(),
            nf.tp.clone(),
            0,
        );
        let val = match self.mem.get(&key) {
            Some(v) => v.clone(),
            None => {
                let v = self.gwhodge_body(g, d, &th, &nf.t1, &nf.ta, &nf.tb, &nf.tp);
                self.mem.insert(key, v.clone());
                v
            }
        };
        val * nf.sign
    }

    fn gwhodge_body(
        &mut self,
        g: u32,
        d: u32,
        th: &[i64],
        t1: &[i64],
        ta: &[(i64, i64)],
        tb: &[(i64, i64)],
        tp: &[i64],
    ) -> Q {
        if has_neg1(t1, ta, tb, tp) {
            return Q::ZERO;
        }
        if th.is_empty() {
            return self.gw_rec(g, d, t1, ta, tb, tp);
        }
        let k = ((th[0] + 1) / 2) as u32;
        let tnew: Vec<i64> = th[1..].to_vec();
        let gh = |s: &mut Self,
                  th: &[i64],
                  t1: &[i64],
                  ta: &[(i64, i64)],
                  tb: &[(i64, i64)],
                  tp: &[i64]| s.gwhodge_rec(g, d, th, t1, ta, tb, tp);
        let mut answer = {
            let mut v = t1.to_vec();
            v.insert(0, 2 * k as i64);
            gh(self, &tnew, &v, ta, tb, tp)
        };
        for i in 0..t1.len() {
            let mut v = t1.to_vec();
            v[i] += 2 * k as i64 - 1;
            answer -= gh(self, &tnew, &v, ta, tb, tp);
        }
        for i in 0..tp.len() {
            let mut v = tp.to_vec();
            v[i] += 2 * k as i64 - 1;
            answer -= gh(self, &tnew, t1, ta, tb, &v);
        }
        for i in 0..ta.len() {
            let mut va = ta.to_vec();
            va[i].0 += 2 * k as i64 - 1;
            answer -= gh(self, &tnew, t1, &va, tb, tp);
        }
        for i in 0..tb.len() {
            let mut vb = tb.to_vec();
            vb[i].0 += 2 * k as i64 - 1;
            answer -= gh(self, &tnew, t1, ta, &vb, tp);
        }
        for l in 0..2 * k as i64 - 1 {
            let mut sum = {
                let mut v1 = t1.to_vec();
                v1.push(l);
                let mut vp = tp.to_vec();
                vp.push(2 * k as i64 - 2 - l);
                gh(self, &tnew, &v1, ta, tb, &vp)
            };
            for ia in 1..=g as i64 {
                let mut va = ta.to_vec();
                va.insert(0, (l, ia));
                let mut vb = tb.to_vec();
                vb.insert(0, (2 * k as i64 - 2 - l, ia));
                sum -= sgn(ta.len() as i64) * gh(self, &tnew, t1, &va, &vb, tp);
            }
            answer += sgn(l) * sum;
        }
        answer * (bernoulli(2 * k) / fact(2 * k))
    }

    /// GWomega：H_k 递归（ω = ch(Rπf*L)，deg = deg L）。
    pub fn gwomega_rec(
        &mut self,
        g: u32,
        d: u32,
        to: &[i64],
        th: &[i64],
        t1: &[i64],
        ta: &[(i64, i64)],
        tb: &[(i64, i64)],
        tp: &[i64],
        deg: i64,
    ) -> Q {
        if to.iter().any(|&w| w <= 0) || th.iter().any(|&h| h <= 0 || h % 2 == 0) {
            return Q::ZERO;
        }
        let Some(nf) = Normal::new(t1, ta, tb, tp) else {
            return Q::ZERO;
        };
        let (to, th) = (sorted(to), sorted(th));
        let key = (
            2u8,
            g,
            d,
            to.clone(),
            th.clone(),
            nf.t1.clone(),
            nf.ta.clone(),
            nf.tb.clone(),
            nf.tp.clone(),
            deg,
        );
        let val = match self.mem.get(&key) {
            Some(v) => v.clone(),
            None => {
                let v = self.gwomega_body(g, d, &to, &th, &nf.t1, &nf.ta, &nf.tb, &nf.tp, deg);
                self.mem.insert(key, v.clone());
                v
            }
        };
        val * nf.sign
    }

    fn gwomega_body(
        &mut self,
        g: u32,
        d: u32,
        to: &[i64],
        th: &[i64],
        t1: &[i64],
        ta: &[(i64, i64)],
        tb: &[(i64, i64)],
        tp: &[i64],
        deg: i64,
    ) -> Q {
        if has_neg1(t1, ta, tb, tp) {
            return Q::ZERO;
        }
        if to.is_empty() {
            return self.gwhodge_rec(g, d, th, t1, ta, tb, tp);
        }
        let w = to[0];
        let tnew: Vec<i64> = to[1..].to_vec();
        let om = |s: &mut Self,
                  to: &[i64],
                  th: &[i64],
                  t1: &[i64],
                  ta: &[(i64, i64)],
                  tb: &[(i64, i64)],
                  tp: &[i64]| s.gwomega_rec(g, d, to, th, t1, ta, tb, tp, deg);
        let dq = Q::from(deg);
        let mut answer = Q::ZERO;
        if w == 1 {
            // H_1
            let mut vp = tp.to_vec();
            vp.insert(0, 1);
            // 源码 answer := −deg/2·(...)——首项为负
            answer -= (dq.clone() / Q::from(2)) * om(self, &tnew, th, t1, ta, tb, &vp);
            for i in 0..t1.len() {
                let l = t1[i];
                let mut v = t1.to_vec();
                v.remove(i);
                let mut vp = tp.to_vec();
                vp.insert(0, l);
                answer += (dq.clone() / Q::from(2)) * om(self, &tnew, th, &v, ta, tb, &vp);
            }
            {
                let mut v = t1.to_vec();
                v.insert(0, 2);
                answer += Q::from(1) / Q::from(12) * om(self, &tnew, th, &v, ta, tb, tp);
            }
            for i in 0..t1.len() {
                let mut v = t1.to_vec();
                v[i] += 1;
                answer -= Q::from(1) / Q::from(12) * om(self, &tnew, th, &v, ta, tb, tp);
            }
            for i in 0..tp.len() {
                let mut vp = tp.to_vec();
                vp[i] += 1;
                answer -= Q::from(1) / Q::from(12) * om(self, &tnew, th, t1, ta, tb, &vp);
            }
            for i in 0..ta.len() {
                let mut va = ta.to_vec();
                va[i].0 += 1;
                let mut vb = tb.to_vec();
                vb[i].0 += 1;
                answer -= Q::from(1) / Q::from(12) * om(self, &tnew, th, t1, &va, tb, tp);
                answer -= Q::from(1) / Q::from(12) * om(self, &tnew, th, t1, ta, &vb, tp);
            }
            {
                let mut v = t1.to_vec();
                v.insert(0, 0);
                let mut vp = tp.to_vec();
                vp.insert(0, 0);
                let tv = om(self, &tnew, th, &v, ta, tb, &vp);
                answer += Q::from(1) / Q::from(12) * tv;
            }
            for i in 1..=g as i64 {
                let mut va = ta.to_vec();
                va.insert(0, (0, i));
                let mut vb = tb.to_vec();
                vb.insert(0, (0, i));
                answer += Q::from(1) / Q::from(12)
                    * sgn(ta.len() as i64 + 1)
                    * om(self, &tnew, th, t1, &va, &vb, tp);
            }
            for i in 0..tnew.len() {
                if tnew[i] == 1 {
                    let mut tr = tnew.clone();
                    tr.remove(i);
                    let mut vp = tp.to_vec();
                    vp.insert(0, 0);
                    vp.insert(0, 0);
                    answer -= (dq.clone() / Q::from(24)) * om(self, &tr, th, t1, ta, tb, &vp);
                }
            }
        } else if w % 2 == 1 {
            // H_{2k−1}
            let k = ((w + 1) / 2) as u32;
            {
                let mut v = t1.to_vec();
                v.insert(0, 2 * k as i64);
                answer = om(self, &tnew, th, &v, ta, tb, tp);
            }
            for i in 0..t1.len() {
                let mut v = t1.to_vec();
                v[i] += 2 * k as i64 - 1;
                answer -= om(self, &tnew, th, &v, ta, tb, tp);
            }
            for i in 0..tp.len() {
                let mut vp = tp.to_vec();
                vp[i] += 2 * k as i64 - 1;
                answer -= om(self, &tnew, th, t1, ta, tb, &vp);
            }
            for i in 0..ta.len() {
                let mut va = ta.to_vec();
                va[i].0 += 2 * k as i64 - 1;
                let mut vb = tb.to_vec();
                vb[i].0 += 2 * k as i64 - 1;
                answer -= om(self, &tnew, th, t1, &va, tb, tp);
                answer -= om(self, &tnew, th, t1, ta, &vb, tp);
            }
            for l in 0..2 * k as i64 - 1 {
                let mut v = t1.to_vec();
                v.insert(0, l);
                let mut vp = tp.to_vec();
                vp.insert(0, 2 * k as i64 - 2 - l);
                answer += sgn(l) * om(self, &tnew, th, &v, ta, tb, &vp);
                for i in 1..=g as i64 {
                    let mut va = ta.to_vec();
                    va.insert(0, (l, i));
                    let mut vb = tb.to_vec();
                    vb.insert(0, (2 * k as i64 - 2 - l, i));
                    answer -= sgn(l + ta.len() as i64) * om(self, &tnew, th, t1, &va, &vb, tp);
                }
            }
            for i in 0..tnew.len() {
                if tnew[i] == 1 {
                    let mut tr = tnew.clone();
                    tr.remove(i);
                    for l in 0..2 * k as i64 - 1 {
                        let mut vp = tp.to_vec();
                        vp.insert(0, 2 * k as i64 - 2 - l);
                        vp.insert(0, l);
                        answer -=
                            (dq.clone() / Q::from(2)) * sgn(l) * om(self, &tr, th, t1, ta, tb, &vp);
                    }
                }
            }
            answer *= bernoulli(2 * k) / fact(2 * k);
        } else {
            // H_{2k}
            let k = (w / 2) as u32;
            {
                let mut vp = tp.to_vec();
                vp.insert(0, 2 * k as i64);
                answer = om(self, &tnew, th, t1, ta, tb, &vp);
            }
            for i in 0..t1.len() {
                let l = t1[i];
                let mut v = t1.to_vec();
                v.remove(i);
                let mut vp = tp.to_vec();
                vp.insert(0, l + 2 * k as i64 - 1);
                answer -= om(self, &tnew, th, &v, ta, tb, &vp);
            }
            for l in 0..2 * k as i64 - 1 {
                let mut vp = tp.to_vec();
                vp.insert(0, 2 * k as i64 - 2 - l);
                vp.insert(0, l);
                answer += (Q::from(1) / Q::from(2)) * sgn(l) * om(self, &tnew, th, t1, ta, tb, &vp);
            }
            answer *= dq * bernoulli(2 * k) / fact(2 * k);
        }
        answer
    }
}

// ===================== gw 前端 =====================

/// 插入类型（对应 GWall 文档的 v_i 编码）：
/// `T1(k)`=τ_k(1)（文档值 k）；`Pt(k)`=τ_k(pt)（值 −k，k=0 用 −1/2）；
/// `A(k,i)`/`B(k,i)`=τ_k(a_i)/τ_k(b_i)（值 [k,±i]）；
/// `ChO(k)`=ch_k(Rπf*O)（值 [k]，须奇）；`ChL(k)`=ch_k(Rπf*L)（值 [−k]）；
/// `LDeg(deg)` 改变 ChL 的线丛度数（值 [[deg]]，默认 2g−2 即 K）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Ins {
    T1(u32),
    Pt(u32),
    A(u32, u32),
    B(u32, u32),
    ChO(u32),
    ChL(u32),
    LDeg(i64),
}

/// 主入口：⟨Π insert⟩_{g,d}（不连通 local curve 不变量）。
pub fn gw(g: u32, d: u32, insert: &[Ins]) -> Q {
    let mut ctx = Ctx::new();
    gw_ctx(&mut ctx, g, d, insert)
}

/// 环路插入可按任意顺序给出：递归族的规范序是"全部 a 环路在前、b 环路在后"
/// （GWall 文档的约定），环路后代反交换，交错输入按把 a 逐个移过其前面的
/// b 所需的对换数取符号；偶类插入与一切交换，不影响符号。
pub fn gw_ctx(ctx: &mut Ctx, g: u32, d: u32, insert: &[Ins]) -> Q {
    let mut omega_degree = 2 * g as i64 - 2;
    let (mut t1, mut ta, mut tb, mut tp) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut to, mut th) = (Vec::new(), Vec::new());
    let mut swaps = 0usize;
    for ins in insert {
        match *ins {
            Ins::T1(k) => t1.push(k as i64),
            Ins::Pt(k) => tp.push(k as i64),
            Ins::A(k, i) => {
                swaps += tb.len();
                ta.push((k as i64, i as i64));
            }
            Ins::B(k, i) => tb.push((k as i64, i as i64)),
            Ins::ChO(k) => th.push(k as i64),
            Ins::ChL(k) => to.push(k as i64),
            Ins::LDeg(deg) => omega_degree = deg,
        }
    }
    let v = if !to.is_empty() {
        ctx.gwomega_rec(g, d, &to, &th, &t1, &ta, &tb, &tp, omega_degree)
    } else if !th.is_empty() {
        ctx.gwhodge_rec(g, d, &th, &t1, &ta, &tb, &tp)
    } else {
        ctx.gw_rec(g, d, &t1, &ta, &tb, &tp)
    };
    if swaps % 2 == 1 { -v } else { v }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// String 方程（L_{−1} 的数学内容）强独立检验：
    /// ⟨τ_0(1)·X⟩ = Σ_x ⟨X|τ_x 降一阶⟩——右侧不含 τ(1)，直达闭式层。
    #[test]
    fn string_equation() {
        let mut ctx = Ctx::new();
        // 例 1：⟨τ_0(1) τ_1(pt) τ_0(pt)⟩ = ⟨τ_0(pt)τ_0(pt)⟩ + ⟨τ_1(pt)·(τ_{−1})=0 侧⟩
        let lhs = gw_ctx(&mut ctx, 2, 1, &[Ins::T1(0), Ins::Pt(1), Ins::Pt(0)]);
        let rhs = gw_ctx(&mut ctx, 2, 1, &[Ins::Pt(0), Ins::Pt(0)]);
        assert_eq!(lhs.to_string(), rhs.to_string());
        // 例 2：⟨τ_0(1) τ_1(pt) τ_1(pt)⟩ = ⟨τ_0(pt)τ_1(pt)⟩ + ⟨τ_1(pt)τ_0(pt)⟩
        let lhs = gw_ctx(&mut ctx, 1, 2, &[Ins::T1(0), Ins::Pt(1), Ins::Pt(1)]);
        let rhs = gw_ctx(&mut ctx, 1, 2, &[Ins::Pt(0), Ins::Pt(1)])
            + gw_ctx(&mut ctx, 1, 2, &[Ins::Pt(1), Ins::Pt(0)]);
        assert_eq!(lhs.to_string(), rhs.to_string());
        // 例 3：τ_0(1) 与 loop 后代：⟨τ_0(1) τ_1(a)τ_0(b)⟩ = ⟨τ_0(a)τ_0(b)⟩
        let lhs = gw_ctx(&mut ctx, 2, 1, &[Ins::T1(0), Ins::A(1, 1), Ins::B(0, 1)]);
        let rhs = gw_ctx(&mut ctx, 2, 1, &[Ins::A(0, 1), Ins::B(0, 1)]);
        assert_eq!(lhs.to_string(), rhs.to_string());
        // 例 4：τ_0(1) 与 τ_1(1)：⟨τ_0(1)τ_1(1)X⟩ = ⟨τ_0(1)X⟩（τ_1(1) 降为 τ_0(1)）
        let lhs = gw_ctx(&mut ctx, 2, 1, &[Ins::T1(0), Ins::T1(1), Ins::Pt(0)]);
        let rhs = gw_ctx(&mut ctx, 2, 1, &[Ins::T1(0), Ins::Pt(0)]);
        assert_eq!(lhs.to_string(), rhs.to_string());
    }

    /// Dilaton 方程的线性结构：⟨τ_1(1)·X⟩ = (a + n)·⟨X⟩，n = X 的点插入数。
    /// 不连通理论含常映射分量（τ_1(1) 落在收缩分量上另有 Hodge 积分贡献），
    /// a 不是连通情形的 2h−2，这里只检验 a 与 n 无关（n = 2, 3）。
    #[test]
    fn dilaton_equation() {
        let mut ctx = Ctx::new();
        let mut a_of = |n: usize| {
            let x: Vec<Ins> = vec![Ins::Pt(0); n];
            let mut l = vec![Ins::T1(1)];
            l.extend_from_slice(&x);
            let vx = gw_ctx(&mut ctx, 2, 1, &x);
            assert!(vx != Q::ZERO);
            gw_ctx(&mut ctx, 2, 1, &l) / vx - Q::from(n as i64)
        };
        assert_eq!(a_of(2), a_of(3), "dilaton 系数应与 n 无关");
    }

    /// 记忆化键规范化：同一上下文中先算打乱顺序的插入、再算规范顺序的（后者
    /// 命中前者存下的规范形），两者只差 a、b 各自的置换符号；偶类插入顺序无关。
    #[test]
    fn memo_normalization_signs() {
        // 值与 v0.1.0（无规范化）逐字一致
        let shuffled_ins = [
            Ins::ChL(1),
            Ins::A(1, 2),
            Ins::A(0, 2),
            Ins::B(0, 2),
            Ins::ChO(5),
            Ins::B(1, 2),
        ];
        let mut ctx = Ctx::new();
        let shuffled = gw_ctx(&mut ctx, 2, 3, &shuffled_ins);
        assert_eq!(shuffled.to_string(), "41833459363/111476736000");
        let canon = gw_ctx(
            &mut ctx,
            2,
            3,
            &[
                Ins::A(0, 2),
                Ins::A(1, 2),
                Ins::B(0, 2),
                Ins::B(1, 2),
                Ins::ChO(5),
                Ins::ChL(1),
            ],
        );
        // a 侧差一次对换：符号相反
        assert_eq!(canon, -shuffled.clone());
        // 新上下文先算规范序、再算打乱序（命中顺序反过来）
        let mut ctx2 = Ctx::new();
        assert_eq!(
            gw_ctx(
                &mut ctx2,
                2,
                3,
                &[
                    Ins::A(0, 2),
                    Ins::A(1, 2),
                    Ins::B(0, 2),
                    Ins::B(1, 2),
                    Ins::ChO(5),
                    Ins::ChL(1)
                ]
            ),
            canon
        );
        assert_eq!(gw_ctx(&mut ctx2, 2, 3, &shuffled_ins), shuffled);
        // 两个相同的奇类插入：反交换 ⇒ 零
        assert_eq!(
            gw(
                2,
                1,
                &[Ins::A(0, 1), Ins::A(0, 1), Ins::B(0, 1), Ins::B(1, 1)]
            ),
            Q::ZERO
        );
    }

    /// 回归：a、b 环路个数不等时不变量为零。GWomega 的 H 分支按同一下标同时
    /// 改 a、b，原先在 a 多于 b 时越界 panic（批处理整个进程退出）。
    #[test]
    fn unequal_loop_counts_vanish() {
        assert_eq!(gw(1, 1, &[Ins::A(1, 1), Ins::ChL(3)]), Q::ZERO);
        assert_eq!(
            gw(
                3,
                2,
                &[
                    Ins::T1(2),
                    Ins::A(2, 2),
                    Ins::B(2, 2),
                    Ins::Pt(2),
                    Ins::A(1, 2),
                    Ins::ChL(3)
                ]
            ),
            Q::ZERO
        );
        assert_eq!(gw(2, 1, &[Ins::B(0, 1), Ins::ChL(1)]), Q::ZERO);
    }

    /// 环路反交换：交错给出 a/b 环路与规范序（a 在前）只差对换符号；
    /// 偶类插入的位置不影响符号。
    #[test]
    fn loop_order_sign() {
        let mut ctx = Ctx::new();
        let canon = gw_ctx(
            &mut ctx,
            2,
            1,
            &[Ins::A(0, 1), Ins::A(1, 2), Ins::B(0, 2), Ins::B(1, 1)],
        );
        assert!(canon != Q::ZERO);
        // b(·,2) 在 a(1,2) 前：一次对换
        let one = gw_ctx(
            &mut ctx,
            2,
            1,
            &[Ins::A(0, 1), Ins::B(0, 2), Ins::A(1, 2), Ins::B(1, 1)],
        );
        assert_eq!(one, -canon.clone());
        // 两个 b 都在 a 前：四次对换
        let four = gw_ctx(
            &mut ctx,
            2,
            1,
            &[Ins::B(0, 2), Ins::B(1, 1), Ins::A(0, 1), Ins::A(1, 2)],
        );
        assert_eq!(four, canon.clone());
        // 偶类插在中间不改变符号
        let mixed = gw_ctx(
            &mut ctx,
            2,
            1,
            &[
                Ins::A(0, 1),
                Ins::Pt(0),
                Ins::A(1, 2),
                Ins::T1(0),
                Ins::B(0, 2),
                Ins::B(1, 1),
            ],
        );
        let mixed_canon = gw_ctx(
            &mut ctx,
            2,
            1,
            &[
                Ins::Pt(0),
                Ins::T1(0),
                Ins::A(0, 1),
                Ins::A(1, 2),
                Ins::B(0, 2),
                Ins::B(1, 1),
            ],
        );
        assert_eq!(mixed, mixed_canon);
    }

    /// 前端分流与消没：非奇 ChO → 0；无 τ(1)/ch 插入直达闭式层。
    #[test]
    fn frontend_dispatch() {
        let mut ctx = Ctx::new();
        // 纯点后代：gw 前端 == gwpoint 直达
        let a = gw_ctx(&mut ctx, 2, 2, &[Ins::Pt(0), Ins::Pt(1)]);
        let b = crate::gwpoint::gwpoint(2, 2, &[0, 1]);
        assert_eq!(a.to_string(), b.to_string());
        // 偶 ch → 0
        assert_eq!(
            gw_ctx(&mut ctx, 2, 1, &[Ins::ChO(2), Ins::Pt(0), Ins::Pt(0)]).to_string(),
            "0"
        );
        // 环路成对：⟨τ_1(a₁)τ_2(b₁)⟩ 走 GW→loop 层
        let c = gw_ctx(&mut ctx, 3, 1, &[Ins::A(1, 1), Ins::B(2, 1)]);
        let direct = crate::gwloop::gwloop(
            3,
            1,
            &[],
            &[LoopIns { tau: 1, idx: 1 }],
            &[LoopIns { tau: 2, idx: 1 }],
        );
        assert_eq!(c.to_string(), direct.to_string());
    }

    /// GWomega 的 deg 约定：[[deg]] 插入改变 ChL 的线丛度。
    #[test]
    fn omega_degree_flag() {
        let mut ctx = Ctx::new();
        // 默认 deg = 2g−2（K）；显式同值应一致
        let a = gw_ctx(
            &mut ctx,
            2,
            1,
            &[Ins::ChL(1), Ins::Pt(0), Ins::Pt(0), Ins::Pt(0)],
        );
        let b = gw_ctx(
            &mut ctx,
            2,
            1,
            &[
                Ins::ChL(1),
                Ins::Pt(0),
                Ins::Pt(0),
                Ins::Pt(0),
                Ins::LDeg(2),
            ],
        );
        assert_eq!(a.to_string(), b.to_string());
    }
}
