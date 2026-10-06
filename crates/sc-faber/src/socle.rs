//! socle 公式：strata 顶点上 κ-ψ 单项式顶度积分的闭式。
//!
//! 移植自 admcycles `DR/evaluation.py`（Pixton 的公式）。这是 M3 的核心：
//! **求值不需要 FZ 关系矩阵**——每个稳定图顶点独立求值：
//! - g=0：CT 公式（Kontsevich 多项系数 + Bernoulli 常数）；
//! - g≥1：STsum = κ 指标集合分划求和 × STrecur（纯 ψ 顶度递归：
//!   string + 两项 dilaton 型修正 + 分离边界分裂项，除以 (2g+n−1)(2g+n−2)）。
//!
//! 带记忆化（L1 进程内；M3c 接 SQLite L2——STrecur 的子问题跨查询/跨会话
//! 复用极多，正是持久化设计的用武之地）。

use sc_num::{One, Q, Zero, is_zero};
use std::cell::RefCell;
use std::collections::HashMap;

thread_local! {
    static MEMO: RefCell<HashMap<Vec<u32>, Q>> = RefCell::new(HashMap::new());
}

/// STrecur 结果的 SQLite L2（跨会话持久；STrecur 子问题复用率高，
/// 是 Faber TAU 表的现代等价物）。
static L2: std::sync::OnceLock<std::sync::Mutex<Option<rusqlite::Connection>>> =
    std::sync::OnceLock::new();

const ALGO_VERSION: &str = "sc-faber 0.1.0";

fn l2() -> &'static std::sync::Mutex<Option<rusqlite::Connection>> {
    L2.get_or_init(|| std::sync::Mutex::new(None))
}

/// 卸下持久缓存（之后只用进程内记忆化）。
pub fn detach_cache() {
    *l2().lock().unwrap() = None;
}

/// 当前是否挂接了持久缓存（进程全局）。
pub fn cache_attached() -> bool {
    l2().lock().unwrap().is_some()
}

/// 挂接持久缓存文件（不存在则创建）。L2 是进程全局的，重复调用替换连接。
pub fn attach_cache(path: &str) -> Result<(), String> {
    let conn = rusqlite::Connection::open(path).map_err(|e| format!("打开缓存失败：{e}"))?;
    conn.pragma_update(None, "journal_mode", "WAL")
        .map_err(|e| e.to_string())?;
    conn.pragma_update(None, "synchronous", "NORMAL")
        .map_err(|e| e.to_string())?;
    conn.execute_batch(
        "CREATE TABLE IF NOT EXISTS st (
            k BLOB PRIMARY KEY,
            v TEXT NOT NULL
        ) WITHOUT ROWID;",
    )
    .map_err(|e| e.to_string())?;
    *l2().lock().unwrap() = Some(conn);
    Ok(())
}

/// LEB128 编码（含版本与种类前缀）。
fn enc_key(kind: &str, nums: &[u32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(nums.len() + 8);
    out.extend_from_slice(ALGO_VERSION.as_bytes());
    out.push(b'|');
    out.extend_from_slice(kind.as_bytes());
    out.push(b'|');
    for &x in nums {
        let mut u = x as u64;
        loop {
            let b = (u & 0x7f) as u8;
            u >>= 7;
            if u == 0 {
                out.push(b);
                break;
            }
            out.push(b | 0x80);
        }
    }
    out
}

fn l2_get_raw(key: &[u8]) -> Option<Q> {
    let guard = l2().lock().unwrap();
    let conn = guard.as_ref()?;
    let v: String = conn
        .query_row("SELECT v FROM st WHERE k = ?1", [key], |r| r.get(0))
        .ok()?;
    let (n, d) = v.split_once('/')?;
    // 分子分母可能超出 i64：按任意精度有理数解析
    let (n, d): (Q, Q) = (n.parse().ok()?, d.parse().ok()?);
    if d == Q::ZERO {
        return None;
    }
    Some(n / d)
}

fn l2_put_raw(key: Vec<u8>, val: &Q) {
    let mut guard = l2().lock().unwrap();
    let Some(conn) = guard.as_mut() else { return };
    let s = val.to_string();
    let (n, d) = s.split_once('/').unwrap_or((s.as_str(), "1"));
    let _ = conn.execute(
        "INSERT OR REPLACE INTO st (k, v) VALUES (?1, ?2)",
        rusqlite::params![key, format!("{n}/{d}")],
    );
}

/// 通用 L2 读取（跨模块共用同一张表，按 kind 前缀区分）。
pub fn cache_get(kind: &str, nums: &[u32]) -> Option<Q> {
    l2_get_raw(&enc_key(kind, nums))
}

/// 通用 L2 写入。
pub fn cache_put(kind: &str, nums: &[u32], val: &Q) {
    l2_put_raw(enc_key(kind, nums), val);
}

fn fact(n: u64) -> Q {
    let mut r = Q::ONE;
    for i in 2..=n {
        r *= Q::from(i as i64);
    }
    r
}

/// 二项系数 C(n, k)。
fn binom(n: u64, k: u64) -> Q {
    if k > n {
        return Q::ZERO;
    }
    let mut r = Q::ONE;
    for i in 0..k {
        r = r * Q::from(n - i) / Q::from(i + 1);
    }
    r
}

/// Bernoulli B_n（有理数；B_1 = −1/2，与 Sage 一致）。
/// 生成函数递推：Σ_{k=0}^{m} C(m+1,k) B_k = 0（m ≥ 1）。
pub fn bernoulli(n: u64) -> Q {
    let n = n as usize;
    let mut b = vec![Q::ZERO; n + 1];
    b[0] = Q::ONE;
    for m in 1..=n {
        let mut acc = Q::ZERO;
        for (k, bk) in b[..m].iter().enumerate() {
            acc += binom(m as u64 + 1, k as u64) * bk;
        }
        b[m] = -acc / Q::from(m as i64 + 1);
    }
    b[n].clone()
}

/// 有理数幂（指数可为负）。
fn qpow(base: i64, e: i32) -> Q {
    let mut r = Q::ONE;
    let b = Q::from(base);
    let mut i = e.abs();
    while i > 0 {
        r = if e < 0 { r / &b } else { r * &b };
        i -= 1;
    }
    r
}

/// n! 的精确值。
fn fact_q(n: u32) -> Q {
    fact(n as u64)
}

/// 多重集的集合分划（规范形去重），返回 (块表, 轨道因子)：
/// 因子 = |Aut(全体)| / (|Aut(块多重集)| · ∏|Aut(块)|)
/// = 给出该形状的"位置集合分划"个数。对齐 admcycles `DR/utils.setparts_with_auts`。
///
/// 直接枚举多重集分划（块以计数向量表示、按字典序非增排列，每个分划
/// 恰好生成一次），代价与分划个数成正比——不经位置分划（Bell 数）再去重。
pub fn setparts_with_auts(symlist: &[u32]) -> Vec<(Vec<Vec<u32>>, Q)> {
    let mut vals: Vec<u32> = symlist.to_vec();
    vals.sort_unstable();
    let mut distinct: Vec<u32> = Vec::new();
    let mut counts: Vec<u32> = Vec::new();
    for v in vals {
        if distinct.last() == Some(&v) {
            *counts.last_mut().unwrap() += 1;
        } else {
            distinct.push(v);
            counts.push(1);
        }
    }
    let whole_aut: Q = counts.iter().fold(Q::ONE, |acc, &c| acc * fact_q(c));
    let mut out = Vec::new();
    let max_block = counts.clone();
    let mut blocks: Vec<Vec<u32>> = Vec::new();
    rec_multiset_parts(
        &mut counts.clone(),
        &max_block,
        &mut blocks,
        &mut |bs: &[Vec<u32>]| {
            // 块多重集自同构：相同块相邻（非增序）
            let mut denom = Q::ONE;
            let mut run = 1u32;
            for w in bs.windows(2) {
                if w[0] == w[1] {
                    run += 1;
                } else {
                    denom *= fact_q(run);
                    run = 1;
                }
            }
            denom *= fact_q(run);
            for b in bs {
                for &c in b {
                    denom *= fact_q(c);
                }
            }
            let mut shape: Vec<Vec<u32>> = bs
                .iter()
                .map(|b| {
                    b.iter()
                        .zip(&distinct)
                        .flat_map(|(&c, &v)| std::iter::repeat_n(v, c as usize))
                        .collect()
                })
                .collect();
            shape.sort_unstable();
            out.push((shape, &whole_aut / denom));
        },
    );
    out
}

/// 把计数向量 `rem` 拆成若干非零块，块按字典序非增（每块 ≤ `bound`）。
fn rec_multiset_parts(
    rem: &mut Vec<u32>,
    bound: &[u32],
    blocks: &mut Vec<Vec<u32>>,
    emit: &mut dyn FnMut(&[Vec<u32>]),
) {
    if rem.iter().all(|&c| c == 0) {
        emit(blocks);
        return;
    }
    // 逐分量枚举 b ≤ rem（逐分量）、b ≤ bound（字典序）、b ≠ 0
    let k = rem.len();
    let mut b = vec![0u32; k];
    fn walk(
        i: usize,
        tight: bool,
        b: &mut Vec<u32>,
        rem: &mut Vec<u32>,
        bound: &[u32],
        blocks: &mut Vec<Vec<u32>>,
        emit: &mut dyn FnMut(&[Vec<u32>]),
    ) {
        if i == b.len() {
            if b.iter().any(|&c| c > 0) {
                for (r, &x) in rem.iter_mut().zip(b.iter()) {
                    *r -= x;
                }
                blocks.push(b.clone());
                let nb = b.clone();
                rec_multiset_parts(rem, &nb, blocks, emit);
                blocks.pop();
                for (r, &x) in rem.iter_mut().zip(b.iter()) {
                    *r += x;
                }
            }
            return;
        }
        let hi = if tight { rem[i].min(bound[i]) } else { rem[i] };
        for x in (0..=hi).rev() {
            b[i] = x;
            walk(i + 1, tight && x == bound[i], b, rem, bound, blocks, emit);
        }
        b[i] = 0;
    }
    walk(0, true, &mut b, rem, bound, blocks, emit);
}

/// ∫纯 ψ 顶度递归（admcycles `STrecur`）。`psi` 内部会排序。
pub fn st_recur(psi: &[u32]) -> Q {
    let mut p = psi.to_vec();
    p.sort_unstable();
    if let Some(v) = MEMO.with(|m| m.borrow().get(&p).cloned()) {
        return v;
    }
    if let Some(v) = cache_get("st", &p) {
        MEMO.with(|m| m.borrow_mut().insert(p.clone(), v.clone()));
        return v;
    }
    let v = st_recur_calc(&p);
    MEMO.with(|m| m.borrow_mut().insert(p.clone(), v.clone()));
    cache_put("st", &p, &v);
    v
}

fn st_recur_calc(psi: &[u32]) -> Q {
    if psi.is_empty() {
        return Q::ONE;
    }
    let n = psi.len() as i64;
    let s: i64 = psi.iter().map(|&x| x as i64).sum();
    if (s - n) % 3 != 0 {
        return Q::ZERO;
    }
    if psi[0] == 0 {
        // string 递归（与 tau.rs 一致的全亏格 string）
        if s == 0 && n == 3 {
            return Q::ONE;
        }
        let mut total = Q::ZERO;
        let psi_end = &psi[1..];
        for i in 0..(n - 1) as usize {
            if psi_end[i] > 0 {
                let mut c = psi_end.to_vec();
                c[i] -= 1;
                total += st_recur(&c);
            }
        }
        return total;
    }
    let g = (s - n) / 3 + 1;
    let d = psi[psi.len() - 1] as i64;
    let mut total = Q::ZERO;
    // 修正项 1：(2d+3)/12 · STrecur([0,0,0,0]+psi, d+1)
    let mut p1 = vec![0u32; 4];
    p1.extend_from_slice(psi);
    let last = p1.len() - 1;
    p1[last] += 1;
    total += Q::from_signeds(2 * d + 3, 12) * st_recur(&p1);
    // 修正项 2：−(2g+n−1)/6 · STrecur([0,0,0]+psi)
    let mut p2 = vec![0u32; 3];
    p2.extend_from_slice(psi);
    total -= Q::from_signeds(2 * g + n - 1, 6) * st_recur(&p2);
    // 分离分裂项：I ⊆ {0..n−2}。各项只依赖 I 与补集的多重集，故按多重集
    // 枚举：值 v 在前 n−1 个指标中出现 m 次，取其中 c 个进 I，权重 ∏C(m,c)。
    // 代价 ∏(m+1)，而非 2^{n−1}（如 ∫κ₁^{3g−3} 的叶子 [2;3g−3]）。
    let mut vals: Vec<u32> = Vec::new();
    let mut mult: Vec<u32> = Vec::new();
    for &x in &psi[..(n - 1) as usize] {
        if vals.last() == Some(&x) {
            *mult.last_mut().unwrap() += 1;
        } else {
            vals.push(x);
            mult.push(1);
        }
    }
    let mut take = vec![0u32; vals.len()];
    loop {
        let n_i = take.iter().sum::<u32>() as i64 + 2;
        let deg_i: i64 = vals.iter().zip(&take).map(|(&v, &c)| (v * c) as i64).sum();
        if (deg_i - n_i) % 3 == 0 && 2 * ((deg_i - n_i) / 3 + 1) - 2 + n_i > 0 {
            let mut psi3 = vec![0u32; 2];
            let mut psi1 = vec![0u32; 2];
            let mut weight = Q::ONE;
            for ((&v, &m), &c) in vals.iter().zip(&mult).zip(&take) {
                psi3.extend(std::iter::repeat_n(v, c as usize));
                psi1.extend(std::iter::repeat_n(v, (m - c) as usize));
                weight *= binom(m as u64, c as u64);
            }
            let x = st_recur(&psi3);
            // x = 0 时贡献为 0（admcycles 此处 raise），跳过
            if !is_zero(&x) {
                psi1.push(d as u32 + 1);
                let mut psi2 = psi1[1..].to_vec();
                let l2 = psi2.len() - 1;
                psi2[l2] = d as u32;
                total += (Q::from(2 * d + 3) * st_recur(&psi1)
                    - Q::from(2 * g + n - 1) * st_recur(&psi2))
                    * x
                    * weight;
            }
        }
        // 计数向量的里程表式递增
        let mut j = 0;
        while j < take.len() && take[j] == mult[j] {
            take[j] = 0;
            j += 1;
        }
        if j == take.len() {
            break;
        }
        take[j] += 1;
    }
    total /= Q::from((2 * g + n - 1) * (2 * g + n - 2));
    total
}

/// κ→ψ 的分划和封装（admcycles `STsum`）。
pub fn st_sum(psilist: &[u32], kappalist: &[u32]) -> Q {
    let mut kappa = kappalist.to_vec();
    kappa.sort_unstable();
    let mut total = Q::ZERO;
    for (blocks, aut) in setparts_with_auts(&kappa) {
        let mut arg: Vec<u32> = blocks.iter().map(|b| 1 + b.iter().sum::<u32>()).collect();
        arg.extend_from_slice(psilist);
        let sign = if blocks.len() % 2 == 0 { 1 } else { -1 };
        total += Q::from(sign) * aut * st_recur(&arg);
    }
    if kappalist.len() % 2 == 1 {
        -total
    } else {
        total
    }
}

/// 多项系数 multinomial(list) = (Σ)!/∏(各元素!)。
fn multinomial(list: &[u32]) -> Q {
    let s: u64 = list.iter().map(|&x| x as u64).sum();
    let mut r = fact(s);
    for &x in list {
        r /= fact(x as u64);
    }
    r
}

/// CT 公式的求和部分。
fn ct_sum(psilist: &[u32], kappalist: &[u32]) -> Q {
    let mut kappa = kappalist.to_vec();
    kappa.sort_unstable();
    let mut total = Q::ZERO;
    for (blocks, aut) in setparts_with_auts(&kappa) {
        let mut arg: Vec<u32> = blocks.iter().map(|b| 1 + b.iter().sum::<u32>()).collect();
        arg.extend_from_slice(psilist);
        let sign = if blocks.len() % 2 == 0 { 1 } else { -1 };
        total += Q::from(sign) * aut * multinomial(&arg);
    }
    if kappalist.len() % 2 == 1 {
        -total
    } else {
        total
    }
}

/// CT 常数：|(2^{2g−1} − 1) B_{2g}| / (2^{2g−1}·(2g)!)，g=0 时 2^{-1}=1/2。
fn ct_const(g: u32) -> Q {
    let gg = 2 * g as u64;
    let power = qpow(2, gg as i32 - 1); // g=0 → 1/2
    let num = (&power - Q::ONE) * bernoulli(gg);
    (if num < Q::ZERO { -&num } else { num }) / (power * fact(gg))
}

/// ∫_{M̄_{g,n}} Πκ·Πψ（顶度；psilist 可含 0）。维度不符返回 Err。
pub fn socle_formula(g: u32, psilist: &[u32], kappalist: &[u32]) -> Result<Q, String> {
    let degree: i64 = psilist.iter().map(|&x| x as i64).sum::<i64>()
        + kappalist.iter().map(|&x| x as i64).sum::<i64>();
    let n = psilist.len() as i64;
    if degree != 3 * g as i64 - 3 + n {
        return Err(format!(
            "socle 维度不符：g={g}, n={n}, 度={degree} ≠ {}",
            3 * g as i64 - 3 + n
        ));
    }
    if g == 0 {
        return Ok(ct_const(0) * ct_sum(psilist, kappalist));
    }
    Ok(st_sum(psilist, kappalist))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bernoulli_values() {
        // B_0=1, B_1=−1/2, B_2=1/6, B_4=−1/30, B_6=1/42, B_8=−1/30, B_10=5/66
        assert_eq!(bernoulli(0), Q::ONE);
        assert_eq!(bernoulli(1), Q::from_signeds(-1, 2));
        assert_eq!(bernoulli(2), Q::from_signeds(1, 6));
        assert_eq!(bernoulli(4), Q::from_signeds(-1, 30));
        assert_eq!(bernoulli(6), Q::from_signeds(1, 42));
        assert_eq!(bernoulli(8), Q::from_signeds(-1, 30));
        assert_eq!(bernoulli(10), Q::from_signeds(5, 66));
    }

    #[test]
    fn ct_const_sequence() {
        // admcycles 文档序列
        let expect = ["1", "1/24", "7/5760", "31/967680", "127/154828800"];
        for (g, e) in expect.iter().enumerate() {
            assert_eq!(ct_const(g as u32).to_string(), *e, "g={g}");
        }
    }

    #[test]
    fn setparts_examples() {
        // admcycles 文档例
        let r = setparts_with_auts(&[1, 1]);
        assert_eq!(r.len(), 2);
        let r = setparts_with_auts(&[1, 1, 1]);
        // ([[1],[1],[1]],1), ([[1],[1,1]],3), ([[1,1,1]],1)
        assert_eq!(r.len(), 3);
        assert!(r.contains(&(vec![vec![1], vec![1, 1]], Q::from(3))));
        assert!(r.contains(&(vec![vec![1, 1, 1]], Q::ONE)));
        assert_eq!(setparts_with_auts(&[]), vec![(vec![], Q::ONE)]);
    }

    /// 参照实现：枚举全部位置集合分划（Bell 数个）后按形状计数。
    fn setparts_bruteforce(symlist: &[u32]) -> HashMap<Vec<Vec<u32>>, u64> {
        fn rec(
            i: usize,
            s: &[u32],
            blocks: &mut Vec<Vec<u32>>,
            out: &mut HashMap<Vec<Vec<u32>>, u64>,
        ) {
            if i == s.len() {
                let mut shape: Vec<Vec<u32>> = blocks
                    .iter()
                    .map(|b| {
                        let mut v = b.clone();
                        v.sort_unstable();
                        v
                    })
                    .collect();
                shape.sort_unstable();
                *out.entry(shape).or_insert(0) += 1;
                return;
            }
            for bi in 0..=blocks.len() {
                if bi == blocks.len() {
                    blocks.push(vec![s[i]]);
                } else {
                    blocks[bi].push(s[i]);
                }
                rec(i + 1, s, blocks, out);
                if blocks[bi].len() == 1 {
                    blocks.pop();
                } else {
                    blocks[bi].pop();
                }
            }
        }
        let mut out = HashMap::new();
        rec(0, symlist, &mut Vec::new(), &mut out);
        out
    }

    #[test]
    fn setparts_match_bruteforce() {
        for s in [
            vec![1u32, 1, 1, 1],
            vec![1, 1, 2, 2],
            vec![1, 2, 3, 4],
            vec![3, 1, 1, 2, 1, 3],
            vec![1, 1, 1, 1, 2, 2, 3],
            vec![2; 8],
        ] {
            let fast = setparts_with_auts(&s);
            let slow = setparts_bruteforce(&s);
            assert_eq!(fast.len(), slow.len(), "{s:?}");
            for (shape, f) in &fast {
                assert_eq!(*f, Q::from(slow[shape]), "{s:?} {shape:?}");
            }
        }
        // 长 κ 表：分划数 p(20) = 627，旧实现需 Bell(20) ≈ 5×10¹³ 次枚举
        let big = setparts_with_auts(&[1; 20]);
        assert_eq!(big.len(), 627);
        let total = big.iter().fold(Q::ZERO, |acc, (_, f)| acc + f);
        assert_eq!(total.to_string(), "51724158235372"); // Bell(20)
    }

    #[test]
    fn st_recur_doc_values() {
        assert_eq!(st_recur(&[0, 0, 0]), Q::ONE);
        assert_eq!(st_recur(&[1]), Q::from_signeds(1, 24));
        assert_eq!(st_recur(&[4]), Q::from_signeds(1, 1152)); // M3 真值
        assert_eq!(st_recur(&[7]), Q::from_signeds(1, 82944));
        assert_eq!(st_recur(&[2, 2, 2]), Q::from_signeds(7, 240));
        assert_eq!(st_recur(&[2; 6]), Q::from_signeds(1225, 144));
        assert_eq!(st_recur(&[1, 1, 1]), Q::from_signeds(1, 12));
    }

    /// 回归：分裂项按多重集枚举（重复指数带二项式权重）须与 dilaton 闭式一致。
    #[test]
    fn st_recur_repeated_entries() {
        // ⟨τ₁ⁿ⟩₁ = (n−1)!/24
        for n in 1..=16usize {
            assert_eq!(
                st_recur(&vec![1; n]),
                fact(n as u64 - 1) / Q::from(24),
                "n={n}"
            );
        }
        // ⟨τ₁ᵏτ₄⟩₂ = (k+2)!/2 · ⟨τ₄⟩₂ = (k+2)!/2304
        for k in 0..=16usize {
            let mut p = vec![1; k];
            p.push(4);
            assert_eq!(st_recur(&p), fact(k as u64 + 2) / Q::from(2304), "k={k}");
        }
    }

    #[test]
    fn socle_kappa_only_g3() {
        // admcycles 文档（MODULI_ST）
        assert_eq!(
            socle_formula(3, &[], &[6]).unwrap(),
            Q::from_signeds(1, 82944)
        );
        assert_eq!(
            socle_formula(3, &[], &[1, 5]).unwrap(),
            Q::from_signeds(1, 5760)
        );
        assert_eq!(
            socle_formula(3, &[], &[2, 4]).unwrap(),
            Q::from_signeds(971, 2903040)
        );
        assert_eq!(
            socle_formula(3, &[], &[1, 1, 4]).unwrap(),
            Q::from_signeds(2173, 967680)
        );
        assert_eq!(
            socle_formula(3, &[], &[1, 1, 1, 1, 1, 1]).unwrap(),
            Q::from_signeds(176557, 107520)
        );
    }

    /// κ₁ 最高幂（Weil–Petersson 体积）：长 κ 表曾因 Bell 数枚举不可算
    /// （g=6 需 Bell(15) ≈ 1.4×10⁹ 个位置分划）。
    #[test]
    fn kappa1_top_powers() {
        assert_eq!(
            socle_formula(4, &[], &[1; 9]).unwrap().to_string(),
            "1959225867017/696729600"
        );
        assert_eq!(
            socle_formula(6, &[], &[1; 15]).unwrap().to_string(),
            "2516292682076619940682627/2522559283200"
        );
    }

    #[test]
    fn socle_g0_kontsevich() {
        // 与 Kontsevich 一致：(n−3)!/∏d!
        assert!(socle_formula(0, &[1, 1, 0, 0], &[]).is_err()); // Σ=2≠1
        assert_eq!(
            socle_formula(0, &[2, 1, 0, 0, 0, 0], &[]).unwrap(),
            Q::from(3)
        );
        // κ_1 on M̄_{0,4}（顶度 1）：CT 路线 = 多项系数
        assert_eq!(socle_formula(0, &[0, 0, 0, 0], &[1]).unwrap(), Q::ONE);
    }

    #[test]
    fn m3_targets() {
        // FABER_REDESIGN 真值表
        assert_eq!(
            socle_formula(2, &[4], &[]).unwrap(),
            Q::from_signeds(1, 1152)
        ); // ∫ψ₁⁴ M̄_{2,1}
        assert_eq!(
            socle_formula(2, &[], &[1, 1, 1]).unwrap(),
            Q::from_signeds(43, 2880)
        ); // ∫κ₁³ M̄_2
    }
}
