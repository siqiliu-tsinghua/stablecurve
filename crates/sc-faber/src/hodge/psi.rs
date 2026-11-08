//! Mumford GRR 的后代形式（Faber–Pandharipande, Proposition 2）：
//! https://arxiv.org/abs/math/9810173
//!
//! 令 C_m = (m+1)! ch_m / B_{m+1}（m 奇数），则插入一个 C_m 等于：
//! 新增 τ_{m+1} − Σ_i 把 τ_{d_i} 换为 τ_{d_i+m}，加上两类边界项
//! ½Σ_{a+b=m−1}(−1)^a（非分离项 + 分离项）。其余 C 在分离边界上
//! 按二项式分配。这一递归始终只保留 ψ 和 C，不生成 κ 多项式。
//!
//! Bernoulli 因子在入口一次性提出，避免递归中反复乘大分母。相同 ψ 次数
//! 和 C 指标按重数分配；string/dilaton 先消去 0/1 次标记。

use super::{ChMono, ch_mul, fact, lambda_expansion};
use crate::cache_source::Registry;
use crate::socle::bernoulli;
use sc_num::{One, Q, Zero};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::{Arc, Weak};

#[cfg(test)]
mod measure_serial;
mod parallel;
mod serial_source;

use serial_source::{SerialSource, SerialValues};

#[derive(Default)]
struct Descendents {
    // ch 与 ψ 的总次数已决定亏格。分两层存储，使命中时可直接借用 ch 切片。
    values: SerialValues,
    source: Option<Arc<SerialSource>>,
    #[cfg(test)]
    shared_values: Option<std::sync::Arc<std::sync::RwLock<measure_serial::SerialValues>>>,
    grr: Vec<Q>,
    factorials: Vec<Q>,
    odd_factorials: Vec<Q>,
    binomials: Vec<Vec<Q>>,
    one_point: Vec<Q>,
}

thread_local! {
    static DESCENDENTS: RefCell<Descendents> = RefCell::new(Descendents::default());
    static CACHED_DESCENDENTS: RefCell<BTreeMap<u64, CachedDescendents>> = const { RefCell::new(BTreeMap::new()) };
}

struct CachedDescendents {
    registry: Weak<Registry>,
    evaluator: Descendents,
}

pub(super) fn clear_cache_context(id: u64) {
    CACHED_DESCENDENTS.with(|cell| {
        cell.borrow_mut().remove(&id);
    });
    parallel::clear_cache_context(id);
}

fn with_descendents(
    registry: Option<&Arc<Registry>>,
    calculate: impl FnOnce(&mut Descendents) -> Q,
) -> Q {
    let Some(registry) = registry else {
        return DESCENDENTS.with(|cell| calculate(&mut cell.borrow_mut()));
    };
    CACHED_DESCENDENTS.with(|cell| {
        let mut contexts = cell.borrow_mut();
        if let Some(context) = contexts.get_mut(&registry.id()) {
            return calculate(&mut context.evaluator);
        }
        // Do not scan every historical memo or evict another live instance's
        // values when a query switches its cache attachment.
        contexts.retain(|_, context| context.registry.strong_count() != 0);
        let context = contexts
            .entry(registry.id())
            .or_insert_with(|| CachedDescendents {
                registry: Arc::downgrade(registry),
                evaluator: Descendents::with_registry(Some(registry)),
            });
        calculate(&mut context.evaluator)
    })
}

fn multinomial(d: &[u32]) -> Q {
    let mut value = fact(d.iter().map(|&x| u64::from(x)).sum());
    for &x in d {
        value /= fact(u64::from(x));
    }
    value
}

pub(super) fn integrate(g: u32, n: u32, psi: &[(u32, u32)], lambda: &[(u32, u32)]) -> Q {
    let threads = thread_count(g, if lambda.is_empty() { 16 } else { 24 });
    integrate_with_threads(g, n, psi, lambda, threads)
}

fn thread_count(g: u32, maximum: usize) -> usize {
    sc_num::parallel::thread_count("SC_FABER_THREADS", if g < 7 { 1 } else { maximum })
}

fn integrate_with_threads(
    g: u32,
    n: u32,
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    threads: usize,
) -> Q {
    integrate_ch_with_threads(g, n, psi, lambda, &[], threads)
}

fn integrate_ch_with_threads(
    g: u32,
    n: u32,
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    explicit_ch: &[(u32, u32)],
    threads: usize,
) -> Q {
    let weight = lambda
        .iter()
        .chain(explicit_ch)
        .map(|&(l, p)| u64::from(l) * u64::from(p))
        .sum::<u64>();
    let psi_weight = psi.iter().map(|&(_, p)| u64::from(p)).sum::<u64>();
    if (weight + psi_weight) as i64 != 3 * i64::from(g) - 3 + i64::from(n) {
        return Q::ZERO;
    }
    let mut d = vec![0; n as usize];
    for &(leg, power) in psi {
        d[leg as usize - 1] = power;
    }
    d.sort_unstable();

    integrate_ch_roots(g, &[(d, Q::ONE)], lambda, explicit_ch, threads)
}

fn integrate_ch_roots(
    g: u32,
    roots: &[(Vec<u32>, Q)],
    lambda: &[(u32, u32)],
    explicit_ch: &[(u32, u32)],
    threads: usize,
) -> Q {
    let weight: u64 = lambda
        .iter()
        .chain(explicit_ch)
        .map(|&(l, p)| u64::from(l) * u64::from(p))
        .sum();
    let factors: Vec<_> = lambda
        .iter()
        .copied()
        .filter(|&(l, p)| l > 0 && p > 0)
        .collect();
    // Mumford 关系 c(E)c(E*)=1 的最高次项：λ_g²=0。
    if g > 0
        && factors
            .iter()
            .filter(|&&(l, _)| l == g)
            .map(|&(_, p)| u64::from(p))
            .sum::<u64>()
            >= 2
    {
        return Q::ZERO;
    }
    // λ_g 定理：∫ λ_g Πψ_i^{d_i} = multinomial(d) · (1−2^{1−2g})|B_{2g}|/(2g)!.
    // https://arxiv.org/abs/math/9908052, Theorem 1 and equation (5).
    if g > 0 && factors == [(g, 1)] && explicit_ch.is_empty() {
        let mut power = Q::ONE;
        for _ in 0..2 * g - 1 {
            power *= Q::from(2);
        }
        let b = bernoulli(2 * u64::from(g));
        let b = if b < Q::ZERO { -b } else { b };
        return roots
            .iter()
            .fold(Q::ZERO, |sum, (d, c)| sum + c * multinomial(d))
            * (Q::ONE - Q::ONE / power)
            * b
            / fact(2 * u64::from(g));
    }
    // Hodge 类从 M_g（g≥2）或 M_{1,1} 拉回；超过其维数时为零。
    if weight
        > if g > 1 {
            3 * u64::from(g) - 3
        } else {
            u64::from(g)
        }
    {
        return Q::ZERO;
    }
    let mut poly = vec![(explicit_ch.to_vec(), Q::ONE)];
    for (l, p) in factors {
        let expansion = lambda_expansion(l);
        for _ in 0..p {
            poly = collect_terms(ch_mul(&poly, &expansion));
        }
    }
    crate::cache_source::with_current_registry(|registry| {
        with_descendents(registry, |evaluator| {
            for (ch, coefficient) in &mut poly {
                for &(m, p) in ch.iter() {
                    let factor = evaluator.grr_factor(m);
                    for _ in 0..p {
                        *coefficient *= &factor;
                    }
                }
            }
            if threads > 1 {
                return parallel::integrate(g, roots, &poly, evaluator, threads, registry);
            }
            if registry
                .is_some_and(|registry| registry.has_history() && registry.check_error().is_ok())
            {
                let mut sum = Q::ZERO;
                for (d, c) in roots {
                    for (ch, coefficient) in &poly {
                        let value = if let Some(value) = evaluator.memoized_root(g, d, ch) {
                            value
                        } else if direct_root(g, d, ch) {
                            evaluator.eval(g, d.clone(), ch)
                        } else {
                            // Reuse the existing layered planner only for the
                            // historical L2 miss path. One evaluation thread
                            // still permits frontier-sized storage requests.
                            return parallel::integrate(g, roots, &poly, evaluator, 1, registry);
                        };
                        sum += c * coefficient * value;
                    }
                }
                return sum;
            }
            let has_parallel_cache =
                registry.is_some_and(|registry| parallel::has_cached_engine(registry));
            if !has_parallel_cache {
                let mut sum = Q::ZERO;
                for (d, c) in roots {
                    for (ch, coefficient) in &poly {
                        let value = evaluator
                            .memoized_root(g, d, ch)
                            .unwrap_or_else(|| evaluator.eval(g, d.clone(), ch));
                        sum += c * coefficient * value;
                    }
                }
                return sum;
            }
            let mut sum = Q::ZERO;
            for (d, c) in roots {
                for (ch, coefficient) in &poly {
                    if let Some(value) = evaluator.memoized_root(g, d, ch) {
                        sum += c * coefficient * value;
                        continue;
                    }
                    let mut descendants = d.clone();
                    let cached = registry.and_then(|registry| {
                        descendants.sort_unstable();
                        let n = descendants.len() as i64;
                        let weight = ch
                            .iter()
                            .map(|&(m, p)| i64::from(m) * i64::from(p))
                            .sum::<i64>();
                        if 2 * i64::from(g) - 2 + n > 0
                            && weight
                                + descendants
                                    .iter()
                                    .map(|&power| i64::from(power))
                                    .sum::<i64>()
                                == 3 * i64::from(g) - 3 + n
                        {
                            parallel::cached(g, &descendants, ch, registry)
                        } else {
                            None
                        }
                    });
                    let value = cached.unwrap_or_else(|| evaluator.eval(g, descendants, ch));
                    sum += c * coefficient * value;
                }
            }
            sum
        })
    })
}

pub(super) fn pure(g: u32, d: &[u32]) -> Q {
    if thread_count(g, 16) == 1 {
        crate::cache_source::with_current_registry(|registry| {
            with_descendents(registry, |evaluator| {
                if let Some(value) = evaluator.memoized_root(g, d, &[]) {
                    return value;
                }
                let mut descendants = d.to_vec();
                if let Some(registry) =
                    registry.filter(|registry| parallel::has_cached_engine(registry))
                {
                    descendants.sort_unstable();
                    let n = descendants.len() as i64;
                    if 2 * i64::from(g) - 2 + n > 0
                        && descendants
                            .iter()
                            .map(|&power| i64::from(power))
                            .sum::<i64>()
                            == 3 * i64::from(g) - 3 + n
                    {
                        if let Some(value) = parallel::cached(g, &descendants, &[], registry) {
                            return value;
                        }
                    }
                }
                if registry.is_some_and(|registry| {
                    registry.has_history() && registry.check_error().is_ok()
                }) && !direct_root(g, &descendants, &[])
                {
                    return parallel::integrate(
                        g,
                        &[(descendants, Q::ONE)],
                        &[(vec![], Q::ONE)],
                        evaluator,
                        1,
                        registry,
                    );
                }
                evaluator.eval(g, descendants, &[])
            })
        })
    } else {
        pure_sum(g, &[(d.to_vec(), Q::ONE)])
    }
}

pub(super) fn pure_sum(g: u32, roots: &[(Vec<u32>, Q)]) -> Q {
    integrate_ch_roots(g, roots, &[], &[], thread_count(g, 16))
}

/// κ 的逆推前公式：对位置集合分划求和，每块 B 新增 τ_{1+Σ_{i∈B}κ_i}，
/// 系数为 (−1)^{κ 的个数−块数}。Hodge 类随忘点映射拉回。
pub(super) fn integrate_mixed(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Q {
    integrate_mixed_with_threads(g, n, kappa, psi, lambda, ch, thread_count(g, 16))
}

fn integrate_mixed_with_threads(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
    threads: usize,
) -> Q {
    if kappa.is_empty() {
        return integrate_ch_with_threads(g, n, psi, lambda, ch, threads);
    }
    let mut kappalist = Vec::new();
    for &(index, power) in kappa {
        kappalist.extend(std::iter::repeat_n(index, power as usize));
    }
    // 不同分划可产生相同的后代幂列表；先合并再递归。
    let mut roots = BTreeMap::<Vec<u32>, Q>::new();
    for (blocks, multiplicity) in crate::socle::setparts_with_auts(&kappalist) {
        let mut powers: Vec<u32> = blocks.iter().map(|b| 1 + b.iter().sum::<u32>()).collect();
        powers.sort_unstable();
        let coefficient = if (kappalist.len() - blocks.len()) % 2 == 0 {
            multiplicity
        } else {
            -multiplicity
        };
        *roots.entry(powers).or_insert(Q::ZERO) += coefficient;
    }
    let mut descendants = vec![0; n as usize];
    for &(leg, power) in psi {
        descendants[leg as usize - 1] = power;
    }
    let roots: Vec<_> = roots
        .into_iter()
        .filter(|(_, c)| *c != Q::ZERO)
        .map(|(powers, c)| {
            let mut d = descendants.clone();
            d.extend(powers);
            d.sort_unstable();
            (d, c)
        })
        .collect();
    integrate_ch_roots(g, &roots, lambda, ch, threads)
}

fn collect_terms(mut terms: Vec<ChMono>) -> Vec<ChMono> {
    terms.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut out: Vec<ChMono> = Vec::with_capacity(terms.len());
    for (mono, coefficient) in terms {
        if let Some(last) = out.last_mut().filter(|last| last.0 == mono) {
            last.1 += coefficient;
        } else {
            out.push((mono, coefficient));
        }
    }
    out.retain(|(_, c)| *c != Q::ZERO);
    out
}

// These are the original inexpensive eval guards and closed formulas. They
// keep their local path instead of waiting for storage during recovery.
fn direct_root(g: u32, d: &[u32], ch: &[(u32, u32)]) -> bool {
    let n = d.len() as i64;
    let weight: u64 = ch.iter().map(|&(m, p)| u64::from(m) * u64::from(p)).sum();
    if 2 * i64::from(g) - 2 + n <= 0
        || weight as i128 + d.iter().map(|&power| i128::from(power)).sum::<i128>()
            != 3 * i128::from(g) - 3 + i128::from(n)
    {
        return true;
    }
    if ch.is_empty() {
        g == 0 || d.len() == 1
    } else {
        g == 0 || ch.last().unwrap().0 >= 2 * g || weight > u64::from((3 * g - 3).max(1)) || g == 1
    }
}

impl Descendents {
    pub(super) fn with_registry(registry: Option<&Arc<Registry>>) -> Self {
        let source = registry.map(|registry| {
            let source = Arc::new(SerialSource::default());
            registry.register(source.clone());
            source
        });
        Self {
            source,
            ..Self::default()
        }
    }

    pub(super) fn with_memo<R>(&self, read: impl FnOnce(&SerialValues) -> R) -> R {
        if let Some(source) = &self.source {
            return read(source.read().values());
        }
        #[cfg(test)]
        if let Some(shared) = &self.shared_values {
            return read(&shared.read().unwrap());
        }
        read(&self.values)
    }

    fn memoized_root(&self, g: u32, d: &[u32], ch: &[(u32, u32)]) -> Option<Q> {
        let source = self.source.as_ref()?;
        let n = d.len() as i64;
        let weight: u32 = ch.iter().map(|&(m, p)| m * p).sum();
        // The original memo omits g because a valid top-degree state fixes it.
        // Check that invariant before borrowing a caller's key without eval.
        if !d.is_sorted()
            || 2 * i64::from(g) - 2 + n <= 0
            || i64::from(weight) + d.iter().map(|&power| i64::from(power)).sum::<i64>()
                != 3 * i64::from(g) - 3 + n
        {
            return None;
        }
        source.get(d, ch)
    }

    fn finish_formula(&self, g: u32, d: Vec<u32>, ch: &[(u32, u32)], value: Q) -> Q {
        if let Some(source) = &self.source {
            source.insert(g, d, ch, &value);
        }
        value
    }

    fn factorial(&mut self, n: usize) -> Q {
        while self.factorials.len() <= n {
            let value = self
                .factorials
                .last()
                .map_or(Q::ONE, |previous| previous * Q::from(self.factorials.len()));
            self.factorials.push(value);
        }
        self.factorials[n].clone()
    }

    fn multinomial(&mut self, d: &[u32]) -> Q {
        let mut value = self.factorial(d.iter().map(|&x| x as usize).sum());
        for &x in d {
            value /= self.factorial(x as usize);
        }
        value
    }

    /// (2n+1)!!，供 DVV 结点权重和左侧系数复用。
    fn odd_factorial(&mut self, n: u32) -> Q {
        while self.odd_factorials.len() <= n as usize {
            let value = self.odd_factorials.last().map_or(Q::ONE, |previous| {
                previous * Q::from(2 * self.odd_factorials.len() + 1)
            });
            self.odd_factorials.push(value);
        }
        self.odd_factorials[n as usize].clone()
    }

    fn binomial(&mut self, n: u32, k: u32) -> Q {
        while self.binomials.len() <= n as usize {
            let mut row = vec![Q::ONE; self.binomials.len() + 1];
            if let Some(previous) = self.binomials.last() {
                for (value, pair) in row.iter_mut().skip(1).zip(previous.windows(2)) {
                    *value = &pair[0] + &pair[1];
                }
            }
            self.binomials.push(row);
        }
        self.binomials[n as usize][k as usize].clone()
    }

    fn grr_factor(&mut self, m: u32) -> Q {
        while self.grr.len() <= (m / 2) as usize {
            let k = 2 * self.grr.len() as u64 + 2;
            self.grr.push(bernoulli(k) / fact(k));
        }
        self.grr[(m / 2) as usize].clone()
    }

    fn eval(&mut self, g: u32, mut d: Vec<u32>, ch: &[(u32, u32)]) -> Q {
        let n = d.len() as i64;
        let weight: u32 = ch.iter().map(|&(m, p)| m * p).sum();
        if 2 * i64::from(g) - 2 + n <= 0
            || i64::from(weight) + d.iter().map(|&x| i64::from(x)).sum::<i64>()
                != 3 * i64::from(g) - 3 + n
        {
            return Q::ZERO;
        }
        d.sort_unstable();
        if let Some(source) = &self.source {
            if let Some(value) = source.get(&d, ch) {
                return value;
            }
        }
        if ch.is_empty() {
            if g == 0 {
                let value = self.multinomial(&d);
                return self.finish_formula(g, d, ch, value);
            }
            if d.len() == 1 {
                while self.one_point.len() <= g as usize {
                    let value = self.one_point.last().map_or(Q::ONE, |previous| {
                        previous / Q::from(24 * self.one_point.len())
                    });
                    self.one_point.push(value);
                }
                return self.finish_formula(g, d, ch, self.one_point[g as usize].clone());
            }
        } else if g == 0 || ch.last().unwrap().0 >= 2 * g || weight > (3 * g - 3).max(1) {
            return self.finish_formula(g, d, ch, Q::ZERO);
        } else if g == 1 {
            // C₁ = 12λ₁；∫λ₁Πψ^d = multinomial(d)/24。
            let value = self.multinomial(&d) / Q::from(2);
            return self.finish_formula(g, d, ch, value);
        }
        if self.source.is_none() {
            #[cfg(test)]
            if let Some(value) = self.measurement_memo_get(&d, ch) {
                return value;
            }
            #[cfg(not(test))]
            if let Some(value) = self.values.get(ch).and_then(|values| values.get(&d)) {
                return value.clone();
            }
        }
        let value = self.compute(g, &d, ch);
        if let Some(source) = &self.source {
            source.insert(g, d, ch, &value);
        } else {
            #[cfg(test)]
            self.measurement_memo_insert(d, ch, &value);
            #[cfg(not(test))]
            self.values
                .entry(ch.to_vec())
                .or_default()
                .insert(d, value.clone());
        }
        value
    }

    fn compute(&mut self, g: u32, d: &[u32], ch: &[(u32, u32)]) -> Q {
        if d.first() == Some(&0) {
            let mut sum = Q::ZERO;
            let mut i = 1;
            while i < d.len() {
                let end = i + d[i..].partition_point(|&x| x == d[i]);
                if d[i] > 0 {
                    let mut next = d[1..].to_vec();
                    next[i - 1] -= 1;
                    sum += Q::from(end - i) * self.eval(g, next, ch);
                }
                i = end;
            }
            return sum;
        }
        if d.first() == Some(&1) {
            // g=1 的纯 ψ 最终保留 ⟨τ₁⟩₁=1/24，不能忘掉最后一个标记。
            let count =
                d.partition_point(|&x| x == 1)
                    .min(if g == 1 { d.len() - 1 } else { d.len() });
            let mut factor = Q::ONE;
            for n in d.len() - count..d.len() {
                factor *= Q::from(2 * u64::from(g) - 2 + n as u64);
            }
            return factor * self.eval(g, d[count..].to_vec(), ch);
        }

        if ch.is_empty() {
            return self.pure_psi(g, d);
        }

        let mut rest = ch.to_vec();
        let last = rest.last_mut().unwrap();
        let m = last.0;
        last.1 -= 1;
        if last.1 == 0 {
            rest.pop();
        }
        let mut added = d.to_vec();
        added.push(m + 1);
        let mut sum = self.eval(g, added, &rest);

        let mut groups = Vec::new();
        let mut i = 0;
        while i < d.len() {
            let end = i + d[i..].partition_point(|&x| x == d[i]);
            let mut shifted = d.to_vec();
            shifted[i] += m;
            sum -= Q::from(end - i) * self.eval(g, shifted, &rest);
            groups.push((d[i], (end - i) as u32));
            i = end;
        }
        let weights: Vec<Q> = (0..m)
            .map(|a| Q::from(if a % 2 == 0 { 1 } else { -1 }))
            .collect();
        sum + self.boundary(g, d, &rest, &groups, &weights)
    }

    /// DVV：选最小 ψ 次数，结点分裂只需 d_min−1 项。
    /// https://msp.org/pjm/2009/241-2/pjm-v241-n2-p05-p.pdf, equation (4).
    /// (2k+1)!!⟨τ_k Πτ_d⟩ = Σ_i (2k+2d_i−1)!!/(2d_i−1)!!⟨τ_{k+d_i−1}…⟩
    /// + ½Σ_{a+b=k−2}(2a+1)!!(2b+1)!!（两类边界项）。
    fn pure_psi(&mut self, g: u32, d: &[u32]) -> Q {
        let k = d[0];
        let rest = &d[1..];
        let mut groups = Vec::new();
        let mut sum = Q::ZERO;
        let mut i = 0;
        while i < rest.len() {
            let end = i + rest[i..].partition_point(|&x| x == rest[i]);
            let mut merged = rest.to_vec();
            merged[i] += k - 1;
            let mut factor = Q::from(end - i);
            for j in 1..=k {
                factor *= Q::from(2 * rest[i] + 2 * j - 1);
            }
            sum += factor * self.eval(g, merged, &[]);
            groups.push((rest[i], (end - i) as u32));
            i = end;
        }
        let weights: Vec<Q> = (0..k - 1)
            .map(|a| self.odd_factorial(a) * self.odd_factorial(k - 2 - a))
            .collect();
        (sum + self.boundary(g, rest, &[], &groups, &weights)) / self.odd_factorial(k)
    }

    fn boundary(
        &mut self,
        g: u32,
        d: &[u32],
        ch: &[(u32, u32)],
        groups: &[(u32, u32)],
        weights: &[Q],
    ) -> Q {
        let m = weights.len() as u32;
        let mut sum = Q::ZERO;
        // 结点权重关于两侧对称，仅中间项保留 1/2。
        for a in 0..m.div_ceil(2) {
            let mut node = d.to_vec();
            node.extend([a, m - 1 - a]);
            let mut value = &weights[a as usize] * self.eval(g - 1, node, ch);
            if 2 * a == m - 1 {
                value /= Q::from(2);
            }
            sum += value;
        }

        let mut split = Split {
            evaluator: self,
            g,
            m,
            groups,
            ch,
            weights,
            left_d: Vec::new(),
            right_d: Vec::new(),
            left_ch: Vec::new(),
            right_ch: Vec::new(),
        };
        sum += split.distribute(0, 0, Q::ONE);
        sum
    }
}

struct Split<'a> {
    evaluator: &'a mut Descendents,
    g: u32,
    m: u32,
    groups: &'a [(u32, u32)],
    ch: &'a [(u32, u32)],
    weights: &'a [Q],
    left_d: Vec<u32>,
    right_d: Vec<u32>,
    left_ch: Vec<(u32, u32)>,
    right_ch: Vec<(u32, u32)>,
}

impl Split<'_> {
    fn distribute(&mut self, index: usize, left_defect: u32, coefficient: Q) -> Q {
        if index == self.groups.len() + self.ch.len() {
            return self.evaluate(left_defect, coefficient);
        }
        let is_psi = index < self.groups.len();
        let (degree, count) = if is_psi {
            self.groups[index]
        } else {
            self.ch[index - self.groups.len()]
        };
        let left_len = self.left_d.len();
        let right_len = self.right_d.len();
        let mut sum = Q::ZERO;
        // defect = Σ(d_i−1) + Σm·p_m；非负结点次数要求 defect≤3h−2。
        // 只枚举 h≤g/2，超过上界的整个分配子树都没有贡献。
        let degree = if is_psi { degree - 1 } else { degree };
        let bound = 3 * (self.g / 2) - 2;
        let max_taken = count.min((bound - left_defect) / degree);
        for taken in 0..=max_taken {
            if is_psi {
                self.left_d.resize(left_len + taken as usize, degree + 1);
                self.right_d
                    .resize(right_len + (count - taken) as usize, degree + 1);
            } else {
                if taken > 0 {
                    self.left_ch.push((degree, taken));
                }
                if taken < count {
                    self.right_ch.push((degree, count - taken));
                }
            }
            let weight = &coefficient * self.evaluator.binomial(count, taken);
            sum += self.distribute(index + 1, left_defect + degree * taken, weight);
            if !is_psi {
                if taken > 0 {
                    self.left_ch.pop();
                }
                if taken < count {
                    self.right_ch.pop();
                }
            }
        }
        self.left_d.truncate(left_len);
        self.right_d.truncate(right_len);
        sum
    }

    fn evaluate(&mut self, left_defect: u32, coefficient: Q) -> Q {
        let mut sum = Q::ZERO;
        // 所有原标记次数≥2，故亏格 0 的分离分支不可能顶维。
        // 维数直接确定结点次数 a，无须另枚举 a=0..m−1。
        let min_h = ((left_defect + 4) / 3).max(1);
        let max_h = ((left_defect + self.m + 1) / 3).min(self.g / 2);
        let left_ch_weight: u32 = self.left_ch.iter().map(|&(m, p)| m * p).sum();
        let right_ch_weight: u32 = self.right_ch.iter().map(|&(m, p)| m * p).sum();
        for h in min_h..=max_h {
            if left_ch_weight > (3 * h - 3).max(1)
                || right_ch_weight > (3 * (self.g - h) - 3).max(1)
                || self.left_ch.last().is_some_and(|&(m, _)| m >= 2 * h)
                || self
                    .right_ch
                    .last()
                    .is_some_and(|&(m, _)| m >= 2 * (self.g - h))
            {
                continue;
            }
            let a = 3 * h - 2 - left_defect;
            let b = self.m - 1 - a;
            let mut diagonal = false;
            if 2 * h == self.g {
                let order = self
                    .left_d
                    .cmp(&self.right_d)
                    .then(self.left_ch.cmp(&self.right_ch))
                    .then(a.cmp(&b));
                if order.is_gt() {
                    continue;
                }
                diagonal = order.is_eq();
            }
            let mut left = self.left_d.clone();
            left.push(a);
            let lv = self.evaluator.eval(h, left, &self.left_ch);
            if lv == Q::ZERO {
                continue;
            }
            let mut right = self.right_d.clone();
            right.push(b);
            let mut value =
                &coefficient * lv * self.evaluator.eval(self.g - h, right, &self.right_ch);
            if diagonal {
                value /= Q::from(2);
            }
            sum += &self.weights[a as usize] * value;
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hodge::{int_class, klp};

    #[test]
    fn borrowed_root_memo_validates_genus_dimension_and_c_semantics() {
        let registry = Arc::new(Registry::new());
        let mut evaluator = Descendents::with_registry(Some(&registry));
        let expected = Q::from_signeds(1, 384);
        assert_eq!(evaluator.eval(2, vec![1, 4], &[]), expected);
        assert_eq!(evaluator.memoized_root(2, &[1, 4], &[]), Some(expected));
        assert_eq!(evaluator.memoized_root(3, &[1, 4], &[]), None);
        assert_eq!(evaluator.memoized_root(2, &[1, 3], &[]), None);
        assert_eq!(evaluator.memoized_root(0, &[], &[]), None);
        assert_eq!(evaluator.memoized_root(2, &[4, 1], &[]), None);

        // The descendant source holds C_1 = 12 lambda_1, without the outer
        // Bernoulli factor used to recover true ch/lambda values.
        assert_eq!(evaluator.eval(1, vec![0], &[(1, 1)]), Q::from_signeds(1, 2));
        assert_eq!(
            evaluator.memoized_root(1, &[0], &[(1, 1)]),
            Some(Q::from_signeds(1, 2)),
        );
        assert_eq!(evaluator.memoized_root(1, &[0], &[]), None);
        assert_eq!(evaluator.eval(0, vec![0; 4], &[(1, 1)]), Q::ZERO);
        assert_eq!(
            evaluator.memoized_root(0, &[0; 4], &[(1, 1)]),
            Some(Q::ZERO),
        );

        let mut original = Descendents::default();
        assert_eq!(original.eval(2, vec![1, 4], &[]), Q::from_signeds(1, 384));
        assert_eq!(original.memoized_root(2, &[1, 4], &[]), None);
    }

    #[test]
    fn cached_serial_roots_preserve_unsorted_aliases_and_reject_wrong_genus() {
        use crate::cache_source::with_registry;

        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            let calculate = |g, d| integrate_ch_roots(g, &[(d, Q::ONE)], &[], &[], 1);
            let expected = Q::from_signeds(1, 384);
            assert_eq!(calculate(2, vec![1, 4]), expected);
            let source = registry.sources().pop().unwrap();
            let ready = source.ready();
            assert_eq!(calculate(2, vec![1, 4]), expected);
            assert_eq!(calculate(2, vec![4, 1]), expected);
            assert_eq!(calculate(3, vec![1, 4]), Q::ZERO);
            assert_eq!(calculate(0, vec![]), Q::ZERO);
            assert_eq!(source.ready(), ready);
        });
    }

    fn legacy(g: u32, d: &[u32], lambda: &[(u32, u32)]) -> Q {
        let psi: Vec<_> = d
            .iter()
            .enumerate()
            .map(|(i, &p)| (i as u32 + 1, p))
            .collect();
        let mut poly = vec![(vec![], Q::ONE)];
        for &(l, p) in lambda {
            let expansion = lambda_expansion(l);
            for _ in 0..p {
                poly = ch_mul(&poly, &expansion);
            }
        }
        poly.into_iter().fold(Q::ZERO, |sum, (ch, c)| {
            sum + c * int_class(g, d.len() as u32, &[], &psi, &ch).unwrap()
        })
    }

    fn direct(g: u32, d: &[u32], lambda: &[(u32, u32)]) -> Q {
        let psi: Vec<_> = d
            .iter()
            .enumerate()
            .map(|(i, &p)| (i as u32 + 1, p))
            .collect();
        klp(g, d.len() as u32, &[], &psi, lambda, &[]).unwrap()
    }

    fn partitions(
        sum: u32,
        min: u32,
        slots: usize,
        prefix: &mut Vec<u32>,
        emit: &mut dyn FnMut(&[u32]),
    ) {
        if slots == 0 {
            if sum == 0 {
                emit(prefix);
            }
            return;
        }
        for value in min..=sum / slots as u32 {
            prefix.push(value);
            partitions(sum - value, value, slots - 1, prefix, emit);
            prefix.pop();
        }
    }

    #[test]
    fn descendants_match_kappa_recursion() {
        // 包含 τ₀/τ₁、相同次数的标记、多个奇 ch、λ_g 闭式和 Hodge 过维消没。
        for g in 1..=4 {
            let mut lambdas = vec![vec![], vec![(1, 2)], vec![(1, 3)], vec![(1, 5)]];
            for l in 1..=g {
                lambdas.push(vec![(l, 1)]);
                lambdas.push(vec![(1, 1), (l, 1)]);
                lambdas.push(vec![(l, 2)]);
            }
            for n in 1..=4 {
                for lambda in &lambdas {
                    let w: u32 = lambda.iter().map(|&(l, p)| l * p).sum();
                    let dimension = 3 * g - 3 + n;
                    if w > dimension {
                        continue;
                    }
                    partitions(dimension - w, 0, n as usize, &mut vec![], &mut |d| {
                        assert_eq!(
                            direct(g, d, lambda),
                            legacy(g, d, lambda),
                            "g={g}, ψ={d:?}, λ={lambda:?}"
                        );
                    });
                }
            }
        }
        for (g, d, lambda) in [
            (5, vec![2, 2, 2], vec![(3, 3)]),
            (5, vec![0, 2, 4], vec![(1, 2), (3, 1), (4, 1)]),
            (6, vec![2, 2], vec![(3, 1), (5, 2)]),
        ] {
            assert_eq!(
                direct(g, &d, &lambda),
                legacy(g, &d, &lambda),
                "g={g}, ψ={d:?}, λ={lambda:?}"
            );
        }
    }

    #[test]
    fn duplicate_lambda_factors_and_zero_powers() {
        assert_eq!(
            direct(3, &[2, 2], &[(1, 1), (1, 3), (0, 7), (3, 0)]),
            legacy(3, &[2, 2], &[(1, 4)])
        );
        assert_eq!(direct(0, &[0, 0, 0, 1], &[(0, 2)]), Q::ONE);
        assert_eq!(direct(2, &[3], &[(1, 2)]), Q::ZERO);
        assert_eq!(direct(3, &[1], &[(3, 1), (3, 1)]), Q::ZERO);
    }

    #[test]
    fn parallel_matches_serial_descendants() {
        for (g, d, lambda) in [
            (1, vec![0, 1, 2], vec![]),
            (3, vec![0, 0, 2], vec![(1, 7)]),
            (4, vec![0, 1, 2, 3], vec![(3, 1), (4, 1)]),
            (5, vec![2, 2, 2], vec![(3, 3)]),
            (6, vec![2, 2], vec![(3, 1), (5, 2)]),
            (6, vec![2, 2, 2, 2], vec![(3, 1), (4, 2)]),
            (6, vec![5, 6, 7], vec![]),
            (7, vec![], vec![(1, 18)]),
        ] {
            let psi: Vec<_> = d
                .iter()
                .enumerate()
                .map(|(i, &p)| (i as u32 + 1, p))
                .collect();
            let serial = integrate_with_threads(g, d.len() as u32, &psi, &lambda, 1);
            for threads in [2, 4] {
                let psi = psi.clone();
                let lambda = lambda.clone();
                let expected = serial.clone();
                std::thread::spawn(move || {
                    for _ in 0..2 {
                        assert_eq!(
                            integrate_with_threads(g, psi.len() as u32, &psi, &lambda, threads),
                            expected,
                            "g={g}, ψ={psi:?}, λ={lambda:?}, threads={threads}"
                        );
                    }
                })
                .join()
                .unwrap();
            }
        }
    }

    #[test]
    fn mixed_and_pure_batches_match_serial_with_fresh_workers() {
        for g in [3, 5, 7] {
            for (kappa, lambda, ch) in [
                (vec![(1, 3)], vec![(1, 3 * g - 8)], vec![(3, 1)]),
                (vec![(1, 1), (2, 2)], vec![(1, 3 * g - 7)], vec![]),
                (vec![(1, 2)], vec![], vec![(1, 3 * g - 4)]),
            ] {
                let expected = integrate_mixed_with_threads(g, 1, &kappa, &[], &lambda, &ch, 1);
                for threads in [2, 4] {
                    std::thread::scope(|s| {
                        s.spawn(|| {
                            let value = integrate_mixed_with_threads(
                                g,
                                1,
                                &kappa,
                                &[],
                                &lambda,
                                &ch,
                                threads,
                            );
                            assert_eq!(
                                value, expected,
                                "g={g}, κ={kappa:?}, λ={lambda:?}, ch={ch:?}"
                            );
                        })
                        .join()
                        .unwrap()
                    });
                }
            }
        }
        std::thread::spawn(|| {
            let roots = vec![
                (vec![2; 12], Q::ONE),
                (vec![3; 6], Q::from(-2)),
                (vec![2, 2, 11], Q::from(3)),
            ];
            let mut serial = Descendents::default();
            let expected = roots.iter().fold(Q::ZERO, |sum, (d, c)| {
                sum + c * serial.eval(5, d.clone(), &[])
            });
            // Use an empty serial table so the parallel DAG must actually be constructed.
            let result = parallel::integrate(
                5,
                &roots,
                &[(vec![], Q::ONE)],
                &mut Descendents::default(),
                4,
                None,
            );
            assert_eq!(result, expected);
        })
        .join()
        .unwrap();
    }

    #[test]
    fn parallel_reuses_results_across_queries() {
        std::thread::spawn(|| {
            for g in [8, 6, 9, 7, 8] {
                let parallel =
                    integrate_with_threads(g, 0, &[], &[(1, 3 * g - 3)], 2 + g as usize % 3);
                let serial = integrate_with_threads(g, 0, &[], &[(1, 3 * g - 3)], 1);
                assert_eq!(parallel, serial, "g={g}");
            }
        })
        .join()
        .unwrap();
    }

    #[test]
    #[ignore = "手动比较串行与不同线程数的后代求值"]
    fn benchmark_parallel_threads() {
        for threads in [1, 2, 4, 8, 16, 24, 32] {
            let (value, elapsed, profile) = std::thread::spawn(move || {
                let start = std::time::Instant::now();
                let value = integrate_with_threads(11, 0, &[], &[(1, 30)], threads);
                (value, start.elapsed(), parallel::profile())
            })
            .join()
            .unwrap();
            assert_eq!(
                value.to_string(),
                "434910980875003666624744/494737211863125"
            );
            eprintln!(
                "threads={threads}: {elapsed:?}; planning={:?}, evaluation={:?}, new_states={}, terms={}",
                profile.0, profile.1, profile.2, profile.3
            );
        }
    }

    /// 每次在新线程计时：两种算法均从空的线程内递归表开始，不计线程启动。
    /// cargo test -p sc-faber --release benchmark_hodge_descendants -- --ignored --nocapture
    #[test]
    #[ignore = "手动比较 Hodge 递归性能"]
    fn benchmark_hodge_descendants() {
        use std::time::Instant;
        for (g, d, lambda) in [
            (5, vec![2, 2, 2], vec![(1, 9)]),
            (6, vec![2, 2, 2, 2], vec![(1, 11)]),
            (5, vec![2, 2, 2], vec![(3, 3)]),
            (6, vec![2, 2], vec![(3, 1), (5, 2)]),
            (7, vec![], vec![(1, 18)]),
            (8, vec![], vec![(1, 21)]),
            (6, vec![3, 3, 3, 3], vec![(6, 1), (1, 1)]),
            (8, vec![5, 5, 5], vec![(8, 1), (1, 1)]),
            (8, vec![5, 5, 6], vec![(8, 1)]),
        ] {
            let mut timings = [Vec::new(), Vec::new()];
            for _ in 0..3 {
                let mut values = Vec::new();
                for (algorithm, times) in timings.iter_mut().enumerate() {
                    let d = d.clone();
                    let lambda = lambda.clone();
                    let (value, elapsed) = std::thread::spawn(move || {
                        let start = Instant::now();
                        let value = if algorithm == 0 {
                            legacy(g, &d, &lambda)
                        } else {
                            direct(g, &d, &lambda)
                        };
                        (value, start.elapsed())
                    })
                    .join()
                    .unwrap();
                    times.push(elapsed);
                    values.push(value);
                }
                assert_eq!(values[0], values[1]);
            }
            timings.iter_mut().for_each(|times| times.sort_unstable());
            eprintln!(
                "g={g} ψ={d:?} λ={lambda:?}: old={:?}, direct={:?}, speedup={:.2}x",
                timings[0][1],
                timings[1][1],
                timings[0][1].as_secs_f64() / timings[1][1].as_secs_f64()
            );
        }
    }

    /// 同一线程顺序计算 ∫_{M_g} λ₁^{3g−3}，复用前面亏格积累的递归结果。
    /// cargo test -p sc-faber --release benchmark_lambda1_series -- --ignored --nocapture
    #[test]
    #[ignore = "手动计时 g=2..14 的连续 λ₁ 积分"]
    fn benchmark_lambda1_series() {
        use std::time::{Duration, Instant};

        let mut total = Duration::ZERO;
        for g in 2..=14 {
            let start = Instant::now();
            let value = direct(g, &[], &[(1, 3 * g - 3)]);
            let elapsed = start.elapsed();
            total += elapsed;
            eprintln!("g={g}: {elapsed:?}, cumulative={total:?}, integral={value}");
        }
        eprintln!("g=2..14 first sequence: {total:?}");

        let start = Instant::now();
        for g in 2..=14 {
            std::hint::black_box(direct(g, &[], &[(1, 3 * g - 3)]));
        }
        eprintln!("g=2..14 repeated sequence: {:?}", start.elapsed());
    }

    /// 分别以 SC_FABER_THREADS=1 和 24 启动新进程，对比 g=15 的首次计算。
    #[test]
    #[ignore = "手动计时 g=15 的 λ₁ 顶次积分"]
    fn benchmark_lambda1_g15() {
        let start = std::time::Instant::now();
        let value = direct(15, &[], &[(1, 42)]);
        let elapsed = start.elapsed();
        eprintln!("g=15, elapsed={elapsed:?}, integral={value}");
    }

    /// cargo test -p sc-faber --release benchmark_lambda1_to_g20 -- --ignored --nocapture
    #[test]
    #[ignore = "手动计算到 g=20 的连续 λ₁ 顶次积分"]
    fn benchmark_lambda1_to_g20() {
        let mut total = std::time::Duration::ZERO;
        for g in 2..=20 {
            eprintln!("starting g={g}");
            let start = std::time::Instant::now();
            let value = direct(g, &[], &[(1, 3 * g - 3)]);
            let elapsed = start.elapsed();
            total += elapsed;
            eprintln!("g={g}: {elapsed:?}, cumulative={total:?}, integral={value}");
            let (planning, evaluation, states, terms) = parallel::profile();
            eprintln!(
                "planning={planning:?}, evaluation={evaluation:?}, new_states={states}, terms={terms}"
            );
        }
        eprintln!("g=2..20 first sequence: {total:?}");
    }
}
