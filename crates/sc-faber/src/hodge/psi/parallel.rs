//! 按 (ch 总次数, 2g−2+n) 分层的后代递归。
//!
//! 只展开本次查询缺少的子问题，以有序表合并重复状态。同一层只读取低层的
//! 有理数结果；工作线程各自累加，整层完成后统一写回。跨查询只保留索引和
//! 结果，求和式随本次求值结束释放。串行递归仍由父模块保留。

use super::{ChMono, Descendents};
use sc_num::{One, Q, Zero};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock, mpsc};

type Index = BTreeMap<Vec<(u32, u32)>, BTreeMap<Vec<u32>, usize>>;
type SerialValues = BTreeMap<Vec<(u32, u32)>, BTreeMap<Vec<u32>, Q>>;

#[derive(Default)]
struct Engine {
    // 直接按两个数学次数分开索引，避免所有新状态争用同一把写锁。
    index: Vec<Vec<RwLock<Index>>>,
    values: Vec<Q>,
    #[cfg(test)]
    profile: (std::time::Duration, std::time::Duration, usize, usize),
}

thread_local! {
    static ENGINE: RefCell<Engine> = RefCell::new(Engine::default());
}

enum Operand {
    Constant(Q),
    Node(usize),
    Scaled(Q, usize),
}

impl Operand {
    fn scale(self, factor: Q) -> Self {
        match self {
            Self::Constant(value) => Self::Constant(factor * value),
            Self::Node(index) => Self::Scaled(factor, index),
            Self::Scaled(value, index) => Self::Scaled(factor * value, index),
        }
    }

    fn factor(self, coefficient: &mut Q) -> Option<usize> {
        match self {
            Self::Constant(value) => {
                *coefficient *= value;
                None
            }
            Self::Node(index) => Some(index),
            Self::Scaled(value, index) => {
                *coefficient *= value;
                Some(index)
            }
        }
    }
}

// GRR 的组合系数通常是小整数或半整数。直接存储它们，避免每条边携带
// 完整的大有理数；超出范围的系数仍保留任意精度，绝不截断。
enum Coefficient {
    Integer(i64),
    Half(i64),
    Rational(Box<Q>),
}

impl Coefficient {
    fn new(value: Q) -> Self {
        if let Ok(integer) = i64::try_from(&value) {
            return Self::Integer(integer);
        }
        if *value.denominator_ref() == 2u32 {
            if let Ok(numerator) = i64::try_from(&(&value * Q::from(2))) {
                return Self::Half(numerator);
            }
        }
        Self::Rational(Box::new(value))
    }

    fn rational(&self) -> Q {
        match self {
            Self::Integer(value) => Q::from(*value),
            Self::Half(value) => Q::from_signeds(*value, 2),
            Self::Rational(value) => value.as_ref().clone(),
        }
    }

    fn add(&mut self, other: &Self) {
        if let (Self::Integer(left), Self::Integer(right)) = (&mut *self, other) {
            if let Some(sum) = left.checked_add(*right) {
                *left = sum;
                return;
            }
        }
        *self = Self::new(self.rational() + other.rational());
    }

    fn multiply(&self, value: &Q) -> Q {
        match self {
            Self::Integer(1) => value.clone(),
            Self::Integer(-1) => -value,
            Self::Integer(factor) => Q::from(*factor) * value,
            Self::Half(factor) => Q::from_signeds(*factor, 2) * value,
            Self::Rational(factor) => factor.as_ref() * value,
        }
    }
}

struct Term {
    factor: Coefficient,
    left: usize,
    // usize::MAX 表示线性项，否则为两个子积分的乘积。
    right: usize,
}

// 常见系数是半整数，节点编号也远小于 2^32。半整数分子与两个 u32
// 共占 16 字节；超过任一范围的项放到 wide，保持任意精度和 usize 索引。
struct CompactTerm {
    numerator: i64,
    left: u32,
    right: u32,
}

const LINEAR: u32 = u32::MAX;
const WIDE: u32 = u32::MAX - 1;

#[derive(Default)]
struct TermBlock {
    compact: Vec<CompactTerm>,
    wide: Vec<Term>,
}

impl TermBlock {
    fn len(&self) -> usize {
        self.compact.len()
    }

    fn push(&mut self, term: Term) {
        let numerator = match &term.factor {
            Coefficient::Integer(value) => value.checked_mul(2),
            Coefficient::Half(value) => Some(*value),
            Coefficient::Rational(_) => None,
        };
        if let Some(numerator) = numerator.filter(|_| {
            term.left < WIDE as usize && (term.right == usize::MAX || term.right < WIDE as usize)
        }) {
            self.compact.push(CompactTerm {
                numerator,
                left: term.left as u32,
                right: if term.right == usize::MAX {
                    LINEAR
                } else {
                    term.right as u32
                },
            });
        } else {
            self.compact.push(CompactTerm {
                // Vec<Term> 长度受 isize::MAX / size_of::<Term>() 限制。
                numerator: self.wide.len() as i64,
                left: 0,
                right: WIDE,
            });
            self.wide.push(term);
        }
    }

    fn shrink_to_fit(&mut self) {
        self.compact.shrink_to_fit();
        self.wide.shrink_to_fit();
    }

    fn evaluate(
        &self,
        range: std::ops::Range<usize>,
        constant: &Q,
        divisor: Option<&Q>,
        values: &[Q],
    ) -> Q {
        let mut sum = constant.clone();
        for term in &self.compact[range] {
            let value = if term.right == WIDE {
                let wide = &self.wide[term.numerator as usize];
                let mut value = wide.factor.multiply(&values[wide.left]);
                if wide.right != usize::MAX {
                    value *= &values[wide.right];
                }
                value
            } else {
                let left = &values[term.left as usize];
                let mut value = match term.numerator {
                    2 => left.clone(),
                    -2 => -left,
                    n if n % 2 == 0 => Q::from(n / 2) * left,
                    n => Q::from_signeds(n, 2) * left,
                };
                if term.right != LINEAR {
                    value *= &values[term.right as usize];
                }
                value
            };
            sum += value;
        }
        if let Some(divisor) = divisor {
            sum /= divisor;
        }
        sum
    }
}

#[derive(Default)]
struct Equation {
    constant: Q,
    terms: Vec<Term>,
    divisor: Option<Q>,
}

impl Equation {
    fn add(&mut self, mut factor: Q, left: Operand, right: Operand) {
        let left = left.factor(&mut factor);
        let right = right.factor(&mut factor);
        if factor == Q::ZERO {
            return;
        }
        match (left, right) {
            (None, None) => self.constant += factor,
            (Some(left), None) | (None, Some(left)) => self.terms.push(Term {
                factor: Coefficient::new(factor),
                left,
                right: usize::MAX,
            }),
            (Some(a), Some(b)) => self.terms.push(Term {
                factor: Coefficient::new(factor),
                left: a.min(b),
                right: a.max(b),
            }),
        }
    }

    fn linear(&mut self, factor: Q, value: Operand) {
        self.add(factor, value, Operand::Constant(Q::ONE));
    }

    fn collect(&mut self) {
        self.terms
            .sort_unstable_by_key(|term| (term.left, term.right));
        self.terms.dedup_by(|term, previous| {
            if previous.left == term.left && previous.right == term.right {
                previous.factor.add(&term.factor);
                true
            } else {
                false
            }
        });
        self.terms
            .retain(|term| !matches!(term.factor, Coefficient::Integer(0)));
    }

    #[cfg(test)]
    fn evaluate(&self, values: &[Q]) -> Q {
        evaluate_terms(&self.constant, &self.terms, self.divisor.as_ref(), values)
    }
}

#[cfg(test)]
fn evaluate_terms(constant: &Q, terms: &[Term], divisor: Option<&Q>, values: &[Q]) -> Q {
    let mut sum = constant.clone();
    for term in terms {
        let mut value = term.factor.multiply(&values[term.left]);
        if term.right != usize::MAX {
            value *= &values[term.right];
        }
        sum += value;
    }
    if let Some(divisor) = divisor {
        sum /= divisor;
    }
    sum
}

#[derive(Default)]
struct PackedEquation {
    constant: Q,
    divisor: Option<Q>,
    block: usize,
    terms: std::ops::Range<usize>,
}

#[derive(Default)]
struct Plan {
    equations: Vec<PackedEquation>,
    // 每个展开批次只保留一个连续块；避免为数百万个求和式单独分配和释放。
    terms: Vec<TermBlock>,
    layers: Vec<Vec<usize>>,
}

impl Plan {
    fn evaluate(&self, id: usize, values: &[Q]) -> Q {
        let equation = &self.equations[id];
        self.terms[equation.block].evaluate(
            equation.terms.clone(),
            &equation.constant,
            equation.divisor.as_ref(),
            values,
        )
    }
}

pub(super) fn integrate(
    g: u32,
    descendants: &[(Vec<u32>, Q)],
    poly: &[ChMono],
    tables: &mut Descendents,
    threads: usize,
) -> Q {
    ENGINE.with(|cell| {
        #[cfg(test)]
        let start = std::time::Instant::now();
        let mut engine = cell.borrow_mut();
        let base = engine.values.len();
        let weight = poly
            .iter()
            .map(|(ch, _)| ch.iter().map(|&(m, p)| m * p).sum::<u32>())
            .max()
            .unwrap_or(0) as usize;
        // W_ch + (2g−2+n) 沿每条递归边不增加。
        let bound =
            (weight + 2 * g as usize + descendants.iter().map(|(d, _)| d.len()).max().unwrap_or(0))
                .saturating_sub(2);
        if engine.index.len() <= weight {
            engine.index.resize_with(weight + 1, Vec::new);
        }
        for row in &mut engine.index[..=weight] {
            if row.len() <= bound {
                row.resize_with(bound + 1, || RwLock::new(Index::new()));
            }
        }
        let next = AtomicUsize::new(base);
        let (roots, mut plan) = {
            let mut planner = Planner::new(&engine, &tables.values, &next);
            let mut roots = Vec::new();
            for (d, factor) in descendants {
                for (ch, c) in poly {
                    roots.push((planner.node(g, d.clone(), ch), factor * c));
                }
            }
            let plan = build_plan(&mut planner, threads);
            (roots, plan)
        };
        #[cfg(test)]
        let planning = start.elapsed();
        engine.values.resize(base + plan.equations.len(), Q::ZERO);
        #[cfg(test)]
        let start = std::time::Instant::now();
        evaluate_layers(&plan, &mut engine.values, base, threads);
        #[cfg(test)]
        {
            engine.profile = (
                planning,
                start.elapsed(),
                plan.equations.len(),
                plan.terms.iter().map(TermBlock::len).sum(),
            );
        }
        // 连续块分摊给工作线程释放，大有理数系数也在这里析构。
        if !plan.terms.is_empty() {
            let chunk = plan.terms.len().div_ceil(threads.max(1));
            std::thread::scope(|scope| {
                for blocks in plan.terms.chunks_mut(chunk) {
                    scope.spawn(move || {
                        for block in blocks {
                            *block = TermBlock::default();
                        }
                    });
                }
            });
        }
        roots.into_iter().fold(Q::ZERO, |sum, (root, coefficient)| {
            sum + coefficient
                * match root {
                    Operand::Constant(value) => value,
                    Operand::Node(index) => engine.values[index].clone(),
                    Operand::Scaled(factor, index) => factor * &engine.values[index],
                }
        })
    })
}

#[cfg(test)]
pub(super) fn profile() -> (std::time::Duration, std::time::Duration, usize, usize) {
    ENGINE.with(|cell| cell.borrow().profile)
}

struct Planner<'a> {
    engine: &'a Engine,
    serial: &'a SerialValues,
    tables: Descendents,
    next: &'a AtomicUsize,
    discovered: Vec<State>,
}

struct State {
    g: u32,
    d: Vec<u32>,
    ch: Vec<(u32, u32)>,
    rank: (u32, u32),
    index: usize,
}

impl<'a> Planner<'a> {
    fn new(engine: &'a Engine, serial: &'a SerialValues, next: &'a AtomicUsize) -> Self {
        Self {
            engine,
            serial,
            tables: Descendents::default(),
            next,
            discovered: Vec::new(),
        }
    }

    fn node(&mut self, g: u32, mut d: Vec<u32>, ch: &[(u32, u32)]) -> Operand {
        let n = d.len() as i64;
        let weight: u32 = ch.iter().map(|&(m, p)| m * p).sum();
        let euler = 2 * i64::from(g) - 2 + n;
        if euler <= 0
            || i64::from(weight) + d.iter().map(|&v| i64::from(v)).sum::<i64>()
                != 3 * i64::from(g) - 3 + n
        {
            return Operand::Constant(Q::ZERO);
        }
        d.sort_unstable();
        if ch.is_empty() {
            if g == 0 || d.len() == 1 {
                return Operand::Constant(self.tables.eval(g, d, ch));
            }
        } else if g == 0 || ch.last().unwrap().0 >= 2 * g || weight > (3 * g - 3).max(1) {
            return Operand::Constant(Q::ZERO);
        } else if g == 1 {
            return Operand::Constant(self.tables.multinomial(&d) / Q::from(2));
        }
        // Dilaton 可在 string 之前消掉任意位置的 τ₁；把纯倍数关系直接吸收
        // 到边系数中，不为它们建立独立任务或反复经过同步层。
        let zeros = d.partition_point(|&x| x == 0);
        let ones = d[zeros..].partition_point(|&x| x == 1);
        let remove = ones.min(if g == 1 { d.len() - 1 } else { d.len() });
        if remove > 0 {
            let mut factor = Q::ONE;
            for n in d.len() - remove..d.len() {
                factor *= Q::from(2 * u64::from(g) - 2 + n as u64);
            }
            d.drain(zeros..zeros + remove);
            return self.node(g, d, ch).scale(factor);
        }
        // String 只有一种正次数时，所有求和项相同，也只保留一个依赖。
        if zeros > 0 && d.get(zeros).is_some_and(|value| d.last() == Some(value)) {
            let multiplicity = d.len() - zeros;
            d.remove(0);
            d[zeros - 1] -= 1;
            return self.node(g, d, ch).scale(Q::from(multiplicity));
        }
        let shard = &self.engine.index[weight as usize][euler as usize];
        let cached = shard
            .read()
            .unwrap()
            .get(ch)
            .and_then(|values| values.get(&d))
            .copied();
        if let Some(index) = cached {
            if index < self.engine.values.len() && self.engine.values[index] == Q::ZERO {
                return Operand::Constant(Q::ZERO);
            }
            return Operand::Node(index);
        }
        // 小查询的串行结果也可用于之后的大查询。
        if let Some(value) = self.serial.get(ch).and_then(|values| values.get(&d)) {
            return Operand::Constant(value.clone());
        }
        let mut registry = shard.write().unwrap();
        let entry = registry.entry(ch.to_vec()).or_default();
        if let Some(&index) = entry.get(&d) {
            return Operand::Node(index);
        }
        let index = self.next.fetch_add(1, Ordering::Relaxed);
        entry.insert(d.clone(), index);
        drop(registry);
        self.discovered.push(State {
            g,
            d,
            ch: ch.to_vec(),
            rank: (weight, euler as u32),
            index,
        });
        Operand::Node(index)
    }

    fn expand(&mut self, g: u32, d: &[u32], ch: &[(u32, u32)], terms: Vec<Term>) -> Equation {
        let mut equation = Equation {
            terms,
            ..Equation::default()
        };
        if d.first() == Some(&0) {
            let mut i = 1;
            while i < d.len() {
                let end = i + d[i..].partition_point(|&x| x == d[i]);
                if d[i] > 0 {
                    let mut next = d[1..].to_vec();
                    next[i - 1] -= 1;
                    equation.linear(Q::from(end - i), self.node(g, next, ch));
                }
                i = end;
            }
            return equation;
        }
        if d.first() == Some(&1) {
            let count =
                d.partition_point(|&x| x == 1)
                    .min(if g == 1 { d.len() - 1 } else { d.len() });
            let mut factor = Q::ONE;
            for n in d.len() - count..d.len() {
                factor *= Q::from(2 * u64::from(g) - 2 + n as u64);
            }
            equation.linear(factor, self.node(g, d[count..].to_vec(), ch));
            return equation;
        }
        if ch.is_empty() {
            let k = d[0];
            let rest = &d[1..];
            let mut groups = Vec::new();
            let mut i = 0;
            while i < rest.len() {
                let end = i + rest[i..].partition_point(|&x| x == rest[i]);
                let mut merged = rest.to_vec();
                merged[i] += k - 1;
                let mut factor = Q::from(end - i);
                for j in 1..=k {
                    factor *= Q::from(2 * rest[i] + 2 * j - 1);
                }
                equation.linear(factor, self.node(g, merged, &[]));
                groups.push((rest[i], (end - i) as u32));
                i = end;
            }
            let weights: Vec<Q> = (0..k - 1)
                .map(|a| self.tables.odd_factorial(a) * self.tables.odd_factorial(k - 2 - a))
                .collect();
            self.boundary(g, rest, &[], &groups, &weights, &mut equation);
            equation.divisor = Some(self.tables.odd_factorial(k));
            return equation;
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
        equation.linear(Q::ONE, self.node(g, added, &rest));
        let mut groups = Vec::new();
        let mut i = 0;
        while i < d.len() {
            let end = i + d[i..].partition_point(|&x| x == d[i]);
            let mut shifted = d.to_vec();
            shifted[i] += m;
            equation.linear(-Q::from(end - i), self.node(g, shifted, &rest));
            groups.push((d[i], (end - i) as u32));
            i = end;
        }
        let weights: Vec<Q> = (0..m)
            .map(|a| Q::from(if a % 2 == 0 { 1 } else { -1 }))
            .collect();
        self.boundary(g, d, &rest, &groups, &weights, &mut equation);
        equation
    }

    fn boundary(
        &mut self,
        g: u32,
        d: &[u32],
        ch: &[(u32, u32)],
        groups: &[(u32, u32)],
        weights: &[Q],
        equation: &mut Equation,
    ) {
        let m = weights.len() as u32;
        for a in 0..m.div_ceil(2) {
            let mut node = d.to_vec();
            node.extend([a, m - 1 - a]);
            let mut factor = weights[a as usize].clone();
            if 2 * a == m - 1 {
                factor /= Q::from(2);
            }
            equation.linear(factor, self.node(g - 1, node, ch));
        }
        Distribution {
            planner: self,
            equation,
            g,
            m,
            groups,
            ch,
            weights,
            left_d: Vec::new(),
            right_d: Vec::new(),
            left_ch: Vec::new(),
            right_ch: Vec::new(),
        }
        .distribute(0, 0, Q::ONE);
    }
}

#[derive(Clone)]
struct BuildJob {
    states: Arc<Vec<State>>,
    start: usize,
    end: usize,
}

struct Built {
    discovered: Vec<State>,
    equations: Vec<(usize, (u32, u32), PackedEquation)>,
    terms: TermBlock,
}

impl Built {
    fn store(
        self,
        plan: &mut Plan,
        layers: &mut BTreeMap<(u32, u32), Vec<usize>>,
        base: usize,
        count: usize,
    ) -> Vec<State> {
        plan.equations.resize_with(count, PackedEquation::default);
        let block = plan.terms.len();
        plan.terms.push(self.terms);
        for (id, rank, mut equation) in self.equations {
            equation.block = block;
            plan.equations[id - base] = equation;
            layers.entry(rank).or_default().push(id - base);
        }
        self.discovered
    }
}

impl Planner<'_> {
    fn build(&mut self, job: BuildJob) -> Built {
        let mut equations = Vec::with_capacity(job.end - job.start);
        let mut terms = TermBlock::default();
        let mut scratch = Vec::new();
        for state in &job.states[job.start..job.end] {
            let mut equation = self.expand(state.g, &state.d, &state.ch, scratch);
            equation.collect();
            let start = terms.len();
            for term in equation.terms.drain(..) {
                terms.push(term);
            }
            scratch = equation.terms;
            equations.push((
                state.index,
                state.rank,
                PackedEquation {
                    constant: equation.constant,
                    divisor: equation.divisor,
                    block: 0,
                    terms: start..terms.len(),
                },
            ));
        }
        terms.shrink_to_fit();
        Built {
            discovered: std::mem::take(&mut self.discovered),
            equations,
            terms,
        }
    }
}

fn build_plan(planner: &mut Planner<'_>, threads: usize) -> Plan {
    let mut pending = std::mem::take(&mut planner.discovered);
    if pending.is_empty() {
        return Plan::default();
    }
    let base = planner.engine.values.len();
    let mut plan = Plan::default();
    let mut layers: BTreeMap<(u32, u32), Vec<usize>> = BTreeMap::new();
    // 缓存命中后的少量补算不启动线程池；依赖前沿足够宽时再并行展开。
    while !pending.is_empty() && pending.len() < 128 {
        let states = Arc::new(std::mem::take(&mut pending));
        let built = planner.build(BuildJob {
            end: states.len(),
            states,
            start: 0,
        });
        pending = built.store(
            &mut plan,
            &mut layers,
            base,
            planner.next.load(Ordering::Relaxed) - base,
        );
    }
    if pending.is_empty() {
        plan.layers = layers.into_values().collect();
        return plan;
    }
    let (jobs_tx, jobs_rx) = mpsc::channel::<BuildJob>();
    let jobs_rx = Mutex::new(jobs_rx);
    let (reply_tx, reply_rx) = mpsc::channel();
    let outcome = std::thread::scope(|scope| {
        for _ in 1..threads {
            let mut worker = Planner::new(planner.engine, planner.serial, planner.next);
            let reply_tx = reply_tx.clone();
            let jobs_rx = &jobs_rx;
            scope.spawn(move || {
                loop {
                    let Ok(job) = jobs_rx.lock().unwrap().recv() else {
                        break;
                    };
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        worker.build(job)
                    }));
                    if reply_tx.send(result).is_err() {
                        break;
                    }
                }
            });
        }
        drop(reply_tx);
        let result = (|| {
            while !pending.is_empty() {
                let states = Arc::new(std::mem::take(&mut pending));
                let chunk = if states.len() < 128 {
                    states.len()
                } else {
                    states.len().div_ceil(threads * 4).max(32)
                };
                let mut jobs = 0;
                for start in (chunk..states.len()).step_by(chunk) {
                    jobs_tx
                        .send(BuildJob {
                            states: Arc::clone(&states),
                            start,
                            end: (start + chunk).min(states.len()),
                        })
                        .unwrap();
                    jobs += 1;
                }
                let built = planner.build(BuildJob {
                    states: Arc::clone(&states),
                    start: 0,
                    end: chunk.min(states.len()),
                });
                pending.extend(built.store(
                    &mut plan,
                    &mut layers,
                    base,
                    planner.next.load(Ordering::Relaxed) - base,
                ));
                for _ in 0..jobs {
                    let built = reply_rx.recv().unwrap()?;
                    pending.extend(built.store(
                        &mut plan,
                        &mut layers,
                        base,
                        planner.next.load(Ordering::Relaxed) - base,
                    ));
                }
            }
            Ok(())
        })();
        drop(jobs_tx);
        result
    });
    if let Err(payload) = outcome {
        std::panic::resume_unwind(payload);
    }
    plan.layers = layers.into_values().collect();
    plan
}

struct Distribution<'a, 'b> {
    planner: &'a mut Planner<'b>,
    equation: &'a mut Equation,
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

impl Distribution<'_, '_> {
    fn distribute(&mut self, index: usize, defect: u32, coefficient: Q) {
        if index == self.groups.len() + self.ch.len() {
            self.emit(defect, coefficient);
            return;
        }
        let is_psi = index < self.groups.len();
        let (degree, count) = if is_psi {
            self.groups[index]
        } else {
            self.ch[index - self.groups.len()]
        };
        let left_len = self.left_d.len();
        let right_len = self.right_d.len();
        let degree = if is_psi { degree - 1 } else { degree };
        let bound = 3 * (self.g / 2) - 2;
        let max_taken = count.min((bound - defect) / degree);
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
            let weight = &coefficient * self.planner.tables.binomial(count, taken);
            self.distribute(index + 1, defect + degree * taken, weight);
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
    }

    fn emit(&mut self, defect: u32, coefficient: Q) {
        let min_h = ((defect + 4) / 3).max(1);
        let max_h = ((defect + self.m + 1) / 3).min(self.g / 2);
        let left_weight: u32 = self.left_ch.iter().map(|&(m, p)| m * p).sum();
        let right_weight: u32 = self.right_ch.iter().map(|&(m, p)| m * p).sum();
        for h in min_h..=max_h {
            if left_weight > (3 * h - 3).max(1)
                || right_weight > (3 * (self.g - h) - 3).max(1)
                || self.left_ch.last().is_some_and(|&(m, _)| m >= 2 * h)
                || self
                    .right_ch
                    .last()
                    .is_some_and(|&(m, _)| m >= 2 * (self.g - h))
            {
                continue;
            }
            let a = 3 * h - 2 - defect;
            let b = self.m - 1 - a;
            let mut factor = &coefficient * &self.weights[a as usize];
            if 2 * h == self.g {
                let order = self
                    .left_d
                    .cmp(&self.right_d)
                    .then(self.left_ch.cmp(&self.right_ch))
                    .then(a.cmp(&b));
                if order.is_gt() {
                    continue;
                }
                if order.is_eq() {
                    factor /= Q::from(2);
                }
            }
            let mut left = self.left_d.clone();
            left.push(a);
            let left = self.planner.node(h, left, &self.left_ch);
            if matches!(&left, Operand::Constant(value) if *value == Q::ZERO) {
                continue;
            }
            let mut right = self.right_d.clone();
            right.push(b);
            let right = self.planner.node(self.g - h, right, &self.right_ch);
            self.equation.add(factor, left, right);
        }
    }
}

#[derive(Clone, Copy)]
struct Job {
    layer: usize,
    start: usize,
    end: usize,
}

fn evaluate_layers(plan: &Plan, values: &mut Vec<Q>, base: usize, threads: usize) {
    if plan.equations.is_empty() {
        return;
    }
    let layers = &plan.layers;
    let width = layers.iter().map(Vec::len).max().unwrap_or(0);
    let threads = if width < 128 {
        1
    } else {
        threads.min(width.div_ceil(32)).max(1)
    };
    if threads == 1 {
        for layer in layers {
            for &id in layer {
                values[base + id] = plan.evaluate(id, values);
            }
        }
        return;
    }
    let stored = RwLock::new(values);
    let (jobs_tx, jobs_rx) = mpsc::channel::<Job>();
    let jobs_rx = Mutex::new(jobs_rx);
    let (reply_tx, reply_rx) = mpsc::channel();
    let calculate = |job: Job| {
        let values = stored.read().unwrap();
        layers[job.layer][job.start..job.end]
            .iter()
            .map(|&id| (id, plan.evaluate(id, &values)))
            .collect::<Vec<_>>()
    };
    let outcome = std::thread::scope(|scope| {
        for _ in 1..threads {
            let reply_tx = reply_tx.clone();
            let jobs_rx = &jobs_rx;
            let calculate = &calculate;
            scope.spawn(move || {
                loop {
                    let Ok(job) = jobs_rx.lock().unwrap().recv() else {
                        break;
                    };
                    let result =
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| calculate(job)));
                    if reply_tx.send(result).is_err() {
                        break;
                    }
                }
            });
        }
        drop(reply_tx);
        let result = (|| {
            for (layer_index, layer) in layers.iter().enumerate() {
                if layer.len() < 128 {
                    let mut values = stored.write().unwrap();
                    for &id in layer {
                        let result = plan.evaluate(id, &values);
                        values[base + id] = result;
                    }
                    continue;
                }
                let chunk = layer.len().div_ceil(threads * 4).max(32);
                let mut jobs = 0;
                for start in (chunk..layer.len()).step_by(chunk) {
                    jobs_tx
                        .send(Job {
                            layer: layer_index,
                            start,
                            end: (start + chunk).min(layer.len()),
                        })
                        .unwrap();
                    jobs += 1;
                }
                let mut completed = calculate(Job {
                    layer: layer_index,
                    start: 0,
                    end: chunk.min(layer.len()),
                });
                for _ in 0..jobs {
                    completed.extend(reply_rx.recv().unwrap()?);
                }
                let mut values = stored.write().unwrap();
                for (id, value) in completed {
                    values[base + id] = value;
                }
            }
            Ok(())
        })();
        drop(jobs_tx);
        result
    });
    if let Err(payload) = outcome {
        std::panic::resume_unwind(payload);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_coefficients_preserve_large_and_fractional_values() {
        let value = Q::from_signeds(7, 11);
        for factor in [
            Q::from(i64::MIN),
            Q::from(i64::MAX),
            Q::from(i64::MAX) * Q::from(2),
            Q::from_signeds(3, 2),
            Q::from_signeds(-3, 2),
            Q::from_signeds(5, 7),
            Q::from(i64::MAX) * Q::from_signeds(i64::MAX, 3),
        ] {
            let coefficient = Coefficient::new(factor.clone());
            assert_eq!(coefficient.rational(), factor);
            assert_eq!(coefficient.multiply(&value), &factor * &value);
        }
    }

    #[test]
    fn packed_terms_keep_exact_coefficients_and_wide_indices() {
        assert_eq!(std::mem::size_of::<CompactTerm>(), 16);
        let values = [Q::from_signeds(7, 11), Q::from_signeds(-13, 17)];
        let mut block = TermBlock::default();
        let mut expected = Q::from(3);
        for factor in [
            Q::ONE,
            -Q::ONE,
            Q::from_signeds(-3, 2),
            Q::from_signeds(5, 7),
            Q::from(i64::MIN),
            Q::from(i64::MAX) * Q::from(2),
        ] {
            for right in [usize::MAX, 1] {
                let mut value = &factor * &values[0];
                if right != usize::MAX {
                    value *= &values[right];
                }
                expected += value;
                block.push(Term {
                    factor: Coefficient::new(factor.clone()),
                    left: 0,
                    right,
                });
            }
        }
        assert_eq!(
            block.evaluate(0..block.len(), &Q::from(3), Some(&Q::from(23)), &values),
            expected / Q::from(23)
        );
        for (left, right) in [
            (WIDE as usize, usize::MAX),
            (0, WIDE as usize),
            (usize::MAX, usize::MAX),
        ] {
            block.push(Term {
                factor: Coefficient::Integer(1),
                left,
                right,
            });
            let packed = block.compact.last().unwrap();
            assert_eq!(packed.right, WIDE);
            let wide = &block.wide[packed.numerator as usize];
            assert_eq!((wide.left, wide.right), (left, right));
        }
    }

    #[test]
    fn collecting_terms_handles_overflow_and_cancellation() {
        let mut equation = Equation::default();
        equation.linear(Q::from(i64::MAX), Operand::Node(0));
        equation.linear(Q::ONE, Operand::Node(0));
        equation.linear(Q::from_signeds(3, 2), Operand::Node(1));
        equation.linear(Q::from_signeds(-3, 2), Operand::Node(1));
        equation.collect();
        assert_eq!(equation.terms.len(), 1);
        let value = Q::from_signeds(7, 11);
        assert_eq!(
            equation.evaluate(&[value.clone(), Q::from(5)]),
            (Q::from(i64::MAX) + Q::ONE) * value
        );
        equation.linear(Q::from(i64::MIN), Operand::Node(0));
        equation.collect();
        assert!(equation.terms.is_empty());
    }
}
