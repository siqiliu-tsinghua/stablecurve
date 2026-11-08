//! 按 (ch 总次数, 2g−2+n) 分层的后代递归。
//!
//! 只展开本次查询缺少的子问题，以有序表合并重复状态。同一层只读取低层的
//! 有理数结果；工作线程各自累加，整层完成后统一写回。跨查询只保留索引和
//! 结果，求和式随本次求值结束释放。串行递归仍由父模块保留。

use super::{ChMono, Descendents, SerialSource};
use crate::cache_source::Registry;
use sc_num::{One, Q, Zero};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock, mpsc};

mod pressure;
mod recovery;
mod source;
mod workspace;
use source::{Batch, ParallelSource};

type Index = BTreeMap<Vec<(u32, u32)>, BTreeMap<Vec<u32>, usize>>;
type SerialValues = BTreeMap<Vec<(u32, u32)>, BTreeMap<Vec<u32>, Q>>;
static EMPTY_SERIAL: SerialValues = BTreeMap::new();
static REPORT_STAGES: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

#[cfg(test)]
mod measure;
#[cfg(test)]
mod measure_identity;
#[cfg(test)]
mod measure_memory;

#[derive(Default)]
struct Engine {
    // 直接按两个数学次数分开索引，避免所有新状态争用同一把写锁。
    index: Vec<Vec<Arc<RwLock<Index>>>>,
    values: Vec<Q>,
    source: Option<Arc<ParallelSource>>,
    #[cfg(test)]
    profile: (std::time::Duration, std::time::Duration, usize, usize),
}

thread_local! {
    static ENGINE: RefCell<Engine> = RefCell::new(Engine::default());
    static SHARED_ENGINES: RefCell<BTreeMap<u64, (std::sync::Weak<Registry>, Engine)>> = const { RefCell::new(BTreeMap::new()) };
}

pub(super) fn clear_cache_context(id: u64) {
    SHARED_ENGINES.with(|cell| {
        cell.borrow_mut().remove(&id);
    });
}

fn with_engine<R>(registry: Option<&Arc<Registry>>, work: impl FnOnce(&mut Engine) -> R) -> R {
    let Some(registry) = registry else {
        return ENGINE.with(|cell| work(&mut cell.borrow_mut()));
    };
    SHARED_ENGINES.with(|cell| {
        let mut engines = cell.borrow_mut();
        if let Some((_, engine)) = engines.get_mut(&registry.id()) {
            let source = Arc::clone(engine.source.as_ref().unwrap());
            let _use = source.begin_use();
            return work(engine);
        }
        engines.retain(|_, (owner, _)| owner.strong_count() > 0);
        let (_, engine) = engines.entry(registry.id()).or_insert_with(|| {
            let source = Arc::new(ParallelSource::default());
            registry.register(source.clone());
            (
                Arc::downgrade(registry),
                Engine {
                    source: Some(source),
                    ..Engine::default()
                },
            )
        });
        let source = Arc::clone(engine.source.as_ref().unwrap());
        let _use = source.begin_use();
        work(engine)
    })
}

pub(super) fn cached(g: u32, d: &[u32], ch: &[(u32, u32)], registry: &Arc<Registry>) -> Option<Q> {
    let weight = ch
        .iter()
        .map(|&(m, p)| u64::from(m) * u64::from(p))
        .sum::<u64>();
    let euler: usize = (2 * i64::from(g) - 2 + d.len() as i64).try_into().ok()?;
    SHARED_ENGINES.with(|cell| {
        let engines = cell.borrow();
        let (_, engine) = engines.get(&registry.id())?;
        let source = engine.source.as_ref()?;
        // Match the reclaimer's order: use permit before any shard/value lock.
        // This also makes its history publication visible before a new miss.
        let _use = source.begin_use();
        let shard = engine
            .index
            .get(usize::try_from(weight).ok()?)?
            .get(euler)?;
        // Keep the original key-to-ID access until the value has been copied.
        // Pressure reclamation clears these same shard locks and the vector.
        let index = shard.read().unwrap();
        let id = *index.get(ch)?.get(d)?;
        let values = source.read();
        let value = values.get(id).cloned();
        drop(values);
        drop(index);
        value
    })
}

pub(super) fn has_cached_engine(registry: &Registry) -> bool {
    SHARED_ENGINES.with(|cell| cell.borrow().contains_key(&registry.id()))
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
        if self.compact.len() == self.compact.capacity() {
            let needed = self.compact.len() + 1;
            pressure::capacity(&mut self.compact, needed);
        }
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
            if self.wide.len() == self.wide.capacity() {
                let needed = self.wide.len() + 1;
                pressure::capacity(&mut self.wide, needed);
            }
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
        if self.compact.len() < self.compact.capacity() {
            pressure::reallocation(
                self.compact
                    .len()
                    .saturating_mul(std::mem::size_of::<CompactTerm>()),
            );
        }
        if self.wide.len() < self.wide.capacity() {
            pressure::reallocation(self.wide.len().saturating_mul(std::mem::size_of::<Term>()));
        }
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
        if (left.is_some() || right.is_some()) && self.terms.len() == self.terms.capacity() {
            let needed = self.terms.len() + 1;
            pressure::capacity(&mut self.terms, needed);
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
    // Only the historical-recovery path populates this list. These IDs are
    // already committed and must not become new dirty persistence records.
    restored: Vec<usize>,
    // Only pressure moves sealed blocks here. Ordinary evaluation retains its
    // original contiguous, directly indexed term vectors.
    workspace: Option<Arc<workspace::PlanWorkspace>>,
    // Updated only for failed waves when stage diagnostics are enabled.
    rollback_count: usize,
    rollback_states: usize,
    rollback_elapsed: std::time::Duration,
}

impl Drop for Plan {
    fn drop(&mut self) {
        if let Some(workspace) = &self.workspace {
            workspace.retire();
        }
    }
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
    registry: Option<&Arc<Registry>>,
) -> Q {
    with_engine(registry, |engine| {
        let base = engine
            .source
            .as_ref()
            .map_or_else(|| engine.values.len(), |source| source.read().len());
        let control = registry.map(|_| {
            let weight = poly
                .iter()
                .map(|(ch, _)| {
                    ch.iter()
                        .map(|&(degree, power)| degree as usize * power as usize)
                        .sum::<usize>()
                })
                .max()
                .unwrap_or(0);
            let marks = descendants.iter().map(|(d, _)| d.len()).max().unwrap_or(0);
            let max_marks = weight
                .saturating_add(2usize.saturating_mul(g as usize))
                .saturating_add(marks);
            let key_bytes = max_marks
                .saturating_mul(std::mem::size_of::<u32>())
                .saturating_add(weight.saturating_mul(std::mem::size_of::<(u32, u32)>()));
            let capacity = engine
                .source
                .as_ref()
                .map_or(0, |source| source.read().capacity());
            pressure::Control::new(base, capacity, key_bytes)
        });
        let next = AtomicUsize::new(base);
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            integrate_engine(
                engine,
                g,
                descendants,
                poly,
                tables,
                threads,
                registry,
                &next,
                base,
                control.as_ref(),
            )
        }));
        match outcome {
            Ok(value) => {
                if let Some(control) = &control {
                    control.finish();
                }
                value
            }
            Err(payload) => {
                let assigned = next.load(Ordering::Relaxed);
                if let Some(control) = &control {
                    control.start_cleanup();
                }
                let completed = engine
                    .source
                    .as_ref()
                    .map(|source| source.completed_since(base, assigned));
                for row in &engine.index {
                    for shard in row {
                        shard
                            .write()
                            .unwrap_or_else(|error| error.into_inner())
                            .retain(|_, values| {
                                values.retain(|_, id| {
                                    *id < base
                                        || completed
                                            .as_ref()
                                            .is_some_and(|complete| complete[*id - base])
                                });
                                !values.is_empty()
                            });
                        shard.clear_poison();
                    }
                }
                if let Some(source) = &engine.source {
                    // Reserve every assigned ID, including unfinished holes;
                    // the backend may still be reading this generation's keys.
                    if let Some(control) = &control {
                        control.resize_values(source, assigned);
                    } else {
                        source.write().resize(assigned, Q::ZERO);
                    }
                } else {
                    engine.values.truncate(base);
                }
                pressure::propagate(payload)
            }
        }
    })
}

fn integrate_engine(
    engine: &mut Engine,
    g: u32,
    descendants: &[(Vec<u32>, Q)],
    poly: &[ChMono],
    tables: &mut Descendents,
    threads: usize,
    registry: Option<&Arc<Registry>>,
    next: &AtomicUsize,
    base: usize,
    control: Option<&pressure::Control>,
) -> Q {
    #[cfg(test)]
    let mut measurement = measure::Run::start();
    let stage_start = (registry.is_some()
        && *REPORT_STAGES
            .get_or_init(|| std::env::var("SC_FABER_CACHE_MEASURE_REPORT").as_deref() == Ok("1")))
    .then(std::time::Instant::now);
    #[cfg(test)]
    let start = std::time::Instant::now();
    let source = engine.source.clone();
    let recovering = registry.is_some_and(|registry| {
        registry
            .query_descendant_history()
            .unwrap_or_else(|error| crate::cache_source::resource_failure(error))
    });
    // The historical path waits for a frontier read. It must not retain a
    // result-table or serial-memo guard while waiting for the storage owner.
    let stored = source
        .as_ref()
        .filter(|_| !recovering)
        .map(|source| source.read());
    let values = stored.as_deref().unwrap_or(&engine.values);
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
    let mut new_shards = Vec::new();
    for row in &mut engine.index[..=weight] {
        if row.len() <= bound {
            row.resize_with(bound + 1, || {
                let shard = Arc::new(RwLock::new(Index::new()));
                if source.is_some() {
                    new_shards.push(Arc::clone(&shard));
                }
                shard
            });
        }
    }
    if let Some(source) = &source {
        source.observe_shards(new_shards);
    }
    let (roots, mut planned) = if recovering {
        let mut planner =
            Planner::recovering(engine, tables.source.as_deref(), base, next, registry);
        planner.control = control;
        let roots = build_roots(&mut planner, g, descendants, poly);
        let plan = build_plan(&mut planner, threads);
        (roots, plan)
    } else {
        tables.with_memo(|serial| {
            let mut planner = Planner::new(engine, serial, values, next, registry);
            planner.control = control;
            let roots = build_roots(&mut planner, g, descendants, poly);
            let plan = build_plan(&mut planner, threads);
            (roots, plan)
        })
    };
    drop(stored);
    let mut plan = loop {
        match planned {
            Planning::Complete(plan) => break plan,
            Planning::Paused(mut draft) => {
                let registry = registry.expect("only cache-backed planning can pause");
                // All borrowed result/serial guards and the old worker scope have
                // ended before this first workspace registration or native sample.
                if draft.plan.workspace.is_none() && !draft.refreshed && pressure::refresh_normal()
                {
                    // An expired observation alone does not require disk backing.
                    // Re-enter the original borrowed path after the safe refresh;
                    // this also keeps source-only diagnostic owners usable.
                    draft.refreshed = true;
                    planned = if recovering {
                        let mut planner = Planner::recovering(
                            engine,
                            tables.source.as_deref(),
                            base,
                            next,
                            Some(registry),
                        );
                        planner.control = control;
                        resume_plan(&mut planner, threads, draft)
                    } else {
                        let stored = source.as_ref().map(|source| source.read());
                        let values = stored.as_deref().unwrap_or(&engine.values);
                        tables.with_memo(|serial| {
                            let mut planner =
                                Planner::new(engine, serial, values, next, Some(registry));
                            planner.control = control;
                            resume_plan(&mut planner, threads, draft)
                        })
                    };
                    continue;
                }
                pressure::back_plan(&mut draft.plan, registry);
                let mut planner = Planner::recovering(
                    engine,
                    tables.source.as_deref(),
                    base,
                    next,
                    Some(registry),
                );
                planner.lookup_history = recovering;
                planner.control = control;
                planned = resume_plan(&mut planner, threads, draft);
            }
        }
    };
    if let Some(source) = &source {
        let reused = roots
            .iter()
            .filter(|(root, _)| match root {
                Operand::Node(id) | Operand::Scaled(_, id) => *id < base,
                Operand::Constant(_) => false,
            })
            .count();
        source.note_query(reused, roots.len());
    }
    #[cfg(test)]
    let planning = start.elapsed();
    #[cfg(test)]
    measurement.planned(engine, &plan);
    let planning_elapsed = stage_start.map(|start| start.elapsed());
    if let Some(planning) = planning_elapsed {
        eprintln!(
            "sc-faber: planning finished genus={g} planning_s={:.6} equations={} blocks={} backed={} failed_waves={} rollback_states={} rollback_s={:.6}",
            planning.as_secs_f64(),
            plan.equations.len(),
            plan.terms.len(),
            plan.workspace.is_some(),
            plan.rollback_count,
            plan.rollback_states,
            plan.rollback_elapsed.as_secs_f64(),
        );
    }
    if let Some(source) = &source {
        if let Some(control) = control {
            control.resize_values(source, base + plan.equations.len());
        } else {
            source.write().resize(base + plan.equations.len(), Q::ZERO);
        }
    } else {
        engine.values.resize(base + plan.equations.len(), Q::ZERO);
    }
    if control.is_some() && plan.workspace.is_none() && !pressure::normal_at_boundary() {
        pressure::back_plan(&mut plan, registry.unwrap());
    }
    #[cfg(test)]
    measurement.resized();
    #[cfg(test)]
    let start = std::time::Instant::now();
    if let Some(source) = &source {
        if !plan.equations.is_empty() {
            let batch = if plan.restored.is_empty() {
                source.begin(base, std::mem::take(&mut plan.layers))
            } else {
                source.begin_restored(
                    base,
                    std::mem::take(&mut plan.layers),
                    std::mem::take(&mut plan.restored),
                )
            };
            if plan.workspace.is_some() {
                pressure::evaluate(&plan, source, &batch, base, threads);
            } else {
                match evaluate_shared(&plan, source, &batch, base, threads) {
                    Evaluation::Complete => {}
                    Evaluation::PausedEvaluation {
                        next_layer,
                        outputs,
                    } => {
                        pressure::back_plan(&mut plan, registry.unwrap());
                        pressure::evaluate_from(
                            &plan, source, &batch, base, threads, next_layer, outputs,
                        );
                    }
                }
            }
        }
    } else {
        #[cfg(not(test))]
        evaluate_layers(&plan, &mut engine.values, base, threads);
        #[cfg(test)]
        measure::evaluate(&mut plan, &mut engine.values, base, threads);
    }
    if let (Some(start), Some(planning)) = (stage_start, planning_elapsed) {
        eprintln!(
            "sc-faber: stages genus={g} planning_s={:.6} evaluation_s={:.6} equations={} blocks={} backed={}",
            planning.as_secs_f64(),
            start.elapsed().saturating_sub(planning).as_secs_f64(),
            plan.equations.len(),
            plan.terms.len(),
            plan.workspace.is_some(),
        );
    }
    #[cfg(test)]
    {
        engine.profile = (
            planning,
            start.elapsed(),
            plan.equations.len(),
            plan.terms.iter().map(TermBlock::len).sum(),
        );
    }
    #[cfg(test)]
    if source.is_none() {
        measure::reevaluate_same_plan(&mut plan, &mut engine.values, base, threads);
    }
    #[cfg(test)]
    measurement.evaluated(engine);
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
    #[cfg(test)]
    measurement.terms_released();
    let stored = source.as_ref().map(|source| source.read());
    let values = stored.as_deref().unwrap_or(&engine.values);
    let result = roots.into_iter().fold(Q::ZERO, |sum, (root, coefficient)| {
        sum + coefficient
            * match root {
                Operand::Constant(value) => value,
                Operand::Node(index) => values[index].clone(),
                Operand::Scaled(factor, index) => factor * &values[index],
            }
    });
    #[cfg(test)]
    {
        drop(plan);
        measurement.finish(engine);
    }
    result
}

#[cfg(test)]
pub(super) fn profile() -> (std::time::Duration, std::time::Duration, usize, usize) {
    ENGINE.with(|cell| cell.borrow().profile)
}

fn build_roots(
    planner: &mut Planner<'_>,
    g: u32,
    descendants: &[(Vec<u32>, Q)],
    poly: &[ChMono],
) -> Vec<(Operand, Q)> {
    let wave = planner.control.map(|_| {
        pressure::Wave::begin_root(descendants.len().saturating_mul(poly.len()), 1)
            .unwrap_or_else(|error| crate::cache_source::resource_failure(error))
    });
    let roots = pressure::Wave::with(wave.as_ref(), || {
        let mut roots = Vec::new();
        for (d, factor) in descendants {
            for (ch, c) in poly {
                if roots.len() == roots.capacity() {
                    let needed = roots.len() + 1;
                    pressure::capacity(&mut roots, needed);
                }
                roots.push((planner.node(g, d.clone(), ch), factor * c));
            }
        }
        roots
    });
    if let Some(wave) = wave {
        wave.finish();
    }
    roots
}

struct Planner<'a> {
    engine: &'a Engine,
    serial: &'a SerialValues,
    values: &'a [Q],
    shared_serial: Option<&'a SerialSource>,
    base: usize,
    pressure_view: bool,
    lookup_history: bool,
    control: Option<&'a pressure::Control>,
    registry: Option<&'a Arc<Registry>>,
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
    fn new(
        engine: &'a Engine,
        serial: &'a SerialValues,
        values: &'a [Q],
        next: &'a AtomicUsize,
        registry: Option<&'a Arc<Registry>>,
    ) -> Self {
        Self {
            engine,
            serial,
            values,
            shared_serial: None,
            base: values.len(),
            pressure_view: false,
            lookup_history: false,
            control: None,
            registry,
            tables: Descendents::with_registry(registry),
            next,
            discovered: Vec::new(),
        }
    }

    fn recovering(
        engine: &'a Engine,
        serial: Option<&'a SerialSource>,
        base: usize,
        next: &'a AtomicUsize,
        registry: Option<&'a Arc<Registry>>,
    ) -> Self {
        let mut planner = Self::new(engine, &EMPTY_SERIAL, &[], next, registry);
        planner.shared_serial = serial;
        planner.base = base;
        planner.pressure_view = true;
        planner.lookup_history = true;
        planner
    }

    fn worker(&self) -> Self {
        let mut worker = if self.pressure_view {
            Self::recovering(
                self.engine,
                self.shared_serial,
                self.base,
                self.next,
                self.registry,
            )
        } else {
            Self::new(
                self.engine,
                self.serial,
                self.values,
                self.next,
                self.registry,
            )
        };
        worker.lookup_history = self.lookup_history;
        worker.control = self.control;
        worker
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
            return Operand::Constant(if self.registry.is_some() {
                self.tables.eval(g, d, ch)
            } else {
                Q::ZERO
            });
        } else if g == 1 {
            return Operand::Constant(if self.registry.is_some() {
                self.tables.eval(g, d, ch)
            } else {
                self.tables.multinomial(&d) / Q::from(2)
            });
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
            let zero = if self.pressure_view && index < self.base {
                self.engine.source.as_ref().is_some_and(|source| {
                    source
                        .read()
                        .get(index)
                        .is_some_and(|value| *value == Q::ZERO)
                })
            } else {
                index < self.values.len() && self.values[index] == Q::ZERO
            };
            if zero {
                return Operand::Constant(Q::ZERO);
            }
            return Operand::Node(index);
        }
        // 小查询的串行结果也可用于之后的大查询。
        if let Some(value) = self.serial.get(ch).and_then(|values| values.get(&d)) {
            return Operand::Constant(value.clone());
        }
        if let Some(value) = self.shared_serial.and_then(|serial| serial.get(&d, ch)) {
            return Operand::Constant(value);
        }
        // Capacity and future-value credit are obtained before locking the
        // original mathematical index and before publishing any new ID.
        pressure::discovery(&mut self.discovered, self.control);
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

struct BuildJob {
    states: Arc<Vec<State>>,
    start: usize,
    end: usize,
    restored: Option<Vec<Option<Q>>>,
    wave: Option<pressure::Wave>,
}

struct Built {
    discovered: Vec<State>,
    equations: Vec<(usize, (u32, u32), PackedEquation)>,
    terms: TermBlock,
    restored: Vec<usize>,
}

impl Built {
    fn store(
        self,
        plan: &mut Plan,
        layers: &mut BTreeMap<(u32, u32), Vec<usize>>,
        base: usize,
        count: usize,
    ) -> Vec<State> {
        pressure::capacity(&mut plan.equations, count);
        plan.equations.resize_with(count, PackedEquation::default);
        let block = plan.terms.len();
        pressure::capacity(&mut plan.terms, block + 1);
        if let Some(workspace) = &plan.workspace {
            let stored = workspace
                .push_block(self.terms)
                .unwrap_or_else(|error| crate::cache_source::resource_failure(error));
            assert_eq!(stored, block);
            plan.terms.push(TermBlock::default());
        } else {
            plan.terms.push(self.terms);
        }
        let restored_len = plan.restored.len() + self.restored.len();
        pressure::capacity(&mut plan.restored, restored_len);
        plan.restored.extend(self.restored);
        for (id, rank, mut equation) in self.equations {
            equation.block = block;
            plan.equations[id - base] = equation;
            let layer = layers.entry(rank).or_default();
            if layer.len() == layer.capacity() {
                let needed = layer.len() + 1;
                pressure::capacity(layer, needed);
            }
            layer.push(id - base);
        }
        self.discovered
    }
}

impl Planner<'_> {
    fn build(&mut self, mut job: BuildJob) -> Built {
        let wave = job.wave.take();
        pressure::Wave::with(wave.as_ref(), || self.build_unscoped(job))
    }

    fn build_unscoped(&mut self, job: BuildJob) -> Built {
        if job.restored.is_some() {
            return self.build_restored(job);
        }
        let mut equations = Vec::new();
        pressure::capacity(&mut equations, job.end - job.start);
        if equations.capacity() < job.end - job.start {
            equations.reserve_exact(job.end - job.start);
        }
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
            restored: Vec::new(),
        }
    }

    fn build_restored(&mut self, mut job: BuildJob) -> Built {
        let values = job.restored.take().unwrap();
        assert_eq!(values.len(), job.end - job.start);
        let mut equations = Vec::new();
        pressure::capacity(&mut equations, values.len());
        if equations.capacity() < values.len() {
            equations.reserve_exact(values.len());
        }
        let mut restored = Vec::new();
        let mut terms = TermBlock::default();
        let mut scratch = Vec::new();
        for (state, value) in job.states[job.start..job.end].iter().zip(values) {
            let mut equation = if let Some(value) = value {
                if restored.len() == restored.capacity() {
                    let needed = restored.len() + 1;
                    pressure::capacity(&mut restored, needed);
                }
                restored.push(state.index);
                Equation {
                    constant: value,
                    ..Equation::default()
                }
            } else {
                self.expand(state.g, &state.d, &state.ch, scratch)
            };
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
            restored,
        }
    }
}

struct FrontierWork {
    states: Arc<Vec<State>>,
    start: usize,
    observed: bool,
}

#[derive(Default)]
struct Draft {
    plan: Plan,
    layers: BTreeMap<(u32, u32), Vec<usize>>,
    pending: Vec<State>,
    frontier: Option<FrontierWork>,
    // At most one native-refresh retry for a refused boundary. A genuinely
    // oversized wave must progress to backing/shrinking, not refresh forever.
    refreshed: bool,
    wave_cap: usize,
    successful_full_waves: u8,
    retry_since: Option<std::time::Instant>,
}

#[cfg(test)]
thread_local! {
    static RECOVERY_CAP_PEAK: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static REJECT_PREFLIGHT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static ROLLBACK_VISITS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static POOL_UNWIND: std::cell::Cell<Option<PoolUnwind>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
#[derive(Clone, Copy)]
enum PoolUnwind {
    Ordinary,
    Resource,
}

enum Planning {
    Complete(Plan),
    Paused(Draft),
}

impl Draft {
    fn wave_count(&self, remaining: usize) -> usize {
        // Historical lookup has its own bounded read batches. Only actual
        // backing or a rejected wave changes the mathematical frontier size.
        let count = pressure::wave_limit(self.plan.workspace.is_some(), remaining);
        if self.wave_cap == 0 {
            count
        } else {
            count.min(self.wave_cap)
        }
    }

    fn prepare(&mut self) -> usize {
        if self
            .frontier
            .as_ref()
            .is_some_and(|frontier| frontier.start == frontier.states.len())
        {
            self.frontier = None;
        }
        if self.frontier.is_none() && !self.pending.is_empty() {
            if self.plan.workspace.is_some() {
                // A Built block packs contiguous frontier states. Group the
                // backed frontier by the same rank used by evaluation layers,
                // so one block is not reread across many unrelated ranks.
                // Sorting moves existing states in place; IDs stay unchanged.
                self.pending.sort_unstable_by_key(|state| state.rank);
            }
            self.frontier = Some(FrontierWork {
                states: Arc::new(std::mem::take(&mut self.pending)),
                start: 0,
                observed: false,
            });
        }
        self.frontier
            .as_ref()
            .map_or(0, |frontier| frontier.states.len() - frontier.start)
    }

    fn observe(&mut self, planner: &Planner<'_>) {
        let frontier = self.frontier.as_mut().unwrap();
        if frontier.observed {
            return;
        }
        if let Some(source) = &planner.engine.source {
            source.observe(&frontier.states);
        }
        #[cfg(test)]
        measure_identity::observe_states(&frontier.states);
        #[cfg(test)]
        measure::check_memory_limit();
        frontier.observed = true;
    }

    fn store(&mut self, built: Built, planner: &Planner<'_>, wave: Option<&pressure::Wave>) {
        pressure::Wave::with(wave, || {
            let discovered = built.store(
                &mut self.plan,
                &mut self.layers,
                planner.base,
                planner.next.load(Ordering::Relaxed) - planner.base,
            );
            let count = self.pending.len().saturating_add(discovered.len());
            pressure::capacity(&mut self.pending, count);
            self.pending.extend(discovered);
        });
    }

    // No equation, child frontier, block or layer becomes part of this Plan
    // until every job succeeded and every fallible capacity grant succeeded.
    // Preflight may grow capacities/add empty rank entries, which preserves
    // the preceding mathematical Plan if its credit is refused partway.
    fn preflight(&mut self, built: &[Built], planner: &Planner<'_>, wave: Option<&pressure::Wave>) {
        #[cfg(test)]
        if REJECT_PREFLIGHT.with(|reject| reject.replace(false)) {
            pressure::exhaustion(
                crate::cache_source::budget::ReservationError::Busy,
                "test preflight capacity rejection",
            );
        }
        if wave.is_none() {
            return;
        }
        pressure::Wave::with(wave, || {
            let count = planner.next.load(Ordering::Relaxed) - planner.base;
            pressure::capacity(&mut self.plan.equations, count);
            let blocks = self
                .plan
                .terms
                .len()
                .checked_add(built.len())
                .unwrap_or_else(|| {
                    crate::cache_source::resource_failure("plan block count overflow")
                });
            pressure::capacity(&mut self.plan.terms, blocks);
            let mut restored = self.plan.restored.len();
            let mut discovered = self.pending.len();
            let mut ranks = BTreeMap::<(u32, u32), usize>::new();
            pressure::reallocation(built.len().saturating_mul(512));
            if let Some(workspace) = &self.plan.workspace {
                workspace
                    .reserve_blocks(built.len(), pressure::reallocation)
                    .unwrap_or_else(|error| crate::cache_source::resource_failure(error));
            }
            let rank_limit = planner.engine.index.iter().map(Vec::len).sum::<usize>();
            let equation_count = built
                .iter()
                .map(|built| built.equations.len())
                .sum::<usize>();
            pressure::reallocation(rank_limit.min(equation_count).saturating_mul(128));
            for built in built {
                restored = restored
                    .checked_add(built.restored.len())
                    .unwrap_or_else(|| {
                        crate::cache_source::resource_failure("restored state count overflow")
                    });
                discovered = discovered
                    .checked_add(built.discovered.len())
                    .unwrap_or_else(|| {
                        crate::cache_source::resource_failure("frontier count overflow")
                    });
                for equation in &built.equations {
                    *ranks.entry(equation.1).or_default() += 1;
                }
            }
            pressure::capacity(&mut self.plan.restored, restored);
            pressure::capacity(&mut self.pending, discovered);
            for (rank, count) in ranks {
                let layer = self.layers.entry(rank).or_default();
                let needed = layer.len().checked_add(count).unwrap_or_else(|| {
                    crate::cache_source::resource_failure("plan layer length overflow")
                });
                pressure::capacity(layer, needed);
            }
        });
    }

    fn commit(
        &mut self,
        built: &mut Vec<Built>,
        planner: &Planner<'_>,
        wave: Option<&pressure::Wave>,
    ) -> Result<(), Box<dyn std::any::Any + Send>> {
        let preflight = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            self.preflight(built, planner, wave);
        }));
        preflight?;
        for built in built.drain(..) {
            self.store(built, planner, wave);
        }
        Ok(())
    }

    fn failed_wave(
        &mut self,
        planner: &mut Planner<'_>,
        cutoff: usize,
        checkpoint: usize,
        count: usize,
        failure: FailedWave,
    ) -> bool {
        let FailedWave {
            payload,
            replies,
            wave,
        } = failure;
        let error = match payload.downcast::<pressure::WaveExhausted>() {
            Ok(error) => error,
            Err(payload) => std::panic::resume_unwind(payload),
        };
        #[cfg(test)]
        if error.message == "test wave capacity rejection" {
            assert!(
                planner.next.load(Ordering::Relaxed) > cutoff,
                "the regression must reject after assigning a child ID"
            );
        }
        // Every reply has arrived and no job can retain a key/ID from this
        // failed wave. Only the immutable input frontier was observed by the
        // result source; its IDs precede cutoff and remain untouched.
        let rollback_start = REPORT_STAGES
            .get()
            .copied()
            .unwrap_or(false)
            .then(std::time::Instant::now);
        let removed = rollback_discovered(
            planner.engine,
            cutoff,
            replies
                .into_discovered()
                .chain(std::mem::take(&mut planner.discovered)),
        );
        planner.next.store(cutoff, Ordering::Relaxed);
        planner
            .control
            .unwrap()
            .rollback(checkpoint)
            .unwrap_or_else(|error| crate::cache_source::resource_failure(error));
        // All partial equations/keys and their capacities have now died. Keep
        // the wave credit alive until then, before sampling/admitting a retry.
        drop(wave);
        if let Some(start) = rollback_start {
            self.plan.rollback_count += 1;
            self.plan.rollback_states += removed;
            self.plan.rollback_elapsed += start.elapsed();
        }
        // Expired observations are not evidence of physical memory pressure.
        // Leave the borrowed scope, refresh, and retry a smaller resident wave.
        // Capacity rejection still goes directly to the backing slow path.
        self.refreshed = !error.stale;
        self.wave_cap = count.div_ceil(2).max(1);
        self.successful_full_waves = 0;
        let started = *self.retry_since.get_or_insert_with(std::time::Instant::now);
        let Some(workspace) = &self.plan.workspace else {
            return false;
        };
        // We are already using the short, unborrowed view. A smaller wave can
        // retry immediately after sampling; the minimum wave waits only for
        // independent sealed-block spill progress, with a finite deadline.
        if count > 1 {
            pressure::refresh_normal();
            true
        } else if pressure::retry(workspace, started) {
            true
        } else {
            crate::cache_source::resource_failure(format!(
                "minimum plan wave cannot fit: {}",
                error.message
            ));
        }
    }

    fn advance(&mut self, count: usize) {
        self.frontier.as_mut().unwrap().start += count;
        // Refresh deduplication belongs to one boundary without progress,
        // rather than the entire query. Later stale observations may refresh.
        self.refreshed = false;
        self.retry_since = None;
        if self.plan.workspace.is_some() {
            self.recover_wave_cap(count);
        } else {
            // A stale observation can shrink a resident retry as well. Once
            // that retry succeeds, restore the original frontier scheduler;
            // the 4096-state pressure ceiling does not apply to resident work.
            self.wave_cap = 0;
            self.successful_full_waves = 0;
        }
    }

    fn recover_wave_cap(&mut self, count: usize) {
        let ceiling = pressure::wave_limit(true, usize::MAX);
        if self.wave_cap == 0 || self.wave_cap >= ceiling || count < self.wave_cap {
            return;
        }
        // A short frontier tail does not demonstrate room for a larger wave.
        // Four fully committed waves permit a bounded retry at twice the size;
        // any later rejection still halves the cap and resets this evidence.
        self.successful_full_waves += 1;
        if self.successful_full_waves == 4 {
            self.wave_cap = self.wave_cap.saturating_mul(2).min(ceiling);
            self.successful_full_waves = 0;
            #[cfg(test)]
            RECOVERY_CAP_PEAK.with(|peak| peak.set(peak.get().max(self.wave_cap)));
        }
    }

    fn finish(mut self) -> Plan {
        self.plan.layers = self
            .layers
            .into_values()
            .filter(|layer| !layer.is_empty())
            .collect();
        self.plan
    }
}

fn build_plan(planner: &mut Planner<'_>, threads: usize) -> Planning {
    let draft = Draft {
        pending: std::mem::take(&mut planner.discovered),
        ..Draft::default()
    };
    resume_plan(planner, threads, draft)
}

// Failure here occurs only at a frontier/wave boundary with no active jobs.
// The borrowed fast view returns Paused; its caller releases guards and joins
// the worker scope before doing any sampling, spill registration, or waiting.
fn admit_wave(
    planner: &Planner<'_>,
    draft: &Draft,
    remaining: usize,
    threads: usize,
) -> Option<(usize, Option<pressure::Wave>)> {
    if planner.control.is_none() {
        return Some((remaining, None));
    }
    let backed = draft.plan.workspace.is_some();
    let mut count = draft.wave_count(remaining);
    let started = std::time::Instant::now();
    loop {
        let jobs = if count < 128 {
            1
        } else {
            (threads.max(1) * 4).min(count.div_ceil(32))
        };
        let admitted = pressure::Wave::begin(count, jobs, backed).and_then(|wave| {
            if planner.lookup_history {
                // The original response Vec and all per-job Option<Q> Vecs
                // coexist during dispatch. One bounded backend reply may also
                // coexist with the original Vec. Charge all three containers
                // once before lookup, so subsequent term growth cannot consume
                // their credit. Rational payload is moved between these Vecs.
                let spare = jobs.checked_add(1).and_then(|jobs| jobs.checked_mul(4));
                let bytes = count
                    .checked_mul(3)
                    .and_then(|slots| spare.and_then(|spare| slots.checked_add(spare)))
                    .and_then(|slots| slots.checked_mul(std::mem::size_of::<Option<Q>>()))
                    .unwrap_or_else(|| {
                        crate::cache_source::resource_failure("recovery response buffer overflow")
                    });
                let admitted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    pressure::Wave::with(Some(&wave), || pressure::reallocation(bytes));
                }));
                if let Err(payload) = admitted {
                    return match payload.downcast::<pressure::WaveExhausted>() {
                        Ok(error) => Err(error.message),
                        Err(payload) => std::panic::resume_unwind(payload),
                    };
                }
            }
            Ok(wave)
        });
        match admitted {
            Ok(wave) => return Some((count, Some(wave))),
            Err(_) if !backed => return None,
            Err(error) => {
                if count > 1 {
                    count = count.div_ceil(2);
                    continue;
                }
                if pressure::retry(draft.plan.workspace.as_ref().unwrap(), started) {
                    continue;
                }
                crate::cache_source::resource_failure(format!(
                    "plan frontier cannot make progress: {error}"
                ));
            }
        }
    }
}

fn restored_wave(
    planner: &Planner<'_>,
    states: &Arc<Vec<State>>,
    start: usize,
    count: usize,
) -> Option<std::vec::IntoIter<Option<Q>>> {
    planner
        .lookup_history
        .then(|| recovery::frontier_range(planner.registry.unwrap(), states, start..start + count))
}

struct FailedBuild {
    payload: Box<dyn std::any::Any + Send>,
    discovered: Vec<State>,
}

type BuildResult = Result<Built, FailedBuild>;

#[derive(Default)]
struct BuildReplies {
    staged: Vec<Built>,
    // These are the original buffers moved out of failed workers, not a second
    // per-node log. The outer container allocates only when a job has failed.
    partial: Vec<Vec<State>>,
    failure: Option<Box<dyn std::any::Any + Send>>,
}

impl BuildReplies {
    fn with_capacity(jobs: usize) -> Self {
        Self {
            staged: Vec::with_capacity(jobs),
            ..Self::default()
        }
    }

    fn commit(
        &mut self,
        draft: &mut Draft,
        planner: &Planner<'_>,
        wave: Option<&pressure::Wave>,
    ) -> Result<(), Box<dyn std::any::Any + Send>> {
        if let Some(payload) = self.failure.take() {
            return Err(payload);
        }
        draft.commit(&mut self.staged, planner, wave)
    }

    fn into_discovered(self) -> impl Iterator<Item = State> {
        self.staged
            .into_iter()
            .flat_map(|built| built.discovered)
            .chain(self.partial.into_iter().flatten())
    }
}

struct FailedWave {
    payload: Box<dyn std::any::Any + Send>,
    replies: BuildReplies,
    wave: Option<pressure::Wave>,
}

// Jobs have all replied before this is called. Their existing discovery lists
// identify every insertion, so old mathematical states are never traversed.
fn rollback_discovered(
    engine: &Engine,
    cutoff: usize,
    discovered: impl Iterator<Item = State>,
) -> usize {
    let mut removed = 0;
    for state in discovered {
        #[cfg(test)]
        ROLLBACK_VISITS.with(|visits| visits.set(visits.get() + 1));
        if state.index < cutoff {
            continue;
        }
        let shard = &engine.index[state.rank.0 as usize][state.rank.1 as usize];
        let mut index = shard.write().unwrap_or_else(|error| error.into_inner());
        if let Some(values) = index.get_mut(&state.ch) {
            // Checking both the full identity and its exact assigned ID keeps
            // even a stale/duplicate record from removing a different state.
            if values.get(&state.d) == Some(&state.index) {
                values.remove(&state.d);
                removed += 1;
            }
            if values.is_empty() {
                index.remove(&state.ch);
            }
        }
        drop(index);
        shard.clear_poison();
    }
    removed
}

fn build_caught(planner: &mut Planner<'_>, job: BuildJob) -> BuildResult {
    catch_build(planner, |planner| planner.build(job))
}

fn catch_build<'a>(
    planner: &mut Planner<'a>,
    build: impl FnOnce(&mut Planner<'a>) -> Built,
) -> BuildResult {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| build(planner))).map_err(|payload| {
        FailedBuild {
            payload,
            // Preserve partial work for exact rollback, while leaving the
            // worker no buffer whose capacity outlives returned ID credit.
            discovered: std::mem::take(&mut planner.discovered),
        }
    })
}

fn stage_reply(reply: BuildResult, replies: &mut BuildReplies) {
    match reply {
        Ok(built) => replies.staged.push(built),
        Err(FailedBuild {
            payload,
            discovered,
        }) => {
            replies.partial.push(discovered);
            // Never hide an ordinary mathematical/programming panic behind a
            // simultaneous recoverable capacity rejection from another job.
            let replace = replies.failure.as_ref().is_none_or(|old| {
                old.is::<pressure::WaveExhausted>() && !payload.is::<pressure::WaveExhausted>()
            });
            if replace {
                replies.failure = Some(payload);
            }
        }
    }
}

fn resume_plan(planner: &mut Planner<'_>, threads: usize, mut draft: Draft) -> Planning {
    // Preserve the original no-pool path for small frontiers and warm misses.
    loop {
        let remaining = draft.prepare();
        if remaining == 0 {
            return Planning::Complete(draft.finish());
        }
        if remaining >= 128 {
            break;
        }
        let Some((count, wave)) = admit_wave(planner, &draft, remaining, 1) else {
            return Planning::Paused(draft);
        };
        draft.observe(planner);
        let frontier = draft.frontier.as_ref().unwrap();
        let start = frontier.start;
        let states = Arc::clone(&frontier.states);
        let cutoff = planner.next.load(Ordering::Relaxed);
        let checkpoint = planner.control.map_or(0, pressure::Control::checkpoint);
        let mut restored = restored_wave(planner, &states, start, count);
        let result = build_caught(
            planner,
            BuildJob {
                states,
                start,
                end: start + count,
                restored: recovery::take(&mut restored, count),
                wave: wave.clone(),
            },
        );
        let mut replies = BuildReplies::with_capacity(1);
        stage_reply(result, &mut replies);
        if let Err(payload) = replies.commit(&mut draft, planner, wave.as_ref()) {
            if draft.failed_wave(
                planner,
                cutoff,
                checkpoint,
                count,
                FailedWave {
                    payload,
                    replies,
                    wave,
                },
            ) {
                continue;
            }
            return Planning::Paused(draft);
        }
        draft.advance(count);
        if let Some(wave) = wave {
            wave.finish();
        }
    }
    let threads = threads.max(1);
    let (jobs_tx, jobs_rx) = mpsc::channel::<BuildJob>();
    let jobs_rx = Mutex::new(jobs_rx);
    let (reply_tx, reply_rx) = mpsc::channel();
    let outcome = std::thread::scope(|scope| {
        // Own the only job sender in this callback's stack frame. Unwinding
        // must disconnect waiting workers before scope attempts to join them.
        let jobs_tx = jobs_tx;
        for _ in 1..threads {
            let mut worker = planner.worker();
            let reply_tx = reply_tx.clone();
            let jobs_rx = &jobs_rx;
            scope.spawn(move || {
                loop {
                    let Ok(job) = jobs_rx.lock().unwrap().recv() else {
                        break;
                    };
                    let result = build_caught(&mut worker, job);
                    if reply_tx.send(result).is_err() {
                        break;
                    }
                }
            });
        }
        drop(reply_tx);
        #[cfg(test)]
        if let Some(unwind) = POOL_UNWIND.with(|unwind| unwind.take()) {
            match unwind {
                PoolUnwind::Ordinary => {
                    std::panic::resume_unwind(Box::new("test planning controller panic"));
                }
                PoolUnwind::Resource => {
                    crate::cache_source::resource_failure(
                        "test planning controller resource error",
                    );
                }
            }
        }
        let result: Result<bool, Box<dyn std::any::Any + Send>> = (|| {
            loop {
                let remaining = draft.prepare();
                if remaining == 0 {
                    return Ok(true);
                }
                let Some((count, wave)) = admit_wave(planner, &draft, remaining, threads) else {
                    return Ok(false);
                };
                draft.observe(planner);
                let frontier = draft.frontier.as_ref().unwrap();
                let states = Arc::clone(&frontier.states);
                let offset = frontier.start;
                let cutoff = planner.next.load(Ordering::Relaxed);
                let checkpoint = planner.control.map_or(0, pressure::Control::checkpoint);
                let mut restored = restored_wave(planner, &states, offset, count);
                let chunk = if count < 128 {
                    count
                } else {
                    count.div_ceil(threads * 4).max(32)
                };
                let first_restored = recovery::take(&mut restored, chunk.min(count));
                let mut jobs = 0;
                for start in (chunk..count).step_by(chunk) {
                    jobs_tx
                        .send(BuildJob {
                            states: Arc::clone(&states),
                            start: offset + start,
                            end: offset + (start + chunk).min(count),
                            restored: recovery::take(&mut restored, chunk.min(count - start)),
                            wave: wave.clone(),
                        })
                        .unwrap();
                    jobs += 1;
                }
                let main = build_caught(
                    planner,
                    BuildJob {
                        states,
                        start: offset,
                        end: offset + chunk.min(count),
                        restored: first_restored,
                        wave: wave.clone(),
                    },
                );
                let mut replies = BuildReplies::with_capacity(jobs + 1);
                stage_reply(main, &mut replies);
                // Drain every dispatched job before touching next/index or
                // returning Paused. No worker can assign an ID during rollback.
                for _ in 0..jobs {
                    stage_reply(reply_rx.recv().unwrap(), &mut replies);
                }
                if let Err(payload) = replies.commit(&mut draft, planner, wave.as_ref()) {
                    if draft.failed_wave(
                        planner,
                        cutoff,
                        checkpoint,
                        count,
                        FailedWave {
                            payload,
                            replies,
                            wave,
                        },
                    ) {
                        continue;
                    }
                    return Ok(false);
                }
                draft.advance(count);
                if let Some(wave) = wave {
                    wave.finish();
                }
            }
        })();
        drop(jobs_tx);
        result
    });
    match outcome {
        Ok(true) => Planning::Complete(draft.finish()),
        Ok(false) => Planning::Paused(draft),
        Err(payload) => std::panic::resume_unwind(payload),
    }
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

// Same scheduling and arithmetic as evaluate_layers. The existing table lock
// belongs to a long-lived source, and publication occurs once per whole layer.
enum Evaluation {
    Complete,
    PausedEvaluation {
        next_layer: usize,
        outputs: pressure::OutputMemory,
    },
}

#[inline]
fn accumulate_completed_heap(value: &Q, total: &mut u64, peak: &mut u64) {
    let heap = pressure::OutputMemory::heap(value);
    *total = (*total)
        .checked_add(heap)
        .unwrap_or_else(|| crate::cache_source::resource_failure("layer output heap overflow"));
    *peak = (*peak).max(heap);
}

fn evaluate_shared(
    plan: &Plan,
    source: &ParallelSource,
    batch: &Batch,
    base: usize,
    threads: usize,
) -> Evaluation {
    let layers = &batch.layers;
    let mut outputs = pressure::OutputMemory::default();
    let width = layers.iter().map(Vec::len).max().unwrap_or(0);
    let threads = if width < 128 {
        1
    } else {
        threads.min(width.div_ceil(32)).max(1)
    };
    if threads == 1 {
        for (layer_index, layer) in layers.iter().enumerate() {
            let Some(credit) = pressure::admit_layer(outputs, layer.len(), 1, false) else {
                return Evaluation::PausedEvaluation {
                    next_layer: layer_index,
                    outputs,
                };
            };
            let (heap_bytes, peak_heap) = {
                let mut values = source.write();
                let (mut heap_bytes, mut peak_heap) = (0, 0);
                for &id in layer {
                    let value = plan.evaluate(id, &values);
                    accumulate_completed_heap(&value, &mut heap_bytes, &mut peak_heap);
                    values[base + id] = value;
                }
                (heap_bytes, peak_heap)
            };
            source.publish(batch, layer_index);
            outputs.finish(credit, layer.len(), 1, false, heap_bytes, peak_heap, 0);
        }
        return Evaluation::Complete;
    }
    let (jobs_tx, jobs_rx) = mpsc::channel::<Job>();
    let jobs_rx = Mutex::new(jobs_rx);
    let (reply_tx, reply_rx) = mpsc::channel();
    let calculate = |job: Job| {
        let values = source.read();
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
                let buffered = layer.len() >= 128;
                let workers = if buffered { threads } else { 1 };
                // The previous layer's replies and source guards have ended.
                // Rejection leaves this entire layer untouched; the worker
                // scope is joined before the caller can migrate Plan.terms.
                let Some(credit) = pressure::admit_layer(outputs, layer.len(), workers, buffered)
                else {
                    return Ok(Evaluation::PausedEvaluation {
                        next_layer: layer_index,
                        outputs,
                    });
                };
                if layer.len() < 128 {
                    let (heap_bytes, peak_heap) = {
                        let mut values = source.write();
                        let (mut heap_bytes, mut peak_heap) = (0, 0);
                        for &id in layer {
                            let value = plan.evaluate(id, &values);
                            accumulate_completed_heap(&value, &mut heap_bytes, &mut peak_heap);
                            values[base + id] = value;
                        }
                        (heap_bytes, peak_heap)
                    };
                    source.publish(batch, layer_index);
                    outputs.finish(credit, layer.len(), 1, false, heap_bytes, peak_heap, 0);
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
                let (heap_bytes, peak_heap) = {
                    let mut values = source.write();
                    let (mut heap_bytes, mut peak_heap) = (0, 0);
                    for (id, value) in completed {
                        accumulate_completed_heap(&value, &mut heap_bytes, &mut peak_heap);
                        values[base + id] = value;
                    }
                    (heap_bytes, peak_heap)
                };
                source.publish(batch, layer_index);
                outputs.finish(credit, layer.len(), workers, true, heap_bytes, peak_heap, 0);
            }
            Ok(Evaluation::Complete)
        })();
        drop(jobs_tx);
        result
    });
    match outcome {
        Ok(evaluation) => evaluation,
        Err(payload) => std::panic::resume_unwind(payload),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn engine_use_permit_survives_table_borrows_and_releases_on_unwind() {
        use crate::cache_source::ResultSource;
        for panic in [false, true] {
            let registry = Arc::new(Registry::new());
            let source = with_engine(Some(&registry), |engine| {
                let source = Arc::clone(engine.source.as_ref().unwrap());
                engine.index = vec![
                    (0..5)
                        .map(|_| Arc::new(RwLock::new(Index::new())))
                        .collect(),
                ];
                let shard = Arc::clone(&engine.index[0][4]);
                shard
                    .write()
                    .unwrap()
                    .entry(vec![])
                    .or_default()
                    .insert(vec![2, 3], 0);
                source.observe_shards(vec![shard]);
                source.write().push(Q::from(7));
                source.observe(&Arc::new(vec![State {
                    g: 2,
                    d: vec![2, 3],
                    ch: Vec::new(),
                    rank: (0, 4),
                    index: 0,
                }]));
                let batch = source.begin(0, vec![vec![0]]);
                source.publish(&batch, 0);
                source.durable(1);
                source
            });
            assert_eq!(cached(2, &[2, 3], &[], &registry), Some(Q::from(7)));
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                with_engine(Some(&registry), |engine| {
                    let id = engine.index[0][4].read().unwrap()[&vec![]][&vec![2, 3]];
                    let reclaimed = std::thread::scope(|scope| {
                        scope
                            .spawn(|| source.reclaim_during_query(&registry))
                            .join()
                            .unwrap()
                    });
                    assert_eq!(reclaimed, 0);
                    assert!(!registry.has_history());
                    if panic {
                        std::panic::resume_unwind(Box::new("test engine use panic"));
                    }
                    // The final owned Q is produced before with_engine drops
                    // its permit. No raw ID survives the mathematical scope.
                    source.read()[id].clone()
                })
            }));
            let value = if panic {
                assert_eq!(
                    result.unwrap_err().downcast_ref::<&str>(),
                    Some(&"test engine use panic")
                );
                None
            } else {
                Some(result.unwrap())
            };
            assert!(source.reclaim_during_query(&registry) > 0);
            assert!(registry.has_history());
            assert!(source.read().is_empty());
            assert!(cached(2, &[2, 3], &[], &registry).is_none());
            if let Some(value) = value {
                assert_eq!(value, Q::from(7));
            }
        }
    }

    fn rollback_engine() -> Engine {
        Engine {
            index: (0..7)
                .map(|_| {
                    (0..9)
                        .map(|_| Arc::new(RwLock::new(Index::new())))
                        .collect()
                })
                .collect(),
            ..Engine::default()
        }
    }

    fn discovered_build(planner: &mut Planner<'_>) -> Built {
        Built {
            discovered: std::mem::take(&mut planner.discovered),
            equations: Vec::new(),
            terms: TermBlock::default(),
            restored: Vec::new(),
        }
    }

    fn reject_test_build() -> ! {
        pressure::exhaustion(
            crate::cache_source::budget::ReservationError::Busy,
            "test worker capacity rejection",
        )
    }

    #[test]
    fn rejected_worker_returns_its_partial_discoveries_for_exact_rollback() {
        let engine = rollback_engine();
        let next = AtomicUsize::new(0);
        let reply = std::thread::scope(|scope| {
            scope
                .spawn(|| {
                    let mut worker = Planner::new(&engine, &EMPTY_SERIAL, &[], &next, None);
                    let reply = catch_build(&mut worker, |worker| {
                        assert!(matches!(worker.node(2, vec![2, 3], &[]), Operand::Node(0)));
                        reject_test_build();
                    });
                    assert_eq!(worker.discovered.capacity(), 0);
                    reply
                })
                .join()
                .unwrap()
        });
        let Err(failed) = reply else {
            panic!("the worker must reject after publishing a child ID");
        };
        assert!(failed.payload.is::<pressure::WaveExhausted>());
        assert_eq!(failed.discovered.len(), 1);
        assert_eq!(next.load(Ordering::Relaxed), 1);
        assert_eq!(
            rollback_discovered(&engine, 0, failed.discovered.into_iter()),
            1
        );
        assert!(engine.index[0][4].read().unwrap().is_empty());
    }

    #[test]
    fn failed_wave_keeps_successful_and_failed_worker_discoveries() {
        let engine = rollback_engine();
        let next = AtomicUsize::new(0);
        let mut planner = Planner::new(&engine, &EMPTY_SERIAL, &[], &next, None);
        assert!(matches!(
            planner.node(3, vec![], &[(1, 6)]),
            Operand::Node(0)
        ));
        let history = std::mem::take(&mut planner.discovered);
        let cutoff = next.load(Ordering::Relaxed);
        let (successful, failed) = std::thread::scope(|scope| {
            let successful = scope.spawn(|| {
                let mut worker = planner.worker();
                catch_build(&mut worker, |worker| {
                    worker.node(2, vec![2, 3], &[]);
                    discovered_build(worker)
                })
            });
            let failed = scope.spawn(|| {
                let mut worker = planner.worker();
                catch_build(&mut worker, |worker| {
                    worker.node(2, vec![2, 2, 2], &[]);
                    reject_test_build();
                })
            });
            (successful.join().unwrap(), failed.join().unwrap())
        });
        let mut replies = BuildReplies::with_capacity(2);
        stage_reply(successful, &mut replies);
        stage_reply(failed, &mut replies);
        assert_eq!(replies.staged.len(), 1);
        assert_eq!(replies.staged[0].discovered.len(), 1);
        assert_eq!(replies.partial.len(), 1);
        assert_eq!(replies.partial[0].len(), 1);
        let mut draft = Draft::default();
        let error = replies.commit(&mut draft, &planner, None).unwrap_err();
        assert!(error.is::<pressure::WaveExhausted>());
        assert!(draft.plan.equations.is_empty());
        assert_eq!(
            rollback_discovered(&engine, cutoff, replies.into_discovered()),
            2
        );
        assert_eq!(
            engine.index[6][4].read().unwrap()[&history[0].ch][&history[0].d],
            0
        );
        assert!(engine.index[0][4].read().unwrap().is_empty());
        assert!(engine.index[0][5].read().unwrap().is_empty());
    }

    #[test]
    fn preflight_rejection_retains_the_original_staged_discovery_buffer() {
        let engine = rollback_engine();
        let next = AtomicUsize::new(0);
        let mut planner = Planner::new(&engine, &EMPTY_SERIAL, &[], &next, None);
        let reply = catch_build(&mut planner, |planner| {
            planner.node(2, vec![2, 3], &[]);
            discovered_build(planner)
        });
        let mut replies = BuildReplies::with_capacity(1);
        stage_reply(reply, &mut replies);
        let allocation = replies.staged[0].discovered.as_ptr();
        REJECT_PREFLIGHT.with(|reject| reject.set(true));
        let mut draft = Draft::default();
        let error = replies.commit(&mut draft, &planner, None).unwrap_err();
        assert!(error.is::<pressure::WaveExhausted>());
        assert_eq!(replies.staged[0].discovered.as_ptr(), allocation);
        assert_eq!(replies.staged[0].discovered.len(), 1);
        assert!(draft.pending.is_empty());
        assert_eq!(
            rollback_discovered(&engine, 0, replies.into_discovered()),
            1
        );
        assert!(engine.index[0][4].read().unwrap().is_empty());
    }

    #[test]
    fn rollback_visits_only_failed_discoveries_and_checks_the_exact_id() {
        let engine = Engine {
            index: vec![
                (0..343)
                    .map(|_| Arc::new(RwLock::new(Index::new())))
                    .collect(),
            ],
            ..Engine::default()
        };
        let shard = &engine.index[0][342];
        {
            let mut index = shard.write().unwrap();
            let values = index.entry(vec![]).or_default();
            for id in 0..=256 {
                values.insert(vec![id as u32, 512 - id as u32], id);
            }
        }
        let records = [(256, 256), (100, 100), (101, 999)]
            .into_iter()
            .map(|(degree, index)| State {
                g: 171,
                d: vec![degree, 512 - degree],
                ch: vec![],
                rank: (0, 342),
                index,
            });
        ROLLBACK_VISITS.with(|visits| visits.set(0));
        assert_eq!(rollback_discovered(&engine, 256, records), 1);
        ROLLBACK_VISITS.with(|visits| assert_eq!(visits.get(), 3));
        let index = shard.read().unwrap();
        let values = index.get(&Vec::new()).unwrap();
        assert_eq!(values.len(), 256);
        for id in 0..256 {
            assert_eq!(values[&vec![id as u32, 512 - id as u32]], id);
        }
    }

    #[test]
    fn ordinary_worker_panic_takes_priority_and_is_not_retried() {
        let engine = rollback_engine();
        let next = AtomicUsize::new(7);
        let mut planner = Planner::new(&engine, &EMPTY_SERIAL, &[], &next, None);
        for ordinary_first in [false, true] {
            let mut replies = BuildReplies::default();
            for ordinary in [ordinary_first, !ordinary_first] {
                let payload: Box<dyn std::any::Any + Send> = if ordinary {
                    Box::new("ordinary mathematical panic")
                } else {
                    Box::new(pressure::WaveExhausted {
                        message: "test capacity rejection".into(),
                        stale: false,
                    })
                };
                stage_reply(
                    Err(FailedBuild {
                        payload,
                        discovered: Vec::new(),
                    }),
                    &mut replies,
                );
            }
            let mut draft = Draft::default();
            let payload = replies.commit(&mut draft, &planner, None).unwrap_err();
            let failure = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                draft.failed_wave(
                    &mut planner,
                    0,
                    0,
                    2,
                    FailedWave {
                        payload,
                        replies,
                        wave: None,
                    },
                )
            }))
            .unwrap_err();
            assert_eq!(
                failure.downcast_ref::<&str>(),
                Some(&"ordinary mathematical panic")
            );
            assert_eq!(next.load(Ordering::Relaxed), 7);
        }
    }

    #[test]
    fn planning_pool_unwind_closes_jobs_before_join() {
        const CHILD: &str = "SC_FABER_TEST_PLANNING_POOL_UNWIND";
        if let Some(kind) = std::env::var_os(CHILD) {
            let ordinary = match kind.to_str() {
                Some("ordinary") => true,
                Some("resource") => false,
                _ => panic!("invalid planning pool test mode"),
            };
            let engine = rollback_engine();
            let next = AtomicUsize::new(128);
            let mut planner = Planner::new(&engine, &EMPTY_SERIAL, &[], &next, None);
            // Only establish a pooled frontier. The controller unwinds after
            // spawning its real workers, before expanding any of these states.
            let draft = Draft {
                pending: (0..128)
                    .map(|index| State {
                        g: 2,
                        d: vec![2, 3],
                        ch: Vec::new(),
                        rank: (0, 4),
                        index,
                    })
                    .collect(),
                ..Draft::default()
            };
            POOL_UNWIND.with(|unwind| {
                unwind.set(Some(if ordinary {
                    PoolUnwind::Ordinary
                } else {
                    PoolUnwind::Resource
                }));
            });
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                resume_plan(&mut planner, 4, draft)
            }));
            let Err(payload) = result else {
                panic!("the planning controller must unwind");
            };
            assert!(POOL_UNWIND.with(|unwind| unwind.get().is_none()));
            if ordinary {
                assert_eq!(
                    payload.downcast_ref::<&str>(),
                    Some(&"test planning controller panic")
                );
            } else {
                let error = payload
                    .downcast_ref::<crate::cache_source::ResourceFailure>()
                    .expect("resource errors must keep their public error type");
                assert_eq!(error.0, "test planning controller resource error");
            }
            return;
        }

        // A deadlock regression must not leave a blocked thread in this test
        // process. Each scenario owns one child and kills/reaps only that child
        // if the pool fails to exit within the bounded wait.
        struct OwnedChild(std::process::Child);
        impl Drop for OwnedChild {
            fn drop(&mut self) {
                let _ = self.0.kill();
                let _ = self.0.wait();
            }
        }
        for kind in ["ordinary", "resource"] {
            let mut child = OwnedChild(
                std::process::Command::new(std::env::current_exe().unwrap())
                    .args([
                        "--exact",
                        "hodge::psi::parallel::tests::planning_pool_unwind_closes_jobs_before_join",
                        "--test-threads=1",
                    ])
                    .env(CHILD, kind)
                    .stdout(std::process::Stdio::null())
                    .spawn()
                    .unwrap(),
            );
            let started = std::time::Instant::now();
            loop {
                if let Some(status) = child.0.try_wait().unwrap() {
                    assert!(
                        status.success(),
                        "{kind} planning pool test failed: {status}"
                    );
                    break;
                }
                assert!(
                    started.elapsed() < std::time::Duration::from_secs(10),
                    "{kind} planning pool did not disconnect and join its workers"
                );
                std::thread::sleep(std::time::Duration::from_millis(10));
            }
        }
    }

    #[test]
    fn resident_history_keeps_full_frontiers_and_backed_waves_keep_their_limit() {
        let registry = Arc::new(Registry::new());
        let engine = rollback_engine();
        let next = AtomicUsize::new(0);
        let planner = Planner::recovering(&engine, None, 0, &next, Some(&registry));
        assert!(planner.lookup_history);
        let backed_limit = pressure::wave_limit(true, usize::MAX);
        let remaining = backed_limit * 3 + 17;
        let mut draft = Draft::default();
        assert_eq!(draft.wave_count(remaining), remaining);
        draft.wave_cap = 137;
        assert_eq!(draft.wave_count(remaining), 137);
        draft.wave_cap = 0;
        draft.plan.workspace = Some(workspace::PlanWorkspace::new(Vec::new()));
        assert_eq!(draft.wave_count(remaining), backed_limit);
        draft.wave_cap = 137;
        assert_eq!(draft.wave_count(remaining), 137);
        assert_eq!(draft.wave_count(19), 19);
    }

    #[test]
    fn resident_history_preserves_exact_c_hits_zero_and_normal_miss_expansion() {
        let registry = Arc::new(Registry::new());
        let engine = rollback_engine();
        let next = AtomicUsize::new(3);
        let mut planner = Planner::recovering(&engine, None, 0, &next, Some(&registry));
        let states = Arc::new(vec![
            State {
                g: 2,
                d: vec![],
                ch: vec![(1, 3)],
                rank: (3, 2),
                index: 0,
            },
            State {
                g: 2,
                d: vec![0, 0],
                ch: vec![(5, 1)],
                rank: (5, 4),
                index: 1,
            },
            State {
                g: 2,
                d: vec![0, 2, 4],
                ch: vec![],
                rank: (0, 5),
                index: 2,
            },
        ]);
        let expected = Descendents::default().eval(2, vec![], &[(1, 3)]);
        let reply = build_caught(
            &mut planner,
            BuildJob {
                states,
                start: 0,
                end: 3,
                restored: Some(vec![Some(expected.clone()), Some(Q::ZERO), None]),
                wave: None,
            },
        );
        let Ok(built) = reply else {
            panic!("resident historical equations must build normally");
        };
        assert_eq!(built.restored, vec![0, 1]);
        for (position, expected) in [expected, Q::ZERO].into_iter().enumerate() {
            let equation = &built.equations[position].2;
            assert_eq!(
                built.terms.evaluate(
                    equation.terms.clone(),
                    &equation.constant,
                    equation.divisor.as_ref(),
                    &[],
                ),
                expected
            );
        }
        assert!(
            !built.equations[2].2.terms.is_empty(),
            "a miss must still expand its ordinary mathematical dependencies"
        );
        assert!(!built.discovered.is_empty());
    }

    #[test]
    fn resident_progress_restores_original_frontier_scheduling() {
        let mut draft = Draft {
            wave_cap: 1,
            successful_full_waves: 3,
            pending: vec![State {
                g: 2,
                d: vec![4],
                ch: vec![],
                rank: (0, 3),
                index: 0,
            }],
            ..Default::default()
        };
        assert_eq!(draft.prepare(), 1);
        draft.advance(1);
        assert_eq!(draft.wave_cap, 0);
        assert_eq!(draft.successful_full_waves, 0);
    }

    #[test]
    fn failed_wave_cap_recovers_from_full_successes_without_unbounded_growth() {
        let mut draft = Draft {
            wave_cap: 8,
            ..Draft::default()
        };
        for _ in 0..4 {
            draft.recover_wave_cap(4);
        }
        assert_eq!(draft.wave_cap, 8, "short tails are not full-wave progress");
        for _ in 0..3 {
            draft.recover_wave_cap(8);
        }
        assert_eq!(draft.wave_cap, 8);
        draft.recover_wave_cap(8);
        assert_eq!(draft.wave_cap, 16);
        let ceiling = pressure::wave_limit(true, usize::MAX);
        while draft.wave_cap < ceiling {
            let old = draft.wave_cap;
            for _ in 0..4 {
                draft.recover_wave_cap(old);
            }
            assert!(draft.wave_cap > old && draft.wave_cap <= ceiling);
        }
        for _ in 0..8 {
            draft.recover_wave_cap(ceiling);
        }
        assert_eq!(draft.wave_cap, ceiling);
        let mut normal = Draft::default();
        for _ in 0..8 {
            normal.recover_wave_cap(ceiling);
        }
        assert_eq!(
            normal.wave_cap, 0,
            "ordinary unrestricted waves stay unchanged"
        );
    }

    #[test]
    fn only_new_backed_frontiers_are_sorted_by_the_evaluation_rank() {
        fn pending() -> Vec<State> {
            [
                (3, vec![], vec![(1, 6)]),
                (2, vec![2, 3], vec![]),
                (2, vec![1], vec![(1, 3)]),
                (2, vec![0, 2], vec![(1, 3)]),
                (2, vec![2, 2, 2], vec![]),
            ]
            .into_iter()
            .enumerate()
            .map(|(index, (g, d, ch))| State {
                g,
                rank: (
                    ch.iter().map(|&(m, power)| m * power).sum(),
                    2 * g - 2 + d.len() as u32,
                ),
                d,
                ch,
                index,
            })
            .collect()
        }

        for backed in [false, true] {
            let mut draft = Draft {
                pending: pending(),
                ..Default::default()
            };
            let allocation = draft.pending.as_ptr();
            if backed {
                draft.plan.workspace = Some(workspace::PlanWorkspace::new(Vec::new()));
            }
            assert_eq!(draft.prepare(), 5);
            let frontier = draft.frontier.as_ref().unwrap();
            assert_eq!(
                frontier.states.as_ptr(),
                allocation,
                "frontier must move without cloning"
            );
            let ids: Vec<_> = frontier.states.iter().map(|state| state.index).collect();
            assert_eq!(
                ids,
                if backed {
                    vec![1, 4, 2, 3, 0]
                } else {
                    vec![0, 1, 2, 3, 4]
                }
            );
            for state in frontier.states.iter() {
                assert_eq!(
                    state.rank,
                    (
                        state.ch.iter().map(|&(m, power)| m * power).sum(),
                        2 * state.g - 2 + state.d.len() as u32,
                    )
                );
            }
            // Entering pressure while an Arc frontier is partially consumed
            // never rearranges IDs already observed by the result consumer.
            draft.frontier.as_mut().unwrap().start = 2;
            if !backed {
                draft.plan.workspace = Some(workspace::PlanWorkspace::new(Vec::new()));
            }
            assert_eq!(draft.prepare(), 3);
            let frontier = draft.frontier.as_ref().unwrap();
            assert_eq!(frontier.states.as_ptr(), allocation);
            assert_eq!(
                frontier
                    .states
                    .iter()
                    .map(|state| state.index)
                    .collect::<Vec<_>>(),
                ids
            );
        }
    }

    #[test]
    fn repeated_stale_observations_after_progress_keep_resident_planning() {
        use crate::cache_source::ResultSource;
        let poly = [(vec![(1, 6)], Q::ONE)];
        let expected = Descendents::default().eval(3, vec![], &poly[0].0);
        // Exercise both boundary admission and an in-wave rejection after a
        // child ID was assigned. No storage exists: accidental backing fails.
        for at_capacity in [false, true] {
            let registry = Arc::new(Registry::new());
            let mut tables = Descendents::with_registry(Some(&registry));
            let (actual, injected, backed) = pressure::stale_across_waves(at_capacity, || {
                integrate(
                    3,
                    &[(vec![], Q::ONE)],
                    &poly,
                    &mut tables,
                    4,
                    Some(&registry),
                )
            });
            assert_eq!(actual, expected);
            assert_eq!(injected, 2, "a committed wave separates the stale events");
            assert_eq!(backed, 0, "stale observations are not memory pressure");
            let source = with_engine(Some(&registry), |engine| {
                engine.source.as_ref().unwrap().clone()
            });
            assert!(source.ready() > 0);
            assert_eq!(
                source.ready(),
                source.read().len(),
                "stale rollback must not leave ID holes"
            );
        }
    }

    #[test]
    fn backed_plan_keeps_exact_descendants_and_fresh_mode_avoids_history_reads() {
        use crate::cache_source::{ResultSource, storage::Storage};
        struct TempCache(std::path::PathBuf);
        impl Drop for TempCache {
            fn drop(&mut self) {
                for stem in ["cache", "capacity"] {
                    for suffix in [".sqlite", ".sqlite-wal", ".sqlite-shm"] {
                        let _ = std::fs::remove_file(self.0.join(format!("{stem}{suffix}")));
                    }
                }
                // This is deliberately nonrecursive: unrelated files in the
                // uniquely created test directory are never removed.
                let _ = std::fs::remove_dir(&self.0);
            }
        }
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "sc-faber-pressure-test-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&directory).unwrap();
        let temporary = TempCache(directory);
        let poly = [(vec![(1, 6)], Q::ONE)];
        let expected = Descendents::default().eval(3, vec![], &poly[0].0);
        for (stem, reject_capacity) in [("cache", false), ("capacity", true)] {
            let path = temporary.0.join(format!("{stem}.sqlite"));
            let registry = Arc::new(Registry::new());
            let storage = Storage::open(path.to_str().unwrap(), Arc::clone(&registry)).unwrap();
            registry.bind_storage(storage.handle());
            let mut tables = Descendents::with_registry(Some(&registry));
            RECOVERY_CAP_PEAK.with(|peak| peak.set(0));
            let calculate = || {
                integrate(
                    3,
                    &[(vec![], Q::ONE)],
                    &poly,
                    &mut tables,
                    4,
                    Some(&registry),
                )
            };
            let (actual, backed) = if reject_capacity {
                pressure::exhaust_next_wave(calculate)
            } else {
                pressure::force_backing(calculate)
            };
            assert_eq!(actual, expected);
            assert_eq!(backed, 1);
            if reject_capacity {
                RECOVERY_CAP_PEAK.with(|peak| {
                    assert!(
                        peak.get() > 1,
                        "real committed equations must recover from cap one"
                    );
                    assert!(peak.get() <= pressure::wave_limit(true, usize::MAX));
                });
            }
            let report = storage.flush().unwrap();
            assert_eq!(
                report.read_keys, 0,
                "fresh pressure is not historical recovery"
            );
            let source = with_engine(Some(&registry), |engine| {
                engine.source.as_ref().unwrap().clone()
            });
            assert!(source.ready() > 0);
            assert_eq!(
                source.ready(),
                source.read().len(),
                "rollback must not leave assigned-ID holes"
            );
            let mut cursor = source.floor();
            let mut seen = BTreeMap::new();
            while cursor < source.ready() {
                if let Some(batch) = source.pull(cursor, 256) {
                    for row in batch.rows {
                        let crate::cache_source::IntegralKey::Descendant { g, d, ch } = row.key
                        else {
                            unreachable!()
                        };
                        assert_eq!(row.value, Descendents::default().eval(g, d.clone(), &ch));
                        assert!(seen.insert((g, d, ch), row.value).is_none());
                    }
                    cursor = batch.next;
                }
            }
            assert_eq!(seen.len(), source.ready() - source.floor());
            storage.close().unwrap();
        }
    }

    #[test]
    fn unavailable_spill_is_resource_error_and_normal_retry_remains_exact() {
        let registry = Arc::new(Registry::new());
        let mut tables = Descendents::with_registry(Some(&registry));
        let poly = [(vec![(1, 6)], Q::ONE)];
        let expected = Descendents::default().eval(3, vec![], &poly[0].0);
        let error = crate::cache_source::catch_resource(|| {
            let (value, _) = pressure::force_backing(|| {
                integrate(
                    3,
                    &[(vec![], Q::ONE)],
                    &poly,
                    &mut tables,
                    2,
                    Some(&registry),
                )
            });
            Ok(value)
        })
        .unwrap_err();
        assert!(error.contains("plan spill unavailable"), "{error}");
        assert!(
            registry.check_error().is_ok(),
            "resource failure must permit retry"
        );
        assert_eq!(
            integrate(
                3,
                &[(vec![], Q::ONE)],
                &poly,
                &mut tables,
                2,
                Some(&registry),
            ),
            expected
        );
    }

    #[test]
    fn shared_parallel_results_survive_queries_and_match_the_serial_oracle() {
        use crate::cache_source::ResultSource;
        let registry = Arc::new(Registry::new());
        let ch = vec![(1, 6)];
        let expected = Descendents::default().eval(3, Vec::new(), &ch);
        let mut tables = Descendents::with_registry(Some(&registry));
        assert_eq!(
            integrate(
                3,
                &[(vec![], Q::ONE)],
                &[(ch.clone(), Q::ONE)],
                &mut tables,
                4,
                Some(&registry)
            ),
            expected
        );
        let source = with_engine(Some(&registry), |engine| {
            assert!(
                engine.values.is_empty(),
                "one mathematical result table only"
            );
            engine.source.as_ref().unwrap().clone()
        });
        let ready = source.ready();
        let slots = source.read().len();
        assert!(ready > 0);
        let mut cursor = 0;
        let mut seen = BTreeMap::new();
        while cursor < ready {
            if let Some(batch) = source.pull(cursor, 256) {
                assert!(batch.next > cursor);
                for row in batch.rows {
                    let crate::cache_source::IntegralKey::Descendant { g, d, ch } = row.key else {
                        unreachable!()
                    };
                    assert_eq!(row.value, Descendents::default().eval(g, d.clone(), &ch));
                    assert!(seen.insert((g, d, ch), row.value).is_none());
                }
                cursor = batch.next;
            }
        }
        assert_eq!(seen.len(), ready);
        assert_eq!(
            integrate(
                3,
                &[(vec![], Q::ONE)],
                &[(ch, Q::ONE)],
                &mut tables,
                4,
                Some(&registry)
            ),
            expected
        );
        assert_eq!(source.ready(), ready);
        assert_eq!(source.read().len(), slots);
        let pure = vec![2; 6];
        let value = cached(3, &pure, &[], &registry).expect("pure descendant of C1^6");
        assert_eq!(value, Descendents::default().eval(3, pure.clone(), &[]));
        assert_eq!(
            integrate(
                3,
                &[(pure, Q::ONE)],
                &[(vec![], Q::ONE)],
                &mut tables,
                1,
                Some(&registry)
            ),
            value
        );
        assert_eq!(source.ready(), ready);
    }

    #[test]
    fn registry_keeps_parallel_values_after_calling_thread_exits() {
        use crate::cache_source::ResultSource;
        let registry = Arc::new(Registry::new());
        let worker_registry = Arc::clone(&registry);
        let source = std::thread::spawn(move || {
            let mut tables = Descendents::with_registry(Some(&worker_registry));
            integrate(
                3,
                &[(vec![], Q::ONE)],
                &[(vec![(1, 6)], Q::ONE)],
                &mut tables,
                4,
                Some(&worker_registry),
            );
            with_engine(Some(&worker_registry), |engine| {
                engine.source.as_ref().unwrap().clone()
            })
        })
        .join()
        .unwrap();
        let ready = source.ready();
        assert!(ready > 0);
        let mut cursor = 0;
        while cursor < ready {
            if let Some(batch) = source.pull(cursor, 128) {
                cursor = batch.next;
            }
        }
        assert_eq!(cursor, ready);
    }

    #[test]
    fn background_reader_covers_real_parallel_and_formula_sources() {
        use crate::cache_source::ResultSource;
        let registry = Arc::new(Registry::new());
        let reader = registry.start_reader().unwrap();
        let mut tables = Descendents::with_registry(Some(&registry));
        let ch = vec![(1, 6)];
        let expected = Descendents::default().eval(3, vec![], &ch);
        assert_eq!(
            integrate(
                3,
                &[(vec![], Q::ONE)],
                &[(ch, Q::ONE)],
                &mut tables,
                4,
                Some(&registry)
            ),
            expected
        );
        let total = registry
            .sources()
            .iter()
            .map(|source| source.ready())
            .sum::<usize>();
        let report = reader.finish();
        assert!(report.reached_target);
        assert!(!report.worker_panicked);
        assert!(!report.invalid_cursor);
        assert_eq!(report.rows, total);
        assert!(total > 0);
        // The read service observes copies; it never consumes the mathematical
        // table or makes the completed descendants unavailable to later calls.
        let source = with_engine(Some(&registry), |engine| {
            engine.source.as_ref().unwrap().clone()
        });
        assert!(source.ready() > 0);
        assert!(cached(3, &[2; 6], &[], &registry).is_some());
    }

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
