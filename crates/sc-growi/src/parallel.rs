//! 分裂公式中的独立数值项。变量项及增量消元仍由原来的调用线程处理。
//! 每个工作线程有独立递归上下文，共享已完成结果，保留乘积的零项短路顺序。
use crate::num::{One, Zero};
use crate::{Ctx, Inv, Q};

// 小批次直接递归；限制待算项数，避免分裂枚举一次性展开占用大量内存。
const PARALLEL_THRESHOLD: usize = 32;
const MAX_PENDING: usize = 512;
const CHUNK_SIZE: usize = 4;

struct Product {
    coefficient: Q,
    factors: Vec<Inv>,
    variable: Option<i32>,
}

enum Contribution {
    Number(Q),
    Variable(i32, Q),
}

impl Product {
    fn evaluate(&self, ctx: &mut Ctx) -> Option<Contribution> {
        let mut value = Q::ONE;
        for (i, inv) in self.factors.iter().enumerate() {
            if let Some(name) = self.variable.filter(|_| i + 1 == self.factors.len()) {
                let mut probe = inv.clone();
                let Some(known) = crate::simplify::simplify(ctx, &mut probe) else {
                    return Some(Contribution::Variable(name, &self.coefficient * value));
                };
                value *= known;
            } else {
                value *= ctx.compute(inv)?;
            }
            if value == Q::ZERO {
                break;
            }
        }
        Some(Contribution::Number(&self.coefficient * value))
    }
}

pub(crate) struct Sum {
    value: Q,
    pending: Vec<Product>,
    threads: usize,
    variables: Vec<(i32, Q)>,
}

impl Sum {
    pub(crate) fn new(ctx: &Ctx) -> Self {
        Self {
            value: Q::ZERO,
            pending: Vec::new(),
            threads: ctx.parallelism(),
            variables: Vec::new(),
        }
    }

    pub(crate) fn pair(
        &mut self,
        ctx: &mut Ctx,
        coefficient: &Q,
        first: Inv,
        second: Inv,
    ) -> Option<()> {
        if self.threads == 1 {
            let left = ctx.compute(&first)?;
            if left != Q::ZERO {
                let right = ctx.compute(&second)?;
                if right != Q::ZERO {
                    self.value += coefficient * left * right;
                }
            }
            return Some(());
        }
        self.push(
            ctx,
            Product {
                coefficient: coefficient.clone(),
                factors: vec![first, second],
                variable: None,
            },
        )
    }

    pub(crate) fn single(&mut self, ctx: &mut Ctx, coefficient: &Q, inv: &Inv) -> Option<()> {
        if self.threads == 1 {
            self.value += coefficient * ctx.compute(inv)?;
            return Some(());
        }
        self.push(
            ctx,
            Product {
                coefficient: coefficient.clone(),
                factors: vec![inv.clone()],
                variable: None,
            },
        )
    }

    pub(crate) fn product(
        &mut self,
        ctx: &mut Ctx,
        coefficient: Q,
        factors: impl Iterator<Item = Inv>,
    ) -> Option<()> {
        if self.threads == 1 {
            let mut value = Q::ONE;
            for inv in factors {
                value *= ctx.compute(&inv)?;
                if value == Q::ZERO {
                    break;
                }
            }
            self.value += coefficient * value;
            return Some(());
        }
        self.push(
            ctx,
            Product {
                coefficient,
                factors: factors.collect(),
                variable: None,
            },
        )
    }

    /// TRR 的两点分量：末因子只用 simplify，尚未知时保留线性变量。
    /// 变量名由未排序的骨架决定；工作线程只返回系数，不接触方程组。
    pub(crate) fn conditional_pair(
        &mut self,
        ctx: &mut Ctx,
        coefficient: &Q,
        first: Inv,
        mut second: Inv,
        name: i32,
    ) -> Option<()> {
        if self.threads == 1 {
            let left = ctx.compute(&first)?;
            if left != Q::ZERO {
                match crate::simplify::simplify(ctx, &mut second) {
                    Some(right) => self.value += coefficient * left * right,
                    None => self.variables.push((name, coefficient * left)),
                }
            }
            return Some(());
        }
        self.push(
            ctx,
            Product {
                coefficient: coefficient.clone(),
                factors: vec![first, second],
                variable: Some(name),
            },
        )
    }

    fn add(&mut self, value: Contribution) {
        match value {
            Contribution::Number(value) => self.value += value,
            Contribution::Variable(name, value) => self.variables.push((name, value)),
        }
    }

    fn push(&mut self, ctx: &mut Ctx, mut product: Product) -> Option<()> {
        let mut known = 0;
        for inv in &mut product.factors {
            let point = if inv.cdeg == 0 && inv.genus >= 0 {
                crate::simplify::mapping_to_point(inv)
            } else {
                None
            };
            let value = point.or_else(|| {
                if inv.deg == 1 {
                    inv.canonicalize();
                    ctx.cache.get(inv)
                } else {
                    // 原始类归约会重命名；查询未命中时保留原输入，交给
                    // simplify 按原顺序排序、归约一次。
                    let mut key = inv.clone();
                    key.canonicalize();
                    key.reduce();
                    ctx.cache.get(&key)
                }
            });
            let Some(value) = value else {
                break;
            };
            if value == Q::ZERO {
                return Some(());
            }
            product.coefficient *= value;
            known += 1;
        }
        if known == product.factors.len() {
            self.value += product.coefficient;
            return Some(());
        }
        product.factors.drain(..known);
        self.pending.push(product);
        if self.pending.len() >= MAX_PENDING {
            self.flush(ctx)?;
        }
        Some(())
    }

    fn flush(&mut self, ctx: &mut Ctx) -> Option<()> {
        if self.pending.len() < PARALLEL_THRESHOLD || self.threads == 1 {
            for product in std::mem::take(&mut self.pending) {
                self.add(product.evaluate(ctx)?);
            }
            return Some(());
        }
        let shared = ctx.cache.start_parallel();
        #[cfg(test)]
        let start = std::time::Instant::now();
        let depth = ctx.depth();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            sc_num::parallel::chunks(
                self.pending.len(),
                CHUNK_SIZE,
                self.threads,
                || Ctx::worker(crate::cache::Cache::worker(shared.clone()), depth),
                |worker, range| {
                    range
                        .map(|i| self.pending[i].evaluate(worker))
                        .collect::<Option<Vec<_>>>()
                },
            )
        }));
        // 工作线程退出后才归并和写回 SQLite，包括求和过程中出现 None 的情况。
        ctx.cache.finish_parallel(shared);
        #[cfg(test)]
        {
            ctx.parallel_profile.0 += 1;
            ctx.parallel_profile.1 += self.pending.len();
            ctx.parallel_profile.2 += start.elapsed();
        }
        let results = match result {
            Ok(results) => results,
            Err(payload) => std::panic::resume_unwind(payload),
        };
        for (chunk, values) in results.into_iter().enumerate() {
            // 可解性可能依赖同批次其它项刚算出的值。归并后顺序重试，
            // 不把工作线程的临时 None 当作最终失败。
            let values = match values {
                Some(values) => values,
                None => self.pending
                    [chunk * CHUNK_SIZE..((chunk + 1) * CHUNK_SIZE).min(self.pending.len())]
                    .iter()
                    .map(|product| product.evaluate(ctx))
                    .collect::<Option<Vec<_>>>()?,
            };
            for value in values {
                self.add(value);
            }
        }
        self.pending.clear();
        Some(())
    }

    pub(crate) fn finish(mut self, ctx: &mut Ctx) -> Option<Q> {
        self.flush(ctx)?;
        debug_assert!(self.variables.is_empty());
        Some(self.value)
    }

    pub(crate) fn finish_linear(
        mut self,
        ctx: &mut Ctx,
        sys: &mut crate::linsys::Linsys,
        coefficient: &Q,
    ) -> Option<Q> {
        self.flush(ctx)?;
        for (name, value) in self.variables {
            if !sys.add_var_term(name, -(coefficient * value)) {
                return None;
            }
        }
        Some(coefficient * self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Pdata;

    #[test]
    fn single_queries_match_serial_and_reuse_results() {
        let mut batches = 0;
        for (text, expected) in [
            ("G=3,D=7 in P^2 thru H^2:23", "153796445095"),
            ("G=2,D=3 in P^3 thru H^2:12", "5930"),
            ("G=2,D=5 in P^4 thru H^3*psi^22", "-41369/110075314176"),
            ("G=0,D=2,E=5,N=1 thru [1,2]:5", "238200"),
            ("G=0,D=2,E=5,N=1 thru [H*psi,2],[H,2]", "25/2"),
            ("G=0,D=1,E=3,N=2 thru pc(1),pc(1)", "-6"),
        ] {
            let query = crate::query::parse(text).unwrap();
            for threads in [1, 4] {
                let mut ctx = Ctx::memory();
                ctx.threads = threads;
                let value = crate::gw(&mut ctx, &query.settings, &query.classes).unwrap();
                assert_eq!(value.to_string(), expected, "{text}, threads={threads}");
                let count = ctx.cache.len();
                assert_eq!(
                    crate::gw(&mut ctx, &query.settings, &query.classes).unwrap(),
                    value
                );
                assert_eq!(ctx.cache.len(), count);
                batches += ctx.parallel_profile.0;
            }
        }
        assert!(batches > 0);
    }

    #[test]
    fn conditional_coefficients_keep_variable_order_and_zero_shortcuts() {
        let one = Inv::new(0, 1, 0, 0, [Pdata::default(); 3]);
        let zero = Inv::new(0, 1, 0, 0, [Pdata::new(1, 0, 0); 3]);
        let unsupported = Inv::new(2, 3, 2, 1, [Pdata::new(2, 1, 0)]);
        let mut rows = Vec::new();
        for threads in [1, 4] {
            let mut ctx = Ctx::memory();
            ctx.threads = threads;
            let mut sum = Sum::new(&ctx);
            for i in 0..40 {
                let h = 1 + i % 2;
                let inv = Inv::new(2, 1, 3, 1, [Pdata::new(h, 6 - h, 0)]);
                let coefficient = Q::from(if i % 3 == 1 { -i - 1 } else { i + 1 });
                sum.conditional_pair(&mut ctx, &coefficient, one.clone(), inv, 3 + i32::from(h))
                    .unwrap();
            }
            sum.pair(&mut ctx, &Q::ONE, zero.clone(), unsupported.clone())
                .unwrap();
            sum.single(&mut ctx, &Q::from_signeds(2, 5), &one).unwrap();
            let mut sys = crate::linsys::Linsys::new(2);
            let constant = sum
                .finish_linear(&mut ctx, &mut sys, &Q::from_signeds(3, 7))
                .unwrap();
            assert_eq!(constant, Q::from_signeds(6, 35));
            assert_eq!((sys.numvars(), sys.varname(0), sys.varname(1)), (2, 4, 5));
            sys.add_const_term(constant);
            rows.push(sys.debug_current_row());
            if threads > 1 {
                assert!(ctx.parallel_profile.0 > 0);
            }
        }
        assert_eq!(rows[0], rows[1]);
    }

    #[test]
    fn unresolved_batch_restores_the_callers_cache() {
        let mut ctx = Ctx::memory();
        ctx.threads = 4;
        let known = Inv::new(2, 1, 0, 1, [Pdata::new(2, 0, 0); 2]);
        ctx.cache.put(known.clone(), Q::ONE);
        let unsupported = Inv::new(2, 3, 2, 1, [Pdata::new(2, 1, 0)]);
        let mut sum = Sum::new(&ctx);
        for _ in 0..32 {
            sum.single(&mut ctx, &Q::ONE, &unsupported).unwrap();
        }
        assert_eq!(sum.finish(&mut ctx), None);
        assert_eq!(ctx.cache.get(&known), Some(Q::ONE));
        assert_eq!(ctx.compute(&known), Some(Q::ONE));
    }

    #[test]
    #[ignore = "手动比较 GROWI 单查询的并行度与批次耗时"]
    fn benchmark_single_query() {
        let query = crate::query::parse("G=3,D=10 in P^2 thru H^2:32").unwrap();
        for threads in [1, 4, 16, 32] {
            let mut ctx = Ctx::memory();
            ctx.threads = threads;
            let start = std::time::Instant::now();
            let result = crate::gw(&mut ctx, &query.settings, &query.classes).unwrap();
            assert_eq!(result.to_string(), "38747510483053595091600");
            eprintln!(
                "threads={threads}, elapsed={:?}, batches={}, terms={}, parallel_time={:?}, states={}",
                start.elapsed(),
                ctx.parallel_profile.0,
                ctx.parallel_profile.1,
                ctx.parallel_profile.2,
                ctx.cache.len()
            );
        }
    }
}
