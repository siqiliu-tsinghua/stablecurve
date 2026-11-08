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

use crate::cache_source::storage::{Storage, StorageHandle};
use crate::cache_source::{Reader, Registry};
use sc_num::{One, Q, Zero, is_zero};
use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Arc;

thread_local! {
    static MEMO: RefCell<HashMap<Vec<u32>, Q>> = RefCell::new(HashMap::new());
}

/// STrecur 结果的 SQLite L2（跨会话持久；STrecur 子问题复用率高，
/// 是 Faber TAU 表的现代等价物）。
static L2: std::sync::OnceLock<std::sync::Mutex<Option<Arc<CacheSession>>>> =
    std::sync::OnceLock::new();

struct Service {
    reader: Option<Reader>,
    storage: Option<Storage>,
    _sampler: crate::cache_source::sampler::SamplerLease,
}

/// Set the optional total-process memory ceiling shared by all FABER cache
/// sessions. `None` derives the budget from currently available system memory.
/// The ceiling includes computation and storage progress space, not just L1.
pub fn set_cache_memory_limit(maximum_bytes: Option<u64>) -> Result<(), String> {
    let config = crate::cache_source::budget::BudgetConfig {
        max_bytes: maximum_bytes,
        ..Default::default()
    };
    crate::cache_source::budget::process_budget().configure(config)
}

enum ServiceState {
    Open(Service),
    Closing,
    Closed(Result<(), String>),
}

/// An independent cache owner. Mathematical calls run concurrently; explicit
/// close rejects new calls, waits for existing calls, then confirms durability.
pub struct CacheSession {
    registry: Arc<Registry>,
    handle: Option<StorageHandle>,
    service: std::sync::Mutex<ServiceState>,
    finished: std::sync::Condvar,
}

impl CacheSession {
    pub fn open(path: &str) -> Result<Self, String> {
        let sampler = crate::cache_source::sampler::start()?;
        let registry = Arc::new(Registry::new());
        let pull = std::env::var("SC_FABER_CACHE_MEASURE_PULL").as_deref() == Ok("1");
        let source_only = std::env::var("SC_FABER_CACHE_MEASURE_SOURCE_ONLY").as_deref() == Ok("1");
        let reader = if pull {
            Some(
                registry
                    .start_reader()
                    .map_err(|error| format!("启动 C1 读取线程失败：{error}"))?,
            )
        } else {
            None
        };
        let storage = if pull || source_only {
            None
        } else {
            let storage = Storage::open(path, Arc::clone(&registry))?;
            registry.bind_storage(storage.handle());
            Some(storage)
        };
        let handle = storage.as_ref().map(Storage::handle);
        Ok(Self {
            registry,
            handle,
            service: std::sync::Mutex::new(ServiceState::Open(Service {
                reader,
                storage,
                _sampler: sampler,
            })),
            finished: std::sync::Condvar::new(),
        })
    }

    /// Bind this owner for one query, including nested mathematical calls.
    /// The session mutex is never held while the query runs.
    pub fn run<T>(&self, calculate: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
        self.registry.check_error()?;
        let result = crate::cache_source::try_with_registry(Arc::clone(&self.registry), || {
            crate::cache_source::catch_resource(calculate)
        })?;
        if result.is_ok() {
            self.registry.check_error()?;
        }
        result
    }

    /// Run explicitly without any persistent owner, including the global one.
    pub fn run_uncached<T>(calculate: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
        crate::cache_source::with_optional_registry(None, calculate)
    }

    pub fn error(&self) -> Option<String> {
        self.registry.check_error().err()
    }

    pub fn is_open(&self) -> bool {
        !self.registry.is_closing() && self.error().is_none()
    }

    /// Commit a finite snapshot of published results. Queries may continue.
    pub fn flush(&self) -> Result<(), String> {
        self.registry.check_error()?;
        if self.registry.is_closing() && !crate::cache_source::has_active_scope(self.registry.id())
        {
            return Err("cache session is closing or closed".to_owned());
        }
        if let Some(storage) = &self.handle {
            let report = storage.flush()?;
            if std::env::var("SC_FABER_CACHE_MEASURE_REPORT").as_deref() == Ok("1") {
                eprintln!("sc-faber: C2 storage flush {report:?}");
            }
        }
        Ok(())
    }

    pub fn close(&self) -> Result<(), String> {
        if crate::cache_source::has_active_scope(self.registry.id()) {
            return Err("cannot close a cache session inside its active query".to_owned());
        }
        self.registry.begin_close();
        let service = {
            let mut state = self
                .service
                .lock()
                .unwrap_or_else(|error| error.into_inner());
            loop {
                match &*state {
                    ServiceState::Open(_) => {
                        let ServiceState::Open(service) =
                            std::mem::replace(&mut *state, ServiceState::Closing)
                        else {
                            unreachable!()
                        };
                        break service;
                    }
                    ServiceState::Closing => {
                        state = self
                            .finished
                            .wait(state)
                            .unwrap_or_else(|error| error.into_inner());
                    }
                    ServiceState::Closed(result) => return result.clone(),
                }
            }
        };
        // Neither the global attachment nor service mutex is held while old
        // queries flush/read the backend or finish their mathematical work.
        self.registry.wait_idle();
        let result =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| finish_service(service)))
                .unwrap_or_else(|_| Err("cache session close panicked".to_owned()));
        let result = result.and_then(|()| self.registry.check_error());
        if let Err(error) = &result {
            self.registry.fail(error.clone());
        }
        if result.is_ok() {
            // Every query has left and the storage worker has committed and
            // joined. Clear the primary tables themselves, including sources
            // retained in another live thread's TLS context.
            for source in self.registry.sources() {
                source.release();
            }
        }
        self.registry.release_sources();
        crate::hodge::clear_cache_context(self.registry.id());
        ST_SOURCES.with(|sources| {
            if let Ok(mut sources) = sources.try_borrow_mut() {
                sources.remove(&self.registry.id());
            }
        });
        let mut state = self
            .service
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        *state = ServiceState::Closed(result.clone());
        self.finished.notify_all();
        result
    }
}

impl Drop for CacheSession {
    fn drop(&mut self) {
        if matches!(
            self.service
                .get_mut()
                .unwrap_or_else(|error| error.into_inner()),
            ServiceState::Closed(_)
        ) {
            return;
        }
        if let Err(error) = self.close() {
            eprintln!("sc-faber: 关闭持久缓存失败：{error}");
        }
    }
}

fn l2() -> &'static std::sync::Mutex<Option<Arc<CacheSession>>> {
    L2.get_or_init(|| std::sync::Mutex::new(None))
}

/// Legacy process-global owner. Resolve it once at a query entry; the captured
/// Arc remains valid if a later attachment replaces or detaches this one.
pub(crate) fn result_registry() -> Option<Arc<Registry>> {
    l2().lock()
        .unwrap()
        .as_ref()
        .map(|session| Arc::clone(&session.registry))
}

/// 卸下持久缓存（之后只用进程内记忆化）。
pub fn detach_cache() {
    if let Err(error) = close_cache() {
        eprintln!("sc-faber: 关闭持久缓存失败：{error}");
    }
}

/// 等待已完成成果提交，然后卸下进程全局兼容缓存。
/// 需要获取异步存储错误时使用本方法，而非无法返回错误的析构路径。
pub fn close_cache() -> Result<(), String> {
    let attached = l2().lock().unwrap().take();
    if let Some(attached) = attached {
        attached.close()?;
    }
    Ok(())
}

/// Diagnostic teardown only. The global attachment lock is never held while
/// a reader captures source watermarks, drains copied rows, or joins.
fn finish_service(mut service: Service) -> Result<(), String> {
    if let Some(reader) = service.reader.take() {
        let started = std::time::Instant::now();
        let report = reader.finish();
        eprintln!(
            "sc-faber: C1 pull reader {report:?}, closing_seconds={:.6}",
            started.elapsed().as_secs_f64()
        );
        if !report.reached_target || report.worker_panicked || report.invalid_cursor {
            return Err("C1 measurement did not finish copying its published target".to_owned());
        }
    }
    if let Some(storage) = service.storage.take() {
        let report = storage.close()?;
        if std::env::var("SC_FABER_CACHE_MEASURE_REPORT").as_deref() == Ok("1") {
            eprintln!("sc-faber: C2 storage close {report:?}");
        }
    }
    Ok(())
}

/// 当前是否挂接了持久缓存（进程全局）。
pub fn cache_attached() -> bool {
    l2().lock()
        .unwrap()
        .as_ref()
        .is_some_and(|session| session.is_open())
}

/// 最近已知的异步缓存错误；数学结果不等同于持久化成功。
pub fn cache_error() -> Option<String> {
    l2().lock()
        .unwrap()
        .as_ref()
        .and_then(|session| session.error())
}

/// 等待调用前已完成的积分落盘；普通积分调用不隐式等待。
pub fn flush_cache() -> Result<(), String> {
    // The handle has a fixed owner and no global lock is held during the wait.
    let session = {
        let guard = l2().lock().unwrap();
        let Some(attached) = guard.as_ref() else {
            return Ok(());
        };
        Arc::clone(attached)
    };
    session.flush()
}

/// 挂接持久缓存文件（不存在则创建）。L2 是进程全局的，重复调用替换连接。
pub fn attach_cache(path: &str) -> Result<(), String> {
    let session = Arc::new(CacheSession::open(path)?);
    let previous = l2().lock().unwrap().replace(session);
    if let Some(previous) = previous {
        previous.close()?;
    }
    Ok(())
}

/// 通用 L2 读取（跨模块共用同一张表，按 kind 前缀区分）。
pub fn cache_get(kind: &str, nums: &[u32]) -> Option<Q> {
    crate::cache_source::with_current_registry(|registry| {
        let registry = registry?;
        if let Some(value) = registry.manual_get(kind, nums) {
            return Some(value);
        }
        if !registry.has_history() {
            return None;
        }
        registry
            .lookup(crate::cache_source::IntegralKey::Legacy {
                kind: kind.into(),
                nums: nums.to_vec(),
            })
            .ok()
            .flatten()
    })
}

/// 通用 L2 写入。
pub fn cache_put(kind: &str, nums: &[u32], val: &Q) {
    crate::cache_source::with_current_registry(|registry| {
        if let Some(registry) = registry {
            registry.manual_put(kind, nums, val);
        }
    });
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
/// 生成函数递推：Σ_{k=0}^{m} C(m+1,k) B_k = 0（m ≥ 1）。表按需延长并缓存：
/// Bernoulli 多项式等调用方逐个取 B_0..B_n，每次从头递推会叠成 O(n⁴)。
pub fn bernoulli(n: u64) -> Q {
    sc_num::bernoulli(n)
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

type StSource = crate::cache_source::legacy_source::LegacySource<Vec<u32>>;

thread_local! {
    static ST_SOURCES: RefCell<std::collections::BTreeMap<u64, std::sync::Weak<StSource>>> =
        const { RefCell::new(std::collections::BTreeMap::new()) };
}

fn st_nums<K: AsRef<[u32]>>(key: &K) -> Vec<u32> {
    key.as_ref().to_vec()
}

fn st_source(registry: &Registry) -> Arc<StSource> {
    ST_SOURCES.with(|cell| {
        let mut sources = cell.borrow_mut();
        if let Some(source) = sources
            .get(&registry.id())
            .and_then(std::sync::Weak::upgrade)
        {
            return source;
        }
        sources.retain(|_, source| source.strong_count() != 0);
        let source = Arc::new(StSource::new("st", st_nums));
        registry.register(source.clone());
        sources.insert(registry.id(), Arc::downgrade(&source));
        source
    })
}

/// ∫纯 ψ 顶度递归（admcycles `STrecur`）。`psi` 内部会排序。
pub fn st_recur(psi: &[u32]) -> Q {
    crate::cache_source::with_current_registry(|registry| {
        let Some(registry) = registry else {
            let mut p = psi.to_vec();
            p.sort_unstable();
            if let Some(value) = MEMO.with(|memo| memo.borrow().get(&p).cloned()) {
                return value;
            }
            let value = st_recur_calc(&p);
            MEMO.with(|memo| memo.borrow_mut().insert(p, value.clone()));
            return value;
        };
        let source = st_source(registry);
        let sorted = psi.is_sorted();
        if sorted {
            if let Some(value) = source.get_slice(psi) {
                return value;
            }
        }
        let mut p = psi.to_vec();
        if !sorted {
            p.sort_unstable();
            if let Some(value) = source.get_slice(&p) {
                return value;
            }
        }
        // A storage error cannot become a false mathematical zero. Once the
        // owner has failed, skip further L2 attempts and retain exact results.
        if registry.has_history() && registry.check_error().is_ok() {
            match registry.lookup(crate::cache_source::IntegralKey::Legacy {
                kind: "st".to_owned(),
                nums: p.clone(),
            }) {
                Ok(Some(value)) => {
                    source.insert(p, &value);
                    return value;
                }
                Ok(None) => {}
                Err(error) => registry.fail(error),
            }
        }
        let value = st_recur_calc(&p);
        source.insert(p, &value);
        value
    })
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
    // dilaton：⟨τ₁X⟩_{g,n} = (2g−2+n−1)⟨X⟩_{g,n−1}（升序表的表头；要求剩余
    // 空间稳定，⟨τ₁⟩₁ 等基例仍走下面的一般递推）。省掉一整轮分裂求和。
    if psi[0] == 1 && 2 * g - 2 + (n - 1) > 0 {
        return Q::from(2 * g - 2 + n - 1) * st_recur(&psi[1..]);
    }
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
    crate::cache_source::with_current_registry(|registry| {
        let Some(registry) = registry else {
            return st_sum_calc(psilist, kappalist);
        };
        let degree = psilist
            .iter()
            .chain(kappalist)
            .map(|&power| i64::from(power))
            .sum::<i64>();
        let numerator = degree - psilist.len() as i64 + 3;
        let genus = if numerator >= 0 && numerator % 3 == 0 {
            u32::try_from(numerator / 3).ok()
        } else {
            None
        };
        let Some(g) = genus else {
            return st_sum_calc(psilist, kappalist);
        };
        match with_socle_root(registry, g, psilist, kappalist, || {
            Ok(st_sum_calc(psilist, kappalist))
        }) {
            Ok(value) => value,
            Err(error) => {
                registry.fail(error);
                st_sum_calc(psilist, kappalist)
            }
        }
    })
}

fn st_sum_calc(psilist: &[u32], kappalist: &[u32]) -> Q {
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
    let cache_enabled = crate::cache_source::with_current_registry(|registry| registry.is_some());
    let mut kappa = kappalist.to_vec();
    kappa.sort_unstable();
    let mut total = Q::ZERO;
    for (blocks, aut) in setparts_with_auts(&kappa) {
        let mut arg: Vec<u32> = blocks.iter().map(|b| 1 + b.iter().sum::<u32>()).collect();
        arg.extend_from_slice(psilist);
        let sign = if blocks.len() % 2 == 0 { 1 } else { -1 };
        let value = if cache_enabled {
            crate::hodge::psi_integral(0, &arg)
        } else {
            multinomial(&arg)
        };
        total += Q::from(sign) * aut * value;
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
    crate::cache_source::with_current_registry(|registry| {
        if let Some(registry) = registry {
            with_socle_root(registry, g, psilist, kappalist, || {
                socle_formula_calc(g, psilist, kappalist)
            })
        } else {
            socle_formula_calc(g, psilist, kappalist)
        }
    })
}

fn with_socle_root(
    registry: &Arc<Registry>,
    g: u32,
    psilist: &[u32],
    kappalist: &[u32],
    calculate: impl FnOnce() -> Result<Q, String>,
) -> Result<Q, String> {
    let Ok(n) = u32::try_from(psilist.len()) else {
        return calculate();
    };
    if 2 * i64::from(g) - 2 + i64::from(n) <= 0 {
        return calculate();
    }
    let psi: Vec<_> = psilist
        .iter()
        .enumerate()
        .filter_map(|(leg, &power)| (power != 0).then_some((leg as u32 + 1, power)))
        .collect();
    let mut indices = kappalist.to_vec();
    indices.sort_unstable();
    let mut kappa: Vec<(u32, u32)> = Vec::new();
    for index in indices {
        if let Some(last) = kappa.last_mut().filter(|last| last.0 == index) {
            last.1 += 1;
        } else {
            kappa.push((index, 1));
        }
    }
    crate::hodge::root_source::with_registry(
        registry.clone(),
        g,
        n,
        &kappa,
        &psi,
        &[],
        &[],
        calculate,
    )
}

fn socle_formula_calc(g: u32, psilist: &[u32], kappalist: &[u32]) -> Result<Q, String> {
    if g == 0 {
        return Ok(ct_const(0) * ct_sum(psilist, kappalist));
    }
    if kappalist.is_empty() {
        return Ok(crate::hodge::psi_integral(g, psilist));
    }
    let mut roots = std::collections::BTreeMap::<Vec<u32>, Q>::new();
    for (blocks, multiplicity) in setparts_with_auts(kappalist) {
        let mut arg = psilist.to_vec();
        arg.extend(blocks.iter().map(|b| 1 + b.iter().sum::<u32>()));
        arg.sort_unstable();
        let coefficient = roots.entry(arg).or_insert(Q::ZERO);
        if (kappalist.len() - blocks.len()) % 2 == 0 {
            *coefficient += multiplicity;
        } else {
            *coefficient -= multiplicity;
        }
    }
    Ok(crate::hodge::psi_sum(
        g,
        &roots
            .into_iter()
            .filter(|(_, c)| *c != Q::ZERO)
            .collect::<Vec<_>>(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn successful_close_reclaims_sources_held_in_another_live_threads_context() {
        let session = Arc::new(CacheSession::open(":memory:").unwrap());
        let (completed, observed) = std::sync::mpsc::channel();
        let (release, released) = std::sync::mpsc::channel();
        let owner = session.clone();
        let worker = std::thread::spawn(move || {
            let value = owner
                .run(|| crate::hodge::klp(7, 0, &[], &[], &[(1, 18)], &[]))
                .unwrap();
            assert_eq!(value, Q::from_signeds(32017001, 638512875));
            assert!(completed.send(owner.registry.sources()).is_ok());
            // Keep this thread and its original Engine/Descendents TLS alive
            // while another thread closes the shared cache session.
            released.recv().unwrap();
        });
        let sources = observed.recv().unwrap();
        assert!(sources.iter().any(|source| source.ready() > source.floor()));
        session.close().unwrap();
        for source in &sources {
            assert_eq!(source.floor(), source.ready());
        }
        release.send(()).unwrap();
        worker.join().unwrap();
    }

    #[test]
    fn independent_sessions_and_explicit_no_cache_keep_their_own_values() {
        let first = CacheSession::open(":memory:").unwrap();
        let second = CacheSession::open(":memory:").unwrap();
        first
            .run(|| {
                cache_put("owner-test", &[1], &Q::from(11));
                second.run(|| {
                    assert_eq!(cache_get("owner-test", &[1]), None);
                    cache_put("owner-test", &[1], &Q::from(22));
                    Ok(())
                })?;
                CacheSession::run_uncached(|| {
                    assert_eq!(cache_get("owner-test", &[1]), None);
                    Ok(())
                })?;
                assert_eq!(cache_get("owner-test", &[1]), Some(Q::from(11)));
                Ok(())
            })
            .unwrap();
        first.close().unwrap();
        second
            .run(|| {
                assert_eq!(cache_get("owner-test", &[1]), Some(Q::from(22)));
                Ok(())
            })
            .unwrap();
        second.close().unwrap();
    }

    #[test]
    fn close_waits_for_concurrent_queries_and_their_final_publications() {
        use std::sync::mpsc;

        let session = Arc::new(CacheSession::open(":memory:").unwrap());
        let (entered, observed) = mpsc::channel();
        let mut workers = Vec::new();
        let mut releases = Vec::new();
        for number in [1, 2] {
            let session = Arc::clone(&session);
            let entered = entered.clone();
            let (release, wait) = mpsc::channel();
            releases.push(release);
            workers.push(std::thread::spawn(move || {
                session.run(|| {
                    entered.send(()).unwrap();
                    wait.recv().unwrap();
                    cache_put("closing-test", &[number], &Q::from(number));
                    // The closer must not own a service/global mutex while it
                    // waits for this already entered query to finish and flush.
                    session.flush()?;
                    Ok(number)
                })
            }));
        }
        observed
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap();
        observed
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap();
        session.registry.begin_close();
        let (closed, closing) = mpsc::channel();
        let owner = Arc::clone(&session);
        let closer = std::thread::spawn(move || {
            closed.send(owner.close()).unwrap();
        });
        assert!(closing.try_recv().is_err());
        assert!(session.run(|| Ok(())).is_err());
        for release in releases {
            release.send(()).unwrap();
        }
        for worker in workers {
            worker.join().unwrap().unwrap();
        }
        closing
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap()
            .unwrap();
        closer.join().unwrap();
        assert!(!session.is_open());
        assert!(session.registry.sources().is_empty());
        assert!(session.registry.manual_get("closing-test", &[1]).is_none());
        session.close().unwrap();
    }

    #[test]
    fn query_unwind_releases_its_permit_and_self_close_is_rejected() {
        let session = CacheSession::open(":memory:").unwrap();
        session
            .run(|| {
                assert!(session.close().is_err());
                Ok(())
            })
            .unwrap();
        assert!(session.is_open());
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = session.run(|| -> Result<(), String> { panic!("test query unwind") });
        }));
        assert!(panic.is_err());
        session.close().unwrap();
    }

    #[test]
    fn masked_outer_owner_cannot_close_itself_and_keeps_its_original_permit() {
        let first = CacheSession::open(":memory:").unwrap();
        let second = CacheSession::open(":memory:").unwrap();
        first
            .run(|| {
                assert!(crate::cache_source::has_active_scope(first.registry.id()));
                CacheSession::run_uncached(|| {
                    assert!(first.close().is_err());
                    Ok(())
                })?;
                second.run(|| {
                    assert!(first.close().is_err());
                    // A pre-existing outer query can rebind its owner while a
                    // closer is waiting; this must neither add a permit nor fail.
                    first.registry.begin_close();
                    first.run(|| {
                        cache_put("nested-owner", &[1], &Q::from(3));
                        Ok(())
                    })
                })?;
                Ok(())
            })
            .unwrap();
        assert!(!crate::cache_source::has_active_scope(first.registry.id()));
        first.close().unwrap();
        second.close().unwrap();
    }

    #[test]
    fn st_cache_uses_its_original_math_table_without_a_manual_value_mirror() {
        use crate::cache_source::{
            IntegralKey, ResultSource, with_optional_registry, with_registry,
        };

        let expected = with_optional_registry(None, || {
            MEMO.with(|memo| memo.borrow_mut().clear());
            let value = st_recur(&[4]);
            MEMO.with(|memo| memo.borrow_mut().clear());
            value
        });
        assert_eq!(expected, Q::from_signeds(1, 1152));
        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            assert_eq!(st_recur(&[4]), expected);
            let source = st_source(&registry);
            assert_eq!(source.get_slice(&[4]), Some(expected.clone()));
            let ready = source.ready();
            assert_eq!(st_recur(&[4]), expected);
            assert_eq!(source.ready(), ready);
            assert_eq!(st_recur(&[4, 1]), Q::from_signeds(1, 384));
            let ready = source.ready();
            assert_eq!(st_recur(&[1, 4]), Q::from_signeds(1, 384));
            assert_eq!(source.ready(), ready);
            assert!(MEMO.with(|memo| memo.borrow().is_empty()));
            assert_eq!(registry.manual_get("st", &[4]), None);
            assert!(source.pull(0, usize::MAX).unwrap().rows.iter().any(|row| {
                row.key
                    == IntegralKey::Legacy {
                        kind: "st".to_owned(),
                        nums: vec![4],
                    }
                    && row.value == expected
            }));
        });
    }

    #[test]
    fn st_infallible_queries_keep_exact_values_after_storage_failure() {
        use crate::cache_source::with_registry;

        let registry = Arc::new(Registry::new());
        registry.set_history(true);
        registry.fail("test ST storage failure".to_owned());
        with_registry(registry.clone(), || {
            assert_eq!(st_recur(&[4]), Q::from_signeds(1, 1152));
            assert_eq!(st_sum(&[1], &[]), Q::from_signeds(1, 24));
            assert_eq!(
                registry.check_error().unwrap_err(),
                "test ST storage failure"
            );
        });
    }

    #[test]
    fn empty_st_convention_is_not_published_as_a_stable_integral() {
        use crate::cache_source::{IntegralKey, with_registry};

        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            assert_eq!(st_recur(&[]), Q::ONE);
            assert_eq!(st_sum(&[], &[]), Q::ONE);
        });
        let rows: Vec<_> = registry
            .sources()
            .iter()
            .flat_map(|source| source.pull(0, usize::MAX).unwrap().rows)
            .collect();
        assert_eq!(rows.len(), 1);
        assert_eq!(
            rows[0].key,
            IntegralKey::Legacy {
                kind: "st".into(),
                nums: vec![]
            }
        );
        assert_eq!(rows[0].value, Q::ONE);
    }

    #[test]
    fn genus_zero_socle_publishes_the_root_and_actual_ct_descendants() {
        use crate::cache_source::{IntegralKey, with_registry};

        let registry = Arc::new(Registry::new());
        with_registry(registry.clone(), || {
            assert!(socle_formula(0, &[1, 1, 0, 0], &[]).is_err());
            assert!(registry.sources().is_empty());
            assert_eq!(socle_formula(0, &[0; 4], &[1]).unwrap(), Q::ONE);
            assert_eq!(st_sum(&[0; 4], &[1]), Q::ONE);
            // This helper convention may have a stable root, while klp still
            // rejects kappa_0 before looking at another API's root cache.
            assert_eq!(socle_formula(0, &[0; 3], &[0]).unwrap(), Q::ONE);
            assert!(crate::hodge::klp(0, 3, &[(0, 1)], &[], &[], &[]).is_err());
        });
        let rows: Vec<_> = registry
            .sources()
            .iter()
            .flat_map(|source| source.pull(0, usize::MAX).unwrap().rows)
            .collect();
        assert!(rows.iter().any(|row| {
            row.key
                == IntegralKey::Integral {
                    g: 0,
                    n: 4,
                    kappa: vec![(1, 1)],
                    psi: vec![],
                    lambda: vec![],
                    ch: vec![],
                }
                && row.value == Q::ONE
        }));
        assert!(rows.iter().any(|row| {
            row.key
                == IntegralKey::Descendant {
                    g: 0,
                    d: vec![0, 0, 0, 0, 2],
                    ch: vec![],
                }
                && row.value == Q::ONE
        }));
    }

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
