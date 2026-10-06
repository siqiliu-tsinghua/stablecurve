//! 计算上下文与策略链驱动——对应原版 `libgw/compute.{h,cc}`。
//!
//! 收编原版全部全局状态（GWdbm、缺省 N/D/G/E）。M1 阶段策略链为
//! simplify → split2 → relsplit（亏格 0 全覆盖）；trr / virasoro 在 M2/M3
//! 插到 relsplit 之前，顺序与原版一致（REDESIGN §1.4）。

use std::path::Path;

use crate::cache::Cache;
use crate::inv::Inv;
use crate::num::Q;

pub struct Ctx {
    pub cache: Cache,
    /// 递归深度保护（调试用；正常计算远达不到）
    depth: u32,
}

impl Ctx {
    pub fn new(cache: Cache) -> Self {
        Ctx { cache, depth: 0 }
    }

    pub fn memory() -> Self {
        Ctx::new(Cache::memory())
    }

    pub fn with_db(path: &Path) -> rusqlite::Result<Self> {
        Ok(Ctx::new(Cache::with_db(path)?))
    }

    pub(crate) fn depth(&self) -> u32 {
        self.depth
    }

    pub fn flush(&mut self) {
        self.cache.flush();
    }

    /// 递归计算一个不变量。`None` = 当前策略链算不出（原版 UNDEFINED）。
    pub fn compute(&mut self, inv: &Inv) -> Option<Q> {
        if crate::env_flag("SC_GROWI_CHAIN") {
            eprintln!("[chain {}] {inv}", self.depth);
        }
        self.depth += 1;
        if self.depth > 10_000 {
            self.depth -= 1;
            if crate::env_flag("SC_GROWI_CYCLE") {
                eprintln!("[cycle d={}] {inv}", self.depth);
            }
            return None;
        }
        let mut iv = inv.clone();
        let r = self.compute_sorted(&mut iv);
        self.depth -= 1;
        r
    }

    /// 内部：输入可能未排序；simplify 会先规范形化。
    /// 策略链与原版 compute 一致：simplify → split2 → trr → [virasoro(M3)] → relsplit。
    fn compute_sorted(&mut self, iv: &mut Inv) -> Option<Q> {
        if let Some(r) = crate::simplify::simplify(self, iv) {
            return Some(r);
        }
        if let Some(r) = crate::split2::split2(self, iv) {
            return Some(r);
        }
        if let Some(r) = crate::trr::trr_genusg(self, iv) {
            return Some(r);
        }
        if let Some(r) = crate::virasoro::virasoro(self, iv) {
            return Some(r);
        }
        if let Some(r) = crate::relsplit::relsplit(self, iv) {
            return Some(r);
        }
        eprintln!("*error* don't know how to compute {iv}");
        None
    }
}
