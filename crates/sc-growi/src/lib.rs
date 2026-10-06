//! # sc-growi
//!
//! Gromov-Witten 不变量计算引擎——Andreas Gathmann 的 GROWI (1.0.4) 的现代 Rust
//! 重实现核心库。设计与决策见工作区根目录 `docs/GROWI_REDESIGN.md`。
//!
//! 本 crate 只包含纯 Rust 的算法核心；CLI 与 pyo3 绑定是后续独立的 `sc-*` crate。
//!
//! 模块布局（对应原版文件，见 REDESIGN §2.3 映射表）：
//! - [`num`]     —— 有理数后端（malachite），对应 `math/rational.cc`
//! - [`poly`]    —— 变量/单项式/多项式，平坦 `Vec` 实现，对应 `math/{var,monomial,poly}.cc`
//! - [`inv`]     —— 不变量数据与规范形，对应 `libgw/inv.{h,cc}`
//! - [`intern`]  —— 不变量 hash-consing 仓，取代原版“排序 + 位打包 + gdbm”缓存键
//! - [`cache`]   —— 两级缓存（L1 内存 + L2 SQLite），取代 `database/database.cc`
//! - [`ctx`]     —— 计算上下文与策略链，对应 `libgw/compute.{h,cc}`
//! - [`simplify`] —— 标准公理链，对应 `libgw/simplify.cc`
//! - [`split2`]  —— 二分量分裂（WDVV / g=0 TRR），对应 `libgw/split2.cc`
//! - [`relinv`]  —— 相对分裂枚举状态机，对应 `libgw/relinv.{h,cc}`
//! - [`relsplit`]—— 相对分裂入口（超曲面绝对不变量），对应 `libgw/relsplit.cc`
//! - [`gw`]      —— GW() 入口（校验 + 多线性展开），对应 `libgw/gwinv.cc`
//!
//! M2/M3 将加入 `trr`（任意亏格 TRR，负亏格链不变量）与 `virasoro`。

pub(crate) fn env_flag(name: &str) -> bool {
    use std::sync::OnceLock;
    static CACHE: OnceLock<Vec<String>> = OnceLock::new();
    let vars = CACHE.get_or_init(|| {
        std::env::vars()
            .filter(|(k, _)| k.starts_with("SC_GROWI_"))
            .map(|(k, _)| k)
            .collect()
    });
    vars.iter().any(|k| k == name)
}

pub mod cache;
pub mod ctx;
pub mod gw;
pub mod intern;
pub mod inv;
pub mod linsys;
pub mod num;
pub mod poly;
pub mod query;
pub mod relinv;
pub mod relsplit;
pub mod simplify;
pub mod split2;
pub mod symfunc;
pub mod trr;
pub mod virasoro;

pub use ctx::Ctx;
pub use gw::{GwClass, GwError, GwSettings, gw};
pub use intern::{InvId, InvStore};
pub use inv::{Inv, Pdata};
pub use num::Q;
pub use poly::{Mono, Poly, Var};
