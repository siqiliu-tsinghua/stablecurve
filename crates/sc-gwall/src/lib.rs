//! sc-gwall：光滑曲线的（不连通）Gromov-Witten 不变量。
//!
//! Jim Bryan《GWall》（Maple，2002）的 Rust 移植（GROWI 作者是 Andreas
//! Gathmann——两个项目勿混）。算法来自 Pandharipande–Okounkov 与
//! Faber–Pandharipande（local curve 理论）。
//!
//! 模块：
//! - [`comb`]：combinat 子集（permute/cartprod/partition）+ 置换奇偶
//! - [`qseries`]：截断有理级数 + Bernoulli（Z/Y 配分函数组装）
//! - [`gwpoint`]：纯点后代闭式（Faber 幂和 + 钩长维数）
//! - [`gwloop`]：环路后代双实现（外代数符号）
//! - [`gwl`]：GW/GWhodge/GWomega 递归族 + gw 前端（Ins 七种插入）
//! - [`pf`]：branch/Ztil/Z/连通化 L/BPS 反演（BP 论文管线）

pub mod comb;
pub mod gwl;
pub mod gwloop;
pub mod gwpoint;
pub mod pf;
pub mod qseries;
