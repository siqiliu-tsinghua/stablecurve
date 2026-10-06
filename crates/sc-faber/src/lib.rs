//! sc-faber：曲线模空间 M̄_{g,n} 上重言代数的计算内核。
//!
//! 移植自 Carel Faber 的 Maple 程序《Numbers》（2001–2004），但以算法生成的
//! Faber–Zagier/Pixton 关系取代其静态预计算表（原版因此受限于 3g−3+n < 20）。
//! 功能对标 admcycles 的重言代数内核（见 docs/FABER_REDESIGN.md）。
//!
//! 分层：
//! - [`graph`] / [`kppoly`] / [`taut`]：稳定图、κ-ψ 多项式、重言类骨架；
//! - [`socle`] / [`hodge`] / [`push`] / [`int`]：socle 求值、λ/ch 类、
//!   Mumford 前推、积分（Faber 原程序核心）；
//! - [`dr`] / [`classes`]：DR cycle（Pixton/Chiodo）、Θ/Hain/超椭圆等
//!   （admcycles 对标增量）；
//! - [`witten`] / [`grr`]：Witten r-旋、广义 Hodge chern（GRRcomp）。

pub mod classes;
pub mod dr;
mod flow;
pub mod generators;
pub mod graph;
pub mod grr;
pub mod hodge;
pub mod int;
pub mod kppoly;
pub mod push;
pub mod socle;
pub mod tau;
pub mod taut;
pub mod witten;

pub use graph::StableGraph;
pub use int::{int_fund, int_fund_monomial};
pub use kppoly::KappaPsiPolynomial;
pub use sc_num::Q;
pub use tau::Tau;
pub use taut::{DecStratum, TautClass};
