//! 与 admcycles 收录的 Hodge 积分数据逐条对拍：M̄_g 上全部 λ 单项式的积分
//! （J. Schmitt, Z. Sun 2026，admcycles/data/hodge_data.sage）。夹具
//! tests/data/admcycles_hodge.txt 由 tools/admcycles_hodge.py 生成（g ≤ 9）；
//! 更高亏格（至 g = 14）用该工具的 check 子命令手动对拍。
//!
//! debug 构建只跑 g ≤ 7（保持 `cargo test` 快速），release 跑全部。

use sc_faber::hodge::klp;

#[test]
fn hodge_integrals_match_admcycles_data() {
    let data = include_str!("data/admcycles_hodge.txt");
    let max_g = if cfg!(debug_assertions) { 7 } else { u32::MAX };
    let (mut checked, mut nonzero) = (0, 0);
    for line in data.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut cols = line.split('\t');
        let g: u32 = cols.next().unwrap().parse().unwrap();
        if g > max_g {
            continue;
        }
        let mut lam: Vec<(u32, u32)> = Vec::new();
        for i in cols.next().unwrap().split(',') {
            let i: u32 = i.parse().unwrap();
            match lam.last_mut() {
                Some(e) if e.0 == i => e.1 += 1,
                _ => lam.push((i, 1)),
            }
        }
        let expect = cols.next().unwrap();
        let got = klp(g, 0, &[], &[], &lam, &[])
            .unwrap_or_else(|e| panic!("g={g} λ={lam:?}：{e}"))
            .to_string();
        assert_eq!(got, expect, "g={g} λ={lam:?}");
        checked += 1;
        if expect != "0" {
            nonzero += 1;
        }
    }
    let expected_rows = if cfg!(debug_assertions) { 432 } else { 2033 };
    assert_eq!(checked, expected_rows, "夹具条数与亏格上限不符");
    assert!(nonzero > checked * 9 / 10);
}
