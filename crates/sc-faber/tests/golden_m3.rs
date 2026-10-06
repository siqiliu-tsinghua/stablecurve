//! M3 黄金测试：34 项与 admcycles 逐项精确一致的值（2026-10-05 差分对拍）。
//! 生成自 tools/faber_diff.py 的首轮全绿运行。

use sc_faber::graph::StableGraph;
use sc_faber::kppoly::KappaPsiPolynomial as KP;
use sc_faber::taut::{DecStratum, TautClass};

fn taut_monomial(g: u32, n: u32, kappa: &[(u32, u32)], psi: &[(u32, u32)]) -> TautClass {
    let mut p = KP::one();
    for &(i, c) in kappa {
        p = p.mul(&KP::kappa_power(0, i, c));
    }
    for &(l, c) in psi {
        p = p.mul(&KP::psi_power(l, c));
    }
    let legs: Vec<u32> = (1..=n).collect();
    let graph = StableGraph::new(&[g], &[legs], &[]).unwrap();
    TautClass::new(vec![DecStratum { graph, poly: p }])
}

#[test]
fn fund_monomials_vs_admcycles() {
    let cases: Vec<(u32, u32, &[(u32, u32)], &[(u32, u32)], &str)> = vec![
        (0, 4, &[(1, 1)], &[], "1"),
        (0, 5, &[(1, 2), (2, 1)], &[], "0"),
        (0, 4, &[], &[(1, 1)], "1"),
        (0, 5, &[], &[(1, 2)], "5"),
        (0, 5, &[(1, 1)], &[(1, 1)], "3"),
        (1, 1, &[(1, 1)], &[], "1/24"),
        (1, 2, &[(1, 1), (2, 1)], &[], "1/24"),
        (1, 2, &[(1, 2)], &[], "1/24"),
        (1, 3, &[(2, 1), (3, 2)], &[], "1/12"),
        (1, 4, &[(1, 2), (2, 2)], &[], "1/6"),
        (1, 4, &[(2, 1), (4, 3)], &[], "1/8"),
        (1, 3, &[(1, 3)], &[], "1/24"),
        (1, 5, &[(4, 2), (5, 3)], &[], "7/24"),
        (2, 1, &[(1, 4)], &[], "1/1152"),
        (2, 0, &[], &[(1, 3)], "43/2880"),
        (2, 2, &[(2, 5)], &[], "1/1152"),
        (2, 2, &[(1, 1), (2, 4)], &[], "1/384"),
        (2, 2, &[(1, 2), (2, 3)], &[], "29/5760"),
        (2, 3, &[(2, 1), (3, 5)], &[], "1/288"),
        (2, 0, &[], &[(3, 1)], "1/1152"),
        (2, 0, &[], &[(1, 1), (2, 1)], "1/240"),
        (2, 1, &[(1, 3)], &[(1, 1)], "29/5760"),
        (2, 2, &[(1, 1), (2, 1)], &[(1, 1)], "0"),
        (3, 0, &[], &[(6, 1)], "1/82944"),
        (3, 0, &[], &[(1, 1), (5, 1)], "1/5760"),
        (3, 0, &[], &[(2, 1), (4, 1)], "971/2903040"),
        (3, 0, &[], &[(1, 2), (4, 1)], "2173/967680"),
        (3, 0, &[], &[(1, 6)], "176557/107520"),
        (3, 1, &[(1, 3)], &[(1, 1), (2, 1)], "0"),
    ];
    for (g, n, psi, kappa, expect) in &cases {
        let v = taut_monomial(*g, *n, kappa, psi)
            .evaluate()
            .unwrap()
            .to_string();
        assert_eq!(&v, expect, "g={g} n={n} psi={psi:?} kappa={kappa:?}");
    }
}

#[test]
fn strata_vs_admcycles() {
    let cases: Vec<(&str, Vec<u32>, Vec<Vec<u32>>, Vec<(u32, u32)>, &str)> = vec![
        ("d0", vec![0], vec![vec![1, 2, 3]], vec![(2, 3)], "1/2"),
        ("d0", vec![1], vec![vec![1, 2]], vec![(1, 2)], "0"),
        ("d0", vec![1], vec![vec![1, 2, 3, 4]], vec![(3, 4)], "0"),
        ("d0", vec![0], vec![vec![1, 2, 3, 4, 5]], vec![(4, 5)], "0"),
        (
            "sep11",
            vec![1, 1],
            vec![vec![1, 3], vec![2, 4]],
            vec![(3, 4)],
            "0",
        ),
    ];
    for (name, genera, legs, edges, expect) in &cases {
        let graph = StableGraph::new(genera, legs, edges).unwrap();
        let cls = TautClass::new(vec![DecStratum {
            graph,
            poly: KP::one(),
        }]);
        let v = cls.evaluate().unwrap().to_string();
        assert_eq!(&v, expect, "stratum {name}");
    }
}
