//! 对拍导出：把一批查询与本引擎的值打成 JSONL（供 tools/faber_diff.py 对照 admcycles）。
//! 用法：cargo run -p sc-faber --example diff_dump [-- 缓存.sqlite]

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

fn main() {
    if let Some(path) = std::env::args().nth(1) {
        let _ = sc_faber::socle::attach_cache(&path);
    }
    let mut out = String::new();
    // 基本类上的 κ/ψ 单项式（顶权）：(g, n, psi, kappa)
    // λ 差分批（M4：低次幂、g≤4——oracle legacy 路径秒级可算）
    let lam: Vec<(u32, u32, Vec<(u32, u32)>, Vec<(u32, u32)>, Vec<(u32, u32)>)> = vec![
        (1, 1, vec![], vec![(1, 1)], vec![(1, 1)]),
        (1, 2, vec![], vec![(1, 1), (2, 1)], vec![(1, 1)]),
        (1, 2, vec![(1, 1)], vec![], vec![(1, 1)]),
        (2, 0, vec![(1, 1)], vec![], vec![(2, 1)]),
        (2, 0, vec![(2, 1)], vec![], vec![(1, 1)]),
        (2, 0, vec![], vec![], vec![(1, 1), (2, 1)]),
        (2, 1, vec![], vec![(1, 3)], vec![(1, 1)]),
        (2, 1, vec![], vec![(1, 2)], vec![(2, 1)]),
        (2, 1, vec![(1, 1)], vec![(1, 1)], vec![(1, 1)]),
        (2, 2, vec![], vec![(1, 2), (2, 1)], vec![(2, 1)]),
        (2, 2, vec![], vec![(1, 2)], vec![(1, 1), (2, 1)]),
        (2, 2, vec![(2, 1)], vec![], vec![(1, 3)]),
        (3, 0, vec![(3, 1)], vec![], vec![(3, 1)]),
        (3, 0, vec![(1, 1)], vec![], vec![(2, 1), (3, 1)]),
        (3, 0, vec![], vec![], vec![(1, 1), (2, 1), (3, 1)]),
        (3, 0, vec![(3, 1)], vec![], vec![(1, 3)]),
        (3, 1, vec![], vec![(1, 3)], vec![(3, 1)]),
        (3, 1, vec![], vec![(1, 1)], vec![(2, 1), (3, 1)]),
        (3, 1, vec![(2, 1)], vec![(1, 1)], vec![(1, 1), (2, 1)]),
    ];
    for (g, n, kappa, psi, lambda) in &lam {
        let v = sc_faber::hodge::klp(*g, *n, kappa, psi, lambda, &[])
            .unwrap_or_else(|e| panic!("g={g} n={n}: {e}"));
        let psi_s: Vec<String> = psi.iter().map(|(l, p)| format!("[{l},{p}]")).collect();
        let kap_s: Vec<String> = kappa.iter().map(|(i, p)| format!("[{i},{p}]")).collect();
        let lam_s: Vec<String> = lambda.iter().map(|(i, p)| format!("[{i},{p}]")).collect();
        out.push_str(&format!(
            "{{\"kind\":\"lam\",\"g\":{g},\"n\":{n},\"psi\":[{}],\"kappa\":[{}],\"lambda\":[{}],\"value\":\"{v}\"}}
",
            psi_s.join(","),
            kap_s.join(","),
            lam_s.join(",")
        ));
    }
    let fund: Vec<(u32, u32, Vec<(u32, u32)>, Vec<(u32, u32)>)> = vec![
        (0, 4, vec![(1, 1)], vec![]),
        (0, 5, vec![(1, 2), (2, 1)], vec![]),
        (0, 4, vec![], vec![(1, 1)]),
        (0, 5, vec![], vec![(1, 2)]),
        (0, 5, vec![(1, 1)], vec![(1, 1)]),
        (1, 1, vec![(1, 1)], vec![]),
        (1, 2, vec![(1, 1), (2, 1)], vec![]),
        (1, 2, vec![(1, 2)], vec![]),
        (1, 3, vec![(2, 1), (3, 2)], vec![]),
        (1, 4, vec![(1, 2), (2, 2)], vec![]),
        (1, 4, vec![(2, 1), (4, 3)], vec![]),
        (1, 3, vec![(1, 3)], vec![]),
        (1, 5, vec![(4, 2), (5, 3)], vec![]),
        (2, 1, vec![(1, 4)], vec![]),
        (2, 0, vec![], vec![(1, 3)]),
        (2, 2, vec![(2, 5)], vec![]),
        (2, 2, vec![(1, 1), (2, 4)], vec![]),
        (2, 2, vec![(1, 2), (2, 3)], vec![]),
        (2, 3, vec![(2, 1), (3, 5)], vec![]),
        (2, 0, vec![], vec![(3, 1)]),
        (2, 0, vec![], vec![(1, 1), (2, 1)]),
        (2, 1, vec![(1, 3)], vec![(1, 1)]),
        (2, 2, vec![(1, 1), (2, 1)], vec![(1, 1)]),
        (3, 0, vec![], vec![(6, 1)]),
        (3, 0, vec![], vec![(1, 1), (5, 1)]),
        (3, 0, vec![], vec![(2, 1), (4, 1)]),
        (3, 0, vec![], vec![(1, 2), (4, 1)]),
        (3, 0, vec![], vec![(1, 6)]),
        (3, 1, vec![(1, 3)], vec![(1, 1), (2, 1)]),
    ];
    for (g, n, psi, kappa) in &fund {
        let v = taut_monomial(*g, *n, kappa, psi)
            .evaluate()
            .unwrap_or_else(|e| panic!("g={g} n={n}: {e}"));
        let psi_s: Vec<String> = psi.iter().map(|(l, p)| format!("[{l},{p}]")).collect();
        let kap_s: Vec<String> = kappa.iter().map(|(i, p)| format!("[{i},{p}]")).collect();
        out.push_str(&format!(
            "{{\"kind\":\"fund\",\"g\":{g},\"n\":{n},\"psi\":[{}],\"kappa\":[{}],\"value\":\"{v}\"}}\n",
            psi_s.join(","),
            kap_s.join(",")
        ));
    }
    // strata：不可约边界 Δ₀ 与分离边界
    let strata: Vec<(&str, u32, u32, Vec<u32>, Vec<Vec<u32>>, Vec<(u32, u32)>)> = vec![
        ("d0", 1, 1, vec![0], vec![vec![1, 2, 3]], vec![(2, 3)]),
        ("d0", 2, 0, vec![1], vec![vec![1, 2]], vec![(1, 2)]),
        ("d0", 2, 2, vec![1], vec![vec![1, 2, 3, 4]], vec![(3, 4)]),
        ("d0", 1, 3, vec![0], vec![vec![1, 2, 3, 4, 5]], vec![(4, 5)]),
        (
            "sep11",
            2,
            2,
            vec![1, 1],
            vec![vec![1, 3], vec![2, 4]],
            vec![(3, 4)],
        ),
    ];
    for (name, g, n, genera, legs, edges) in &strata {
        let graph = StableGraph::new(genera, legs, edges).unwrap();
        let cls = TautClass::new(vec![DecStratum {
            graph,
            poly: KP::one(),
        }]);
        let v = cls.evaluate().unwrap();
        out.push_str(&format!(
            "{{\"kind\":\"stratum\",\"name\":\"{name}\",\"g\":{g},\"n\":{n},\"value\":\"{v}\"}}\n"
        ));
    }
    print!("{out}");
}
