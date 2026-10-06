//! 探针：逐个算已知小不变量，定位数值 bug。
use sc_growi::{GwClass, GwSettings, Poly, Var, gw};

fn p2(d: i16, classes: Vec<GwClass>) -> Result<sc_growi::Q, sc_growi::GwError> {
    let mut ctx = sc_growi::Ctx::memory();
    let s = GwSettings {
        dim: 2,
        deg: 1,
        genus: 0,
        cdeg: d,
    };
    gw(&mut ctx, &s, &classes)
}

fn cls(poly: Poly, count: u32) -> GwClass {
    GwClass {
        poly,
        mult: 0,
        count,
    }
}

fn h(k: u32) -> Poly {
    Poly::from(Var::H).pow(k)
}

fn main() {
    let cases: Vec<(&str, i16, Vec<GwClass>, &str)> = vec![
        ("lines thru 2 pts", 1, vec![cls(h(2), 2)], "1"),
        ("conics thru 5 pts", 2, vec![cls(h(2), 5)], "1"),
        ("cubics thru 8 pts", 3, vec![cls(h(2), 8)], "12"),
        // 除子公理：<H, H2^5> d=2 = 2 * 1
        (
            "divisor H + 5pts d2",
            2,
            vec![cls(h(1), 1), cls(h(2), 5)],
            "2",
        ),
        // 两层除子：<H2^2, H^2> d=1 = 1
        ("2pts+2lines d1", 1, vec![cls(h(2), 2), cls(h(1), 2)], "1"),
        // <H^3> d=1 = 1（单点，projspace_1point）
        ("H3 1pt d1", 1, vec![cls(h(3), 1)], "1"),
        // <H^6> d=2 = 1
        ("H6 1pt d2", 2, vec![cls(h(6), 1)], "1"),
        // <H^9> d=3 = 1（3d 公式）
        ("H9 1pt d3", 3, vec![cls(h(9), 1)], "1"),
    ];
    for (name, d, classes, expect) in cases {
        match p2(d, classes) {
            Ok(v) => println!(
                "{name:<24} d={d}: {} (期望 {expect}) {}",
                v,
                if v.to_string() == expect {
                    "✓"
                } else {
                    "✗"
                }
            ),
            Err(e) => println!("{name:<24} d={d}: ERR {e}"),
        }
    }
}
