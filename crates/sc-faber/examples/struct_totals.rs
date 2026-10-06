//! 逐结构贡献对拍：同一键函数 + 同一积分器，只有系数来源不同。
//! Part A：解析 oracle 系数表行 → 结构键 → Σ coeff × I(oracle 装饰) × aut。
//! Part B：我们的 71 生成元 → 同键 → Σ coeff' × I'(我们的装饰求值)。
//! 输出两列总额，python 侧合并对比。

use sc_faber::classes::{dr_coeff_one, int_against};
use sc_faber::generators::tautgens;
use sc_faber::graph::StableGraph;
use sc_faber::kppoly::{KappaPsiPolynomial as KP, KpTerm};
use sc_faber::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero, is_zero};

/// 结构键（不含装饰）：排序顶点 (genus, 标记集) + 边对多重集（按排序后位序）。
fn struct_key(graph: &StableGraph) -> String {
    let n = graph.num_markings();
    let plain: Vec<(u32, Vec<u32>)> = graph
        .vertices()
        .iter()
        .map(|v| {
            let mut m: Vec<u32> = v.legs.iter().copied().filter(|l| *l <= n).collect();
            m.sort_unstable();
            (v.genus, m)
        })
        .collect();
    let order: Vec<usize> = (0..plain.len()).collect();
    let order = {
        let mut o = order;
        o.sort_by_key(|&i| format!("{:?}", plain[i]));
        o
    };
    let rank: std::collections::HashMap<usize, usize> =
        order.iter().enumerate().map(|(r, &v)| (v, r)).collect();
    let pos = |l: u32| {
        graph
            .vertices()
            .iter()
            .position(|v| v.legs.contains(&l))
            .unwrap()
    };
    let mut epairs: Vec<(usize, usize)> = graph
        .edges()
        .iter()
        .map(|&(a, b)| {
            let (x, y) = (rank[&pos(a)], rank[&pos(b)]);
            (x.min(y), y.max(x))
        })
        .collect();
    epairs.sort_unstable();
    let mut vs: Vec<String> = plain
        .iter()
        .map(|(g, m)| format!("(g{g},m{m:?})"))
        .collect();
    vs.sort();
    format!("{};e{:?}", vs.join(";"), epairs)
}

fn parse_line(
    l: &str,
) -> Option<(
    Vec<u32>,
    Vec<Vec<u32>>,
    Vec<(u32, u32)>,
    Vec<Vec<u32>>,
    Vec<(u32, u32)>,
    String,
)> {
    let l = l.strip_prefix("G ")?.trim();
    fn take_group(s: &str) -> Option<(&str, &str)> {
        let s = s.trim_start();
        let open = s.chars().next()?;
        let close = match open {
            '[' => ']',
            '{' => '}',
            _ => return None,
        };
        let mut depth = 0;
        for (i, c) in s.char_indices() {
            if c == open {
                depth += 1;
            } else if c == close {
                depth -= 1;
                if depth == 0 {
                    return Some((&s[..=i], s[i + 1..].trim_start()));
                }
            }
        }
        None
    }
    let (g_s, r) = take_group(l)?;
    let (legs_s, r) = take_group(r)?;
    let (edges_s, r) = take_group(r)?;
    let (kappa_s, r) = take_group(r)?;
    let (psi_s, coeff_s) = take_group(r)?;
    let numlist = |s: &str| -> Vec<Vec<u32>> {
        s[1..s.len() - 1]
            .split("], [")
            .map(|p| {
                p.replace(['[', ']'], "")
                    .split(",")
                    .map(|x| x.trim())
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect()
            })
            .collect()
    };
    let genera: Vec<u32> = g_s[1..g_s.len() - 1]
        .split(",")
        .map(|x| x.trim())
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().unwrap())
        .collect();
    let legs = numlist(legs_s);
    let edges: Vec<(u32, u32)> = if edges_s == "[]" {
        vec![]
    } else {
        edges_s[1..edges_s.len() - 1]
            .split("), (")
            .map(|p| {
                let toks: Vec<u32> = p
                    .replace(['(', ')'], "")
                    .split(",")
                    .map(|x| x.trim())
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect();
                (toks[0], toks[1])
            })
            .collect()
    };
    let kappa_counts = numlist(kappa_s);
    let kappa_part: Vec<Vec<u32>> = kappa_counts
        .iter()
        .map(|counts| {
            let mut p = Vec::new();
            for (i, &c) in counts.iter().enumerate() {
                for _ in 0..c {
                    p.push(i as u32 + 1);
                }
            }
            p
        })
        .collect();
    let psi: Vec<(u32, u32)> = if psi_s == "{}" {
        vec![]
    } else {
        psi_s[1..psi_s.len() - 1]
            .split(", ")
            .map(|kv| {
                let (a, b) = kv.split_once(": ").unwrap();
                (a.trim().parse().unwrap(), b.trim().parse().unwrap())
            })
            .collect()
    };
    Some((genera, legs, edges, kappa_part, psi, coeff_s.to_string()))
}

fn main() {
    let mode = std::env::args().nth(1).expect("模式：oracle|ours [表路径]");
    let mut totals: std::collections::BTreeMap<String, Q> = std::collections::BTreeMap::new();
    if mode == "oracle" {
        let path = std::env::args().nth(2).expect("表路径");
        for line in std::fs::read_to_string(&path).unwrap().lines() {
            let Some((genera, legs, edges, kappa, psi, coeff_s)) = parse_line(line) else {
                continue;
            };
            let coeff = {
                let (a, b) = coeff_s.split_once('/').unwrap_or((coeff_s.as_str(), "1"));
                Q::from_signeds(a.parse::<i64>().unwrap(), b.parse::<i64>().unwrap())
            };
            let Ok(graph) = StableGraph::new_raw(&genera, &legs, &edges) else {
                continue;
            };
            let nv = genera.len();
            let mut kk = kappa;
            kk.resize(nv, Vec::new());
            let cls = TautClass::new(vec![DecStratum {
                graph: graph.clone(),
                poly: KP::new(vec![KpTerm {
                    kappa: kk,
                    psi,
                    coeff: Q::ONE,
                }]),
            }]);
            // oracle 口径：ξ_*dec 不除 aut → 我们的 int_against(已除 aut) 乘回 aut
            let aut = graph.automorphism_number();
            let i_val = int_against(&cls, &[], &[(1, 1)]).unwrap_or_default() * Q::from(aut as i64);
            let key = struct_key(&graph);
            *totals.entry(key).or_insert(Q::ZERO) += coeff * i_val;
        }
    } else {
        for generator in tautgens(2, 1, 3) {
            let Some(c) = dr_coeff_one(&generator, 2, &[1], 3, -1, 2) else {
                continue;
            };
            if is_zero(&c) {
                continue;
            }
            // dr_coeff_one 已含 1/aut（oracle 口径）→ 积分乘回 aut 到 ξ_* 口径
            let cls = TautClass::new(vec![generator.clone()]);
            let aut = generator.graph.automorphism_number();
            let i_val = int_against(&cls, &[], &[(1, 1)]).unwrap_or_default() * Q::from(aut as i64);
            let key = struct_key(&generator.graph);
            *totals.entry(key).or_insert(Q::ZERO) += c * i_val;
        }
    }
    for (k, v) in &totals {
        if !is_zero(v) {
            println!("{k}\t{v}");
        }
    }
}
// （main 之外）Δ₀ 结构逐装饰三向：MODE=delta0 时输出每装饰的 (oracle 系数, 我们积分, 我们系数)
