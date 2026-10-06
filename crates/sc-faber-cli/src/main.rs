//! sc-faber 引擎二进制。
//!
//! 用法：
//!   faber --batch [--cache PATH] < queries.jsonl     批处理（stdin 逐行 JSON，stdout 逐行 JSON）
//!   faber '{"op":"klp",...}' [--cache PATH]          单查询（值直接打印到 stdout）
//!
//! 查询协议（每行一个 JSON 对象；κ/ψ/λ/ch 均为 [[指标或腿号, 幂],..]）：
//!   {"op":"klp","g":G,"n":N,"kappa":..,"psi":..,"lambda":..,"ch":..}   κ/ψ/λ/ch 单项式积分
//!   {"op":"dr","g":G,"a":[..],"r":R,"kappa":..,"psi":..}               Pixton DR（r 缺省 = g）
//!   {"op":"drpoly","g":G,"a":[..],"k":K,...}                          Chiodo 插值主路径
//!   {"op":"drspin","g":G,"a":[..],"k":K,...}                          spin DR（A 全奇、k 奇）
//!   {"op":"witten","g":G,"a3":[[cn,cd,d],..],...}                     Witten r-旋极限类
//!   {"op":"hyperell"|"biell","g":G,"n":N,"m":M,...}                   超椭圆/双椭圆轨迹
//!   {"op":"grr","g":G,"n":N,"t":T,"l":L,"dvec":[..],"a3":[[h,mask,v],..],...}  广义 Hodge chern
//! drpoly/drspin 的 k 缺省为 ΣA/(2g−2+n)（须整除）。
//! 批处理应答：{"line":k,"status":"ok","value":"p/q"} 或 {"line":k,"status":"error","error":"..."}

use sc_faber::kppoly::check_marking_legs;
use serde::Deserialize;
use serde_json::json;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Query {
    op: String,
    g: u32,
    #[serde(default)]
    n: Option<u32>,
    #[serde(default)]
    m: Option<u32>,
    #[serde(default)]
    k: Option<i64>,
    #[serde(default)]
    a3: Option<Vec<Vec<i64>>>,
    #[serde(default)]
    t: Option<u32>,
    #[serde(default)]
    l: Option<i64>,
    #[serde(default)]
    dvec: Option<Vec<i64>>,
    #[serde(default)]
    a: Option<Vec<i64>>,
    #[serde(default)]
    r: Option<u32>,
    #[serde(default)]
    kappa: Vec<(u32, u32)>,
    #[serde(default)]
    psi: Vec<(u32, u32)>,
    #[serde(default)]
    lambda: Vec<(u32, u32)>,
    #[serde(default)]
    ch: Vec<(u32, u32)>,
}

const OPS: &str = "klp / dr / drpoly / drspin / witten / hyperell / biell / grr";

/// drpoly/drspin 的 k：显式给出，或 ΣA/(2g−2+n)（须整除）。
fn twist(q: &Query, a: &[i64]) -> Result<i64, String> {
    let chi = 2 * q.g as i64 - 2 + a.len() as i64;
    if chi <= 0 {
        return Err(format!("M̄_({},{}) 不稳定", q.g, a.len()));
    }
    if let Some(k) = q.k {
        return Ok(k);
    }
    let s: i64 = a.iter().sum();
    if s % chi != 0 {
        return Err(format!("ΣA={s} 不是 2g−2+n={chi} 的倍数，请显式给出 k"));
    }
    Ok(s / chi)
}

fn run(q: &Query) -> Result<String, String> {
    use sc_faber::classes::{dr_chiodo_poly, int_against};
    match q.op.as_str() {
        "klp" => {
            let n = q.n.ok_or("klp 需要 n")?;
            sc_faber::hodge::klp(q.g, n, &q.kappa, &q.psi, &q.lambda, &q.ch).map(|v| v.to_string())
        }
        "dr" => {
            let a = q.a.as_ref().ok_or("dr 需要 a")?;
            sc_faber::dr::int_dr(q.g, a, q.r, &q.kappa, &q.psi).map(|v| v.to_string())
        }
        "drpoly" | "drspin" => {
            let a = q.a.as_ref().ok_or_else(|| format!("{} 需要 a", q.op))?;
            check_marking_legs(&q.psi, a.len() as u32)?;
            let k = twist(q, a)?;
            let cls = dr_chiodo_poly(q.g, a, q.g, k, q.op == "drspin")?;
            int_against(&cls, &q.kappa, &q.psi).map(|v| v.to_string())
        }
        "grr" => {
            let n = q.n.ok_or("grr 需要 n")?;
            let t = q.t.ok_or("grr 需要 t")?;
            let l = q.l.unwrap_or(0);
            let d = q.dvec.as_ref().ok_or("grr 需要 dvec")?;
            check_marking_legs(&q.psi, n)?;
            let mut a: Vec<(u32, Vec<u32>, i64)> = Vec::new();
            for tr in q.a3.iter().flatten() {
                let [h, mask, v] = tr[..] else {
                    return Err("grr 的 a3 项须为 [h, s_mask, v]".into());
                };
                let h = u32::try_from(h).map_err(|_| format!("grr 的 h={h} 须为非负"))?;
                let mask =
                    u64::try_from(mask).map_err(|_| format!("grr 的 s_mask={mask} 须为非负"))?;
                let s: Vec<u32> = (0..64)
                    .filter(|i| mask >> i & 1 == 1)
                    .map(|i| i + 1)
                    .collect();
                a.push((h, s, v));
            }
            let cls = sc_faber::grr::generalized_hodge_chern_single(t, l, d, &a, q.g, n)?;
            int_against(&cls, &q.kappa, &q.psi).map(|v| v.to_string())
        }
        "witten" => {
            let a3 = q.a3.as_ref().ok_or("witten 需要 a3=[[cn,cd,d],..]")?;
            check_marking_legs(&q.psi, a3.len() as u32)?;
            let mut a: Vec<(i64, i64, i64)> = Vec::new();
            for tr in a3 {
                let [cn, cd, d] = tr[..] else {
                    return Err("witten 的 a3 项须为 [cn, cd, d]".into());
                };
                a.push((cn, cd, d));
            }
            let cls = sc_faber::witten::witten_class(q.g, &a)?;
            int_against(&cls, &q.kappa, &q.psi).map(|v| v.to_string())
        }
        "hyperell" => {
            let n = q.n.ok_or("hyperell 需要 n")?;
            sc_faber::classes::hyperell_int(q.g, n, q.m.unwrap_or(0), &q.kappa, &q.psi)
                .map(|v| v.to_string())
        }
        "biell" => {
            let n = q.n.ok_or("biell 需要 n")?;
            sc_faber::classes::biell_int(q.g, n, q.m.unwrap_or(0), &q.kappa, &q.psi)
                .map(|v| v.to_string())
        }
        other => Err(format!("未知操作 '{other}'（可用：{OPS}）")),
    }
}

fn answer(line: &str) -> Result<String, String> {
    let q: Query = serde_json::from_str(line).map_err(|e| format!("JSON 解析失败：{e}"))?;
    run(&q)
}

const USAGE: &str = "用法：faber --batch [--cache PATH] < queries.jsonl\n      faber '{\"op\":\"klp\",\"g\":2,\"n\":1,\"psi\":[[1,4]]}' [--cache PATH]";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut batch = false;
    let mut cache: Option<String> = None;
    let mut single: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--batch" => batch = true,
            "--cache" => {
                i += 1;
                let Some(p) = args.get(i) else {
                    eprintln!("--cache 需要路径\n{USAGE}");
                    std::process::exit(2);
                };
                cache = Some(p.clone());
            }
            "--help" | "-h" => {
                println!("{USAGE}");
                return;
            }
            s if s.starts_with("--") => {
                eprintln!("未知选项 {s}\n{USAGE}");
                std::process::exit(2);
            }
            s => single = Some(s.to_string()),
        }
        i += 1;
    }
    if let Some(c) = &cache {
        if let Err(e) = sc_faber::socle::attach_cache(c) {
            eprintln!("warning: 缓存不可用（{e}），继续无 L2");
        }
    }
    if let (false, Some(line)) = (batch, &single) {
        match answer(line) {
            Ok(v) => println!("{v}"),
            Err(e) => {
                eprintln!("error: {e}");
                std::process::exit(2);
            }
        }
        return;
    }
    use std::io::BufRead;
    for (k, line) in std::io::stdin().lock().lines().enumerate() {
        let resp = match line {
            Err(e) => json!({"line": k + 1, "status": "error", "error": format!("读取失败：{e}")}),
            Ok(line) if line.trim().is_empty() => continue,
            Ok(line) => match answer(&line) {
                Ok(v) => json!({"line": k + 1, "status": "ok", "value": v}),
                Err(e) => json!({"line": k + 1, "status": "error", "error": e}),
            },
        };
        println!("{resp}");
    }
}
