//! sc-gwall 引擎二进制。
//!
//! 用法：
//!   gwall --batch < queries.jsonl          批处理（stdin 逐行 JSON，stdout 逐行 JSON）
//!   gwall '{"op":"gw",...}'                单查询（值直接打印到 stdout）
//!
//! 查询协议：
//!   {"op":"gw","g":G,"d":D,"insert":[..],"ldeg":L}
//!       insert 元素两种写法等价：
//!         数组  ["t1",K] | ["pt",K] | ["a",K,I] | ["b",K,I] | ["chO",K] | ["chL",K] | ["ldeg",V]
//!         对象  {"t":"t1","k":K} | {"t":"a","k":K,"i":I} | {"t":"ldeg","v":V} | …
//!       （对应 Maple 编码 k / −k / [k,i] / [k,−i] / [k] / [−k] / [[deg]]；
//!         a/b 环路可按任意顺序给出，符号自动计入）
//!   {"op":"ztil","g":G,"d":D,"n":N}        t 级数系数（t^0..t^N）
//!   {"op":"z","g":G,"d":D,"n":N}           去常映射（exp(−N₀) 修正）
//!   {"op":"l","g":G,"d":D,"n":N}           度变量 log 连通化 L_d
//!   {"op":"bps","g":G,"d":D,"h":H}         local BPS n_d^h(g)（GV 反演）
//! 批处理应答：{"line":k,"status":"ok","value":"p/q"}（级数 op 的 value 为字符串数组）
//!             或 {"line":k,"status":"error","error":"..."}

use sc_gwall::gwl::{Ctx, Ins, gw_ctx};
use sc_gwall::pf::{bps_local, l_series, z, ztil};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Query {
    op: String,
    g: u32,
    #[serde(default)]
    d: Option<u32>,
    #[serde(default)]
    n: Option<u32>,
    #[serde(default)]
    h: Option<u32>,
    #[serde(default)]
    insert: Vec<Value>,
    #[serde(default)]
    ldeg: Option<i64>,
}

/// 计算结果：单值或级数系数表。
enum Answer {
    Scalar(String),
    Series(Vec<String>),
}

impl Answer {
    fn to_json(&self) -> Value {
        match self {
            Answer::Scalar(s) => json!(s),
            Answer::Series(v) => json!(v),
        }
    }

    fn to_text(&self) -> String {
        match self {
            Answer::Scalar(s) => s.clone(),
            Answer::Series(v) => format!("[{}]", v.join(",")),
        }
    }
}

fn get_u32(v: Option<&Value>, what: &str) -> Result<u32, String> {
    v.and_then(Value::as_u64)
        .and_then(|x| u32::try_from(x).ok())
        .ok_or_else(|| format!("{what} 须为非负整数"))
}

fn get_i64(v: Option<&Value>, what: &str) -> Result<i64, String> {
    v.and_then(Value::as_i64)
        .ok_or_else(|| format!("{what} 须为整数"))
}

/// 一个插入项：数组 `["a", K, I]` 或对象 `{"t":"a","k":K,"i":I}`。
fn parse_insert(x: &Value) -> Result<Ins, String> {
    let (t, a1, a2): (&str, Option<&Value>, Option<&Value>) = match x {
        Value::Array(arr) => {
            let t = arr
                .first()
                .and_then(Value::as_str)
                .ok_or("插入项数组须以类型字符串开头")?;
            let want = match t {
                "a" | "b" => 3,
                _ => 2,
            };
            if arr.len() != want {
                return Err(format!("插入项 {x} 应有 {want} 个元素"));
            }
            (t, arr.get(1), arr.get(2))
        }
        Value::Object(obj) => {
            let t = obj
                .get("t")
                .and_then(Value::as_str)
                .ok_or("插入项对象缺少字符串字段 t")?;
            let first = if t == "ldeg" {
                obj.get("v")
            } else {
                obj.get("k")
            };
            (t, first, obj.get("i"))
        }
        _ => return Err(format!("插入项 {x} 须为数组或对象")),
    };
    Ok(match t {
        "t1" => Ins::T1(get_u32(a1, "t1 的 k")?),
        "pt" => Ins::Pt(get_u32(a1, "pt 的 k")?),
        "a" => Ins::A(get_u32(a1, "a 的 k")?, get_u32(a2, "a 的 i")?),
        "b" => Ins::B(get_u32(a1, "b 的 k")?, get_u32(a2, "b 的 i")?),
        "chO" => Ins::ChO(get_u32(a1, "chO 的 k")?),
        "chL" => Ins::ChL(get_u32(a1, "chL 的 k")?),
        "ldeg" => Ins::LDeg(get_i64(a1, "ldeg 的 v")?),
        other => {
            return Err(format!(
                "未知插入类型 '{other}'（可用：t1/pt/a/b/chO/chL/ldeg）"
            ));
        }
    })
}

fn run(ctx: &mut Ctx, q: &Query) -> Result<Answer, String> {
    let d = q.d.ok_or("缺少 d（曲线次数）")?;
    let series = |c: &[sc_num::Q], n: u32| -> Answer {
        Answer::Series(
            (0..=n as usize)
                .map(|i| {
                    c.get(i)
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "0".into())
                })
                .collect(),
        )
    };
    match q.op.as_str() {
        "gw" => {
            let mut ins = q
                .insert
                .iter()
                .map(parse_insert)
                .collect::<Result<Vec<_>, _>>()?;
            if let Some(l) = q.ldeg {
                ins.push(Ins::LDeg(l));
            }
            Ok(Answer::Scalar(gw_ctx(ctx, q.g, d, &ins).to_string()))
        }
        "ztil" => {
            let n = q.n.ok_or("ztil 需要 n")?;
            Ok(series(&ztil(ctx, q.g, d, n)?.c, n))
        }
        "z" => {
            let n = q.n.ok_or("z 需要 n")?;
            Ok(series(&z(ctx, q.g, d, n)?.c, n))
        }
        "l" => {
            let n = q.n.ok_or("l 需要 n")?;
            Ok(series(&l_series(ctx, q.g, d, n)?.c, n))
        }
        "bps" => {
            let h = q.h.ok_or("bps 需要 h")?;
            Ok(Answer::Scalar(bps_local(ctx, q.g, d, h)?.to_string()))
        }
        other => Err(format!("未知操作 '{other}'（可用：gw/ztil/z/l/bps）")),
    }
}

fn answer(ctx: &mut Ctx, line: &str) -> Result<Answer, String> {
    let q: Query = serde_json::from_str(line).map_err(|e| format!("JSON 解析失败：{e}"))?;
    run(ctx, &q)
}

const USAGE: &str =
    "用法：gwall --batch < queries.jsonl\n      gwall '{\"op\":\"bps\",\"g\":3,\"d\":2,\"h\":1}'";

fn main() {
    let mut batch = false;
    let mut single: Option<String> = None;
    for a in std::env::args().skip(1) {
        match a.as_str() {
            "--batch" => batch = true,
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
    }
    let mut ctx = Ctx::new();
    if batch {
        use std::io::BufRead;
        for (k, line) in std::io::stdin().lock().lines().enumerate() {
            let resp = match line {
                Err(e) => {
                    json!({"line": k + 1, "status": "error", "error": format!("读取失败：{e}")})
                }
                Ok(line) if line.trim().is_empty() => continue,
                Ok(line) => match answer(&mut ctx, &line) {
                    Ok(v) => json!({"line": k + 1, "status": "ok", "value": v.to_json()}),
                    Err(e) => json!({"line": k + 1, "status": "error", "error": e}),
                },
            };
            println!("{resp}");
        }
    } else if let Some(s) = single {
        match answer(&mut ctx, &s) {
            Ok(v) => println!("{}", v.to_text()),
            Err(e) => {
                eprintln!("错误：{e}");
                std::process::exit(2);
            }
        }
    } else {
        eprintln!("{USAGE}");
        std::process::exit(2);
    }
}
