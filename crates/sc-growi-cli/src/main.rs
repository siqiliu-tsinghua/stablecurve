//! sc-growi 引擎二进制。
//!
//! 用法：
//!   growi "<查询表达式>" [--cache PATH]          单查询，值打印到 stdout
//!   growi --batch [--cache PATH] < queries.txt  批处理：stdin 每行一个表达式，stdout 逐行 JSON
//!   growi "lines in quintic threefold"        → 2875
//!   growi "G=2,D=5 in P^4 thru H^3*psi^22"    → -41369/110075314176
//!
//! 缓存缺省为 ~/.growi/cache.sqlite；`--no-cache` 只用内存。
//! 表达式语法是原版 `expr.cc` 关键字表的直译（见 `sc_growi::query`）。
//! 批处理应答：{"line":k,"status":"ok","value":"p/q"}、{"line":k,"status":"unknown"}
//! （当前策略链算不出，原版 UNDEFINED）或 {"line":k,"status":"invalid","message":"..."}。

use sc_growi::query::parse;
use sc_growi::{Ctx, GwError, gw};
use serde_json::json;

const USAGE: &str = "用法: growi \"<查询表达式>\" [--cache PATH | --no-cache]\n      growi --batch [--cache PATH | --no-cache] < queries.txt";

fn open_ctx(cache: Option<&str>, no_cache: bool) -> Result<Ctx, String> {
    if no_cache {
        return Ok(Ctx::memory());
    }
    match cache {
        Some(p) => {
            Ctx::with_db(std::path::Path::new(p)).map_err(|e| format!("无法打开缓存 {p}: {e}"))
        }
        None => {
            // 默认持久化：~/.growi/cache.sqlite（两级缓存的 L2）
            let home = std::env::var_os("HOME")
                .or_else(|| std::env::var_os("USERPROFILE"))
                .ok_or("无法确定主目录")?;
            let dir = std::path::PathBuf::from(home).join(".growi");
            std::fs::create_dir_all(&dir)
                .map_err(|e| format!("无法创建 {}: {e}", dir.display()))?;
            Ctx::with_db(&dir.join("cache.sqlite")).map_err(|e| format!("无法打开默认缓存: {e}"))
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut expr: Option<&str> = None;
    let mut cache: Option<&str> = None;
    let mut no_cache = false;
    let mut batch = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--cache" => {
                i += 1;
                let Some(p) = args.get(i) else {
                    eprintln!("--cache 需要路径\n{USAGE}");
                    std::process::exit(2);
                };
                cache = Some(p);
            }
            "--no-cache" => no_cache = true,
            "--batch" => batch = true,
            "--help" | "-h" => {
                println!("{USAGE}");
                return;
            }
            a if a.starts_with("--") => {
                eprintln!("未知选项 {a}\n{USAGE}");
                std::process::exit(2);
            }
            a => expr = Some(a),
        }
        i += 1;
    }
    let mut ctx = match open_ctx(cache, no_cache) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(2);
        }
    };
    if batch {
        run_batch(&mut ctx);
        return;
    }
    let Some(expr) = expr else {
        eprintln!("{USAGE}");
        std::process::exit(2);
    };
    let query = match parse(expr) {
        Ok(q) => q,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(2);
        }
    };
    match gw(&mut ctx, &query.settings, &query.classes) {
        Ok(v) => println!("{v}"),
        Err(GwError::Unknown(iv)) => println!("UNDEFINED ({iv})"),
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(2);
        }
    }
}

/// 批处理：stdin 每行一个表达式，stdout 逐行 JSON。
fn run_batch(ctx: &mut Ctx) {
    use std::io::BufRead;
    for (i, line) in std::io::stdin().lock().lines().enumerate() {
        let line_no = i + 1;
        let out = match line {
            Err(e) => {
                json!({"line": line_no, "status": "invalid", "message": format!("读取失败：{e}")})
            }
            Ok(line) if line.trim().is_empty() => continue,
            Ok(line) => match parse(line.trim()) {
                Err(e) => json!({"line": line_no, "status": "invalid", "message": e}),
                Ok(q) => match gw(ctx, &q.settings, &q.classes) {
                    Ok(v) => json!({"line": line_no, "status": "ok", "value": v.to_string()}),
                    Err(GwError::Unknown(_)) => json!({"line": line_no, "status": "unknown"}),
                    Err(e) => {
                        json!({"line": line_no, "status": "invalid", "message": e.to_string()})
                    }
                },
            },
        };
        println!("{out}");
    }
}
