//! growi 二进制的协议测试：批处理每行都是合法 JSON。

use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

fn batch(input: &str) -> Vec<Value> {
    let mut child = Command::new(env!("CARGO_BIN_EXE_growi"))
        .args(["--batch", "--no-cache"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success());
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("非法 JSON 行 {l:?}: {e}")))
        .collect()
}

#[test]
fn batch_statuses() {
    let r = batch(concat!(
        "lines in cubic surface\n",
        "\n",
        "lines in quintic 3fold\n",
        "lines in P^2 thru H^2 \"quoted\" fold\n",
        "elliptic D=2 in quintic threefold\n",
        "D=65537 in P^2 thru H^2:2\n",
    ));
    assert_eq!(r.len(), 5);
    assert_eq!(
        (r[0]["line"].as_u64(), r[0]["value"].as_str()),
        (Some(1), Some("27"))
    );
    assert_eq!(
        (r[1]["line"].as_u64(), r[1]["value"].as_str()),
        (Some(3), Some("2875"))
    );
    assert_eq!(r[2]["status"], "invalid");
    assert_eq!(r[3]["status"], "unknown");
    assert_eq!(r[4]["status"], "invalid");
}

#[test]
fn single_and_flags() {
    let out = Command::new(env!("CARGO_BIN_EXE_growi"))
        .args(["--no-cache", "G=2,D=5 in P^4 thru H^3*psi^22"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap().trim(),
        "-41369/110075314176"
    );
    let out = Command::new(env!("CARGO_BIN_EXE_growi"))
        .args(["lines", "--cache"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let out = Command::new(env!("CARGO_BIN_EXE_growi"))
        .args(["--bogus", "lines"])
        .output()
        .unwrap();
    assert!(!out.status.success());
}
