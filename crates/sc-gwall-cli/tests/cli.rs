//! gwall 二进制的协议测试：批处理每行都是合法 JSON，两种插入写法等价。

use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

fn batch(input: &str) -> Vec<Value> {
    let mut child = Command::new(env!("CARGO_BIN_EXE_gwall"))
        .arg("--batch")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
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
fn array_and_object_inserts_agree() {
    let r = batch(concat!(
        r#"{"op":"gw","g":2,"d":1,"insert":[["t1",0],["pt",1],["pt",0]]}"#,
        "\n",
        r#"{"op":"gw","g":2,"d":1,"insert":[{"t":"t1","k":0},{"t":"pt","k":1},{"t":"pt","k":0}]}"#,
        "\n",
        r#"{"op":"gw","g":2,"d":1,"insert":[["chL",1],["pt",0],["pt",0],["pt",0],["ldeg",2]]}"#,
        "\n",
        r#"{"op":"gw","g":2,"d":1,"insert":[["chL",1],["pt",0],["pt",0],["pt",0],{"t":"ldeg","v":2}]}"#,
        "\n",
        r#"{"op":"gw","g":2,"d":1,"insert":[["chL",1],["pt",0],["pt",0],["pt",0]],"ldeg":2}"#,
        "\n",
    ));
    assert_eq!(r.len(), 5);
    for x in &r {
        assert_eq!(x["status"], "ok", "{x}");
    }
    assert_eq!(r[0]["value"], r[1]["value"]);
    assert_eq!(r[2]["value"], r[3]["value"]);
    assert_eq!(r[2]["value"], r[4]["value"]);
}

#[test]
fn values_and_series() {
    let r = batch(concat!(
        r#"{"op":"bps","g":3,"d":2,"h":1}"#,
        "\n",
        "\n",
        r#"{"op":"l","g":2,"d":2,"n":2}"#,
        "\n",
    ));
    assert_eq!(r.len(), 2);
    assert_eq!(r[0]["line"], 1);
    assert_eq!(r[0]["value"], "4");
    assert_eq!(r[1]["line"], 3);
    assert_eq!(r[1]["value"], serde_json::json!(["15/2", "0", "-5/4"]));
}

#[test]
fn errors_are_valid_json() {
    let r = batch(concat!(
        r#"{"op":"bps","g":"3","d":2,"h":1}"#,
        "\n",
        r#"{"op":"bps","g":3,"h":1}"#,
        "\n",
        r#"{"op":"gw","g":2,"d":1,"insert":[["zz",1]]}"#,
        "\n",
        "not json \"quoted\" \\ text\n",
    ));
    assert_eq!(r.len(), 4);
    for x in &r {
        assert_eq!(x["status"], "error", "{x}");
        assert!(x["error"].is_string());
    }
}

#[test]
fn single_query() {
    let out = Command::new(env!("CARGO_BIN_EXE_gwall"))
        .arg(r#"{"op":"l","g":2,"d":2,"n":2}"#)
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap().trim(),
        "[15/2,0,-5/4]"
    );
}
