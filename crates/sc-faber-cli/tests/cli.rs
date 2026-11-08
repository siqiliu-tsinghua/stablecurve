//! faber 二进制的协议测试：批处理每行都是合法 JSON；USER_GUIDE 例子成立。

use serde_json::Value;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

struct CacheDirectory(PathBuf);

impl CacheDirectory {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "sc-faber-cli-{}-{timestamp}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for CacheDirectory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn batch(input: &str) -> Vec<Value> {
    let mut child = Command::new(env!("CARGO_BIN_EXE_faber"))
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
fn user_guide_examples() {
    let r = batch(concat!(
        r#"{"op":"klp","g":2,"n":2,"psi":[[1,2]],"lambda":[[1,3]]}"#,
        "\n",
        r#"{"op":"dr","g":2,"a":[1,1,-2],"r":2,"kappa":[[1,1]],"psi":[[1,3]]}"#,
        "\n",
        r#"{"op":"hyperell","g":2,"n":0,"m":0,"kappa":[[3,1]],"psi":[]}"#,
        "\n",
        r#"{"op":"biell","g":2,"n":0,"m":0,"kappa":[[2,1]],"psi":[]}"#,
        "\n",
        r#"{"op":"grr","g":2,"n":1,"t":1,"l":1,"dvec":[0],"kappa":[[1,1]],"psi":[[1,2]]}"#,
        "\n",
        r#"{"op":"drpoly","g":1,"a":[2,3,-5],"kappa":[[2,1]]}"#,
        "\n",
        r#"{"op":"drspin","g":1,"a":[1,1],"kappa":[[1,1]]}"#,
        "\n",
    ));
    let want = [
        "1/1440", "13/144", "1/1152", "1/32", "5/576", "37/24", "1/24",
    ];
    assert_eq!(r.len(), want.len());
    for (x, w) in r.iter().zip(want) {
        assert_eq!(x["status"], "ok", "{x}");
        assert_eq!(x["value"], w, "{x}");
    }
}

#[test]
fn errors_are_valid_json_and_helpful() {
    let r = batch(concat!(
        r#"{"op":"klp","g":"2","n":1}"#,
        "\n",
        r#"{"op":"dr","g":2,"A":[1,1,-2]}"#,
        "\n",
        r#"{"op":"dr","g":2,"a":[1,1,-2],"psi":[[4,1]]}"#,
        "\n",
        r#"{"op":"drpoly","g":1,"a":[]}"#,
        "\n",
        r#"{"op":"drpoly","g":1,"a":[1,2]}"#,
        "\n",
        r#"{"op":"drpoly","g":1,"a":[1,2],"k":1}"#,
        "\n",
        r#"{"op":"nope","g":1}"#,
        "\n",
        "not json \"quoted\" \\ text\n",
    ));
    assert_eq!(r.len(), 8);
    for x in &r {
        assert_eq!(x["status"], "error", "{x}");
    }
    // 大小写写错的字段名要被指出
    assert!(r[1]["error"].as_str().unwrap().contains('A'), "{}", r[1]);
    // 未知 op 的提示列出全部操作
    assert!(
        r[6]["error"].as_str().unwrap().contains("witten"),
        "{}",
        r[6]
    );
}

#[test]
fn single_query() {
    let out = Command::new(env!("CARGO_BIN_EXE_faber"))
        .arg(r#"{"op":"klp","g":2,"n":1,"psi":[[1,2],[1,2]]}"#)
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(String::from_utf8(out.stdout).unwrap().trim(), "1/1152");
}

#[test]
fn explicit_cache_open_failure_is_a_process_error() {
    let directory = CacheDirectory::new();
    let out = Command::new(env!("CARGO_BIN_EXE_faber"))
        .arg(r#"{"op":"klp","g":0,"n":3}"#)
        .arg("--cache")
        .arg(&directory.0)
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(2));
    assert!(out.stdout.is_empty());
    assert!(
        String::from_utf8(out.stderr)
            .unwrap()
            .contains("无法打开缓存")
    );
}

#[test]
fn successful_exit_has_committed_its_named_integral_root() {
    let directory = CacheDirectory::new();
    let path = directory.0.join("cache.sqlite");
    let out = Command::new(env!("CARGO_BIN_EXE_faber"))
        .arg(r#"{"op":"dr","g":1,"a":[1,-1],"r":1,"kappa":[[1,1]]}"#)
        .arg("--cache")
        .arg(&path)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(String::from_utf8(out.stdout).unwrap().trim(), "1/24");

    let session = sc_faber::socle::CacheSession::open(path.to_str().unwrap()).unwrap();
    // DR-v1: g, degree; i64-slice length and complete signed weights;
    // then length-prefixed kappa and psi pairs. This lookup does no mathematics.
    let key = [1, 1, 2, 0, 1, 0, u32::MAX, u32::MAX, 1, 0, 1, 1, 0, 0];
    let persisted = session
        .run(|| Ok(sc_faber::socle::cache_get("scalar-dr-v1", &key)))
        .unwrap();
    assert_eq!(persisted, Some(sc_faber::Q::from_signeds(1, 24)));
    session.close().unwrap();
}
