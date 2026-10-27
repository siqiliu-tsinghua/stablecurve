//! 持久缓存写入失败须可见：卸下 L2（`cache_attached` 变为 false），不静默
//! 吞掉错误。L2 是进程全局的，故单独成一个测试二进制。

use sc_faber::socle::{attach_cache, cache_attached, cache_put};

#[test]
fn failed_write_detaches_cache() {
    let dir = std::env::temp_dir().join(format!("sc-faber-cache-fail-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("cache.sqlite");
    attach_cache(path.to_str().unwrap()).unwrap();
    assert!(cache_attached());
    // 另一个连接删表，模拟库被破坏
    rusqlite::Connection::open(&path)
        .unwrap()
        .execute("DROP TABLE st", [])
        .unwrap();
    cache_put("st", &[1, 2, 3], &sc_faber::Q::from(1));
    assert!(!cache_attached());
    let _ = std::fs::remove_dir_all(&dir);
}
