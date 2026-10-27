//! 两级缓存：L1 内存 + L2 SQLite（REDESIGN §2.5）。
//!
//! 语义对齐原版 GW.dbm——"算过即落盘、跨会话复用中间不变量"；机制换代：
//! 读穿（L1 → L2 → 计算）+ 写回（dirty 队列批量事务），取代原版逐条同步
//! `gdbm_store`。
//!
//! 键 = 规范形 [`Inv`] 的版本化变长整数编码（LEB128/zigzag），不再用原版的
//! 私有位打包。值 = 分子/分母的十进制字符串（精确、跨平台）。

use std::collections::HashMap;
use std::path::Path;

use rusqlite::Connection;

use crate::inv::Inv;
use crate::num::Q;

/// 键格式版本；编码语义变化时 +1，旧库自动作废。
const KEY_FORMAT_VERSION: u8 = 1;
/// 算法版本；策略链行为变化时更新，防止旧值污染新算法。
const ALGO_VERSION: &str = "sc-growi 0.1.0";
/// dirty 队列达到此长度即落盘。
const DIRTY_LIMIT: usize = 4096;

/// 不变量缓存（L1 + 可选 L2）。
pub struct Cache {
    l1: HashMap<Inv, Q>,
    dirty: Vec<(Inv, Q)>,
    db: Option<Db>,
}

struct Db {
    conn: Connection,
}

impl Cache {
    /// 纯内存缓存（测试用）。
    pub fn memory() -> Self {
        Cache {
            l1: HashMap::new(),
            dirty: Vec::new(),
            db: None,
        }
    }

    /// 带磁盘 L2 的缓存；文件不存在则创建。
    pub fn with_db(path: &Path) -> rusqlite::Result<Self> {
        let mut cache = Cache::memory();
        cache.db = Some(Db::open(path)?);
        Ok(cache)
    }

    /// 读穿查询：L1 → L2 →（都没有）None。L2 命中时提升进 L1。
    pub fn get(&mut self, inv: &Inv) -> Option<Q> {
        if let Some(v) = self.l1.get(inv) {
            return Some(v.clone());
        }
        let db = self.db.as_ref()?;
        let key = encode_key(inv);
        let value = db.lookup(&key)?;
        self.l1.insert(inv.clone(), value.clone());
        Some(value)
    }

    /// 写 L1 并记入 dirty；达到阈值自动落盘。
    pub fn put(&mut self, inv: Inv, value: Q) {
        if crate::env_flag("SC_GROWI_KEYS") {
            eprintln!("[key] {inv} = {value}");
        }
        if self.l1.insert(inv.clone(), value.clone()).is_none() {
            self.dirty.push((inv, value));
        }
        if self.dirty.len() >= DIRTY_LIMIT {
            self.flush();
        }
    }

    /// 批量落盘 dirty 队列（单事务）。
    ///
    /// 写回失败（busy_timeout 之后仍失败，多为磁盘满、只读或库损坏）时打印
    /// 警告并卸下 L2：本会话的结果仍在 L1，只是不再持久化。
    pub fn flush(&mut self) {
        if self.dirty.is_empty() {
            return;
        }
        if let Some(db) = self.db.as_mut() {
            if let Err(e) = db.store_batch(&self.dirty) {
                eprintln!(
                    "sc-growi: 缓存写回失败（{e}），{} 条结果未保存；本会话停用持久缓存",
                    self.dirty.len()
                );
                self.db = None;
            }
        }
        self.dirty.clear();
    }

    pub fn len(&self) -> usize {
        self.l1.len()
    }

    pub fn is_empty(&self) -> bool {
        self.l1.is_empty()
    }
}

impl Drop for Cache {
    fn drop(&mut self) {
        self.flush();
    }
}

impl Db {
    fn open(path: &Path) -> rusqlite::Result<Self> {
        if let Some(dir) = path.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let conn = Connection::open(path)?;
        // WAL + NORMAL：缓存数据可再生，免逐写 fsync，换批量吞吐（§2.5）
        conn.pragma_update(None, "journal_mode", "WAL")?;
        conn.pragma_update(None, "synchronous", "NORMAL")?;
        conn.pragma_update(None, "busy_timeout", 5000)?;
        conn.execute_batch(
            "CREATE TABLE IF NOT EXISTS meta(
                k TEXT PRIMARY KEY, v TEXT NOT NULL);
             CREATE TABLE IF NOT EXISTS invariants(
                key BLOB PRIMARY KEY, val TEXT NOT NULL)
                WITHOUT ROWID;",
        )?;
        let mut db = Db { conn };
        db.check_version()?;
        Ok(db)
    }

    /// 版本校验：格式或算法版本不符 → 整库作废重建（取代原版 CHECK_TIME）。
    fn check_version(&mut self) -> rusqlite::Result<()> {
        let read = |db: &Db, k: &str| -> Option<String> {
            db.conn
                .query_row("SELECT v FROM meta WHERE k=?1", [k], |r| r.get(0))
                .ok()
        };
        let want = |v: &str| v.to_string();
        if read(self, "key_format") != Some(want(&KEY_FORMAT_VERSION.to_string()))
            || read(self, "algo") != Some(want(ALGO_VERSION))
        {
            self.conn.execute("DELETE FROM invariants", [])?;
            let mut stmt = self
                .conn
                .prepare("INSERT OR REPLACE INTO meta(k,v) VALUES(?1,?2)")?;
            stmt.execute(("key_format", KEY_FORMAT_VERSION.to_string()))?;
            stmt.execute(("algo", ALGO_VERSION))?;
        }
        Ok(())
    }

    fn lookup(&self, key: &[u8]) -> Option<Q> {
        // 值 = 有理数的规范文本（"a/b" 或 "n"，含符号），FromStr 精确回读
        let s: String = self
            .conn
            .query_row("SELECT val FROM invariants WHERE key=?1", [key], |r| {
                r.get(0)
            })
            .ok()?;
        s.parse().ok()
    }

    fn store_batch(&mut self, entries: &[(Inv, Q)]) -> rusqlite::Result<()> {
        let tx = self.conn.transaction()?;
        {
            let mut stmt =
                tx.prepare_cached("INSERT OR REPLACE INTO invariants(key,val) VALUES(?1,?2)")?;
            for (inv, q) in entries {
                stmt.execute((encode_key(inv).as_slice(), q.to_string()))?;
            }
        }
        tx.commit()
    }
}

/// 规范形不变量 → 版本化 LEB128 键。
pub fn encode_key(inv: &Inv) -> Vec<u8> {
    let mut out = vec![KEY_FORMAT_VERSION];
    write_zigzag(&mut out, inv.dim as i64);
    write_zigzag(&mut out, inv.deg as i64);
    write_zigzag(&mut out, inv.genus as i64);
    write_zigzag(&mut out, inv.cdeg as i64);
    write_zigzag(&mut out, inv.pts.len() as i64);
    for p in inv.pts.iter() {
        write_zigzag(&mut out, p.mult as i64);
        write_zigzag(&mut out, p.psi as i64);
        write_zigzag(&mut out, p.h as i64);
    }
    out
}

fn write_zigzag(out: &mut Vec<u8>, v: i64) {
    let mut u = ((v << 1) ^ (v >> 63)) as u64;
    loop {
        let byte = (u & 0x7f) as u8;
        u >>= 7;
        if u == 0 {
            out.push(byte);
            break;
        }
        out.push(byte | 0x80);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inv::Pdata;

    fn lines_p2() -> Inv {
        Inv::new(2, 1, 0, 1, [Pdata::new(2, 0, 0); 2])
    }

    #[test]
    fn memory_roundtrip() {
        let mut c = Cache::memory();
        assert!(c.get(&lines_p2()).is_none());
        c.put(lines_p2(), Q::from(1));
        assert_eq!(c.get(&lines_p2()), Some(Q::from(1)));
    }

    #[test]
    fn sqlite_roundtrip() {
        let dir = std::env::temp_dir().join("sc-growi-test-cache");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.sqlite");
        {
            let mut c = Cache::with_db(&path).unwrap();
            assert!(c.get(&lines_p2()).is_none());
            c.put(lines_p2(), crate::num::q_frac(-1, 12));
            c.flush();
        }
        // 新进程（新 Cache 实例）应从 L2 读到
        let mut c2 = Cache::with_db(&path).unwrap();
        assert_eq!(c2.get(&lines_p2()), Some(crate::num::q_frac(-1, 12)));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// 写回失败不静默：卸下 L2（可见），结果仍留在 L1。
    #[test]
    fn failed_flush_detaches_l2() {
        let dir = std::env::temp_dir().join("sc-growi-test-cache-fail");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let mut c = Cache::with_db(&dir.join("cache.sqlite")).unwrap();
        let db = c.db.as_ref().unwrap();
        db.conn.execute("DROP TABLE invariants", []).unwrap();
        c.put(lines_p2(), Q::from(1));
        c.flush();
        assert!(c.db.is_none());
        assert_eq!(c.get(&lines_p2()), Some(Q::from(1)));
        drop(c);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn key_encoding_deterministic() {
        let k1 = encode_key(&lines_p2());
        let k2 = encode_key(&lines_p2().into_canonical());
        assert_eq!(k1, k2);
        assert_ne!(
            k1,
            encode_key(&Inv::new(2, 1, 0, 2, [Pdata::new(2, 0, 0); 2]))
        );
    }
}
