//! # sc_growi
//!
//! sc-growi 计算引擎的 Python 原生模块（pyo3）。
//!
//! ```python
//! from sc_growi import GROWI
//!
//! g = GROWI()                                   # ~/.growi/cache.sqlite（两级缓存）
//! g = GROWI(cache="")                           # 只用内存
//! g.gw("lines in quintic threefold")            # → Fraction(2875)
//! g.gw("G=2,D=5 in P^4 thru H^3*psi^22")        # → Fraction(-41369, 110075314176)
//! ```
//!
//! 语义映射：算出 → `fractions.Fraction`（精确值，大数以字符串过界）；
//! 算不出 → `None`（原版 UNDEFINED）；输入非法 → `ValueError`。
//! `compute` 是长计算——pyo3 侧已释放 GIL，计算期间 Python 其他线程照常运行。

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use sc_growi::gw;

fn to_fraction(py: Python<'_>, v: &sc_growi::Q) -> PyObject {
    let s = v.to_string(); // "num" 或 "num/den"，Fraction 的字符串构造子精确解析
    PyModule::import(py, "fractions")
        .and_then(|m| m.getattr("Fraction"))
        .and_then(|f| f.call1((s,)))
        .expect("构造 Fraction 失败")
        .into()
}

// 类名即 Python 接口名（沿用原程序名），保留全大写
#[allow(clippy::upper_case_acronyms)]
#[pyclass(unsendable)]
struct GROWI {
    ctx: sc_growi::Ctx,
}

#[pymethods]
impl GROWI {
    /// GROWI(cache=None)：`cache` 为 SQLite 路径；None → ~/.growi/cache.sqlite；
    /// "" → 只用内存（不落盘）。
    #[new]
    #[pyo3(signature = (cache=None))]
    fn new(cache: Option<&str>) -> PyResult<Self> {
        let ctx = match cache {
            Some("") => sc_growi::Ctx::memory(),
            Some(p) => sc_growi::Ctx::with_db(std::path::Path::new(p))
                .map_err(|e| PyValueError::new_err(format!("无法打开缓存 {p}: {e}")))?,
            None => {
                let home = std::env::var_os("HOME")
                    .or_else(|| std::env::var_os("USERPROFILE"))
                    .ok_or_else(|| PyValueError::new_err("无法确定主目录"))?;
                let dir = std::path::PathBuf::from(home).join(".growi");
                std::fs::create_dir_all(&dir).map_err(|e| {
                    PyValueError::new_err(format!("无法创建 {}: {e}", dir.display()))
                })?;
                sc_growi::Ctx::with_db(&dir.join("cache.sqlite"))
                    .map_err(|e| PyValueError::new_err(format!("无法打开默认缓存: {e}")))?
            }
        };
        Ok(GROWI { ctx })
    }

    /// gw(expr) → Fraction | None。计算期间释放 GIL。
    fn gw(&mut self, py: Python<'_>, expr: &str) -> PyResult<Option<PyObject>> {
        let query = sc_growi::query::parse(expr).map_err(PyValueError::new_err)?;
        let r = py.allow_threads(|| gw(&mut self.ctx, &query.settings, &query.classes));
        match r {
            Ok(v) => Ok(Some(to_fraction(py, &v))),
            Err(sc_growi::GwError::Unknown(_)) => Ok(None),
            Err(e) => Err(PyValueError::new_err(e.to_string())),
        }
    }

    /// 已缓存的中间不变量个数（L1）。
    #[getter]
    fn cache_size(&self) -> usize {
        self.ctx.cache.len()
    }

    /// 把 dirty 队列立即落盘（析构时也会自动 flush）。
    fn flush(&mut self) {
        self.ctx.flush();
    }
}

/// Python 侧的类表达式解析（复用 CLI 的解析器，便于程序化构造查询）。
#[pyfunction]
fn parse(expr: &str) -> PyResult<String> {
    let q = sc_growi::query::parse(expr).map_err(PyValueError::new_err)?;
    Ok(format!(
        "N={} E={} G={} D={} classes={}",
        q.settings.dim,
        q.settings.deg,
        q.settings.genus,
        q.settings.cdeg,
        q.classes
            .iter()
            .map(|c| format!("{}(mult={},count={})", c.poly, c.mult, c.count))
            .collect::<Vec<_>>()
            .join(",")
    ))
}

#[pymodule]
fn _sc_growi(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GROWI>()?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    Ok(())
}
