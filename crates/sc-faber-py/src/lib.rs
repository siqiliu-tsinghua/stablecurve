//! # sc_faber
//!
//! sc-faber 计算引擎的 Python 原生模块（pyo3）。
//!
//! ```python
//! from sc_faber import FABER
//!
//! f = FABER()                                  # ~/.faber/cache.sqlite（L2）
//! f.klp(2, 2, psi=[(1,2)], lam=[(1,3)])        # → Fraction(1, 1440)
//! f.dr(2, [1,1,-2], r=2, kappa=[(1,1)], psi=[(1,3)])  # → Fraction(13, 144)
//! ```
//!
//! 参数：kappa/psi/lam/ch 均为 (指标, 幂) 对的列表。计算释放 GIL。

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn to_fraction(py: Python<'_>, v: &faber::Q) -> PyObject {
    let s = v.to_string();
    PyModule::import(py, "fractions")
        .and_then(|m| m.getattr("Fraction"))
        .and_then(|f| f.call1((s,)))
        .expect("构造 Fraction 失败")
        .into()
}

// 类名即 Python 接口名（沿用原程序名），保留全大写
#[allow(clippy::upper_case_acronyms)]
#[pyclass(unsendable)]
struct FABER {
    session: Option<faber::socle::CacheSession>,
}

fn with_session<T>(
    session: Option<&faber::socle::CacheSession>,
    calculate: impl FnOnce() -> Result<T, String>,
) -> Result<T, String> {
    match session {
        Some(session) => session.run(calculate),
        None => faber::socle::CacheSession::run_uncached(calculate),
    }
}

#[pymethods]
impl FABER {
    /// Set the total-process memory budget shared by cache-backed calculations.
    /// None selects the automatic budget derived from available system memory.
    #[staticmethod]
    #[pyo3(signature = (maximum_bytes=None))]
    fn set_memory_limit(maximum_bytes: Option<u64>) -> PyResult<()> {
        faber::socle::set_cache_memory_limit(maximum_bytes).map_err(PyValueError::new_err)
    }

    /// FABER(cache=None)：`cache` 为 SQLite L2 路径；None → ~/.faber/cache.sqlite；
    /// "" → 不用持久缓存。每个实例独立持有自己的缓存会话。
    #[new]
    #[pyo3(signature = (cache=None))]
    fn new(cache: Option<&str>) -> PyResult<Self> {
        let session = match cache {
            Some("") => None,
            Some(p) => Some(
                faber::socle::CacheSession::open(p)
                    .map_err(|e| PyValueError::new_err(format!("无法打开缓存 {p}: {e}")))?,
            ),
            None => {
                let home = std::env::var_os("HOME")
                    .or_else(|| std::env::var_os("USERPROFILE"))
                    .ok_or_else(|| PyValueError::new_err("无法确定主目录"))?;
                let dir = std::path::PathBuf::from(home).join(".faber");
                std::fs::create_dir_all(&dir).map_err(|e| {
                    PyValueError::new_err(format!("无法创建 {}: {e}", dir.display()))
                })?;
                let path = dir.join("cache.sqlite");
                Some(
                    faber::socle::CacheSession::open(&path.to_string_lossy()).map_err(|e| {
                        PyValueError::new_err(format!("无法打开缓存 {}: {e}", path.display()))
                    })?,
                )
            }
        };
        Ok(FABER { session })
    }

    /// 当前实例是否有可用的持久缓存。
    #[getter]
    fn cache_attached(&self) -> bool {
        self.session
            .as_ref()
            .is_some_and(|session| session.is_open())
    }

    #[getter]
    fn cache_error(&self) -> Option<String> {
        self.session.as_ref().and_then(|session| session.error())
    }

    fn flush_cache(&self, py: Python<'_>) -> PyResult<()> {
        let session = self.session.as_ref();
        py.allow_threads(|| session.map_or(Ok(()), |session| session.flush()))
            .map_err(PyValueError::new_err)
    }

    fn close_cache(&mut self, py: Python<'_>) -> PyResult<()> {
        let session = self.session.as_ref();
        py.allow_threads(|| session.map_or(Ok(()), |session| session.close()))
            .map_err(PyValueError::new_err)?;
        // A failed close retains its owner and its latched error. Later calls
        // must not silently change from persistence to pure-memory operation.
        self.session.take();
        Ok(())
    }

    /// 清空内存中的稳定图缓存（dr 等按 (g, n) 逐层缓存枚举出的稳定图，
    /// 大的 (g, n, 余维) 可达数十万个图、数百 MB；只影响当前线程）。
    fn clear_graph_cache(&self) {
        faber::generators::clear_graph_cache();
    }

    /// ∫ κ^kappa · ψ^psi · λ^lam · ch^ch（M̄_{g,n} 基本类上单项式）→ Fraction。
    #[pyo3(signature = (g, n, kappa=Vec::new(), psi=Vec::new(), lam=Vec::new(), ch=Vec::new()))]
    fn klp(
        &mut self,
        py: Python<'_>,
        g: u32,
        n: u32,
        kappa: Vec<(u32, u32)>,
        psi: Vec<(u32, u32)>,
        lam: Vec<(u32, u32)>,
        ch: Vec<(u32, u32)>,
    ) -> PyResult<PyObject> {
        let session = self.session.as_ref();
        let v = py.allow_threads(|| {
            with_session(session, || faber::hodge::klp(g, n, &kappa, &psi, &lam, &ch))
        });
        match v {
            Ok(v) => Ok(to_fraction(py, &v)),
            Err(e) => Err(PyValueError::new_err(e)),
        }
    }

    /// ∫ DR_g(A) · κ^kappa · ψ^psi → Fraction（r 省略 = g）。
    #[pyo3(signature = (g, a, r=None, kappa=Vec::new(), psi=Vec::new()))]
    fn dr(
        &mut self,
        py: Python<'_>,
        g: u32,
        a: Vec<i64>,
        r: Option<u32>,
        kappa: Vec<(u32, u32)>,
        psi: Vec<(u32, u32)>,
    ) -> PyResult<PyObject> {
        let session = self.session.as_ref();
        let v = py
            .allow_threads(|| with_session(session, || faber::dr::int_dr(g, &a, r, &kappa, &psi)));
        match v {
            Ok(v) => Ok(to_fraction(py, &v)),
            Err(e) => Err(PyValueError::new_err(e)),
        }
    }
}

#[pymodule]
fn _sc_faber(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FABER>()?;
    Ok(())
}
