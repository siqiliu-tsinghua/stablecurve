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
struct FABER;

#[pymethods]
impl FABER {
    /// FABER(cache=None)：`cache` 为 SQLite L2 路径；None → ~/.faber/cache.sqlite；
    /// "" → 不用持久缓存。注意 L2 是进程全局的：新建实例会替换（或卸下）
    /// 此前实例挂接的缓存。
    #[new]
    #[pyo3(signature = (cache=None))]
    fn new(cache: Option<&str>) -> PyResult<Self> {
        match cache {
            Some("") => faber::socle::detach_cache(),
            Some(p) => faber::socle::attach_cache(p)
                .map_err(|e| PyValueError::new_err(format!("无法打开缓存 {p}: {e}")))?,
            None => {
                let home = std::env::var_os("HOME")
                    .or_else(|| std::env::var_os("USERPROFILE"))
                    .ok_or_else(|| PyValueError::new_err("无法确定主目录"))?;
                let dir = std::path::PathBuf::from(home).join(".faber");
                std::fs::create_dir_all(&dir).map_err(|e| {
                    PyValueError::new_err(format!("无法创建 {}: {e}", dir.display()))
                })?;
                // 缺省缓存打不开时降级为纯内存（cache_attached 反映实际状态）
                if faber::socle::attach_cache(&dir.join("cache.sqlite").to_string_lossy()).is_err()
                {
                    faber::socle::detach_cache();
                }
            }
        }
        Ok(FABER)
    }

    /// L2 缓存当前是否挂接（进程全局状态）。
    #[getter]
    fn cache_attached(&self) -> bool {
        faber::socle::cache_attached()
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
        let v = py.allow_threads(|| faber::hodge::klp(g, n, &kappa, &psi, &lam, &ch));
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
        let v = py.allow_threads(|| faber::dr::int_dr(g, &a, r, &kappa, &psi));
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
