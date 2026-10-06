//! sc-gwall 的 pyo3 绑定。
//!
//! ```python
//! from sc_gwall import GWALL
//! w = GWALL()
//! w.gw(2, 1, [("t1", 0), ("pt", 0), ("pt", 0)])   # Fraction
//! w.bps(3, 2, 1)                                   # Fraction(4) —— 论文表 n_2^1(3)
//! ```
//!
//! 插入项可用元组或列表；a/b 环路可按任意顺序给出（符号自动计入）。
//! 计算期间释放 GIL。

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PySequence};
use sc_gwall::gwl::{Ctx, Ins, gw_ctx};
use sc_gwall::pf::{bps_local, l_series, z, ztil};

// 类名即 Python 接口名（沿用原程序名），保留全大写
#[allow(clippy::upper_case_acronyms)]
#[pyclass]
struct GWALL {
    ctx: Ctx,
}

fn q_to_py(py: Python<'_>, q: &sc_num::Q) -> PyResult<PyObject> {
    let frac = PyModule::import(py, "fractions")?.getattr("Fraction")?;
    Ok(frac.call1((q.to_string(),))?.unbind())
}

fn series_to_py(py: Python<'_>, c: &[sc_num::Q], n: u32) -> PyResult<PyObject> {
    let list = PyList::empty(py);
    let zero = sc_num::Q::from(0);
    for i in 0..=n as usize {
        list.append(q_to_py(py, c.get(i).unwrap_or(&zero))?)?;
    }
    Ok(list.into_any().unbind())
}

fn ins_from(items: &Bound<'_, PyAny>) -> PyResult<Vec<Ins>> {
    let mut out = Vec::new();
    for it in items.try_iter()? {
        let it = it?;
        let seq = it.downcast::<PySequence>().map_err(|_| {
            PyValueError::new_err(format!("插入项 {it} 须为元组或列表，如 (\"pt\", 0)"))
        })?;
        let t: String = seq.get_item(0)?.extract()?;
        let want = if t == "a" || t == "b" { 3 } else { 2 };
        if seq.len()? != want {
            return Err(PyValueError::new_err(format!(
                "插入项 {it} 应有 {want} 个元素"
            )));
        }
        let ins = match t.as_str() {
            "t1" => Ins::T1(seq.get_item(1)?.extract()?),
            "pt" => Ins::Pt(seq.get_item(1)?.extract()?),
            "a" => Ins::A(seq.get_item(1)?.extract()?, seq.get_item(2)?.extract()?),
            "b" => Ins::B(seq.get_item(1)?.extract()?, seq.get_item(2)?.extract()?),
            "chO" => Ins::ChO(seq.get_item(1)?.extract()?),
            "chL" => Ins::ChL(seq.get_item(1)?.extract()?),
            "ldeg" => Ins::LDeg(seq.get_item(1)?.extract()?),
            other => {
                return Err(PyValueError::new_err(format!(
                    "未知插入类型 '{other}'（可用：t1/pt/a/b/chO/chL/ldeg）"
                )));
            }
        };
        out.push(ins);
    }
    Ok(out)
}

#[pymethods]
impl GWALL {
    #[new]
    fn new() -> Self {
        GWALL { ctx: Ctx::new() }
    }

    /// gw(g, d, insert) —— 不连通 local curve GW 不变量。
    /// insert 形如 [("t1", 0), ("pt", 1), ("a", 2, 1), ("chL", 3), ("ldeg", 1)]。
    #[pyo3(signature = (g, d, insert))]
    fn gw(
        &mut self,
        py: Python<'_>,
        g: u32,
        d: u32,
        insert: &Bound<'_, PyAny>,
    ) -> PyResult<PyObject> {
        let ins = ins_from(insert)?;
        let ctx = &mut self.ctx;
        let v = py.allow_threads(|| gw_ctx(ctx, g, d, &ins));
        q_to_py(py, &v)
    }

    /// bps(g, d, h) —— local BPS n_d^h(g)（全量 GW 反演）。
    #[pyo3(signature = (g, d, h))]
    fn bps(&mut self, py: Python<'_>, g: u32, d: u32, h: u32) -> PyResult<PyObject> {
        let ctx = &mut self.ctx;
        let v = py
            .allow_threads(|| bps_local(ctx, g, d, h))
            .map_err(PyValueError::new_err)?;
        q_to_py(py, &v)
    }

    /// ztil(g, d, n) —— t 级数系数列表（含常映射分量）。
    #[pyo3(signature = (g, d, n))]
    fn ztil(&mut self, py: Python<'_>, g: u32, d: u32, n: u32) -> PyResult<PyObject> {
        let ctx = &mut self.ctx;
        let s = py
            .allow_threads(|| ztil(ctx, g, d, n))
            .map_err(PyValueError::new_err)?;
        series_to_py(py, &s.c, n)
    }

    /// z(g, d, n) —— 去常映射（exp(−N₀)）级数。
    #[pyo3(signature = (g, d, n))]
    fn z(&mut self, py: Python<'_>, g: u32, d: u32, n: u32) -> PyResult<PyObject> {
        let ctx = &mut self.ctx;
        let s = py
            .allow_threads(|| z(ctx, g, d, n))
            .map_err(PyValueError::new_err)?;
        series_to_py(py, &s.c, n)
    }

    /// l(g, d, n) —— 度变量 log 连通化 L_d 级数（论文表 N 的母函数）。
    #[pyo3(signature = (g, d, n))]
    fn l(&mut self, py: Python<'_>, g: u32, d: u32, n: u32) -> PyResult<PyObject> {
        let ctx = &mut self.ctx;
        let s = py
            .allow_threads(|| l_series(ctx, g, d, n))
            .map_err(PyValueError::new_err)?;
        series_to_py(py, &s.c, n)
    }
}

#[pymodule]
fn _sc_gwall(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GWALL>()?;
    Ok(())
}
