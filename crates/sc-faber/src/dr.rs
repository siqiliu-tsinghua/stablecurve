//! 双重分支循环 DR_{g}(A)（Pixton 公式，插值形式）。
//!
//! 移植自 admcycles `DR/double_ramification_cycle.py` 的 `DR_compute/DR_coeff_new`
//! （kval=0 的普通 DR）。公式：对每个余维 r 的生成元 Γ（稳定图 + ψ 装饰；
//! 含 κ 的生成元系数为 0），
//!
//! coeff(Γ) = (1/|Aut Γ|) · Π_{标记 m} d_m^{2 e_m}/e_m! · Π_{结点} (−1)/(e₁!e₂!(e₁+e₂+1))
//!            · Q_Γ(0)
//!
//! 其中 e_m = ψ 在标记 m 处的幂、e₁,e₂ = ψ 在结点两半边的幂，
//! Q_Γ(m) = m^{−b₁(Γ)} Σ_{平衡权重 w: E→{0..m−1}} Π_e w_e^{2(e₁+e₂+1)}，
//! 平衡 = 每顶点 (Σ_{v 的标记} d + Σ_{出边} w − Σ_{入边} w) ≡ 0 (mod m)；
//! Q 为次数 ≤ 2Σexp+1 的多项式（Pixton 多项式性定理），在
//! m ∈ {m₀+1..m₀+deg+2}（m₀ = ⌈Σ|d|/2⌉）取样后拉格朗日插值到 0。
//! 总体再乘 2^{−r}。
//!
//! 手验锚点（admcycles 文档值 × 自同构调整后）：DR(1,(2,3,−5)) =
//! 2ψ₁+9/2ψ₂+25/2ψ₃ −(1/12)[Δ₀] −(25/2)[3|12] −…（我们的 [Δ₀] = 物理轨迹类

use crate::kppoly::KappaPsiPolynomial;
use crate::taut::{DecStratum, TautClass};
use sc_num::{One, Q, Zero, is_zero};
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

fn fact(n: u64) -> Q {
    let mut r = Q::ONE;
    for i in 2..=n {
        r *= Q::from(i as i64);
    }
    r
}

/// 边因子表 [w^{2s}]_{w=0..m−1}（按 (s, m) 跨生成元缓存）。
fn weight_table(s: u32, m: i64) -> Rc<Vec<Q>> {
    thread_local! {
        static MEMO: RefCell<HashMap<(u32, i64), Rc<Vec<Q>>>> = RefCell::new(HashMap::new());
    }
    if let Some(t) = MEMO.with(|c| c.borrow().get(&(s, m)).cloned()) {
        return t;
    }
    let t: Rc<Vec<Q>> = Rc::new(
        (0..m)
            .map(|w| {
                let mut t = Q::ONE;
                for _ in 0..2 * s {
                    t *= Q::from(w);
                }
                t
            })
            .collect(),
    );
    MEMO.with(|c| c.borrow_mut().insert((s, m), t.clone()));
    t
}

/// 平衡权重系统的和：Σ Π_e w_e^{2 exp_e}（mod m 平衡，u 出边 +w、v 入边 −w）。
fn balanced_weight_sum(edges: &[(usize, usize)], exps: &[u32], given: &[i64], m: i64) -> Q {
    let tables: Vec<Rc<Vec<Q>>> = exps.iter().map(|&s| weight_table(s, m)).collect();
    let table: Vec<&[Q]> = tables.iter().map(|t| t.as_slice()).collect();
    crate::flow::balanced_flow_sum(edges, given, m, &table)
}

/// 拉格朗日插值：给定点 (xs, ys)，求 p(0)。
fn interpolate_at_zero(xs: &[i64], ys: &[Q]) -> Q {
    let mut acc = Q::ZERO;
    for (j, (&xj, yj)) in xs.iter().zip(ys).enumerate() {
        let mut term = yj.clone();
        for (k, &xk) in xs.iter().enumerate() {
            if k != j {
                // (0 − xk)/(xj − xk)
                term *= Q::from(-xk) / Q::from(xj - xk);
            }
        }
        acc += term;
    }
    acc
}

/// DR_g(A) 余维 r 的类（A 长度 = n，ΣA = 0；r 默认 g）。
pub fn dr_cycle(g: u32, dvector: &[i64], r: Option<u32>) -> Result<TautClass, String> {
    dr_cycle_impl(g, dvector, r, false)
}

/// 生成元口径（含 1/|Aut|，与 admcycles 的 DR_compute 逐系数可比）；
/// 消没定理（r>g ⇒ 0）在此口径下逐项成立。
pub fn dr_cycle_gen(g: u32, dvector: &[i64], r: Option<u32>) -> Result<TautClass, String> {
    dr_cycle_impl(g, dvector, r, true)
}

fn dr_cycle_impl(
    g: u32,
    dvector: &[i64],
    r: Option<u32>,
    gen_norm: bool,
) -> Result<TautClass, String> {
    let n = dvector.len() as u32;
    let r = r.unwrap_or(g);
    if dvector.iter().sum::<i64>() != 0 {
        return Err(format!("DR 要求 ΣA = 0，得到 {dvector:?}"));
    }
    if 2 * g as i64 - 2 + n as i64 <= 0 {
        return Err(format!("M̄_({g},{n}) 不稳定"));
    }
    if r as i64 > 3 * g as i64 - 3 + n as i64 {
        return Err("余维超过模空间维数".into());
    }
    let m0 = (dvector.iter().map(|d| d.abs()).sum::<i64>() + 1) / 2;
    let mut terms: Vec<DecStratum> = Vec::new();
    // 含 κ 的生成元系数为 0（kval = 0），只枚举纯 ψ 装饰（顺序与 tautgens
    // 中 κ 为空的部分一致）。顶点、给定值、圈秩逐图只算一次；Q_Γ(0) 只依赖
    // 图与各结点的指数，图内按结点指数缓存。
    for c in 0..=r {
        for graph in crate::generators::layer(g, n, c).iter() {
            let verts = graph.vertices();
            let nv = verts.len();
            let max_leg = verts
                .iter()
                .flat_map(|v| v.legs.iter().copied())
                .max()
                .unwrap_or(0) as usize;
            let mut owner = vec![0usize; max_leg + 1];
            for (i, v) in verts.iter().enumerate() {
                for &l in &v.legs {
                    owner[l as usize] = i;
                }
            }
            let edges_v: Vec<(usize, usize)> = graph
                .edges()
                .iter()
                .map(|&(a, b)| (owner[a as usize], owner[b as usize]))
                .collect();
            // given[v] = Σ_{v 的标记} d
            let mut given = vec![0i64; nv];
            for l in 1..=n {
                given[owner[l as usize]] += dvector[(l - 1) as usize];
            }
            let b1 = (graph.edges().len() + 1) as i64 - nv as i64;
            // 注意：admcycles 的 DR 系数含 1/autom_count（其生成元口径 ξ_*1）；
            // 我们的 TautClass 口径为 ξ_*/|Aut|（evaluate 时统一除），故此处
            // 不再放 1/aut，积分结果与 admcycles 直接可比。
            let aut_scalar = if gen_norm {
                Q::from_signeds(1, graph.automorphism_number() as i64)
            } else {
                Q::ONE
            };
            let mut q0_cache: HashMap<Vec<u32>, Q> = HashMap::new();
            for poly in crate::generators::psi_monomials(graph, r - c) {
                let t = &poly.terms[0];
                // 标记部分：d^{2e}/e!
                let mut scalar = aut_scalar.clone();
                for &(l, e) in &t.psi {
                    if l <= n {
                        let d = dvector[(l - 1) as usize];
                        let mut p = Q::ONE;
                        for _ in 0..2 * e {
                            p *= Q::from(d);
                        }
                        scalar = scalar * p / fact(e as u64);
                    }
                }
                // 结点部分：exp = e1+e2+1，因子 −1/(e1!e2!(e1+e2+1))
                let mut exps: Vec<u32> = Vec::with_capacity(edges_v.len());
                for &(a, b) in graph.edges() {
                    let e1 = t.psi_power_at(a);
                    let e2 = t.psi_power_at(b);
                    let s = e1 + e2 + 1;
                    exps.push(s);
                    scalar = scalar * Q::from(-1)
                        / (fact(e1 as u64) * fact(e2 as u64) * Q::from(s as i64));
                }
                // Q(0)：插值
                let q0 = q0_cache
                    .entry(exps.clone())
                    .or_insert_with(|| {
                        let deg: u32 = 2 * exps.iter().sum::<u32>();
                        let mut xs = Vec::new();
                        let mut ys = Vec::new();
                        for m in (m0 + 1)..=(m0 + deg as i64 + 2) {
                            let mut q = balanced_weight_sum(&edges_v, &exps, &given, m);
                            // b1 ≥ 0（连通图）
                            for _ in 0..b1.max(0) {
                                q /= Q::from(m);
                            }
                            xs.push(m);
                            ys.push(q);
                        }
                        interpolate_at_zero(&xs, &ys)
                    })
                    .clone();
                let coeff = scalar * q0 * Q::from_signeds(1, 2i64.pow(r));
                if !is_zero(&coeff) {
                    // 展开成单项式列表（装饰中 ψ 只在标记与结点半边上）
                    terms.push(DecStratum {
                        graph: graph.clone(),
                        poly: KappaPsiPolynomial::new(vec![crate::kppoly::KpTerm {
                            kappa: t.kappa.clone(),
                            psi: t.psi.clone(),
                            coeff,
                        }]),
                    });
                }
            }
        }
    }
    Ok(TautClass::new(terms))
}

/// ∫_{M̄_{g,n}} DR · X，X 为基本类上的 κ-ψ 单项式。
/// 限制公式（Arbarello–Cornalba 的 κ，结点半边在各顶点上算作标记）：
/// κ_j|_Γ = Σ_v κ_j^{(v)}；ψ_i|_Γ = ψ_i（其顶点）。求值见 classes::int_against。
pub fn int_dr(
    g: u32,
    dvector: &[i64],
    r: Option<u32>,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
) -> Result<Q, String> {
    crate::kppoly::check_marking_legs(psi, dvector.len() as u32)?;
    let cls = dr_cycle(g, dvector, r)?;
    // DR 的生成元不带 κ 装饰，int_against 的"保留原项 κ"在此退化为空
    crate::classes::int_against(&cls, kappa, psi)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// admcycles 文档锚点：DR(1,(2,3,−5)) =
    /// 2ψ₁ + 9/2ψ₂ + 25/2ψ₃（基本类） −1/24（Δ₀） −25/2 −9/2 −2（三个分离边界）。
    /// 归一化注意：admcycles 的 δ_irr 生成元 = 2×物理轨迹类；我们除以 |Aut|=2，
    /// 故我们的 Δ₀ 系数应为 −1/12。
    #[test]
    fn dr1_anchor() {
        let cls = dr_cycle(1, &[2, 3, -5], None).unwrap();
        let mut found_fund = false;
        let mut found_d0 = false;
        let mut sep_vals = Vec::new();
        for t in &cls.terms {
            let g = &t.graph;
            if g.vertices().len() == 1 && g.edges().is_empty() {
                // 基本类：多项式应为 2ψ₁+9/2ψ₂+25/2ψ₃
                let p = &t.poly;
                assert_eq!(p.terms.len(), 3);
                for term in &p.terms {
                    let (l, e) = term.psi[0];
                    assert_eq!(e, 1);
                    let expect = Q::from(dvector_sq(l)) / Q::from(2);
                    assert_eq!(term.coeff, expect, "ψ_{l}");
                }
                found_fund = true;
            } else if g.vertices().len() == 1 && g.edges().len() == 1 {
                // Δ₀：我们 −1/12（= 文档 −1/24 × 2）
                assert_eq!(t.poly.terms[0].coeff, Q::from_signeds(-1, 12));
                found_d0 = true;
            } else {
                // 分离边界：−d²/2（单标记一侧）
                sep_vals.push(t.poly.terms[0].coeff.clone());
            }
        }
        assert!(found_fund && found_d0);
        sep_vals.sort();
        let expect = [
            Q::from_signeds(-25, 2),
            Q::from_signeds(-9, 2),
            Q::from_signeds(-2, 1),
        ];
        let mut e = expect.to_vec();
        e.sort();
        // 排序比较（malachite Q 的 Ord 是数值序？并非——手动比对集合）
        for v in &sep_vals {
            assert!(e.contains(v), "分离项 {v} 不在期望集");
        }
        assert_eq!(sep_vals.len(), 3);
    }

    fn dvector_sq(l: u32) -> i64 {
        [2i64, 3, -5][(l - 1) as usize] * [2i64, 3, -5][(l - 1) as usize]
    }

    /// r > g 时恒为零（CJ18）。用积分判据：类为零 ⟺ 对一组张成补空间的
    /// X 全部积分为零（逐项形式消没依赖生成元基的选取，不作判据）。
    #[test]
    fn dr_vanishing_above_genus() {
        let xs: Vec<(Vec<(u32, u32)>, Vec<(u32, u32)>)> = vec![
            (vec![], vec![(1, 1)]),
            (vec![], vec![(1, 2), (2, 1), (3, 1)]),
            (vec![(1, 1)], vec![(1, 1), (2, 1)]),
            (vec![(2, 1)], vec![(3, 1)]),
            (vec![(1, 2)], vec![]),
        ];
        for (k, p) in &xs {
            let v = int_dr(1, &[1, 3, -4], Some(2), k, p).unwrap();
            assert!(v == sc_num::Q::ZERO, "∫DR·X 非零（k={k:?} p={p:?}）：{v}");
        }
    }

    /// λ_g 与 DR(g,()) 的关系：admcycles 实测两值（DR_cycle 归一化下
    /// g=2 同号、g=3 反号；其 doctest 恒等式用的是另一入口
    /// double_ramification_cycle 的约定）。我们两条路线与其对应路线逐值一致。
    #[test]
    fn lambda_vs_dr_admcycles_anchored() {
        use crate::hodge::klp;
        // g=3 κ₃：λ 路 = 31/967680（= admcycles lambdaclass）；DR 路 = −31/967680（= admcycles DR_cycle）
        let lam = klp(3, 0, &[(3, 1)], &[], &[(3, 1)], &[]).unwrap();
        assert_eq!(lam.to_string(), "31/967680");
        let dr = int_dr(3, &[], None, &[(3, 1)], &[]).unwrap();
        assert_eq!(dr.to_string(), "-31/967680");
        // g=2 κ₁：两者同号 7/5760
        let lam = klp(2, 0, &[(1, 1)], &[], &[(2, 1)], &[]).unwrap();
        assert_eq!(lam.to_string(), "7/5760");
        let dr = int_dr(2, &[], None, &[(1, 1)], &[]).unwrap();
        assert_eq!(dr.to_string(), "7/5760");
    }

    /// 腿号 n+1.. 是结点半边：用户 ψ 越界必须报错，而非落到结点上。
    #[test]
    fn leg_validation() {
        assert!(int_dr(2, &[1, 1, -2], Some(2), &[], &[(4, 1)]).is_err());
        assert!(int_dr(2, &[1, 1, -2], Some(2), &[], &[(0, 1)]).is_err());
        assert!(int_dr(1, &[], None, &[], &[]).is_err()); // M̄_{1,0} 不稳定
        // 同腿 ψ 合并
        assert_eq!(
            int_dr(2, &[1, 1, -2], Some(2), &[], &[(1, 1), (1, 1), (2, 2)]).unwrap(),
            int_dr(2, &[1, 1, -2], Some(2), &[], &[(1, 2), (2, 2)]).unwrap()
        );
    }

    /// g=2 r=2 的 DR·X 五连（admcycles 2026-10-05 实测值）。
    #[test]
    fn dr2_integrals_vs_admcycles() {
        let cases: Vec<(Vec<(u32, u32)>, Vec<(u32, u32)>, &str)> = vec![
            (vec![(1, 1)], vec![(1, 3)], "13/144"),
            (vec![], vec![(1, 2), (2, 2)], "1/120"),
            (vec![(1, 1)], vec![(1, 1), (2, 1), (3, 1)], "7/48"),
            (vec![(2, 1)], vec![(1, 1), (2, 1)], "5/48"),
            (vec![], vec![(1, 1), (2, 1), (3, 2)], "1/80"),
        ];
        for (k, p, expect) in &cases {
            let v = int_dr(2, &[1, 1, -2], Some(2), k, p).unwrap();
            assert_eq!(&v.to_string(), expect, "k={k:?} p={p:?}");
        }
    }

    /// DR_compute(1,1,2,(2,−2)) 文档值：ψ₁、ψ₂ 系数 2、2；Δ₀（我们口径）−1/12。
    #[test]
    fn dr1_two_points() {
        let cls = dr_cycle(1, &[2, -2], None).unwrap();
        for t in &cls.terms {
            let g = &t.graph;
            if g.vertices().len() == 1 && g.edges().is_empty() {
                assert_eq!(t.poly.terms.len(), 2);
                for term in &t.poly.terms {
                    assert_eq!(term.coeff, Q::from(2));
                }
            } else if g.edges().len() == 1 && g.vertices().len() == 1 {
                assert_eq!(t.poly.terms[0].coeff, Q::from_signeds(-1, 12));
            } else {
                panic!("(1,2) 不应有分离边界项");
            }
        }
    }
}
