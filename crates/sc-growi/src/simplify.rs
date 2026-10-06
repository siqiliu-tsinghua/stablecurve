//! 标准公理化简——原版 `libgw/simplify.cc` 的直译。
//!
//! 公理链的顺序是算法语义的一部分（便宜在前，REDESIGN §1.4）：
//! 度 0 映射到点 → 排序/归约 → 缓存查询 → 余切线公理 → Pⁿ 单点不变量 →
//! 基本类公理 → 除子公理 → g>1 无点加点。
//!
//! 语义对齐注意：原版度 0 的 mapping-to-point 结果**不入库**（便宜）；
//! 其余公理成功后先入库再返回。返回 `None` = 本函数算不出（对应原版
//! UNDEFINED），且此时入参已被排成规范形。

use crate::cache::Cache;
use crate::inv::{Inv, Pdata};
use crate::num::{One, Q, Zero, q, q_frac};

use crate::ctx::Ctx;

/// 度 0 不变量：映射到点的公理。对应原版 `mapping_to_point`。
fn mapping_to_point(iv: &Inv) -> Option<Q> {
    if iv.genus == 0 {
        // 亏格 0：多项式系数 × 普通交积
        let mut result = Q::from(iv.deg as i64); // 顶部 H 幂给因子 E
        let mut k: i64 = 0; // 分子阶乘计数
        let mut l: i64 = 0; // H 幂累积；负值 = 原始类标记
        for p in iv.pts.iter() {
            if p.h > 0 && l < 0 {
                return Some(Q::ZERO); // H × 原始类 → 0
            }
            if p.h < 0 {
                if l == 0 {
                    l = p.h as i64;
                } else if l != p.h as i64 {
                    return Some(Q::ZERO); // 两个不同原始类之积 → 0
                } else {
                    // 同一原始类的平方 = (1/E)·H^top
                    l = iv.dim as i64;
                    result /= Q::from(iv.deg as i64);
                }
            } else {
                l += p.h as i64;
            }
            // ψ 幂给二项式系数（连乘 (k+1)/j）
            for j in 1..=p.psi as i64 {
                k += 1;
                result *= Q::from(k);
                result /= Q::from(j);
            }
        }
        if l != iv.dim as i64 {
            return Some(Q::ZERO);
        }
        return Some(result);
    }

    if iv.genus == 1 {
        // 亏格 1：见 Kontsevich-Manin；维数条件保证可化到单点情形
        if iv.pts.is_empty() {
            return Some(Q::ZERO);
        }
        let mut l: i64 = 0;
        for p in iv.pts.iter() {
            if p.h < 0 {
                return Some(Q::ZERO); // 不出现原始类
            }
            l += p.h as i64;
        }
        if l >= 2 {
            return Some(Q::ZERO);
        }
        if iv.pts.len() > 1 {
            return None; // 算不出（对应原版 UNDEFINED）
        }
        let (n, e) = (iv.dim as i64, iv.deg as i64);
        if l != 0 {
            // <H> = −(top−1 陈类) × <λ>，<λ> on M_{1,1} = 1/24
            let one_minus_e = q(1 - e);
            Some(
                (q_pow(&one_minus_e, n + 2)
                    - Q::ONE
                    - Q::from((n + 2) * e) * (q_frac(e, 2) * Q::from(n + 1) - Q::ONE))
                    / (Q::from(24 * e * e)),
            )
        } else {
            // <psi> = top 陈类 × 1/24
            let one_minus_e = q(1 - e);
            Some((q_pow(&one_minus_e, n + 2) + Q::from((n + 2) * e) - Q::ONE) / (Q::from(24 * e)))
        }
    } else {
        None // 其它亏格的度 0 不变量目前未知
    }
}

/// 平方乘幂（指数可为负）。
fn q_pow(x: &Q, n: i64) -> Q {
    let neg = n < 0;
    let mut e = n.unsigned_abs();
    let mut base = x.clone();
    let mut acc = Q::ONE;
    while e > 0 {
        if e & 1 == 1 {
            acc = &acc * &base;
        }
        base = &base * &base;
        e >>= 1;
    }
    if neg { Q::ONE / &acc } else { acc }
}

/// Pⁿ 的单点有理不变量：J-函数 ∏_{i=1..d} 1/(H+i)^{n+1} 的递推。
/// 对应原版 `projspace_1point`。允许也要求 iw 可变（原版就地递推）。
fn projspace_1point(ctx: &mut Ctx, iw: &mut Inv) -> Option<Q> {
    if iw.pts[0].h < 0 {
        return Some(Q::ZERO); // 无原始类
    }

    if iw.cdeg == 1 {
        // 递推起点：binomial(−n−1, n−h)
        let mut result = Q::ONE;
        let mut j = -iw.dim as i64 - 1;
        for i in 1..=(iw.dim as i64 - iw.pts[0].h as i64) {
            result *= Q::from(j);
            j -= 1;
            result /= Q::from(i);
        }
        return Some(result);
    }

    // 递推步：二项展开 1/(H+d)^{n+1} 直接编入
    let mut result = Q::ZERO;
    let mut j: i64 = 0;
    let mut coeff = q_pow(&Q::from(iw.cdeg as i64), -(iw.dim as i64 + 1));
    iw.cdeg -= 1;
    iw.pts[0].psi -= iw.dim + 1;
    while iw.pts[0].h <= iw.dim {
        result += &coeff * ctx.compute(iw)?.clone();
        j += 1;
        coeff *= Q::from(-(iw.dim as i64 + j));
        coeff /= Q::from(j * (iw.cdeg as i64 + 1));
        iw.pts[0].h += 1;
        iw.pts[0].psi -= 1;
    }
    Some(result)
}

/// 尝试用标准公理计算不变量。成功则入库并返回；失败返回 `None` 且
/// `iv` 已被排成规范形（原始类已归约）。
pub fn simplify(ctx: &mut Ctx, iv: &mut Inv) -> Option<Q> {
    // 度 0 → 映射到点公理（不入库）
    if iv.cdeg == 0 && iv.genus >= 0 {
        if let Some(r) = mapping_to_point(iv) {
            return Some(r);
        }
    }

    // 排序 + 原始类归约
    iv.canonicalize();
    if iv.deg != 1 {
        iv.reduce();
    }

    // 缓存查询
    if let Some(v) = ctx.cache.get(iv) {
        return Some(v);
    }

    // 相对不变量与负亏格不变量没有可用的公理
    if !iv.pts.is_empty() && iv.pts[iv.pts.len() - 1].mult != 0 {
        return None;
    }
    if iv.genus < 0 {
        return None;
    }

    let npts = iv.pts.len();

    // 余切线公理（绝对不变量）
    if !iv.pts.is_empty() && iv.pts[npts - 1].psi != 0 && (iv.genus <= 1 || npts > 1) {
        for i in 0..npts {
            if iv.pts[i].psi > 1 {
                break;
            }
            if iv.pts[i].h != 0 || iv.pts[i].psi != 1 {
                continue;
            }
            let mut iw = iv.clone();
            iw.del(i);
            let n = iw.pts.len() as i64;
            let r = if 2 * iv.genus as i64 - 2 + n != 0 {
                ctx.compute(&iw)
                    .map(|v| v * Q::from(2 * iv.genus as i64 - 2 + n))
            } else {
                Some(Q::ZERO)
            };
            if let Some(r) = r {
                return store(ctx, iv, r);
            }
            break; // 一个余切线类不行，其余也不行
        }
    }

    // Pⁿ 单点不变量
    if npts == 1 && iv.genus == 0 && iv.deg == 1 {
        let mut iw = iv.clone();
        let r = projspace_1point(ctx, &mut iw);
        if let Some(r) = r {
            return store(ctx, iv, r);
        }
    }

    // 基本类公理
    if iv.genus <= 1 || npts > 1 {
        for i in 0..npts {
            if iv.pts[i].h > 0 || iv.pts[i].psi != 0 {
                break;
            }
            if iv.pts[i].h != 0 {
                continue; // 原始类
            }
            let r = if iv.pts[npts - 1].psi == 0 {
                // 无 descendant → 不变量为 0
                Some(Q::ZERO)
            } else {
                let mut iw = iv.clone();
                iw.del(i);
                let mut acc = Q::ZERO;
                let mut ok = true;
                for j in (0..iw.pts.len()).rev() {
                    if iw.pts[j].psi == 0 {
                        break;
                    }
                    iw.pts[j].psi -= 1;
                    match ctx.compute(&iw) {
                        Some(v) => acc += v,
                        None => {
                            ok = false;
                            break;
                        }
                    }
                    iw.pts[j].psi += 1;
                }
                if ok { Some(acc) } else { None }
            };
            if let Some(r) = r {
                return store(ctx, iv, r);
            }
            break;
        }
    }

    // 除子公理（原始类做除子（维数 2 时）不在此处理）
    if iv.genus >= 1 || npts > 1 {
        for i in 0..npts {
            if iv.pts[i].psi != 0 || iv.pts[i].h > 1 {
                break;
            }
            if iv.pts[i].h == 0 {
                continue;
            }
            let d = iv.cdeg as i64;
            let mut iw = iv.clone();
            iw.del(i);
            // 主项：DD × <去掉除子类的不变量>（原版 if (DD) result = DD*compute(iw)）
            let mut base = if d != 0 {
                Q::from(d) * ctx.compute(&iw)?
            } else {
                Q::ZERO
            };
            let mut ok = true;
            for j in (0..iw.pts.len()).rev() {
                if iw.pts[j].psi == 0 {
                    break;
                }
                if iw.pts[j].h < iw.dim && iw.pts[j].h >= 0 {
                    iw.pts[j].psi -= 1;
                    iw.pts[j].h += 1;
                    match ctx.compute(&iw) {
                        Some(v) => base += v,
                        None => {
                            ok = false;
                            break;
                        }
                    }
                    iw.pts[j].h -= 1;
                    iw.pts[j].psi += 1;
                }
            }
            if ok {
                return store(ctx, iv, base);
            }
            break;
        }
    }

    // g>1 无标记点：加一个点（ψ 类）以便 Virasoro 处理
    if iv.pts.is_empty() && iv.genus > 1 {
        let mut iw = iv.clone();
        iw.pts = vec![Pdata::new(0, 1, 0)].into();
        let r = ctx
            .compute(&iw)
            .map(|v| v / Q::from(2 * iv.genus as i64 - 2));
        if let Some(r) = r {
            return store(ctx, iv, r);
        }
    }

    None
}

fn store(ctx: &mut Ctx, iv: &Inv, r: Q) -> Option<Q> {
    ctx.cache.put(iv.clone(), r.clone());
    Some(r)
}

// Cache 仅供 doc 链接引用，避免 unused import 警告
#[allow(unused)]
fn _t(_: &Cache) {}
