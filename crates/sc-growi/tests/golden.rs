//! 黄金测试：全链已知数值（来源：gwcheck.cm、README、经典文献）；
//! M2-M4 后含任意亏格、超曲面、五次三维流形序列（d=5 定案见下方注释）。

use sc_growi::{GwClass, GwSettings, Poly, Var, gw};

fn h2() -> Poly {
    Poly::from(Var::H).pow(2)
}

/// P^2 的 Kontsevich 数：过 3d−1 个点的 d 次有理曲线数。
#[test]
fn kontsevich_p2() {
    let mut ctx = sc_growi::Ctx::memory();
    let known = [
        (1i16, "1"),
        (2, "1"),
        (3, "12"),
        (4, "620"),
        (5, "87304"),
        (6, "26312976"),
        (7, "14616808192"),
    ];
    for &(d, v) in known.iter() {
        let s = GwSettings {
            dim: 2,
            deg: 1,
            genus: 0,
            cdeg: d,
        };
        let cls = vec![GwClass {
            poly: h2(),
            mult: 0,
            count: (3 * d - 1) as u32,
        }];
        let r = gw(&mut ctx, &s, &cls).unwrap_or_else(|e| panic!("d={d}: {e}"));
        assert_eq!(r.to_string(), v, "Kontsevich d={d}");
    }
    // d=8,9 只算不断言（暖缓存），d=10 对齐 gwcheck
    for d in 8..=10i16 {
        let s = GwSettings {
            dim: 2,
            deg: 1,
            genus: 0,
            cdeg: d,
        };
        let cls = vec![GwClass {
            poly: h2(),
            mult: 0,
            count: (3 * d - 1) as u32,
        }];
        let r = gw(&mut ctx, &s, &cls).unwrap_or_else(|e| panic!("d={d}: {e}"));
        if d == 10 {
            assert_eq!(r.to_string(), "40739017561997799680", "Kontsevich d=10");
        }
    }
}

/// 三次曲面上的直线：27。
#[test]
fn lines_on_cubic_surface() {
    let mut ctx = sc_growi::Ctx::memory();
    let s = GwSettings {
        dim: 2,
        deg: 3,
        genus: 0,
        cdeg: 1,
    };
    let r = gw(&mut ctx, &s, &[]).unwrap();
    assert_eq!(r.to_string(), "27");
}

/// 五次三维流形的有理曲线数（镜对称经典序列）。
///
/// d=5 澄清（codex 终局仲裁 + 本方 GV 多覆盖独立复算）：
/// GW 不变量 N₅ = 229305888887648 是**镜对称正确值**——GV 多覆盖公式
/// N_d = Σ_{k|d} n_{d/k}/k³ 给 N₅ = n₅ + 2875/5³ = ...625 + 23 = ...648
/// （d=1..6 六值全部自洽复算）。文献佐证：Dang (arXiv:1409.3712 §5
/// localization 独立计算)、Cooper–Zinger (Michigan Math. J. 63 (2014)
/// Table 1)、Givental 镜像公式级数直接展开。此前记录的"镜对称真值 …624"
/// 系我方转抄误差：枚举数 n₅ = ...625 误作 ...624（624 = 625 − 1），
/// "+24 偏差"实为 23（多覆盖项）+ 1（转抄误差）。gwcheck.cm 断言与
/// 独立数学一致，无任何程序缺陷。
#[test]
fn quintic_threefold() {
    let mut ctx = sc_growi::Ctx::memory();
    let known = [
        (1i16, "2875"),
        (2, "4876875/8"),
        (3, "8564575000/27"),
        (4, "15517926796875/64"),
        (5, "229305888887648"), // 镜对称正确值（GV 多覆盖复算 + 三文献）
        (6, "248249742157695375"),
    ];
    for &(d, v) in known.iter() {
        let s = GwSettings {
            dim: 3,
            deg: 5,
            genus: 0,
            cdeg: d,
        };
        let r = gw(&mut ctx, &s, &[]).unwrap_or_else(|e| panic!("d={d}: {e}"));
        assert_eq!(r.to_string(), v, "quintic d={d}");
    }
}

/// README 相对不变量示例：过 4 点且与给定直线相切的平面二次曲线。
/// GW(N=1,E=1,D=2,G=0,H^2:4,[1,2])
///
/// 原版实际输出 2（经差分预言机确认），与经典特征数（4 点 + 1 切线的
/// 二次曲线数 = 2）一致；README 写的 1 系笔误。
#[test]
fn conics_tangent_to_line() {
    let mut ctx = sc_growi::Ctx::memory();
    let s = GwSettings {
        dim: 1,
        deg: 1,
        genus: 0,
        cdeg: 2,
    };
    let cls = vec![
        GwClass {
            poly: h2(),
            mult: 0,
            count: 4,
        },
        GwClass {
            poly: Poly::from(1), // [1,2]：基本类带重数 2
            mult: 2,
            count: 1,
        },
    ];
    let r = gw(&mut ctx, &s, &cls).unwrap();
    assert_eq!(r.to_string(), "2");
}

/// 缓存行为：磁盘 L2 命中后第二次计算瞬时返回相同结果。
#[test]
fn sqlite_cache_reuse() {
    let dir = std::env::temp_dir().join("sc-growi-golden-cache");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("cache.sqlite");

    let s = GwSettings {
        dim: 2,
        deg: 3,
        genus: 0,
        cdeg: 1,
    };
    let v1 = {
        let mut ctx = sc_growi::Ctx::with_db(&path).unwrap();
        gw(&mut ctx, &s, &[]).unwrap()
    };
    let v2 = {
        let mut ctx = sc_growi::Ctx::with_db(&path).unwrap();
        gw(&mut ctx, &s, &[]).unwrap()
    };
    assert_eq!(v1, v2);
    assert_eq!(v1.to_string(), "27");
    let _ = std::fs::remove_dir_all(&dir);
}

/// gwcheck 的 Getzler 椭圆不变量（g=1）。
#[test]
fn getzler_genus1() {
    let mut ctx = sc_growi::Ctx::memory();
    let cases = [
        // (dim, d, H^dim 点数, 其余 H^k 点数, 期望)
        (3i16, 1i16, 2u32, 0u32, "-1/12"), // <H3:2>
        (3, 1, 0, 4, "-1/6"),              // <H2:4>（H^2 条件）
        (3, 1, 1, 2, "-1/12"),             // <H3,H2:2>
        (3, 5, 10, 0, "-147/4"),
        (3, 5, 5, 10, "-7179606"),
        (3, 5, 0, 20, "-1984020394752"),
        (2, 6, 18, 0, "57435240"),
        (2, 7, 21, 0, "60478511040"),
        (2, 8, 24, 0, "96212546526096"),
    ];
    for &(n, d, top, h2, expect) in cases.iter() {
        let s = GwSettings {
            dim: n,
            deg: 1,
            genus: 1,
            cdeg: d,
        };
        let mut cls = vec![GwClass {
            poly: Poly::from(Var::H).pow(n as u32),
            mult: 0,
            count: top,
        }];
        if h2 > 0 {
            cls.push(GwClass {
                poly: Poly::from(Var::H).pow(2),
                mult: 0,
                count: h2,
            });
        }
        let r = gw(&mut ctx, &s, &cls).unwrap_or_else(|e| panic!("n={n} d={d}: {e}"));
        assert_eq!(r.to_string(), expect, "g=1 n={n} d={d}");
    }
}

/// gwcheck 的 EHX/Caporaso-Harris 亏格 2 不变量。
#[test]
fn ehx_genus2() {
    let mut ctx = sc_growi::Ctx::memory();
    // P^2：d=1..8，点数 3d+1
    let p2 = [
        (1i16, "0"),
        (2, "0"),
        (3, "0"),
        (4, "27"),
        (5, "36855"),
        (6, "58444767"),
        (7, "122824720116"),
        (8, "346860150644700"),
    ];
    for &(d, v) in p2.iter() {
        let s = GwSettings {
            dim: 2,
            deg: 1,
            genus: 2,
            cdeg: d,
        };
        let cls = vec![GwClass {
            poly: Poly::from(Var::H).pow(2),
            mult: 0,
            count: (3 * d + 1) as u32,
        }];
        let r = gw(&mut ctx, &s, &cls).unwrap();
        assert_eq!(r.to_string(), v, "g=2 P2 d={d}");
    }
    // P^3 d=3：H^3:k + H^2:2(6-k)
    let p3 = [
        (6u32, "1/12"),
        (5, "5/12"),
        (4, "5/2"),
        (3, "46/3"),
        (2, "307/3"),
        (1, "747"),
        (0, "5930"),
    ];
    for &(k, v) in p3.iter() {
        let s = GwSettings {
            dim: 3,
            deg: 1,
            genus: 2,
            cdeg: 3,
        };
        let mut cls = vec![GwClass {
            poly: Poly::from(Var::H).pow(3),
            mult: 0,
            count: k,
        }];
        let h2 = 2 * (6 - k);
        if h2 > 0 {
            cls.push(GwClass {
                poly: Poly::from(Var::H).pow(2),
                mult: 0,
                count: h2,
            });
        }
        let r = gw(&mut ctx, &s, &cls).unwrap();
        assert_eq!(r.to_string(), v, "g=2 P3 k={k}");
    }
}

/// gwcheck 的 EHX 亏格 3 不变量（P^2，点数 3d+2）。
#[test]
fn ehx_genus3() {
    let mut ctx = sc_growi::Ctx::memory();
    let p2 = [
        (1i16, "0"),
        (2, "0"),
        (3, "0"),
        (4, "1"),
        (5, "7915"),
        (6, "34435125"),
        (7, "153796445095"),
        (8, "800457740515775"),
    ];
    for &(d, v) in p2.iter() {
        let s = GwSettings {
            dim: 2,
            deg: 1,
            genus: 3,
            cdeg: d,
        };
        let cls = vec![GwClass {
            poly: Poly::from(Var::H).pow(2),
            mult: 0,
            count: (3 * d + 2) as u32,
        }];
        let r = gw(&mut ctx, &s, &cls).unwrap();
        assert_eq!(r.to_string(), v, "g=3 d={d}");
    }
}

/// M4：五重切线问题（habilitation 例 4.1.10/4.1.11 真值）。
/// N_d = d(d−3)(d−4)(d⁷+12d⁶−18d⁵−540d⁴+311d³+5457d²−2133d−12690)/5!
/// 是 d 次平面上五重切二次曲线的虚拟数；GROWI 相对不变量
/// GW(N=1,E=d,D=2,[1,2]:5) = 5!·N_d（五个带标号的切点）。
#[test]
fn conics_fivefold_tangent() {
    let mut ctx = sc_growi::Ctx::memory();
    // 论文表值 N_d × 120
    let known = [
        (5i16, "238200"),
        (6, "8573040"),
        (7, "82547640"),
        (8, "467190720"),
    ];
    for &(e, v) in known.iter() {
        let s = GwSettings {
            dim: 1,
            deg: e,
            genus: 0,
            cdeg: 2,
        };
        let cls = vec![GwClass {
            poly: Poly::from(1),
            mult: 2,
            count: 5,
        }];
        let r = gw(&mut ctx, &s, &cls).unwrap();
        assert_eq!(r.to_string(), v, "5-fold tangent conics d={e}");
    }
}

/// M4：带 descendant 的相对不变量 + 原类平方（与差分预言机核对过）。
#[test]
fn relative_descendants_and_primitives() {
    let mut ctx = sc_growi::Ctx::memory();
    // <[H²ψ,3]>_2（P² rel 五次曲线）
    let s = GwSettings {
        dim: 1,
        deg: 5,
        genus: 0,
        cdeg: 2,
    };
    let r = gw(
        &mut ctx,
        &s,
        &[GwClass {
            poly: &Poly::from(Var::H).pow(2) * &Poly::from(Var::Psi),
            mult: 3,
            count: 1,
        }],
    )
    .unwrap();
    assert_eq!(r.to_string(), "0");

    // <[Hψ,2],[H,2]>_2 = 25/2
    let r = gw(
        &mut ctx,
        &s,
        &[
            GwClass {
                poly: &Poly::from(Var::H) * &Poly::from(Var::Psi),
                mult: 2,
                count: 1,
            },
            GwClass {
                poly: Poly::from(Var::H),
                mult: 2,
                count: 1,
            },
        ],
    )
    .unwrap();
    assert_eq!(r.to_string(), "25/2");

    // <[Hψ,4],[H,1]>_1 五次三维流形相对 = 525
    let s = GwSettings {
        dim: 3,
        deg: 5,
        genus: 0,
        cdeg: 1,
    };
    let r = gw(
        &mut ctx,
        &s,
        &[
            GwClass {
                poly: &Poly::from(Var::H) * &Poly::from(Var::Psi),
                mult: 4,
                count: 1,
            },
            GwClass {
                poly: Poly::from(Var::H),
                mult: 1,
                count: 1,
            },
        ],
    )
    .unwrap();
    assert_eq!(r.to_string(), "525");

    // 原类平方：pc(1)²= H³/E（引擎自动折 1/5 因子），值为 H³ψ 情形 1575/2 ÷ 5
    let s = GwSettings {
        dim: 3,
        deg: 5,
        genus: 0,
        cdeg: 2,
    };
    let r = gw(
        &mut ctx,
        &s,
        &[GwClass {
            poly: &Poly::from(Var::Prim(1)).pow(2) * &Poly::from(Var::Psi),
            mult: 8,
            count: 1,
        }],
    )
    .unwrap();
    assert_eq!(r.to_string(), "315/2");
    // pc(1)² 无 ψ：972 = 4860/5
    let r = gw(
        &mut ctx,
        &s,
        &[GwClass {
            poly: Poly::from(Var::Prim(1)).pow(2),
            mult: 9,
            count: 1,
        }],
    )
    .unwrap();
    assert_eq!(r.to_string(), "972");
}
