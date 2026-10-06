"""sc-faber 全量验收：Instructions、DR 锚点、缓存。"""
from fractions import Fraction

import pytest
from sc_faber import FABER


@pytest.fixture(scope="module")
def f(tmp_path_factory):
    return FABER(cache=str(tmp_path_factory.mktemp("fb") / "cache.sqlite"))


def test_instructions_klp(f):
    # Faber《Instructions》数值例（g<=4 全部）
    assert f.klp(2, 2, psi=[(1, 2)], lam=[(1, 3)]) == Fraction(1, 1440)
    assert f.klp(2, 2, kappa=[(2, 1)], psi=[(2, 1)], lam=[(1, 2)]) == Fraction(7, 720)
    assert f.klp(3, 0, kappa=[(1, 1)], lam=[(3, 1), (2, 1)]) == Fraction(1, 120960)
    assert f.klp(3, 1, psi=[(1, 2)], lam=[(3, 1), (2, 1)]) == Fraction(1, 120960)
    assert f.klp(4, 2, psi=[(1, 2), (2, 2)], lam=[(4, 1), (3, 1)]) == Fraction(1, 276480)
    assert f.klp(4, 1, psi=[(1, 3)], lam=[(4, 1), (3, 1)]) == Fraction(1, 3225600)
    assert f.klp(3, 0, lam=[(1, 6)]) == Fraction(1, 90720)
    assert f.klp(4, 0, lam=[(1, 9)]) == Fraction(1, 113400)


def test_pure_psi_kappa(f):
    # Kontsevich 与 socle（admcycles 差分过的值）
    assert f.klp(0, 4, psi=[(1, 1)]) == Fraction(1)
    assert f.klp(0, 6, psi=[(1, 2), (2, 1)]) == Fraction(3)
    assert f.klp(2, 0, kappa=[(1, 3)]) == Fraction(43, 2880)
    assert f.klp(2, 1, psi=[(1, 4)]) == Fraction(1, 1152)


def test_dr(f):
    # DR 锚点（admcycles 实测）
    assert f.dr(1, [1, -1], r=1, kappa=[(1, 1)]) == Fraction(1, 24)
    assert f.dr(1, [1, -1], r=1, psi=[(1, 1)]) == Fraction(0)
    assert f.dr(2, [1, 1, -2], r=2, kappa=[(1, 1)], psi=[(1, 3)]) == Fraction(13, 144)
    assert f.dr(2, [1, 1, -2], r=2, psi=[(1, 2), (2, 2)]) == Fraction(1, 120)


def test_dr_vanishing(f):
    assert f.dr(1, [1, 3, -4], r=2, psi=[(1, 1)]) == Fraction(0)
    assert f.dr(1, [1, 3, -4], r=2, kappa=[(1, 2)]) == Fraction(0)


def test_lambda_dr_relation(f):
    # λ_g 与 DR(g,())：g=2 同号、g=3 反号（admcycles 双侧实测）
    assert f.klp(2, 0, kappa=[(1, 1)], lam=[(2, 1)]) == Fraction(7, 5760)
    assert f.dr(2, [], r=2, kappa=[(1, 1)]) == Fraction(7, 5760)
    assert f.klp(3, 0, kappa=[(3, 1)], lam=[(3, 1)]) == Fraction(31, 967680)
    assert f.dr(3, [], r=3, kappa=[(3, 1)]) == Fraction(-31, 967680)


def test_errors(f):
    with pytest.raises(ValueError):
        f.klp(2, 2, psi=[(9, 1)])  # 腿越界
    with pytest.raises(ValueError):
        f.dr(2, [1, 1])  # ΣA ≠ 0


def test_cache_roundtrip(tmp_path):
    c = tmp_path / "c.sqlite"
    f1 = FABER(cache=str(c))
    assert f1.klp(2, 1, psi=[(1, 4)]) == Fraction(1, 1152)
    f2 = FABER(cache=str(c))
    assert f2.klp(2, 1, psi=[(1, 4)]) == Fraction(1, 1152)  # L2 暖命中


def test_input_normalization(f):
    # 同腿 ψ 幂合并；ψ 腿号越界（会撞上结点半边编号）报错
    assert f.klp(2, 1, psi=[(1, 2), (1, 2)]) == Fraction(1, 1152)
    with pytest.raises(ValueError):
        f.dr(2, [1, 1, -2], r=2, psi=[(4, 1)])
    # 秩外 λ 为 0
    assert f.klp(2, 2, psi=[(1, 2)], lam=[(3, 1)]) == Fraction(0)


def test_cache_is_process_global(tmp_path):
    c = tmp_path / "g.sqlite"
    a = FABER(cache=str(c))
    assert a.cache_attached
    b = FABER(cache="")
    assert not b.cache_attached
    assert not a.cache_attached  # L2 为进程全局，b 已将其卸下
