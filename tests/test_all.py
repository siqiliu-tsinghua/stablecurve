"""全量验收：README 与 gwcheck 的全部示例。"""
from fractions import Fraction

import pytest
from sc_growi import GROWI


@pytest.fixture(scope="module")
def g(tmp_path_factory):
    return GROWI(cache=str(tmp_path_factory.mktemp("gw") / "cache.sqlite"))


def test_readme_examples(g):
    assert g.gw("lines in cubic surface") == Fraction(27)
    assert g.gw("elliptic cubics in P^2 thru H^2:9") == Fraction(1)
    assert g.gw("G=2,D=5 in P^4 thru H^3*psi^22") == Fraction(-41369, 110075314176)


def test_kontsevich(g):
    for d, v in [(5, 87304), (10, 40739017561997799680)]:
        assert g.gw(f"rational D={d} in P^2 thru H^2:{3*d-1}") == Fraction(v)


def test_getzler_elliptic(g):
    assert g.gw("elliptic sextics in P^2 thru H^2:18") == Fraction(57435240)
    assert g.gw("elliptic septics in P^2 thru H^2:21") == Fraction(60478511040)
    assert g.gw("elliptic octics in P^2 thru H^2:24") == Fraction(96212546526096)


def test_ehx_genus2(g):
    for d, v in [(4, 27), (5, 36855), (6, 58444767), (8, 346860150644700)]:
        assert g.gw(f"G=2,D={d} in P^2 thru H^2:{3*d+1}") == Fraction(v)


def test_relative(g):
    # 过 4 点且切于给定直线的二次曲线（经典特征数 = 2）
    assert g.gw("conics in P^1 thru H^2:4,[1,2]") == Fraction(2)


def test_unknown_returns_none(g):
    # E>1 的高亏格超曲面不变量当前不支持
    assert g.gw("elliptic D=2 in quintic threefold") is None


def test_memory_only_and_zero_class():
    g = GROWI(cache="")
    assert g.gw("lines in quintic 3fold") == Fraction(2875)
    # 零类 ⇒ 不变量为 0（原版语义）
    assert g.gw("lines in P^2 thru H^2:2, H-H") == Fraction(0)
    with pytest.raises(ValueError):
        g.gw("D=65537 in P^2 thru H^2:2")
