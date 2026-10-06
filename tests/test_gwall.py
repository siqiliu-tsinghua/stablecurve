"""sc-gwall 验收：BP 论文表锚点、插入写法、环路符号、错误映射。"""
from fractions import Fraction

import pytest
from sc_gwall import GWALL


@pytest.fixture(scope="module")
def w():
    return GWALL()


def test_bps_paper_values(w):
    # BP 论文 local 表 n_2^h(g) 与 §1 的 g=0/g=1 完整结果
    assert w.bps(3, 2, 1) == Fraction(4)
    assert w.bps(2, 2, 0) == Fraction(-2)
    assert w.bps(0, 1, 0) == Fraction(1)
    assert w.bps(0, 2, 0) == Fraction(0)
    assert w.bps(1, 3, 0) == Fraction(1)
    assert w.bps(1, 3, 1) == Fraction(0)


def test_series(w):
    assert w.l(2, 2, 2) == [Fraction(15, 2), Fraction(0), Fraction(-5, 4)]
    assert w.gw(2, 2, []) == w.ztil(2, 2, 0)[0]


def test_insert_forms(w):
    tup = w.gw(2, 1, [("t1", 0), ("pt", 1), ("pt", 0)])
    lst = w.gw(2, 1, [["t1", 0], ["pt", 1], ["pt", 0]])
    assert tup == lst
    assert w.gw(2, 1, (("t1", 0), ("pt", 1), ("pt", 0))) == tup


def test_loop_order_sign(w):
    canon = w.gw(2, 1, [("a", 0, 1), ("a", 1, 2), ("b", 0, 2), ("b", 1, 1)])
    swapped = w.gw(2, 1, [("a", 0, 1), ("b", 0, 2), ("a", 1, 2), ("b", 1, 1)])
    assert canon != 0
    assert swapped == -canon


def test_errors(w):
    with pytest.raises(ValueError):
        w.gw(2, 1, [("zz", 1)])
    with pytest.raises(ValueError):
        w.gw(2, 1, [("a", 1)])  # 环路需要 (k, i)
    with pytest.raises(ValueError):
        w.bps(2, 3, 9)  # 超出 branch ≤ 10 截断
