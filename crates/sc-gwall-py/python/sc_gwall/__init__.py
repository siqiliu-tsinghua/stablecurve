"""sc-gwall：光滑曲线的（不连通）Gromov-Witten 不变量（Jim Bryan GWall 的 Rust 移植）。

用法：
    from sc_gwall import GWALL
    w = GWALL()
    w.gw(2, 1, [("t1", 0), ("pt", 0), ("pt", 0)])   # Fraction
    w.bps(3, 2, 1)                                   # Fraction(4)
"""

from ._sc_gwall import GWALL

__all__ = ["GWALL"]
__version__ = "0.2.0"
