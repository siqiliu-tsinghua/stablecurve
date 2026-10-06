"""原版 GROWI README 中 test.cm 示例的 Python 复现：P^3 中过 2d 个点的椭圆空间曲线数。

运行：uv run python examples/elliptic_space_curves.py
"""
from fractions import Fraction

from sc_growi import GROWI

g = GROWI(cache="")
for d in range(1, 9):
    n0 = g.gw(f"rational D={d} in P^3 thru H^3:{2*d}")
    n1 = g.gw(f"elliptic D={d} in P^3 thru H^3:{2*d}")
    print(d, n1 + Fraction(2 * d - 1, 12) * n0)
