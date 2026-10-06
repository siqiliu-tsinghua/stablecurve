"""Witten r-spin 移植锚点：Pm 多项式数值、edge factor、极限类积分。"""
from admcycles.witten import Pmpolynomial, rspin_edge_factor, rspin_leg_factor
from admcycles import Wittenrspin, psiclass, kappaclass

print("== Pm 数值（X, a）==")
for (m, x, a) in [(1, 5, 3), (1, 7, 0), (1, 6, -2), (2, 5, 4), (2, 9, 0), (3, 7, 2), (0, 5, 1)]:
    v = Pmpolynomial(m).subs(X=x, a=a) if m else 1
    print(f"PM {m} {x} {a} = {v}")

print("== rspin_leg_factor 数值（X=m 代入）==")
for (d, a, x) in [(1, 0, 9), (2, 0, 9), (1, 3, 9), (2, -1, 9), (3, 2, 11)]:
    aa = x + a if a < 0 else a
    print(f"LEG {d} {a} {x} = {Pmpolynomial(d).subs(X=x, a=aa)}")

print("== rspin_edge_factor(w1, m, d1, d2) ==")
for (w1, m, d1, d2) in [(0, 7, 0, 0), (1, 7, 0, 0), (2, 9, 1, 0), (0, 11, 1, 1), (3, 11, 2, 0)]:
    print(f"EDGE {w1} {m} {d1} {d2} = {rspin_edge_factor(w1, m, d1, d2)}")

print("== Wittenrspin 极限类积分 ==")
def Xterm(g, N, kappa, psi):
    from admcycles import TautologicalRing
    t = TautologicalRing(g, N).one()
    for (i, c) in kappa:
        t = t * kappaclass(i, g, N) ** c
    for (l, p) in psi:
        t = t * psiclass(l, g, N) ** p
    return t
cases = [
    (1, [2, (1, -2)], [], [[1, 1]]),
    (1, [2, (1, -2)], [[1, 1]], []),
    (1, [2, (1, -2)], [], [[2, 1]]),
    (2, [2], [[2, 1]], []),
    (2, [2], [[1, 2]], []),
    (1, [3, (1/2, -1), (1/2, -2)], [], [[1, 1], [2, 1]]),
    (1, [3, (1/2, -1), (1/2, -2)], [[1, 1]], [[3, 1]]),
    (2, [1, 1], [[2, 1]], []),
    (2, [1, 1], [], [[1, 1], [2, 1]]),
]
for (g, A, kap, psi) in cases:
    N = len(A)
    W = Wittenrspin(g, A)
    v = (W * Xterm(g, N, kap, psi)).evaluate()
    print(f"WITT {g} {A} {kap} {psi} = {v}")
