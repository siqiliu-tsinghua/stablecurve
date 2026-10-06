"""Witten r-spin 差分：oracle 侧（g=2 为主）。A_i = (C_num, C_den, D)。"""
import json
import time

from admcycles import Wittenrspin, psiclass, kappaclass
from admcycles import TautologicalRing


def X_term(g, N, kappa, psi):
    R = TautologicalRing(g, N)
    t = R.one()
    for (i, c) in kappa:
        t = t * kappaclass(i, g, N) ** c
    for (l, p) in psi:
        t = t * psiclass(l, g, N) ** p
    return t


CASES = [
    (2, [[0, 1, 2]], [[3, 1]], []),                       # W(2,(2,))·k3
    (2, [[0, 1, 2]], [[1, 1], [2, 1]], []),               # ·k1k2
    (2, [[0, 1, 2]], [[2, 1]], [[1, 1]]),                 # ·k2 psi1
    (2, [[0, 1, 2]], [[1, 1]], [[1, 2]]),                 # ·k1 psi1^2
    (2, [[0, 1, 2]], [], [[1, 3]]),                       # ·psi1^3
    (2, [[0, 1, 3], [0, 1, -1]], [[1, 1]], [[1, 1], [2, 1]]),  # W(2,(3,-1))
    (2, [[0, 1, 3], [0, 1, -1]], [[2, 1]], [[1, 1]]),
    (2, [[0, 1, 3], [0, 1, -1]], [], [[1, 2], [2, 1]]),
    (2, [[0, 1, 1], [0, 1, 1]], [[1, 1]], [[1, 1], [2, 1]]),   # W(2,(1,1))
    (1, [[1, 1, 2], [1, 1, -3], [0, 1, 0]], [[1, 1]], [[1, 1]]),  # 分数 C
]


def main():
    for (g, A, kap, psi) in CASES:
        N = len(A)
        Avec = [(cn / cd) * None if False else ((cn, cd, dd)) for (cn, cd, dd) in A]
        # 构造 (C, D) 元组供 Wittenrspin
        from fractions import Fraction
        avec = []
        for (cn, cd, dd) in A:
            if cn == 0:
                avec.append(dd)
            else:
                avec.append((Fraction(cn, cd), dd))
        t0 = time.time()
        W = Wittenrspin(g, avec)
        tC = time.time() - t0
        t0 = time.time()
        v = (W * X_term(g, N, kap, psi)).evaluate()
        tE = time.time() - t0
        print(json.dumps({"g": g, "A": A, "kappa": [list(x) for x in kap],
                          "psi": [list(x) for x in psi], "value": str(v),
                          "t_class": round(tC, 2), "t_eval": round(tE, 2)}), flush=True)


if __name__ == "__main__":
    main()
