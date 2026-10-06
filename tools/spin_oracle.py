#!/usr/bin/env python3
"""drpoly（Chiodo 插值主路径）/ drspin（spin DR）对拍：oracle 侧。

∫ DR_cycle(g,A,d=g,k,chiodo_coeff=True[,spin])·X，X 为 κ-ψ 单项式。
输出 JSONL：{"op":..,"g":..,"A":[..],"k":..,"kappa":[..],"psi":[..],"value":"p/q"}
"""
import json
import time

from admcycles import DR_cycle, psiclass, kappaclass
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
    # drpoly g=1, n=3, X deg 2
    ("drpoly", 1, (2, 3, -5), None, [], [[1, 1], [2, 1]]),
    ("drpoly", 1, (2, 3, -5), None, [[2, 1]], []),
    ("drpoly", 1, (2, 3, -5), None, [], [[3, 2]]),
    ("drpoly", 1, (2, 3, -5), None, [[1, 1]], [[2, 1]]),
    # drpoly g=2, n=3, k=1, X deg 4
    ("drpoly", 2, (3, 4, -2), 1, [[4, 1]], []),
    ("drpoly", 2, (3, 4, -2), 1, [[2, 2]], []),
    ("drpoly", 2, (3, 4, -2), 1, [[1, 1], [3, 1]], []),
    # drspin g=1, n=2, X deg 2
    ("drspin", 1, (1, 1), 1, [[2, 1]], []),
    ("drspin", 1, (1, 1), 1, [], [[1, 2]]),
    ("drspin", 1, (1, 1), 1, [], [[1, 1], [2, 1]]),
    ("drspin", 1, (3, -1), 1, [[2, 1]], []),
    ("drspin", 1, (3, -1), 1, [], [[1, 1], [2, 1]]),
    # drspin g=2, n=3, k=1, X deg 4
    ("drspin", 2, (9, -3, -1), 1, [[4, 1]], []),
    ("drspin", 2, (9, -3, -1), 1, [[1, 1], [3, 1]], []),
]


def main():
    for (op, g, A, k, kappa, psi) in CASES:
        N = len(A)
        spin = op == "drspin"
        t0 = time.time()
        D = DR_cycle(g, A, g, k, chiodo_coeff=True, spin=spin)
        tC = time.time() - t0
        t0 = time.time()
        v = (D * X_term(g, N, kappa, psi)).evaluate()
        tE = time.time() - t0
        print(json.dumps({
            "op": op, "g": g, "A": list(A), "k": k,
            "kappa": [list(x) for x in kappa], "psi": [list(x) for x in psi],
            "value": str(v), "t_class": round(tC, 2), "t_eval": round(tE, 2),
        }), flush=True)


if __name__ == "__main__":
    main()
