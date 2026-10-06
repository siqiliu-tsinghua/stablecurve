#!/usr/bin/env python3
"""Hyperell/Biell 配对对拍：oracle 侧。

输出 JSONL：{"op":..,"g":..,"n":..,"m":..,"kappa":[..],"psi":[..],"value":"p/q"}
查询集按 batch 参数分组：
  1: Hyperell(1,1,1) deg X=2（4 例）
  2: Hyperell(2,0,0) deg X=3（3 例）
  4: Biell(2,0,0) deg X=0（1 例）
  3: Hyperell(2,2,1) deg X=4（2 例）
  5: Biell(3,0,0) deg X=4（2 例）
用法：python hyb_oracle.py 1 2 4  < /dev/null（参数=要跑的 batch 号）
"""
import json
import sys
import time

from admcycles import Hyperell, Biell, psiclass, kappaclass
from admcycles import TautologicalRing


def X_term(g, N, kappa, psi):
    R = TautologicalRing(g, N)
    t = R.one()
    for (i, c) in kappa:
        t = t * kappaclass(i, g, N) ** c
    for (l, p) in psi:
        t = t * psiclass(l, g, N) ** p
    return t


CASES = {
    1: [("hyperell", 1, 1, 1, [], [(2, 1), (3, 1)]),
        ("hyperell", 1, 1, 1, [], [(2, 2)]),
        ("hyperell", 1, 1, 1, [(1, 1)], [(1, 1)]),
        ("hyperell", 1, 1, 1, [(2, 1)], [])],
    2: [("hyperell", 2, 0, 0, [(3, 1)], []),
        ("hyperell", 2, 0, 0, [(1, 3)], []),
        ("hyperell", 2, 0, 0, [(1, 1), (2, 1)], [])],
    4: [("biell", 2, 0, 0, [[2, 1]], []),
        ("biell", 2, 0, 0, [[1, 2]], [])],
    3: [("hyperell", 2, 2, 1, [], [(1, 1), (2, 1), (3, 1), (4, 1)]),
        ("hyperell", 2, 2, 1, [(1, 1)], [(1, 1), (3, 1)])],
    6: [("hyperell", 1, 1, 1, [], [[1, 2]]),
        ("hyperell", 1, 1, 1, [], [[1, 1], [2, 1]]),
        ("hyperell", 1, 1, 1, [], [[2, 2]]),
        ("hyperell", 1, 1, 1, [], [[2, 1], [3, 1]]),
        ("hyperell", 1, 1, 1, [[2, 1]], []),
        ("hyperell", 1, 1, 1, [[1, 1]], [[1, 1]]),
        ("hyperell", 1, 1, 1, [[1, 1]], [[2, 1]])],
    5: [("biell", 3, 0, 0, [(4, 1)], []),
        ("biell", 3, 0, 0, [(2, 2)], [])],
    7: [("hyperell", 2, 1, 0, [], [[1, 3]]),
        ("hyperell", 2, 2, 0, [], [[1, 2], [2, 1]]),
        ("hyperell", 2, 2, 0, [[1, 1]], [[1, 2]]),
        ("hyperell", 2, 2, 0, [[3, 1]], []),
        ("hyperell", 3, 0, 0, [[5, 1]], []),
        ("hyperell", 2, 1, 1, [[2, 1]], [[1, 1]]),
        ("hyperell", 2, 1, 1, [[1, 1]], [[2, 1], [3, 1]]),
        ("hyperell", 2, 1, 1, [], [[1, 2], [2, 1], [3, 1]]),
        ("biell", 2, 1, 0, [], [[1, 3]]),
        ("biell", 2, 1, 0, [[1, 1]], [[1, 2]]),
        ("biell", 2, 1, 0, [[3, 1]], []),
        ("biell", 3, 1, 0, [[1, 1], [4, 1]], [])],
    8: [("biell", 2, 1, 0, [], [[1, 2]]),
        ("biell", 2, 1, 0, [[2, 1]], []),
        ("biell", 2, 1, 0, [[1, 1]], [[1, 1]]),
        ("hyperell", 2, 1, 1, [[2, 1]], [[2, 1]]),
        ("hyperell", 2, 1, 1, [[1, 2]], [[3, 1]])],
}

CLS = {"hyperell": Hyperell, "biell": Biell}


def main():
    batches = [int(a) for a in sys.argv[1:]] or [1]
    for b in batches:
        for (op, g, n, m, kappa, psi) in CASES[b]:
            N = n + 2 * m
            t0 = time.time()
            H = CLS[op](g, n, m)
            tH = time.time() - t0
            t0 = time.time()
            v = (H * X_term(g, N, kappa, psi)).evaluate()
            tE = time.time() - t0
            rec = {"op": op, "g": g, "n": n, "m": m,
                   "kappa": [list(k) for k in kappa], "psi": [list(p) for p in psi],
                   "value": str(v), "t_class": round(tH, 2), "t_eval": round(tE, 2)}
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
