"""GRRcomp 差分：chern_single 各 (t,l,d,a) 的积分 vs oracle。"""
import json
import time

from admcycles import psiclass, kappaclass
from admcycles.GRRcomp import generalized_hodge_chern_single, twist_div_to_zero
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
    # (t, l, d, a, g, n, X_kappa, X_psi)
    (1, 0, [1, -1], [[1, [], -1]], 1, 2, [], [[1, 1]]),          # 主图 psi
    (1, 0, [1, -1], [[1, [], -1]], 1, 2, [[1, 1]], []),          # kappa
    (2, 1, [0, 0], [], 2, 2, [[1, 1]], [[1, 1]]),                # ch2 · k1 psi1
    (2, 1, [0, 0], [], 2, 2, [[2, 1]], []),                      # ch2 · k2
    (2, 0, [1, -1], [[1, [], -1]], 1, 2, [], [[1, 2]]),          # ch2 g=1
    (1, 2, [0, 0], [], 2, 1, [[1, 1]], []),                      # l=2 ch1 · k1
    (1, 1, [1, 0], [], 1, 2, [], [[1, 1]]),                      # d=[1,0]
    (2, 1, [0], [], 3, 1, [[1, 1], [2, 1]], []),                 # g=3 ch2
    (1, 0, [2, -1], [[1, [1], -1], [0, [2], 1]], 1, 2, [], [[2, 1]]),  # 多 a 分量
]


def main():
    for (t, l, d, a, g, n, kap, psi) in CASES:
        t0 = time.time()
        ch = generalized_hodge_chern_single(t, l, d, a, g, n)
        tc = time.time() - t0
        t0 = time.time()
        v = (ch * X_term(g, n, kap, psi)).evaluate()
        te = time.time() - t0
        print(json.dumps({
            "t": t, "l": l, "d": d, "a": a, "g": g, "n": n,
            "kappa": [list(x) for x in kap], "psi": [list(x) for x in psi],
            "value": str(v), "t_class": round(tc, 2), "t_eval": round(te, 2),
        }), flush=True)


if __name__ == "__main__":
    main()
