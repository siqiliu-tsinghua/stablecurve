"""GRRcomp 差分第二批：t=2 与 g=2 强例（度数严格匹配）。"""
import json
from admcycles import psiclass, kappaclass
from admcycles.GRRcomp import generalized_hodge_chern_single
from admcycles import TautologicalRing

CASES = [
    (2, 1, [0, 0], [], 2, 2, [[1, 1]], [[1, 2]]),   # k1 psi1^2
    (2, 1, [0, 0], [], 2, 2, [[2, 1]], [[1, 1]]),   # k2 psi1
    (2, 1, [0, 0], [], 2, 2, [[1, 1], [1, 1]], []), # k1^2? kappa=[(1,2)]
    (2, 1, [0, 0], [], 2, 2, [[3, 1]], []),         # k3
    (1, 1, [0], [], 2, 1, [[1, 1]], [[1, 2]]),      # g=2 n=1 ch1 · k1 psi1^2
    (1, 1, [0], [], 2, 1, [[3, 1]], []),            # ch1 · k3
    (2, 0, [1, -1], [[1, [], -1]], 1, 2, [[1, 1]], []),  # g=1 t=2 · k1
    (2, 0, [1, -1], [[1, [], -1]], 1, 2, [], [[1, 1]]),  # g=1 t=2 · psi1
    (1, 0, [1, 1, -2], [[1, [1], -1]], 1, 3, [], [[1, 1], [2, 1]]),
]

for (t, l, d, a, g, n, kap, psi) in CASES:
    ch = generalized_hodge_chern_single(t, l, d, a, g, n)
    R = TautologicalRing(g, n)
    x = R.one()
    for (i, c) in kap:
        x = x * kappaclass(i, g, n) ** c
    for (ll, p) in psi:
        x = x * psiclass(ll, g, n) ** p
    print(json.dumps({"t": t, "l": l, "d": d, "a": a, "g": g, "n": n,
                      "kappa": [list(z) for z in kap], "psi": [list(z) for z in psi],
                      "value": str((ch * x).evaluate())}), flush=True)
