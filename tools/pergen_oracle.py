"""oracle 逐生成元 evaluate（完整不截断），与混合实验逐项对照。"""
from admcycles import psiclass, tautgens
from admcycles.double_ramification_cycle import DR_coeff_new
from sage.all import QQ

gens = tautgens(2, 1, 3, decst=True)
psis = psiclass(1, 2, 1)
for i, t in enumerate(gens):
    c = DR_coeff_new(t, 2, 1, 3, (1,), -1, False, True, 2)
    if c == 0:
        continue
    v = (c * t * psis).evaluate()
    if v != 0:
        print("PG", i, v, t.gamma, t.poly)
tot = QQ(0)
for t in gens:
    c = DR_coeff_new(t, 2, 1, 3, (1,), -1, False, True, 2)
    if c == 0:
        continue
    tot += (c * t * psis).evaluate()
print("TOTAL", tot)
