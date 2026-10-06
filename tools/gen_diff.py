"""Per-generator chiodo DR coefficients from oracle (Theta debug aid)."""
from admcycles import tautgens
from admcycles.double_ramification_cycle import DR_coeff_new

gens = tautgens(2, 1, 3, decst=True)

n = 0
for t in gens:
    c = DR_coeff_new(t, 2, 1, 3, (1, 1), -1, False, True, 2)
    if c != 0:
        g = t.gamma
        kappa, psi = t.poly.monom[0]
        print("G", list(g._genera), [list(l) for l in g._legs], list(g._edges), kappa, dict(psi), c)
        n += 1
print("NONZERO", n, "of", len(gens))
