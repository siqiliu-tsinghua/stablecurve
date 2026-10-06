"""GRRcomp 移植锚点：基础件数值 + todd/chern/lambda 类结构。"""
from admcycles import *
from admcycles.GRRcomp import *
from sage.all import bernoulli_polynomial, QQ

print("== 基础件 ==")
print("DB12 =", divided_bernoulli(12))
for (l, n) in [(0, 4), (2, 4), (1, 2), (-1, 3), (3, 1)]:
    print(f"DBP({l},{n}) =", bernoulli_polynomial(l, n) / ZZ(n).factorial())

print("== todd_inv_exp(a=[], t=1, g=1, n=2) ==")
print(todd_inv_exp([], 1, 1, 2))

print("== todd_inv_exp_graphs 图数 ==")
print("len =", len(todd_inv_exp_graphs([[1, [1], -1]], 1, [1], 3, 2)._terms))
print(todd_inv_exp_graphs([[1, [1], -1]], 1, [1], 3, 2))

print("== chern_single(1, l=0, d=[1,-1], a=[[1,[],-1]], g=1, n=2) ==")
print(generalized_hodge_chern_single(1, 0, [1, -1], [[1, [], -1]], 1, 2))

print("== generalized_lambda(1, 1, [0], [], 2, 1) ==")
print(generalized_lambda(1, 1, [0], [], 2, 1))

print("== twist_div_to_zero(3,[3,4],2,2) ==")
print(twist_div_to_zero(3, [3, 4], 2, 2))
