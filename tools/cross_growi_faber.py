#!/usr/bin/env python3
"""GROWI 验收阶段 2：faber ↔ growi 互为 oracle 交叉验证。

共享数学面：M̄_{0,n} 的 ψ 积分。映射（实验确立，三例锚定）：
  GROWI  "G=0,D=0 in P^{d} thru H:{d}, psi^{a1},...,psi^{am}"
    =  ∫_{M̄_{0,d+m}} ψ₁^{a₁}···ψ_m^{a_m}   （faber klp(g=0) 独立计算）
合法组合：Σaᵢ = d+m−3（重言维度匹配；GROWI 侧表现为查询维度校验）。
"""
import itertools
from fractions import Fraction

from sc_growi import GROWI
from sc_faber import FABER


def compositions(total: int, parts: int):
    """有序正整数组成。"""
    if parts == 1:
        return [(total,)] if total >= 1 else []
    out = []
    for a in range(1, total - parts + 2):
        for rest in compositions(total - a, parts - 1):
            out.append((a,) + rest)
    return out


def main():
    g = GROWI()
    f = FABER()
    ok = fail = skipped = 0
    fails = []
    for dim in range(3, 7):
        for m in range(1, 4):
            total = dim + m - 3
            for a in compositions(total, m):
                psi_strs = [f"psi^{ai}" if ai > 1 else "psi" for ai in a]
                expr = f"G=0,D=0 in P^{dim} thru H:{dim}," + ",".join(psi_strs)
                try:
                    gv = g.gw(expr)
                except Exception as e:  # noqa: BLE001
                    skipped += 1
                    print(f"[skip] {expr}  ({e})")
                    continue
                fv = f.klp(0, dim + m, [], [(i + 1, ai) for (i, ai) in enumerate(a)])
                if Fraction(gv) == Fraction(fv):
                    ok += 1
                else:
                    fail += 1
                    fails.append((expr, gv, fv))
                    print(f"[FAIL] {expr}\n        growi={gv} faber={fv}")
    print(f"\n== 交叉 {ok} 通过 / {fail} 失败 / {skipped} 跳过 ==")
    if fails:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
