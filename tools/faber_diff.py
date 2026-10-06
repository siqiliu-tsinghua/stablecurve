# -*- coding: utf-8 -*-
"""sc-faber vs admcycles 差分对拍。
用法（WSL 内）：python3 faber_diff.py <rust_jsonl 路径>
Rust 侧由 examples/diff_dump.rs 导出；admcycles 侧逐项同构造比对。
strata 语义对齐：Rust 的 evaluate 已除 |Aut|；admcycles 裸 decstratum 不除，
脚本对 stratum 项同时计算 raw 与 raw/aut，按匹配者计。"""
import json
import sys

from admcycles import psiclass, kappaclass, stgraph, lambdaclass, fundclass


def fund_value(g, n, psi, kappa):
    return fund_value_raw(g, n, psi, kappa).evaluate()


def fund_value_raw(g, n, psi, kappa):
    cls = None
    for (l, p) in psi:
        t = psiclass(l, g, n) ** p
        cls = t if cls is None else cls * t
    for (i, p) in kappa:
        t = kappaclass(i, g, n) ** p
        cls = t if cls is None else cls * t
    if cls is None:
        cls = fundclass(g, n)
    return cls


STRATA = {
    ("d0", 1, 1): ([0], [[1, 2, 3]], [(2, 3)]),
    ("d0", 2, 0): ([1], [[1, 2]], [(1, 2)]),
    ("d0", 2, 2): ([1], [[1, 2, 3, 4]], [(3, 4)]),
    ("d0", 1, 3): ([0], [[1, 2, 3, 4, 5]], [(4, 5)]),
    ("sep11", 2, 2): ([1, 1], [[1, 3], [2, 4]], [(3, 4)]),
}


def main(path):
    ok = bad = 0
    for line in open(path, encoding="utf-8"):
        d = json.loads(line)
        ours = d["value"]
        if d["kind"] == "fund":
            theirs = str(fund_value(d["g"], d["n"], d["psi"], d["kappa"]))
        elif d["kind"] == "lam":
            cls = None
            for (l, p) in d["lambda"]:
                t = lambdaclass(l, d["g"], d["n"]) ** p
                cls = t if cls is None else cls * t
            base = fund_value_raw(d["g"], d["n"], d["psi"], d["kappa"])
            cls = base if cls is None else cls * base
            theirs = str(cls.evaluate())
        else:
            genera, legs, edges = STRATA[(d["name"], d["g"], d["n"])]
            tc = stgraph(genera, legs, edges).to_tautclass()
            raw = tc.evaluate()
            aut = stgraph(genera, legs, edges).automorphism_number()
            cands = {str(raw), str(raw / aut)}
            theirs = ours if ours in cands else f"{raw}|{raw / aut}"
        if ours == theirs:
            ok += 1
        else:
            bad += 1
            print(f"MISMATCH g={d.get('g')} n={d.get('n')} psi={d.get('psi')} "
                  f"kappa={d.get('kappa')} name={d.get('name')}: "
                  f"ours={ours} admcycles={theirs}")
    print(f"一致 {ok} / 不一致 {bad}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
