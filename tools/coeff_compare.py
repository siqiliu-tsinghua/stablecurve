"""逐生成元系数比对：oracle 表（gen_diff.py 输出）vs 我们的 coeff_table 输出。
两侧各自生成规范键（边数 + 顶点(亏格,标记集,κ,标记ψ,结点ψ) 多重集），
多重集配对比较系数。"""
import ast
import re
import sys
from fractions import Fraction as F


def oracle_keys(path):
    out = []
    for line in open(path, encoding="utf-8"):
        m = re.match(
            r"G (\[.*?\]) (\[\[.*?\]\]|\[\]) (\[\(.*?\)\]|\[\]) (\[\[.*?\]\]|\[\]) (\{.*?\}) (-?\d+/\d+|-?\d+)",
            line.strip(),
        )
        if not m:
            continue
        genera = ast.literal_eval(m.group(1))
        legs = ast.literal_eval(m.group(2))
        edges = ast.literal_eval(m.group(3))
        kappa = ast.literal_eval(m.group(4))
        psi = ast.literal_eval(m.group(5))
        coeff = F(m.group(6))
        edge_legs = {l for e in edges for l in e}
        vks = []
        plain = []
        for i in range(len(genera)):
            marks = sorted(l for l in legs[i] if l not in edge_legs)
            kap = []
            if i < len(kappa):
                for j, c in enumerate(kappa[i]):
                    kap += [j + 1] * c
            mpsi = sorted((l, p) for l, p in psi.items() if l not in edge_legs and l in legs[i])
            npsi = sorted(p for l, p in psi.items() if l in edge_legs and l in legs[i])
            vks.append(f"(g{genera[i]},m{marks},k{kap},p{mpsi},n{npsi})")
            plain.append((genera[i], marks, kap))
        vks.sort()
        # 边对按【vks 排序后】的顶点序（与 Rust 侧一致）
        rank = {vi: r for r, vi in enumerate(sorted(range(len(genera)), key=lambda i: str(plain[i])))}

        def vtx_of(l):
            for i, vs in enumerate(legs):
                if l in vs:
                    return i
            raise ValueError

        epairs = sorted(
            tuple(sorted((rank[vtx_of(a)], rank[vtx_of(b)]))) for a, b in edges
        )
        out.append((f"E{len(edges)};{';'.join(vks)};e{epairs}", coeff))
    return out


def our_keys(path):
    out = []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        key, coeff = line.split("\t")
        out.append((key, F(coeff)))
    return out


def main():
    o = oracle_keys(sys.argv[1])
    u = our_keys(sys.argv[2])
    od = {}
    for k, c in o:
        od[k] = od.get(k, F(0)) + c
    print(f"oracle {len(o)} 行 / 我们 {len(u)} 行")
    print(f"我们独有键 {len(set(k for k, _ in u) - set(od))}：")
    for k in sorted(set(k for k, _ in u) - set(od)):
        cs = [c for kk, c in u if kk == k]
        print("  +", k, cs)
    only_o = set(od) - set(k for k, _ in u)
    print(f"oracle 独有键 {len(only_o)}：")
    for k in sorted(only_o):
        print("  -", k, od[k])
    both = 0
    bad = 0
    ud = {}
    for k, c in u:
        ud[k] = ud.get(k, F(0)) + c
    for k in sorted(set(ud) & set(od)):
        both += 1
        if ud[k] != od[k]:
            bad += 1
            print("≠", k, "oracle=", od[k], "ours=", ud[k])
    print(f"共同键 {both}，其中系数不同 {bad}")


if __name__ == "__main__":
    main()
