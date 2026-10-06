"""按 Gathmann (math.AG/9908054) 独立实现的五次曲面亏格 0 递归参考。"""
from fractions import Fraction as F
from math import factorial
from collections import Counter
import subprocess

E, XDIM = 5, 3
# 注：5 次 = 229305888887648 是 GW 不变量正确值（GV 多覆盖 N₅ = n₅+2875/125；
# 枚举数 n₅ = ...625，此前 ...624 为转抄误差）
TRUE = {1: F(2875), 2: F(4876875,8), 3: F(8564575000,27), 4: F(15517926796875,64), 5: F(229305888887648)}

proc = subprocess.Popen(["./refs/oracle/oracle.exe"], stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE, text=True)
p4cache = {}
def absP4(classes, d):
    key = (tuple(sorted(classes)), d)
    if key in p4cache: return p4cache[key]
    inp = [f"4 1 0 {d}", str(len(classes))] + [f"{h} {p} 0" for h,p in classes]
    proc.stdin.write("\n".join(inp) + "\n"); proc.stdin.flush()
    line = proc.stdout.readline().strip()
    v = F(0) if line == "UNDEFINED" else F(line.split()[-1])
    p4cache[key] = v
    return v

def map_pt(classes):
    l = sum(h for h,_ in classes)
    if l != XDIM: return F(0)
    tq = sum(p for _,p in classes)
    m = F(factorial(tq))
    for _,p in classes: m /= factorial(p)
    return F(E) * m

def expand(cnt):
    out = []
    for (b,q), n in cnt.items(): out += [(b,q)] * n
    return out

def submultisets(cnt):
    items = list(cnt.items())
    acc = [({}, dict(cnt))]
    for v, n in items:
        nxt = []
        for S, rest in acc:
            for take in range(0, n+1):
                S2 = dict(S); r2 = dict(rest)
                if take: S2[v] = take
                r2[v] = n - take
                if r2[v] == 0: del r2[v]
                nxt.append((S2, r2))
        acc = nxt
    for S, r in acc:
        yield S, r

BAD = []
relcache = {}
def rel(classes, relpt, m, d):
    key = (tuple(sorted(classes)), relpt, m, d)
    if key in relcache: return relcache[key]
    if m == 0:
        v = absP4(tuple(sorted(classes)) + (relpt,), d)
        relcache[key] = v; return v
    mp = m - 1
    a, psi = relpt
    total = F(0)
    for cof, int_classes, d0, exts in enum_configs(classes, relpt, mp, d):
        val = cof * absQ(int_classes, d0)
        for (scls, sa, sm, sd) in exts:
            if val == 0: break
            val *= rel(scls, (sa, 0), sm, sd)
        total -= val
    if mp == 0:
        if 0 <= a <= XDIM:
            total += E * absP4(tuple(sorted(classes)) + ((a+1, psi),), d)
    else:
        if a >= 0 and a < XDIM:
            total += E * rel(classes, (a+1, psi), mp, d)
        total += mp * rel(classes, (a, psi+1), mp, d)
    relcache[key] = total
    return total

def absQ(classes, d):
    cs = sorted(classes)
    if d == 0: return map_pt(cs)
    if not cs:
        return rel((), (1,0), E*d+1, d) / d
    return rel(tuple(cs[:-1]), cs[-1], E*d+1, d)

def enum_configs(classes, relpt, mp, d):
    a, psi = relpt
    M = Counter(classes)
    out = []
    for d0 in range(0, d):
        rem = d - d0
        need = mp - E * d0
        if need < 1: continue

        def rec2(rem_d, need_m, maxd, specs, left):
            if rem_d == 0:
                if need_m != 0: return
                # 无序性断言：specs 按规范序非严格递减
                keys = [(sp[0], sp[1], tuple(sorted(sp[3]))) for sp in specs]
                assert all(keys[i] >= keys[i+1] for i in range(len(keys)-1)), keys
                int_classes = [(a, psi)] + expand(left) + [(XDIM - sp[1], 0) for sp in specs]
                s = sum(h_ + q_ - 1 for h_, q_ in int_classes)
                if s != 0: BAD.append((classes, relpt, mp, d, int_classes, d0, specs)); return
                cof = F(1)
                for sp in specs: cof *= sp[2]
                cof /= F(E) ** len(specs)
                for k in Counter(specs).values(): cof /= factorial(k)
                for v, n in M.items():
                    denom = factorial(left.get(v, 0))
                    for sp in specs: denom *= factorial(dict(sp[3]).get(v, 0))
                    cof *= F(factorial(n), denom)
                exts = [(expand(dict(sp[3])), sp[1], sp[2], sp[0]) for sp in specs]
                out.append((cof, tuple(sorted(int_classes)), d0, exts))
                return
            for di in range(1, min(rem_d, maxd) + 1):
                for ai in range(0, XDIM + 1):
                    for S, rest in submultisets(left):
                        sdim = sum((b + q - 1) * n for (b, q), n in S.items())
                        mi = E * di + 2 - ai - sdim
                        if mi < 1 or mi > E * di or mi > need_m: continue
                        Sf = frozenset(S.items())
                        key = (di, ai, tuple(sorted(S.items())))
                        if specs and key > (specs[-1][0], specs[-1][1], tuple(sorted(specs[-1][3]))):
                            continue  # 规范序：避免同构多重集重复枚举
                        rec2(rem_d - di, need_m - mi, di, specs + [(di, ai, mi, Sf)], rest)

        rec2(rem, need, rem, [], dict(M))
    return out

if __name__ == "__main__":
    import sys
    ds = [int(x) for x in sys.argv[1:]] or [1,2,3,4,5]
    for d in ds:
        v = absQ((), d)
        t = TRUE.get(d)
        ok = "?" if t is None else ("YES" if v == t else "NO")
        print(f"n_{d} 参考 = {v}  真值匹配: {ok}")
    print(f"P4 叶子数: {len(p4cache)}  rel 状态数: {len(relcache)}")
