# -*- coding: utf-8 -*-
# 注意：λ 高次幂（如 λ1^9）是 admcycles 的纯 Python 死区（1200s 超时），按
# 规矩不纳入对拍；我们 2.1s 算 g=5 λ1^12，供参考。
"""sc-faber（CLI）vs admcycles（WSL）速度对拍。只测本机两边都能算的集合。
用法（Windows 侧）：uv run --no-project python tools/faber_bench.py
"""
import json
import subprocess
import time

FABER = "./target/release/faber.exe"

QUERIES = [
    # (名称, faber JSON, admcycles 构造代码)
    ("Kontsevich d=10", '{"op":"klp","g":0,"n":10,"psi":[[1,7]]}',
     "psiclass(1,0,10)**7"),
    ("κ1^3 M̄2", '{"op":"klp","g":2,"n":0,"kappa":[[1,3]]}',
     "kappaclass(1,2,0)**3"),
    ("ψ^4 M̄2,1", '{"op":"klp","g":2,"n":1,"psi":[[1,4]]}',
     "psiclass(1,2,1)**4"),
    ("κ6 M̄3", '{"op":"klp","g":3,"n":0,"kappa":[[6,1]]}',
     "kappaclass(6,3,0)"),
    ("λ3λ2ψ² M̄3,1", '{"op":"klp","g":3,"n":1,"psi":[[1,2]],"lambda":[[3,1],[2,1]]}',
     "lambdaclass(3,3,1)*lambdaclass(2,3,1)*psiclass(1,3,1)**2"),
    ("DR(2)·κψ", '{"op":"dr","g":2,"a":[1,1,-2],"r":2,"kappa":[[1,1]],"psi":[[1,3]]}',
     "DR_cycle(2,(1,1,-2),d=2)*kappaclass(1,2,3)*psiclass(1,2,3)**3"),
]


def bench_faber(q):
    best = None
    for _ in range(3):
        t0 = time.perf_counter()
        r = subprocess.run([FABER, "--batch"], input=q + "\n",
                           capture_output=True, text=True)
        dt = time.perf_counter() - t0
        v = json.loads(r.stdout.strip())["value"]
        best = dt if best is None else min(best, dt)
    return best, v


ADR = r'''
import time, sys
from admcycles import psiclass, kappaclass, lambdaclass, DR_cycle
cases = %s
out = []
for name, expr in cases:
    t0 = time.perf_counter()
    v = eval(expr).evaluate()
    out.append((name, time.perf_counter()-t0, str(v)))
print("RESULT", out)
'''

def main():
    adr_cases = [(n, e) for n, _, e in QUERIES]
    r = subprocess.run(
        ["wsl", "-d", "Debian", "-u", "oracle", "--",
         "/home/oracle/miniforge3/envs/sage/bin/python", "-u", "-c",
         ADR % repr(adr_cases)],
        capture_output=True, text=True, timeout=1200)
    line = [l for l in r.stdout.splitlines() if l.startswith("RESULT")][-1]
    adr = {n: (t, v) for n, t, v in eval(line[len("RESULT "):])}
    print(f"{'查询':<18} {'faber':>9} {'admcycles':>10}  值一致")
    for name, q, _ in QUERIES:
        dt, v = bench_faber(q)
        at, av = adr[name]
        ok = "✓" if v == av else f"✗ ({v} vs {av})"
        print(f"{name:<18} {dt*1000:7.1f}ms {at*1000:8.1f}ms  {ok}")


if __name__ == "__main__":
    main()
