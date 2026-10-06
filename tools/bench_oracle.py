"""最终基准：oracle（i128+内存库，原版算法的加强版）vs Rust 移植。
head-to-head：oracle 安全范围；Rust-only：oracle 崩溃区间（gwcheck.cm 真值校验）。"""
import subprocess, json, time, statistics, os, sys
from fractions import Fraction as F

ORACLE = "./refs/oracle/oracle.exe"
RUST = "./target/release/growi.exe"
TMP = os.environ.get("TEMP", "/tmp").replace("\\", "/")

def rm_cache(c):
    for suf in ("", "-wal", "-shm"):
        try: os.remove(c + suf)
        except FileNotFoundError: pass

def oracle_run(queries):
    p = subprocess.Popen([ORACLE], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    t0 = time.perf_counter()
    out = {}
    for q in queries:
        p.stdin.write(q + "\n"); p.stdin.flush()
        while True:
            line = p.stdout.readline()
            if not line: out[q] = None; break
            s = line.strip()
            if s.startswith("RESULT "): out[q] = F(s.split()[-1]); break
            if s == "UNDEFINED": out[q] = None; break
    dt = time.perf_counter() - t0
    p.stdin.close(); p.wait()
    return dt, out

def rust_run(exprs, cache):
    p = subprocess.Popen([RUST, "--batch", "--cache", cache],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    t0 = time.perf_counter()
    out = {}
    for e in exprs:
        p.stdin.write(e + "\n"); p.stdin.flush()
        j = json.loads(p.stdout.readline().strip())
        out[e] = None if j["status"] != "ok" else F(j["value"])
    dt = time.perf_counter() - t0
    p.stdin.close(); p.wait()
    return dt, out

def pts_line(n, e, g, d, pts):
    return f"{n} {e} {g} {d}\n{len(pts)}\n" + "\n".join(f"{h} {p} {m}" for h, p, m in pts)

# ---- head-to-head 工作负载 ----
H2H = []
H2H.append(("启动开销(1平凡查询)", [
    ("rational D=1 in P^2 thru H^2:2", pts_line(2,1,0,1,[(2,0,0)]*2))]))
H2H.append(("Kontsevich P^2 d=8..10", [
    (f"rational D={d} in P^2 thru H^2:{3*d-1}", pts_line(2,1,0,d,[(2,0,0)]*(3*d-1))) for d in (8,9,10)]))
H2H.append(("椭圆 P^2 d=6..8", [
    (f"elliptic D={d} in P^2 thru H^2:{3*d}", pts_line(2,1,1,d,[(2,0,0)]*(3*d))) for d in (6,7,8)]))
H2H.append(("亏格2 P^2 d=5..6", [
    (f"G=2,D={d} in P^2 thru H^2:{3*d+1}", pts_line(2,1,2,d,[(2,0,0)]*(3*d+1))) for d in (5,6)]))
H2H.append(("五次三维流形 d=1..6", [
    (f"D={d} in quintic threefold", pts_line(3,5,0,d,[])) for d in (1,2,3,4,5,6)]))
H2H.append(("五重切二次曲线 e=5..8", [
    (f"N=1,E={e},D=2,[1,2]:5", pts_line(1,e,0,2,[(0,0,2)]*5)) for e in (5,6,7,8)]))
H2H.append(("Getzler P^3 g=1", [
    ("elliptic lines in P^3 thru H^3:2", pts_line(3,1,1,1,[(3,0,0)]*2)),
    ("G=1,D=1 in P^3 thru H^2:4", pts_line(3,1,1,1,[(2,0,0)]*4)),
    ("G=1,D=5 in P^3 thru H^3:10", pts_line(3,1,1,5,[(3,0,0)]*10)),
    ("G=1,D=5 in P^3 thru H^3:5,H^2:10", pts_line(3,1,1,5,[(3,0,0)]*5+[(2,0,0)]*10)),
    ("G=1,D=5 in P^3 thru H^2:20", pts_line(3,1,1,5,[(2,0,0)]*20))]))

# ---- Rust-only（oracle 崩溃区间，gwcheck.cm 真值校验）----
GWCHK = {"G=2,D=7 in P^2 thru H^2:22": "122824720116",
         "G=2,D=8 in P^2 thru H^2:25": "346860150644700",
         "G=3,D=4 in P^2 thru H^2:14": "1",
         "G=3,D=5 in P^2 thru H^2:17": "7915",
         "G=3,D=6 in P^2 thru H^2:20": "34435125",
         "G=3,D=7 in P^2 thru H^2:23": "153796445095",
         "G=3,D=8 in P^2 thru H^2:26": "800457740515775"}
RO = [("亏格2 P^2 d=7..8", [e for e in GWCHK if e.startswith("G=2")]),
      ("亏格3 P^2 d=4..8", [e for e in GWCHK if e.startswith("G=3")])]

print(f"{'工作负载':<24} {'oracle':>10} {'rust冷':>10} {'rust暖':>9} {'加速(冷)':>8}  一致  抽查值")
rows = []
for name, pairs in H2H:
    exprs = [e for e,_ in pairs]; queries = [q for _,q in pairs]
    ots, rts, wts, ovals, rvals = [], [], [], None, None
    for _ in range(3):
        t, ov = oracle_run(queries); ots.append(t); ovals = ov
    for k in range(3):
        c = f"{TMP}/bench_c{k}.sqlite"; rm_cache(c)
        t, rv = rust_run(exprs, c); rts.append(t); rvals = rv
    c = f"{TMP}/bench_warm.sqlite"; rm_cache(c)
    rust_run(exprs, c)
    for _ in range(3):
        t, _ = rust_run(exprs, c); wts.append(t)
    ok = all(ovals[q] is not None and rvals[e] is not None and ovals[q] == rvals[e]
             for (e,q) in pairs)
    om, rm_, wm = statistics.median(ots), statistics.median(rts), statistics.median(wts)
    last = pairs[-1][0]
    rows.append((name, om, rm_, wm, ok))
    print(f"{name:<24} {om:9.3f}s {rm_:9.3f}s {wm:8.3f}s {om/max(rm_,1e-9):7.1f}x  "
          f"{'OK' if ok else 'FAIL':<4}  {last.split()[-1]}={rvals[last]}")
    sys.stdout.flush()

print()
for name, exprs in RO:
    rts, wts, rvals = [], [], None
    for k in range(3):
        c = f"{TMP}/bench_r{k}.sqlite"; rm_cache(c)
        t, rv = rust_run(exprs, c); rts.append(t); rvals = rv
    c = f"{TMP}/bench_rw.sqlite"; rm_cache(c)
    rust_run(exprs, c)
    for _ in range(3):
        t, _ = rust_run(exprs, c); wts.append(t)
    ok = all(rvals[e] is not None and str(rvals[e]) == GWCHK[e] for e in exprs)
    rm_, wm = statistics.median(rts), statistics.median(wts)
    rows.append((name + " [Rust单侧]", None, rm_, wm, ok))
    print(f"{name+' [Rust单侧]':<24} {'崩溃':>10} {rm_:9.3f}s {wm:8.3f}s  {'—':>8}  "
          f"{'gwcheck✓' if ok else 'FAIL':<4}  {exprs[-1].split()[-1]}={rvals[exprs[-1]]}")
    sys.stdout.flush()

ho = sum(r[1] for r in rows if r[1]); hr = sum(r[2] for r in rows if r[1])
print(f"\n对拍部分合计: oracle {ho:.2f}s  rust冷 {hr:.2f}s  加速 {ho/hr:.2f}x")
print("全部一致" if all(r[4] for r in rows) else "存在不一致！")
