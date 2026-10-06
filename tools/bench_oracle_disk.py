"""oracle 落盘版基准：模拟原版 gdbm 一算一存的代价（在仓库根目录运行）。"""
import subprocess, time, statistics, os, re
ORACLE = "./refs/oracle/oracle.exe"  # 自行编译的原版 GROWI（i128 版），见 tools/README.md
TMP = os.environ.get("TEMP", "/tmp")

def pts(n,e,g,d,pts): return f"{n} {e} {g} {d}\n{len(pts)}\n" + "\n".join(f"{h} {p} {m}" for h,p,m in pts)

W = {
 "Kontsevich P^2 d=8..10": [pts(2,1,0,d,[(2,0,0)]*(3*d-1)) for d in (8,9,10)],
 "椭圆 P^2 d=6..8": [pts(2,1,1,d,[(2,0,0)]*(3*d)) for d in (6,7,8)],
 "亏格2 P^2 d=5..6": [pts(2,1,2,d,[(2,0,0)]*(3*d+1)) for d in (5,6)],
 "五次三维流形 d=1..6": [pts(3,5,0,d,[]) for d in (1,2,3,4,5,6)],
 "五重切 e=5..8": [pts(1,e,0,2,[(0,0,2)]*5) for e in (5,6,7,8)],
 "Getzler P^3 g=1": [pts(3,1,1,1,[(3,0,0)]*2), pts(3,1,1,1,[(2,0,0)]*4),
                      pts(3,1,1,5,[(3,0,0)]*10), pts(3,1,1,5,[(3,0,0)]*5+[(2,0,0)]*10),
                      pts(3,1,1,5,[(2,0,0)]*20)],
}

def run(queries, disk=None):
    env = dict(os.environ)
    if disk: env["ORACLE_DB_DISK"] = disk
    p = subprocess.Popen([ORACLE], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, env=env)
    t0 = time.perf_counter()
    for q in queries:
        p.stdin.write(q+"\n"); p.stdin.flush()
        p.stdout.readline()
    dt = time.perf_counter()-t0
    p.stdin.close(); err = p.stderr.read(); p.wait()
    rd=wr=hits=0
    for m in re.finditer(r"reads=(\d+) hits=(\d+) writes=(\d+)", err):
        rd+=int(m.group(1)); hits+=int(m.group(2)); wr+=int(m.group(3))
    return dt, rd, hits, wr

print(f"{'工作负载':<22} {'内存库':>9} {'磁盘库(一算一存)':>12} {'写次数':>7} {'命中率':>6} {'每写代价':>9}")
for name, qs in W.items():
    mem = statistics.median([run(qs)[0] for _ in range(3)])
    disk_ts, rd, hits, wr = [], 0, 0, 0
    for k in range(3):
        d = f"{TMP}/dbtest{k}.bin"
        try: os.remove(d)
        except: pass
        t, r, h, w = run(qs, d); disk_ts.append(t); rd, hits, wr = r, h, w
    disk = statistics.median(disk_ts)
    per = (disk-mem)/wr*1000 if wr else 0
    print(f"{name:<22} {mem*1000:7.1f}ms {disk*1000:10.1f}ms {wr:7d} {hits/max(rd,1)*100:5.0f}% {per:7.3f}ms")
