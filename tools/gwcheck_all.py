#!/usr/bin/env python3
"""GROWI 验收阶段 1：gwcheck.cm 全部 41 条断言的 parity。

语法映射：GW(D=5,G=0,P^2,H^2:14) → "G=0,D=5 in P^2 thru H^2:14"。
多插入用停用词分隔（thru 等价逗号）。
d=5 定案（codex 仲裁 + GV 复算）：…648 即镜对称 GW 正确值
（N₅ = n₅ + 2875/5³ = ...625 + 23；此前"真值 …624"为转抄误差）。
"""
import sys
import time
from fractions import Fraction

from sc_growi import GROWI

CASES = [
    # —— Kontsevich P²（3）——
    ("G=0,D=5 in P^2 thru H^2:14", "87304"),
    ("G=0,D=10 in P^2 thru H^2:29", "40739017561997799680"),
    ("G=0,D=15 in P^2 thru H^2:44", "150246278745658335777587625061177835520"),
    # —— Getzler 椭圆 P²（3）——
    ("G=1,D=6 in P^2 thru H^2:18", "57435240"),
    ("G=1,D=7 in P^2 thru H^2:21", "60478511040"),
    ("G=1,D=8 in P^2 thru H^2:24", "96212546526096"),
    # —— Getzler 椭圆 P³（6）——
    ("G=1,D=1 in P^3 thru H^3:2", "-1/12"),
    ("G=1,D=1 in P^3 thru H^3 thru H^2:2", "-1/12"),
    ("G=1,D=1 in P^3 thru H^2:4", "-1/6"),
    ("G=1,D=5 in P^3 thru H^3:10", "-147/4"),
    ("G=1,D=5 in P^3 thru H^3:5 thru H^2:10", "-7179606"),
    ("G=1,D=5 in P^3 thru H^2:20", "-1984020394752"),
    # —— EHX g=2 P²（8）——
    ("G=2,D=1 in P^2 thru H^2:4", "0"),
    ("G=2,D=2 in P^2 thru H^2:7", "0"),
    ("G=2,D=3 in P^2 thru H^2:10", "0"),
    ("G=2,D=4 in P^2 thru H^2:13", "27"),
    ("G=2,D=5 in P^2 thru H^2:16", "36855"),
    ("G=2,D=6 in P^2 thru H^2:19", "58444767"),
    ("G=2,D=7 in P^2 thru H^2:22", "122824720116"),
    ("G=2,D=8 in P^2 thru H^2:25", "346860150644700"),
    # —— EHX g=2 P³（7）——
    ("G=2,D=3 in P^3 thru H^3:6", "1/12"),
    ("G=2,D=3 in P^3 thru H^3:5 thru H^2:2", "5/12"),
    ("G=2,D=3 in P^3 thru H^3:4 thru H^2:4", "5/2"),
    ("G=2,D=3 in P^3 thru H^3:3 thru H^2:6", "46/3"),
    ("G=2,D=3 in P^3 thru H^3:2 thru H^2:8", "307/3"),
    ("G=2,D=3 in P^3 thru H^3 thru H^2:10", "747"),
    ("G=2,D=3 in P^3 thru H^2:12", "5930"),
    # —— EHX g=3 P²（8）——
    ("G=3,D=1 in P^2 thru H^2:5", "0"),
    ("G=3,D=2 in P^2 thru H^2:8", "0"),
    ("G=3,D=3 in P^2 thru H^2:11", "0"),
    ("G=3,D=4 in P^2 thru H^2:14", "1"),
    ("G=3,D=5 in P^2 thru H^2:17", "7915"),
    ("G=3,D=6 in P^2 thru H^2:20", "34435125"),
    ("G=3,D=7 in P^2 thru H^2:23", "153796445095"),
    ("G=3,D=8 in P^2 thru H^2:26", "800457740515775"),
    # —— 五次三维流形有理（6）——
    ("G=0,D=1,E=5,N=3", "2875"),
    ("G=0,D=2,E=5,N=3", "4876875/8"),
    ("G=0,D=3,E=5,N=3", "8564575000/27"),
    ("G=0,D=4,E=5,N=3", "15517926796875/64"),
    ("G=0,D=5,E=5,N=3", "229305888887648"),  # 镜对称正确值（GV 多覆盖复算）
    ("G=0,D=6,E=5,N=3", "248249742157695375"),
]


def parse_rational(s: str) -> Fraction:
    s = s.strip().strip('"')
    if "/" in s:
        a, b = s.split("/")
        return Fraction(int(a), int(b))
    return Fraction(int(s))


def main():
    cache = sys.argv[1] if len(sys.argv) > 1 else None
    g = GROWI(cache=cache) if cache else GROWI()
    ok = fail = 0
    t_all = time.time()
    times = []
    for (expr, expected) in CASES:
        t0 = time.time()
        err = None
        try:
            got = g.gw(expr)
        except Exception as e:  # noqa: BLE001
            got = None
            err = str(e)
        dt = time.time() - t0
        times.append((dt, expr))
        want = parse_rational(expected)
        mark = "ok" if got == want else "FAIL"
        if got == want:
            ok += 1
        else:
            fail += 1
        print(f"[{mark:>4}] {dt:7.2f}s  {expr}")
        print(f"          got {got}  want {want}" + (f"  ERR {err}" if err else ""))
    total = time.time() - t_all
    print(f"\n== {ok}/{len(CASES)} 通过，{fail} 失败，总耗时 {total:.1f}s ==")
    times.sort(reverse=True)
    print("最慢 5 条：")
    for (dt, expr) in times[:5]:
        print(f"  {dt:7.2f}s  {expr}")
    if fail:
        raise SystemExit(1)  # 供 CI 判定


if __name__ == "__main__":
    main()
