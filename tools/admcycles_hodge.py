#!/usr/bin/env python3
"""sc-faber ↔ admcycles 收录的 Hodge 积分数据（无须本地 Sage）。

数据：admcycles/data/hodge_data.sage 的 hodge_number——M̄_g（g = 2..14）上
全部 λ 单项式的积分，共 53159 条（J. Schmitt, Z. Sun, "Computed data for
Hodge integrals and Jacobian locus coefficients", 2026；admcycles 以 GPLv2+
发布）。默认下载固定提交的版本，也可用 --source 指定本地文件。

  uv run python tools/admcycles_hodge.py fixture            # 生成测试夹具（g ≤ 9）
  uv run python tools/admcycles_hodge.py check --max-genus 11
                                                             # 逐条对拍（经 Python 绑定）
"""
import argparse
import re
import sys
import time
import urllib.request
from collections import Counter
from fractions import Fraction
from pathlib import Path

COMMIT = "e1e53753"  # 2026-01-07，hodge_number 扩充到 g ≤ 14
URL = f"https://gitlab.com/modulispaces/admcycles/-/raw/{COMMIT}/admcycles/data/hodge_data.sage"
FIXTURE = Path(__file__).resolve().parent.parent / "crates/sc-faber/tests/data/admcycles_hodge.txt"
LINE = re.compile(r"hodge_number\[\(([\d, ]+)\)\]\s*=\s*(-?[\d/]+)\s*$")


def load(source):
    """→ [(g, λ 指标升序, 值字符串)]，按 (g, 指标) 排序。"""
    if source:
        text = Path(source).read_text(encoding="utf-8")
    else:
        with urllib.request.urlopen(URL, timeout=120) as r:
            text = r.read().decode("utf-8")
    rows = []
    for line in text.splitlines():
        m = LINE.match(line)
        if m:
            idx = sorted(int(x) for x in m.group(1).replace(" ", "").split(",") if x)
            assert sum(idx) % 3 == 0, line
            rows.append((sum(idx) // 3 + 1, idx, m.group(2)))
    rows.sort()
    return rows


def write_fixture(rows, max_genus):
    sel = [r for r in rows if r[0] <= max_genus]
    head = [
        "# admcycles 收录的 Hodge 积分：M̄_g 上 λ 单项式的积分（g ≤ %d，共 %d 条）。" % (max_genus, len(sel)),
        "# 来源：admcycles/data/hodge_data.sage（提交 %s），J. Schmitt, Z. Sun," % COMMIT,
        "#   \"Computed data for Hodge integrals and Jacobian locus coefficients\" (2026)，",
        "#   https://gitlab.com/modulispaces/admcycles，GPLv2+。",
        "# 由 tools/admcycles_hodge.py fixture 生成。每行：g<TAB>λ 指标（升序，逗号分隔）<TAB>积分值",
    ]
    body = ["%d\t%s\t%s" % (g, ",".join(map(str, idx)), v) for g, idx, v in sel]
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text("\n".join(head + body) + "\n", encoding="utf-8", newline="\n")
    print(f"写入 {FIXTURE}（{len(sel)} 条）")


def check(rows, min_genus, max_genus):
    from sc_faber import FABER

    f = FABER(cache="")
    ok = True
    for g in range(min_genus, max_genus + 1):
        sel = [r for r in rows if r[0] == g]
        t0 = time.time()
        bad = []
        for _, idx, v in sel:
            lam = sorted(Counter(idx).items())
            try:
                got = f.klp(g, 0, lam=lam)
            except ValueError as e:
                bad.append((idx, v, f"错误：{e}"))
                continue
            if got != Fraction(v):
                bad.append((idx, v, got))
        print(f"g={g}: {len(sel)} 条，不一致 {len(bad)}（{time.time() - t0:.1f}s）")
        for idx, v, got in bad[:5]:
            print(f"    λ{idx}: admcycles {v}，sc-faber {got}")
        ok &= not bad
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["fixture", "check"])
    ap.add_argument("--source", help="本地 hodge_data.sage（缺省下载提交 %s 的版本）" % COMMIT)
    ap.add_argument("--min-genus", type=int, default=2)
    ap.add_argument("--max-genus", type=int, default=None)
    a = ap.parse_args()
    rows = load(a.source)
    if a.cmd == "fixture":
        write_fixture(rows, a.max_genus or 9)
    else:
        sys.exit(0 if check(rows, a.min_genus, a.max_genus or 10) else 1)


if __name__ == "__main__":
    main()
