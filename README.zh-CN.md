中文 | [English](README.md)

# stablecurve

[![CI](https://github.com/siqiliu-tsinghua/stablecurve/actions/workflows/ci.yml/badge.svg)](https://github.com/siqiliu-tsinghua/stablecurve/actions/workflows/ci.yml)

三个代数几何计算引擎的 Rust 工作区，共用有理数层（`sc-num`）：

| crate | 内容 | 前身 |
|---|---|---|
| `sc-growi` | Gromov-Witten 不变量：P^n（无固定亏格上限）、超曲面（E > 1 的正亏格只部分支持，与原版相同）、相对不变量 | Andreas Gathmann《GROWI》1.0.4（C++） |
| `sc-faber` | 曲线模空间 M̄_{g,n} 的重言代数：ψ/κ/λ/ch 交数、DR（Pixton/Chiodo/k-twist/spin）、Θ 类、Hain 除子、Hyperell/Biell 配对、Witten r-旋、广义 Hodge chern（GRRcomp） | Carel Faber《Numbers》（Maple）+ admcycles 选定的七项功能 |
| `sc-gwall` | 光滑曲线的 local GW 不变量（后代/环路/ch 插入）与 local BPS（BP 论文表全部 121 个数值格，含 4 格论文勘误）；全部亏格 g ≥ 0，度 d ≥ 2 时受原版 t 截断 ≤ 10 的限制 | Jim Bryan《GWall》（Maple）+ BP 论文 |
| `sc-num` | malachite 有理数包装（三引擎共用） | — |

- **用法（CLI + Python）**：[USER_GUIDE.zh-CN.md](USER_GUIDE.zh-CN.md)
- **分析与设计**：[docs/GROWI_REDESIGN.md](docs/GROWI_REDESIGN.md)、
  [docs/FABER_REDESIGN.md](docs/FABER_REDESIGN.md)、
  [docs/GWALL_REDESIGN.md](docs/GWALL_REDESIGN.md)

**接口。** 每个引擎都是 Rust 库，另有 CLI 与 Python 模块，结果均为精确
有理数（Python 中为 `Fraction`）。`sc-growi` 与 `sc-gwall` 三种接口功能
相同；`sc-faber` 在 Rust 之外只覆盖一部分：

| sc-faber | Rust | CLI | Python |
|---|---|---|---|
| `klp`（κ/ψ/λ/ch 单项式）、`dr`（Pixton DR） | ✓ | ✓ | ✓ |
| `drpoly`、`drspin`、`witten`、`hyperell`、`biell`、`grr` | ✓ | ✓ | — |
| Θ 类、Hain 除子 | ✓ | — | — |

## 快速开始

```bash
uv sync                          # 构建并安装 Python 模块 sc_growi / sc_faber / sc_gwall
cargo build --release            # CLI 二进制 target/release/{growi,faber,gwall}
uv run pytest                    # Python 测试
# Rust 工作区同时构建 pyo3 crate，需要 Python ≥ 3.12：
export PYO3_PYTHON="$PWD/.venv/Scripts/python.exe"   # Windows；其他平台 .venv/bin/python
cargo test --workspace --release # Rust 测试
```

三引擎均已通过各自"黄金数字"验收：GROWI 对照 gwcheck.cm 全部 41 条
断言；faber 对照 Faber《Instructions》、七项选定功能的 admcycles 差分，
以及 admcycles 数据集中 M̄_g（g ≤ 14）上全部 53159 个 λ 单项式 Hodge
积分；gwall 对照 Bryan–Pandharipande 论文 local BPS 表（n₂–n₅）全部
121 个数值格（其中 4 格订正论文印刷错误），以及论文给出的 g = 0、
g = 1 完整结果。详见各设计文档的验收记录。

**性能。** gwcheck 全套在现代硬件上冷缓存启动用时 4.3s。同一台机器上，
计算主导的负载比加强版原版 GROWI 快 1.3–1.5 倍（毫秒级查询以持久缓存的
固定开销为主），且能算该版本会崩溃的亏格 2、度 ≥ 7 与亏格 3；基准见
[docs/GROWI_REDESIGN.md](docs/GROWI_REDESIGN.md)。

**持续集成。** 每次推送都在 GitHub Actions 上测试
（[`ci.yml`](.github/workflows/ci.yml)），覆盖 Linux（x86_64）、macOS
（Apple Silicon）与 Windows（x86_64）：全部 Rust 测试（Rust stable 与最低
支持版本 1.85；含 g ≤ 9 的 2033 个 admcycles Hodge 积分与各独立计算路径
之间的一致性检验）、rustfmt 与 clippy、Python 测试、gwcheck 验收与
growi↔faber 交叉验证。与 admcycles 的 g ≤ 14 全量对拍是一次性运行
`tools/admcycles_hodge.py check --max-genus 14`（约一小时）。v0.1.0 发布前
曾在 Windows 11、Debian 13（rustc 1.85.1）与 macOS 26.7（Intel）上人工
验证。

## 许可证

Copyright (C) 2026 Si-Qi Liu。本项目是自由软件，采用
[GNU 通用公共许可证第 3 版或更新版本](LICENSE)（GPL-3.0-or-later）：
允许使用、修改与再分发，但衍生作品须以同一许可证连同源代码一并分发。

`sc-growi` 译自 GROWI（GPL），`sc-faber` 的部分模块对应 admcycles
（GPL-2.0-or-later）；`sc-faber` 的 κ/ψ/λ/ch 内核移植自 C. Faber 的
《Numbers》，`sc-gwall` 移植自 J. Bryan 的《GWall》，二者均为作者公开
发布的研究程序。上游作品与署名见 [NOTICE](NOTICE)。

## 开发者信息

本项目在 GLM-5.3、GPT-6-Astra 与 Claude Opus 5.5 协助下完成。

## 参考文献

- A. Gathmann, *GROWI 1.0.4*（原程序与 gwcheck.cm）：
  https://personal.math.uni-sb.de/~agathmann/growi/
- C. Faber, *Numbers*（Maple 程序与 Instructions）：
  https://math.stanford.edu/~vakil/programs/Numbers.tar
- J. Bryan, *GWall*（原程序）：
  https://personal.math.ubc.ca/~jbryan/gwall.html ；
  J. Bryan, R. Pandharipande, *BPS states of curves in Calabi–Yau
  threefolds*： https://arxiv.org/abs/math/0009025
- admcycles（Sage 库，faber 的差分 oracle），开发者：Aaron Pixton、
  Johannes Schmitt、Vincent Delecroix、Jason van Zelm、
  Jonathan Zachhuber、Johannes Schwab：
  https://gitlab.com/modulispaces/admcycles
