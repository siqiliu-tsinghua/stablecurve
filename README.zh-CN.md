中文 | [English](README.md)

# stablecurve

[![CI](https://github.com/siqiliu-tsinghua/stablecurve/actions/workflows/ci.yml/badge.svg)](https://github.com/siqiliu-tsinghua/stablecurve/actions/workflows/ci.yml)

三个代数几何计算引擎的 Rust 工作区，共用有理数层（`sc-num`）：

| crate | 内容 | 前身 |
|---|---|---|
| `sc-growi` | Gromov-Witten 不变量（P^n 任意亏格、超曲面、相对不变量） | Andreas Gathmann《GROWI》1.0.4（C++） |
| `sc-faber` | 曲线模空间 M̄_{g,n} 的重言代数：ψ/κ/λ/ch 交数、DR（Pixton/Chiodo/k-twist/spin）、Θ 类、Hain 除子、Hyperell/Biell 配对、Witten r-旋、广义 Hodge chern（GRRcomp） | Carel Faber《Numbers》（Maple）+ admcycles 增量 |
| `sc-gwall` | 光滑曲线的 local GW 不变量（后代/环路/ch 插入）与 local BPS（BP 论文表全部 121 个数值格，含 4 格论文勘误；任意亏格 g ≥ 0） | Jim Bryan《GWall》（Maple）+ BP 论文 |
| `sc-num` | malachite 有理数包装（三引擎共用） | — |

- **用法（CLI + Python）**：[USER_GUIDE.zh-CN.md](USER_GUIDE.zh-CN.md)
- **分析与设计**：[docs/GROWI_REDESIGN.md](docs/GROWI_REDESIGN.md)、
  [docs/FABER_REDESIGN.md](docs/FABER_REDESIGN.md)、
  [docs/GWALL_REDESIGN.md](docs/GWALL_REDESIGN.md)

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
断言（冷启动 4.3s，约 40× 原版）；faber 对照 Faber《Instructions》与
admcycles 差分（对标增量一档七项全部完成），以及 admcycles 数据集中
M̄_g（g ≤ 14）上全部 53159 个 λ 单项式 Hodge 积分；gwall 对照 Bryan–
Pandharipande 论文 local BPS 表（n₂–n₅）全部 121 个数值格（其中 4 格
订正论文印刷错误），以及论文给出的 g = 0、g = 1 完整结果。详见各设计
文档的验收记录。

**持续集成。** 每次推送都在 GitHub Actions 上测试
（[`ci.yml`](.github/workflows/ci.yml)），覆盖 Linux（x86_64）、macOS
（Apple Silicon）与 Windows（x86_64）：全部 Rust 测试（Rust stable 与最低
支持版本 1.85）、rustfmt 与 clippy、Python 测试、gwcheck 验收与
growi↔faber 交叉验证。v0.1.0 发布前曾在 Windows 11、Debian 13
（rustc 1.85.1）与 macOS 26.7（Intel）上人工验证。

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
