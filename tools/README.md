# tools/ —— 开发用脚本

这些脚本是开发期的差分对拍与基准工具，**不是**用户接口，也不随 Python 包
分发。多数依赖本仓库之外、需自行准备的参照环境，放在被 git 忽略的 `refs/`
下；在仓库根目录运行。

| 脚本 | 用途 | 依赖 |
|---|---|---|
| `gwcheck_all.py` | GROWI 原版 gwcheck.cm 全部 41 条断言（已内嵌）对拍 | Python 模块（`uv sync`） |
| `cross_growi_faber.py` | growi ↔ faber 在 M̄_{0,n} ψ 积分上的交叉验证 | Python 模块（`uv sync`） |
| `bench_oracle.py`、`bench_oracle_disk.py` | 与原版 GROWI 的速度对比 | `refs/oracle/oracle.exe`（自行编译的原版，i128 版） |
| `ref_quintic.py` | 五次三维流形亏格 0 的独立递归参考 | `refs/oracle/oracle.exe` |
| `faber_diff.py`、`faber_bench.py` | sc-faber ↔ admcycles 差分与速度对比 | WSL 中装有 Sage + admcycles 的环境；Rust 侧由 `cargo run -p sc-faber --example diff_dump` 导出 |
| `gen_diff.py`、`pergen_oracle.py`、`coeff_compare.py`、`hyb_oracle.py`、`spin_oracle.py` | DR/Θ/spin 的逐生成元系数对拍（Θ 排障期） | admcycles |
| `witten_*.py`、`grr_*.py` | Witten r-旋、GRRcomp 的锚点与差分 | admcycles |

admcycles 环境的搭建记录见 `docs/FABER_REDESIGN.md`（第三部分 oracle 环境）。
`faber_bench.py` 中的 WSL 发行版名与用户名是作者本机的设置，使用前需改成自己的。
