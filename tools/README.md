# tools/ —— 开发用脚本

这些脚本是开发期的差分对拍与基准工具，**不是**用户接口，也不随 Python 包
分发。多数依赖本仓库之外、需自行准备的参照环境，放在被 git 忽略的 `refs/`
下；在仓库根目录运行。

| 脚本 | 用途 | 依赖 |
|---|---|---|
| `gwcheck_all.py` | GROWI 原版 gwcheck.cm 全部 41 条断言（已内嵌）对拍 | Python 模块（`uv sync`） |
| `cross_growi_faber.py` | growi ↔ faber 在 M̄_{0,n} ψ 积分上的交叉验证 | Python 模块（`uv sync`） |
| `admcycles_hodge.py` | 与 admcycles 收录的 Hodge 积分数据（λ 单项式，g ≤ 14，53159 条）逐条对拍；`fixture` 子命令生成 `crates/sc-faber/tests/data/` 的测试夹具 | Python 模块（`uv sync`）；数据自动下载，或 `--source` 指定本地文件 |
| `bench_oracle.py`、`bench_oracle_disk.py` | 与原版 GROWI 的速度对比 | `refs/oracle/oracle.exe`（自行编译的原版，i128 版） |
| `ref_quintic.py` | 五次三维流形亏格 0 的独立递归参考 | `refs/oracle/oracle.exe` |
| `faber_diff.py`、`faber_bench.py` | sc-faber ↔ admcycles 差分与速度对比 | WSL 中装有 Sage + admcycles 的环境；Rust 侧由 `cargo run -p sc-faber --example diff_dump` 导出 |
| `gen_diff.py`、`pergen_oracle.py`、`coeff_compare.py`、`hyb_oracle.py`、`spin_oracle.py` | DR/Θ/spin 的逐生成元系数对拍（Θ 排障期） | admcycles |
| `witten_*.py`、`grr_*.py` | Witten r-旋、GRRcomp 的锚点与差分 | admcycles |

不需要本地环境的差分测试有两条：`admcycles_hodge.py`（及由它生成、随
`cargo test` 运行的夹具测试 `crates/sc-faber/tests/admcycles_hodge.rs`），
以及 `crates/sc-faber` 各测试中固化的 admcycles 实测值。需要活 oracle 的
脚本（下表其余各项）只在补新对拍时使用；admcycles 环境的搭建记录见
`docs/FABER_REDESIGN.md`（第三部分 oracle 环境）。
`faber_bench.py` 中的 WSL 发行版名与用户名是作者本机的设置，使用前需改成自己的。
