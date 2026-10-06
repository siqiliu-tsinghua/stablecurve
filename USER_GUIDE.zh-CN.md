中文 | [English](USER_GUIDE.md)

# stablecurve 用户指南

三个引擎各自提供 **CLI 二进制**（批处理友好）与 **Python 模块**（pyo3
直连，`fractions.Fraction` 出值）。本文是完整的用法参考；各引擎的
算法、验证与设计决策见 `docs/` 下的设计文档。

前置：Rust 工具链（≥ 1.85；Debian 13 stable 自带的 `rustc` 即可）与
[uv](https://docs.astral.sh/uv/)。

- CLI 二进制：`cargo build --release` 生成 `target/release/growi`、
  `faber`、`gwall`。
- Python 模块：`uv sync` 构建 `sc_growi`、`sc_faber`、`sc_gwall` 三个
  扩展并以可编辑方式装入 `.venv`；Rust 源码变动后会自动重建。要打独立
  wheel：`uvx maturin build --release -m crates/sc-growi-py/Cargo.toml`
  （`sc-faber-py`、`sc-gwall-py` 同理）。

三个 CLI 的批处理协议一致：输入每行一个查询，输出每行一个 JSON 对象，
`"line"` 为从 1 起的输入行号（空行跳过）。

---

## 1. sc-growi——GW 不变量（GROWI 移植）

计算 P^n、超曲面（含五次 Calabi–Yau 三维流形）与相对条件下任意亏格的
Gromov-Witten 不变量。

### 1.1 表达式语法（原版 `expr.cc` 直译）

- 赋值：`N=3`（环境空间维数）、`E=5`（超曲面次数）、`G=2`（亏格）、
  `D=5`（曲线次数）
- 关键字：`rational`/`elliptic`（=G0/G1）、`lines`/`conics`/`cubics`…
  （=D1/D2/D3…）、`point`/`curve`/`surface`/`threefold`（或 `3fold`）…、
  `quadric`/`cubic`/`quintic`…（=E2/E3/E5…）、`P^2`/`P^3`/…
- 停用词（等价逗号）：`in`/`on`/`to`/`thru`/`through`/`with`/`and`
- 类表达式：`H`、`psi` 的 `+ * ^ ( )` 多项式；`H:9` = H 重复 9 次；
  `[类,重数]` = 相对条件；`pc(k)` = 超曲面的第 k 号原始类
- 恒为零的类使不变量为 0（与原版一致）。数值有范围检查（N/E/G/D ≤ 32767，
  类的总次数 ≤ 255），越界报错，绝不截断。
- 例：`lines in cubic surface`、`elliptic D=6 in P^2 thru H^2:18`、
  `G=2,D=5 in P^4 thru H^3*psi^22`、`conics in P^1 thru H^2:4,[1,2]`

### 1.2 CLI

```bash
./target/release/growi "lines in quintic threefold"        # 2875
./target/release/growi --batch --cache gw.sqlite < q.txt   # 每行一个表达式
```

缓存缺省为 `~/.growi/cache.sqlite`；`--no-cache` 只用内存。批处理应答：

```json
{"line":1,"status":"ok","value":"2875"}
{"line":2,"status":"unknown"}
{"line":3,"status":"invalid","message":"..."}
```

`unknown` 表示当前策略链算不出（原版 UNDEFINED；单查询模式打印
`UNDEFINED (...)`）。

### 1.3 Python

```python
from sc_growi import GROWI
g = GROWI()                                   # 两级缓存，落盘到 ~/.growi/cache.sqlite
g = GROWI(cache="gw.sqlite")                  # 指定缓存文件
g = GROWI(cache="")                           # 只用内存
g.gw("G=2,D=5 in P^4 thru H^3*psi^22")       # Fraction(-41369, 110075314176)
```

超出支持范围（如 E>1 的高亏格超曲面不变量）返回 `None`；输入非法抛
`ValueError`。计算期间释放 GIL。

---

## 2. sc-faber——重言代数积分（Faber Numbers 移植 + admcycles 对标）

计算 M̄_{g,n} 上的重言类积分：κ/ψ/λ/ch 混合单项式、DR cycle
（Pixton/Chiodo、k-twist、spin）、Θ 类、Hain 除子、超椭圆/双椭圆轨迹、
Witten r-旋、广义 Hodge chern。返回值均为**精确有理数**。

### 2.1 CLI（JSONL 协议）

```bash
./target/release/faber '{"op":"klp","g":2,"n":2,"psi":[[1,2]],"lambda":[[1,3]]}'
./target/release/faber --batch --cache fb.sqlite < queries.jsonl
```

应答：`{"line":k,"status":"ok","value":"p/q"}` 或
`{"line":k,"status":"error","error":"..."}`。公共字段：`g`（亏格）、
`kappa`（[[指标,幂],…]）、`psi`（[[腿号,幂],…]）。ψ 腿号须为标记
1..n（同一腿的多项相乘合并）；未知字段名会被拒绝，字段拼错时报错而不是
被静默忽略。

操作（均含 `g`，均返回 ∫类·X）：

| op | 专有字段 | 含义 |
|---|---|---|
| `klp` | `n`；`lambda`/`ch` 同格式 | κ^κ ψ^ψ λ^λ ch^ch 单项式积分 |
| `dr` | `a`（= A 向量，ΣA = 0）、`r`（缺省 g） | Pixton DR 积分 ∫DR(g,A,r)·X |
| `drpoly` | `a`、`k` | Chiodo 插值主路径（= DR_cycle chiodo_coeff=True） |
| `drspin` | `a`、`k` | spin DR（CSS21；A 全奇、k 奇） |
| `witten` | `a3`=[[C分子,C分母,D],…] | Witten r-旋极限类（A_i = C_i·r + D_i；ΣD = 2g−2） |
| `hyperell` / `biell` | `n`（固定标记）、`m`（共轭对） | 超椭圆 / 双椭圆轨迹积分 |
| `grr` | `n`、`t`（ch 次数）、`l`、`dvec`=[d₁..dₙ]、`a3`=[[h, s掩码, v],…] | 广义 Hodge chern ch_t（PRvZ20）；S 按位掩码 Σ2^{i−1} 编码 |

`drpoly`/`drspin` 要求 ΣA = k(2g−2+n)；ΣA 能被 2g−2+n 整除时 `k` 可省。

例：

```json
{"op":"dr","g":2,"a":[1,1,-2],"r":2,"kappa":[[1,1]],"psi":[[1,3]]}   -> 13/144
{"op":"hyperell","g":2,"n":0,"m":0,"kappa":[[3,1]],"psi":[]}          -> 1/1152
{"op":"biell","g":2,"n":0,"m":0,"kappa":[[2,1]],"psi":[]}             -> 1/32
{"op":"grr","g":2,"n":1,"t":1,"l":1,"dvec":[0],"kappa":[[1,1]],"psi":[[1,2]]} -> 5/576
```

### 2.2 Python

```python
from sc_faber import FABER
f = FABER()                                  # 两级缓存，落盘到 ~/.faber/cache.sqlite
f.klp(2, 2, psi=[(1,2)], lam=[(1,3)])        # Fraction(1, 1440)
f.klp(7, 0, lam=[(1,18)])                    # Fraction(32017001, 638512875)
f.dr(2, [1,1,-2], r=2, kappa=[(1,1)], psi=[(1,3)])   # Fraction(13, 144)
```

`FABER(cache="路径")` 使用指定文件，`FABER(cache="")` 不用持久缓存。
持久缓存是**进程全局**的：新建 `FABER` 会替换（或卸下）此前实例挂接的
缓存，`cache_attached` 反映当前状态。

`dr` 等枚举出的稳定图也按 (g, n) 逐层留在内存里，重复调用直接复用。
大的情形可达数十万个图（数百 MB），`f.clear_graph_cache()` 可将其释放。

注：pyo3 层目前只暴露 `klp`/`dr`，其余 op 走 CLI。

---

## 3. sc-gwall——光滑曲线的 local GW 与 BPS（GWall 移植）

计算到光滑曲线 C_g 的（不连通）GW 不变量（点/环路后代、ch 插入）及
local BPS 数。

### 3.1 插入编码

对应原版 Maple 的 v_i 编码。Python 中插入项为元组或列表；CLI 中为 JSON
数组，或等价的对象 `{"t":类型,"k":K,"i":I}`（`ldeg` 用 `{"t":"ldeg","v":V}`）。

| 编码 | 含义 |
|---|---|
| `("t1", K)` | 后代 τ_K(1) |
| `("pt", K)` | 后代 τ_K(pt) |
| `("a", K, I)` / `("b", K, I)` | 环路后代 τ_K(a_I) / τ_K(b_I) |
| `("chO", K)` | ch_K(Rπ*f*O) |
| `("chL", K)` | ch_K(Rπ*f*L)（L 缺省为 2g−2 次即 K） |
| `("ldeg", V)` | 把 L 的次数改为 V（原版 `[[deg]]`） |

环路插入可按任意顺序给出：环路后代反交换，换到规范序（全部 a 环路在
b 环路之前，原程序的要求）所需的符号自动计入。

### 3.2 CLI（JSONL 协议）

```bash
./target/release/gwall '{"op":"gw","g":2,"d":1,"insert":[["t1",0],["pt",1],["pt",0]]}'
./target/release/gwall --batch < queries.jsonl
```

| op | 字段 | 含义 |
|---|---|---|
| `gw` | `g`、`d`、`insert`、可选 `ldeg` | 不变量 ⟨Π insert⟩ |
| `ztil` / `z` | `g`、`d`、`n` | t 级数系数 t^0..t^n（ztil 含常映射，z 去除） |
| `l` | `g`、`d`、`n` | 度变量 log 连通化 L_d（论文 N 的母函数） |
| `bps` | `g`、`d`、`h` | local BPS 数 n_d^h(g)（全量 GW 反演） |

`d` 必填。级数类 op 返回系数数组，如
`{"line":1,"status":"ok","value":["15/2","0","-5/4"]}`。

### 3.3 Python

```python
from sc_gwall import GWALL
w = GWALL()
w.bps(3, 2, 1)                               # Fraction(4) —— 论文表
w.bps(0, 1, 0)                               # Fraction(1) —— 支持 g = 0
w.gw(2, 2, [])                               # Fraction(8) —— Ztil 常数项
w.l(2, 2, 2)                                 # [15/2, 0, -5/4] —— L_d 系数
w.z(2, 1, 4)                                 # 度 1 闭式系数
```

支持全部亏格 g ≥ 0。度 d ≥ 2 时与原版一样受 t 截断 ≤ 10 的限制
（"branch number 10"），超限查询立即抛 `ValueError`。计算期间释放 GIL。

---

## 4. 跨引擎一致性（选读）

faber 与 growi 在共享的数学面（M̄_{0,n} 上的 ψ 积分）互为 oracle：
`growi "G=0,D=0 in P^d thru H:d,psi^a1,..,psi^am"` = `faber klp(0, d+m,
psi=[(1,a1),..])`，34 例全部精确一致（`tools/cross_growi_faber.py`）；
gwall 的 N₀ 闭式与 faber 的 λ³ 积分在三个独立点上咬合。
