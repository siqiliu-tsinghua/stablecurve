# GROWI 现代化重做方案（Rust）

基于对 `growi-1.0.4`（约 10,200 行 C++，Andreas Gathmann，2003 年首版）源码的逐文件分析。

---

## 第一部分：原程序技术分析

### 1.1 总体架构

```
用户输入（.cm 文件或命令行表达式）
        │
        ▼
cm 预处理器（cm/*.cc）            —— .cm → .cc 源码到源码翻译
        │
        ▼
g++ 运行时编译 + 链接（growi 脚本）—— 链接 libgw / math / database / -lgmp -lgdbm
        │
        ▼
GW(gwlist)  [libgw/gwinv.cc]      —— 输入校验 + 多线性展开
        │
        ▼
compute(inv) [libgw/compute.cc]   —— 递归主循环，策略链分派：
        │   1. simplify  —— 标准公理 + 数据库查询
        │   2. split2    —— 二分量分裂（WDVV / 亏格 0,1 TRR）
        │   3. trr       —— 任意亏格拓扑递归（"负亏格"技巧）
        │   4. virasoro  —— Virasoro 约束 + 线性方程组
        │   5. relsplit  —— 相对不变量分裂
        ▼
math 库（rational/poly/monomial/linsys）+ database 库（gdbm）
```

### 1.2 各模块职责

**math 层（无算法，纯数据结构）**

| 文件 | 内容 |
|---|---|
| `rational.cc` (583 行) | GMP `mpq_t` 的包装。每个有理数是一个堆分配的 `rational_data`（内含 mpq + 引用计数），`rational` 类本身只是裸指针。运算符用宏批量生成（每个二元运算 3 个变体处理 int 混合）。"UNDEFINED" 用空指针表示，除零不报错而是传染性地产生 UNDEFINED。 |
| `monomial.cc` (363 行) | 单项式 = 排序的 (变量,指数) 数组，同样引用计数。变量编码进单个 int（名字 6 bit + 索引移位）。 |
| `poly.cc` (707 行) | 多项式 = **排序单链表**，带头结点（头结点存引用计数）。每加一个单项式就 `new` 一个节点，节点里嵌着引用计数的 rational。提供 `mon_first/mon_next/mon_coeff/mon_vars` 顺序读取 API，游标 `cur_ptr` **存在对象内部**。 |
| `linsys.cc` (215 行) | 稠密高斯消元，增量式加方程，随时报告 矛盾/不定/唯一可解。给 Virasoro 用。 |

**database 层**

- `database.cc` (485 行)：gdbm 键值库的包装。**位级打包的 key**（`set_key_bits` 按位写入整数，省字节），不变量的 key 还做了游程编码（相邻相同点合并，见 `inv.cc:249`）。写一条就 `gdbm_store` 一次。删库靠 `system("rm -rf ...")`。
- 两个库：`GW.dbm`（不变量）、`symfunc.dbm`（Virasoro 系数）。

**libgw 算法层（数学核心）**

- `inv.h/cc`：不变量数据 `inv { xdim, xdeg, cgenus, cdeg, points: pdata[] }`，每个点 `pdata { Hpow, psipow, mult }`（Hpow 为负表示原始上同调类）。点的**排序规范形**就是数据库 key（原始的 hash-consing 思想）。`vdim()` 给出虚维数贡献。
- `gwinv.cc`：`GW()` 入口。校验 N/D/G/E，把多项式插入类按**多线性**展开——一个里程计式循环遍历所有单项式组合，增量维护系数和虚维数。
- `simplify.cc`：降度优先的公理链——度 0 映射到点公理（g=0 多项式系数乘拓扑交积；g=1 用陈类公式 /24）；排序、原始类归约；**数据库查询**；余切线公理；Pⁿ 单点不变量（J-函数 `∏ 1/(H+i)^{n+1}` 的递推）；基本类公理；除子公理（含修正项）；g>1 无点不变量的加点技巧。
- `split2.cc`：主力。`next_pair()` 枚举两类点集的所有分裂——**4 指针双缓冲交换**（把新分裂写到 ar[2]/ar[3] 再换指针，避免来回拷贝），按块移动排序类、增量维护二项式系数和维数。`split2_sum()` 在其外层再枚举度分裂（有维数可行性界 d1..d2，对角类带 1/E 因子，先算低亏格/低度的一侧、为零即跳过）。其上实现 WDVV 重建阶梯：两点单/双 descendant → Kontsevich–Manin 首重建定理 → g=0 TRR → g=1 TRR（带环路项 /24）。
- `trr.cc`：**负亏格不变量**（Gathmann 的技巧，见 `trr.h` 头注）：内部记 g=−k，表示亏格 0 分量链图之和（权 (−1)^节点数），从而任意亏格 TRR 只用二分量分裂即可递推。`trr_reduce()` 在 ψ≥3g−1 的点上展开；`trr_variable()` 把未知不变量留作线性方程组的变量（供 Virasoro 用）。
- `virasoro.cc`：把 L_k 约束 (A)–(E) 项组装成线性方程，`linsys` 求解。仅支持 Pⁿ（Fano）。
- `relinv.cc` (1,177 行，全库最大) + `relsplit.cc`：相对理论的分裂公式枚举器。`next()` 嵌套枚举：度分配 → 类分组 → 节点数 → 节点重数 → 节点类，自动去重对称构型并把组合因子和重数因子折进 `cof()`。是全程序最难懂、最难移植的部分。
- `symfunc.cc`：Virasoro 需要的初等对称函数系数，同样入库缓存。

**cm 前端**

- `cm.cc` + `expr.cc` + `inout.cc`：手写词法/语法分析，做四类翻译：`a/b → divconv(a)/b`（有理数除法）；`a^b → pow(a,b)`（带优先级爬升的正确加括号）；`GW(x,y,z) → GW(gwlist(x)+gwlist(y)+...)`；自然语言关键字表（lines/conics/cubic/threefold/stopword in/on/thru…，`expr.cc` 的 op 表）。
- `growi` 脚本：翻译 → 调 g++ 编译 → 链接 → 执行。**每次运行用户程序都要起一次完整 C++ 编译**。

### 1.3 "时代特色"低效/危险做法清单

这是重做要清除的东西（括号内是现代替代）：

1. **万物皆引用计数堆对象**：rational/monomial/poly 全是"句柄→堆数据+引用计数"，修改时写时复制（`rational.cc` 的 `RATIONAL_OP1` 宏）。2003 年为省内存；如今指针追踪 + 引用计数流量 + 分配器压力彻底摧毁缓存局部性。（→ 值语义 + 小向量）
2. **多项式是逐节点 `new` 的链表**（`poly.cc:42 paddsub`）：每次加法每个单项式一次堆分配，系数又是嵌套的堆对象。本程序的类多项式实际上只有 2–3 个变量、几个单项式。（→ 排序 `Vec<(系数,单项式)>`，快 1–2 个数量级）
3. **游标存在多项式对象内部，且拷贝时"翻倍"共享**（`gwinv.h:93` 的注释自述）：靠 gwentry 拷贝构造函数复制 cur_ptr 实现两份拷贝独立迭代——极聪明的别名戏法，也极脆弱。（→ `iter()`）
4. **`next_pair` 的 4 指针双缓冲交换**（`split2.cc:92`）：为避免 1:1 回拷而设计的 API，正确性依赖调用方严格约。（→ 迭代器拥有自己的 scratch 空间）
5. **内层递归每次查磁盘数据库**：`simplify` 对每个不变量先做 gdbm 查询，`split2` 每次成功后再 `gdbm_store`——最热路径上是同步磁盘 I/O。（→ 进程内哈希表为主，持久化为辅）
6. **位打包 key + 游程编码 + gdbm**：字节抠到 bit 级（`database.cc set_key_bits`），gdbm 已是遗产依赖（最后一次大更新 2011 年前后，Windows 支持差）；删库用 `system("rm -rf")`，不可移植且有注入风险。（→ 结构化二进制格式 / redb）
7. **用"翻译成 C++ + 现场编译"实现一个小 DSL**：.cm 实际是一门极小的 DSL（查询表达式关键字 + GW 调用），但实现是五段流水线：写 .cm → cm 预处理成 .cc → 调 g++ 编译链接 → 运行 → 清理临时文件。每次运行付一次完整编译，还要求用户机器上装有编译器，临时文件管理也烦人。这是 DSL 技术不发达年代的凑合做法。（→ 查询表达式由引擎二进制直接解析；用户侧程序用 Python 拼 CLI，见 2.4 节）
8. **单线程**：分裂项之间天然独立却串行计算。（→ rayon 并行）
9. UNDEFINED 用空指针+cout 报错，错误沿运算传染。（→ `Result<Q, ComputeError>`）

### 1.4 必须原样保留的聪明设计

这些是程序的价值所在，重做时数学语义必须逐条对齐：

- **策略链顺序**：便宜公理在前（simplify 的公理是有序的，度 0 → 排序 → 查库 → 余切线 → …），失败才走贵的分裂/Virasoro。
- **排序规范形 + 入库缓存**：不变量排序后作为缓存键（原始的 hash-consing）。
- **负亏格链不变量**技巧（`trr.h` 长注释）。
- **增量维护**：枚举分裂时只做相对更新（维数 dw、二项式系数 coeff），配合维数可行性界 d1/d2 提前剪枝、"先算小侧、为零即跳过"的求值顺序启发式。
- **多线性展开**：多项式插入类逐单项式笛卡尔展开（`gwinv.cc` 的里程计循环——语义对，实现换成迭代器）。
- **relinv 的对称因子簿记**：组合数学正确性核心，重写时逻辑必须等价。
- 数学依据：Gathmann, Duke Math. J. 115 (2002) 与 math.AG/0305361。

---

## 第二部分：Rust 重做方案

### 2.0 设计原则

1. **算法忠实**：策略链、公式、规范形语义与原版逐条等价；用原版数字做黄金测试。
2. **数据结构换代**：链表→数组，堆句柄→值语义，磁盘缓存→内存 cache + 可选持久化。
3. **库优先**：核心是 crate，CLI/解析器/脚本都只是它的前端；不再运行时编译任何东西。
4. **先单线程正确，再并行**。

### 2.1 依赖选型

| 用途 | 原版 | 新选择 | 理由 |
|---|---|---|---|
| 大有理数 | GMP (`mpq_t`) | **malachite**（`malachite-q`） | 纯 Rust、活跃维护、对小额数值有优化；留 `rug`（GMP 绑定）作 feature 开关，要极限性能时一键切换 |
| 持久缓存 | gdbm | **SQLite**（`rusqlite` + `bundled` 特性，随包编译、零系统依赖） | 跨会话累积中间不变量是原版的核心价值，必须保留并强化；SQLite 成熟、单文件、WAL 并发读、Windows 原生支持，还能用 sqlite3 CLI 直接检视积累的库 |
| 表达式解析 | 手写 cm 预处理器 | **winnow**（在引擎二进制内解析查询表达式） | 现代、零拷贝、错误信息好；只解析"那一门小 DSL"，不再做源码翻译 |
| 用户脚本 | .cm + 运行时 g++ 编译 | **Python + pyo3 原生模块**（`sc-growi-py` crate，maturin 打包成 `import growi`） | 进程内直连，无子进程/序列化；原版 DSL 承担的循环/算术/打印由成熟脚本生态承接，`fractions.Fraction` 精确承接有理数结果 |
| 并行 | 无 | **rayon**（后期，引擎内）+ Python 多进程分片 | 分裂项 fan-out 天然适合 rayon；跨机/跨批次分片用多进程 |
| CLI | growi shell 脚本 | **clap**（引擎二进制） | 单发 / 批量 / 常驻三种模式 + 缓存子命令 |

### 2.2 核心数据结构（最重要的架构决策）

**hash-consing 不变量仓**是进程内那一级的支点，取代原版的"排序+位打包"编码，并与磁盘 SQLite 组成两级缓存（见 2.5 节）：

```rust
// 每个标记点 6~12 字节，纯值语义
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Pdata { h: i16, psi: i16, mult: i16 }   // h<0 表示第 |h| 号原始类

// 不变量：pts 恒为排序规范形
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Inv {
    dim: u8,        // xdim
    deg: u16,       // xdeg
    genus: i8,      // cgenus，负数 = 链不变量 g=-k
    cdeg: u16,
    pts: Box<[Pdata]>,
}

// 全局唯一化：canonical Inv -> InvId(u32)
// 相等性变成整数比较；哈希在插入时算好存槽里
pub struct InvStore { map: HashMap<Inv, InvId>, slots: Vec<Inv>, }

// 两级缓存：L1 内存（本进程热数据），L2 SQLite（跨会话中间结果库）
pub struct Cache {
    l1: HashMap<InvId, Q>,     // 命中即返回
    dirty: VecDeque<InvId>,    // 新算出、待批量落盘
    db: InvDb,                 // rusqlite 句柄，WAL 模式
}
// trait 抽象，后续并行/换实现不动算法层
pub trait InvariantCache { fn get(&mut self, inv: &Inv) -> Option<Q>; fn put(&mut self, inv: &Inv, v: &Q); }
```

递归全程传 `InvId`：比较 O(1)、克隆免费、缓存键就是 id。原版每次 `compute()` 的"gdbm 磁盘往返"变成一次内存哈希查询——这是单项收益最大的改动。

**多项式**只存在于用户边界（类是 H、ψ、Pᵢ 的多项式，实际最多几个单项式）：

```rust
pub struct Poly(Vec<Term>);                 // 排序的平坦数组
struct Term { c: Q, m: Mono }
type Mono = SmallVec<[(VarId, u8); 4]>;     // smallvec，栈上常见情形零分配
```

GW() 的多线性展开改成对各多项式 `iter()` 的笛卡尔积迭代器，配 `itertools::multi_cartesian_product` 语义的生成器——彻底删除游标别名戏法。

**有理数**：`type Q = malachite::Rational` 薄封装（便于以后切 rug 或加小额快速路径）。中高次数的不变量数值会到 10³⁸⁺，必须任意精度，没有捷径。

**计算主循环**：

```rust
pub fn compute(ctx: &mut Ctx, id: InvId) -> Result<Q, Unknown> {
    if let Some(v) = ctx.cache.get(id) { return Ok(v); }
    let r = simplify::apply(ctx, id)
        .or_else(|_| split2::apply(ctx, id))
        .or_else(|_| trr::apply(ctx, id))
        .or_else(|_| virasoro::apply(ctx, id))
        .or_else(|_| relsplit::apply(ctx, id))?;
    ctx.cache.insert(id, r.clone());
    Ok(r)
}
```

UNDEFINED 传染变成 `Result` 短路；`Ctx` 收编所有原全局状态（N/D/G/E 缺省值、两个缓存）。

**分裂枚举器**：`next_pair` 重写为拥有 scratch 空间的迭代器：

```rust
pub struct PairSplitter { /* 两块 scratch: Box<[Pdata]> + 游标 */ }
impl Iterator for PairSplitter {
    type Item = (InvId, InvId, Q);   // 两侧 + 增量二项式系数
}
```
维数界剪枝（d1/d2）、对角类 1/E、先算小侧等启发式全部保留，只是不再靠双缓冲指针交换。

**线性方程组**：稠密高斯消元照搬（系统很小），`Vec<Vec<Q>>` + `HashMap<i64, usize>` 变量名表。

### 2.3 模块映射表

| 原文件 | 新位置（crate `sc-growi`）| 改动性质 |
|---|---|---|
| `math/rational.cc` | 删除，`num.rs` 用 malachite | 整层退役 |
| `math/poly.cc`、`monomial.cc`、`var.cc` | `poly.rs`（~200 行）| 链表→Vec |
| `math/linsys.cc` | `linsys.rs` | 直译 + Result 化 |
| `database/database.cc` | `cache.rs`（L1 内存）+ `invdb.rs`（L2 SQLite）| 重写：语义保留（算过即落盘、跨会话复用），机制换代（批量事务代替逐条同步写）|
| `libgw/inv.cc` | `inv.rs` + `intern.rs` | 规范形保留，序列化重写 |
| `libgw/gwinv.cc` | `gw.rs` | 校验保留，里程计循环→笛卡尔积迭代器 |
| `libgw/simplify.cc` | `simplify.rs` | 公式直译（全部公式有注释出处）|
| `libgw/split2.cc` | `split2.rs` | `next_pair`→迭代器，数学不变 |
| `libgw/trr.cc` | `trr.rs` | 负亏格技巧原样保留 |
| `libgw/virasoro.cc` | `virasoro.rs` | 直译，仅 Pⁿ 限制照旧 |
| `libgw/symfunc.cc` | `symfunc.rs` | 直译，缓存改内存 |
| `libgw/relinv.cc`、`relsplit.cc` | `relinv.rs`、`relsplit.rs` | **最高风险**：组合枚举逐函数对照移植 |
| `cm/*` + growi 脚本 | `sc-growi-cli` 内的查询解析器（winnow，即 `expr.cc` 关键字表的直译，CLI 与 pyo3 模块共用）+ `sc-growi-py`（pyo3 绑定）| 整层退役：DSL 不再做源码翻译，查询表达式由引擎直接解析；用户"程序"就是普通 Python 脚本 |

### 2.4 前端：三种接入面——pyo3 原生模块（主）/ CLI（辅）/ Rust crate（彻底消灭源码翻译）

原版前端是一条五段流水线：写 .cm → 预处理成 .cc → 调 g++ 编译链接 → 运行 → 清理临时文件。它实际设计的是一门很小的 DSL，只是在 DSL 工具不发达的年代用"生成 C++ 再现场编译"来实现。重做把这件事拆成三种接入面，各归各位：

**① pyo3 原生模块——Python 的主接入方式。** `sc-growi-py` crate 用 pyo3 包装 `sc-growi`，maturin 打包成 pip 可装的 `growi` 模块，Python 进程内直连引擎，无子进程、无序列化层。原版 test.cm 的椭圆空间曲线示例变成：

```python
# elliptic_space_curves.py
from growi import GROWI          # pyo3 原生模块

g = GROWI()                      # 打开 ~/.growi/cache.sqlite（两级缓存）
for d in range(1, 9):
    n0 = g.gw(f"rational D={d} in P^3 thru H^3:{2*d}")      # → Fraction
    n1 = g.gw(f"elliptic D={d} in P^3 thru H^3:{2*d}")
    print(d, n1 + (2*d - 1) * n0 / 12)
```

绑定层的关键设计：

- **语义映射**：算出 → `fractions.Fraction`（大数以字符串过界，Python 侧精确构造）；算不出 → `None`（原版 UNDEFINED 的直译）；输入非法 → 抛 `ValueError`——取代原版"UNDEFINED 空指针传染 + cout 报错"。
- **GIL 释放**：`compute` 是长计算，绑定层用 `py.allow_threads(...)` 把 GIL 放掉再进引擎——计算期间 Python 其他线程照常跑；内部并行届时用 rayon，不受 GIL 约束。
- **线程模型**：`GROWI` 对象内部状态（Ctx/缓存）以 `Mutex` 包裹（rusqlite 连接 `Send` 不 `Sync`，正合适）；同进程多线程真并行用"每线程各开一个 `GROWI` 实例 + WAL 多连接"或进程内 rayon，跨机器/跨批次用多进程分片。
- **边界纪律**：pyo3 类型只出现在 `sc-growi-py`，`sc-growi` 保持纯 Rust——CLI 与模块共享同一实现，不会长出第二套语义。

**② CLI——留给命令行场景自己玩。** 单发 `growi "lines in quintic threefold"` → `2875`；批量 `growi --batch`（stdin 每行一个表达式，stdout 逐行 JSONL）；缓存运维 `growi cache stats / merge / vacuum / export`。语法即 `expr.cc` 关键字表的直译（`P^n` / `lines` / `cubic` / `thru` / `H^2:9` / `[class,mult]:n`），winnow 实现——这套解析器 CLI 与 pyo3 模块共用同一个（`g.gw("...")` 吃的就是同样的字符串）。原版 `.cm` 里 `1/2` 是有理数除法、`^` 是幂，在这里是自然语义，无需任何翻译技巧。机器语义同上：JSON 输出 `{"status":"ok","value":"87304"}` / `"unknown"` / `"invalid"`，退出码 0 = ok/unknown、2 = 非法。批量/行协议服务模式留给非 Python 的编排（shell 管道等）。

**③ Rust crate。** `sc-growi` 照常发布，愿意直接写 Rust 的研究者照用。

**多进程并行**：Python 起 N 个分片（按度数切）依然是最省心的横向扩展。SQLite WAL 天然支持多进程并发读 + 串行写，配 `busy_timeout`；不想争用写锁的长批次，各分片写临时库后 `growi cache merge`（ATTACH + INSERT OR IGNORE）合并。

**项目结构**：

```
stablecurve/             # cargo workspace 根 + uv 工作区根（pyproject.toml）
├── crates/
│   ├── sc-num/          # 有理数层（三引擎共用）
│   ├── sc-growi/        # 算法库（纯 Rust，不含 pyo3）——crate 命名规范 sc-*
│   ├── sc-growi-cli/    # 引擎二进制 growi：单查询 + JSONL 批处理
│   ├── sc-growi-py/     # pyo3 绑定（cdylib）+ python/sc_growi，maturin 打包为 `import sc_growi`
│   └── sc-faber*/、sc-gwall*/   # 另两个引擎，同样三件套
├── tests/               # pytest（三个 Python 模块）
└── examples/            # 示例脚本
```

REPL 不单做：`python -i` 配原生模块天然就是交互环境。

### 2.5 持久化：磁盘 SQLite + 内存的二级缓存（核心特性，非可选）

原版 GW.dbm 的本质是**跨会话的中间结果知识库**：算一个高次/高亏格不变量会递归产生成千上万个更低阶的中间不变量，而不同计算任务之间的中间结果大量重叠——上次算完的库让这次计算近乎免费。原版 README 明说库文件"可以变得非常大"，且删除后"下次运行重建"。这个语义必须原样保留；要改的只是机制：原版**每算出一个值就同步写一次盘**（gdbm_store 逐条、写放大严重），重做改为**读穿 + 批量写回**。

```
compute(id)
  ├─ L1（内存 HashMap）命中 ──────────────► 返回
  ├─ L1 未命中 → SELECT L2（SQLite）─► 命中：填入 L1，返回
  └─ L2 未命中 → 递归计算 ─────────► 写 L1 + 记入 dirty
dirty 队列：每 N 条（如 4096）/ T 秒 / 进程退出时，单事务批量 INSERT
```

**Schema**（`WITHOUT ROWID` 让 BLOB 主键直接做 B 树键，正对键值负载）：

```sql
CREATE TABLE meta (k TEXT PRIMARY KEY, v TEXT);   -- 格式版本 + 算法版本
CREATE TABLE invariants (
  key BLOB PRIMARY KEY,   -- 规范形 Inv 的版本化紧凑编码（LEB128 变长整数）
  num BLOB NOT NULL,      -- 分子（大端字节）
  den BLOB NOT NULL       -- 分母
) WITHOUT ROWID;
-- symfunc 同构，另建一张表
```

**关键决策：**

- **WAL + `synchronous=NORMAL`**：缓存数据可再生，无需逐写 fsync，批量事务吞吐高两三个数量级；崩溃最多丢最后一批 dirty，语义与原版一致（丢了重算即可）。
- **启动版本校验**：meta 表存格式版本 + 算法版本，不符则整库作废重建——取代原版 CHECK_TIME 时间戳技巧，语义更准（算法变了旧值才作废，单纯过期不作废）。
- **不做启动全量预载**：库可能很大，按 miss 逐条读穿即可；L1 每项几十字节、无上限，真到内存压力再加 LRU。
- **key 编码用 LEB128 变长整数**：保留原版"小值占字节少"的洞察，但用标准编码而非私有位打包，跨平台无歧义。
- **并行阶段预留**：WAL 天然支持"多读者 + 单写者"——每工作线程一个只读连接，专用写线程消费 dirty 队列；接口现在就抽成 trait（`InvariantCache`），算法层不感知。

### 2.6 测试策略

1. **黄金数字**：`gwcheck.cm` 全部断言直接移植——Kontsevich 数（P² d≤15）、Getzler 椭圆（P²/P³）、EHX/Caporaso–Harris g=2,3、五次三维流形镜对称序列（2875, 4876875/8, …）。
2. **差分测试**：✅ 已建成（`refs/oracle/`，未入库）——原版算法源码 + MSVC 编译 + 精确 128 位有理数桩 + 内存缓存桩，绕开 GMP/gdbm 依赖；stdin 喂不变量、stdout 出值。M1 期间用它仲裁了三起疑案（见下）。
3. **性质测试**（proptest）：维数条件守恒、除子公理恒等式、结果与插入类排序无关性。
4. **回归对照**：`SC_GROWI_TERMS`（split2 项级转储）、`SC_GROWI_REL`（relsplit 候选转储）、`SC_GROWI_DEBUG`（compute 计数）——env 开关，零开销。

**M1 差分发现（重要）**：
- ✅ **已定案（2026-10-06，见下方"definitive"段）**：原版 quintic d=5 输出 2293058888876**48** 是镜对称正确值（GV 多覆盖 N₅ = n₅ + 2875/125；此前"真值 …624"为本方转抄误差——枚举数 ...625 误记 ...624）。原版无缺陷，移植如实复现即为正确。
- README 的相对不变量示例 `GW(N=1,E=1,D=2,G=0,H^2:4,[1,2])`：原版实际输出 **2**（经典特征数也是 2），README 写 1 系笔误。

**d=5 +24 法医进展（M5 期间第一轮排查，已排除的假设）**：
1. ❌ 伪候选：oracle 切除实验（顶层 + 全层级各跳过第 N 个候选，N=0..200）无一例给出真值 …624——不存在"多算了一个整体配置"。
2. ❌ 缓存中间量被污染：d=5 运行写入的全部 300 个缓存键中，仅最终答案本身为 d=5 独有；其余全部与（正确的）d≤4 或 d=6 运行共享，且 d=6 对这些量的依赖系数非零 ⇒ 共享值全对。
3. ❌ 空间切换后的 1 点 P⁴ 值：用镜对称闭式 `⟨H^hψ^p⟩_{0,d}(P⁴) = [H^{4−h}]∏_{m=1}^d (H+m)^{−5}` 独立验证 d=1..6 全部命中（含 d=5 的 −137/298598400000、16519/5971968000000）。
4. ❌ 度 0 内部分量：<H³>₀=5（=∫_Y H³）等全部正确。
5. ❌ 顶层主项系数（5、25）与 rs_genus0_abs 的除子补偿（/d）核对无误。
- **剩余嫌疑**：relinv 的 `cof()` 组合簿记在某个 d=5 特有配置上偏差恰为 −24（=4!，或 /5 后 −120=5!）——重点查 same 值跟踪（same(i,ncls+1)）与原始类因子（prims=204 的 `(prims+i)` 连乘与 `factorial(cpt_prim[i]−j)`）。

**d=5 定案（2026-10-06，codex 终局仲裁 + 本方独立复算）——原版无缺陷，+24 是我方转抄误差**：
1. GW 不变量 N₅ = 2293058888876**48** 是镜对称正确值。GV 多覆盖公式
   N_d = Σ_{k|d} n_{d/k}/k³ 给 N₅ = n₅ + 2875/5³ = ...62**25** + 23 = ...648
   （本方对 d=1..6 六值全部独立复算自洽；d=5 处 2875/125 = 23）。
2. 文献佐证（codex 联网检索）：Dang (arXiv:1409.3712 §5，localization
   独立计算)、Cooper–Zinger (Michigan Math. J. 63 (2014) Table 1)、
   Givental 镜像公式级数展开（codex 用精确有理数算至 d=6）三处独立
   均为 ...648；未找到任何支持 ...624 的文献或 erratum。
3. "+24"分解：648 − 624 = 23 + 1——23 恰为多覆盖项 2875/125；余 1 是
   **我方此前的转抄误差**（COGP 枚举数 n₅ = ...625 误记 ...624）。
   上文"镜对称真值 …624"的所有记载均源于此误记；gwcheck.cm 断言与
   独立数学一致，1.0.3/1.0.4/编译补丁/相对点更新全部无嫌疑。
4. 处置：移植值 ...648 保持（无任何 override）；golden.rs / 工具脚本
   / 本文档案的 624 引用全部更正。法医第一轮的五个"已排除假设"
   结论仍有效（它们本来就没有发现任何异常）。
- **下一步方案**：按 Gathmann Duke 论文 §4 的分裂系数公式**独立实现**一份 `cof_reference()`，对 d=5 全部候选逐个对拍 cof；这是比继续切除实验更强的定位手段。

### 2.7 里程碑

| 里程碑 | 内容 | 验收 |
|---|---|---|
| M0 | workspace + `Q/poly/inv/intern` + 单测 | 单元测试过 |
| M1 | simplify + split2 + **两级缓存（含 SQLite L2）** + 最小 CLI（`growi "expr"`）+ Pⁿ 亏格 0 全链 | `growi "lines in quintic threefold"` 输出 2875、`quintic d=1..6`、Kontsevich d≤10；重跑同题秒回（库命中）|
| M2 | g=1（Getzler）+ 负亏格 TRR → 任意亏格 | gwcheck 的 g=1,2,3 全部 |
| M3 | virasoro + linsys | 高亏格 P³ 系数与原版对拍 |
| M4 | **relsplit/relinv**（相对不变量 + 超曲面 E>1）| README 例：`conics thru H^2:4,[1,2]`、椭圆空间曲线计数表 |
| M5 | 查询表达式全量（关键字表全覆盖）+ pyo3 模块（maturin，GIL 释放）+ CLI 批量模式 + 缓存运维（stats/merge/VACUUM）| `pip install` 后 `import growi` 可用；README 与 test.cm 全部示例以 Python 脚本复现 |
| M6（可选）| rayon 并行、rug 后端开关、wasm demo | 基准：gwcheck 用时对比原版 |

### 2.8 风险与对策

| 风险 | 对策 |
|---|---|
| relinv 组合枚举移植出错（1,177 行状态机）| 逐函数对照 + 差分测试为主力；放 M4 单独成期 |
| 原版本身有 bug（NEWS 提到 2006 年才修的编译期错误）| 以数学文献和 gwcheck 数字为准，不以原版输出为准 |
| 大数性能不如 GMP | malachite 先行，`rug` feature 兜底；基准在 M1 就建好 |
| Virasoro 猜想假设（仅 Pⁿ）| 照原版限制，文档写明 |
| 缓存内存膨胀 | L1 每项几十字节且 L2 已持久化，可主动收缩 L1（换 LRU）而不丢数据；SQLite 单文件可 VACUUM 压实 |

### 2.9 性能预期

- 首次计算：内层"逐条同步写盘"→"内存哈希 + 批量事务写回"，缓存密集的计算（gwcheck 这类）预期 1–2 个数量级提升；链表多项式 → 平坦数组再给边界代码 10–100×。
- **再次计算同类问题：中间不变量全部从 SQLite 命中，近乎免费**——这是对原版核心价值（跨会话中间结果复用）的保留与放大。
- 前端：没有预处理/编译/清理环节；Python 走 pyo3 进程内直连，零子进程、零序列化开销；CLI 留给命令行场景；多进程分片可进一步扩展。
- rayon 分裂项并行：额外线性于核数的空间（WAL 多读者单写者结构与 dirty 队列天然衔接）。
- 单纯移植而不做并行，目标就已是从"3 分钟 gwcheck（2003 年硬件）"到秒级。

## 性能基准（2026-10，M5 后与原版对拍）

对照物：`refs/oracle`（原版 1.0.4 算法逐字拷贝 + MSVC /O2 + i128 精确有理数 + 进程内内存库）——
即**原版的加强版**：比真机 GROWI（GMP 大分数 + gdbm 磁盘库）更快。Rust 若不输它，即全面胜过原版。
协议：每工作负载固定顺序批量查询；oracle 每次重复开新进程（恒冷），Rust 每次重复换全新 SQLite 缓存文件（冷），
另测同缓存二跑（暖）；各 3 次取中位。脚本 `tools/bench_oracle.py`。

| 工作负载 | oracle 冷 | Rust 冷 | Rust 暖 | 冷加速 |
|---|---|---|---|---|
| 启动开销（1 平凡查询） | 2 ms | 10 ms | 4 ms | — |
| Kontsevich P² d=8..10 | 2 ms | 10 ms | 4 ms | — |
| 椭圆 P² d=6..8 | 10 ms | 16 ms | 4 ms | 0.7× |
| 亏格 2 P² d=5..6 | 185 ms | 125 ms | 4 ms | 1.5× |
| 五次三维流形 d=1..6 | 5 ms | 11 ms | 5 ms | 0.5× |
| 五重切二次曲线 e=5..8 | 11 ms | 18 ms | 4 ms | 0.6× |
| Getzler P³ g=1（5 查询） | 37 ms | 28 ms | 4 ms | 1.3× |
| 亏格 2 P² d=7..8 | 崩溃 | 192 ms | 5 ms | — |
| 亏格 3 P² d=4..8 | 崩溃 | 3.83 s | 4 ms | — |

补充实验（同一 oracle，`ORACLE_DB_DISK` 开关模拟原版 gdbm 一算一存：每条记录追加落盘 + `_commit`
物理刷盘；读仍走内存，等价 gdbm 桶缓存；2026 年 NVMe 实测）：

| 工作负载 | 内存库（对拍用） | 磁盘库（一算一存） | 写次数 | 命中率 | 放大 |
|---|---|---|---|---|---|
| Kontsevich P² d=8..10 | 3 ms | 50 ms | 66 | 86% | 16× |
| 椭圆 P² d=6..8 | 12 ms | 376 ms | 680 | 92% | 32× |
| 亏格 2 P² d=5..6 | 194 ms | 3457 ms | 9366 | 95% | 18× |
| 五次三维流形 d=1..6 | 7 ms | 188 ms | 443 | 73% | 29× |
| 五重切二次曲线 e=5..8 | 15 ms | 526 ms | 1324 | 70% | 36× |
| Getzler P³ g=1 | 45 ms | 835 ms | 2288 | 95% | 19× |

即：**同机同算法，仅"每写一条刷一次盘"这一项就把原版拖慢 16–36 倍**（每写 0.35–0.7 ms）。
真机原版还要叠加：gdbm 随机桶写入比顺序追加更伤、GMP 大分数比 i128 慢、2003 年 HDD 每次同步
5–15 ms（再乘 15–40×）。这与"2003 年跑一遍 gwcheck 约 3 分钟"的量级完全吻合。
Rust 的 L2 用 WAL + synchronous=NORMAL（成组提交、检查点时刷盘）——保住持久性的同时
免掉每写一刷，这正是冷启动仍能快 1–2 个数量级的设计原因。70–95% 的读命中率也印证了
"不同不变量共享中间数据"这一持久化核心前提的量级。

结论：
- **所有对拍数值逐查询精确一致**；Rust 单侧两行与 gwcheck.cm 文献真值精确一致（g=2 d=8=346860150644700、
  g=3 d=8=800457740515775）。
- **毫秒级负载被固定开销主导**：Rust 冷启动多出 ~6 ms（进程 + SQLite 建库/WAL），这是持久化设计的固定成本，
  随负载增大被摊薄；原版纯内存进程无此成本。
- **计算主导区间 Rust 反超**（g=2 1.5×、Getzler 1.3×），且算得动 oracle 算不动的区间（g=2 d≥7、g=3）。
- **暖缓存近乎免费**（4–5 ms 与负载无关）：SQLite 二级缓存兑现了"跨查询/跨会话复用中间数据"的核心设计目标。
- oracle 在 g=2 d≥7、g=3 崩溃（栈溢出）是**预言机桩自身的残余缺陷**（已排除 i128 溢出——计数器为 0；
  已修复有理数键文本编码歧义与 GWdbm/SYMdbm 共库两处硬化；成因未完全定位，留待后续），与算法无关。



## GROWI 全量验收（2026-10-05，与 sc-faber 一档完成联动）

**阶段 1：gwcheck.cm 全部 41 条 parity**（`tools/gwcheck_all.py`）——
41/41 精确通过。三轴：
- 正确性：Kontsevich P² d=15（39 位数）、Getzler P²/P³、EHX g=2/g=3 P²P³
  d≤8、五次三维流形 d≤6 全部与 gwcheck.cm 断言精确一致；
- 耗时：**冷启动（全新 SQLite 缓存）总 4.3s**，最慢单条 EHX g=3 d=8 = 0.85s
  （原版 2003 年硬件约 3 分钟 → ~40×）；
- 内存：重载探针（g=3 d=8 + d=15 + quintic d=6）峰值 python 进程 RSS ~53 MB，
  全量 41 条缓存持久足迹 10.6 MB。
quintic d=5：输出 …648 = gwcheck.cm 断言值 = 镜对称 GW 正确值
（GV 多覆盖复算 + 三文献佐证，见"d=5 定案"段）——原版无缺陷，
疑案撤销。

**阶段 2：faber ↔ growi 互为 oracle 交叉**（`tools/cross_growi_faber.py`）——
两引擎的共享数学面 = M̄_{0,n} ψ 积分（GROWI 度 0 点公理的内部地基，
faber klp 的 g=0 socle）。映射（实验确立）：
  "G=0,D=0 in P^{d} thru H:{d}, psi^{a1},...,psi^{am}" = ∫_{M̄_{0,d+m}}Πψ_i^{a_i}
合法组合（Σa = d+m−3）全量 34 例 **34/34 精确一致、零跳过**（d≤6、m≤3）。
g=1 侧的共享常数为 GROWI 公理硬编码的 1/24 = ∫_{M̄_{1,1}}ψ（faber klp 同值，
常数级交叉）。由此两机器在最深的公共地基上互相背书：GROWI 的全部 g=0
递归根落到该面，faber 的 g=0 socle 被 GROWI 十余年实战验证过的公理链独立复算。

**回归**：验收全程零代码改动；当时 Rust workspace 117 测试 + pytest 13 全绿
（终态计数见 README——后续引擎扩展有增加，此处为验收时点快照）。

## 发布前修订（2026-10-06）

发布前审查发现并修复的问题（均配回归测试）：

- **零类语义**：插入类恒为零（如 `H-H`）时原版返回 0（gwinv.cc 逐点先查零类），
  移植版曾把该点删掉，变成另一个不变量；现与原版一致。
- **数值截断**：解析器把 N/E/G/D、重数、重复次数 `as i16` 截断（`D=65537` 被当作
  D=1），整数累加可回绕；现全部做范围检查，超限报错。单项式指数以 u8 存储，
  类的总次数 > 255 时原先静默丢掉该变量，现报错。
- **关键字 `3fold`…`9fold`**：词法器先把数字吃掉，这些原版关键字从未生效；现数字
  开头的词整体作为词。
- **CLI**：JSON 应答改由 serde 生成（原手工拼接，错误信息含引号时输出非法 JSON）；
  `--cache` 缺值、未知选项报错；新增 `--no-cache`。Python `GROWI(cache="")` 只用内存。
- 缓存算法版本号更新为 `sc-growi 0.1.0`（旧库自动作废）。
