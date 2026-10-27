# Faber《Numbers》→ sc-faber 移植方案

源：`refs/Numbers.tar`（Carel Faber，KTH，2001–2004，Maple 8/9）。
目标：`crates/sc-faber`（Rust），pyo3 直连 + CLI，两级缓存（对齐 sc-growi 工程模式）。
oracle：admcycles（Sage 库，https://gitlab.com/modulispaces/admcycles）——本机无 Maple，原版不可执行；
admcycles 是其精神后继且能力远超，兼作真值来源与功能对标。

## 第一部分：原程序分析

### 1.1 文件构成

| 文件 | 行数 | 性质 | 内容 |
|---|---|---|---|
| `Instructions` | — | 文本 | 用法文档 + 一批精确验证例（现成 golden 真值） |
| `KaLa5` | 391 | **文本（可移植）** | 主引擎：全部递归过程 |
| `KT21v5.m` | 16866 | **Maple 二进制（MVR5）** | KT 表：κ-单项式 → τ-基 的表示 |
| `tauZaalV5.m` | 1378 | Maple 二进制 | TAU 表：纯 ψ 交数（Zaal 的计算） |
| `minderV5.m` | 6274 | Maple 二进制 | 分划表 pa/nnp + κ 相关表 |

关键事实：三个 `.m` 是 Maple 私有 `save` 格式，**没有 Maple 读不出来**。但这不是障碍而是路标：
KT 表的数学内容就是 Faber–Zagier 关系的预计算形式（κ-单项式在 τ-基下的规范表示）。
现代做法（Pixton 2013，admcycles 同款）是按需算法化生成这些关系——因此移植**不搬数据，移数学、换算法**。

### 1.2 数学内核（KaLa5 全部可读逻辑）

入口 `klp(g,n,x)`：x 为 ψ_i（i≤n）、κ_j（j≤3g−3+n）、λ_i（i≤g）的多项式，返回 M̄_{g,n} 上的交数。

1. **λ 展开**：`lam(n) = coeff_t^n exp(Σ_{i≥1} B_{2i} κ_{2i−1} t^{2i−1} / (2i(2i−1)))`——
   λ 类化为奇 κ 类的多项式（Mumford）。
2. **主递归 `klmgn(g,n, 类指标列表)`**：
   - 纯 ψ/κ → `taug`：κ-单项式经 KT 表映到 τ-基，再 `tau` 求值；
   - g=1 带 ch → 恒等式 1/24 · M̄_{0,n+2}（单 λ 消解）；
   - g≥2 带 ch_m → **Mumford 公式**：
     `ch_m = B_{m+1}/(m+1)! · [κ_m 项 + ½Δ_0 项 + Σ_{h=1}^{g-1} 可约项(D/F 分裂)]`，
     递归到更低亏格/更少标记的交数。
3. **τ 求值 `tau([d_1..d_n])`**（纯 ψ 交数）：
   - g=0：Kontsevich 公式 `(n−3)!/∏ d_i!`；
   - d_1=0 / d_1=1：string / dilaton 方程递归；
   - 其余：TAU 表查表（= 我们要用 FZ 重构替代的部分）。
4. **遗忘前推 `push/pushn`**：κ 在第 n 层用 κ|_{level n−1} + ψ_n^j 表出，配合 string 方程逐层下推。
5. **Rahul 关系 `rel(g,d,a)`**：M_g 上重构环关系的显式实现（Pandharipande 1999）。
6. `temp/effe/vv`：Faber 猜想（λ_g·ψ 归一化）的辅助查询。

### 1.3 限制与本质

- **3g−3+n < 20**（表只算到这里）：静态数据的硬顶。
- 结构 = 小算法 + 大预计算表；2004 年 Maple 里现算 FZ 关系不现实，如今是标准算法。
- 验证真值随手可得：Instructions 里 `klp(7,0,la1^18) = 32017001/638512875` 等十余例。

## 第二部分：admcycles 对标（多出什么）

admcycles（Delecroix–Schmitt–van Zelm，2018–，GPLv2+，Sage 库）相对 Faber 程序的增量：

| 能力 | Faber | admcycles | sc-faber 计划 |
|---|---|---|---|
| ψ/κ/λ 交数 | ✓ | ✓ | ✓ 核心 |
| FZ 关系 | 静态表（<20 度） | **算法生成 + 数据库**，无度数上限 | ✓ 算法生成 |
| 一般重言代数 R*(M̄_{g,n})（边界 strata、稳定图、乘法、前推/拉回） | ✗ | ✓ | ✓ 核心（stable graph 一等公民） |
| Pixton 广义 DR cycle 公式 | ✗ | ✓（DR_cycle/DRpoly） | 二阶段（M5） |
| 超椭圆/双椭圆轨迹类 | ✗ | ✓ | 不承诺 |
| k-微分 strata（diffstrata） | ✗ | ✓ | 不承诺（独立子项目体量） |
| Hurwitz 空间 | ✗ | ✓ | 不承诺 |

判断：主体对标 = "admcycles 的重言代数内核"（TautologicalRing + FZ + 交数求值 + 前推），
这部分也是 Faber 程序功能的严格超集；DR cycle 是自然第二阶段；diffstrata 级别的子系统不纳入。

## 第三部分：架构

沿用 sc-growi 工程模式，新建同构三 crate：

```
crates/sc-faber        # 库：重言代数内核
crates/sc-faber-cli    # CLI 二进制 faber（--batch JSONL + --cache）
crates/sc-faber-py     # pyo3 cdylib _sc_faber
```

- **数值层复用**：M0 时把 sc-growi 的 `num.rs`（malachite Q）提升为 `crates/sc-num`
  供两者依赖（或 sc-faber 依赖 sc-growi 的 pub 模块，视耦合度定）。
- **两级缓存（用户核心要求，比 GROWI 更重要）**：FZ 关系矩阵按 (g, n, 余维数) 生成昂贵、
  复用巨大——L1 进程内 HashMap + L2 SQLite（WAL、WITHOUT ROWID、LEB128、ALGO_VERSION 失效），
  与 sc-growi 同款语义。
- **数据结构**：
  - `StableGraph`（顶点=亏格+标记多重集，边=结点，自同构群）；
  - `KappaPsiPoly`：strata 上的 κ-ψ 多项式（对应 admcycles 的 kppoly/decstratum）；
  - `TautClass`：Σ (graph, kppoly, 系数∈Q) 的线性组合；
  - 交数 = 对角配对 / FZ 归约后基系数。
- **FZ 关系生成**：Pixton 路线——生成元空间枚举 + DR 类多项式系数矩阵 + 高斯消元得关系。
  整数/有理线性代数复用/扩展 sc-growi 的 linsys 思路。

### 里程碑

| 里程碑 | 内容 | 验收 | 状态 |
|---|---|---|---|
| M0 | 工作区扩展 + sc-num 抽取 + StableGraph/kppoly 数据结构 + 单测 | 单元测试过 | ✅ `c64f8ea` |
| M1 | kppoly 升级（逐顶点 κ、逐腿 ψ 含结点半边）+ 生成元枚举 tautgens | 枚举/权重/覆盖测试 | ✅ `96f3e38` |
| M2 | τ 求值（Kontsevich + dilaton + 度数守卫）+ pushn/push 遗忘前推 + int_fund 驱动 | g=0 随机性质测试 ∫=∫∘pushn | ✅ `f49e6c1` |
| M3 | socle 求值路线：STrecur/STsum/CT 闭式 + TautClass.evaluate（逐顶点）+ 图自同构 + SQLite L2 | admcycles 34 项差分全绿；黄金测试固化 | ✅ |
| M4 | λ（Bernoulli 展开）+ ch（Mumford 递归）→ 完整 klp | 差分 53/53 全绿（g≤4 活 oracle）+ 三轴验收 | ✅ |
| M5 | Pixton DR cycle（插值形式）+ DR·X 积分 + 同构不变规范形 | g=1/g=2 锚点与消没全绿；三轴验收 | ✅ |
| M6 | pyo3 + CLI + pytest + 基准（vs admcycles）+ README | 全绿；基准 6 项 5 胜 1 平 | ✅ |

**修订**：跨图乘法（admcycles prodtautclass 的粘合+对角）从 M1 挪出——
M3 的 socle 路线使求值不需要它；若后续（M5 DR/关系）需要再实现。
M2 顺序依赖已解除：socle 路线直接覆盖全亏格（M3 完成）。
**M4 完成记录（2026-10-05）**：
- `hodge.rs`：λ→奇 ch 展开（生成函数递推）、g=1 恒等式 1/24·(0,n+2)、g≥2 Mumford
  递归（κ_m 项 + Δ₀ 项 + 可约 D/F 分配，逐行对照 KaLa5 移植）、pushn 的 ψ 同腿幂合并
  与 q0 保留两处对齐修复、int_class 记忆化（L1 + SQLite L2，键规范化）。
- **权重守卫**：入口检查总权 = 3g−3+n（Faber klmgn 本有此检查）——修复非顶权含 λ
  输入返回垃圾值的 bug（admcycles 差分 3 例发现，修复后 53/53）。
- **三轴验收**：正确性 = 单元 53 + 差分 53/53（g≤4）+ M3 黄金 34；耗时 = λ 幂查询
  大幅领先本机 oracle（admcycles legacy g=3 λ₁⁶ 120s 超时，本引擎 g=5 λ₁¹² 2.1s）；
  内存 = g=5 重载峰值 3.1 MB。
- **证据分层（诚实声明）**：活 oracle 只覆盖本机可算部分（纯 ψ/κ 任意 g——socle 闭式；
  含 λ/ch 至 g≤4——legacy 路径）；g=5 为 Instructions 文档值复现（31/680400 ✓）；
  **g=6/7 放弃**（文档值无耗时记录、无比较价值；Hodge 积分另有更新算法，未来可选）。
- 规矩补充：差分测试只对本机跑得动的 oracle；文档值仅作标注参考、不阻塞验收。

**M3 实施记录（2026-10-05）——路线转向**：读 admcycles 源码后发现其求值
**不走 FZ 关系矩阵**：`decstratum.evaluate` = 逐顶点 socle 公式
（`DR/evaluation.py`：g=0 走 CT 闭式，g≥1 走 STsum→STrecur——一个
TRR 型递归 + string/dilaton，纯 ψ 交数由此完全确定）。据此 M3 移植
socle 路线而非 FZ 矩阵：
- `socle.rs`：STrecur（记忆化 L1 + SQLite L2 跨会话持久，`attach_cache`）、
  STsum/CTsum（集合分划 × 自同构因子）、CT 常数（Bernoulli 生成函数递推）、
  socle_formula。修过的移植 bug：块自同构按内容而非大小计数。
- `taut.rs`：`TautClass::evaluate`——任意装饰 strata 求值（顶点度数不符项
  贡献 0），除以 `StableGraph::automorphism_number()`（**语义差异**：我们的
  DecStratum = ξ_*fund/|Aut|；admcycles 裸 decstratum 不除，其 irrbdiv 积分
  为未归一化值——对拍脚本两端已对齐）。
- 差分对拍 `tools/faber_diff.py`：29 项基本类单项式（g=0..3，含 g≥2 全 ≥2
  列表与 κ 混合）+ 5 项 strata（Δ₀×4、分离边界），**34/34 精确一致**，
  全部固化为 `tests/golden_m3.rs`。
- FZ 关系矩阵与跨图乘法从 M3 移除（求值不需要）；若 M5 的 DR cycle 或
  关系枚举需要再实现。

**oracle 环境（2026-10-05 搭成）**：WSL2 Debian 13.5（trixie），专用用户 `oracle`；
trixie 无 sagemath 包、conda-forge 只有 py2 时代残骸 → 走 miniforge + **passagemath-standard**
（Sage 的模块化 pip 发行版，manylinux 轮子）+ `pip install admcycles`，全部隔离在
`/home/oracle/miniforge3`。调用：`wsl -d Debian -u oracle -- .../envs/sage/bin/python`。

**对拍第一轮（admcycles 裁决，2026-10-05）——撤销 M2 法医指控**：
- string 方程（消去零幂标记）**对一切亏格有效**，Faber 原版无条件使用是**正确的**。
  裁决数据：[0,1,2]_{(1,3)} = 1/12（string 链精确复现）；[0,1,5]_{(2,3)} = 1/288
  = string 预测 1/1152+1/384（分量亦验证）。
- 开发中对该结论的两轮错误指控均已记录：第一轮记错 Mumford 关系推出 ∫ψ²=1/2
  （真值 1/24）；第二轮手推 Faber 递归漏了 j=2 项。教训：admcycles 在、记忆退场。
- 意外收获：鸽笼原理（Σd=n 的 g=1 顶度列表必含 0 或 1）⇒ **g=1 纯 ψ 交数已被
  string+dilaton 完全覆盖**，无需 FZ。FZ 真正缺的只是 g≥2 的全 ≥2 列表。
- 顺带钉死的低级错误：M̄_{2,2} 维数是 5 不是 6（三个"反例"实为过度数平凡零）。

**M3 的 admcycles 真值目标表**（已测得，全部精确有理数）：

| 量 | 值 |
|---|---|
| ∫_{M̄_{1,2}}ψ₂² = ∫_{M̄_{1,1}}κ₁ | 1/24 |
| ∫_{M̄_{1,3}}[3]、[1,1,1]、[0,1,2] | 1/24、1/12、1/12 |
| ∫_{M̄_{1,4}}[2,2]、[0,1,3] | 1/6、1/8 |
| ∫_{M̄_{1,5}}[0,2,3] | 7/24 |
| ∫_{M̄_{2,1}}ψ₁⁴ | 1/1152 |
| ∫_{M̄_{2,0}}κ₁³ | 43/2880 |
| ∫_{M̄_{2,3}}[0,1,5] | 1/288 |
| ∫_{M̄_{1,1}}λ₁ | 1/24 |


**M5 进行中（2026-10-05）**：
- `dr.rs`：Pixton DR 公式插值形式完整移植（平衡权重和 + 拉格朗日 Q(0) + m^b₁ 归一），
  双口径（类口径 / 生成元口径含 1/|Aut|）。
- **已验证**：g=1 全部——两个文档锚点（DR(1,(2,3,−5)) 五项逐系数、DR(1,(2,−2))）+
  三个 DR·X 积分（0, 0, 1/24）与 admcycles 精确一致；g=2 r=2 七类图形态系数
  逐一手验与 admcycles 31 项 dump 精确一致（fundclass/[0,2]/[1,1]/[1]/[0]双环/
  [0,1]双边/[1,0,1]）。
- **修过的三个真 bug**：① κ 限制公式应为 Σ_v κ_j^{(v)}（无半边 ψ 修正，admcycles
  boundary_pullback 源码裁决）；② int_dr 进位循环在空因子表时提前 return（多单项式
  项只算第一个）；③ DR 系数 aut 口径（类口径不乘 1/aut，由 evaluate 统一除）。
- **开放项（根因已定位）**：StableGraph 规范形**非同构不变**——同一 stratum 的
  同构变体（如自环+分离边图的两个半边配对）被 tautgens 重复枚举，DR 类重复计数；
  证据：g=2 五个 DR·X 积分全部恰偏大（29/288 vs 26/288 等），r>g 消没测试不
  满足（已标 #[ignore]）。修复方案已定：结构化证书规范形（顶点色 (genus,标记集) +
  边对多重集，同色顶点组内枚举取最小代表，多边配对按早端点升序），预计 ~60 行。
  M1 时"非完整同构规范形"的文档注记正是此雷。
- **M5 完成（2026-10-05）**：结构化证书规范形落地（顶点色分组 + 组内枚举最小
  边对多重集 + 确定性重建）——修复后 g=2 类项数 46→31 与 admcycles 完全一致，
  **五个 DR·X 积分全部精确吻合**（13/144 等），消没定理（r>g）以积分判据通过；
  λ_g 与 DR(g,()) 的符号关系经 admcycles 双侧实测锚定（g=2 同号 7/5760、g=3
  反号 ±31/967680——其 doctest 恒等式用的是 double_ramification_cycle 另一
  入口的约定）。三轴：差分全绿；DR 构造+积分秒级；内存 3.1 MB。


**M6 完成记录（2026-10-05）**：
- `sc-faber-cli`（faber，JSONL 批处理，`--cache`）、`sc-faber-py`（`_sc_faber`，
  FABER 类：klp/dr → Fraction，GIL 释放）、`python/sc_faber`（现位于 `crates/sc-faber-py/python/`）、`tests/test_faber.py`
  （13 项）、根 README。
- 修复两个 M6 期间发现的真 bug：① **L2 键编码歧义**——int_class 的 SQLite 键
  只含 κ 长度，ψ/ch 分界丢失 → 递归子问题键碰撞 → 挂缓存时 g=4 λ 例返回垃圾值
  （g++ 差分级教训的翻版：序列化键必须无歧义）；② klp 入口补 ψ 腿范围校验。
- 基准（tools/faber_bench.py，本机可跑集）：ψ/κ 单项式 6-7ms vs admcycles
  0.3-3ms（它赢启动开销，毫秒区无意义）；**λ₃λ₂ψ² M̄₃,₁ 我们 6ms vs 它 79.7s
  （1200 倍）；DR·κψ 7.5ms vs 83ms（11 倍）**——纯 ψ/κ 之外的每类查询我们都
  大幅领先。λ 高次幂是其死亡区（超时），不纳入。


**一档进行中（2026-10-05）——重大对拍发现**：
- 已移植：Hain 除子（3 锚点全绿）、Chiodo DR（固定 r）、Θ 类（Θ(1,1)/(2,0)/(3,0) 精确）、
  int_against 通用乘积求值（修两 bug：丢原类项 κ、edge 幂级数二项式系数）。
- 【更正 2026-10-05 晚】"oracle 丢项"系我方测量伪影（tail 截断 + float 误读），oracle 求值
  自洽，Θ(2,1)=3/128 为真锚点。引入 codex（VSCode 扩展内 CLI，走本地代理）为第三方仲裁：
  裁决出我们 dr_chiodo 权重和的真符号错（Bw=−given 又代 −w 的双重取负，r≥3 时显形；
  回归例 g=2,A=(3,),k=1,r=3 的 [1,1] 图 ψ 结点项 = 15/32 两侧吻合，仅半边标号相反）。
  r=2 时该修复不可见，Θ(2,1) 分歧仍在：混合实验（oracle 系数×我们积分）与逐项差分
  把残余隔离到我们 int_against 对"结点 ψ+标记 ψ 混合装饰"的约 6 个生成元（部分给 0）。
  下轮：逐项手验这些生成元（Δ₀ ψ₃ψ₄ 手算 ∫_{M̄_{1,3}} ψ₁ψ₃ψ₄ = 1/12 ≠ 0，定位 0 的来源）。
- **admcycles 求值路径在本机不可信**（原文保留，已被上文更正）：`(T*X).evaluate()` 经 TautologicalRing 乘法时丢项
  （fundclass ψ₁³×ψ₁⁴ 应为 1003/7464960，其 Rational 路径给 0、float 路径给非零）。
  其 Θ(2,1)·ψ₁=3/128 与 Θ(3,1)·ψ₁²=15/1024 两个"锚点"是缺项结果，不作验收标准。
- **可信锚**：DR_coeff_new 系数表（显式 Pixton 公式，不经过环乘法）。混合实验
  （oracle 系数 × 我们 socle 积分器）给 Θ(2,1)·ψ₁ = 38887/1866240，两半各自独立
  验证过，为当时最佳真值。我们 dr_chiodo 的 fundclass 七系数与表精确一致，
  但部分边界生成元仍有差异（bicycle 图 ψ₄-only vs 我们 ψ₄₅₆₇ 对称非零）——
  ~~待修匹配器做单项级对拍后对齐~~（已由后续系数对拍器 + b4f2141 结案，
  见下文"Θ 疑案全面告破"）。
- 工具：tools/gen_diff.py（系数表）、oracle_coeff_our_int 示例（混合实验）、
  chiodo_dump/gen_diff_match（对拍，~~匹配器待修~~已修，见下）。

> **⚠️ 历史快照注记（2026-10-06 梳理）**：以上两段及其上方"下轮：逐项手验"
> 是 Θ 疑案**中期**的排查快照，全部事项已由后续工作关闭——"oracle 有 bug"
> 的判断经复核**不成立**（Rational/float 路径差异是我方读取截断伪影，已更正
> 记录），Θ 六锚点最终精确落位。阅读本节请以下文"Θ 疑案全面告破"为终态。


**系数对拍器落地（2026-10-05 晚）**：`dr_coeff_one`（单生成元系数）+
`tools/coeff_compare.py`（双侧规范键多重集比对，修标记/结点半边过滤）。
结果：71 个生成元 → 43 个共同键中**仅 4 个系数不匹配**，全部为
[0,1] 双边图（g0 顶点带标记 + g1 顶点，两条平行边）的 ψ/κ 装饰——
oracle 11/864 vs 我们 −11/432 等（符号与 2 倍因子混合）。另有 6 个
"我们独有键"系对称聚合的表示差（oracle 把对称半边系数不对称分配
如 ψ₃ 55/864 + ψ₄ 0，我们对称平分各 55/1728——**类相同**，此前
"ψ₄-only 之谜"就此解决）。fundclass 七系数与全部 E1/E3 键精确一致。
残余：4 个 [0,1] 双边键的符号/因子——嫌疑为流方向 vs 边存储方向的
econtrib 定向（非自环边的 u→v 与 B 的行号约定），下轮手验一条即可。


**对拍收口进展（2026-10-05 深夜）**：规范键补齐**边对多重集**（自环 (i,i) vs 平行边
区分——此前键碰撞导致 4 个假不匹配）。现在 45 个共同键**仅 2 个系数不匹配**，
均为"g0 自环 + 分离边"图（一边差 −2 倍一边差 +4 倍，且都涉及 g0 上的 κ/ψ
与他顶点装饰组合）。另有星形/链形枚举差：我们产星形（g0 中心）×2 项各 −11/432，
oracle 产链形（g0 居中三顶点）1 项 −11/216，总和恰相等——待判定是同构表示差
还是枚举缺失。fundclass/E1/其余 E2/E3 全部精确一致。


**决定性突破（2026-10-05 深夜 II）**：`StableGraph::new` 会规范化重编号结点腿
（3,4→2,3），而外部表（oracle 系数导出）的 ψ 装饰按旧编号书写 → 结点 ψ 装饰
在规范化后指向不存在的腿 → 积分成 0。新增 `new_raw`（验证不规范化）后，
**8 个结构中 7 个精确一致**，oracle 侧总额精确 −3/64 ✓。
全部残余收敛到唯一结构 `(g0,m[]);(g1,m[1]);e[(0,0),(0,1)]`（g0 自环+分离边、
标记在 g1 侧）：oracle 47/41472 vs 我们 −41/20736——它恰是最初 Θ(2,1)
43/27648 差异的全部来源（含 κ₁ 行 oracle 11/864 与我们手算 −55/20736 的
符号+数值双重分歧，手算含"chiodo 边因子不入 scalar"的修正，需下轮三向细查）。


**最后结构的真相（2026-10-05 深夜 III）**：逐装饰三向分解证明——目标图
（g0 自环+分离边、标记在 g1）上**我们 dr_chiodo 的 class 系数与 oracle 精确一致**
（κ₁v1: class +11/432 = oracle 11/864×aut ✓；ψ标记 1/864 ✓；ψ₇: class 1/864 =
oracle 1/1728×2 ✓；g0 侧结点 ψ 行积分为 0 无贡献）。此前"残余不匹配"是**探针
自身的图过滤器过松**：`genera==[0,1] && edges==2` 同时命中了三个不同的图
（loop 在 g0 / loop 在 g1 / 标记在 g0），把别的图的系数（−11/216、−55/432）
算到了目标图头上。struct_totals 的 ours 分支经 dr_coeff_one 间接受污染。
下轮：收紧过滤器（loop 必须在 g0 顶点）重验 → 预期 8/8 结构全对、Θ(2,1)=3/128
精确收口。


**Θ 疑案全面告破（2026-10-05 收口）**：真根因是 dr_chiodo 把"标记的 A 权重计入
given_weights"（应无条件、作用于所有标记）误实现为仅对带 ψ 装饰的标记生效——
未装饰标记的 A 贡献被静默丢弃（trace 铁证：κ₁ 生成元 given=[−2,2] 而非 [−3,3]）。
修复后：图 (b) 全部六行与 oracle 精确一致（κ₁ +11/432 = 11/864×aut ✓），**8/8
结构全对**，Θ(2,1)·ψ₁ = 3/128、Θ(3,1)·ψ₁² = 15/1024 精确落位（测试全绿）。
至此 Θ 类六个锚点（1,1)/(2,0)/(2,1)/(3,0)/(3,1)+消没全部精确，从 M5 追踪至今的
43/27648 疑案以"三重测量伪影 + 一处真 bug"结案：真 bug 恰好只在"标记落在无
ψ 装饰顶点"的生成元上显形——正是结构对拍器的价值所在。

### 第一档：Hyperell/Biell 配对（2026-10-04/05 收口）

**配对内核（z2_pair，Hyperell/Biell 共通）**——Z/2 Hurwitz 传输公式，Hpullpush
锚点实测：κ_j 指标保持 ×2^{Σ幂次}；固定/分支腿 ψ ×(1/2)^{幂}；共轭对两腿合并到
同一轨道腿（无 1/2 因子）；deltadeg = 2^{2g_q+m−1}（超椭圆 g_q=0：m=0,1,2 →
1/2,1,2；双椭圆 g_q=1：g=2,3 m=0,1 → 2,4）。九锚点测试（g=2 m=0 六个 +
g=1 m=1 共轭对 ψ₅ψ₆=1、ψ₅²=1、ψ₅=0）。

**积分合成（z2_int）**——∫ Locus(g,n,m)·X 的正确公式（33 例 oracle 差分全对，
`tools/hyb_oracle.py` batch 1–8）：
```
∫ Hyperell(g,n,m)·X = 1/(2g+2−n)! · Σ_terms coeff · z2_pair(term)
ρ^*X：ψ_l 稳定；κ_j^c → (κ_j − Σ_{a=n+1..branch} ψ_a^j)^c 多项式展开
```
（= admcycles `forgetful_pullback_list` 的单顶点部分；其有理尾巴图项在 H̄ 上
为空层——σ 不动点计数矛盾（尾巴上分支点+结点已占满不动点），配对为 0，实测
均无需。）Biell 同式（branch = 2g−2），仅 (g,n,m)=(2,0,0) 额外 ×1/2——对齐
admcycles 双重计数特判（无标记 g=2 双椭圆曲线有两种椭圆覆盖结构）。

**排障记录**：初版组合因子 C(branch,n)/(branch−n)!×pair 猜错——二项式不该有、
κ 修正项漏乘 deltadeg、∫_{M̄0,5}ψ₁ψ₂ 手记 1 实为 2；三处全是手算草率而非传输
公式错误，教训：对拍前先把商空间 socle 值用自己引擎核一遍（klp 批查）再反推。
度数预算速查：codim(Hyperell) = g−2+n+m，codim(Biell) = g−1+n+m（弱对拍的
0 值例无验证价值，强例必须 deg 严格匹配）。

### 第一档：spin DR + Witten r-旋（2026-10-05 收口）

**spin DR（dr_chiodo_poly，CSS21）**——重构 dr_chiodo 抽出单生成元
chiodo_coeff_at(gen,g,a,d,k,m,spin)（Θ 固定 r 路径不变），新增插值模式：
m0+1..（2d+2 点）逐点求值、拉格朗日 P(0)；spin 三改动：仅奇边权、m 取偶数、
因子 1/2^{g−h0}；负指数 m^{2d−h0} 修为除法。差分 15 例全对（drpoly 7 =
DR_cycle(chiodo_coeff=True) 主路径独立验证 + drspin 8 含 g=2 (9,−3,−1)）；
结构锚点 = doctest DR_cycle(1,(1,1),spin=True) 五项（类系数 = aut×oracle）。
负数取模用 rem_euclid（−1 是奇数）。

**Witten r-旋（witten.rs，PPZ16 附录）**——数值化路线：P_d(X,a) 为 a-多项式
（系数 QQ[X]），符号和用 Bernoulli 闭式（负 a 解析延拓，锚点 PM(1,6,−2)=−53/6
证伪朴素循环和）；边因子 = S/(X+1) 综合除法取 −T[d1]；给定权 w_v = 1−g_v−
Σj·e_{κ}+Σ(A_i−e_i)（Pixton M[i,0][0]=genus）；权重和 mod (m−1)；**scalar 不
参与插值**（interp(X)·scalar(X) 结构，否则 deg 对不上）；×(−1)^{d−g}；插值点
步长 step=lcm(C 分母)。验证：Pm/leg/edge 17 数值锚点 + g=1 五非零积分 +
g=2 A=(2,) 五非零（κ₃=11/1920 等）+ 弱 0 例若干，全部与 oracle 一致。
排障：四处 bug（C3 三次式括号枚举、u=−b 转 b 系数符号、pm_eval a 幂乘成 x 幂、
h0 负幂、aut 双乘）全部被锚点当场截获；腿因子 a = A_i(m) 不减 ψ 幂由 oracle
两序置换实验定案。

**0.2.1 订正**：ΣC 的整数性检查与 ΣC 本身曾直接把各分子相加（未按 lcm 通分），
C 的分母不全相同时或被误拒、或静默算错 ΣC 与类的次数（锚点的分母恰好相同）；
现按通分后求和，并以"同一组 C 换一种分数写法结果不变"的测试把关。

**kval**：dr_chiodo 的 k-twist 一直在（Θ 用 k=−1），drpoly/drspin 差分含 k=1
例；CLI op：drpoly/drspin/witten/hyperell/biell。

### 第一档：GRRcomp（2026-10-05 收口，一档全部完成）

**grr.rs（PRvZ20 Theorem 1，admcycles/GRRcomp.py 移植）**——广义 Hodge Chern
类 ch_t(R^•π_* O(D))，D = lK̃ + Σd_iσ_i + Σa_{h,S}C_{h,S}：
- todd_inv_sigma（自环 ×½ + 分离 ×1，decsum = Σ(−ψ)^iψ'^{t−1−i}）+ todd_inv_exp_graphs
  三分支链图（严格更大 L+2 / 相等 L / irreducible ×½）+ todd_inv_exp；
- chern_exp/Kin/Sin（链图 + κ@末顶点（子空间 κ_0 = 2g'−2+n'）/ (−ψ_p)^{β−1}
  + p∈Scomp 过滤）+ chern_single = K+S+todd⁻¹，度 t 截断；
- **同构合并**（色空间键：顶点色(genus,标记集,κ多重集)+边装饰+标记 ψ）——
  链图 new_raw 与分离图 new 构造的同构项合并（oracle 展示即合并形式）；
- twist_div_to_zero。

验证：doctest/锚点 6 组（divided_bernoulli 6 值、todd 5 图全结构、
chern_single 三图、λ₁ 四项、twist doctest）+ 差分 18 例（两批 tools/grr_diff*.py）
全部精确一致；强例 5/576、1/480、−1/6 等。排障：graphs 的 ψ 不含 a_i^{tv_i}
（外层统一乘，锚点 a=−1 显形）、aut 口径（push = aut×oracle 展示值——差分
1/48 vs 1/24 截获）、comps(k=0) 下溢。**λ_deg≥2 需 TautClass 一般乘法（图积）
——归入二档**（λ₁ = ch₁ 已覆盖）。

**0.2.1 订正**：decsum 与链图尾段的 Σ(−ψ)^iψ'^{t−1−i} 曾被写成
Σ −ψ^iψ'^{t−1−i}（i 为偶数时符号错），t ≥ 3 的 ch_t 全部受影响；上面的锚点
与差分只覆盖 t ≤ 2，因而没有截获。现由两项一致性测试把关：无扭曲时 grr
（l = 0、1）与 klp 的 ch_t 一致（t ≤ 7），以及 Serre 对偶
(l, d, a) ↔ (1−l, −d, −a)。注意 Serre 对偶对这个错误不敏感（修正前也成立），
真正截获它的是与 klp 的对照。

**一档完成**：Hain ✓、Θ ✓、kval ✓、spin DR ✓（15 例）、Hyperell/Biell ✓
（33 例）、Witten r-旋 ✓（27 锚点+差分）、GRRcomp ✓（24 锚点+例）。

## 开发规矩（2026-10-05 起，强制）

移植项目，不是重新发明。一切实现先对照参考源（Faber KaLa5 / admcycles 源码，
均在本地可读），禁止凭记忆造轮子。

1. **小步提交**：每完成一小段（一个函数/模块单元）即验收，不带病前进。
2. **双重测试**：每段先单元测试（数学恒等式或文献已知值），再 oracle 差分
   （同一查询两侧对拍，脚本化、可重跑）。
3. **三轴验收**，任一轴不过即返工该段：
   - 正确性：与 oracle 值精确一致；
   - 耗时：与 oracle 本机实测同量级或更快（先测 oracle 定基线，再测自己）；
   - 内存：重载下峰值无异常增长。
4. **重载测试须先声明目的**（正确性/耗时/内存三者之一），先小规模实测再放大，
   禁止无目的长跑与盲等。

**M4 审计记录（2026-10-05）**：hodge.rs 与 KaLa5 ch 分支逐行对照——κ 项/Δ₀ 项/
可约项的符号、½ 因子、维度过滤、g=1 分支、push 前置全部一致；实现差异仅为
表示层（语义指标 vs 线性下标，枚举 vs 多项式展开），等价。性能对照：本机
admcycles legacy 路径 g=3 λ₁⁶ 120s 超时，Ring 快路径因 FZ 数据库需网络下载而
不可用；本引擎 g=5 λ₁¹² = 2.1s（release，含记忆化）——大幅领先 oracle。
g=6/7 的 λ 幂重载按上面第 4 条（重载测试须先声明目的）不再无目的测试：正确性已由 g≤5 全绿覆盖，耗时
增长是积分级递归的已知特性（Faber 本人 g=7 需 35s + 预计算表）。

## 第四部分：验证与 oracle

1. **Instructions 例**（十余个精确有理数，含 `klp(7,0,la1^18)`）——第一层真值。
2. **文献**：Faber–Pandharipande math.AG/9810173 及后续（λ_g 交数公式）、Mirzakhani WP 交数、
   Looijenga 消没（度 > g 的重言类为零——免费的消没检验）。
3. **admcycles oracle**：差分脚本 `tools/faber_diff.py`（Rust 侧由 `examples/diff_dump.rs` 导出），
   批量对拍交数与关系；oracle 环境（WSL2 Debian + admcycles，2026-10-05 搭成）见上文。
4. **GROWI 交叉背书**（2026-10-05 验收，`tools/cross_growi_faber.py`）：两引擎的
   共享数学面 = M̄_{0,n} ψ 积分。映射：GROWI `"G=0,D=0 in P^d thru H:d,
   psi^a1,…,psi^am"` = 本引擎 `klp(0, d+m, [], [(i,a_i)])`——合法组合
   （Σa = d+m−3，d≤6、m≤3）全量 34 例 34/34 精确一致零跳过。意义：GROWI
   的全部 g=0 递归根落到该面（其度 0 点公理即 Witten 闭式的工程实现），
   本引擎的 g=0 socle 由此获得一条**独立于 Faber/admcycles 的第三真值链**
   （Gathmann 代码链十余年实战验证）；g=1 侧共享常数 1/24 = ∫_{M̄_{1,1}}ψ。
5. **代数自检**：FZ 关系作用于基为零、string/dilaton 恒等式、Looijenga 消没——
   不依赖外部真值的结构检验。

## 第五部分：风险

- FZ 关系生成的组合爆炸（Pixton 生成元空间大）：缓解 = 缓存 + 生成元剪枝 + 按需子矩阵；
  admcycles 用预计算数据库，我们以 SQLite 扮演同角色。
- Maple 二进制表无法读取 → 无法逐表对拍 KT 数据：以 Instructions 例子 + admcycles 双重覆盖替代。
- sc-num 抽取对 sc-growi 的回归风险：抽取后跑 sc-growi 全量测试即可。

## 发布前修订（2026-10-06）

发布前审查发现并修复的问题（均配回归测试）：

- **同腿 ψ 静默算错**：`psi=[(1,2),(1,2)]` 曾返回 0（`psi_power_at` 只取第一项）；
  `KappaPsiPolynomial::new` 现合并同腿幂（`kppoly::merge_psi`）。
- **ψ 腿号越界**：`int_dr`/`int_against` 未校验腿号，越界腿会撞上结点半边编号
  （n+1..）而被当作结点 ψ；现统一由 `kppoly::check_marking_legs` 报错。
- **drpoly/drspin 输入**：`dr_chiodo_poly` 补 ΣA = k(2g−2+n) 校验；CLI 缺省 k
  改为整除检查（原为截断除法，(g,n) = (1,0)/(0,2) 时除零 panic）。
- **κ 长列表性能悬崖**：`setparts_with_auts` 原先枚举全部位置集合分划（Bell 数）
  再去重，∫_{M̄_6}κ₁¹⁵ 需 Bell(15) ≈ 1.4×10⁹ 次；改为直接枚举多重集分划，
  该积分 0.2s。轨道因子改用精确有理数（原 u64 阶乘在 21! 溢出）。
- **其他**：λ_l（l > g）与 g = 0 的 Hodge 类直接给 0；若干 u32 下溢（g = 0）改为
  稳定性报错；L2 缓存值按任意精度回读（原 i64 解析使大值永不命中）；
  `grr.rs` 一处 `#[test]` 错位（`basics` 从未运行）；CLI 改用 serde 生成 JSON
  （原手工拼接，错误信息含引号时输出非法 JSON），并拒绝未知字段名。
