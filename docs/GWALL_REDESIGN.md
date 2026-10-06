# GWall → sc-gwall 移植方案

Jim Bryan《GWall》（Maple，2002，refs/gwall/GWall.maple）的 Rust 移植。
数学：Pandharipande–Okounkov / Faber–Pandharipande 的 local curve GW 理论
（CY3 中曲线的 BPS 态，BP 论文 arXiv math/0009025，LaTeX 源码
refs/gwall/bp/2001-9.tex）。

## 架构

| 模块 | 内容 | Maple 对应 |
|---|---|---|
| `comb` | combinat 子集（permute/cartprod/partition）+ 置换奇偶 | with(combinat) |
| `qseries` | 截断有理级数 + Bernoulli + z/(2sinh(z/2)) 系数 | series/coeftayl |
| `gwpoint` | 纯点后代闭式（Faber 幂和 + 钩长维数） | GWpoint/power |
| `gwloop` | 环路后代双实现（外代数符号） | GWloop/fastGWloop |
| `gwl` | GW/GWhodge/GWomega 递归族 + gw 前端 | 同名过程 |
| `pf` | branch/Ztil/Z/连通化 L/BPS 反演/Y | branch*/Ztil/Z/Y + 论文公式 |

依赖注记：`with(group)` 在原版中零使用（纯残留）；combinat 只用三个
函数，自实现 <150 行。BPS 反演层（log 连通化 + GV 反演）**原版 Maple
无对应物**——原程序止于 Y（disconnected），论文表由 étale/GRR 公式
另算（见里程碑）。

## 关键移植决策

- **evils 魔数弃用**：原版符号公式依赖 Maple permute 枚举顺序且作者
  自注 m<9 有效；改用 parity(σ)（逆序数 = parity(disjcyc)）。
- **branch 公式化**：不抄 branch2..10 手写展开，运行时按 chernnum 分拆
  和生成（与手写版 16 项系数逐一对照）。
- **度 1 闭式特判**：Z₁ = (4sin²(t/2)/t²)^{g−1}（Pandharipande），绕过
  branch 上限（原版 t^10 / "No can do more than branch number 10"）。
- **BPS 管线**（M4 悬案破案后定型）：
  Z_d = exp(−ΣN₀t^{2j})·Ztil（去常映射）→
  L_d = [q^d]log(1+ΣZ_e q^e)（**度变量连通化**）→
  N^h_d(g) = [t^{2(h−(d−1)(g−1))}]L_d →
  n_d^h(g) = Σ_{h'≤h}Σ_{m|d} μ(m)m^{2(g+h')−3}·α_{g+h,g+h'}·N^{h'}_{d/m}(g)，
  α_{G,G'} = [r^{G−G'}](arcsin(√r/2)/(√r/2))^{2G'−2}（**注意 /2**）。

## 里程碑

- **M0**（cfef491）：组合件 + 截断级数。排障：z/(2sinh) 系数 = B_l(1/2)/l!
  （母函数取 x=½——初版 B_l/l! 被 z² 系数 −1/24 锚点截获）。
- **M1**（77b1cc9）：GWpoint/power/N₀。跨引擎交叉：N₀(g,2g−4)/(3g−3)
  = ±∫_{M̄_{g−1}}λ_{g−2}³（sc-faber klp 独立积分，h=2,3,4 三点咬合）。
- **M2**：GWloop 双实现互验（200+ 组合）+ fast 版跨组场景 5.7×。
- **M3**：GW 递归族 + 前端。验证：string 方程四例（左走 L₋₁ 递归、
  右侧直达闭式——两条独立代码路径同值）。
- **M4**：branch/Ztil/Z/Y/BPS。三个真 bug 被 BPS 表对拍截获：H₁ 首项
  符号、branch 阶乘重复乘、截断自适应公式错误（内存爆炸的根因——
  补上原版 branch≤10 守卫）。
- **M4 悬案破案**（codex 仲裁 + 本方数值验收，见下）：连通化与 α 归一化。
- **M5**：n_{3,4,5} 表 58/58、CLI（gwall）/pyo3（sc_gwall）绑定、本文档。

## M4/M5 悬案档案（codex 仲裁两案）

**工作流**：论文级数学问题交给 OpenAI Codex（读
refs/gwall/bp/2001-9.tex + 我们的代码）出方案，本方逐分量手核后实施
——比自行通读论文高效得多。

1. **BPS 管线两处隐性错误**（M4）：① Z 系数须经度变量 log 连通化
   （论文 N 定义在连通模空间，tex 345 行；GWall 的 Z 是 disconnected，
   Maple 无此步骤——原程序到 Y 为止）；② α 级数变量漏除 4 的幂
   （tex 737 行 arcsin(√r/2)/(...)）。此前 n₂ 表 h=0 行全对是两处误差
   相消/不可见的假象。附带教训：测试表曾误抄 étale 表行尾（论文两块
   表 local 比 étale 每行多一个非零格）。
2. **论文勘误四格**（M5）：n₃⁶(4)（论文表印 692352）违反其自身 étale
   界引理（tex 491-502：h ≤ (d_min−1)(g−1) 时 local = étale = 93913）；
   n₄ 三个行尾格源于 tex 622 行 GRR 展开 ε 指数印刷错误（正确
   μ(m)m^{2d_min(g−1)+1}，印成 m^{d_min(g−1)+2}；d=2 时 m=1 两指数相同
   ——n₂ 表全对的原因）。四格以我们的全量 GW 反演值为准（测试
   bps_paper_erratum）。

## 验证体系

1. **BP 论文 local BPS 表**（黄金数字）：n₂–n₅ 四张表全部 121 个数值格
   （n₂ 33、n₃ 42、n₄ 18、n₅ 28；含 4 格勘误；n₅(g=2) 含未知参数的格除外）
   全部精确一致——GWpoint 的 dim 假设、branch 层、递归族、连通化、反演链的
   全链背书。另有论文 §1 引述的 g = 0、g = 1 完整结果（n_d^h(0) = δ_{d1}δ_{h0}、
   n_d^h(1) = δ_{h0}）。
2. **度 1 全系锚**（Pandharipande 闭式）：g=2..4 全系数命中。
3. **跨引擎交叉**：N₀ ↔ sc-faber klp 的 λ³ 积分。
4. **内部互验**：GWloop 双实现、string/dilaton 方程、branch 与手写
   展开对照、度 1 闭式 ↔ branch 路径（ztil(·,1) 双路径一致）。

## 发布前修订（2026-10-06）

发布前审查发现并修复的问题（均配回归测试）：

- **g = 0**：α 中 A 的幂 2G′−2 与度 1 闭式的幂 g−1 在 g = 0 时为负，原用 u32
  计算——release 下回绕成约 4×10⁹ 次级数乘法（等于卡死）。现以级数求逆支持负幂
  （`QSeries::inverse`/`pow_signed`），g = 0 结果与论文（Faber–Pandharipande）
  一致：仅 n_1^0(0) = 1。`y_series` 在 g = 0 时为 Laurent 级数，明确报错。
- **度 1 闭式的截断长度**：g = 1 时 Z_1 曾返回长度 1 的级数，而 `QSeries::mul`
  按较短因子截断。g = 1 的连通级数在数学上恰为常数，故未产生错值（论文 g = 1
  结果全部吻合），但属潜在隐患，现保持完整长度。
- **环路顺序**：原先要求 a 环路全部在 b 环路之前（仅文档约定，交错输入静默得到
  错误符号）；现 gw 前端自动计入换序符号，任意顺序均可。
- **截断上限**：`bps` 超出 branch ≤ 10 时原先算完低阶项才报错（可达数分钟），
  现入口处立即报错。
- **测试**：`pf.rs` 一处 `#[test]` 错位（`arcsin_series` 从未运行，
  `bps_d345_table` 被注册两次）；补全论文表 n₂(g=7)、n₃(g=7)、n₅(g=5) 三行；
  `y_series` 改为以度 1 精确闭式 Y = (−y)^{g−1} 检验；dilaton 测试清理。
- **接口**：CLI 插入项支持数组与对象两种写法（原 `["ldeg",V]` 会把 V 解析到 k），
  `d` 改为必填，级数结果为 JSON 数组，批处理共享记忆化；Python 绑定接受元组或
  列表并在计算期间释放 GIL。
