# Hodge integral enumeration / Hodge 积分枚举

## English

This pure Rust example computes every integral of a monomial in
`lambda_1, ..., lambda_g` of weighted degree `3g-3` on the moduli space of
unmarked stable genus-g curves, for `2 <= g <= 20`. Results are exact rationals,
including zeros. No Python bindings or Sage installation are needed to compute.

Run from the repository root:

```sh
cargo run --release -p hodge-data
cargo run --release -p hodge-data -- --genus 6 --threads 8 --memory-gib 4
```

The default computes only g=2 through 4, with a 24-thread upper bound and a
2 GiB process-memory budget. Output defaults to `target/hodge-data/`, which is
ignored by Git. Use `--data-dir PATH` to choose another output and cache directory.

For the full range on a machine with sufficient memory:

```sh
cargo run --release -p hodge-data -- --min-genus 2 --max-genus 20 --threads 24 --memory-gib 100
```

The example calls the library's streaming batch API:

```rust
sc_faber::hodge::lambda_integrals(g, |indices, value| {
    // Write this monomial's exact value and update progress.
    Ok(())
})?;
```

The library owns partition enumeration, Mumford relations, closed formulas and
the GRR fallback. It reuses earlier values within each genus, retaining the
batch optimization rather than solving every monomial independently. The example
only handles arguments, cache ownership, output and progress. The thread setting
controls backend computation; monomials are visited sequentially. Values needed
by the batch remain in memory until the genus finishes.

The example's current scope stops at g=20; higher genera have not been validated.
The memory budget must fit the workload's remaining working set. SQLite
persistence and work-block spilling do not guarantee completion under
arbitrarily low budgets.

Each completed genus produces `hodge_data_gG.sage`. A descending tuple lists
lambda indices, so `(4, 3, 2)` represents `lambda_4 * lambda_3 * lambda_2`.
To read a result in Sage:

```python
hodge_number = globals().get("hodge_number", {})
load("target/hodge-data/hodge_data_g4.sage")
hodge_number[(4, 3, 2)]  # 1/87091200
```

This matches the `hodge_number` table format documented by
[admcycles](https://gitlab.com/modulispaces/admcycles/-/blob/master/README.rst).
Loading additional genus files adds entries to the same table; existing keys
are replaced with the loaded values. The files contain only Hodge-number
assignments, so use them to extend a table rather than replacing admcycles'
complete bundled data file, which also contains Jacobian locus coefficients.
The SQLite database belongs to `sc-faber`; it is not an admcycles cache.
Loading the Sage table does not automatically seed the separate
[`hodge_recursion` cache](https://gitlab.com/modulispaces/admcycles/-/blob/master/admcycles/hodge_recursion.py).

SQLite stores completed integral values and the library's persistent intermediate
results. Reuse the same data directory to resume: completed genus files are
skipped, and an unfinished genus is re-enumerated using committed cached values.
Unfinished query planning is rebuilt. A `.sage.part` file is renamed only after
the genus finishes and the cache and output have been flushed.

Small exact comparisons with the original evaluator run in the library's CI:

```sh
cargo test --release -p sc-faber hodge::lambda
```

## 中文

这个纯 Rust 示例计算给定亏格 g 的全部纯 lambda 顶次积分：单项式由
`lambda_1, ..., lambda_g` 组成，加权次数为 `3g-3`，支持 `2 <= g <= 20`。
结果包含零值，全部使用精确有理数；计算无需 Python 绑定或 Sage。

在仓库根目录运行上面的命令。默认只计算 g=2～4，线程上限为 24，进程内存
预算为 2 GiB；完整的 g≤20 计算需要显式指定范围及合适的资源预算。
`--genus G` 指定单个亏格，`--threads N`、`--memory-gib N` 分别指定后端线程
上限和进程内存预算。单项式按顺序枚举，各积分使用后端的并行计算能力。

示例直接调用 `sc_faber::hodge::lambda_integrals`，只处理参数、缓存会话、输出
和进度。分拆枚举、Mumford 关系、闭式与 GRR 后端全部由主库负责，并保留批量
计算中对同一亏格较早结果的复用，避免逐条独立计算造成回归。批量复用的数值
会保留到该亏格完成；回调逐条输出不代表内存用量恒定。小亏格下与原算法的精确
对照纳入主库 CI。

示例当前范围止于 g=20，更高亏格尚未验证。内存预算必须能容纳任务仍需常驻的
工作集；SQLite 持久化和工作块溢写不保证任意任务都能在极低预算下完成。

默认输出与缓存目录为 Git 忽略的 `target/hodge-data/`，可用 `--data-dir PATH`
修改。每个亏格输出 `hodge_data_gG.sage`，降序元组表示 lambda 下标列表；
上面的 Sage 例子展示了读取方式。

输出与 [admcycles 的 `hodge_number` 表](https://gitlab.com/modulispaces/admcycles/-/blob/master/README.rst)
格式兼容，可在同一 SageMath 会话中直接加载多个亏格文件，追加到已有表中；
同名键会被加载值覆盖。我们的文件只包含 Hodge 积分赋值；上游完整数据文件还
包含 Jacobian locus 系数，因此应追加加载数据。SQLite 属于 `sc-faber`，
加载 Sage 数据表也不会自动填充 admcycles 的 `hodge_recursion` 算法缓存。

使用相同目录重启即可复用已提交的 SQLite 缓存；已完成的亏格文件直接跳过，
未完成的亏格重新枚举，并复用已持久化的数值。未完成查询的工作图仍需重建。
每个亏格完成并刷新缓存、输出后，才将 `.sage.part` 重命名为正式结果文件。
