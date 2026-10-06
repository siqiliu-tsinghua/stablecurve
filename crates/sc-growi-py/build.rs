// macOS 上用纯 cargo 构建扩展模块时需要 `-undefined dynamic_lookup` 链接参数
// （maturin 会自动加；这里让 `cargo build --workspace` 在各平台都能链接）。
fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
