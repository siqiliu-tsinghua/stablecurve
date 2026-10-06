English | [中文](README.zh-CN.md)

# stablecurve

A Rust workspace of three computational algebraic-geometry engines sharing a
common rational-number layer (`sc-num`):

| crate | contents | origin |
|---|---|---|
| `sc-growi` | Gromov-Witten invariants (P^n at any genus, hypersurfaces, relative invariants) | Andreas Gathmann, *GROWI* 1.0.4 (C++) |
| `sc-faber` | Tautological algebra of the moduli spaces M̄_{g,n}: ψ/κ/λ/ch intersection numbers, DR cycles (Pixton/Chiodo, k-twist, spin), Θ class, Hain divisor, Hyperell/Biell pairings, Witten r-spin, generalized Hodge chern (GRRcomp) | Carel Faber, *Numbers* (Maple) + admcycles feature parity |
| `sc-gwall` | Local GW invariants of a smooth curve (point/loop descendants, ch insertions) and local BPS numbers (all 121 numeric cells of the BP paper tables, 4 correcting paper misprints; any genus g ≥ 0) | Jim Bryan, *GWall* (Maple) + the BP paper |
| `sc-num` | malachite rational wrapper (shared by all three engines) | — |

- **Usage (CLI + Python)**: [USER_GUIDE.md](USER_GUIDE.md)
- **Analysis & design**: [docs/GROWI_REDESIGN.md](docs/GROWI_REDESIGN.md),
  [docs/FABER_REDESIGN.md](docs/FABER_REDESIGN.md),
  [docs/GWALL_REDESIGN.md](docs/GWALL_REDESIGN.md)

## Quick start

```bash
uv sync                          # builds and installs the Python modules sc_growi / sc_faber / sc_gwall
cargo build --release            # CLI binaries target/release/{growi,faber,gwall}
uv run pytest                    # 21 Python tests
# The Rust workspace also builds the pyo3 crates, which need a Python >= 3.12:
export PYO3_PYTHON="$PWD/.venv/Scripts/python.exe"   # Windows; elsewhere .venv/bin/python
cargo test --workspace --release # 164 Rust tests
```

All three engines have passed their "golden-number" acceptance: GROWI
against every one of the 41 assertions of gwcheck.cm (cold start 4.3 s,
about 40× the original); faber against Faber's *Instructions* examples and
admcycles differentials (all seven parity items completed); gwall against
all 121 numeric cells of the local BPS tables (n₂–n₅) in the
Bryan–Pandharipande paper, 4 of which correct misprints in the paper, and
against the paper's closed results for g = 0 and g = 1.
See the acceptance records in each design document for details.

**Tested platforms (v0.1.0).** The full Rust test suite, the Python tests
(both via wheel installation and via `uv sync`), the gwcheck acceptance
suite and the growi↔faber cross-check all pass on:

| platform | Rust | Python |
|---|---|---|
| Windows 11, x86_64 | 1.99 | 3.12 |
| Debian 13 (Linux), x86_64 | 1.85.1 (Debian's packaged `rustc`) | 3.12, 3.13 |
| macOS 26.7, x86_64 (Intel) | 1.85.1 and 1.99 | 3.12 |

Apple Silicon (arm64) has not been tested. The project has no CI; these
results were obtained manually before release.

## License

Copyright (C) 2026 Si-Qi Liu. This project is free software, licensed
under the [GNU General Public License, version 3 or later](LICENSE)
(GPL-3.0-or-later): you may use, modify and redistribute it, provided
derived works are distributed under the same license with their source
code.

`sc-growi` is a translation of GROWI (GPL) and parts of `sc-faber` follow
admcycles (GPL-2.0-or-later). The κ/ψ/λ/ch kernel of `sc-faber` is ported
from C. Faber's *Numbers* and `sc-gwall` from J. Bryan's *GWall*, both
research programs distributed publicly by their authors. See
[NOTICE](NOTICE) for the upstream works and attributions.

## Acknowledgements

This project was completed with the assistance of GLM-5.3, GPT-6-Astra and
Claude Opus 5.5.

## References

- A. Gathmann, *GROWI 1.0.4* (original program and gwcheck.cm):
  https://personal.math.uni-sb.de/~agathmann/growi/
- C. Faber, *Numbers* (Maple program and Instructions):
  https://math.stanford.edu/~vakil/programs/Numbers.tar
- J. Bryan, *GWall* (original program):
  https://personal.math.ubc.ca/~jbryan/gwall.html ;
  J. Bryan, R. Pandharipande, *BPS states of curves in Calabi–Yau
  threefolds*: https://arxiv.org/abs/math/0009025
- admcycles (Sage library, the differential-test oracle for faber), by
  Aaron Pixton, Johannes Schmitt, Vincent Delecroix, Jason van Zelm,
  Jonathan Zachhuber and Johannes Schwab:
  https://gitlab.com/modulispaces/admcycles
