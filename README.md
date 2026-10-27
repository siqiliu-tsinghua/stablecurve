English | [中文](README.zh-CN.md)

# stablecurve

[![CI](https://github.com/siqiliu-tsinghua/stablecurve/actions/workflows/ci.yml/badge.svg)](https://github.com/siqiliu-tsinghua/stablecurve/actions/workflows/ci.yml)

A Rust workspace of three computational algebraic-geometry engines sharing a
common rational-number layer (`sc-num`):

| crate | contents | origin |
|---|---|---|
| `sc-growi` | Gromov-Witten invariants: P^n with no fixed genus bound, hypersurfaces (positive genus with E > 1 only partly supported, as in the original), relative invariants | Andreas Gathmann, *GROWI* 1.0.4 (C++) |
| `sc-faber` | Tautological algebra of the moduli spaces M̄_{g,n}: ψ/κ/λ/ch intersection numbers, DR cycles (Pixton/Chiodo, k-twist, spin), Θ class, Hain divisor, Hyperell/Biell pairings, Witten r-spin, generalized Hodge chern (GRRcomp) | Carel Faber, *Numbers* (Maple) + seven selected admcycles features |
| `sc-gwall` | Local GW invariants of a smooth curve (point/loop descendants, ch insertions) and local BPS numbers (all 121 numeric cells of the BP paper tables, 4 correcting paper misprints); every genus g ≥ 0, degree d ≥ 2 limited by the original's t-truncation ≤ 10 | Jim Bryan, *GWall* (Maple) + the BP paper |
| `sc-num` | malachite rational wrapper (shared by all three engines) | — |

- **Usage (CLI + Python)**: [USER_GUIDE.md](USER_GUIDE.md)
- **Analysis & design**: [docs/GROWI_REDESIGN.md](docs/GROWI_REDESIGN.md),
  [docs/FABER_REDESIGN.md](docs/FABER_REDESIGN.md),
  [docs/GWALL_REDESIGN.md](docs/GWALL_REDESIGN.md)

**Interfaces.** Each engine is a Rust library with a CLI and a Python
module; results are exact rationals (`Fraction` in Python). `sc-growi` and
`sc-gwall` offer the same functionality in all three. `sc-faber` covers
less outside Rust:

| sc-faber | Rust | CLI | Python |
|---|---|---|---|
| `klp` (κ/ψ/λ/ch monomials), `dr` (Pixton DR) | ✓ | ✓ | ✓ |
| `drpoly`, `drspin`, `witten`, `hyperell`, `biell`, `grr` | ✓ | ✓ | — |
| Θ class, Hain divisor | ✓ | — | — |

## Quick start

```bash
uv sync                          # builds and installs the Python modules sc_growi / sc_faber / sc_gwall
cargo build --release            # CLI binaries target/release/{growi,faber,gwall}
uv run pytest                    # Python tests
# The Rust workspace also builds the pyo3 crates, which need a Python >= 3.12:
export PYO3_PYTHON="$PWD/.venv/Scripts/python.exe"   # Windows; elsewhere .venv/bin/python
cargo test --workspace --release # Rust tests
```

All three engines have passed their "golden-number" acceptance: GROWI
against every one of the 41 assertions of gwcheck.cm; faber against Faber's
*Instructions* examples, admcycles differential tests of the seven selected
features, and all 53159 λ-monomial Hodge integrals on M̄_g (g ≤ 14) in the
admcycles data set; gwall against all 121 numeric cells of the local BPS
tables (n₂–n₅) in the Bryan–Pandharipande paper, 4 of which correct
misprints in the paper, and against the paper's closed results for g = 0
and g = 1. See the acceptance records in each design document for details.

**Performance.** The gwcheck suite takes 4.3 s from a cold cache on current
hardware. On the same machine, compute-bound workloads run 1.3–1.5× faster
than a strengthened build of the original GROWI (millisecond queries are
dominated by the fixed cost of the persistent cache), and genus 2 with
degree ≥ 7 and genus 3 are computed where that build crashes; see the
benchmarks in [docs/GROWI_REDESIGN.md](docs/GROWI_REDESIGN.md).

**Continuous integration.** Every push is tested on GitHub Actions
([`ci.yml`](.github/workflows/ci.yml)) on Linux (x86_64), macOS (Apple
Silicon) and Windows (x86_64): the full Rust test suite with Rust stable and
with the minimum supported version 1.85 (including the 2033 admcycles Hodge
integrals with g ≤ 9 and consistency checks between independent code
paths), rustfmt and clippy, the Python tests, the gwcheck acceptance suite
and the growi↔faber cross-check. The complete g ≤ 14 comparison with
admcycles is a one-off run of `tools/admcycles_hodge.py check --max-genus 14`
(about an hour). v0.1.0 was tested manually before release on Windows 11, Debian 13
(rustc 1.85.1) and macOS 26.7 (Intel).

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
