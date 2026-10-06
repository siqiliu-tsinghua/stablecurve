English | [中文](USER_GUIDE.zh-CN.md)

# stablecurve User Guide

Each engine ships a **CLI binary** (batch-friendly) and a **Python module**
(pyo3, returning `fractions.Fraction`). This is the complete usage
reference; for algorithms, verification and design decisions see the design
documents under `docs/`.

Prerequisites: a Rust toolchain (≥ 1.85; e.g. the `rustc` packaged in Debian 13
stable works) and [uv](https://docs.astral.sh/uv/).

- CLI binaries: `cargo build --release` produces `target/release/growi`,
  `faber` and `gwall`.
- Python modules: `uv sync` builds the three extensions `sc_growi`,
  `sc_faber`, `sc_gwall` and installs them (editable) into `.venv`; it
  rebuilds them automatically when the Rust sources change. To build
  standalone wheels: `uvx maturin build --release -m crates/sc-growi-py/Cargo.toml`
  (likewise `sc-faber-py`, `sc-gwall-py`).

All three CLIs speak the same batch protocol: one query per input line,
one JSON object per output line, `"line"` being the 1-based input line
number (blank lines are skipped).

---

## 1. sc-growi — GW invariants (GROWI port)

Computes Gromov-Witten invariants of P^n and hypersurfaces (including the
quintic Calabi–Yau threefold), at any genus, with relative conditions.

### 1.1 Expression syntax (direct port of the original `expr.cc`)

- Assignments: `N=3` (ambient dimension), `E=5` (hypersurface degree),
  `G=2` (genus), `D=5` (curve degree)
- Keywords: `rational`/`elliptic` (= G0/G1), `lines`/`conics`/`cubics`…
  (= D1/D2/D3…), `point`/`curve`/`surface`/`threefold` (or `3fold`)…,
  `quadric`/`cubic`/`quintic`… (= E2/E3/E5…), `P^2`/`P^3`/…
- Stop words (equivalent to commas): `in`/`on`/`to`/`thru`/`through`/
  `with`/`and`
- Class expressions: polynomials in `H` and `psi` with `+ * ^ ( )`;
  `H:9` = H repeated 9 times; `[class,multiplicity]` = a relative condition;
  `pc(k)` = the k-th primitive class of a hypersurface
- A class that is identically zero makes the invariant 0 (as in the
  original). Numbers are range-checked (N/E/G/D ≤ 32767, total degree of a
  class ≤ 255); out-of-range input is an error, never truncated.
- Examples: `lines in cubic surface`, `elliptic D=6 in P^2 thru H^2:18`,
  `G=2,D=5 in P^4 thru H^3*psi^22`, `conics in P^1 thru H^2:4,[1,2]`

### 1.2 CLI

```bash
./target/release/growi "lines in quintic threefold"        # 2875
./target/release/growi --batch --cache gw.sqlite < q.txt   # one expression per line
```

The cache defaults to `~/.growi/cache.sqlite`; `--no-cache` keeps it in
memory only. Batch responses:

```json
{"line":1,"status":"ok","value":"2875"}
{"line":2,"status":"unknown"}
{"line":3,"status":"invalid","message":"..."}
```

`unknown` means the strategy chain cannot compute the invariant (the
original's UNDEFINED; single-query mode prints `UNDEFINED (...)`).

### 1.3 Python

```python
from sc_growi import GROWI
g = GROWI()                                   # two-level cache at ~/.growi/cache.sqlite
g = GROWI(cache="gw.sqlite")                  # cache in a given file
g = GROWI(cache="")                           # in-memory only
g.gw("G=2,D=5 in P^4 thru H^3*psi^22")       # Fraction(-41369, 110075314176)
```

Unsupported queries (e.g. high-genus hypersurface invariants with E>1)
return `None`; invalid input raises `ValueError`. The GIL is released
during computation.

---

## 2. sc-faber — tautological-algebra integrals (Faber Numbers port + admcycles parity)

Computes integrals of tautological classes over M̄_{g,n}: mixed κ/ψ/λ/ch
monomials, DR cycles (Pixton/Chiodo, k-twist, spin), the Θ class, the Hain
divisor, hyperelliptic/bielliptic loci, Witten r-spin classes and the
generalized Hodge chern. All results are **exact rationals**.

### 2.1 CLI (JSONL protocol)

```bash
./target/release/faber '{"op":"klp","g":2,"n":2,"psi":[[1,2]],"lambda":[[1,3]]}'
./target/release/faber --batch --cache fb.sqlite < queries.jsonl
```

Responses: `{"line":k,"status":"ok","value":"p/q"}` or
`{"line":k,"status":"error","error":"..."}`. Common fields: `g` (genus),
`kappa` ([[index,power],…]), `psi` ([[leg,power],…]). ψ legs must be
markings 1..n (entries for the same leg are multiplied together); unknown
field names are rejected, so a misspelt field is reported rather than
ignored.

Operations (all take `g`; all return ∫class·X):

| op | own fields | meaning |
|---|---|---|
| `klp` | `n`; `lambda`/`ch` same format | ∫ of a κ^κ ψ^ψ λ^λ ch^ch monomial |
| `dr` | `a` (= A vector, ΣA = 0), `r` (defaults to g) | Pixton DR integral ∫DR(g,A,r)·X |
| `drpoly` | `a`, `k` | Chiodo interpolation main path (= DR_cycle chiodo_coeff=True) |
| `drspin` | `a`, `k` | spin DR (CSS21; A all odd, k odd) |
| `witten` | `a3`=[[C-num,C-den,D],…] | Witten r-spin limit class (A_i = C_i·r + D_i; ΣD = 2g−2) |
| `hyperell` / `biell` | `n` (fixed marks), `m` (conjugate pairs) | hyperelliptic / bielliptic locus integral |
| `grr` | `n`, `t` (ch degree), `l`, `dvec`=[d₁..dₙ], `a3`=[[h, s-mask, v],…] | generalized Hodge chern ch_t (PRvZ20); S encoded as bitmask Σ2^{i−1} |

For `drpoly`/`drspin`, ΣA = k(2g−2+n) is required; `k` may be omitted when
ΣA is divisible by 2g−2+n.

Examples:

```json
{"op":"dr","g":2,"a":[1,1,-2],"r":2,"kappa":[[1,1]],"psi":[[1,3]]}   -> 13/144
{"op":"hyperell","g":2,"n":0,"m":0,"kappa":[[3,1]],"psi":[]}          -> 1/1152
{"op":"biell","g":2,"n":0,"m":0,"kappa":[[2,1]],"psi":[]}             -> 1/32
{"op":"grr","g":2,"n":1,"t":1,"l":1,"dvec":[0],"kappa":[[1,1]],"psi":[[1,2]]} -> 5/576
```

### 2.2 Python

```python
from sc_faber import FABER
f = FABER()                                  # two-level cache at ~/.faber/cache.sqlite
f.klp(2, 2, psi=[(1,2)], lam=[(1,3)])        # Fraction(1, 1440)
f.klp(7, 0, lam=[(1,18)])                    # Fraction(32017001, 638512875)
f.dr(2, [1,1,-2], r=2, kappa=[(1,1)], psi=[(1,3)])   # Fraction(13, 144)
```

`FABER(cache="path")` uses the given file and `FABER(cache="")` disables
the persistent cache. The persistent cache is **process-wide**: creating a
new `FABER` replaces (or detaches) the cache of earlier instances, and
`cache_attached` reports the current state.

Note: the pyo3 layer currently exposes `klp`/`dr`; use the CLI for the
remaining ops.

---

## 3. sc-gwall — local-curve GW and BPS (GWall port)

Computes (disconnected) GW invariants of maps to a smooth curve C_g
(point/loop descendants, ch insertions) and local BPS numbers.

### 3.1 Insertion encoding

Corresponds to the v_i encoding of the original Maple program. In Python
an insertion is a tuple or list; in the CLI a JSON array, or equivalently
an object `{"t":kind,"k":K,"i":I}` (`{"t":"ldeg","v":V}` for `ldeg`).

| encoding | meaning |
|---|---|
| `("t1", K)` | descendant τ_K(1) |
| `("pt", K)` | descendant τ_K(pt) |
| `("a", K, I)` / `("b", K, I)` | loop descendants τ_K(a_I) / τ_K(b_I) |
| `("chO", K)` | ch_K(Rπ*f*O) |
| `("chL", K)` | ch_K(Rπ*f*L) (L of degree 2g−2 = K by default) |
| `("ldeg", V)` | change the degree of L to V (original `[[deg]]`) |

Loop insertions may be given in any order: loop descendants anticommute,
and the sign of reordering to the canonical order (all a-loops before all
b-loops, as the original program requires) is applied automatically.

### 3.2 CLI (JSONL protocol)

```bash
./target/release/gwall '{"op":"gw","g":2,"d":1,"insert":[["t1",0],["pt",1],["pt",0]]}'
./target/release/gwall --batch < queries.jsonl
```

| op | fields | meaning |
|---|---|---|
| `gw` | `g`, `d`, `insert`, optional `ldeg` | the invariant ⟨Π insert⟩ |
| `ztil` / `z` | `g`, `d`, `n` | t-series coefficients t^0..t^n (ztil includes constant maps; z does not) |
| `l` | `g`, `d`, `n` | degree-log connected series L_d (generating function of the paper's N) |
| `bps` | `g`, `d`, `h` | local BPS number n_d^h(g) (full GW inversion) |

`d` is required. Series ops return an array of coefficients, e.g.
`{"line":1,"status":"ok","value":["15/2","0","-5/4"]}`.

### 3.3 Python

```python
from sc_gwall import GWALL
w = GWALL()
w.bps(3, 2, 1)                               # Fraction(4) — BP paper table
w.bps(0, 1, 0)                               # Fraction(1) — g = 0 is supported
w.gw(2, 2, [])                               # Fraction(8) — Ztil constant term
w.l(2, 2, 2)                                 # [15/2, 0, -5/4] — L_d coefficients
w.z(2, 1, 4)                                 # degree-1 closed-form coefficients
```

All genera g ≥ 0 are supported. For degree d ≥ 2 the computation is
limited, as in the original, to t-truncation ≤ 10 ("branch number 10");
queries beyond it raise `ValueError` immediately. The GIL is released
during computation.

---

## 4. Cross-engine consistency (optional reading)

faber and growi serve as oracles for each other on their shared
mathematical surface (ψ-integrals over M̄_{0,n}):
`growi "G=0,D=0 in P^d thru H:d,psi^a1,..,psi^am"` = `faber klp(0, d+m,
psi=[(1,a1),..])`, agreeing exactly on all 34 cases
(`tools/cross_growi_faber.py`); gwall's closed form for N₀ and faber's λ³
integrals agree at three independent points.
