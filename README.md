# WIEN2k Parameter Optimizer

**Author:** Majid Haddad — majidhaddad@gmail.com

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

**Comprehensive, automatic, scientifically-documented optimization of all WIEN2k preprocessing parameters from a single `case.struct` file — now with convergence-verified Aitken Δ² extrapolation.**

---

## Quick Start

```bash
# Install
pip install -e /path/to/wien2k-parameter-optimizer

# Run
opt_wien2k                                    # auto-detects *.struct in current dir
opt_wien2k -i                                 # interactive wizard
opt_wien2k BaTiO3.struct --precision high
opt_wien2k case.struct --converge rkmax kmesh  # convergence-verified (Aitken)
```

Full documentation: [USER_GUIDE.md](USER_GUIDE.md)

---

## Features

| Parameter | Optimization Method |
|-----------|-------------------|
| **RMT** | 4 strict conditions: non-overlap (0.90 SCF / 0.85 relax), core leakage, ratio ≤ 1.5, structural margin |
| **RKMAX** | Blaha 2020 Table I for C/N/O/Fe/Cu; all others default to 7.0 (WIEN2k default). Convergence-verified via Aitken Δ² extrapolation |
| **GMAX** | Adaptive per precision level; boosted for H (<0.8), Li (<1.2), halogens (<2.0), f-elements (<2.5) |
| **LMAX/LVNS** | Dynamic per atom: sp=6, d=10, f=12. LVNS: 4/6/8 |
| **k-mesh** | Adaptive Monkhorst-Pack; hexagonal gamma-detection; system-type-aware density; bandgap override |
| **Mixing** | PRATT (insulator 0.30, semiconductor 0.25) / MSR1a (metal 0.20) / MSEC1 (magnetic 0.15) |
| **TEMP** | Fermi smearing: insulator 0.0001, semiconductor 0.001, metal 0.001, magnetic 0.002 Ry |
| **Core/Valence** | Ecut per precision level (-4 to -8 Ry); HDLO for RMT > 2.2 (d/f) or > 2.5 (any) |
| **HSE hybrid** | in2c file generation with FOCK operator parameters |

All WIEN2k input files generated: `in0`, `in1`, `in1c`, `in2`, `in2c` (HSE), `inm`, `klist`, `machines`, `struct_optimized`

### New: Submit Script Generation

```bash
# Generate SLURM script
opt_wien2k Si.struct --machines 16 --submit slurm

# Generate PBS script
opt_wien2k Fe.struct --magnetic --machines 32 --submit pbs
```

### New: Strict-FAQ Mode

```bash
# No precision offsets — uses table values directly (Blaha defaults)
opt_wien2k Si.struct --strict-faq

# With convergence verification
opt_wien2k Si.struct --strict-faq --converge rkmax kmesh
```

### Convergence-Verified Optimization with Aitken Δ² Extrapolation

The `--converge` flag runs **real SCF cycles** and uses Aitken's acceleration method to extrapolate converged parameters with minimal computational cost — typically **~2× fewer SCF runs** than brute-force linear sweeps.

```
E_inf = E3 - (E3 - E2)² / (E3 - 2·E2 + E1)      [Aitken, 1926]
```

- **3 seed points** (loose criteria) → extrapolate → **1 confirm run** (production criteria)
- Non-monotonic fallback: 4+ point exponential least-squares regression (pure Python, no scipy)
- Separately tunable tolerance: `--etol 1.0` (mRy/atom) — auto-tightened 10× for forces
- Full markdown convergence report with raw data, fit parameters, and confirmation residuals

```bash
opt_wien2k BaTiO3.struct --converge rkmax kmesh --etol 1.0
opt_wien2k Si.struct --calc-type forces --converge rkmax kmesh    # auto-tightens etol
opt_wien2k Fe.struct --magnetic --converge rkmax kmesh --cluster-submit
```

---

## Usage Examples

```bash
# Basic — auto-detect struct file
opt_wien2k

# Single step only
opt_wien2k BaTiO3.struct --only rmt

# Magnetic metal
opt_wien2k Fe.struct --magnetic --precision high

# Publication quality
opt_wien2k BaTiO3.struct --precision very_high --refinement fine

# Strict FAQ mode (no precision offsets)
opt_wien2k Si.struct --strict-faq

# Interactive wizard
opt_wien2k -i

# Generate SLURM submit script with machines file
opt_wien2k Si.struct --machines 16 --submit slurm

# HSE hybrid functional
opt_wien2k Si.struct --vxc hse

# Convergence-verified (Aitken — recommended)
opt_wien2k BaTiO3.struct --converge rkmax kmesh --etol 1.0

# Band gap override
opt_wien2k case.struct --bandgap 0       # force metal treatment
opt_wien2k case.struct --bandgap 1.5     # force insulator treatment

# Override default init_lapw -prec flag
opt_wien2k Si.struct --prec 2

# Show all options
opt_wien2k --help
```

---

## The Four RMT Conditions (Blaha)

The optimizer enforces all four RMT constraints from the WIEN2k FAQ:

1. **Non-overlap:** `RMT(i) + RMT(j) ≤ 0.90 × NN_distance` (SCF) or `≤ 0.85` (relaxation)
2. **Core leakage:** `:NEC01` checked from SCF output; critical > 0.01 e⁻, warning > 0.002 e⁻
3. **Ratio balance:** max(RMT)/min(RMT) ≤ 1.5; strict ≤ 1.3 for sp-only systems
4. **Structural margin:** RMT reduced by 7% for relaxation/optimization/EOS calculations

Nearest-neighbour distances are calculated using all symmetry-equivalent positions (including MULT>1 sites in real WIEN2k struct files).

---

## RKMAX Reference Table

RKMAX table values and their provenance:

| Element | RKMAX | Source |
|---------|-------|--------|
| C | 5.5 | Blaha 2020, Table I |
| N | 6.5 | Blaha 2020, Table I |
| O | 7.0 | Blaha 2020, Table I |
| Fe | 8.0 | Blaha 2020, Table I |
| Cu | 8.0 | Blaha 2020, Table I |
| All others | 7.0 | WIEN2k default |

Effective RKMAX per atom: `RKMAX_eff(i) = RKMAX × RMT_min / RMT(i)`

---

## Precision Levels

| Level | Ecut (Ry) | GMAX | RKMAX Offset | init_lapw -prec |
|-------|----------|------|-------------|-----------------|
| screening | -4.0 | 10.0 | -1.0 | 0 |
| coarse | -5.0 | 12.0 | -0.5 | 0 |
| medium | -6.0 | 14.0 | 0.0 | 1 |
| high | -7.0 | 16.0 | +0.5 | 2 |
| very_high | -8.0 | 20.0 | +1.5 | 3 |

---

## Installation

```bash
# Editable install (recommended for development)
cd /path/to/wien2k-parameter-optimizer
pip install -e .

# System-wide
pip install /path/to/wien2k-parameter-optimizer
```

No external dependencies for core functionality — Python 3.8+ standard library only. Convergence modes require WIEN2k in PATH.

### Run tests

```bash
make test
# or
python3 -m pytest tests/ -v
```

---

## Project Structure

```
wien2k-parameter-optimizer/
├── optimize_wien2k.py          # Main CLI entry point
├── pyproject.toml              # Modern packaging config
├── Makefile                    # Install/test/clean targets
├── setup.py                    # Legacy pip setup (fallback)
├── optim_wien/                 # Package
│   ├── __init__.py             # Version
│   ├── cli.py                  # Interactive TUI (menus, colors, progress)
│   ├── constants.py            # All reference tables & constants
│   ├── struct_parser.py        # WIEN2k struct parser — MULT>1 equiv positions
│   ├── rmt.py                  # 4-condition RMT optimization
│   ├── rkmax.py                # Blaha RKMAX table lookup (seed for converge)
│   ├── gmax.py                 # GMAX optimization (halogen/f-element aware)
│   ├── lmax.py                 # LMAX / LVNS / HDLO
│   ├── kmesh.py                # Adaptive Monkhorst-Pack (bandgap override)
│   ├── mixing.py               # PRATT/MSR1a/MSEC1 mixing, TEMP, convergence
│   ├── core_valence.py         # Ecut & HDLO recommendations
│   ├── input_generator.py      # in0, in1, in1c, in2, in2c, inm, klist, machines
│   ├── submit.py               # SLURM & PBS script generators
│   ├── convergence.py          # Legacy auto-convergence (linear sweep)
│   ├── converge.py             # Aitken Δ² engine, SCF runner, markdown report
│   └── report.py               # Scientific report generation
├── tests/                      # Unit tests (31 tests, stdlib-only)
│   ├── test_converge.py        # Aitken extrapolation tests
│   ├── test_klist.py           # k-point list format tests
│   ├── test_kmesh.py           # k-mesh volume factor tests
│   └── test_struct_parser.py   # Struct parser tests
├── Si.struct                   # Example struct file
├── BaTiO3.struct               # Example struct file
├── USER_GUIDE.md               # Comprehensive user guide
└── README.md                   # This file
```

---

## Changelog (Bug Fixes)

| Issue | Before | After |
|-------|--------|-------|
| MULT>1 equivalent positions | Parser treated equiv position lines as elements | Skips equiv lines, stores all positions for NN min |
| numk encoding | `n1*100 + n2*10 + n3` (failed for multi-digit) | `n1*n2*n3` (correct) |
| in1 header format | 2 values | 3 values (RKMAX, global LMAX, V-NMT) |
| GMAX halogen check | Missing | Added: GMAX+2 for halogens (RMT<2.0), f-elements (<2.5) |
| GMAX method | Always used precision table | Strict-FAQ mode returns 12.0 (WIEN2k default) |
| NN fallback | `(i, i+1)` crash for last pair | Wrap-around check |
| subprocess calls | `shell=True` | `shlex.split` + list args |
| MSR1a mixing value | 0.40 | **0.20** (correct Blaha value) |
| MSEC1 mixing value | 0.30 | **0.15** (correct Blaha value) |
| klist format | Missing for auto-converge path | Full MP coords in all paths |
| Exponential fit sign | Wrong denominator sign | `1-exp(-α·h)` (correct) |
| RKMAX_TABLE | ~80 fabricated values | Only 5 from Blaha 2020 Table I |
| ELECTRONEGATIVITY | 350-line dead dict | Removed |
| RMT overwrite | Silent when struct RMT < 1.0 | Warning message |
| Same-element RMT | Forced equal (breaks surfaces) | Warning + manual review |
| Interactive config | Missing 8 keys | All keys present |
| klist inversion flag | Hardcoded `add INV=1` | Reads actual value |
| `run_lapw` in help | Missing `-p` flag | `run_lapw -p` |

---

## Requirements

- **Python 3.8+** (stdlib only — no pip packages required for core)
- **WIEN2k** (only for `--converge` or `--auto-converge` modes)
- **pytest** (only for running unit tests)

---

## References

- Blaha et al., *J. Chem. Phys.* **152**, 074101 (2020)
- A.C. Aitken, *Proc. Royal Soc. Edinburgh* **46**, 289–305 (1926)
- WIEN2k User's Guide & FAQ: http://www.wien2k.at/
- D.J. Singh, L. Nordstrom, *Planewaves, Pseudopotentials, and the LAPW Method*, 2nd Ed., Springer (2006)

---

## License

MIT — see LICENSE file.

## Citation

If you use this tool in your research:

> M. Haddad, *WIEN2k Parameter Optimizer*, https://github.com/majidhaddad02/wien2k-parameter-optimizer (2024)

---

## Acknowledgments

Based on the work of P. Blaha, K. Schwarz, and the WIEN2k development team at TU Wien. Aitken extrapolation method from A.C. Aitken (1926).
