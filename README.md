# WIEN2k Parameter Optimizer

**Author:** Dr. Majid Haddad — [dr.majidhaddad@gmail.com](mailto:dr.majidhaddad@gmail.com)

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

**Comprehensive, automatic, scientifically-documented optimization of all WIEN2k preprocessing parameters from a single `case.struct` file — now with convergence-verified Aitken Δ² extrapolation.**

---

## Quick Start

```bash
./install.sh --here && export PATH="$(pwd):$PATH"
opt_wien2k                                    # auto-detects *.struct in current dir
opt_wien2k -i                                 # interactive wizard
opt_wien2k BaTiO3.struct --precision high
opt_wien2k case.struct --converge rkmax,kmesh # convergence-verified (Aitken)
```

📖 **Full documentation:** [USER_GUIDE.md](USER_GUIDE.md)

---

## Features

| Parameter | Optimization Method |
|-----------|-------------------|
| **RMT** | 4 strict Blaha conditions: non-overlap, core leakage, ratio ≤ 1.5, structural margin |
| **RKMAX** | Element-specific reference table; convergence-verified via Aitken Δ² extrapolation |
| **GMAX** | Adaptive — depends on H/Li/halogens/f-elements; verification via lapw2 output parsing |
| **LMAX/LVNS** | Dynamic per atom type (sp=6, d=10, f=12) |
| **k-mesh** | Correct Monkhorst-Pack formula; hexagonal/gamma detection; metallicity-aware |
| **Mixing** | PRATT (insulators, 0.25) / MSR1a (metals, 0.20) / MSEC1 (large metals, 0.15) |
| **TEMP** | Adaptive Fermi smearing (0.001–0.002 Ry) |
| **Core/Valence** | Ecut per precision level; HDLO for large RMT |

All WIEN2k input files generated: `in0`, `in1`, `in1c`, `in2`, `inm`, `klist`, `struct_optimized`

### New: Convergence-Verified Optimization with Aitken Δ² Extrapolation

The `--converge` flag runs **real SCF cycles** and uses Aitken's acceleration method to extrapolate converged parameters with minimal computational cost — typically **~2× fewer SCF runs** than brute-force linear sweeps.

```
E_inf = E3 - (E3 - E2)² / (E3 - 2·E2 + E1)      [Aitken, 1926]
```

- **3 seed points** (loose criteria) → extrapolate → **1 confirm run** (production criteria)
- Non-monotonic fallback: 4+ point exponential least-squares regression (pure Python, no scipy)
- Separately tunable tolerance: `--etol 1.0` (mRy/atom) — auto-tightened 10× for forces
- Full markdown convergence report with raw data, fit parameters, and confirmation residuals

```bash
opt_wien2k BaTiO3.struct --converge rkmax,kmesh --etol 1.0
opt_wien2k Si.struct --calc-type forces --converge rkmax,kmesh    # auto-tightens etol
opt_wien2k Fe.struct --magnetic --converge rkmax,kmesh --cluster-submit
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

# Interactive wizard
opt_wien2k -i

# Legacy full convergence with WIEN2k
opt_wien2k BaTiO3.struct --auto-converge --precision high

# Convergence-verified (Aitken Δ² — recommended)
opt_wien2k BaTiO3.struct --converge rkmax,kmesh --etol 1.0

# Band gap override
opt_wien2k case.struct --bandgap 0       # force metal treatment
opt_wien2k case.struct --bandgap 1.5     # force insulator treatment

# Show all options
opt_wien2k --help
```

---

## The Four RMT Conditions (Blaha)

The optimizer enforces all four RMT constraints from the WIEN2k FAQ:

1. **Non-overlap:** `RMT(i) + RMT(j) ≤ 0.90 × NN_distance` — computed from all equivalent positions
2. **Core leakage:** `:NEC01 < 0.002` (verified via SCF in converge modes)
3. **Ratio balance:** 1.3 for sp-systems, 1.5 for mixed sp/d/f
4. **Structural margin:** extra 10–15% reduction with free internal coordinates

No existing tool enforces all four simultaneously. Nearest-neighbour distances are calculated using all symmetry-equivalent positions (including MULT>1 sites in real WIEN2k struct files).

---

## RKMAX Reference Table (Seeds for Convergence)

| Element Type | RKMAX | Notes |
|-------------|-------|-------|
| H, He | 3.5 | Very small spheres |
| Li, Be, B | 4.5 | Small atoms |
| C | 6.5 | sp light |
| N, S | 6.5 | per Blaha 2020 |
| O, F, Cl | 7.0 | Electronegative |
| Si, P | 5.0 | sp-elements |
| Na–Ar | 5.5 | Alkali |
| K, Ca | 6.0 | Large sp |
| 3d TM (Sc–Zn) | 7.5–8.0 | Transition metals |
| 4d TM (Y–Cd) | 8.0 | Heavy TM |
| 5d TM (Hf–Hg) | 8.0–8.5 | Very heavy |
| Lanthanides | 8.0–8.5 | f-elements |
| Actinides | 8.5 | 5f-elements |

*These are seeds for the convergence engine — never the final answer. Every converged value is traceable to `.scf` output.*

---

## Convergence Hierarchy

**Table-based (default):**
```
RMT (4 geometric constraints)
  → RKMAX (Blaha reference table)
    → GMAX (precision-dependent)
```

**Convergence-verified (`--converge`):**
```
RMT robustness (nudge ~5%, verify ΔE < etol)
  → k-mesh (3-point Aitken + confirmation)
    → RKMAX (3-point Aitken + confirmation)
      → GMAX verification (parse lapw2 output)
```

**Blaha's sequence:**
1. RMT robustness check
2. k-mesh convergence
3. RKMAX convergence
4. GMAX verification
5. Final confirmation

---

## Installation

```bash
./install.sh --here        # alias in current directory (recommended)
./install.sh --user        # pip install to ~/.local
sudo ./install.sh --system # system-wide
./install.sh --uninstall   # remove
```

Zero external dependencies — Python 3.8+ standard library only. Convergence modes require WIEN2k in PATH.

---

## Project Structure

```
wien2k-parameter-optimizer/
├── optimize_wien2k.py          # Main CLI entry point
├── optim_wien/                 # Package
│   ├── cli.py                  # Interactive TUI (menus, colors, progress)
│   ├── constants.py            # All reference tables & constants
│   ├── struct_parser.py        # WIEN2k struct parser — MULT>1 equiv positions
│   ├── rmt.py                  # 4-condition RMT optimization
│   ├── rkmax.py                # Blaha RKMAX table lookup (seed for converge)
│   ├── gmax.py                 # GMAX optimization (halogen/f-element aware)
│   ├── lmax.py                 # LMAX / LVNS / HDLO
│   ├── kmesh.py                # Adaptive Monkhorst-Pack (bandgap override)
│   ├── mixing.py               # MSR1a=0.20, MSEC1=0.15, TEMP
│   ├── core_valence.py         # Ecut & HDLO recommendations
│   ├── input_generator.py      # in0, in1, in1c, in2, inm, klist (full MP coords)
│   ├── convergence.py          # Legacy auto-convergence (linear sweep)
│   ├── converge.py             # NEW: Aitken Δ² engine, SCF runner, markdown report
│   └── report.py               # Scientific report generation
├── install.sh                  # Multi-mode installer
├── setup.py                    # pip package setup
├── USER_GUIDE.md               # Comprehensive user guide
└── README.md                   # This file
```

---

## Bug Fixes (Changelog)

| Issue | Before | After |
|-------|--------|-------|
| MULT>1 equivalent positions | Parser treated equiv position lines as elements → wrong NN distances | Skips MULT-1 lines, stores equiv positions, expands all for NN min |
| numk encoding | `n1*100 + n2*10 + n3` (failed for multi-digit like 12×12×12 → 1332) | `n1*n2*n3` (correct: 1728) |
| in1 header format | 2 values (RKMAX, V-NMT) | 3 values (RKMAX, global LMAX, V-NMT) |
| GMAX halogen check | Missing — halogens with RMT<2.0 not boosted | Added: GMAX+2 for halogens (RMT<2.0), f-elements (RMT<2.5) |
| NN fallback for last atom pair | `(i, i+1)` with i=last → crash | Added wrap-around check |
| subprocess calls | `shell=True` or string commands | `shlex.split` + list args + `shell=False` |
| MSR1a mixing value | 0.40 | **0.20** (correct Blaha value) |
| MSEC1 mixing value | 0.30 | **0.15** (correct Blaha value) |
| klist generation | Missing for legacy auto-converge path | Full Monkhorst-Pack coords + weights in all paths |
| Exponential fit sign | `exp(-α·h)-1` in denominator | `1-exp(-α·h)` (correct sign for positive α) |

---

## Requirements

- **Python 3.8+** (stdlib only — no pip packages required for core)
- **WIEN2k** (only for `--converge` or `--auto-converge` modes)
- **scipy** (optional — used only in 4-point exponential fit fallback; pure-Python fallback included)

---

## References

- Blaha et al., *J. Chem. Phys.* **152**, 074101 (2020)
- A.C. Aitken, *Proc. Royal Soc. Edinburgh* **46**, 289–305 (1926)
- WIEN2k User's Guide & FAQ: http://www.wien2k.at/
- D.J. Singh, L. Nordstrom, *Planewaves, Pseudopotentials, and the LAPW Method*, 2nd Ed., Springer (2006)

---

## License

MIT — see [LICENSE](LICENSE) file.

## Citation

If you use this tool in your research:

> M. Haddad, *WIEN2k Parameter Optimizer*, https://github.com/majidhaddad02/wien2k-parameter-optimizer (2024)

---

## Acknowledgments

Based on the work of P. Blaha, K. Schwarz, and the WIEN2k development team at TU Wien. Aitken extrapolation method from A.C. Aitken (1926).
