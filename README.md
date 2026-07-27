# WIEN2k Parameter Optimizer

**Author:** Dr. Majid Haddad — [dr.majidhaddad@gmail.com](mailto:dr.majidhaddad@gmail.com)

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

**Comprehensive, automatic, scientifically-documented optimization of all WIEN2k preprocessing parameters from a single `case.struct` file.**

---

## Quick Start

```bash
./install.sh --here && export PATH="$(pwd):$PATH"
opt_wien2k                        # auto-detects *.struct in current dir
opt_wien2k -i                     # interactive wizard
opt_wien2k BaTiO3.struct --precision high
```

📖 **Full documentation:** [USER_GUIDE.md](USER_GUIDE.md)

---

## Features

| Parameter | Optimization Method |
|-----------|-------------------|
| **RMT** | 4 strict Blaha conditions: non-overlap, core leakage, ratio ≤ 1.5, structural margin |
| **RKMAX** | Element-specific reference table; effective RKMAX per atom |
| **GMAX** | Adaptive — depends on H/Li/halogens/f-elements |
| **LMAX/LVNS** | Dynamic per atom type (sp=6, d=10, f=12) |
| **k-mesh** | Correct Monkhorst-Pack formula; hexagonal/gamma detection; magnetic ×1.5 density |
| **Mixing** | PRATT (insulators) / MSR1a (metals) / MSEC1 (large metals) |
| **TEMP** | Adaptive Fermi smearing (0.001–0.002 Ry) |
| **Core/Valence** | Ecut per precision level; HDLO for large RMT |

All WIEN2k input files generated: `in0`, `in1`, `in2`, `inm`, `klist`, `struct_optimized`

Optional auto-convergence: runs WIEN2k SCF cycles to converge RMT (via `:NEC01`), k-mesh, and RKMAX following Blaha's hierarchy.

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

# Full convergence with WIEN2k
opt_wien2k BaTiO3.struct --auto-converge --precision high

# Show all options
opt_wien2k --help
```

---

## The Four RMT Conditions (Blaha)

The optimizer enforces all four RMT constraints from the WIEN2k FAQ:

1. **Non-overlap:** `RMT(i) + RMT(j) ≤ 0.90 × NN_distance`
2. **Core leakage:** `:NEC01 < 0.002` (verified via SCF in auto-converge mode)
3. **Ratio balance:** 1.3 for sp-systems, 1.5 for mixed sp/d/f
4. **Structural margin:** extra 10–15% reduction with free internal coordinates

No existing tool enforces all four simultaneously.

---

## Installation

```bash
./install.sh --here        # alias in current directory (recommended)
./install.sh --user        # pip install to ~/.local
sudo ./install.sh --system # system-wide
./install.sh --uninstall   # remove
```

Zero external dependencies — Python 3.8+ standard library only.

---

## RMT / RKMAX Reference Table

| Element Type | RMT (bohr) | RKMAX | Notes |
|-------------|-----------|-------|-------|
| H, He | 0.80–1.00 | 3.5 | Very small spheres |
| Li, Be, B | 1.60–1.80 | 4.5 | Small atoms |
| C, N, O, F | 1.00–1.40 | 6.5 | Electronegative |
| Si, P, S | 1.80–2.10 | 5.0 | sp-elements |
| Na–Ar | 1.80–2.20 | 5.5 | Alkali |
| K, Ca | 2.20–2.50 | 6.0 | Large sp |
| 3d TM | 1.90–2.20 | 7.0–8.0 | Transition metals |
| 4d TM | 2.10–2.40 | 8.0 | Heavy TM |
| 5d TM | 2.20–2.60 | 8.5 | Very heavy |
| Lanthanides | 2.30–2.70 | 8.0 | f-elements |
| Actinides | 2.40–2.80 | 8.5 | 5f-elements |

---

## Convergence Hierarchy (Blaha)

```
RMT core leakage (< 0.002 e⁻)
    ↓
k-mesh convergence (ΔE < 0.1 mRy)
    ↓
RKMAX convergence (ΔE < 0.1 mRy)
    ↓
Forces convergence (ΔF < 1 mRy/bohr)
    ↓
Property convergence (gap, EFG, DOS, optics)
```

**Key principle:** total energy converges faster than forces. Always converge the **property of interest**, not just energy.

---

## Project Structure

```
wien2k-parameter-optimizer/
├── optimize_wien2k.py          # Main CLI entry point
├── optim_wien/                 # Package
│   ├── cli.py                  # Interactive TUI (menus, colors, progress)
│   ├── constants.py            # All reference tables & constants
│   ├── struct_parser.py        # WIEN2k struct file parser + NN distance
│   ├── rmt.py                  # 4-condition RMT optimization
│   ├── rkmax.py                # Blaha RKMAX table lookup
│   ├── gmax.py                 # GMAX optimization
│   ├── lmax.py                 # LMAX / LVNS / HDLO
│   ├── kmesh.py                # Adaptive k-mesh (incl. system detection)
│   ├── mixing.py               # SCF mixing scheme & TEMP
│   ├── core_valence.py         # Ecut & HDLO recommendations
│   ├── input_generator.py      # WIEN2k file generation (in0, in1, in2, inm, klist)
│   ├── convergence.py          # Auto-convergence engine
│   └── report.py               # Scientific report generation
├── install.sh                  # Multi-mode installer
├── setup.py                    # pip package setup
├── USER_GUIDE.md               # Comprehensive user guide
└── README.md                   # This file
```

---

## Requirements

- **Python 3.8+** (stdlib only — no pip packages)
- **WIEN2k** (only for `--auto-converge`)

---

## References

- Blaha et al., *J. Chem. Phys.* **152**, 074101 (2020)
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

Based on the work of P. Blaha, K. Schwarz, and the WIEN2k development team at TU Wien.
