# WIEN2k Parameter Optimizer

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

**Comprehensive, automatic, scientifically-documented optimization of all WIEN2k preprocessing parameters from a single `case.struct` file.**

---

## Overview

This script takes a WIEN2k `case.struct` file and automatically determines the optimal combination of **all** preprocessing parameters based on the official WIEN2k documentation, Peter Blaha's reference tables, and the four strict RMT conditions. It generates ready-to-use WIEN2k input files and a comprehensive scientific report.

### What It Optimizes

| Parameter | Description | Reference |
|-----------|-------------|-----------|
| **RMT** | Muffin-Tin radii (4 strict conditions) | [FAQ RMT](http://www.wien2k.at/reg_user/faq/rmt.html) |
| **RKMAX** | Plane-wave cutoff (Blaha table) | [FAQ RKMAX](http://www.wien2k.at/reg_user/faq/rkmax.html) |
| **GMAX** | Fourier expansion for density/potential | User's Guide §4.2 |
| **LMAX / LVNS** | Angular momentum cutoffs | User's Guide §4.4 |
| **k-mesh** | Adaptive Monkhorst-Pack grid | [FAQ kgen](http://www.wien2k.at/reg_user/faq/kgen.html) |
| **Mixing** | SCF convergence scheme (PRATT/MSR1a/MSEC1) | User's Guide §4.3 |
| **TEMP** | Fermi smearing for metallic systems | User's Guide §4.3 |
| **Core/Valence** | Ecut and HDLO recommendations | User's Guide §4.4 |

---

## Quick Start

### Requirements

- Python 3.8+
- No external Python dependencies (standard library only)
- WIEN2k (optional — only needed for `--auto-converge`)

### Basic Usage

```bash
# Simple parameter recommendation
python optimize_wien2k.py BaTiO3.struct

# High-precision relaxation
python optimize_wien2k.py BaTiO3.struct --calc-type relaxation --precision high

# Magnetic system
python optimize_wien2k.py Fe.struct --magnetic --precision very_high

# Skip generating input files (report only)
python optimize_wien2k.py BaTiO3.struct --no-input-files
```

### Auto-Convergence (requires WIEN2k)

```bash
# Full Blaha hierarchy: RMT → k-mesh → RKMAX
python optimize_wien2k.py BaTiO3.struct --auto-converge
```

## Generated Output

```
optim_results/
├── BaTiO3_optimization_report.txt   # 10-section scientific report
├── BaTiO3.struct_optimized           # Struct with optimized RMTs
├── BaTiO3.in0                        # GMAX, VXCTYPE, SPIN
├── BaTiO3.in1                        # RKMAX, LMAX, E-params
├── BaTiO3.in2                        # Mixing, TEMP, convergence criteria
├── BaTiO3.inm                        # Magnetization settings
└── BaTiO3.klist                      # k-point mesh
```

---

## The Four RMT Conditions

The script enforces four strict conditions for Muffin-Tin radii optimization:

### Condition 1: Non-Overlapping Spheres
RMT<sub>i</sub> + RMT<sub>j</sub> ≤ safety_margin × NN<sub>ij</sub>

- **SCF**: safety_margin = 0.90
- **Relaxation/Optimization**: safety_margin = 0.85

> *"You can save a lot of CPU-time by changing RMT to almost touching spheres."* — Blaha

### Condition 2: Core Leakage Control (`:NEC01`)
Core charge must remain inside the muffin-tin sphere:

- **Acceptable**: leakage < 0.001 e⁻
- **Warning**: 0.001 < leakage < 0.002 e⁻
- **Critical**: leakage > 0.01 e⁻ → must fix

### Condition 3: Relative Size Ratio
max(RMT) / min(RMT) ≤ 1.5

- d-elements: 20% larger than sp
- f-elements: 30% larger than sp
- Ratio > 1.5 → ghostbands and QTL-B error

### Condition 4: Structural Change Margin
For relaxation/EOS calculations, RMTs are reduced by 7% to allow atomic motion.

---

## The Blaha Convergence Hierarchy

When using `--auto-converge`, the script follows Peter Blaha's recommended order:

```
┌──────────────────────────────────────┐
│ 1. RMT (geometric + core leakage)    │
│    • 4 strict conditions             │
│    • Run SCF, parse :NEC01           │
│    • Iterate until leakage < 0.002   │
└──────────────────┬───────────────────┘
                   ▼
┌──────────────────────────────────────┐
│ 2. k-mesh convergence                │
│    • Start with mesh/2               │
│    • Increase until ΔE < 0.1 mRy     │
└──────────────────┬───────────────────┘
                   ▼
┌──────────────────────────────────────┐
│ 3. RKMAX convergence                 │
│    • Start with RKMAX−1.0            │
│    • Increase 0.5 until ΔE < 0.1 mRy │
└──────────────────┬───────────────────┘
                   ▼
              Final Report
```

---

## CLI Options

| Flag | Values | Default | Description |
|------|--------|---------|-------------|
| `--calc-type` | `scf`, `relaxation`, `optimization`, `eos`, `forces`, `efg` | `scf` | Calculation type |
| `--precision` | `screening`, `coarse`, `medium`, `high`, `very_high` | `medium` | RKMAX & Ecut precision |
| `--refinement` | `coarse`, `medium`, `fine`, `very_fine` | `medium` | k-mesh refinement factor |
| `--system-type` | `metal_small`, `semiconductor`, `insulator`, `metal_large`, `insulator_large`, `surface`, `molecule` | auto | Override system detection |
| `--vxc` | `pbe`, `lda`, `wc`, `pbesol`, `scan`, `hse` | `pbe` | Exchange-correlation functional |
| `--magnetic` | flag | off | Spin-polarized calculation |
| `--auto-converge` | flag | off | Run full convergence hierarchy |
| `--output` | path | `./optim_results` | Output directory |
| `--no-input-files` | flag | off | Skip generating WIEN2k input files |
| `--quiet` | flag | off | Suppress console output |

---

## System Type Detection

The script automatically classifies the system based on element composition:

| Has TM | Has O/N/F/S | Atoms | Classification | k-density |
|--------|-------------|-------|----------------|-----------|
| ✗ | ✗ | ≤10 | semiconductor | 500 /bohr³ |
| ✓ | ✓ | ≤20 | semiconductor | 500 /bohr³ |
| ✓ | ✗ | ≤10 | metal_small | 5000 /bohr³ |
| any | any | ≥40 | insulator_large | 10 /bohr³ |

TM = transition metal (d or f block)

Hexagonal cells (γ = 120°) are automatically detected and use Gamma-centered k-meshes.

---

## Element-Specific Rules

### Hydrogen
- Short bonds (C–H, O–H): RMT(H) ≈ 0.5 × RMT(partner)
- GMAX increased to 20 for density expansion
- RKMAX capped at 3.5 for very small H spheres

### Same Element
- All atoms of the same element get identical RMT values
- `setrmt_lapw` behavior replicated automatically

### Large RMT (> 2.2 bohr)
- HDLOs (High-Derivative Local Orbitals) recommended for d/f elements
- LVNS increased to 6 (RMT > 2.2) or 8 (RMT > 2.5)

### Phase Comparison
For cohesive energy calculations between phases, ensure identical effective RKMAX across phases:
```
RKMAX_eff = RKMAX_input × RMT / R_min
```

---

## Mixing Schemes

| System | Scheme | Mixing Factor | TEMP (Ry) |
|--------|--------|---------------|-----------|
| Insulator | PRATT | 0.30 | 0.0001 |
| Semiconductor | PRATT | 0.25 | 0.001 |
| Metal | MSR1a | 0.40 | 0.003 |
| Magnetic Metal | MSEC1 | 0.30 | 0.002 |

---

## Project Structure

```
wien2k-parameter-optimizer/
├── optimize_wien2k.py              # Main CLI entry point
├── optim_wien/                     # Core package
│   ├── __init__.py
│   ├── constants.py                # Element tables, reference data
│   ├── struct_parser.py            # case.struct parser + NN computation
│   ├── rmt.py                      # 4-condition RMT optimization
│   ├── rkmax.py                    # Blaha RKMAX reference table
│   ├── gmax.py                     # GMAX for density expansion
│   ├── lmax.py                     # LMAX/LVNS/HDLO optimization
│   ├── kmesh.py                    # Adaptive k-mesh generation
│   ├── mixing.py                   # Mixing & TEMP optimization
│   ├── core_valence.py             # Core/valence separation
│   ├── input_generator.py          # WIEN2k file generation
│   ├── convergence.py              # Auto-convergence engine
│   └── report.py                   # Comprehensive report generation
├── BaTiO3.struct                   # Example structure
├── Si.struct                       # Example structure
└── README.md
```

---

## Scientific References

| Source | Link |
|--------|------|
| P. Blaha et al., *J. Chem. Phys.* **152**, 074101 (2020) | [DOI](https://doi.org/10.1063/1.5143061) |
| WIEN2k User's Guide | [usersguide.pdf](http://susi.theochem.tuwien.ac.at/reg_user/textbooks/usersguide.pdf) |
| WIEN2k FAQ — RMT | [rmt.html](http://www.wien2k.at/reg_user/faq/rmt.html) |
| WIEN2k FAQ — RKMAX | [rkmax.html](http://www.wien2k.at/reg_user/faq/rkmax.html) |
| WIEN2k FAQ — kgen | [kgen.html](http://www.wien2k.at/reg_user/faq/kgen.html) |
| L.D. Marks, Optimization Notes | [WIEN2k Workshop 2006](http://www.wien2k.at/events/ws2006/Exercises.pdf) |
| P. Blaha, Lecture Slides | [WIEN2k Workshop 2015](http://www.wien2k.at/events/ws2015/WIEN2k-Blaha-lectures.pdf) |
| P. Blaha, Getting Started | [WIEN2k Workshop 2024](http://susi.theochem.tuwien.ac.at/events/ws24/Blaha-getting_started.pdf) |
| WIEN2k Mailing List Archive | [mail-archive](https://www.mail-archive.com/wien@zeus.theochem.tuwien.ac.at/) |

---

## Limitations

1. **NN distances for MULT > 1**: Nearest-neighbor computation is exact only when each non-equivalent atom has MULT=1. For space groups with MULT > 1 positions (e.g., diamond Si), NN distances are approximate. Use WIEN2k's `x nn` for precise values.

2. **Core leakage requires SCF**: The `:NEC01` core leakage check in auto-converge mode requires running a full SCF cycle (WIEN2k must be installed).

3. **No spin-orbit coupling**: `.inso` file generation and spin-orbit specific parameter adjustments are not yet supported.

4. **GGA+U / mBJ**: `.inorb` and `.in0_grr` file generation for GGA+U and modified Becke-Johnson potentials is not yet supported.

---

## License

MIT License — See [LICENSE](LICENSE) file.

---

## Citation

If you use this tool in your research, please cite:

```
P. Blaha, K. Schwarz, F. Tran, R. Laskowski, G. K. H. Madsen, and L. D. Marks,
J. Chem. Phys. 152, 074101 (2020).
```

And mention this optimizer:

```
WIEN2k Parameter Optimizer, https://github.com/majidhaddad02/wien2k-parameter-optimizer
```

---

*Based on the official WIEN2k documentation by Peter Blaha, Karlheinz Schwarz, and the WIEN2k team at TU Wien.*
