# WIEN2k Parameter Optimizer — User Guide

**Author:** Dr. Majid Haddad  
**Email:** dr.majidhaddad@gmail.com  
**License:** MIT  
**Repository:** https://github.com/majidhaddad02/wien2k-parameter-optimizer

---

## Table of Contents

1. [Overview](#1-overview)
2. [Quick Start](#2-quick-start)
3. [Installation](#3-installation)
4. [Basic Usage](#4-basic-usage)
5. [Interactive Mode](#5-interactive-mode)
6. [Command-Line Options](#6-command-line-options)
7. [Single-Step Optimization](#7-single-step-optimization)
8. [Parameter Reference](#8-parameter-reference)
9. [Calculation Types](#9-calculation-types)
10. [Precision Levels](#10-precision-levels)
11. [Auto-Convergence](#11-auto-convergence)
12. [Output Files](#12-output-files)
13. [Examples](#13-examples)
14. [Troubleshooting](#14-troubleshooting)
15. [References](#15-references)

---

## 1. Overview

The **WIEN2k Parameter Optimizer** automatically determines all preprocessing parameters for the WIEN2k density-functional-theory (DFT) code based on a `case.struct` file. It follows the strict convergence hierarchy established by Blaha et al. and documented in the official WIEN2k user guides.

**What it does:**
- Reads your `case.struct` file
- Computes pairwise nearest-neighbor distances
- Optimizes **RMT** (muffin-tin radii) through four strict physical constraints
- Determines **RKMAX** (plane-wave cutoff) using element-specific reference tables
- Computes **GMAX** (Fourier expansion cutoff for charge density)
- Sets **LMAX/LVNS** (angular momentum cutoffs per atom)
- Generates adaptive **k-point mesh** with correct Monkhors-Pack formula
- Configures **SCF mixing scheme** (PRATT/MSR1a/MSEC1) and Fermi smearing
- Sets **core-valence separation** (Ecut and HDLO recommendations)
- Generates all WIEN2k input files (`in0`, `in1`, `in2`, `inm`, `klist`)
- Optionally **auto-converges** parameters by running WIEN2k itself

**Scientific foundation:** Based on the WIEN2k FAQ (Blaha, Schwarz, Luitz), *J. Chem. Phys. 152, 074101 (2020)*, and the official WIEN2k user's guide.

---

## 2. Quick Start

```bash
# Install (one command)
./install.sh --here
export PATH="$(pwd):$PATH"

# Run — auto-detects *.struct in current directory
opt_wien2k

# Or with a file
opt_wien2k BaTiO3.struct

# Interactive wizard
opt_wien2k -i
```

That's it. The optimizer reads your struct file, optimizes all parameters, and writes input files to `./optim_results/`.

---

## 3. Installation

### Method A: Local alias (recommended for single user)

```bash
cd /path/to/wien2k-parameter-optimizer
./install.sh --here
echo 'export PATH="/path/to/wien2k-parameter-optimizer:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Method B: User install via pip

```bash
./install.sh --user
# Binary installed to ~/.local/bin/opt_wien2k
```

### Method C: System-wide (requires root)

```bash
sudo ./install.sh --system
# Binary installed to /usr/local/bin/opt_wien2k
```

### Method D: Custom prefix

```bash
./install.sh --prefix /opt/local
# Binary installed to /opt/local/bin/opt_wien2k
```

### Dependencies

- **Python 3.8+** (standard library only — no pip packages required)
- **WIEN2k** (only required for `--auto-converge` mode)

### Uninstall

```bash
./install.sh --uninstall
```

---

## 4. Basic Usage

### Run on a single struct file

```bash
opt_wien2k BaTiO3.struct
```

### Auto-detect struct file (simplest)

```bash
# If only one *.struct file exists in the current directory:
opt_wien2k

# If multiple exist, a menu appears:
#   [1] BaTiO3.struct
#   [2] Si.struct
#   Pick one [1-2]:
```

### Quick run with specific precision

```bash
opt_wien2k BaTiO3.struct --precision high
```

### Specify output directory

```bash
opt_wien2k BaTiO3.struct --output ./my_results
```

### Quiet mode (minimal output)

```bash
opt_wien2k BaTiO3.struct --quiet
```

---

## 5. Interactive Mode

Launch the step-by-step wizard with:

```bash
opt_wien2k -i
```

The wizard guides you through **8 steps**:

| Step | Description | Options |
|------|------------|---------|
| 1 | Struct file path | Auto-detected, type path, or `back` to cancel |
| 2 | Calculation type | SCF / Forces / Relaxation / Optimization / EOS / EFG |
| 3 | Precision level | Screening / Coarse / Medium / High / Very High |
| 4 | k-mesh refinement | Coarse / Medium / Fine / Very Fine |
| 5 | XC functional | PBE / PBEsol / LDA / WC / SCAN / HSE |
| 6 | Magnetic | Yes / No |
| 7 | Auto-convergence | Yes / No |
| 8 | Output directory | Default: `./optim_results` |

After confirming, all steps run automatically with a colored progress tracker.

### Post-run menu

After optimization completes in interactive mode, you can:

1. **View Full Report** — display the complete scientific report
2. **View Generated Files** — list all output files with sizes
3. **Re-run with Different Settings** — go back to parameter selection
4. **Edit RMT Manually** — adjust specific muffin-tin radii
5. **Edit RKMAX** — change the plane-wave cutoff
6. **Exit** — done

---

## 6. Command-Line Options

| Flag | Values | Default | Description |
|------|--------|---------|-------------|
| `struct_file` | path | auto-detected | Path to `case.struct` file (positional, optional) |
| `-i`, `--interactive` | — | off | Launch interactive wizard |
| `--calc-type` | scf, forces, relaxation, optimization, eos, efg | scf | WIEN2k calculation type |
| `--precision` | screening, coarse, medium, high, very_high | medium | Accuracy level |
| `--refinement` | coarse, medium, fine, very_fine | medium | k-mesh density multiplier |
| `--system-type` | metal_small, semiconductor, insulator, ... | auto | Override auto-detected system type |
| `--vxc` | pbe, lda, wc, pbesol, scan, hse | pbe | Exchange-correlation functional |
| `--magnetic` | — | off | Enable spin-polarized calculation |
| `--auto-converge` | — | off | Run WIEN2k SCF to converge parameters |
| `--output` | directory path | ./optim_results | Output directory |
| `--no-input-files` | — | off | Skip generating input files |
| `--quiet` | — | off | Minimal console output |
| `--only` | rmt, rkmax, gmax, lmax, kmesh, mixing, core | all | Run only specified step(s) |

---

## 7. Single-Step Optimization

The `--only` flag runs one or more specific optimization steps, automatically resolving dependencies.

```bash
# Only optimize muffin-tin radii
opt_wien2k BaTiO3.struct --only rmt

# Only optimize plane-wave cutoff (RMT runs silently as dependency)
opt_wien2k --only rkmax

# Multiple steps
opt_wien2k --only rmt kmesh

# Core-valence separation only (RMT as silent dependency)
opt_wien2k --only core
```

**Available steps:** `rmt` `rkmax` `gmax` `lmax` `kmesh` `mixing` `core`

**Dependency resolution:**
- `rkmax`, `gmax`, `lmax`, `core` → need `rmt` (runs silently)
- `mixing` → needs `kmesh` (runs silently)
- `rmt`, `kmesh` → no dependencies

Even in single-step mode, a full report and input files are generated using recommended defaults for non-selected steps.

---

## 8. Parameter Reference

### RMT (Muffin-Tin Radii)

**Optimized through four strict conditions** from the WIEN2k FAQ:

1. **Non-overlap:** `RMT(i) + RMT(j) ≤ 0.90 × NN_distance(i,j)`
2. **Core leakage:** `:NEC01 < 0.002` electrons (verified via SCF in auto-converge mode)
3. **Ratio balance:** `max(RMT) / min(RMT) ≤ 1.5` (~1.3 for sp-elements)
4. **Structural margin:** RMT reduced for cells with free internal coordinates

### RKMAX (Plane-Wave Cutoff)

`RKMAX = Min(RMT) × K_max`

Reference values from Blaha's table:

| Element Type | RKMAX | Notes |
|-------------|-------|-------|
| H, He | 3.5 | Very small RMT |
| Li, Be, B | 4.5 | Small atoms |
| C, N, O, F | 6.5 | Electronegative, high electron density between spheres |
| Si, P, S | 5.0 | sp-elements |
| 3d TM (Ti-Cu) | 7.0–8.0 | d-elements |
| 4d TM (Zr-Ag) | 8.0 | Heavy d-elements |
| 5d TM (Hf-Au) | 8.5 | Very heavy d-elements |
| Lanthanides | 8.0 | f-elements |
| Actinides | 8.5 | 5f-elements |

Effective RKMAX per atom: `RKMAX_eff(i) = RKMAX × Min(RMT) / RMT(i)`

### GMAX (Fourier Expansion)

| Precision | GMAX | Description |
|-----------|------|-------------|
| screening | 10.0 | Quick, rough estimates |
| coarse | 12.0 | Exploratory calculations |
| medium | 14.0 | Production quality |
| high | 16.0 | Publication quality |
| very_high | 18.0 | Benchmark quality |

GMAX increases for systems containing H, Li, halogens, or f-elements.

### LMAX/LVNS (Angular Momentum)

| Element Type | LMAX | Description |
|-------------|------|-------------|
| sp (H, O, Si, ...) | 6 | Standard sp-elements |
| d (Ti, Fe, Co, ...) | 10 | Transition metals |
| f (Ce, Eu, U, ...) | 12 | f-block elements |

**LVNS** (L-max for non-spherical contributions):
- 4 for sp-only systems
- 6 for systems with d-elements or large RMTs (> 2.2)
- 8 for systems with f-elements

### k-point Mesh

Adaptive Monkhors-Pack mesh based on:

```
k_density = target points per bohr⁻³ (depends on system type)
n_i = round(|b_i| × (k_density × N)^(1/3))
```

| System Type | k-Density | Description |
|-------------|-----------|-------------|
| metal_small | 3000 | Small metallic unit cells |
| semiconductor | 500 | Semiconductors and insulators |
| insulator | 200 | Large-gap insulators |
| metal_large | 300 | Large metallic supercells |
| insulator_large | 10 | Very large insulating supercells |

Magnetic systems get ×1.5 density multiplier.

### Mixing Scheme

| System Type | Scheme | Mixing Factor | Notes |
|-------------|--------|---------------|-------|
| Insulator | PRATT | 0.25 | Fast convergence, no smearing |
| Semiconductor | PRATT | 0.25 | Same as insulator |
| Metal (small) | MSR1a | 0.20 | 0.001 Ry Fermi smearing |
| Metal (large) | MSEC1 | 0.15 | 0.002 Ry Fermi smearing |

### TEMP (Fermi Smearing)

| System | TEMP (Ry) |
|--------|-----------|
| Metals (small cell) | 0.001 |
| Metals (large cell) | 0.002 |
| Semiconductors/insulators | 0.001 |

### Core-Valence (Ecut)

| Precision | Ecut (Ry) | HDLO |
|-----------|-----------|------|
| screening | -4.0 | No |
| coarse | -5.0 | No |
| medium | -6.0 | For RMT > 2.5 |
| high | -7.0 | For RMT > 2.2 |
| very_high | -8.0 | For RMT > 2.0 |

**HDLO** (high-derivative local orbitals) recommended for atoms with large RMT to improve the linearization.

---

## 9. Calculation Types

| Type | Description | Key Settings |
|------|------------|-------------|
| `scf` | Self-consistent field (default) | TOT mode in in2 |
| `forces` | Forces and geometry relaxation | FOR mode in in2, TETRA=NO |
| `relaxation` | Volume + positions / positions only | — |
| `optimization` | Internal positions only | — |
| `eos` | Equation of state (E vs V) | — |
| `efg` | Electric field gradient | — |

---

## 10. Precision Levels

| Level | Ecut | GMAX | ECONV | Typical Use |
|-------|------|------|-------|-------------|
| screening | -4.0 | 10.0 | 1e-5 | Quick exploratory scans |
| coarse | -5.0 | 12.0 | 1e-5 | Rough estimates |
| medium | -6.0 | 14.0 | 5e-6 | Production quality |
| high | -7.0 | 16.0 | 1e-5 | Publication quality |
| very_high | -8.0 | 18.0 | 5e-6 | Benchmark quality |

---

## 11. Auto-Convergence

When `--auto-converge` is enabled and WIEN2k is in PATH, the tool runs an automatic convergence cycle following Blaha's hierarchy:

### RMT Core Leakage Verification

```
1. Run SCF with current RMTs
2. Parse :NEC01 core leakage per atom
3. If leakage > 0.002 e⁻ for any atom → increase RMT by 5%
4. Recheck non-overlap constraint
5. Repeat until converged or max 5 iterations
```

### k-mesh Convergence

```
1. Start with 1/2 of recommended mesh
2. Compare ΔE between successive meshes
3. Refine: ×1.0 → ×1.5 → ×2.0 → ×3.0 → ×4.0
4. Stop when |ΔE| < 0.1 mRy
```

### RKMAX Convergence

```
1. Start at base_rkmax - 1.0
2. Increment by 0.5 until max_rkmax
3. Compare ΔE between successive RKMAX values
4. Stop when |ΔE| < 0.1 mRy
```

**Important:** Total energy converges faster than forces.
- Converge the **property of interest**, not just total energy.

---

## 12. Output Files

After optimization, the output directory contains:

| File | Description |
|------|-------------|
| `*_optimization_report.txt` | Full scientific report (10 sections) |
| `*.struct_optimized` | Optimized struct file with new RMTs |
| `*.in0` | Core density / exchange-correlation input |
| `*.in1` | Linearization energies per atom |
| `*.in2` | SCF convergence parameters |
| `*.inm` | Density initialization |
| `*.klist` | k-point list |

### Ready for WIEN2k

The tool prints the exact commands:

```bash
init_lapw -b -rkmax 7.0 -numk 777 -ecut 7
run_lapw -p
```

---

## 13. Examples

### Example 1: Simple semiconductor (Si)

```bash
# Battery material — silicon
opt_wien2k Si.struct --precision high

# Output (summary):
#   RMT: Si=2.000
#   RKMAX: 5.5  |  GMAX: 16.0
#   k-mesh: 5×5×5 (semiconductor, 125 pts)
#   Ecut: -7.0 Ry  |  Mixing: PRATT 0.25
```

### Example 2: Perovskite oxide (BaTiO₃)

```bash
opt_wien2k BaTiO3.struct --precision very_high --refinement fine

# Output:
#   RMT: Ba=2.076, Ti=1.927, O=1.459
#   RKMAX: 7.0  |  GMAX: 18.0
#   k-mesh: 14×14×14 (2744 pts)  |  k-density: 2744 pts
#   Ecut: -8.0 Ry  |  Mixing: PRATT 0.25
```

### Example 3: Magnetic metal (Fe)

```bash
opt_wien2k Fe.struct --magnetic --precision high

# Output:
#   RMT: Fe=2.100
#   RKMAX: 8.0  |  GMAX: 16.0
#   k-mesh: 12×12×12 (1728 pts, metal_small)  |  k-density × 1.5 (magnetic)
#   TEMP: 0.001 Ry  |  Mixing: MSR1a 0.20
```

### Example 4: Transition metal nitride (TiN)

```bash
opt_wien2k TiN.struct --precision high

# Auto-detected as metal_small (TM nitride without oxide)
#   RKMAX: 7.0  |  k-mesh: higher density for metal
#   TEMP: 0.001 Ry  |  Mixing: MSR1a 0.20
```

### Example 5: Interactive wizard with forces

```bash
opt_wien2k -i

# Step-by-step wizard:
#   1. Select Si.struct
#   2. Calculation → Forces
#   3. Precision → High
#   4. k-mesh → Fine
#   5. XC → PBE
#   6. Magnetic → No
#   7. Auto-converge → Yes
#   8. Output → ./si_forces
#   → Confirmation → Run
```

### Example 6: Only optimize RMT and k-mesh

```bash
# Quick RMT + k-mesh check without changing other params
opt_wien2k BaTiO3.struct --only rmt kmesh --no-input-files
```

### Example 7: Auto-convergence

```bash
# Full convergence with WIEN2k
opt_wien2k BaTiO3.struct --precision high --auto-converge --output ./converged

# Output includes:
#   - RMT core leakage iteration history
#   - k-mesh convergence table (ΔE values)
#   - RKMAX convergence table (ΔE values)
#   - Final converged parameters
```

### Example 8: EOS calculation

```bash
opt_wien2k Si.struct --calc-type eos --precision high

# Generates input files suitable for equation-of-state calculations
```

### Example 9: Quick screening

```bash
# Fast screening with lower precision
opt_wien2k --precision screening --no-input-files

# Saves report only — useful for quick parameter estimates
```

### Example 10: Multiple struct files in directory

```bash
ls *.struct
# BaTiO3.struct  Fe.struct  Si.struct

opt_wien2k
# ━━ MULTIPLE STRUCT FILES FOUND ━━
#   [1] BaTiO3.struct
#   [2] Fe.struct
#   [3] Si.struct
#   Pick one [1-3]: 1
```

---

## 14. Troubleshooting

### "No struct file found in current directory"

The tool couldn't find any `*.struct` file in your current directory.

**Solution:** Either `cd` to the directory containing your struct file, or specify it explicitly:
```bash
opt_wien2k /path/to/my/structure.struct
```

### "WIEN2k not found in PATH" (during --auto-converge)

The tool needs WIEN2k executables (`run_lapw`, `init_lapw`) for auto-convergence.

**Solution:** Source the WIEN2k environment first:
```bash
source /path/to/wien2k/wien2k_environment.sh
opt_wien2k BaTiO3.struct --auto-converge
```

### "RMT ratio exceeds 1.5 — CRITICAL WARNING"

The RMT optimization couldn't satisfy all four constraints simultaneously.

**Solution:**
1. Check your struct file — are the positions correct?
2. Try `--only rmt` to inspect the RMT values
3. Use the interactive post-run menu to manually adjust RMTs
4. For very tight structures, consider reducing RMT manually

### "Core leakage not converged after 5 iterations"

Some atoms have persistent core leakage > 0.002 e⁻.

**Solution:**
1. Consider transferring semi-core states to valence (more negative Ecut)
2. Use `--precision very_high` for stricter Ecut
3. Manually increase the problematic atom's RMT via interactive mode

### "SCF failed" during auto-convergence

**Solution:**
1. Check WIEN2k installation with `run_lapw --help`
2. Try reducing precision or k-mesh density first
3. Increase SCF timeout in convergence.py if needed

---

## 15. References

1. **P. Blaha, K. Schwarz, F. Tran, R. Laskowski, G.K.H. Madsen, L.D. Marks**  
   *WIEN2k: An APW+lo Program for Calculating the Properties of Solids*  
   J. Chem. Phys. **152**, 074101 (2020)

2. **P. Blaha, K. Schwarz, G.K.H. Madsen, D. Kvasnicka, J. Luitz**  
   *WIEN2k User's Guide, Version 19.1*  
   http://susi.theochem.tuwien.ac.at/

3. **WIEN2k FAQ — RMT**  
   http://www.wien2k.at/reg_user/faq/rmt.html

4. **WIEN2k FAQ — RKMAX**  
   http://www.wien2k.at/reg_user/faq/rkmax.html

5. **WIEN2k FAQ — k-point Generation**  
   http://www.wien2k.at/reg_user/faq/kgen.html

6. **L.D. Marks**  
   *Optimization Notes*, WIEN2k Workshop 2006

7. **Blaha Lectures**, WIEN2k Workshop 2015

8. **P. Blaha, K. Schwarz, P. Sorantin, S.B. Trickey**  
   *Full-potential, linearized augmented plane wave programs for crystalline systems*  
   Comput. Phys. Commun. **59**, 399 (1990)

9. **D.J. Singh, L. Nordstrom**  
   *Planewaves, Pseudopotentials, and the LAPW Method*, 2nd Ed., Springer (2006)

---

## Citation

If you use this tool in your research, please cite:

> M. Haddad, *WIEN2k Parameter Optimizer*, https://github.com/majidhaddad02/wien2k-parameter-optimizer (2024)

---

*Last updated: July 2026*
