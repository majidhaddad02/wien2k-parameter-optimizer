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
8. [Convergence-Verified Optimization](#8-convergence-verified-optimization)
9. [Parameter Reference](#9-parameter-reference)
10. [Calculation Types](#10-calculation-types)
11. [Precision Levels](#11-precision-levels)
12. [Auto-Convergence (Legacy)](#12-auto-convergence-legacy)
13. [Output Files](#13-output-files)
14. [Examples](#14-examples)
15. [Troubleshooting](#15-troubleshooting)
16. [References](#16-references)

---

## 1. Overview

The **WIEN2k Parameter Optimizer** automatically determines all preprocessing parameters for the WIEN2k density-functional-theory (DFT) code from a `case.struct` file. It follows the strict convergence hierarchy established by Blaha et al. and documented in the official WIEN2k user guides.

**What it does:**

- Reads your `case.struct` file (including MULT>1 equivalent positions)
- Computes pairwise nearest-neighbour distances using all symmetry-equivalent sites
- Optimizes **RMT** (muffin-tin radii) through four strict physical constraints
- Determines **RKMAX** (plane-wave cutoff) using element-specific tables as seed
- Computes **GMAX** (Fourier expansion cutoff for charge density)
- Sets **LMAX/LVNS** (angular momentum cutoffs per atom)
- Generates adaptive **k-point mesh** with correct Monkhorst-Pack formula
- Configures **SCF mixing scheme** (PRATT/MSR1a/MSEC1) and Fermi smearing
- Sets **core-valence separation** (Ecut and HDLO recommendations)
- Generates all WIEN2k input files (`in0`, `in1`, `in1c`, `in2`, `inm`, `klist`)
- **NEW:** Convergence-verified optimization — Aitken Δ² extrapolation with SCF confirmation

**Two optimization approaches:**

| Mode | Flag | Description |
|------|------|-------------|
| Table-based (fast) | default | Uses element-specific reference tables — runs instantly, no WIEN2k needed |
| Convergence-verified | `--converge rkmax,kmesh` | Runs real SCF, fits exponential convergence model, confirms via extra SCF — requires WIEN2k |
| Legacy auto-converge | `--auto-converge` | Linear sweep — slower, kept for backwards compatibility |

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

# Convergence-verified (Aitken extrapolation — fast)
opt_wien2k BaTiO3.struct --converge rkmax,kmesh --etol 1.0
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
- **WIEN2k** (only required for `--converge` or `--auto-converge` modes)

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

### Known band gap override

```bash
# Force metal treatment (useful for small-gap semiconductors)
opt_wien2k BaTiO3.struct --bandgap 0

# Force insulator treatment
opt_wien2k Si.struct --bandgap 1.5
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

### Table-based optimization flags

| Flag | Values | Default | Description |
|------|--------|---------|-------------|
| `struct_file` | path | auto-detected | Path to `case.struct` file (positional, optional) |
| `-i`, `--interactive` | — | off | Launch interactive wizard |
| `--calc-type` | scf, forces, relaxation, optimization, eos, efg | scf | WIEN2k calculation type |
| `--precision` | screening, coarse, medium, high, very_high | medium | Accuracy level |
| `--refinement` | coarse, medium, fine, very_fine | medium | k-mesh density multiplier |
| `--system-type` | metal_small, semiconductor, insulator, ... | auto | Override auto-detected system type |
| `--bandgap` | float (eV) | none | Known band gap — 0=metal, <1=semiconductor, ≥1=insulator |
| `--vxc` | pbe, lda, wc, pbesol, scan, hse | pbe | Exchange-correlation functional |
| `--magnetic` | — | off | Enable spin-polarized calculation |
| `--output` | directory path | ./optim_results | Output directory |
| `--no-input-files` | — | off | Skip generating input files |
| `--quiet` | — | off | Minimal console output |
| `--only` | rmt, rkmax, gmax, lmax, kmesh, mixing, core | all | Run only specified step(s) |

### Convergence-verified optimization flags (NEW)

| Flag | Values | Default | Description |
|------|--------|---------|-------------|
| `--converge` | rmt, rkmax, kmesh, gmax | rkmax kmesh | Parameters to converge via Aitken extrapolation |
| `--etol` | float | 1.0 | Convergence tolerance in mRy/atom |
| `--cluster-submit` | — | off | Submit convergence jobs to HPC scheduler |
| `--converge-report` | path | auto | Path for markdown convergence report |

### Legacy auto-convergence flag

| Flag | Description |
|------|-------------|
| `--auto-converge` | Linear-sweep k-mesh & RKMAX convergence (slower, kept for compatibility) |

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

## 8. Convergence-Verified Optimization

The `--converge` flag replaces blind trust in lookup tables with a minimal number of real WIEN2k SCF runs, Aitken Δ² extrapolation, and explicit confirmation steps.

### Method

**Governing principle:** every "converged" number in the final report must be traceable to actual `.scf` output — never to a table or curve fit alone.

The Aitken Δ² process (Aitken, 1926) accelerates linearly convergent sequences:
```
E_inf = E3 - (E3 - E2)² / (E3 - 2·E2 + E1)
```
for three equally-spaced parameter values (x, x+h, x+2h).

### Usage

```bash
# Converge RKMAX and k-mesh with 1 mRy/atom tolerance
opt_wien2k BaTiO3.struct --converge rkmax,kmesh --etol 1.0

# Converge all parameters, tighter tolerance for forces
opt_wien2k Si.struct --calc-type forces --converge rmt,rkmax,kmesh,gmax --etol 0.1

# Submit convergence jobs to HPC
opt_wien2k Fe.struct --magnetic --converge rkmax,kmesh --cluster-submit

# Custom report path
opt_wien2k case.struct --converge rkmax,kmesh --converge-report converge_report.md
```

### RKMAX convergence

1. **Seed:** 3 points from `RKMAX_TABLE` — x0, x0+1, x0+2 (equally spaced)
2. **SCF:** runs all 3 in parallel (loose criteria: `-ec 0.001 -cc 0.005`)
3. **Fit:** Aitken Δ² → extrapolated E_inf, decay rate α, amplitude A
4. **Predict:** solve `|A|·exp(-α·x) < etol` → predicted RKMAX (rounded UP to 0.5)
5. **Confirm:** one SCF at predicted RKMAX (production criteria: `-ec 0.0001 -cc 0.001`)
6. **Fallback:** if confirmation fails, add point + refit via exponential least-squares (up to 6 points)
7. **Report:** all raw data, fit parameters, residual, SCF-count comparison

### k-mesh convergence

**Insulators/semiconductors:** same 3-point Aitken approach as RKMAX, parameterized by effective linear mesh dimension n.

**Metals:** 4-point bracket + smearing cross-check at the finest mesh. Fermi-surface integration can converge non-monotonically — the tool checks for this and flags non-monotonic behaviour for manual review.

### RMT robustness (not a sweep)

RMT is a geometric packing constraint — there is no "RMT → ∞" limit. Instead, the tool runs a single robustness check: nudges the two hardest-constrained atoms' RMTs down ~5%, reruns SCF, and confirms ΔE is negligible.

### GMAX verification (not a sweep)

GMAX is checked by parsing WIEN2k's own `.output2` warnings after the confirmation run. If WIEN2k reports an insufficient GMAX, the tool increases it by the suggested amount and reruns.

### Efficiency

A typical RKMAX convergence uses 4 SCF runs (3 seed + 1 confirm) instead of ~8 for a linear sweep — roughly **2× fewer** with Aitken extrapolation.

### Report

A markdown report is generated containing:

- Raw (parameter, energy) tables for each swept parameter
- Fit used, its parameters (E_inf, α, A), and the confirmation residual
- SCF-count comparison vs. naive brute-force baseline
- All anomaly flags (non-monotonic, failed fit, smearing sensitivity)
- Total wall-clock time

```
# WIEN2k Convergence-Verified Optimization Report

## Final Parameters
| Parameter | Value |
|-----------|-------|
| RKMAX     | 9.5   |
| k-mesh    | 10×10×10 (1000 pts) |

## RKMAX Convergence
| RKMAX | Energy (Ry) | ΔE (Ry) |
|-------|-------------|---------|
| 6.0   | -100.9004   | —       |
| 7.0   | -100.9396   | -0.0392 |
| 8.0   | -100.9634   | -0.0238 |
| 9.5   | -100.9986   | -0.0352 |  ← confirmation

**Aitken Δ² Extrapolation:**
- E_inf = -101.0000 Ry, α = 0.500, Confirmation residual: 0.285 mRy/atom
- Efficiency: 4 SCF runs vs. ~8 linear sweep — 2.0× fewer
```

### Handling forces/phonon calculations

For `--calc-type forces` or `relaxation`, the tolerance is automatically tightened 10× (e.g., 0.1 mRy/atom instead of 1.0). This is documented in the output with a citation to Blaha 2020, Sec. III.B: "forces converge ~10× more slowly with basis size than total energy."

---

## 9. Parameter Reference

### RMT (Muffin-Tin Radii)

**Optimized through four strict conditions** from the WIEN2k FAQ:

1. **Non-overlap:** `RMT(i) + RMT(j) ≤ 0.90 × NN_distance(i,j)`
2. **Core leakage:** `:NEC01 < 0.002` electrons
3. **Ratio balance:** `max(RMT) / min(RMT) ≤ 1.5` (~1.3 for sp-elements)
4. **Structural margin:** RMT reduced for cells with free internal coordinates

Nearest-neighbour distances are computed using **all symmetry-equivalent positions** (including MULT>1 equivalent sites in real WIEN2k struct files), not just the representative position.

### RKMAX (Plane-Wave Cutoff)

`RKMAX = Min(RMT) × K_max`

Reference values from Blaha's table (used as seed in `--converge` mode):

| Element Type | RKMAX | Notes |
|-------------|-------|-------|
| H, Li, Be, B | 3.0–5.0 | Very small atoms |
| C, N, O, F | 6.5–7.0 | Electronegative, high electron density between spheres |
| Si, P, S | 5.0–5.5 | sp-elements |
| Na–Ca | 6.5 | Alkali/alkaline earth |
| 3d TM (Sc–Zn) | 7.5–8.0 | d-elements |
| 4d TM (Y–Cd) | 7.5–8.0 | Heavy d-elements |
| 5d TM (Hf–Hg) | 8.0–8.5 | Very heavy d-elements |
| Lanthanides | 8.0–8.5 | f-elements |
| Actinides | 8.5 | 5f-elements |

Effective RKMAX per atom: `RKMAX_eff(i) = RKMAX × Min(RMT) / RMT(i)`

### GMAX (Fourier Expansion)

| Precision | GMAX | Description |
|-----------|------|-------------|
| screening | 10.0 | Quick, rough estimates |
| coarse | 12.0 | Exploratory calculations |
| medium | 14.0 | Production quality |
| high | 16.0 | Publication quality |
| very_high | 20.0 | Benchmark quality |

GMAX increases for systems containing H (RMT<0.8), Li (RMT<1.2), halogens (RMT<2.0), or f-elements (RMT<2.5).

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

Adaptive Monkhorst-Pack mesh based on reciprocal-space density:

```
k_density = target points per bohr⁻³ (depends on system type)
n_i = round(|b_i| × (k_density × V_BZ)^(1/3))
```

| System Type | k-Density | Description |
|-------------|-----------|-------------|
| metal_small | 3000 | Small metallic unit cells |
| semiconductor | 500 | Semiconductors and insulators |
| insulator | 200 | Large-gap insulators |
| metal_large | 300 | Large metallic supercells |
| insulator_large | 10 | Very large insulating supercells |

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

## 10. Calculation Types

| Type | Description | Convergence Note |
|------|------------|-----------------|
| `scf` | Self-consistent field (default) | Standard `--etol 1.0` |
| `forces` | Forces and geometry relaxation | Auto-tightens etol 10×: `--etol 0.1` |
| `relaxation` | Volume + positions / positions only | Auto-tightens etol 10× |
| `optimization` | Internal positions only | Auto-tightens etol 10× |
| `eos` | Equation of state (E vs V) | — |
| `efg` | Electric field gradient | — |

---

## 11. Precision Levels

| Level | Ecut | GMAX | RKMAX Offset | Typical Use |
|-------|------|------|-------------|-------------|
| screening | -4.0 | 10.0 | -1.0 | Quick exploratory scans |
| coarse | -5.0 | 12.0 | -0.5 | Rough estimates |
| medium | -6.0 | 14.0 | 0.0 | Production quality |
| high | -7.0 | 16.0 | +0.5 | Publication quality |
| very_high | -8.0 | 20.0 | +1.5 | Benchmark quality |

---

## 12. Auto-Convergence (Legacy)

The `--auto-converge` flag runs a linear-sweep convergence cycle. It is kept for backwards compatibility but `--converge` is recommended for new work.

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

---

## 13. Output Files

After optimization, the output directory contains:

| File | Description |
|------|-------------|
| `*_optimization_report.txt` | Full scientific report (10 sections) |
| `*_convergence_report.md` | Convergence-verified markdown report (when `--converge` is used) |
| `*.struct_optimized` | Optimized struct file — valid WIEN2k format with all equivalent positions |
| `*.in0` | Core density / exchange-correlation input — valid 5-line WIEN2k format |
| `*.in1` | Linearization energies per atom — RKMAX, LMAX, V-NMT header |
| `*.in1c` | Core state linearization (WFFIL format) |
| `*.in2` | SCF convergence parameters — TETRA 101, mixing, TEMP |
| `*.inm` | Density initialization |
| `*.klist` | k-point list — full Monkhorst-Pack coordinates with weights |
| `*.lcore` | (only if core leakage detected) |

### Ready for WIEN2k

The tool prints the exact commands:

```
$ init_lapw -b -rkmax 7.0 -numk 343 -ecut 7
$ run_lapw -p
```

The `-numk` value = n1×n2×n3 (total k-points in the full BZ), correct for multi-digit meshes.

---

## 14. Examples

### Example 1: Simple semiconductor (Si)

```bash
opt_wien2k Si.struct --precision high

# Output:
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
#   RKMAX: 7.0  |  GMAX: 20.0
#   k-mesh: 14×14×14 (2744 pts)
#   Ecut: -8.0 Ry  |  Mixing: PRATT 0.25
```

### Example 3: Magnetic metal (Fe)

```bash
opt_wien2k Fe.struct --magnetic --precision high

# Output:
#   RMT: Fe=2.100
#   RKMAX: 8.0  |  GMAX: 16.0
#   k-mesh: 12×12×12 (1728 pts, metal_small)
#   TEMP: 0.001 Ry  |  Mixing: MSR1a 0.20
```

### Example 4: Convergence-verified RKMAX + k-mesh

```bash
opt_wien2k BaTiO3.struct --converge rkmax,kmesh --etol 1.0 --output ./converged

# Runs 4 RKMAX SCF + 4 k-mesh SCF + 1 confirmation = ~9 SCF total
# Generates:
#   ./converged/BaTiO3_convergence_report.md
#   ./converged/BaTiO3_optimization_report.txt
#   + all WIEN2k input files with converged parameters
```

### Example 5: Forces calculation with tight tolerance

```bash
opt_wien2k Si.struct --calc-type forces --converge rkmax,kmesh

# etol auto-tightened to 0.1 mRy/atom for forces
```

### Example 6: Only optimize RMT and k-mesh

```bash
opt_wien2k BaTiO3.struct --only rmt kmesh --no-input-files
```

### Example 7: Interactive wizard with convergence

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

### Example 8: EOS calculation

```bash
opt_wien2k Si.struct --calc-type eos --precision high
```

### Example 9: Quick screening

```bash
opt_wien2k --precision screening --no-input-files
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

## 15. Troubleshooting

### "No struct file found in current directory"

The tool couldn't find any `*.struct` file in your current directory.

**Solution:** Either `cd` to the directory containing your struct file, or specify it explicitly:
```bash
opt_wien2k /path/to/my/structure.struct
```

### "WIEN2k not found in PATH" (during --converge or --auto-converge)

The tool needs WIEN2k executables (`run_lapw`, `init_lapw`) for convergence modes.

**Solution:** Source the WIEN2k environment first:
```bash
source /path/to/wien2k/wien2k_environment.sh
opt_wien2k BaTiO3.struct --converge rkmax,kmesh
```

### "Aitken denominator ≈ 0" or "Decay ratio outside (0,1)" during convergence

The 3-point extrapolation could not fit a valid exponential convergence model.

**Solution:**
1. Check that your SCF calculations are properly converged
2. Try a wider RKMAX range manually
3. The tool will attempt a 4+ point exponential fit automatically
4. If all attempts fail, it flags for manual review rather than returning a guessed value

### "Energy NOT monotonically decreasing" during convergence

The SCF runs did not produce monotonically decreasing energies with increasing basis size. This usually signals an SCF problem.

**Solution:**
1. Check your mixing scheme — metals need MSR1a/MSEC1, not PRATT
2. Verify the core-valence separation is correct
3. Check for ghostbands or linearization errors
4. The tool aborts extrapolation and flags for manual review — this is intentional

### "RMT ratio exceeds 1.5 — CRITICAL WARNING"

The RMT optimization couldn't satisfy all four constraints simultaneously.

**Solution:**
1. Check your struct file — are the positions correct?
2. Try `--only rmt` to inspect the RMT values
3. Use the interactive post-run menu to manually adjust RMTs
4. For very tight structures, consider reducing RMT manually

### "SCF failed" during auto-convergence or converge mode

**Solution:**
1. Check WIEN2k installation with `run_lapw --help`
2. Try reducing precision or k-mesh density first
3. Verify the struct file is valid WIEN2k format
4. Check that `init_lapw` was run successfully

---

## 16. References

1. **P. Blaha, K. Schwarz, F. Tran, R. Laskowski, G.K.H. Madsen, L.D. Marks**
   *WIEN2k: An APW+lo Program for Calculating the Properties of Solids*
   J. Chem. Phys. **152**, 074101 (2020)

2. **A.C. Aitken**
   *On Bernoulli's Numerical Solution of Algebraic Equations*
   Proc. Royal Soc. Edinburgh, **46**, 289–305 (1926)

3. **P. Blaha, K. Schwarz, G.K.H. Madsen, D. Kvasnicka, J. Luitz**
   *WIEN2k User's Guide, Version 19.1*
   http://susi.theochem.tuwien.ac.at/

4. **WIEN2k FAQ — RMT**
   http://www.wien2k.at/reg_user/faq/rmt.html

5. **WIEN2k FAQ — RKMAX**
   http://www.wien2k.at/reg_user/faq/rkmax.html

6. **WIEN2k FAQ — k-point Generation**
   http://www.wien2k.at/reg_user/faq/kgen.html

7. **L.D. Marks**
   *Optimization Notes*, WIEN2k Workshop 2006

8. **Blaha Lectures**, WIEN2k Workshop 2015

9. **P. Blaha, K. Schwarz, P. Sorantin, S.B. Trickey**
   *Full-potential, linearized augmented plane wave programs for crystalline systems*
   Comput. Phys. Commun. **59**, 399 (1990)

10. **D.J. Singh, L. Nordstrom**
    *Planewaves, Pseudopotentials, and the LAPW Method*, 2nd Ed., Springer (2006)

---

## Citation

If you use this tool in your research, please cite:

> M. Haddad, *WIEN2k Parameter Optimizer*, https://github.com/majidhaddad02/wien2k-parameter-optimizer (2024)

---

*Last updated: July 2026*
