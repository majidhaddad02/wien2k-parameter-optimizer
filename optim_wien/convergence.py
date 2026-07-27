"""
Automatic Convergence Module — Hierarchical Parameter Convergence.

Implements the Blaha convergence hierarchy:
  k-mesh → RKMAX → (forces) → (properties)

Requirements:
  - WIEN2k installed and in PATH
  - Working WIEN2k environment (w2web init'd)
  - The case directory must exist with valid struct file

References:
  - P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)
  - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
"""

import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field


_WIEN2K_AVAILABLE = None


def wien2k_available():
    """Check if WIEN2k commands are accessible."""
    global _WIEN2K_AVAILABLE
    if _WIEN2K_AVAILABLE is not None:
        return _WIEN2K_AVAILABLE
    try:
        result = subprocess.run(
            ["which", "run_lapw"], capture_output=True, text=True, timeout=5
        )
        _WIEN2K_AVAILABLE = result.returncode == 0
    except Exception:
        _WIEN2K_AVAILABLE = False
    return _WIEN2K_AVAILABLE


@dataclass
class ConvergenceResult:
    converged: bool = False
    final_rmts: list = field(default_factory=list)
    final_kmesh: tuple = (0, 0, 0)
    final_rkmax: float = 0.0
    rmt_history: list = field(default_factory=list)
    kmesh_history: list = field(default_factory=list)
    rkmax_history: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    total_runtime: float = 0.0


def auto_converge(
    case_dir,
    basename,
    structure,
    initial_rmt_result,
    initial_rkmax_result,
    initial_kmesh_result,
    initial_mixing_result,
    initial_core_valence_result,
    initial_gmax_result,
    initial_lmax_result,
    vxc_type=13,
    magnetic=False,
    kmesh_threshold=0.0001,
    rkmax_threshold=0.0001,
    max_kmesh_scale=4.0,
    max_rkmax=12.0,
    parallel=True,
):
    """
    Full automatic convergence of k-mesh and RKMAX.

    Hierarchy (Blaha):
      1. Fix RMTs and write optimized struct
      2. Initialize WIEN2k with recommended base parameters
      3. k-mesh convergence loop
      4. RKMAX convergence loop

    Returns ConvergenceResult.
    """
    if not wien2k_available():
        return ConvergenceResult(
            converged=False,
            warnings=["WIEN2k not found in PATH. Cannot auto-converge. "
                      "Using recommended parameters only."]
        )

    result = ConvergenceResult()
    t0 = time.time()
    original_dir = os.getcwd()

    try:
        os.chdir(case_dir)

        _write_all_inputs(
            basename, structure,
            initial_rmt_result, initial_rkmax_result,
            initial_kmesh_result, initial_mixing_result,
            initial_core_valence_result, initial_gmax_result,
            initial_lmax_result, vxc_type, magnetic,
        )

        _run_init(basename, initial_rkmax_result.rkmax,
                   initial_kmesh_result.mesh,
                   initial_core_valence_result.ecut,
                   parallel)

        rmt_converged, final_rmts = _converge_rmt_core_leakage(
            basename, structure,
            initial_rmt_result, initial_rkmax_result,
            initial_kmesh_result, initial_mixing_result,
            initial_core_valence_result, initial_gmax_result,
            initial_lmax_result, vxc_type, magnetic,
            parallel, result,
        )
        result.final_rmts = final_rmts

        k_converged, final_kmesh = _converge_kmesh(
            basename, initial_kmesh_result,
            initial_rkmax_result, initial_mixing_result,
            initial_core_valence_result, initial_gmax_result,
            initial_lmax_result, structure,
            vxc_type, magnetic, kmesh_threshold, max_kmesh_scale,
            parallel, result,
        )
        result.final_kmesh = final_kmesh

        r_converged, final_rkmax = _converge_rkmax(
            basename, initial_rkmax_result,
            final_kmesh, initial_mixing_result,
            initial_core_valence_result, initial_gmax_result,
            initial_lmax_result, structure,
            vxc_type, magnetic, rkmax_threshold, max_rkmax,
            parallel, result,
        )
        result.final_rkmax = final_rkmax

        result.converged = rmt_converged and k_converged and r_converged

    except Exception as e:
        result.warnings.append(f"Convergence failed: {e}")
        result.converged = False
    finally:
        os.chdir(original_dir)

    result.total_runtime = time.time() - t0
    return result


def _write_all_inputs(basename, structure, rmt, rk, km, mix, cv, gm, lm,
                      vxc, mag):
    from .input_generator import generate_all_inputs
    generate_all_inputs(
        ".", basename, structure,
        rmt, rk, gm, lm, km, mix, cv,
        "scf", vxc, mag, mag,
    )
    from .report import write_optimized_struct
    write_optimized_struct(structure, rmt.rmt_values, f"{basename}.struct")


def _run_init(basename, rkmax, kmesh, ecut, parallel):
    n1, n2, n3 = kmesh
    numk = n1 * 100 + n2 * 10 + n3
    ecut_abs = int(abs(ecut))
    pflag = "-p" if parallel else ""
    cmd = (
        f"init_lapw -b -rkmax {rkmax} -numk {numk} "
        f"-ecut {ecut_abs} {pflag}"
    )
    subprocess.run(cmd, shell=True, check=True, timeout=300)


def _run_scf(basename, parallel):
    pflag = "-p" if parallel else ""
    subprocess.run(
        f"run_lapw {pflag} -ec 0.0001 -cc 0.001",
        shell=True, check=True, timeout=3600,
    )


def _read_energy(case_file):
    try:
        with open(case_file, "r") as f:
            for line in f:
                m = re.search(r":ENE\s*:.*=\s*(-?[\d.]+)", line)
                if m:
                    return float(m.group(1))
    except (FileNotFoundError, ValueError):
        pass
    return None


def _read_core_leakage(case_file):
    """
    Parse case.scf (or case.scf2 after lstart) for :NEC01 lines.

    Returns dict: {atom_index: leakage_in_electrons}
    """
    leakage = {}
    current_atom = -1
    try:
        with open(case_file, "r") as f:
            for line in f:
                atm_match = re.match(r"\s*ATOM\s+(\d+)", line)
                if atm_match:
                    current_atom = int(atm_match.group(1)) - 1
                if ":NEC01" in line or "CORE-CHARGE OUTSIDE" in line:
                    m = re.search(
                        r"(\d+\.\d+)\s+(?:CORE|core)\s+(?:electrons|ELECTRONS)",
                        line
                    )
                    if not m:
                        m = re.search(r"OUTSIDE SPHERE:\s*(\d+\.\d+)", line)
                    if m and current_atom >= 0:
                        leakage[current_atom] = float(m.group(1))
    except (FileNotFoundError, ValueError):
        pass
    return leakage


def _read_all_core_leakage(scf_file):
    """More thorough scan for :NEC01 in the full scf output."""
    leakage = {}
    try:
        with open(scf_file, "r") as f:
            content = f.read()
        blocks = re.split(r":NEC01", content)
        for block in blocks[1:]:
            m = re.search(
                r"CORE-CHARGE OUTSIDE SPHERE:\s*(\d+\.\d+)", block
            )
            if not m:
                m = re.search(r"(\d+\.\d+)\s+CORE.*leak", block, re.IGNORECASE)
            if m:
                atm_match = re.search(r"ATOM\s+(\d+)", block)
                if atm_match:
                    idx = int(atm_match.group(1)) - 1
                    leakage[idx] = float(m.group(1))
    except (FileNotFoundError, ValueError):
        pass
    return leakage


def _update_struct_rmt(basename, structure, rmt_values):
    """Rewrite case.struct with updated RMT values."""
    a = structure.a
    b = structure.b
    c = structure.c
    alpha = structure.alpha
    beta = structure.beta
    gamma = structure.gamma

    lines = [structure.title, structure.lattice_type, structure.mode]
    lines.append(
        f"             {a:12.6f} {b:12.6f} {c:12.6f} "
        f"{alpha:6.1f} {beta:6.1f} {gamma:6.1f}"
    )
    for i, atom in enumerate(structure.atoms):
        rm = rmt_values[i] if i < len(rmt_values) else atom.rmt
        lines.append(
            f"ATOM {i+1:3d}: X={atom.position[0]:.10f} "
            f"Y={atom.position[1]:.10f} Z={atom.position[2]:.10f}"
        )
        lines.append(f"          MULT= {atom.mult:2d}          ISPLIT= {atom.isplit:2d}")
        lines.append(
            f"{atom.element:10}        NPT=  781  "
            f"R0=0.00010000 RMT=    {rm:10.5f}   Z: {atom.z:.1f}"
        )
        lines.append("LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000")
        lines.append("                     0.0000000 1.0000000 0.0000000")
        lines.append("                     0.0000000 0.0000000 1.0000000")

    with open(f"{basename}.struct", "w") as f:
        f.write("\n".join(lines) + "\n")


def _converge_rmt_core_leakage(
    basename, structure,
    initial_rmt_result, initial_rkmax_result,
    initial_kmesh_result, initial_mixing_result,
    initial_core_valence_result, initial_gmax_result,
    initial_lmax_result, vxc_type, magnetic,
    parallel, result, max_iterations=5,
):
    """
    RMT convergence via core leakage (:NEC01) check.

    Algorithm:
      1. Run SCF with current RMTs
      2. Parse :NEC01 core leakage per atom
      3. If leakage > 0.002 for any atom, increase RMT by 5%
      4. Recheck non-overlap constraint
      5. If constraint violated, reduce the other atom
      6. Rewrite struct, save density, restart
      7. Repeat until converged or max iterations
    """
    from .struct_parser import compute_pairwise_min_distances
    from .rmt import RMTOptimizationResult

    atoms = structure.atoms
    rmt_values = list(initial_rmt_result.rmt_values)
    pw = compute_pairwise_min_distances(structure)
    core_threshold = 0.002

    result.rmt_history = []

    for iteration in range(max_iterations):
        entry = {"iteration": iteration, "rmts": list(rmt_values), "leakage": {}}

        try:
            _run_scf(basename, parallel)
        except subprocess.CalledProcessError:
            result.warnings.append(
                f"RMT convergence: SCF failed at iteration {iteration}"
            )
            break

        leakage = _read_all_core_leakage(f"{basename}.scf")
        entry["leakage"] = dict(leakage)

        max_leak = max(leakage.values()) if leakage else 0.0
        entry["max_leak"] = max_leak
        result.rmt_history.append(entry)

        if max_leak < core_threshold:
            result.warnings.append(
                f"RMT core leakage converged at iteration {iteration}: "
                f"max leakage = {max_leak:.5f} e⁻ < {core_threshold}"
            )
            return True, rmt_values

        worst_atom = max(leakage, key=leakage.get) if leakage else -1
        if worst_atom < 0:
            result.warnings.append("RMT: no leakage data found. RMT not verified.")
            return True, rmt_values

        old_rmt = rmt_values[worst_atom]
        rmt_values[worst_atom] *= 1.05
        entry["action"] = (
            f"Increased {atoms[worst_atom].element}({worst_atom}) "
            f"RMT: {old_rmt:.4f} → {rmt_values[worst_atom]:.4f}"
        )

        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                nn = pw.get((i, j), 10.0)
                if rmt_values[i] + rmt_values[j] > nn:
                    if i == worst_atom:
                        rmt_values[j] = nn - rmt_values[i] - 0.02
                    else:
                        rmt_values[i] = nn - rmt_values[j] - 0.02
                    entry.setdefault("overlap_fixes", []).append(
                        f"Fixed overlap {atoms[i].element}-{atoms[j].element}: "
                        f"RMTs → {rmt_values[i]:.3f}, {rmt_values[j]:.3f}"
                    )

        _update_struct_rmt(basename, structure, rmt_values)

        _run_init(basename, initial_rkmax_result.rkmax,
                   initial_kmesh_result.mesh,
                   initial_core_valence_result.ecut,
                   parallel)

    result.warnings.append(
        f"RMT core leakage NOT converged after {max_iterations} iterations. "
        f"Max leakage = {max_leak:.5f} e⁻ > {core_threshold}. "
        f"Consider manual review or transferring states to valence."
    )
    return False, rmt_values


def _update_klist(basename, mesh, shift=(0.5, 0.5, 0.5)):
    n1, n2, n3 = mesh
    s1, s2, s3 = shift
    with open(f"{basename}.klist", "w") as f:
        f.write(f" {basename}\n")
        f.write(f" {n1:4d} {n2:4d} {n3:4d}    {n1*n2*n3}  number of k-points\n")
        f.write("   -6  1     add INV\n")
        f.write(f" {s1:5.2f} {s2:5.2f} {s3:5.2f}     shift\n")
        f.write("END\n")


def _update_in1(basename, rkmax, use_hdlo=False):
    vnmt = 10 if use_hdlo else 6
    with open(f"{basename}.in1", "r") as f:
        content = f.read()
    content = re.sub(
        r"^\s*[\d.]+\s+\d+\s+.*R-MT",
        f" {rkmax:5.1f}      {vnmt}     (R-MT*K-MAX; MAX L IN WF; V-NMT)",
        content, flags=re.MULTILINE,
    )
    with open(f"{basename}.in1", "w") as f:
        f.write(content)


def _converge_kmesh(basename, initial_kmesh, initial_rkmax, mix, cv, gm,
                     lm, structure, vxc, mag, threshold, max_scale,
                     parallel, result):
    base = initial_kmesh.mesh
    test_meshes = _generate_kmesh_sequence(base, max_scale)

    result.kmesh_history = []
    prev_energy = None
    converged_mesh = base

    for mesh in test_meshes:
        _update_klist(basename, mesh)

        try:
            _run_scf(basename, parallel)
        except subprocess.CalledProcessError:
            result.warnings.append(
                f"SCF failed for k-mesh {mesh[0]}×{mesh[1]}×{mesh[2]}"
            )
            break

        energy = _read_energy(f"{basename}.scf")
        entry = {"mesh": mesh, "energy": energy}
        result.kmesh_history.append(entry)

        if energy is None:
            result.warnings.append(
                f"Could not read energy for k-mesh {mesh[0]}×{mesh[1]}×{mesh[2]}"
            )
            break

        if prev_energy is not None:
            de = abs(energy - prev_energy)
            if de < threshold:
                converged_mesh = mesh
                result.warnings.append(
                    f"k-mesh converged at {mesh[0]}×{mesh[1]}×{mesh[2]} "
                    f"(ΔE = {de:.6f} Ry < {threshold:.6f})"
                )
                break

        prev_energy = energy
    else:
        result.warnings.append(
            f"k-mesh NOT converged up to {test_meshes[-1][0]}×"
            f"{test_meshes[-1][1]}×{test_meshes[-1][2]}. "
            f"Using finest mesh tested."
        )
        converged_mesh = test_meshes[-1]

    return True, converged_mesh


def _converge_rkmax(basename, initial_rkmax, kmesh, mix, cv, gm, lm,
                     structure, vxc, mag, threshold, max_rkmax,
                     parallel, result):
    base = initial_rkmax.rkmax
    test_rkmax_values = _generate_rkmax_sequence(base, max_rkmax)

    _update_klist(basename, kmesh)

    result.rkmax_history = []
    prev_energy = None
    converged_rkmax = base

    for rk in test_rkmax_values:
        _update_in1(basename, rk)

        try:
            _run_scf(basename, parallel)
        except subprocess.CalledProcessError:
            result.warnings.append(f"SCF failed for RKMAX = {rk}")
            break

        energy = _read_energy(f"{basename}.scf")
        entry = {"rkmax": rk, "energy": energy}
        result.rkmax_history.append(entry)

        if energy is None:
            result.warnings.append(f"Could not read energy for RKMAX = {rk}")
            break

        if prev_energy is not None:
            de = abs(energy - prev_energy)
            if de < threshold:
                converged_rkmax = rk
                result.warnings.append(
                    f"RKMAX converged at {rk} (ΔE = {de:.6f} Ry < {threshold:.6f})"
                )
                break

        prev_energy = energy
    else:
        result.warnings.append(
            f"RKMAX NOT converged up to {test_rkmax_values[-1]}. "
            f"Using finest value tested."
        )
        converged_rkmax = test_rkmax_values[-1]

    return True, converged_rkmax


def _generate_kmesh_sequence(base_mesh, max_scale):
    n1, n2, n3 = base_mesh
    seq = []
    current = (max(1, n1 // 2), max(1, n2 // 2), max(1, n3 // 2))
    seq.append(current)

    steps = [base_mesh]
    for factor in [1.5, 2.0, 3.0, 4.0]:
        candidate = (
            max(1, round(n1 * factor)),
            max(1, round(n2 * factor)),
            max(1, round(n3 * factor)),
        )
        if candidate != steps[-1]:
            steps.append(candidate)

    seq.extend(steps)

    target_max = max(1, round(max(n1, n2, n3) * max_scale))
    if seq[-1][0] > target_max or seq[-1][1] > target_max or seq[-1][2] > target_max:
        pass

    return seq


def _generate_rkmax_sequence(base_rkmax, max_rkmax):
    seq = []
    start = max(2.0, base_rkmax - 1.0)
    current = start
    while current <= max_rkmax:
        seq.append(round(current, 1))
        current += 0.5
    if not seq:
        seq = [base_rkmax]
    return seq
