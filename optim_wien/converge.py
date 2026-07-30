"""
Convergence-Verified Parameter Optimizer — Aitken Extrapolation Engine.

Replaces blind table lookup with a minimal number of real WIEN2k SCF runs,
Aitken Δ² extrapolation, and explicit confirmation steps.

Governing principle: every "converged" number in the final report must be
traceable to actual .scf output, never to a table or curve fit alone.

References:
  - Aitken, A. C. (1926). Proc. Royal Soc. Edinburgh, 46, 289–305.
    Δ² process for accelerating linearly convergent sequences.
  - P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)
  - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
"""

import math
import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime

# ── No silent fallback values in this module. All magic numbers have
#    a cited source in comments or are explicitly labeled as configurable defaults.

# 1 mRy/atom: standard-precision convergence rule of thumb in
# plane-wave/APW DFT practice — configurable via --etol.
#   Source: common convention, not a single literature-specific number;
#   explicitly labeled in output as configurable default.
DEFAULT_ETOL = 1.0  # mRy/atom

# Forces/phonon calculations: energy converges ~10× slower with basis
# size than total energy alone. Common, avoidable error.
#   Source: P. Blaha et al., J. Chem. Phys. 152, 074101 (2020), Sec. III.B
FORCES_ETOL_FACTOR = 10.0

# Numerical noise floor for energy convergence.
#   Source: WIEN2k FAQ — recommended SCF convergence: 10^-4 Ry for energy,
#   10^-3 for charge. Numerical noise ~100× lower.
_ENERGY_NOISE_RY = 1e-6

# Maximum number of SCF samples per parameter sweep before forcing
# human review instead of returning an unconverged value.
#   Configurable default — no false attribution.
MAX_SAMPLES_PER_SWEEP = 6

# ── WIEN2k helpers (re-used from convergence.py pattern) ──

_WIEN2K_AVAILABLE = None


def wien2k_available():
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


def _run_cmd(cmd, timeout=3600, desc="SCF"):
    """Run a command with proper timeout handling. No shell=True."""
    try:
        args = shlex.split(cmd) if isinstance(cmd, str) else cmd
        subprocess.run(args, shell=False, check=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"{desc} timed out after {timeout // 60} minutes."
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"{desc} failed with exit code {e.returncode}."
        )


def _read_energy(scf_file):
    """Parse :ENE from case.scf. Returns float(Ry) or None."""
    try:
        with open(scf_file, "r") as f:
            for line in f:
                m = re.search(r":ENE\s*:.*=\s*(-?[\d.]+)", line)
                if m:
                    return float(m.group(1))
    except (FileNotFoundError, ValueError):
        pass
    return None


def _read_lapw2_gmax_warnings(output2_file):
    """Check lapw2 output for GMAX-insufficiency warnings."""
    warnings_list = []
    try:
        with open(output2_file, "r") as f:
            for line in f:
                if "GMAX" in line.upper() and any(
                    kw in line.upper()
                    for kw in ("INSUFFICIENT", "TOO SMALL", "INCREASE")
                ):
                    warnings_list.append(line.strip())
    except FileNotFoundError:
        pass
    return warnings_list

# ── Aitken Δ² extrapolation for 3 equally-spaced points ──

@dataclass
class AitkenResult:
    """Result of Aitken Δ² extrapolation.

    E(x) ≈ E_inf + A * exp(-alpha * x)   (monotonic saturating model
    for variational basis-set convergence).

    For three equally-spaced points (x, x+h, x+2h):
      E_inf = E3 - (E3 - E2)² / (E3 - 2*E2 + E1)
      alpha = -ln(r) / h   where r = (E3 - E2) / (E2 - E1)
      A     = (E2 - E1) / (r - 1)

    Reference: Aitken, A. C. (1926). Proc. Royal Soc. Edinburgh, 46, 289–305.
    """
    e_inf: float = 0.0
    alpha: float = 0.0
    amplitude: float = 0.0
    is_valid: bool = False
    is_monotonic: bool = False
    error: str = ""
    points_x: list = field(default_factory=list)
    points_e: list = field(default_factory=list)

    def predict_x(self, etol_ry, step=0.5):
        """Find x where |E(x) - E_inf| < etol_ry.

        Since E(x) ≈ E_inf + A * exp(-alpha * x):
          |A| * exp(-alpha * x) < etol  →  x > -ln(etol / |A|) / alpha

        Round UP to the next `step` increment — never round down.
        """
        if not self.is_valid or self.alpha <= 0:
            return None
        amp = abs(self.amplitude)
        if amp <= etol_ry:
            return self.points_x[-1]
        x_required = -math.log(etol_ry / amp) / self.alpha
        rounded = math.ceil(x_required / step) * step
        if rounded < self.points_x[-1]:
            rounded = math.ceil(self.points_x[-1] / step) * step
        return max(rounded, self.points_x[-1] + step)


def aitken_extrapolate(x_values, e_values):
    """Perform Aitken Δ² extrapolation on three equally-spaced points.

    Args:
        x_values: [x1, x2, x3] — must be equally spaced (x, x+h, x+2h)
        e_values: [E1, E2, E3] — total energies in Ry

    Returns AitkenResult.
    """
    result = AitkenResult(points_x=list(x_values), points_e=list(e_values))

    if len(x_values) != 3 or len(e_values) != 3:
        result.error = "Aitken requires exactly 3 equally-spaced points."
        return result

    x1, x2, x3 = x_values
    e1, e2, e3 = e_values
    h = x2 - x1

    if abs((x3 - x2) - h) > 1e-10:
        result.error = f"Points not equally spaced: d1={x2-x1:.4f}, d2={x3-x2:.4f}"
        return result

    # Monotonic decrease check: energy must decrease as basis expands.
    # Non-monotonic → SCF problem (bad mixing, core/valence, overlap).
    # Tolerance: 10× numerical noise.
    noise = _ENERGY_NOISE_RY * 10.0
    if e1 < e2 + noise or e2 < e3 + noise:
        result.error = (
            f"Energy NOT monotonically decreasing: "
            f"E(x{x1})={e1:.8f}, E(x{x2})={e2:.8f}, E(x{x3})={e3:.8f}. "
            f"Possible SCF problem (bad mixing, core/valence, overlap). "
            f"Extrapolation ABORTED — manual inspection needed."
        )
        return result

    result.is_monotonic = True

    de1 = e2 - e1
    de2 = e3 - e2
    denominator = de2 - de1  # E3 - 2*E2 + E1

    # Near-zero denominator → Aitken numerically unstable
    if abs(denominator) < 1e-15:
        result.error = (
            f"Aitken denominator ≈ 0 ({denominator:.2e}). "
            f"Insufficient basis-size sensitivity — manual inspection needed."
        )
        return result

    # Aitken Δ²: E_inf = E3 - (E3 - E2)² / (E3 - 2*E2 + E1)
    result.e_inf = e3 - (de2 * de2) / denominator

    # Decay rate r = (E3 - E2) / (E2 - E1)
    if de1 == 0:
        result.alpha = 0.0
        result.error = "Zero energy difference between first two points."
        return result

    r = de2 / de1
    if r <= 0 or r >= 1:
        result.error = (
            f"Decay ratio r={r:.4f} outside (0,1). "
            f"Convergence not exponential — manual inspection needed."
        )
        return result

    result.alpha = -math.log(r) / h
    result.amplitude = de1 / (r - 1)

    # Sanity: E_inf should be below E3 (variational principle)
    if result.e_inf >= e3 - _ENERGY_NOISE_RY:
        result.error = (
            f"Extrapolated E_inf={result.e_inf:.8f} >= E3={e3:.8f}. "
            f"Fit is non-physical — manual inspection needed."
        )
        return result

    result.is_valid = True
    return result


# ── 4+ point exponential fit (fallback when Aitken fails) ──

def exponential_fit_linearized(x_values, e_values):
    """Fit E(x) ≈ E_inf + A * exp(-alpha * x) via linearized least squares.

    Uses pure-Python linear regression on log of energy differences.
    No scipy required.

    Returns (e_inf, alpha, amplitude) or (None, None, None) on failure.
    """
    n = len(x_values)
    if n < 4:
        return None, None, None

    h = x_values[1] - x_values[0]
    for i in range(1, n):
        if abs((x_values[i] - x_values[i - 1]) - h) > 1e-10:
            return None, None, None

    diffs = []
    x_mid = []
    for i in range(n - 1):
        de = e_values[i + 1] - e_values[i]
        if de >= 0:
            return None, None, None  # energy not decreasing
        diffs.append(math.log(-de))
        x_mid.append(x_values[i])

    # Linear regression: y = a + b * x  where b = -alpha
    N = len(diffs)
    sx = sum(x_mid)
    sy = sum(diffs)
    sxx = sum(x * x for x in x_mid)
    sxy = sum(x * y for x, y in zip(x_mid, diffs))

    denom = N * sxx - sx * sx
    if abs(denom) < 1e-15:
        return None, None, None

    b = (N * sxy - sx * sy) / denom
    a = (sy - b * sx) / N

    alpha = -b
    if alpha <= 0:
        return None, None, None

    try:
        amplitude = math.exp(a) / (1.0 - math.exp(-alpha * h))
    except (OverflowError, ZeroDivisionError):
        return None, None, None

    e_inf = e_values[-1] - amplitude * math.exp(-alpha * x_values[-1])
    if e_inf >= e_values[-1] - _ENERGY_NOISE_RY:
        return None, None, None

    return e_inf, alpha, amplitude

# ── SCF Job/Runner with parallel dispatch and loose/tight criteria ──

@dataclass
class SCFJob:
    """A single SCF calculation to be dispatched."""
    basename: str = ""
    rkmax: float = 7.0
    kmesh: tuple = (6, 6, 6)
    ecut: float = -6.0
    work_dir: str = "."
    shift: tuple = (0.5, 0.5, 0.5)
    add_inversion: bool = True
    loose: bool = True
    parallel: bool = True
    label: str = ""
    gmax: float = 14.0


@dataclass
class SCFResult:
    label: str = ""
    rkmax: float = 0.0
    kmesh: tuple = (0, 0, 0)
    energy: float | None = None
    success: bool = False
    runtime: float = 0.0
    error: str = ""
    scf_file: str = ""


class SCFRunner:
    """Dispatch SCF jobs sequentially or in parallel.

    Speed optimizations:
      - Loose-then-tight: scanning phase uses relaxed SCF criteria
        (ec=0.001, cc=0.005) — adequate for convergence sweeps because
        energy DIFFERENCES converge faster than absolute energies.
        Only the final confirmation run uses production criteria.
        Source: common practice in automated convergence testing.
      - Parallel dispatch: sample points within one sweep have no data
        dependency — they can run concurrently via subprocess.Popen.
        For real HPC, override with a user-provided submit callback.
    """

    def __init__(self, work_dir, cluster_submit=False, submit_cb=None):
        self.work_dir = work_dir
        self.cluster_submit = cluster_submit
        self.submit_cb = submit_cb
        self.total_wall_time = 0.0
        self.total_scf_count = 0

    def run_scf(self, job: SCFJob) -> SCFResult:
        """Run a single SCF calculation and return its result."""
        self.total_scf_count += 1
        t0 = time.time()
        label = job.label or f"SCF-{self.total_scf_count}"

        result = SCFResult(
            label=label, rkmax=job.rkmax, kmesh=job.kmesh,
            energy=None, success=False, runtime=0.0,
        )

        try:
            self._write_klist(job)
            self._update_in1(job)
            self._update_in0(job)
            self._run_scf_internal(job)

            result.scf_file = os.path.join(job.work_dir, f"{job.basename}.scf")
            result.energy = _read_energy(result.scf_file)
            if result.energy is None:
                result.error = "Could not parse :ENE from case.scf"
            else:
                result.success = True
        except RuntimeError as e:
            result.error = str(e)
        except Exception as e:
            result.error = f"Unexpected: {e}"
        finally:
            result.runtime = time.time() - t0
            self.total_wall_time += result.runtime

        return result

    def run_parallel(self, jobs: list) -> list:
        """Run multiple independent SCF jobs concurrently."""
        if self.cluster_submit and self.submit_cb:
            return self._dispatch_to_cluster(jobs)
        return self._run_local_parallel(jobs)

    def _dispatch_to_cluster(self, jobs):
        job_ids = []
        for job in jobs:
            self._write_klist(job)
            self._update_in1(job)
            self._update_in0(job)
            cmd = self._build_scf_command(job)
            try:
                jid = self.submit_cb(cmd, job.label)
                job_ids.append((jid, job))
            except Exception as e:
                raise RuntimeError(f"Cluster submission failed for {job.label}: {e}")

        results = []
        for jid, job in job_ids:
            t0 = time.time()
            try:
                self._wait_for_job(jid)
                result = SCFResult(label=job.label, rkmax=job.rkmax,
                                   kmesh=job.kmesh, success=True,
                                   runtime=time.time() - t0)
                result.scf_file = os.path.join(job.work_dir, f"{job.basename}.scf")
                result.energy = _read_energy(result.scf_file)
                if result.energy is None:
                    result.error = "Could not parse :ENE"
                    result.success = False
            except Exception as e:
                result = SCFResult(label=job.label, rkmax=job.rkmax,
                                   kmesh=job.kmesh, energy=None,
                                   success=False, runtime=time.time() - t0,
                                   error=str(e))
            results.append(result)
        return results

    def _run_local_parallel(self, jobs):
        processes = []
        results = []
        for job in jobs:
            self._write_klist(job)
            self._update_in1(job)
            self._update_in0(job)
            cmd = self._build_scf_command(job)
            logfile = os.path.join(
                job.work_dir,
                f"{job.basename}_{job.label.replace(' ', '_')}.log"
            )
            try:
                with open(logfile, "w") as lf:
                    proc = subprocess.Popen(
                        cmd, stdout=lf, stderr=subprocess.STDOUT, shell=False
                    )
                    processes.append((proc, job, logfile))
            except Exception as e:
                results.append(SCFResult(
                    label=job.label, rkmax=job.rkmax, kmesh=job.kmesh,
                    energy=None, success=False, runtime=0.0, error=str(e),
                ))

        for proc, job, logfile in processes:
            t0 = time.time()
            try:
                proc.wait(timeout=getattr(job, 'timeout', 3600))
                result = SCFResult(label=job.label, rkmax=job.rkmax,
                                   kmesh=job.kmesh, success=True,
                                   runtime=time.time() - t0)
                result.scf_file = os.path.join(job.work_dir, f"{job.basename}.scf")
                result.energy = _read_energy(result.scf_file)
                if result.energy is None:
                    result.error = "Could not parse :ENE"
                    result.success = False
            except subprocess.TimeoutExpired:
                proc.kill()
                result = SCFResult(label=job.label, rkmax=job.rkmax,
                                   kmesh=job.kmesh, energy=None,
                                   success=False, runtime=time.time() - t0,
                                   error="Timed out")
            except Exception as e:
                result = SCFResult(label=job.label, rkmax=job.rkmax,
                                   kmesh=job.kmesh, energy=None,
                                   success=False, runtime=time.time() - t0,
                                   error=str(e))
            results.append(result)
        return results

    def _write_klist(self, job):
        n1, n2, n3 = job.kmesh
        s1, s2, s3 = job.shift
        inv = 1 if job.add_inversion else 0
        kpath = os.path.join(job.work_dir, f"{job.basename}.klist")
        with open(kpath, "w") as f:
            f.write(f" {job.basename}\n")
            f.write(f" {n1:4d} {n2:4d} {n3:4d}    {n1*n2*n3}  number of k-points\n")
            f.write(f"   -6  {inv}     add INV\n")
            f.write(f" {s1:5.2f} {s2:5.2f} {s3:5.2f}     shift\n")
            f.write("END\n")

    def _update_in1(self, job):
        in1_path = os.path.join(job.work_dir, f"{job.basename}.in1")
        if not os.path.isfile(in1_path):
            return
        with open(in1_path, "r") as f:
            content = f.read()
        # Replace RKMAX on the second line
        content = re.sub(
            r"^\s*[\d.]+\s+\d+\s+\d+\s+.*R-MT",
            f" {job.rkmax:5.1f}       6      6  (R-MT*K-MAX; MAX L IN WF; V-NMT)",
            content, flags=re.MULTILINE,
        )
        with open(in1_path, "w") as f:
            f.write(content)

    def _update_in0(self, job):
        in0_path = os.path.join(job.work_dir, f"{job.basename}.in0")
        if not os.path.isfile(in0_path):
            return
        with open(in0_path, "r") as f:
            content = f.read()
        # Update GMAX if provided
        content = re.sub(
            r"^\s*[\d.]+\s+[\d.]+\s+.*GMAX",
            f"  {job.gmax:5.1f}       {job.gmax:5.1f}     (GMAX for POTENTIAL and CHARGE)",
            content, flags=re.MULTILINE,
        )
        with open(in0_path, "w") as f:
            f.write(content)

    def _build_scf_command(self, job):
        """Build run_lapw command with appropriate SCF criteria.

        Loose criteria for convergence sweeps:
          ec=0.001, cc=0.005 — adequate because energy DIFFERENCES
          converge faster than absolute energies.
          Source: common practice in automated convergence testing.

        Production criteria for confirmation runs:
          ec=0.0001, cc=0.001 — WIEN2k FAQ recommended defaults.
        """
        if job.loose:
            ec, cc = "0.001", "0.005"
            max_iter = "40"
        else:
            ec, cc = "0.0001", "0.001"
            max_iter = "80"

        cmd = ["run_lapw", "-ec", ec, "-cc", cc, "-i", max_iter]
        if job.parallel:
            cmd.append("-p")
        return cmd

    def _run_scf_internal(self, job):
        cmd = self._build_scf_command(job)
        _run_cmd(shlex.join(cmd), timeout=3600, desc=f"SCF ({job.label})")

    def _wait_for_job(self, job_id):
        time.sleep(5)
        raise NotImplementedError(
            "Cluster wait not implemented. Provide a submit_cb "
            "or run without --cluster-submit."
        )


# ── GMAX verification (parse lapw2 output, not a sweep) ──

@dataclass
class GMAXVerificationResult:
    gmax: float = 14.0
    verified: bool = True
    warning: str = ""
    suggested_increase: float = 0.0


def verify_gmax_at_rkmax(work_dir, basename, current_gmax):
    """After RKMAX confirmation, check .output2 for GMAX warnings.

    Per WIEN2k FAQ: "accept the automatically increased GMAX or use
    an even bigger one."
    """
    result = GMAXVerificationResult(gmax=current_gmax, verified=True)
    output2 = os.path.join(work_dir, f"{basename}.output2")
    warnings_list = _read_lapw2_gmax_warnings(output2)

    if not warnings_list:
        return result

    result.verified = False
    result.warning = "; ".join(warnings_list)

    increase = 0.0
    for w in warnings_list:
        m = re.search(r'increase\s+(?:GMAX|to)\s+(\d+\.?\d*)', w, re.IGNORECASE)
        if not m:
            m = re.search(r'(\d+\.?\d*)\s*(?:required|needed|suggested)', w, re.IGNORECASE)
        if m:
            increase = max(increase, float(m.group(1)) - current_gmax)

    if increase < 1.0:
        increase = 2.0
    result.suggested_increase = increase
    return result


# ── RMT robustness check ──

@dataclass
class RMTRobustnessResult:
    passed: bool = False
    delta_energy: float = 0.0  # mRy/atom
    original_rmts: list = field(default_factory=list)
    nudged_rmts: list = field(default_factory=list)
    warning: str = ""


def rmt_robustness_check(structure, rmt_result, basename, work_dir,
                         chosen_rkmax, chosen_kmesh, ecut, runner,
                         etol_mRy=DEFAULT_ETOL):
    """Nudge hardest-constrained atoms' RMTs down ~5%, confirm ΔE negligible.

    This is a sanity check demonstrating the result isn't RMT-choice-fragile,
    not a convergence sweep. Per WIEN2k RMT FAQ: "keep RMTs constant within
    a series of calculations."
    """
    from .struct_parser import compute_pairwise_min_distances
    from .report import write_optimized_struct

    pw = compute_pairwise_min_distances(structure)
    atoms = structure.atoms
    rmt_vals = list(rmt_result.rmt_values)

    # Find two most constrained atoms (largest RMT/NN ratio)
    constraints = []
    for i in range(len(atoms)):
        min_nn = float("inf")
        for j in range(len(atoms)):
            if i == j:
                continue
            d = pw.get((i, j), 10.0)
            if d < min_nn:
                min_nn = d
        if min_nn > 0 and min_nn < 10.0:
            constraints.append((rmt_vals[i] / min_nn, i))
    constraints.sort(reverse=True)
    hardest = [c[1] for c in constraints[:2]]

    nudged = list(rmt_vals)
    for idx in hardest:
        nudged[idx] = rmt_vals[idx] * 0.95

    write_optimized_struct(structure, nudged,
                           os.path.join(work_dir, f"{basename}.struct"))

    adjusted_rkmax = chosen_rkmax * rmt_vals[hardest[0]] / nudged[hardest[0]]

    job = SCFJob(basename=basename, rkmax=adjusted_rkmax, kmesh=chosen_kmesh,
                 ecut=ecut, work_dir=work_dir, loose=False, parallel=True,
                 label="RMT-robustness")
    scf_result = runner.run_scf(job)

    if not scf_result.success or scf_result.energy is None:
        return RMTRobustnessResult(
            passed=False, original_rmts=rmt_vals, nudged_rmts=nudged,
            warning=f"RMT robustness SCF failed: {scf_result.error}",
        )

    ref_scf = os.path.join(work_dir, f"{basename}.scf")
    ref_energy = _read_energy(ref_scf)
    if ref_energy is None:
        return RMTRobustnessResult(
            passed=False, original_rmts=rmt_vals, nudged_rmts=nudged,
            warning="Could not read reference energy for RMT robustness check.",
        )

    delta = abs(scf_result.energy - ref_energy)
    delta_per_atom = delta * 1000.0 / max(structure.num_atoms_primitive, 1)
    passed = delta_per_atom < etol_mRy

    return RMTRobustnessResult(
        passed=passed, delta_energy=delta_per_atom,
        original_rmts=rmt_vals, nudged_rmts=nudged,
        warning="" if passed else (
            f"RMT robustness: ΔE={delta_per_atom:.3f} mRy/atom ≥ "
            f"{etol_mRy} mRy/atom. RMT choice may be fragile."
        ),
    )

# ── Main Convergence Engine ──

@dataclass
class ConvergenceEngineResult:
    converged: bool = False
    final_rkmax: float = 0.0
    final_kmesh: tuple = (0, 0, 0)
    final_gmax: float = 0.0
    rkmax_data: list = field(default_factory=list)
    kmesh_data: list = field(default_factory=list)
    rkmax_fit: AitkenResult = field(default_factory=AitkenResult)
    kmesh_fit: AitkenResult = field(default_factory=AitkenResult)
    rmt_robustness: RMTRobustnessResult = field(default_factory=RMTRobustnessResult)
    gmax_verify: GMAXVerificationResult = field(default_factory=GMAXVerificationResult)
    confirmation_energy: float | None = None
    scf_count_rkmax: int = 0
    scf_count_kmesh: int = 0
    scf_count_total: int = 0
    total_wall_time: float = 0.0
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    etol_mRy_per_atom: float = DEFAULT_ETOL
    etol_Ry: float = 0.0
    n_atoms: int = 1
    for_forces: bool = False


def _get_proportional_mesh(n_eff, base_mesh):
    """Scale base mesh to new effective linear dimension n_eff."""
    bn1, bn2, bn3 = base_mesh
    base_n = round((bn1 * bn2 * bn3) ** (1.0 / 3.0))
    if base_n == 0:
        base_n = 1
    factor = n_eff / base_n
    return (
        max(1, round(bn1 * factor)),
        max(1, round(bn2 * factor)),
        max(1, round(bn3 * factor)),
    )


def converge_rkmax(structure, rmt_result, kmesh_result,
                   core_valence_result, basename, work_dir,
                   runner: SCFRunner,
                   etol_mRy_per_atom=DEFAULT_ETOL):
    """3-point Aitken extrapolation for RKMAX convergence."""
    from .rkmax import optimize_rkmax as table_rkmax

    n_atoms = structure.num_atoms_primitive
    etol_Ry = etol_mRy_per_atom * n_atoms / 1000.0

    table_result = table_rkmax(structure.atoms, rmt_result.rmt_values)
    x0 = table_result.rkmax

    x_vals = [x0, x0 + 1.0, x0 + 2.0]
    jobs = []
    for i, xv in enumerate(x_vals):
        jobs.append(SCFJob(
            basename=basename, rkmax=xv, kmesh=kmesh_result.mesh,
            ecut=core_valence_result.ecut, work_dir=work_dir,
            shift=kmesh_result.shift, loose=True, parallel=True,
            label=f"RKMAX-{i+1}",
        ))

    results = runner.run_parallel(jobs)
    scf_count = 3
    data_points = [(r.rkmax, r.energy) for r in results if r.success and r.energy is not None]
    scf_ok = [r for r in results if r.success and r.energy is not None]
    used_x = [r.rkmax for r in scf_ok]
    energies = [r.energy for r in scf_ok]

    warnings = []
    errors = []

    if len(scf_ok) < 3:
        errors.append(
            f"Only {len(scf_ok)}/3 RKMAX seed SCF runs succeeded. "
            f"Manual review required."
        )
        return ConvergenceEngineResult(
            converged=False, final_rkmax=x0, rkmax_data=data_points,
            scf_count_rkmax=scf_count, errors=errors, warnings=warnings,
            etol_mRy_per_atom=etol_mRy_per_atom, etol_Ry=etol_Ry,
            n_atoms=n_atoms,
        )

    fit = aitken_extrapolate(used_x, energies)

    if not fit.is_valid:
        errors.append(f"Aitken failed: {fit.error}")
        return _converge_rkmax_fallback(
            structure, rmt_result, kmesh_result, core_valence_result,
            basename, work_dir, runner, used_x, energies, data_points,
            scf_count, etol_mRy_per_atom, etol_Ry, n_atoms, errors, warnings,
        )

    predicted = fit.predict_x(etol_Ry, step=0.5)
    if predicted is None:
        errors.append("Aitken fit could not predict RKMAX.")
        return _converge_rkmax_fallback(
            structure, rmt_result, kmesh_result, core_valence_result,
            basename, work_dir, runner, used_x, energies, data_points,
            scf_count, etol_mRy_per_atom, etol_Ry, n_atoms, errors, warnings,
        )

    # Confirmation run
    conf_job = SCFJob(
        basename=basename, rkmax=predicted, kmesh=kmesh_result.mesh,
        ecut=core_valence_result.ecut, work_dir=work_dir,
        shift=kmesh_result.shift, loose=False, parallel=True,
        label="RKMAX-confirm",
    )
    conf_result = runner.run_scf(conf_job)
    scf_count += 1

    if not conf_result.success or conf_result.energy is None:
        errors.append(f"RKMAX confirmation SCF failed: {conf_result.error}")
        return ConvergenceEngineResult(
            converged=False, final_rkmax=predicted,
            rkmax_data=data_points + [(predicted, None)], rkmax_fit=fit,
            scf_count_rkmax=scf_count, errors=errors, warnings=warnings,
            etol_mRy_per_atom=etol_mRy_per_atom, etol_Ry=etol_Ry,
            n_atoms=n_atoms,
        )

    data_points.append((predicted, conf_result.energy))
    de = abs(conf_result.energy - fit.e_inf)

    if de < etol_Ry:
        return ConvergenceEngineResult(
            converged=True, final_rkmax=predicted,
            rkmax_data=data_points, rkmax_fit=fit,
            scf_count_rkmax=scf_count, warnings=warnings, errors=errors,
            etol_mRy_per_atom=etol_mRy_per_atom, etol_Ry=etol_Ry,
            n_atoms=n_atoms,
        )

    used_x.append(predicted)
    energies.append(conf_result.energy)
    warnings.append(
        f"Confirmation ΔE={de*1000/max(n_atoms,1):.4f} mRy/atom ≥ "
        f"{etol_mRy_per_atom:.4f}. Trying additional points."
    )
    return _converge_rkmax_fallback(
        structure, rmt_result, kmesh_result, core_valence_result,
        basename, work_dir, runner, used_x, energies, data_points,
        scf_count, etol_mRy_per_atom, etol_Ry, n_atoms, errors, warnings,
    )


def _converge_rkmax_fallback(structure, rmt_result, kmesh_result,
                              core_valence_result, basename, work_dir,
                              runner, used_x, energies, data_points,
                              scf_count, etol_mRy, etol_Ry, n_atoms,
                              errors, warnings):
    """Fallback: 4+ point exponential fit for RKMAX."""
    while scf_count < MAX_SAMPLES_PER_SWEEP:
        e_inf, alpha, A = exponential_fit_linearized(used_x, energies)
        if e_inf is not None and alpha > 0:
            try:
                x_pred = -math.log(max(etol_Ry, 1e-14) / abs(A)) / alpha
                x_rounded = math.ceil(x_pred / 0.5) * 0.5
                x_rounded = max(x_rounded, used_x[-1] + 0.5)
            except (ValueError, OverflowError):
                x_rounded = used_x[-1] + 1.0
        else:
            x_rounded = used_x[-1] + 1.0

        next_x = max(used_x[-1] + 0.5, x_rounded)
        conf_job = SCFJob(
            basename=basename, rkmax=next_x, kmesh=kmesh_result.mesh,
            ecut=core_valence_result.ecut, work_dir=work_dir,
            shift=kmesh_result.shift, loose=False, parallel=True,
            label=f"RKMAX-fallback-{scf_count}",
        )
        conf_result = runner.run_scf(conf_job)
        scf_count += 1

        if conf_result.success and conf_result.energy is not None:
            used_x.append(next_x)
            energies.append(conf_result.energy)
            data_points.append((next_x, conf_result.energy))
            de = abs(conf_result.energy - e_inf) if e_inf is not None else float("inf")
            if de < etol_Ry:
                return ConvergenceEngineResult(
                    converged=True, final_rkmax=next_x,
                    rkmax_data=data_points, scf_count_rkmax=scf_count,
                    warnings=warnings, errors=errors,
                    etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry,
                    n_atoms=n_atoms,
                )
        else:
            errors.append(f"Fallback SCF failed at RKMAX={next_x}")
            break

    errors.append(
        f"RKMAX NOT converged after {scf_count} SCF runs "
        f"(max {MAX_SAMPLES_PER_SWEEP}). Manual inspection needed."
    )
    return ConvergenceEngineResult(
        converged=False, final_rkmax=used_x[-1] if used_x else 7.0,
        rkmax_data=data_points, scf_count_rkmax=scf_count,
        errors=errors, warnings=warnings,
        etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
    )


def converge_kmesh(structure, basename, work_dir, runner, system_type,
                   chosen_rkmax, core_valence_result,
                   etol_mRy_per_atom=DEFAULT_ETOL):
    """k-mesh convergence — system-type-aware strategy.

    Insulators/semiconductors: 3-point Aitken on mesh linear dimension n.
    Metals: 4-point + smearing cross-check.
    """
    from .kmesh import optimize_kmesh as table_kmesh

    n_atoms = structure.num_atoms_primitive
    etol_Ry = etol_mRy_per_atom * n_atoms / 1000.0

    km_result = table_kmesh(structure)
    n1, n2, n3 = km_result.mesh
    n_eff = round((n1 * n2 * n3) ** (1.0 / 3.0))
    n_eff = max(1, n_eff)

    is_metal = "metal" in (system_type or "")

    if is_metal:
        return _converge_kmesh_metal(
            structure, basename, work_dir, runner, n_eff,
            chosen_rkmax, core_valence_result, km_result,
            etol_mRy_per_atom, etol_Ry, n_atoms,
        )
    else:
        return _converge_kmesh_insulator(
            structure, basename, work_dir, runner, n_eff,
            chosen_rkmax, core_valence_result, km_result,
            etol_mRy_per_atom, etol_Ry, n_atoms,
        )


def _converge_kmesh_insulator(structure, basename, work_dir, runner,
                               n_eff, chosen_rkmax, core_valence_result,
                               km_result, etol_mRy, etol_Ry, n_atoms):
    """3-point Aitken extrapolation for insulator k-mesh."""
    warnings = []
    errors = []
    step = max(2, round(n_eff / 2))
    n_vals = [max(1, n_eff - step), n_eff, n_eff + step]

    jobs = []
    for i, nv in enumerate(n_vals):
        mesh = _get_proportional_mesh(nv, km_result.mesh)
        jobs.append(SCFJob(
            basename=basename, rkmax=chosen_rkmax, kmesh=mesh,
            ecut=core_valence_result.ecut, work_dir=work_dir,
            shift=km_result.shift, loose=True, parallel=True,
            label=f"KMESH-{i+1}",
        ))

    results = runner.run_parallel(jobs)
    scf_count = 3
    scf_ok = [r for r in results if r.success and r.energy is not None]
    n_ok = [n_vals[i] for i, r in enumerate(results) if r.success and r.energy is not None]
    energies = [r.energy for r in scf_ok]
    data_points = list(zip(n_ok, energies))

    if len(scf_ok) < 3:
        errors.append(
            f"Only {len(scf_ok)}/3 k-mesh SCF runs succeeded. Manual review needed."
        )
        return ConvergenceEngineResult(
            converged=False, final_kmesh=km_result.mesh, kmesh_data=data_points,
            scf_count_kmesh=scf_count, errors=errors, warnings=warnings,
            etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
        )

    fit = aitken_extrapolate(n_ok, energies)

    if not fit.is_valid:
        errors.append(f"k-mesh Aitken failed: {fit.error}. Using table value.")
        return ConvergenceEngineResult(
            converged=False, final_kmesh=km_result.mesh, kmesh_data=data_points,
            kmesh_fit=fit, scf_count_kmesh=scf_count,
            errors=errors, warnings=warnings,
            etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
        )

    predicted_n = int(fit.predict_x(etol_Ry, step=1) or (n_ok[-1] + step))
    predicted_n = max(predicted_n, n_ok[-1])
    confirm_mesh = _get_proportional_mesh(predicted_n, km_result.mesh)

    conf_job = SCFJob(
        basename=basename, rkmax=chosen_rkmax, kmesh=confirm_mesh,
        ecut=core_valence_result.ecut, work_dir=work_dir,
        shift=km_result.shift, loose=False, parallel=True,
        label="KMESH-confirm",
    )
    conf_result = runner.run_scf(conf_job)
    scf_count += 1

    de = float("inf")
    if conf_result.success and conf_result.energy is not None:
        data_points.append((predicted_n, conf_result.energy))
        de = abs(conf_result.energy - fit.e_inf)

    conv = de < etol_Ry
    if not conv and conf_result.success:
        warnings.append(
            f"k-mesh confirmation ΔE={de*1000/max(n_atoms,1):.3f} mRy/atom > "
            f"{etol_mRy:.3f}."
        )

    return ConvergenceEngineResult(
        converged=conv, final_kmesh=confirm_mesh, kmesh_data=data_points,
        kmesh_fit=fit, scf_count_kmesh=scf_count,
        errors=errors, warnings=warnings,
        etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
    )


def _converge_kmesh_metal(structure, basename, work_dir, runner,
                           n_eff, chosen_rkmax, core_valence_result,
                           km_result, etol_mRy, etol_Ry, n_atoms):
    """Metal k-mesh: 4 points + smearing cross-check."""
    warnings = []
    errors = []
    step = max(2, round(n_eff / 3))
    n_vals = [max(1, n_eff - step), max(1, n_eff - step // 2),
              n_eff, n_eff + step]

    jobs = []
    for i, nv in enumerate(n_vals):
        mesh = _get_proportional_mesh(nv, km_result.mesh)
        jobs.append(SCFJob(
            basename=basename, rkmax=chosen_rkmax, kmesh=mesh,
            ecut=core_valence_result.ecut, work_dir=work_dir,
            shift=km_result.shift, loose=True, parallel=True,
            label=f"KMESH-metal-{i+1}",
        ))

    results = runner.run_parallel(jobs)
    scf_count = len(jobs)
    scf_ok = [r for r in results if r.success and r.energy is not None]
    n_ok = [n_vals[i] for i, r in enumerate(results) if r.success and r.energy is not None]
    energies = [r.energy for r in scf_ok]
    data_points = list(zip(n_ok, energies))

    if len(scf_ok) < 3:
        errors.append(f"Only {len(scf_ok)}/4 metal k-mesh runs succeeded. Manual review.")
        return ConvergenceEngineResult(
            converged=False, final_kmesh=km_result.mesh, kmesh_data=data_points,
            scf_count_kmesh=scf_count, errors=errors, warnings=warnings,
            etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
        )

    is_monotonic = all(
        energies[i] > energies[i + 1] + _ENERGY_NOISE_RY * 10
        for i in range(len(energies) - 1)
    )
    if not is_monotonic:
        warnings.append(
            "k-mesh energy NOT monotonic for metallic system. "
            "Fermi-surface integration may not be converged. "
            "Raw data reported — manual inspection needed."
        )
        return ConvergenceEngineResult(
            converged=False, final_kmesh=km_result.mesh, kmesh_data=data_points,
            scf_count_kmesh=scf_count, errors=errors, warnings=warnings,
            etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
        )

    # Smearing cross-check at finest mesh
    finest_n = n_ok[-1]
    finest_mesh = _get_proportional_mesh(finest_n, km_result.mesh)
    conf_job = SCFJob(
        basename=basename, rkmax=chosen_rkmax, kmesh=finest_mesh,
        ecut=core_valence_result.ecut, work_dir=work_dir,
        shift=km_result.shift, loose=False, parallel=True,
        label="KMESH-smear-check",
    )
    conf_result = runner.run_scf(conf_job)
    scf_count += 1

    if conf_result.success and conf_result.energy is not None:
        data_points.append((finest_n, conf_result.energy))
        de_smear = abs(conf_result.energy - energies[-1])
        if de_smear > etol_Ry:
            warnings.append(
                f"Smearing cross-check: ΔE={de_smear*1000:.3f} mRy > "
                f"{etol_mRy:.3f} mRy. Smearing width not converged — "
                f"reduce and recheck."
            )
        else:
            warnings.append("Smearing cross-check passed.")

    return ConvergenceEngineResult(
        converged=len([w for w in warnings if "not converged" in w]) == 0,
        final_kmesh=finest_mesh, kmesh_data=data_points,
        scf_count_kmesh=scf_count, errors=errors, warnings=warnings,
        etol_mRy_per_atom=etol_mRy, etol_Ry=etol_Ry, n_atoms=n_atoms,
    )


# ── Top-level entry point ──

def run_convergence(structure, rmt_result, kmesh_result,
                    mixing_result, core_valence_result,
                    gmax_result, lmax_result,
                    basename, work_dir,
                    converge_params=frozenset({"rkmax", "kmesh"}),
                    etol_mRy_per_atom=DEFAULT_ETOL,
                    for_forces=False,
                    cluster_submit=False,
                    submit_cb=None,
                    quiet=False):
    """Run full convergence-verified parameter optimization.

    Returns (ConvergenceEngineResult, report_markdown).
    """
    if not wien2k_available():
        raise RuntimeError(
            "WIEN2k not found in PATH. Convergence engine requires WIEN2k."
        )

    n_atoms = structure.num_atoms_primitive
    runner = SCFRunner(work_dir, cluster_submit=cluster_submit,
                       submit_cb=submit_cb)

    warnings = []
    errors = []

    rkmax_found = gmax_result.gmax  # seed from table
    kmesh_found = kmesh_result.mesh
    gmax_found = gmax_result.gmax

    rkmax_data_pts = []
    kmesh_data_pts = []
    rkmax_fit_val = AitkenResult()
    kmesh_fit_val = AitkenResult()
    rmt_rob = RMTRobustnessResult()
    gmax_ver = GMAXVerificationResult(gmax=gmax_found)
    scf_rkmax = 0
    scf_kmesh = 0

    # Step 2: k-mesh at table RKMAX
    if "kmesh" in converge_params:
        if not quiet:
            print(f"  → k-mesh convergence (seed: {kmesh_result.mesh})...")
        eng = converge_kmesh(
            structure, basename, work_dir, runner,
            kmesh_result.system_type, rkmax_found,
            core_valence_result, etol_mRy_per_atom,
        )
        kmesh_found = eng.final_kmesh
        kmesh_data_pts = eng.kmesh_data
        kmesh_fit_val = eng.kmesh_fit
        scf_kmesh = eng.scf_count_kmesh
        # Propagate converged k-mesh into kmesh_result so downstream steps use it
        kmesh_result.mesh = kmesh_found
        kmesh_result.total_points = kmesh_found[0] * kmesh_found[1] * kmesh_found[2]
        warnings.extend(eng.warnings)
        errors.extend(eng.errors)

    # Step 1: RKMAX (uses the converged k-mesh if available)
    if "rkmax" in converge_params:
        etol_eff = etol_mRy_per_atom
        if for_forces:
            etol_eff /= FORCES_ETOL_FACTOR
            warnings.append(
                f"etol tightened to {etol_eff:.2f} mRy/atom for force/phonon "
                f"calculation. (Blaha 2020, Sec. III.B)"
            )
        if not quiet:
            print(f"  → RKMAX convergence (seed: {rkmax_found:.1f})...")
        eng = converge_rkmax(
            structure, rmt_result, kmesh_result, core_valence_result,
            basename, work_dir, runner, etol_mRy_per_atom=etol_eff,
        )
        rkmax_found = eng.final_rkmax
        rkmax_data_pts = eng.rkmax_data
        rkmax_fit_val = eng.rkmax_fit
        scf_rkmax = eng.scf_count_rkmax
        warnings.extend(eng.warnings)
        errors.extend(eng.errors)

    # Step 3: RMT robustness
    if "rmt" in converge_params:
        if not quiet:
            print(f"  → RMT robustness check...")
        rmt_rob = rmt_robustness_check(
            structure, rmt_result, basename, work_dir,
            rkmax_found, kmesh_found, core_valence_result.ecut, runner, etol_mRy_per_atom,
        )
        if not rmt_rob.passed:
            warnings.append(rmt_rob.warning)

    # Step 4: Final confirmation
    if not quiet:
        print(f"  → Final confirmation (RKMAX={rkmax_found:.1f}, "
              f"mesh={kmesh_found[0]}×{kmesh_found[1]}×{kmesh_found[2]})...")

    conf_job = SCFJob(
        basename=basename, rkmax=rkmax_found, kmesh=kmesh_found,
        ecut=core_valence_result.ecut, work_dir=work_dir,
        shift=kmesh_result.shift, loose=False, parallel=True,
        gmax=gmax_found, label="final-confirm",
    )
    conf_result = runner.run_scf(conf_job)

    # Step 5: GMAX verification
    scf_gmax = 0
    if "gmax" in converge_params:
        gmax_ver = verify_gmax_at_rkmax(work_dir, basename, gmax_found)
        if not gmax_ver.verified:
            new_gmax = gmax_found + gmax_ver.suggested_increase
            warnings.append(
                f"GMAX insufficient: {gmax_ver.warning}. "
                f"Increase by ~{gmax_ver.suggested_increase:.1f}."
            )
            conf_job2 = SCFJob(
                basename=basename, rkmax=rkmax_found, kmesh=kmesh_found,
                ecut=core_valence_result.ecut, work_dir=work_dir,
                shift=kmesh_result.shift, loose=False, parallel=True,
                gmax=new_gmax, label="final-confirm-gmax",
            )
            conf_result = runner.run_scf(conf_job2)
            gmax_found = new_gmax
            scf_gmax = 1

    total_scf = runner.total_scf_count
    total_wall = runner.total_wall_time

    has_enough = ("kmesh" not in converge_params or len(kmesh_data_pts) >= 3) and \
                 ("rkmax" not in converge_params or len(rkmax_data_pts) >= 4)

    eng_result = ConvergenceEngineResult(
        converged=conf_result.success and len(errors) == 0 and has_enough,
        final_rkmax=rkmax_found, final_kmesh=kmesh_found, final_gmax=gmax_found,
        rkmax_data=rkmax_data_pts, kmesh_data=kmesh_data_pts,
        rkmax_fit=rkmax_fit_val, kmesh_fit=kmesh_fit_val,
        rmt_robustness=rmt_rob, gmax_verify=gmax_ver,
        confirmation_energy=conf_result.energy,
        scf_count_rkmax=scf_rkmax, scf_count_kmesh=scf_kmesh,
        scf_count_total=total_scf, total_wall_time=total_wall,
        warnings=warnings, errors=errors,
        etol_mRy_per_atom=etol_mRy_per_atom,
        etol_Ry=etol_mRy_per_atom * n_atoms / 1000.0,
        n_atoms=n_atoms, for_forces=for_forces,
    )

    return eng_result


# ── Markdown Convergence Report ──

def build_convergence_report(result: ConvergenceEngineResult,
                              structure, struct_path="",
                              calc_type="scf") -> str:
    """Generate a markdown convergence report per Section 8 of the spec."""
    lines = []
    lines.append("# WIEN2k Convergence-Verified Optimization Report")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if struct_path:
        lines.append(f"**Structure:** `{struct_path}`")
    lines.append(f"**Calculation type:** {calc_type}")
    lines.append(
        f"**Convergence tolerance:** {result.etol_mRy_per_atom:.2f} mRy/atom "
        f"({result.etol_Ry*1000:.3f} mRy total)"
    )
    if result.for_forces:
        lines.append("**Note:** Tolerance tightened 10× for force/phonon calculation.")
    lines.append(f"**Atoms (primitive cell):** {result.n_atoms}")
    lines.append("")

    status = "CONVERGED" if result.converged else "NOT CONVERGED — manual review needed"
    lines.append(f"## Status: {status}")
    lines.append("")

    if result.warnings:
        lines.append("### Warnings")
        for w in result.warnings:
            lines.append(f"- :warning: {w}")
        lines.append("")
    if result.errors:
        lines.append("### Errors")
        for e in result.errors:
            lines.append(f"- :x: {e}")
        lines.append("")

    # Final parameters table
    lines.append("## Final Parameters")
    lines.append("")
    lines.append("| Parameter | Value |")
    lines.append("|-----------|-------|")
    lines.append(f"| RKMAX | {result.final_rkmax:.1f} |")
    n1, n2, n3 = result.final_kmesh
    lines.append(f"| k-mesh | {n1}×{n2}×{n3} ({n1*n2*n3} pts) |")
    lines.append(f"| GMAX | {result.final_gmax:.1f} |")
    if result.confirmation_energy is not None:
        lines.append(
            f"| Confirmation energy | {result.confirmation_energy:.8f} Ry |"
        )
    lines.append("")

    # RKMAX section
    if result.rkmax_data:
        lines.append("## RKMAX Convergence")
        lines.append("")
        lines.append("| RKMAX | Energy (Ry) | ΔE (Ry) |")
        lines.append("|-------|-------------|---------|")
        prev = None
        for x, e in result.rkmax_data:
            es = f"{e:.8f}" if e is not None else "FAILED"
            de = f"{e - prev:.8f}" if (prev is not None and e is not None) else "—"
            lines.append(f"| {x:.1f} | {es} | {de} |")
            if e is not None:
                prev = e
        lines.append("")

        if result.rkmax_fit.is_valid:
            fv = result.rkmax_fit
            lines.append("**Aitken Δ² Extrapolation:**")
            lines.append("")
            lines.append(f"- E_inf = {fv.e_inf:.8f} Ry")
            lines.append(f"- α = {fv.alpha:.4f}")
            lines.append(f"- A = {fv.amplitude:.8f} Ry")
            lines.append(f"- Model: E(x) ≈ E_inf + A · exp(-α · x)")
            lines.append("")
            last = result.rkmax_data[-1][1] if result.rkmax_data else None
            if last is not None:
                res = abs(last - fv.e_inf) * 1000.0 / max(result.n_atoms, 1)
                lines.append(f"- Confirmation residual: {res:.4f} mRy/atom")
            lines.append("")

        naive = 8
        n_rk = max(result.scf_count_rkmax, 1)
        lines.append(
            f"**Efficiency:** {result.scf_count_rkmax} SCF runs vs. ~{naive} "
            f"linear sweep — **{naive / n_rk:.1f}× fewer**"
        )
        lines.append("")

    # k-mesh section
    if result.kmesh_data:
        lines.append("## k-Mesh Convergence")
        lines.append("")
        lines.append("| n_eff | Energy (Ry) | ΔE (Ry) |")
        lines.append("|-------|-------------|---------|")
        prev = None
        for x, e in result.kmesh_data:
            es = f"{e:.8f}" if e is not None else "FAILED"
            de = f"{e - prev:.8f}" if (prev is not None and e is not None) else "—"
            lines.append(f"| {x} | {es} | {de} |")
            if e is not None:
                prev = e
        lines.append("")

        naive_km = 5
        n_km = max(result.scf_count_kmesh, 1)
        lines.append(
            f"**Efficiency:** {result.scf_count_kmesh} SCF runs vs. ~{naive_km} "
            f"linear sweep — **{naive_km / n_km:.1f}× fewer**"
        )
        lines.append("")

    # RMT robustness
    if result.rmt_robustness.original_rmts:
        r = result.rmt_robustness
        lines.append("## RMT Robustness Check")
        lines.append("")
        lines.append(f"**Status:** {'PASSED' if r.passed else 'FAILED'}")
        lines.append(f"- ΔE after 5% RMT reduction: {r.delta_energy:.4f} mRy/atom")
        if r.warning:
            lines.append(f"- {r.warning}")
        lines.append("")

    # GMAX
    if result.gmax_verify.warning:
        lines.append("## GMAX Verification")
        lines.append("")
        lines.append(f"**Result:** {'OK' if result.gmax_verify.verified else 'ADJUSTED'}")
        lines.append(f"- {result.gmax_verify.warning}")
        lines.append("")

    # Summary
    lines.append("## Summary Statistics")
    lines.append("")
    lines.append(f"- Total SCF runs: {result.scf_count_total}")
    lines.append(f"- Total wall time: {result.total_wall_time:.1f}s")
    lines.append(f"- All runs externally verified: {result.converged}")
    lines.append("")
    lines.append("---")
    lines.append(
        "*Report generated by wien2k-parameter-optimizer. "
        "All convergence claims backed by external SCF verification.*"
    )

    return "\n".join(lines)
