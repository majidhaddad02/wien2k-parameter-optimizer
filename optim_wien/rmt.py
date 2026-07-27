"""
RMT Optimization — Four Strict Conditions.

1. Non-overlapping spheres with safety margin
2. Core leakage control (NEC01)
3. Relative size ratio ≤ 1.5
4. Structural change allowance for relaxation
"""

import math
from dataclasses import dataclass, field

from .constants import (
    ELEMENT_TYPE, INITIAL_RMT, SAFETY_MARGIN_SCF, SAFETY_MARGIN_RELAX,
    MAX_RMT_RATIO, MAX_RMT_RATIO_STRICT,
    RMT_REDUCTION_RELAX, CORE_LEAKAGE_OK, CORE_LEAKAGE_WARN, CORE_LEAKAGE_CRITICAL,
    RMT_LARGE_THRESHOLD, CalcType, Precision,
)
from .struct_parser import compute_pairwise_min_distances


@dataclass
class RMTOptimizationResult:
    rmt_values: list = field(default_factory=list)
    pairwise_distances: dict = field(default_factory=dict)
    overlap_warnings: list = field(default_factory=list)
    core_leakage_warnings: list = field(default_factory=list)
    ratio_warnings: list = field(default_factory=list)
    structural_notes: list = field(default_factory=list)
    hydrogen_notes: list = field(default_factory=list)
    hdlo_notes: list = field(default_factory=list)
    same_element_notes: list = field(default_factory=list)


def optimize_rmt(
    structure,
    calc_type=CalcType.SCF,
    precision=Precision.MEDIUM,
    core_leakage_data=None,
):
    if core_leakage_data is None:
        core_leakage_data = {}

    atoms = structure.atoms
    n = len(atoms)
    rmt = [a.rmt for a in atoms]
    pw = compute_pairwise_min_distances(structure)

    result = RMTOptimizationResult()
    safety = SAFETY_MARGIN_SCF if calc_type == CalcType.SCF else SAFETY_MARGIN_RELAX

    _init_from_table(atoms, rmt, result)
    _nonoverlap(rmt, atoms, pw, safety)
    _core_leakage(rmt, atoms, core_leakage_data, result)
    _ratio_balance(rmt, atoms, pw, safety, result)
    _structural_margin(rmt, atoms, calc_type, result)
    _same_elements(rmt, atoms, result)
    _handle_hydrogen(rmt, atoms, pw, result)
    _nonoverlap(rmt, atoms, pw, safety)
    _check_large(rmt, atoms, result)
    _validate(rmt, atoms, pw, safety, result)

    result.rmt_values = [round(r, 4) for r in rmt]
    result.pairwise_distances = pw
    return result


def _init_from_table(atoms, rmt, result):
    for i, a in enumerate(atoms):
        tv = INITIAL_RMT.get(a.element, 1.8)
        if rmt[i] < 1.0:
            rmt[i] = tv


def _nonoverlap(rmt, atoms, pw, safety):
    n = len(atoms)
    for _ in range(50):
        ok = True
        for i in range(n):
            for j in range(i + 1, n):
                nn = pw.get((i, j), 10.0)
                if rmt[i] + rmt[j] > safety * nn + 1e-10:
                    ok = False
                    exc = rmt[i] + rmt[j] - safety * nn + 0.002
                    if rmt[i] >= rmt[j]:
                        rmt[i] = max(0.3, rmt[i] - exc * 0.55)
                        rmt[j] = max(0.3, rmt[j] - exc * 0.45)
                    else:
                        rmt[j] = max(0.3, rmt[j] - exc * 0.55)
                        rmt[i] = max(0.3, rmt[i] - exc * 0.45)
        if ok:
            break


def _core_leakage(rmt, atoms, cld, result):
    for i, a in enumerate(atoms):
        leak = cld.get(a.element, 0.0)
        if leak > CORE_LEAKAGE_CRITICAL:
            result.core_leakage_warnings.append(
                f"CRITICAL: {a.element}({i}) core leakage {leak:.4f}e⁻ "
                f"> {CORE_LEAKAGE_CRITICAL}. Increase RMT or transfer state "
                f"to valence."
            )
            rmt[i] *= 1.05
        elif leak > CORE_LEAKAGE_WARN:
            result.core_leakage_warnings.append(
                f"WARNING: {a.element}({i}) core leakage {leak:.4f}e⁻ "
                f"> {CORE_LEAKAGE_WARN}. Increasing RMT by 5%."
            )
            rmt[i] *= 1.05
        elif leak > CORE_LEAKAGE_OK:
            result.core_leakage_warnings.append(
                f"Note: {a.element}({i}) core leakage {leak:.4f}e⁻ (marginal)."
            )


def _ratio_balance(rmt, atoms, pw, safety, result):
    n = len(atoms)
    for _ in range(30):
        rmax = max(rmt)
        rmin = max(min(rmt), 0.5)
        ratio = rmax / rmin
        if ratio <= MAX_RMT_RATIO:
            break

        max_idx = rmt.index(rmax)
        min_idx = rmt.index(rmin)
        max_type = ELEMENT_TYPE.get(atoms[max_idx].element, "sp")
        min_type = ELEMENT_TYPE.get(atoms[min_idx].element, "sp")

        if max_type in ("d", "f"):
            rmt[max_idx] = max(0.3, rmt[max_idx] * 0.95)
            if min_type == "sp":
                nn = pw.get((max_idx, min_idx), 10.0)
                max_allowed = safety * nn - rmt[max_idx] - 0.01
                rmt[min_idx] = min(rmt[min_idx] * 1.02, max(max_allowed, 0.5))
        else:
            rmt[max_idx] = max(0.3, rmt[max_idx] * 0.93)

    final_ratio = max(rmt) / max(min(rmt), 0.5)
    if final_ratio <= MAX_RMT_RATIO_STRICT:
        result.ratio_warnings.append(
            f"RMT ratio {final_ratio:.2f} ≤ 1.3 (OK)."
        )
    elif final_ratio <= MAX_RMT_RATIO:
        result.ratio_warnings.append(
            f"RMT ratio {final_ratio:.2f} ≤ 1.5 (acceptable)."
        )
    else:
        result.ratio_warnings.append(
            f"WARNING: RMT ratio {final_ratio:.2f} > 1.5. "
            f"Ghostbands/QTL-B possible. Consider HDLOs."
        )


def _structural_margin(rmt, atoms, calc_type, result):
    if calc_type in (CalcType.RELAXATION, CalcType.OPTIMIZATION, CalcType.EOS):
        for i in range(len(atoms)):
            rmt[i] *= RMT_REDUCTION_RELAX
        result.structural_notes.append(
            f"Reduced RMTs by 7% for {calc_type.value}."
        )


def _same_elements(rmt, atoms, result):
    groups = {}
    for i, a in enumerate(atoms):
        groups.setdefault(a.element, []).append((i, rmt[i]))
    for el, entries in groups.items():
        if len(entries) > 1:
            minv = min(e[1] for e in entries)
            for idx, _ in entries:
                rmt[idx] = minv
            result.same_element_notes.append(
                f"Unified RMT for {el} sites: {minv:.4f} bohr."
            )


def _handle_hydrogen(rmt, atoms, pw, result):
    for i, a in enumerate(atoms):
        if a.element != "H":
            continue
        min_d = float("inf")
        min_j = -1
        for j in range(len(atoms)):
            if i == j:
                continue
            d = pw.get((i, j), 10.0)
            if d < min_d:
                min_d = d
                min_j = j
        if min_j >= 0 and min_d < 2.5:
            partner_rmt = rmt[min_j]
            suggested = min(0.5 * partner_rmt, 0.8)
            if rmt[i] > suggested:
                rmt[i] = suggested
                result.hydrogen_notes.append(
                    f"H({i}): short bond to {atoms[min_j].element} "
                    f"({min_d:.2f} bohr). "
                    f"RMT(H) = {suggested:.3f}."
                )
        if rmt[i] < 0.5:
            rmt[i] = 0.6
            result.hydrogen_notes.append(f"H({i}): floor RMT = 0.6 bohr.")


def _check_large(rmt, atoms, result):
    for i, a in enumerate(atoms):
        if rmt[i] > 2.5:
            result.hdlo_notes.append(
                f"{a.element}({i}): RMT={rmt[i]:.2f} > 2.5. "
                f"HDLOs strongly recommended."
            )
        elif rmt[i] > RMT_LARGE_THRESHOLD:
            et = ELEMENT_TYPE.get(a.element, "sp")
            if et in ("d", "f"):
                result.hdlo_notes.append(
                    f"{a.element}({i}): RMT={rmt[i]:.2f} > 2.2 "
                    f"({et}-element). HDLOs recommended."
                )


def _validate(rmt, atoms, pw, safety, result):
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            nn = pw.get((i, j), 10.0)
            if rmt[i] + rmt[j] > nn + 1e-10:
                result.overlap_warnings.append(
                    f"FINAL WARNING: {atoms[i].element}({i}) + "
                    f"{atoms[j].element}({j}) RMT sum "
                    f"{rmt[i]+rmt[j]:.3f} > NN {nn:.3f}. Manual review!"
                )
