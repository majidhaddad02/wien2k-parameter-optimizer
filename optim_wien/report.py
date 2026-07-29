"""
Comprehensive Optimization Report Generator.

Generates detailed reports covering all WIEN2k parameters:
RMT, RKMAX, GMAX, LMAX/LVNS, k-mesh, mixing, core/valence.
"""

import os
from datetime import datetime

from .constants import ELEMENT_TYPE, REFERENCES

_SEP = "=" * 80
_SUB = "-" * 80


def generate_report(
    structure,
    rmt_result,
    rkmax_result,
    gmax_result,
    lmax_result,
    kmesh_result,
    mixing_result,
    core_valence_result,
    calc_type="scf",
    precision="medium",
    struct_path="",
):
    lines = []
    lines.append(_header(struct_path, calc_type, precision))
    lines.append(_rmt_section(structure, rmt_result))
    lines.append(_rkmax_section(structure, rmt_result, rkmax_result))
    lines.append(_gmax_section(gmax_result))
    lines.append(_lmax_section(structure, rmt_result, lmax_result))
    lines.append(_kmesh_section(structure, kmesh_result))
    lines.append(_mixing_section(mixing_result))
    lines.append(_core_valence_section(core_valence_result))
    lines.append(_convergence_section(kmesh_result, rkmax_result))
    lines.append(_summary_section(structure, rmt_result, rkmax_result,
                                   kmesh_result, gmax_result,
                                   core_valence_result, mixing_result))
    lines.append(_references())
    return "\n".join(lines)


def _header(path, ct, pr):
    return "\n".join([
        _SEP,
        "WIEN2k Comprehensive Parameter Optimization Report",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Structure: {path}" if path else "",
        f"Calculation: {ct}  |  Precision: {pr}",
        "",
        "Based on: WIEN2k FAQ (Blaha, Schwarz, Luitz) | J. Chem. Phys. 152, 074101 (2020)",
        _SEP,
    ])


def _rmt_section(structure, rr):
    lines = ["", "1. RMT OPTIMIZATION", _SUB]
    hdr = f"{'#':>3}  {'El':<6} {'Type':<4} {'RMT(bohr)':>10}  {'Notes':<40}"
    lines.append(hdr)
    lines.append("-" * len(hdr))

    for i, a in enumerate(structure.atoms):
        rm = rr.rmt_values[i]
        et = ELEMENT_TYPE.get(a.element, "sp")
        notes = ""
        if rr.hdlo_notes:
            for n in rr.hdlo_notes:
                if n.startswith(f"{a.element}({i})"):
                    notes = "HDLO recommended"
                    break
        lines.append(f"{i+1:>3}  {a.element:<6} {et:<4} {rm:>10.4f}  {notes:<40}")

    ratio = max(rr.rmt_values) / max(min(rr.rmt_values), 0.5)
    lines.append(f"\nRMT ratio (max/min): {ratio:.2f}")

    for w in rr.ratio_warnings:
        lines.append(f"  {w}")
    for w in rr.core_leakage_warnings:
        lines.append(f"  {w}")
    for w in rr.structural_notes:
        lines.append(f"  {w}")
    for w in rr.hydrogen_notes:
        lines.append(f"  {w}")
    for w in rr.same_element_notes:
        lines.append(f"  {w}")
    for w in rr.overlap_warnings:
        lines.append(f"  {w}")

    return "\n".join(lines)


def _rkmax_section(structure, rr, kr):
    lines = ["", "2. RKMAX OPTIMIZATION", _SUB]
    lines.append(f"Smallest atom: {kr.min_element} (RMT = {kr.r_min:.2f} bohr)")
    lines.append(f"Base RKMAX (Blaha table): {kr.base_rkmax}")
    lines.append(f"Recommended RKMAX: {kr.rkmax}")
    lines.append("")
    lines.append(f"  {'Atom':>8} {'RMT':>8}  {'RKMAX_eff':>10}")
    lines.append(f"  {'-'*8:>8} {'-'*8:>8}  {'-'*10:>10}")

    for i, a in enumerate(structure.atoms):
        ek = kr.effective_rkmax[i] if i < len(kr.effective_rkmax) else 0
        rm = rr.rmt_values[i]
        lines.append(f"  {a.element:>8} {rm:>8.3f}  {ek:>10.2f}")

    if kr.hydrogen_note:
        lines.append(f"\n  {kr.hydrogen_note}")
    for w in kr.warnings:
        lines.append(f"  {w}")

    return "\n".join(lines)


def _gmax_section(gr):
    lines = ["", "3. GMAX (Fourier Expansion)", _SUB]
    lines.append(f"GMAX: {gr.gmax}")
    lines.append(f"NR2V: {gr.nr2v}")
    lines.append(f"Rationale: {gr.rationale}")
    return "\n".join(lines)


def _lmax_section(structure, rr, lr):
    lines = ["", "4. LMAX / LVNS", _SUB]
    lines.append(f"Global LVNS: {lr.lvns_global}")
    lines.append("")
    for i, a in enumerate(structure.atoms):
        lmax = lr.lmax_values[i] if i < len(lr.lmax_values) else 6
        rm = rr.rmt_values[i]
        lines.append(f"  {a.element}({i}): LMAX={lmax}, RMT={rm:.3f}")
    for n in lr.notes:
        lines.append(f"  {n}")
    for h in lr.hdlo_recommendations:
        lines.append(f"  {h}")
    return "\n".join(lines)


def _kmesh_section(structure, kr):
    lines = ["", "5. K-MESH", _SUB]
    lines.append(f"Cell volume: {structure.volume:.1f} bohr³")
    lines.append(f"Atoms (primitive): {structure.num_atoms_primitive}")
    lines.append(f"System type: {kr.system_type}")
    n1, n2, n3 = kr.mesh
    lines.append(f"Mesh: {n1}×{n2}×{n3} ({kr.total_points} points)")
    lines.append(f"Shift: {kr.shift}")
    lines.append(f"Inversion: {'YES' if kr.add_inversion else 'NO'}")
    for n in kr.notes:
        lines.append(f"  {n}")
    return "\n".join(lines)


def _mixing_section(mr):
    lines = ["", "6. MIXING & SCF", _SUB]
    lines.append(f"Scheme: {mr.scheme}")
    lines.append(f"Mixing factor: {mr.mixing_factor:.2f}")
    lines.append(f"TEMP: {mr.temp:.4f} Ry")
    lines.append(f"TETRA: {'YES' if mr.tetra else 'NO'}")
    lines.append(f"Energy convergence: {mr.econv:.6f} Ry")
    lines.append(f"Charge convergence: {mr.cconv:.6f}")
    lines.append(f"Max SCF cycles: {mr.max_scf}")
    for n in mr.notes:
        lines.append(f"  {n}")
    return "\n".join(lines)


def _core_valence_section(cr):
    lines = ["", "7. CORE/VALENCE", _SUB]
    lines.append(f"Ecut: {cr.ecut:.1f} Ry")
    lines.append(f"HDLOs: {'YES' if cr.use_hdlo else 'NO'}")
    for n in cr.notes:
        lines.append(f"  {n}")
    return "\n".join(lines)


def _convergence_section(kr, rkr):
    lines = ["", "8. CONVERGENCE STRATEGY (Blaha method)", _SUB]
    steps = [
        ("k-mesh", f"Start with coarser mesh, compare ΔE < 0.1 mRy"),
        ("RKMAX", f"Test RKMAX ±0.5 around {rkr.rkmax}, compare ΔE"),
        ("Forces", "Test RKMAX ±0.5, compare ΔF < 1 mRy/bohr"),
        ("k-mesh (forces)", "Refine k-mesh, compare ΔF"),
        ("Volume", "Test RKMAX ±1.0, compare ΔV/V < 0.1%"),
        ("Properties", "Gap, EFG — converge to target accuracy"),
    ]
    for title, desc in steps:
        lines.append(f"  {title}: {desc}")

    lines.append("\nk-point convergence:")
    for s in kr.convergence_steps:
        lines.append(f"  {s}")

    lines.append("\nIMPORTANT (Blaha): Total energy converges faster than forces!")
    lines.append("  Converge the PROPERTY OF INTEREST, not just total energy.")
    return "\n".join(lines)


def _summary_section(structure, rr, rkr, kr, gr, cr, mr):
    lines = ["", "9. QUICK SUMMARY", _SUB]
    lines.append(f"RMT: {', '.join(f'{a.element}={rr.rmt_values[i]:.3f}' for i, a in enumerate(structure.atoms))}")
    lines.append(f"RKMAX: {rkr.rkmax}  |  GMAX: {gr.gmax}")
    n1, n2, n3 = kr.mesh
    lines.append(f"k-mesh: {n1}×{n2}×{n3} ({kr.system_type}, {kr.total_points} pts)")
    lines.append(f"Ecut: {cr.ecut:.1f} Ry  |  Mixing: {mr.scheme} {mr.mixing_factor:.2f}")
    lines.append("")
    lines.append(f"Recommended WIEN2k init: "
                 f"init_lapw -b -rkmax {rkr.rkmax} "
                 f"-numk {n1*n2*n3} "
                 f"-ecut {int(abs(cr.ecut))}")
    return "\n".join(lines)


def _references():
    lines = ["", "10. REFERENCES", _SUB]
    for r in REFERENCES:
        lines.append(f"  {r}")
    lines.append("")
    lines.append(_SEP)
    return "\n".join(lines)


def write_optimized_struct(structure, rmt_values, filepath):
    lines = [structure.title, structure.lattice_type, structure.mode]
    lines.append(
        f"             {structure.a:12.6f} {structure.b:12.6f} "
        f"{structure.c:12.6f} {structure.alpha:6.1f} "
        f"{structure.beta:6.1f} {structure.gamma:6.1f}"
    )
    for i, a in enumerate(structure.atoms):
        rm = rmt_values[i] if i < len(rmt_values) else a.rmt
        lines.append(f"ATOM {i+1:3d}: X={a.position[0]:.10f} "
                     f"Y={a.position[1]:.10f} Z={a.position[2]:.10f}")
        lines.append(f"          MULT= {a.mult:2d}          ISPLIT= {a.isplit:2d}")
        for nidx, ep in enumerate(a.equivalent_positions, start=2):
            sign = "-" if nidx > 2 else " "
            lines.append(
                f"      -{nidx}: X={ep[0]:.10f} Y={ep[1]:.10f} Z={ep[2]:.10f}"
            )
        lines.append(f"{a.element:10}        NPT=  781  "
                     f"R0=0.00010000 RMT=    {rm:10.5f}   Z: {a.z:.1f}")
        rot_lines = a.local_rot_lines
        if rot_lines and len(rot_lines) >= 3:
            for rl in rot_lines[:3]:
                lines.append(rl)
        else:
            lines.append("LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000")
            lines.append("                     0.0000000 1.0000000 0.0000000")
            lines.append("                     0.0000000 0.0000000 1.0000000")
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")
