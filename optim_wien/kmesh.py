"""
k-mesh Optimization — Adaptive Monkhorst-Pack Grid.

Based on WIEN2k FAQ kgen and Blaha's convergence strategy.
"""

import math
from dataclasses import dataclass, field

from .constants import ELEMENT_TYPE


def _detect_system(structure):
    atoms = structure.atoms
    n = structure.num_atoms_primitive
    elements = [a.element for a in atoms]
    has_tm = any(ELEMENT_TYPE.get(el, "sp") in ("d", "f") for el in elements)
    has_oxide = any(el in ("O", "F", "S", "Cl", "Se", "Br", "Te", "I")
                    for el in elements)
    has_nitride = any(el in ("N", "C", "P", "As", "B") for el in elements)
    is_hex = abs(structure.gamma - 120.0) < 1e-6

    if n >= 40:
        return "metal_large" if has_tm else "insulator_large"

    if has_tm and has_nitride and not has_oxide and n <= 4:
        return "metal_small"

    if has_tm and has_oxide and n <= 20:
        return "semiconductor"

    if has_tm and n <= 10:
        return "metal_small"

    if n <= 10:
        return "semiconductor"

    return "semiconductor"


_K_DENSITY = {
    "metal_small": 3000, "semiconductor": 500, "insulator": 200,
    "metal_large": 300, "insulator_large": 10,
}

_REFINE = {"coarse": 0.5, "medium": 1.0, "fine": 2.0, "very_fine": 4.0}


@dataclass
class KMeshResult:
    mesh: tuple = (6, 6, 6)
    shift: tuple = (0.5, 0.5, 0.5)
    gamma_centered: bool = False
    add_inversion: bool = True
    total_points: int = 0
    system_type: str = "semiconductor"
    convergence_steps: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def _primitive_volume_factor(lattice_type):
    lt = lattice_type.strip().upper()
    if lt in ("F", "FCC", "CF"):
        return 0.25
    if lt in ("B", "BCC", "BC"):
        return 0.50
    if lt in ("CXY", "CYZ", "CXZ", "C"):
        return 0.50
    if lt in ("R",):
        return 1.0 / 3.0
    return 1.0


def optimize_kmesh(structure, refinement="medium", system_type=None,
                   gamma_centered=None, bandgap=None):
    result = KMeshResult()

    st = system_type or _detect_system(structure)

    if bandgap is not None:
        if bandgap <= 0.0:
            st = "metal_small"
            result.notes.append(f"bandgap={bandgap:.1f} eV → treating as metal.")
        elif bandgap < 1.0:
            st = "semiconductor"
            result.notes.append(f"bandgap={bandgap:.1f} eV → treating as semiconductor.")
        else:
            st = "insulator"
            result.notes.append(f"bandgap={bandgap:.1f} eV → treating as insulator.")

    result.system_type = st

    if st in ("surface",):
        result.mesh = (20, 20, 1)
        result.shift = (0.5, 0.5, 0.0)
        result.total_points = 400
        result.gamma_centered = False
        result.notes.append("Surface/slab: 20×20×1.")
        return result
    if st in ("molecule",):
        result.mesh = (1, 1, 1)
        result.shift = (0.0, 0.0, 0.0)
        result.total_points = 1
        result.gamma_centered = True
        result.notes.append("Molecule: Gamma-point only.")
        return result

    density = _K_DENSITY.get(st, 200)
    V = structure.volume * _primitive_volume_factor(structure.lattice_type)
    Vbz = (2 * math.pi)**3 / max(V, 1e-10)
    target_total = max(1, density * Vbz)

    bl = structure.reciprocal_lengths
    b1, b2, b3 = bl
    product = max(b1 * b2 * b3, 1e-10)

    scale = (target_total / product) ** (1.0 / 3.0)
    n1 = max(1, round(b1 * scale))
    n2 = max(1, round(b2 * scale))
    n3 = max(1, round(b3 * scale))

    factor = _REFINE.get(refinement, 1.0)
    n1 = max(1, round(n1 * factor))
    n2 = max(1, round(n2 * factor))
    n3 = max(1, round(n3 * factor))
    if n1 + n2 + n3 == 3:
        n1 = 2

    result.mesh = (n1, n2, n3)
    result.total_points = n1 * n2 * n3

    is_hex = abs(structure.gamma - 120.0) < 1e-6
    if gamma_centered is None:
        gamma_centered = is_hex
    result.gamma_centered = gamma_centered
    result.shift = (0.0, 0.0, 0.0) if gamma_centered else (0.5, 0.5, 0.5)
    result.add_inversion = True

    if is_hex:
        result.notes.append("Hexagonal cell detected — using Gamma-centered mesh.")

    if not gamma_centered and st in ("semiconductor", "insulator"):
        result.notes.append(
            "⚠ Shifted mesh may miss gap at Γ/X for semiconductors. "
            "Use Gamma-centered mesh for accurate band gaps."
        )

    _build_steps(result, st, n1, n2, n3)
    return result


def _build_steps(result, st, n1, n2, n3):
    c1, c2, c3 = max(1, n1//2), max(1, n2//2), max(1, n3//2)
    f1 = n1 + max(1, n1//2)
    f2 = n2 + max(1, n2//2)
    f3 = n3 + max(1, n3//2)
    v1, v2, v3 = n1 * 2, n2 * 2, n3 * 2

    result.convergence_steps = [
        f"1. SCF with {c1}×{c2}×{c3} k-mesh",
        f"2. Refine to {n1}×{n2}×{n3}, compare ΔE < 0.1 mRy",
        f"3. Further refine to {f1}×{f2}×{f3}, compare forces",
        f"4. If needed, use {v1}×{v2}×{v3} for DOS/optics",
    ]
    if "metal" in st:
        result.convergence_steps.append(
            "NOTE: metallic — use Methfessel-Paxton (TETRA 101)."
        )
