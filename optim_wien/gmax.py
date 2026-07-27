"""
GMAX Optimization — Fourier Expansion of Density/Potential.

Based on WIEN2k User's Guide, Section 4.2 (case.in0).
"""

from dataclasses import dataclass, field

from .constants import GMAX_DEFAULT, GMAX_HYDROGEN, GMAX_PRECISION, Precision


@dataclass
class GMAXOptimizationResult:
    gmax: float = 12.0
    nr2v: int = 1
    rationale: str = ""


def optimize_gmax(atoms, rmt_values, precision=Precision.MEDIUM):
    result = GMAXOptimizationResult()

    has_h = any(a.element == "H" and rmt_values[i] < 0.8
                for i, a in enumerate(atoms))
    has_li = any(a.element == "Li" and rmt_values[i] < 1.2
                 for i, a in enumerate(atoms))
    has_halogen = any(a.element in ("F", "Cl", "Br", "I") for a in atoms)
    has_f_el = any(a.element in (
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu",
    ) for a in atoms)
    has_tm = any(a.element in (
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    ) for a in atoms)

    if has_h:
        result.gmax = GMAX_HYDROGEN
        result.rationale = "H with small RMT → GMAX = 20."
    elif has_li:
        result.gmax = 16.0
        result.rationale = "Li with small RMT → GMAX = 16."
    elif has_halogen:
        result.gmax = max(GMAX_PRECISION.get(precision, 14), 16.0)
        result.rationale = "Halogen present → GMAX ≥ 16."
    elif has_f_el:
        result.gmax = 16.0
        result.rationale = "f-elements → GMAX = 16."
    else:
        result.gmax = GMAX_PRECISION.get(precision, GMAX_DEFAULT)
        result.rationale = f"Standard GMAX = {result.gmax} (precision={precision.value})."

    if precision in (Precision.VERY_HIGH,):
        result.nr2v = 4
    elif precision in (Precision.HIGH,):
        result.nr2v = 2
    else:
        result.nr2v = 1

    return result
