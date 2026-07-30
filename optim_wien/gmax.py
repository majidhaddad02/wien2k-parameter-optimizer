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


def optimize_gmax(atoms, rmt_values, precision=Precision.MEDIUM, strict_faq=False):
    result = GMAXOptimizationResult()

    if strict_faq:
        result.gmax = 12.0
        result.rationale = "Strict-FAQ mode: GMAX = 12.0 (WIEN2k default)."
        return result

    has_h_small = any(a.element == "H" and rmt_values[i] < 0.8
                      for i, a in enumerate(atoms))
    has_li_small = any(a.element == "Li" and rmt_values[i] < 1.2
                       for i, a in enumerate(atoms))
    has_halogen_small = any(
        a.element in ("F", "Cl", "Br", "I") and rmt_values[i] < 2.0
        for i, a in enumerate(atoms)
    )
    has_f_small = any(
        a.element in (
            "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
            "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu",
        ) and rmt_values[i] < 2.5
        for i, a in enumerate(atoms)
    )

    if has_h_small:
        result.gmax = GMAX_HYDROGEN
        result.rationale = "H with small RMT (<0.8) → GMAX = 20."
    elif has_li_small:
        result.gmax = 16.0
        result.rationale = "Li with small RMT (<1.2) → GMAX = 16."
    elif has_halogen_small:
        result.gmax = max(GMAX_PRECISION.get(precision, 14), 16.0)
        result.rationale = "Halogen with small RMT (<2.0) → GMAX ≥ 16."
    elif has_f_small:
        result.gmax = 16.0
        result.rationale = "f-element with small RMT (<2.5) → GMAX = 16."
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
