"""
RKMAX Optimization — Blaha's Reference Table.

RKMAX = R_min × K_max
Computational cost ∝ RKMAX³
"""

from dataclasses import dataclass, field

from .constants import RKMAX_TABLE, RKMAX_OFFSET, RKMAX_OFFSET_H_SMALL, Precision


@dataclass
class RKMAXOptimizationResult:
    rkmax: float = 7.0
    r_min: float = 2.0
    min_element: str = ""
    base_rkmax: float = 7.0
    effective_rkmax: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    hydrogen_note: str = ""


def optimize_rkmax(atoms, rmt_values, precision=Precision.MEDIUM, strict_faq=False):
    result = RKMAXOptimizationResult()

    r_min = min(rmt_values)
    idx = rmt_values.index(r_min)
    el = atoms[idx].element

    result.r_min = r_min
    result.min_element = el
    result.base_rkmax = RKMAX_TABLE.get(el, 7.0)

    if strict_faq:
        result.rkmax = result.base_rkmax
        result.warnings.append(
            f"Strict-FAQ mode: RKMAX={result.base_rkmax} (table value, no offset)."
        )
    elif el == "H" and r_min < 0.7:
        offset = RKMAX_OFFSET_H_SMALL.get(precision, 0.0)
        result.rkmax = round(result.base_rkmax + offset, 1)
        result.hydrogen_note = (
            f"H with small RMT ({r_min:.2f}): RKMAX={result.rkmax}. "
            f"GMAX should be 18-24."
        )
    elif el == "Li" and r_min < 1.0:
        offset = RKMAX_OFFSET_H_SMALL.get(precision, 0.0)
        result.rkmax = round(result.base_rkmax + offset, 1)
    else:
        offset = RKMAX_OFFSET.get(precision, 0.0)
        result.rkmax = round(result.base_rkmax + offset, 1)

    eff = []
    for i, a in enumerate(atoms):
        ek = result.rkmax * rmt_values[i] / max(r_min, 0.01)
        eff.append(round(ek, 2))
        tv = RKMAX_TABLE.get(a.element, 7.0)
        if ek > 12:
            result.warnings.append(
                f"{a.element}({i}) RKMAX_eff={ek:.1f} > 12 — "
                f"consider adjusting RMTs."
            )
        elif ek < tv - 1.5:
            result.warnings.append(
                f"{a.element}({i}) RKMAX_eff={ek:.1f} below "
                f"recommended {tv} — consider increasing RKMAX."
            )

    result.effective_rkmax = eff
    return result
