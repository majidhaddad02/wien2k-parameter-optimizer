"""
LMAX / LVNS Optimization — Angular Momentum Cutoffs.

LMAX: maximum l for wavefunction expansion inside spheres.
LVNS: maximum l for non-spherical matrix elements (VNS).
"""

from dataclasses import dataclass, field

from .constants import (
    ELEMENT_TYPE, LMAX_SP, LMAX_D, LMAX_F,
    LVNS_SMALL, LVNS_MEDIUM, LVNS_LARGE, RMT_LARGE_THRESHOLD,
)


@dataclass
class LMAXResult:
    lmax_values: list = field(default_factory=list)
    lvns_global: int = 4
    hdlo_recommendations: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def optimize_lmax(atoms, rmt_values):
    result = LMAXResult()
    lmax_list = []

    for i, a in enumerate(atoms):
        et = ELEMENT_TYPE.get(a.element, "sp")
        rmt = rmt_values[i]

        if et == "f":
            lmax = LMAX_F
        elif et == "d":
            lmax = LMAX_D
        else:
            lmax = LMAX_SP

        if rmt > 2.5:
            result.hdlo_recommendations.append(
                f"{a.element}({i}): RMT={rmt:.2f} > 2.5. HDLOs required. LVNS≥6."
            )
        elif rmt > RMT_LARGE_THRESHOLD and et in ("d", "f"):
            result.hdlo_recommendations.append(
                f"{a.element}({i}): RMT={rmt:.2f} — HDLOs recommended for {et}-states."
            )

        lmax_list.append(lmax)

    result.lmax_values = lmax_list

    max_rmt = max(rmt_values) if rmt_values else 2.0
    has_f = any(ELEMENT_TYPE.get(a.element, "sp") == "f" for a in atoms)
    has_d = any(ELEMENT_TYPE.get(a.element, "sp") == "d" for a in atoms)

    if has_f or max_rmt > 2.5:
        result.lvns_global = LVNS_LARGE
        result.notes.append("LVNS = 8 (f-elements or large RMTs > 2.5).")
    elif has_d or max_rmt > RMT_LARGE_THRESHOLD:
        result.lvns_global = LVNS_MEDIUM
        result.notes.append("LVNS = 6 (d-elements or large RMTs > 2.2).")
    else:
        result.lvns_global = LVNS_SMALL
        result.notes.append("LVNS = 4 (sp-elements).")

    return result
