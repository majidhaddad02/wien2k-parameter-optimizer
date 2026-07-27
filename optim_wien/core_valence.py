"""
Core/Valence Separation — Ecut Optimization.

Controls which states are treated as core vs. valence.
Based on WIEN2k User's Guide, Section 4.4 (lstart).
"""

from dataclasses import dataclass, field

from .constants import (
    ELEMENT_TYPE, ECUT_DEFAULT, ECUT_REDUCED, ECUT_PRECISE,
    ECUT_PRECISION, Precision,
)


@dataclass
class CoreValenceResult:
    ecut: float = -6.0
    use_hdlo: bool = False
    notes: list = field(default_factory=list)


def optimize_core_valence(atoms, rmt_values, precision=Precision.MEDIUM):
    result = CoreValenceResult()

    ecut = ECUT_PRECISION.get(precision, ECUT_DEFAULT)
    result.ecut = ecut

    use_hdlo = False
    for i, a in enumerate(atoms):
        rmt = rmt_values[i]
        et = ELEMENT_TYPE.get(a.element, "sp")
        if rmt > 2.2 and et in ("d", "f"):
            use_hdlo = True
        if rmt > 2.5:
            use_hdlo = True

    result.use_hdlo = use_hdlo

    if use_hdlo:
        result.notes.append("HDLOs activated for large-RMT atoms.")
    result.notes.append(f"Ecut = {ecut:.1f} Ry (precision={precision.value}).")

    return result
