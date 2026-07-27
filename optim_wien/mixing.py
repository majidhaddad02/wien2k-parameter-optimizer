"""
Mixing & TEMP Optimization — SCF Convergence Parameters.

Based on WIEN2k User's Guide, Section 4.3 (case.in2).
"""

from dataclasses import dataclass, field

from .constants import (
    ELEMENT_TYPE,
    MIXING_PRATT_INSULATOR, MIXING_PRATT_SEMICONDUCTOR,
    MIXING_MSR1A_METAL, MIXING_MSEC1_MAGNETIC,
    TEMP_INSULATOR, TEMP_SEMICONDUCTOR, TEMP_METAL, TEMP_MAGNETIC,
    ECONV_DEFAULT, ECONV_FORCES, CCONV_DEFAULT,
    MAX_SCF_CYCLES_DEFAULT, MAX_SCF_CYCLES_METAL,
    CalcType, Precision,
)


@dataclass
class MixingResult:
    scheme: str = "PRATT"
    mixing_factor: float = 0.30
    temp: float = 0.001
    econv: float = 0.0001
    cconv: float = 0.001
    max_scf: int = 40
    tetra: bool = False
    notes: list = field(default_factory=list)


def optimize_mixing(
    structure,
    system_type="semiconductor",
    calc_type=CalcType.SCF,
    precision=Precision.MEDIUM,
    magnetic=False,
):
    result = MixingResult()
    atoms = structure.atoms
    has_tm = any(ELEMENT_TYPE.get(a.element, "sp") in ("d", "f") for a in atoms)
    is_metal = "metal" in system_type

    if is_metal:
        if magnetic or "large" in system_type:
            result.scheme = "MSEC1"
            result.mixing_factor = MIXING_MSEC1_MAGNETIC
            if magnetic:
                result.notes.append("Metallic + magnetic → MSEC1 mixing.")
            else:
                result.notes.append("Large metallic system → MSEC1 mixing.")
        else:
            result.scheme = "MSR1a"
            result.mixing_factor = MIXING_MSR1A_METAL
            result.notes.append("Metallic system → MSR1a mixing.")
        result.temp = TEMP_MAGNETIC if (magnetic or "large" in system_type) else TEMP_METAL
        result.tetra = True
        result.max_scf = MAX_SCF_CYCLES_METAL
    else:
        result.scheme = "PRATT"
        if system_type == "insulator" or system_type == "insulator_large":
            result.mixing_factor = MIXING_PRATT_INSULATOR
            result.temp = TEMP_INSULATOR
        else:
            result.mixing_factor = MIXING_PRATT_SEMICONDUCTOR
            result.temp = TEMP_SEMICONDUCTOR
        result.notes.append("Semiconductor/insulator → PRATT mixing.")
        result.tetra = False
        result.max_scf = MAX_SCF_CYCLES_DEFAULT

    if calc_type in (CalcType.FORCES, CalcType.EFG):
        result.econv = ECONV_FORCES
        result.cconv = CCONV_DEFAULT / 10
        result.notes.append("Forces/EFG → stricter convergence.")
    elif precision in (Precision.HIGH, Precision.VERY_HIGH):
        result.econv = ECONV_FORCES
        result.cconv = CCONV_DEFAULT
        result.notes.append("High precision → stricter energy convergence.")
    else:
        result.econv = ECONV_DEFAULT
        result.cconv = CCONV_DEFAULT

    return result
