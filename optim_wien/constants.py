"""
WIEN2k Parameter Optimization — Constants & Reference Data

References:
  - P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)
  - WIEN2k User's Guide
  - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/

Provenance notes:
  - INITIAL_RMT: heuristic defaults based on common WIEN2k workshop
    practice. These are NOT from any published table. A professional
    user should use RMT values from their own case.struct file, which
    takes priority over these defaults.
  - RKMAX_TABLE: only O=7, N=6.5, C=5.5, Fe=8, Cu=8 are from Blaha
    2020 Table I. All other elements default to 7.0 (WIEN2k default).
    See rkmax.py for the fallback.
"""

from enum import Enum

__version__ = "1.0.0"


class CalcType(Enum):
    SCF = "scf"
    RELAXATION = "relaxation"
    OPTIMIZATION = "optimization"
    EOS = "eos"
    DOS = "dos"
    OPTICS = "optics"
    FORCES = "forces"
    EFG = "efg"


class Precision(Enum):
    SCREENING = "screening"
    COARSE = "coarse"
    MEDIUM = "medium"
    HIGH = "high"
    VERY_HIGH = "very_high"


ELEMENT_TYPE = {
    "H": "sp", "He": "sp", "Li": "sp", "Be": "sp", "B": "sp", "C": "sp",
    "N": "sp", "O": "sp", "F": "sp", "Ne": "sp", "Na": "sp", "Mg": "sp",
    "Al": "sp", "Si": "sp", "P": "sp", "S": "sp", "Cl": "sp", "Ar": "sp",
    "K": "sp", "Ca": "sp", "Ga": "sp", "Ge": "sp", "As": "sp", "Se": "sp",
    "Br": "sp", "Kr": "sp", "Rb": "sp", "Sr": "sp", "In": "sp", "Sn": "sp",
    "Sb": "sp", "Te": "sp", "I": "sp", "Xe": "sp", "Cs": "sp", "Ba": "sp",
    "Tl": "sp", "Pb": "sp", "Bi": "sp", "Po": "sp", "At": "sp", "Rn": "sp",
    "Fr": "sp", "Ra": "sp",
    "Sc": "d", "Ti": "d", "V": "d", "Cr": "d", "Mn": "d", "Fe": "d",
    "Co": "d", "Ni": "d", "Cu": "d", "Zn": "d", "Y": "d", "Zr": "d",
    "Nb": "d", "Mo": "d", "Tc": "d", "Ru": "d", "Rh": "d", "Pd": "d",
    "Ag": "d", "Cd": "d", "Hf": "d", "Ta": "d", "W": "d", "Re": "d",
    "Os": "d", "Ir": "d", "Pt": "d", "Au": "d", "Hg": "d",
    "La": "f", "Ce": "f", "Pr": "f", "Nd": "f", "Pm": "f", "Sm": "f",
    "Eu": "f", "Gd": "f", "Tb": "f", "Dy": "f", "Ho": "f", "Er": "f",
    "Tm": "f", "Yb": "f", "Lu": "f", "Ac": "f", "Th": "f", "Pa": "f",
    "U": "f", "Np": "f", "Pu": "f", "Am": "f", "Cm": "f", "Bk": "f",
    "Cf": "f", "Es": "f", "Fm": "f", "Md": "f", "No": "f", "Lr": "f",
}


INITIAL_RMT = {
    # Heuristic defaults — NOT from a published table.
    # A professional user should use the RMT from their own
    # case.struct file. These are only used when struct RMT < 1.0.
    "H": 0.80, "He": 1.00, "Li": 1.60, "Be": 1.40, "B": 1.30, "C": 1.30,
    "N": 1.30, "O": 1.40, "F": 1.30, "Ne": 1.20, "Na": 1.90, "Mg": 1.90,
    "Al": 2.00, "Si": 1.70, "P": 1.70, "S": 1.70, "Cl": 1.70, "Ar": 1.60,
    "K": 2.20, "Ca": 2.20, "Sc": 2.10, "Ti": 2.10, "V": 2.10, "Cr": 2.10,
    "Mn": 2.10, "Fe": 2.10, "Co": 2.10, "Ni": 2.10, "Cu": 2.10, "Zn": 2.10,
    "Ga": 1.90, "Ge": 1.80, "As": 1.80, "Se": 1.80, "Br": 1.80, "Kr": 1.70,
    "Rb": 1.90, "Sr": 1.90, "Y": 2.30, "Zr": 2.30, "Nb": 2.30, "Mo": 2.30,
    "Tc": 2.30, "Ru": 2.30, "Rh": 2.30, "Pd": 2.30, "Ag": 2.30, "Cd": 2.20,
    "In": 2.00, "Sn": 2.00, "Sb": 2.00, "Te": 2.00, "I": 2.00, "Xe": 1.90,
    "Cs": 2.00, "Ba": 2.30, "La": 2.50, "Ce": 2.50, "Pr": 2.45, "Nd": 2.45,
    "Pm": 2.45, "Sm": 2.45, "Eu": 2.40, "Gd": 2.45, "Tb": 2.45, "Dy": 2.45,
    "Ho": 2.45, "Er": 2.45, "Tm": 2.45, "Yb": 2.45, "Lu": 2.45,
    "Hf": 2.45, "Ta": 2.45, "W": 2.45, "Re": 2.45, "Os": 2.45,
    "Ir": 2.45, "Pt": 2.45, "Au": 2.45, "Hg": 2.45,
    "Tl": 2.20, "Pb": 2.20, "Bi": 2.20, "Po": 2.20, "At": 2.20, "Rn": 2.10,
    "Fr": 2.10, "Ra": 2.30,
    "Ac": 2.60, "Th": 2.60, "Pa": 2.60, "U": 2.60, "Np": 2.60, "Pu": 2.60,
    "Am": 2.60, "Cm": 2.60, "Bk": 2.60, "Cf": 2.60, "Es": 2.60, "Fm": 2.60,
    "Md": 2.60, "No": 2.60, "Lr": 2.60,
}

RKMAX_TABLE = {
    # Blaha 2020 Table I — only these 5 values are published.
    # All other elements default to 7.0 (WIEN2k default) via
    # RKMAX_TABLE.get(element, 7.0) in rkmax.py.
    "C": 5.5,
    "N": 6.5,
    "O": 7.0,
    "Fe": 8.0,
    "Cu": 8.0,
}

RKMAX_OFFSET = {
    Precision.SCREENING: -1.0,
    Precision.COARSE:    -0.5,
    Precision.MEDIUM:     0.0,
    Precision.HIGH:       0.5,
    Precision.VERY_HIGH:  1.5,
}

RKMAX_OFFSET_H_SMALL = {
    Precision.SCREENING: -0.5,
    Precision.COARSE:     0.0,
    Precision.MEDIUM:     0.0,
    Precision.HIGH:       0.0,
    Precision.VERY_HIGH:  0.5,
}

SAFETY_MARGIN_SCF = 0.90
SAFETY_MARGIN_RELAX = 0.85
MAX_RMT_RATIO = 1.50
MAX_RMT_RATIO_STRICT = 1.30
RMT_REDUCTION_RELAX = 0.93
CORE_LEAKAGE_OK = 0.001
CORE_LEAKAGE_WARN = 0.002
CORE_LEAKAGE_CRITICAL = 0.01

GMAX_DEFAULT = 12.0
GMAX_HYDROGEN = 20.0
GMAX_PRECISION = {
    Precision.SCREENING: 10.0,
    Precision.COARSE:    12.0,
    Precision.MEDIUM:    14.0,
    Precision.HIGH:      16.0,
    Precision.VERY_HIGH: 20.0,
}

LMAX_SP = 6
LMAX_D  = 10
LMAX_F  = 12
LVNS_SMALL = 4
LVNS_MEDIUM = 6
LVNS_LARGE = 8
RMT_LARGE_THRESHOLD = 2.20

MIXING_PRATT_INSULATOR = 0.30
MIXING_PRATT_SEMICONDUCTOR = 0.25
MIXING_MSR1A_METAL = 0.20
MIXING_MSEC1_MAGNETIC = 0.15
TEMP_INSULATOR = 0.0001
TEMP_SEMICONDUCTOR = 0.001
TEMP_METAL = 0.001
TEMP_MAGNETIC = 0.002

ECONV_DEFAULT = 0.0001
ECONV_FORCES = 0.00001
ECONV_EFG = 0.000005
CCONV_DEFAULT = 0.001
CCONV_STRICT = 0.0001
MAX_SCF_CYCLES_DEFAULT = 40
MAX_SCF_CYCLES_METAL = 80

ECUT_DEFAULT = -6.0
ECUT_REDUCED = -4.0
ECUT_PRECISE = -8.0
ECUT_PRECISION = {
    Precision.SCREENING: -4.0,
    Precision.COARSE:    -5.0,
    Precision.MEDIUM:    -6.0,
    Precision.HIGH:       -7.0,
    Precision.VERY_HIGH:  -8.0,
}

VXCTYPE_DEFAULT = 13
VXCTYPE_PBE = 13
VXCTYPE_LDA = 5
VXCTYPE_WC = 11
VXCTYPE_PBESOL = 19
VXCTYPE_SCAN = 28
VXCTYPE_HSE = 40

NR2V_DEFAULT = 1

INIT_LAPW_PREC = {
    Precision.SCREENING: 0,
    Precision.COARSE: 0,
    Precision.MEDIUM: 1,
    Precision.HIGH: 2,
    Precision.VERY_HIGH: 3,
}

PRECISION_TO_PREC_FLAG = {
    "screening": 0, "coarse": 0, "medium": 1, "high": 2, "very_high": 3,
}

REFERENCES = [
    "P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)",
    "WIEN2k User's Guide, http://susi.theochem.tuwien.ac.at/reg_user/textbooks/usersguide.pdf",
    "WIEN2k FAQ RMT: http://www.wien2k.at/reg_user/faq/rmt.html",
    "WIEN2k FAQ RKMAX: http://www.wien2k.at/reg_user/faq/rkmax.html",
    "WIEN2k FAQ kgen: http://www.wien2k.at/reg_user/faq/kgen.html",
    "L.D. Marks, Optimization Notes, WIEN2k Workshop 2006",
    "Blaha Lectures, WIEN2k Workshop 2015",
]
