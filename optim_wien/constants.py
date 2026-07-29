"""
WIEN2k Parameter Optimization — Constants & Reference Data

References:
  - P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)
  - WIEN2k User's Guide
  - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
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

ELECTRONEGATIVITY = {
    "H": 2.20, "He": 0, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55,
    "N": 3.04, "O": 3.44, "F": 3.98, "Ne": 0, "Na": 0.93, "Mg": 1.31,
    "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": 0,
    "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66,
    "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
    "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96, "Kr": 3.00,
    "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16,
    "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69,
    "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.10, "I": 2.66, "Xe": 2.60,
    "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13,
    "Nd": 1.14, "Pm": 1.13, "Sm": 1.17, "Eu": 1.20, "Gd": 1.20,
    "Tb": 1.10, "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25,
    "Yb": 1.10, "Lu": 1.27, "Hf": 1.30, "Ta": 1.50, "W": 2.36,
    "Re": 1.90, "Os": 2.20, "Ir": 2.20, "Pt": 2.28, "Au": 2.54,
    "Hg": 2.00, "Tl": 1.62, "Pb": 2.33, "Bi": 2.02, "Po": 2.00,
    "At": 2.20, "Rn": 2.20, "Fr": 0.70, "Ra": 0.90, "Ac": 1.10,
    "Th": 1.30, "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28,
}

INITIAL_RMT = {
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
    "H": 3.0, "He": 3.5, "Li": 4.5, "Be": 5.0, "B": 5.0, "Si": 5.0, "C": 5.5, "P": 5.5,
    "N": 6.5, "S": 6.5,
    "O": 7.0, "F": 7.0, "Cl": 7.0,
    "Na": 6.5, "K": 6.5, "Rb": 6.5, "Cs": 6.5, "Fr": 6.5,
    "Mg": 6.5, "Ca": 6.5, "Sr": 6.5, "Ba": 6.5, "Ra": 6.5,
    "Al": 6.5, "Ga": 6.5, "Ge": 6.5,
    "Ne": 5.5, "Ar": 5.5, "Kr": 6.5, "Xe": 6.5, "Rn": 7.0,
    "Sc": 7.5, "Ti": 7.5, "V": 7.5, "Cr": 7.5, "Mn": 8.0,
    "Fe": 8.0, "Co": 8.0, "Ni": 8.0, "Cu": 8.0, "Zn": 8.0,
    "As": 7.5, "Se": 7.5, "Br": 7.5,
    "Y": 7.5, "Zr": 7.5, "Nb": 7.5, "Mo": 7.5, "Tc": 8.0,
    "Ru": 8.0, "Rh": 8.0, "Pd": 8.0, "Ag": 8.0, "Cd": 8.0,
    "In": 8.0, "Sn": 8.0, "Sb": 8.0, "Te": 8.0, "I": 8.0,
    "La": 8.0, "Ce": 8.0, "Hf": 8.0, "Ta": 8.0, "W": 8.0, "Re": 8.0,
    "Os": 8.5, "Ir": 8.5, "Pt": 8.5, "Au": 8.5, "Hg": 8.5,
    "Tl": 8.5, "Pb": 8.5, "Bi": 8.5, "Po": 8.5, "At": 8.5,
    "Pr": 8.5, "Nd": 8.5, "Pm": 8.5, "Sm": 8.5, "Eu": 8.5,
    "Gd": 8.5, "Tb": 8.5, "Dy": 8.5, "Ho": 8.5, "Er": 8.5,
    "Tm": 8.5, "Yb": 8.5, "Lu": 8.5,
    "Ac": 8.5, "Th": 8.5, "Pa": 8.5, "U": 8.5, "Np": 8.5,
    "Pu": 8.5, "Am": 8.5, "Cm": 8.5, "Bk": 8.5, "Cf": 8.5,
    "Es": 8.5, "Fm": 8.5, "Md": 8.5, "No": 8.5, "Lr": 8.5,
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

REFERENCES = [
    "P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)",
    "WIEN2k User's Guide, http://susi.theochem.tuwien.ac.at/reg_user/textbooks/usersguide.pdf",
    "WIEN2k FAQ RMT: http://www.wien2k.at/reg_user/faq/rmt.html",
    "WIEN2k FAQ RKMAX: http://www.wien2k.at/reg_user/faq/rkmax.html",
    "WIEN2k FAQ kgen: http://www.wien2k.at/reg_user/faq/kgen.html",
    "L.D. Marks, Optimization Notes, WIEN2k Workshop 2006",
    "Blaha Lectures, WIEN2k Workshop 2015",
]
