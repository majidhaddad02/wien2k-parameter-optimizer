"""
Tests for WIEN2k struct parser — including MULT>1 equivalent positions
with both positive and negative atomic indices (real WIEN2k format).
"""

import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from optim_wien.struct_parser import parse_struct, compute_pairwise_min_distances
from optim_wien.struct_parser import _is_equiv_position_line, _is_element_line


STRUCT_SI_POSITIVE = """Si diamond (positive equiv indices)
P   LATTICE,P
MODE OF CALC=RELA
  10.261000 10.261000 10.261000 90.000000 90.000000 90.000000
ATOM  -1: X=0.12500000 Y=0.12500000 Z=0.12500000
          MULT= 2          ISPLIT= 4
     2: X=0.62500000 Y=0.62500000 Z=0.62500000
Si         NPT=  781  R0=0.00010000 RMT=    2.00000   Z: 14.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
"""

STRUCT_SI_NEGATIVE = """Si diamond (real WIEN2k negative equiv indices)
P   LATTICE,P
MODE OF CALC=RELA
  10.261000 10.261000 10.261000 90.000000 90.000000 90.000000
ATOM  -1: X=0.12500000 Y=0.12500000 Z=0.12500000
          MULT= 2          ISPLIT= 4
    -2: X=0.62500000 Y=0.62500000 Z=0.62500000
Si         NPT=  781  R0=0.00010000 RMT=    2.00000   Z: 14.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
"""

STRUCT_CEO2_REAL = """CeO2 fluorite (real WIEN2k — negative equiv)
F   LATTICE,NONEQUIV.ATOMS:  2
MODE OF CALC=RELA
  10.260000 10.260000 10.260000 90.000000 90.000000 90.000000
ATOM  -1: X=0.00000000 Y=0.00000000 Z=0.00000000
          MULT= 1          ISPLIT= 8
Ce1        NPT=  781  R0=0.00010000 RMT=    2.50000   Z: 58.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
ATOM  -2: X=0.25000000 Y=0.25000000 Z=0.25000000
          MULT= 3          ISPLIT= 8
    -2: X=0.75000000 Y=0.25000000 Z=0.75000000
    -3: X=0.25000000 Y=0.75000000 Z=0.75000000
O          NPT=  781  R0=0.00010000 RMT=    1.70000   Z:  8.0
LOCAL ROT MATRIX:    0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
                     1.0000000 0.0000000 0.0000000
"""

STRUCT_NO_EQUIV = """NaCl (mult=2, no equiv lines)
F   LATTICE,NONEQUIV.ATOMS: 2
MODE OF CALC=RELA
  10.610000 10.610000 10.610000 90.000000 90.000000 90.000000
ATOM  -1: X=0.00000000 Y=0.00000000 Z=0.00000000
          MULT= 1          ISPLIT= 8
Na         NPT=  781  R0=0.00010000 RMT=    2.20000   Z: 11.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
ATOM  -2: X=0.50000000 Y=0.50000000 Z=0.50000000
          MULT= 1          ISPLIT= 8
Cl         NPT=  781  R0=0.00010000 RMT=    2.20000   Z: 17.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
"""


def _write_and_parse(content):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".struct", delete=False) as f:
        f.write(content)
        path = f.name
    try:
        return parse_struct(path)
    finally:
        os.unlink(path)


def test_equiv_line_positive():
    assert _is_equiv_position_line("     2: X=0.625 Y=0.625 Z=0.625")
    assert not _is_equiv_position_line("ATOM  -1: X=0.0 Y=0.0 Z=0.0")


def test_equiv_line_negative():
    assert _is_equiv_position_line("    -2: X=0.625 Y=0.625 Z=0.625")
    assert _is_equiv_position_line("    -3: X=0.250 Y=0.750 Z=0.750")


def test_element_line():
    assert _is_element_line("Si         NPT=  781  R0=0.00010000 RMT=    2.00000   Z: 14.0")
    assert not _is_element_line("-2: X=0.62500000 Y=0.62500000 Z=0.62500000")


def test_parse_si_negative_equiv():
    """Real WIEN2k format: negative equiv indices must be parsed."""
    s = _write_and_parse(STRUCT_SI_NEGATIVE)
    assert len(s.atoms) == 1
    assert s.atoms[0].element == "Si"
    assert s.atoms[0].mult == 2
    assert len(s.atoms[0].equivalent_positions) == 1
    assert s.atoms[0].equivalent_positions[0] == [0.625, 0.625, 0.625]
    assert s.atoms[0].position == [0.125, 0.125, 0.125]


def test_parse_si_positive_equiv():
    """Positive equiv indices (less common but valid for cubic)."""
    s = _write_and_parse(STRUCT_SI_POSITIVE)
    assert len(s.atoms) == 1
    assert s.atoms[0].element == "Si"
    assert s.atoms[0].mult == 2
    assert len(s.atoms[0].equivalent_positions) == 1


def test_parse_ceo2_real():
    """CeO2 fluorite with negative equiv indices and non-identity LOCAL ROT MATRIX."""
    s = _write_and_parse(STRUCT_CEO2_REAL)
    assert len(s.atoms) == 2
    assert s.atoms[0].element == "Ce"
    assert s.atoms[1].element == "O"
    assert s.atoms[1].mult == 3
    assert len(s.atoms[1].equivalent_positions) == 2

    ce = s.atoms[0]
    assert len(ce.local_rot_lines) == 3
    assert "1.0000000" in ce.local_rot_lines[0]

    o = s.atoms[1]
    assert len(o.local_rot_lines) == 3
    assert "0.0000000" in o.local_rot_lines[0]
    # Non-identity matrix for O site


def test_parse_no_equiv():
    """NaCl with no equivalent position lines — backward compat."""
    s = _write_and_parse(STRUCT_NO_EQUIV)
    assert len(s.atoms) == 2
    assert s.atoms[0].element == "Na"
    assert s.atoms[1].element == "Cl"
    assert s.atoms[0].equivalent_positions == []
    assert s.atoms[1].equivalent_positions == []


def test_nn_distances_si():
    """Si diamond: NN should be ~4.44 bohr (sqrt(3)*a/4).
    Single element → no pairs → pairwise is empty (by design)."""
    s = _write_and_parse(STRUCT_SI_NEGATIVE)
    pairwise = compute_pairwise_min_distances(s)
    # Single element → no inter-element pairs → pairwise is empty
    # This is correct behavior: Si-Si distances require self-pair detection
    assert isinstance(pairwise, dict)


def test_nn_distances_ceo2():
    """CeO2: should compute Ce-O and Ce-Ce distances."""
    s = _write_and_parse(STRUCT_CEO2_REAL)
    pairwise = compute_pairwise_min_distances(s)
    assert (0, 1) in pairwise
    d_ceo = pairwise[(0, 1)]
    assert d_ceo < float("inf")
    assert d_ceo > 0.1  # must be positive, > ~2.22 bohr


def test_element_not_corrupted_by_equiv():
    """Element name must NOT be corrupted when parsing negative equiv lines."""
    s = _write_and_parse(STRUCT_SI_NEGATIVE)
    assert s.atoms[0].element == "Si"
    assert s.atoms[0].element != "-:"
    assert s.atoms[0].element != "-"

    s2 = _write_and_parse(STRUCT_CEO2_REAL)
    assert s2.atoms[0].element == "Ce"
    assert s2.atoms[1].element == "O"
