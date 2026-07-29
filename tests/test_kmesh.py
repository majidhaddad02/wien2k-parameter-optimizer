"""
Tests for k-mesh volume correction for centered lattices.
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from optim_wien.kmesh import optimize_kmesh, _primitive_volume_factor


def test_primitive_volume_factor_fcc():
    assert _primitive_volume_factor("F") == 0.25
    assert _primitive_volume_factor("FCC") == 0.25
    assert _primitive_volume_factor("CF") == 0.25


def test_primitive_volume_factor_bcc():
    assert _primitive_volume_factor("B") == 0.50
    assert _primitive_volume_factor("BCC") == 0.50
    assert _primitive_volume_factor("BC") == 0.50


def test_primitive_volume_factor_base():
    assert _primitive_volume_factor("CXY") == 0.50
    assert _primitive_volume_factor("CYZ") == 0.50
    assert _primitive_volume_factor("CXZ") == 0.50
    assert _primitive_volume_factor("C") == 0.50


def test_primitive_volume_factor_simple():
    assert _primitive_volume_factor("P") == 1.0
    assert _primitive_volume_factor("H") == 1.0
    assert _primitive_volume_factor("R") == 1.0 / 3.0
    assert _primitive_volume_factor("UNKNOWN") == 1.0
    assert _primitive_volume_factor("") == 1.0


class FakeStructure:
    def __init__(self, lt="P", vol=1000.0, ga=90.0):
        self.lattice_type = lt
        self._volume = vol
        self.gamma = ga
        self.a = (vol) ** (1.0 / 3.0)
        self.b = self.a
        self.c = self.a

    @property
    def volume(self):
        return self._volume

    @property
    def num_atoms_primitive(self):
        return 2

    @property
    def reciprocal_lengths(self):
        import math
        a_star = 2 * math.pi / self.a
        return (a_star, a_star, a_star)

    @property
    def lattice_vectors(self):
        a = self.a
        return ([a, 0, 0], [0, a, 0], [0, 0, a])

    @property
    def atoms(self):
        class A:
            element = "Si"
        return [A()]


def test_kmesh_fcc_uses_correct_volume():
    """FCC with conventional volume V → primitive V/4 → 4× larger BZ."""
    s_p = FakeStructure(lt="P", vol=1000.0)
    s_f = FakeStructure(lt="F", vol=1000.0)

    km_p = optimize_kmesh(s_p, refinement="medium", system_type="semiconductor")
    km_f = optimize_kmesh(s_f, refinement="medium", system_type="semiconductor")

    # FCC should have DENSER mesh than P for same conventional volume
    # because primitive cell is 4× smaller → BZ is 4× larger
    p_total = km_p.total_points
    f_total = km_f.total_points

    assert f_total >= p_total, (
        f"FCC mesh ({f_total} pts) should be >= P mesh ({p_total} pts) "
        f"for same conventional volume"
    )


def test_kmesh_bcc_uses_correct_volume():
    """BCC with conventional volume V → primitive V/2 → 2× larger BZ."""
    s_p = FakeStructure(lt="P", vol=1000.0)
    s_b = FakeStructure(lt="B", vol=1000.0)

    km_p = optimize_kmesh(s_p, refinement="medium", system_type="semiconductor")
    km_b = optimize_kmesh(s_b, refinement="medium", system_type="semiconductor")

    b_total = km_b.total_points
    p_total = km_p.total_points

    assert b_total >= p_total, (
        f"BCC mesh ({b_total} pts) should be >= P mesh ({p_total} pts)"
    )
