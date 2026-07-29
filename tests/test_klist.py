"""
Tests for klist generation — END must be last line.
"""

import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from optim_wien.input_generator import _gen_klist


class FakeKMesh:
    def __init__(self, mesh=(4, 4, 4), shift=(0.5, 0.5, 0.5)):
        self.mesh = mesh
        self.shift = shift
        self.total_points = mesh[0] * mesh[1] * mesh[2]
        self.gamma_centered = False
        self.add_inversion = True


def test_end_is_last_line():
    """END must appear at the very end of klist, after all k-point coordinates."""
    km = FakeKMesh(mesh=(4, 4, 4))
    klist = _gen_klist("test", km)
    lines = klist.strip().split("\n")

    # Total lines: header (4) + k-points (64) + END (1) = 69
    assert len(lines) == 69, f"Expected 69 lines, got {len(lines)}"

    # END must be the last line
    assert lines[-1] == "END", f"Last line should be END, got '{lines[-1]}'"

    # END must NOT be in the first 5 lines
    header_lines = lines[:5]
    assert "END" not in header_lines, f"END found early in klist: {header_lines}"

    # Verify k-point lines exist and have correct format
    first_k = lines[4]  # line 5 (0-indexed), first k-point
    parts = first_k.split()
    assert len(parts) == 4, f"Expected 4 values per k-point, got {len(parts)}: {first_k}"
    assert float(parts[3]) > 0  # weight should be positive


def test_single_kpoint():
    """Gamma-only mesh (1x1x1) still puts END last."""
    km = FakeKMesh(mesh=(1, 1, 1), shift=(0.0, 0.0, 0.0))
    klist = _gen_klist("test", km)
    lines = klist.strip().split("\n")
    assert lines[-1] == "END"
    assert len(lines) == 6  # 4 header + 1 k-point + END + blank


def test_large_mesh():
    """12x12x12 mesh — END must still be last line."""
    km = FakeKMesh(mesh=(12, 12, 12))
    klist = _gen_klist("test", km)
    lines = klist.strip().split("\n")
    assert lines[-1] == "END"
    # 4 header + 1728 k-points + END = 1733 total
    assert len(lines) == 1733


def test_shifted_mesh():
    """Non-Gamma mesh with shift (0.5, 0.5, 0.5)."""
    km = FakeKMesh(mesh=(8, 8, 8), shift=(0.5, 0.5, 0.5))
    klist = _gen_klist("test", km)
    lines = klist.strip().split("\n")
    assert lines[-1] == "END"
    assert "0.50" in lines[3]  # shift line


def test_no_inversion():
    """No inversion gives 'no INV' in header."""
    km = FakeKMesh(mesh=(4, 4, 4))
    km.add_inversion = False
    klist = _gen_klist("test", km)
    assert "no INV" in klist


def test_hex_mesh():
    """Hexagonal mesh 12x12x8 — END still last."""
    km = FakeKMesh(mesh=(12, 12, 8))
    klist = _gen_klist("test", km)
    lines = klist.strip().split("\n")
    assert lines[-1] == "END"
