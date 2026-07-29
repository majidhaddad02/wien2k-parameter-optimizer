"""
Tests for Aitken Δ² extrapolation and exponential fit.
"""

import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from optim_wien.converge import aitken_extrapolate, exponential_fit_linearized


def test_aitken_synthetic_exact():
    """Exact exponential: E = -100 + 2*exp(-0.5*x) at x=0,2,4."""
    x = [0.0, 2.0, 4.0]
    E = [-100 + 2.0 * math.exp(-0.5 * xi) for xi in x]

    result = aitken_extrapolate(x, E)

    assert result.is_valid
    assert abs(result.e_inf - (-100.0)) < 0.01
    assert abs(result.alpha - 0.5) < 0.05
    assert abs(result.amplitude - 2.0) < 0.5


def test_aitken_monotonic_check():
    """Non-monotonically varying energies → extrapolation aborted."""
    result = aitken_extrapolate(
        [0.0, 2.0, 4.0],
        [-100.0, -99.9, -99.95],  # E1 > E2 > E3? No: -100 < -99.9 < -99.95
    )
    assert not result.is_valid
    assert not result.is_monotonic


def test_aitken_flat_converged():
    """When energies barely change across three points, denominator ≈ 0."""
    result = aitken_extrapolate(
        [0.0, 2.0, 4.0],
        [-100.0, -100.00001, -100.00002],
    )
    # Nearly equal energies → denominator ≈ 0 → not valid
    assert not result.is_valid


def test_aitken_denominator_near_zero():
    """Linear spacing → denominator ~0 → fit fails gracefully."""
    result = aitken_extrapolate(
        [0.0, 2.0, 4.0],
        [-100.0, -99.5, -99.0],  # linear ΔE=0.5 → de1==de2 → denominator=0
    )
    assert not result.is_valid


def test_aitken_predict_x():
    """Prediction using AitkenResult.predict_x()."""
    result = aitken_extrapolate(
        [0.0, 2.0, 4.0],
        [-100 + 2.0 * math.exp(-0.5 * xi) for xi in [0.0, 2.0, 4.0]],
    )
    assert result.is_valid

    # Need |A|*exp(-alpha*x) < etol → x > -ln(etol/|A|)/alpha
    x_req = result.predict_x(etol_ry=1.0, step=2.0)
    assert x_req is not None
    assert x_req > 0

    x_tight = result.predict_x(etol_ry=0.1, step=2.0)
    assert x_tight is not None
    assert x_tight >= x_req


def test_aitken_physics():
    """Verify decay ratio r ∈ (0,1) for proper exponential convergence."""
    x = [0.0, 1.0, 2.0]
    E = [-100 + 10 * math.exp(-1.0 * xi) for xi in x]
    result = aitken_extrapolate(x, E)
    assert result.is_valid
    # r = (E3-E2)/(E2-E1): both differences are negative and decreasing
    # so 0 < r < 1 for proper exponential convergence


def test_exponential_fit_linearized():
    """4+ point exponential fit recovers known parameters."""
    x = [0.0, 2.0, 4.0, 6.0, 8.0]
    E = [-100 + 2.0 * math.exp(-0.5 * xi) for xi in x]

    e_inf, alpha, amp = exponential_fit_linearized(x, E)

    assert e_inf is not None
    assert alpha is not None
    assert abs(e_inf - (-100.0)) < 0.5
    assert abs(alpha - 0.5) < 0.2
    assert abs(amp - 2.0) < 1.0


def test_exponential_fit_few_points():
    """Fewer than 4 points → returns None."""
    x = [0.0, 2.0, 4.0]
    E = [-100 + 2.0 * math.exp(-0.5 * xi) for xi in x]
    e_inf, alpha, amp = exponential_fit_linearized(x, E)
    assert e_inf is None


def test_exponential_fit_non_monotonic():
    """Non-monotonic data → returns None."""
    x = [0.0, 2.0, 4.0, 6.0]
    E = [-100.0, -99.9, -100.1, -100.2]  # goes up between 2 and 4
    e_inf, alpha, amp = exponential_fit_linearized(x, E)
    assert e_inf is None
