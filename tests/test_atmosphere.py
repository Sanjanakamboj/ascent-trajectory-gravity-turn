import math

import numpy as np
import pytest

from ascent import atmosphere as atmo


def test_sea_level_reference_value():
    assert atmo.density(0.0) == pytest.approx(atmo.RHO0)


def test_positivity_within_range():
    for h in [0, 1000, 10_000, 50_000, atmo.H_MAX]:
        assert atmo.density(h) > 0.0


def test_monotonic_decrease():
    hs = np.linspace(0, atmo.H_MAX, 50)
    rhos = [atmo.density(h) for h in hs]
    assert all(rhos[i] > rhos[i + 1] for i in range(len(rhos) - 1))


def test_exact_exponential_scaling():
    h = 8500.0  # one scale height
    assert atmo.density(h) == pytest.approx(atmo.RHO0 * math.exp(-1.0))
    h2 = 17000.0  # two scale heights
    assert atmo.density(h2) == pytest.approx(atmo.RHO0 * math.exp(-2.0))


def test_high_altitude_limit_is_exact_vacuum():
    assert atmo.density(atmo.H_MAX + 1.0) == 0.0
    assert atmo.density(1e7) == 0.0


def test_negative_altitude_clamped_to_sea_level():
    assert atmo.density(-100.0) == pytest.approx(atmo.RHO0)
    assert atmo.density(-1e6) == pytest.approx(atmo.RHO0)


def test_vectorized_input():
    hs = np.array([0.0, 8500.0, atmo.H_MAX + 1.0])
    rhos = atmo.density(hs)
    assert isinstance(rhos, np.ndarray)
    np.testing.assert_allclose(rhos, [atmo.RHO0, atmo.RHO0 * math.exp(-1.0), 0.0])
