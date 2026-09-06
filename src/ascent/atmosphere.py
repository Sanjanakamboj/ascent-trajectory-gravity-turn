"""Simplified atmosphere model for M2 verification.

This is a **single-exponential** density model, not the full layered US Standard
Atmosphere 1976 that DESIGN.md S1.6 originally flagged as "planned." That is a deliberate
M2 scope decision (documented in DESIGN.md's M2 section as a discrepancy from the M1
note): a single exponential is deterministic, trivial to verify analytically, and
sufficient for M2's goal of verifying the dynamics/propulsion/drag plumbing. It is
explicitly **not** a high-fidelity launch atmosphere and must not be read as one.

Model
-----
    rho(h) = RHO0 * exp(-h / SCALE_HEIGHT)   for 0 <= h <= H_MAX
    rho(h) = RHO0                            for h < 0   (clamped to sea level)
    rho(h) = 0.0                             for h > H_MAX  (vacuum approximation)

Parameters
----------
RHO0 : float
    Reference (sea-level) density, 1.225 kg/m^3 (standard sea-level air density).
SCALE_HEIGHT : float
    8500 m, a representative single scale height for Earth's lower/middle atmosphere.
H_MAX : float
    100,000 m. Above this altitude density is treated as exactly zero (vacuum). This
    is an explicit, documented cutoff, not a physical claim that the atmosphere ends
    there — it exists so drag deterministically vanishes once altitude is well beyond
    the region this simplified model is meant to represent.

Treatment of h < 0
-------------------
Altitudes below zero (e.g. small integrator overshoot near ground level) are clamped to
h = 0 rather than extrapolated, so density never exceeds the sea-level reference value.

Treatment of h > H_MAX
-----------------------
Density is exactly zero; drag consequently vanishes exactly there (see
``test_atmosphere.py::test_high_altitude_limit``).
"""

import numpy as np

RHO0 = 1.225            # kg/m^3
SCALE_HEIGHT = 8500.0   # m
H_MAX = 100_000.0       # m


def density(h):
    """Return atmospheric density rho(h) [kg/m^3] for altitude(s) h [m].

    Accepts a scalar or array-like ``h``. Returns the same shape (scalar in, scalar
    out; array in, ``numpy.ndarray`` out).
    """
    h_arr = np.asarray(h, dtype=float)
    h_clamped = np.clip(h_arr, 0.0, H_MAX)
    rho = RHO0 * np.exp(-h_clamped / SCALE_HEIGHT)
    rho = np.where(h_arr > H_MAX, 0.0, rho)
    if np.isscalar(h) or (hasattr(h, "shape") and h_arr.shape == ()):
        return float(rho)
    return rho
