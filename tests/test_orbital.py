import math

import pytest

from ascent import orbital
from ascent.constants import MU_EARTH, R_EARTH, TARGET_ALTITUDE


def test_circular_orbit_case():
    r = R_EARTH + TARGET_ALTITUDE
    v = math.sqrt(MU_EARTH / r)
    state = orbital.compute_orbital_state(r, v, gamma=0.0)

    assert state.specific_energy == pytest.approx(-MU_EARTH / (2 * r), rel=1e-9)
    assert state.eccentricity == pytest.approx(0.0, abs=1e-9)
    assert state.is_bound
    assert not state.intersects_earth
    assert state.perigee_altitude == pytest.approx(TARGET_ALTITUDE, rel=1e-6)
    assert state.apogee_altitude == pytest.approx(TARGET_ALTITUDE, rel=1e-6)
    assert state.circular_speed_here == pytest.approx(v, rel=1e-12)


def test_elliptical_orbit_semi_major_axis_matches_vis_viva():
    # Pick an elliptical state and independently verify via vis-viva: v^2 = mu(2/r - 1/a)
    r = R_EARTH + 300_000.0
    v = 7900.0
    gamma = math.radians(5.0)
    state = orbital.compute_orbital_state(r, v, gamma)

    assert state.is_bound
    a_expected = 1.0 / (2.0 / r - v**2 / MU_EARTH)
    assert state.semi_major_axis == pytest.approx(a_expected, rel=1e-9)

    # cross-check energy/angular-momentum reconstruction from (a, e)
    a, e = state.semi_major_axis, state.eccentricity
    energy_from_a = -MU_EARTH / (2 * a)
    assert energy_from_a == pytest.approx(state.specific_energy, rel=1e-9)
    h_from_ae = math.sqrt(MU_EARTH * a * (1 - e**2))
    assert h_from_ae == pytest.approx(state.specific_angular_momentum, rel=1e-6)


def test_hyperbolic_case_has_no_apogee_and_eccentricity_above_one():
    r = R_EARTH + 300_000.0
    v = 12_000.0  # well above local escape speed
    state = orbital.compute_orbital_state(r, v, gamma=0.0)

    assert not state.is_bound
    assert state.eccentricity >= 1.0
    assert state.semi_major_axis is None
    assert state.apogee_altitude is None
    # perigee is still defined for an open trajectory
    assert state.perigee_radius > 0


def test_suborbital_trajectory_flagged_as_intersecting_earth():
    # Apogee well above target, but perigee well below the surface -- e.g. a steep,
    # slow lofted trajectory. Must NOT be reported as a valid/circular orbit.
    r = R_EARTH + 500_000.0
    v = 4000.0
    gamma = math.radians(40.0)
    state = orbital.compute_orbital_state(r, v, gamma)

    assert state.perigee_radius < R_EARTH
    assert state.intersects_earth
    assert not orbital.is_circular_orbit_achieved(state, TARGET_ALTITUDE)


def test_max_altitude_alone_does_not_imply_orbit_achieved():
    # High apogee (well above 400 km) but this instantaneous state's perigee is deep
    # inside the Earth (steep flight path angle) -- must not be classified as orbit.
    r = R_EARTH + 900_000.0
    v = 3500.0
    gamma = math.radians(70.0)
    state = orbital.compute_orbital_state(r, v, gamma)

    assert state.altitude > TARGET_ALTITUDE  # "max altitude > target" is true here...
    assert not orbital.is_circular_orbit_achieved(state, TARGET_ALTITUDE)  # ...but not orbit


def test_orbit_achieved_criterion_true_for_near_circular_state():
    r = R_EARTH + TARGET_ALTITUDE
    v = math.sqrt(MU_EARTH / r) + 2.0  # tiny overspeed, within tolerance
    state = orbital.compute_orbital_state(r, v, gamma=0.0)
    assert orbital.is_circular_orbit_achieved(state, TARGET_ALTITUDE, altitude_tolerance=10_000.0)


def test_velocity_component_round_trip():
    v, gamma = 5000.0, math.radians(23.0)
    vr, vt = orbital.velocity_components(v, gamma)
    v2, gamma2 = orbital.speed_and_gamma_from_components(vr, vt)
    assert v2 == pytest.approx(v, rel=1e-12)
    assert gamma2 == pytest.approx(gamma, rel=1e-12)


def test_state_reconstruction_matches_manual_trig():
    v, gamma = 3000.0, math.radians(-15.0)
    state = orbital.compute_orbital_state(R_EARTH + 50_000.0, v, gamma)
    assert state.v_radial == pytest.approx(v * math.sin(gamma))
    assert state.v_tangential == pytest.approx(v * math.cos(gamma))
