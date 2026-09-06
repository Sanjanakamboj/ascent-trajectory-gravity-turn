"""M4 insertion-criterion tests (DESIGN.md M4 S2)."""

import math

import pytest

from ascent import orbital
from ascent.constants import MU_EARTH, R_EARTH, TARGET_ALTITUDE


def test_circular_state_needs_almost_no_circularization():
    r = R_EARTH + TARGET_ALTITUDE
    v = math.sqrt(MU_EARTH / r)
    state = orbital.compute_orbital_state(r, v, gamma=0.0)
    dv = orbital.circularization_delta_v_at_apogee(state)
    assert dv == pytest.approx(0.0, abs=1e-6)
    ins = orbital.evaluate_orbit_insertion(state, TARGET_ALTITUDE)
    assert ins.meets_insertion_criterion


def test_eccentric_transfer_with_apogee_at_target_needs_positive_circularization_dv():
    # Perigee well below target, apogee exactly at target: this is a classic transfer
    # ellipse: circularizing at apogee should need a clear positive delta-v.
    r_apogee = R_EARTH + TARGET_ALTITUDE
    r_perigee = R_EARTH + 100_000.0
    a = 0.5 * (r_apogee + r_perigee)
    v_apogee = math.sqrt(MU_EARTH * (2.0 / r_apogee - 1.0 / a))
    state = orbital.compute_orbital_state(r_apogee, v_apogee, gamma=0.0)
    assert state.apogee_altitude == pytest.approx(TARGET_ALTITUDE, rel=1e-6)
    dv = orbital.circularization_delta_v_at_apogee(state)
    assert dv > 50.0  # a real, nontrivial circularization burn


def test_hyperbolic_state_has_no_circularization_dv():
    state = orbital.compute_orbital_state(R_EARTH + 300_000.0, 12_000.0, gamma=0.0)
    assert state.apogee_radius is None
    assert orbital.circularization_delta_v_at_apogee(state) is None
    ins = orbital.evaluate_orbit_insertion(state, TARGET_ALTITUDE)
    assert not ins.meets_insertion_criterion
    assert not ins.is_bound


def test_apogee_far_from_target_fails_even_with_small_circularization_dv():
    # Apogee at 100 km (well outside tolerance of 400 km target), near-circular there
    # (tiny circularization dv) -- must still fail: an idealized burn cannot move WHERE
    # the apogee is, only its shape.
    r = R_EARTH + 100_000.0
    v = math.sqrt(MU_EARTH / r) - 5.0  # just below local circular speed
    state = orbital.compute_orbital_state(r, v, gamma=0.0)
    ins = orbital.evaluate_orbit_insertion(state, TARGET_ALTITUDE)
    assert not ins.meets_altitude_target
    assert not ins.meets_insertion_criterion


def test_excessive_circularization_dv_fails_the_cap_even_with_apogee_on_target():
    # Apogee exactly at target, but very low perigee (near-vertical arrival) -> large
    # circularization dv -- must fail the documented delta_v_cap even though apogee
    # matches (DESIGN.md M4 S2 condition 4).
    r_apogee = R_EARTH + TARGET_ALTITUDE
    r_perigee = 500_000.0  # an ABSOLUTE radius deep inside the Earth (not an altitude)
    # -- purely a synthetic orbital-mechanics construction to force high eccentricity
    # and a large circularization delta-v; not meant to represent an achievable state.
    a = 0.5 * (r_apogee + r_perigee)
    v_apogee = math.sqrt(MU_EARTH * (2.0 / r_apogee - 1.0 / a))
    state = orbital.compute_orbital_state(r_apogee, v_apogee, gamma=0.0)
    ins = orbital.evaluate_orbit_insertion(state, TARGET_ALTITUDE, delta_v_cap=1500.0)
    assert ins.meets_altitude_target
    assert ins.circularization_delta_v > 1500.0
    assert not ins.circularization_within_cap
    assert not ins.meets_insertion_criterion


def test_reaches_bound_orbit_is_independent_flag_from_full_criterion():
    # A bound, non-intersecting orbit whose apogee is nowhere near the target: passes
    # the coarse "reaches some bound orbit" sanity floor but not the full criterion.
    r = R_EARTH + TARGET_ALTITUDE
    v = math.sqrt(MU_EARTH / r)
    state = orbital.compute_orbital_state(r, v, gamma=0.0)
    ins_wrong_target = orbital.evaluate_orbit_insertion(state, TARGET_ALTITUDE + 5_000_000.0)
    assert ins_wrong_target.reaches_bound_orbit
    assert not ins_wrong_target.meets_insertion_criterion
