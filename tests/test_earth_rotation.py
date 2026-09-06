"""Tests for Earth-rotation / atmosphere-relative-velocity handling (DESIGN.md S3 Note,
M2 dynamics.py docstring). These are written specifically so that they FAIL if inertial
speed were accidentally substituted for atmosphere-relative speed.
"""

import math

import pytest

from ascent import dynamics as dyn
from ascent.constants import LAUNCH_LATITUDE_DEG, OMEGA_EARTH, R_EARTH

LAT = math.radians(LAUNCH_LATITUDE_DEG)


def test_surface_rotational_speed_matches_m1_hand_calc():
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    assert v_rot == pytest.approx(408.73884297255603, rel=1e-9)


def test_stationary_on_pad_has_zero_air_relative_speed():
    # A vehicle sitting on the pad, co-rotating with the Earth, has inertial tangential
    # velocity exactly equal to the local atmospheric co-rotation speed and zero radial
    # velocity. Represent this directly via (v, gamma) = (v_rot, 0) i.e. purely
    # horizontal motion at the surface -- relative speed must be (approximately) zero.
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    v_rel = dyn.relative_speed(v=v_rot, gamma=0.0, r=R_EARTH, latitude_rad=LAT)
    assert v_rel == pytest.approx(0.0, abs=1e-9)


def test_due_east_inertial_velocity_gets_correct_rotation_subtraction():
    # Purely horizontal (gamma=0) inertial speed v, at radius r: relative speed must be
    # exactly |v - v_atm(r)|, not v itself.
    r = R_EARTH + 50_000.0
    v = 2000.0
    v_atm = dyn.atmosphere_corotation_speed(r, LAT)
    v_rel = dyn.relative_speed(v=v, gamma=0.0, r=r, latitude_rad=LAT)
    assert v_rel == pytest.approx(abs(v - v_atm), rel=1e-12)
    # This is the failure-mode check: relative speed must differ from raw inertial speed
    # whenever v_atm is non-negligible.
    assert v_rel != pytest.approx(v, rel=1e-6)


def test_zero_earth_rotation_limit_reduces_to_inertial_speed():
    r = R_EARTH + 50_000.0
    v = 2000.0
    gamma = math.radians(20.0)
    v_rel = dyn.relative_speed(v=v, gamma=gamma, r=r, latitude_rad=LAT, omega_earth=0.0)
    assert v_rel == pytest.approx(v, rel=1e-12)


def test_relative_speed_uses_full_vector_not_just_tangential_component():
    # With a purely radial (vertical) inertial velocity (gamma=90deg), Earth rotation of
    # the atmosphere (tangential) still leaves a nonzero relative tangential component if
    # the vehicle itself has zero tangential (inertial) speed -- i.e. v_rel should not
    # collapse to the radial component alone unless v_atm is also zero.
    r = R_EARTH
    v = 100.0
    gamma = math.pi / 2  # purely vertical inertial velocity, no tangential component
    v_atm = dyn.atmosphere_corotation_speed(r, LAT)
    v_rel = dyn.relative_speed(v=v, gamma=gamma, r=r, latitude_rad=LAT)
    expected = math.sqrt(v**2 + v_atm**2)  # radial=v, relative tangential = 0 - v_atm
    assert v_rel == pytest.approx(expected, rel=1e-9)
