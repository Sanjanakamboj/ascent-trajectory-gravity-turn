"""M5 direct-ascent geometry tests (DESIGN.md M5 S2, checks B/C/D/E)."""

import math

import pytest

from ascent import inclination as inc
from ascent.constants import LAUNCH_LATITUDE_DEG, OMEGA_EARTH, R_EARTH

LAT = math.radians(LAUNCH_LATITUDE_DEG)


# --- B. Azimuth/inclination round trip ----------------------------------------------

@pytest.mark.parametrize("i_deg", [28.5, 30.0, 35.0, 45.0, 55.0, 70.0, 85.0, 90.0])
def test_inclination_azimuth_round_trip(i_deg):
    i = math.radians(i_deg)
    az = inc.azimuth_from_inclination(i, LAT)
    i_back = inc.inclination_from_azimuth(az, LAT)
    assert i_back == pytest.approx(i, abs=1e-9)


# --- C. Impossible inclinations rejected --------------------------------------------

def test_inclination_below_latitude_raises():
    with pytest.raises(ValueError):
        inc.azimuth_from_inclination(math.radians(20.0), LAT)


def test_retrograde_inclination_raises():
    with pytest.raises(ValueError):
        inc.azimuth_from_inclination(math.radians(120.0), LAT)


# --- D. Polar limit: useful eastward contribution -> 0 ------------------------------

def test_polar_limit_zero_useful_boost_full_cross_track():
    az = inc.azimuth_from_inclination(math.radians(90.0), LAT)
    assert az == pytest.approx(0.0, abs=1e-9)
    useful = inc.useful_rotational_boost(az, LAT)
    cross = inc.cross_track_rotational_component(az, LAT)
    v_rot = inc.site_rotational_speed(LAT)
    assert useful == pytest.approx(0.0, abs=1e-6)
    assert cross == pytest.approx(v_rot, rel=1e-9)


# --- E. Due-east limit reproduces M1/M4 initial condition ---------------------------

def test_due_east_limit_matches_m1_m4_value():
    az = inc.azimuth_from_inclination(LAT, LAT)
    assert az == pytest.approx(math.pi / 2, abs=1e-9)
    useful = inc.useful_rotational_boost(az, LAT)
    assert useful == pytest.approx(408.73884297255603, rel=1e-9)  # M1 S7 hand value
    cross = inc.cross_track_rotational_component(az, LAT)
    assert cross == pytest.approx(0.0, abs=1e-6)


# --- Closed-form cross-check: useful_boost(i) == omega*R*cos(i) ---------------------

@pytest.mark.parametrize("i_deg", [28.5, 35.0, 45.0, 55.0, 70.0, 90.0])
def test_useful_boost_matches_closed_form(i_deg):
    i = math.radians(i_deg)
    az = inc.azimuth_from_inclination(i, LAT)
    useful = inc.useful_rotational_boost(az, LAT)
    closed_form = OMEGA_EARTH * R_EARTH * math.cos(i)
    assert useful == pytest.approx(closed_form, rel=1e-9)


# --- Monotonic decrease of useful boost with increasing inclination -----------------

def test_useful_boost_monotonically_decreases_with_inclination():
    inclinations_deg = [28.5, 35.0, 45.0, 55.0, 70.0, 90.0]
    boosts = []
    for i_deg in inclinations_deg:
        i = math.radians(i_deg)
        az = inc.azimuth_from_inclination(i, LAT)
        boosts.append(inc.useful_rotational_boost(az, LAT))
    assert all(boosts[i] > boosts[i + 1] for i in range(len(boosts) - 1))


# --- Vector decomposition self-consistency: useful^2 + cross^2 == v_rot^2 -----------

@pytest.mark.parametrize("i_deg", [28.5, 45.0, 70.0, 90.0])
def test_useful_and_cross_form_full_rotational_vector(i_deg):
    i = math.radians(i_deg)
    az = inc.azimuth_from_inclination(i, LAT)
    useful = inc.useful_rotational_boost(az, LAT)
    cross = inc.cross_track_rotational_component(az, LAT)
    v_rot = inc.site_rotational_speed(LAT)
    assert math.hypot(useful, cross) == pytest.approx(v_rot, rel=1e-9)
