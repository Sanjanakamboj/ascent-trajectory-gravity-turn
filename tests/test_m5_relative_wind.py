"""M5 atmospheric-relative-wind regression tests (check L, DESIGN.md M5 S3).

``dynamics.relative_speed`` gained an ``azimuth_rad`` parameter for M5 (default 90 deg,
due east). These tests confirm: (1) the default exactly reproduces the pre-M5 M1-M4
formula (already covered by the full M1-M3 regression, test_earth_rotation.py, which
still passes unchanged), and (2) the new cross-track term behaves correctly for
non-due-east azimuths.
"""

import math

import pytest

from ascent import dynamics as dyn
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH

LAT = math.radians(LAUNCH_LATITUDE_DEG)


def test_default_azimuth_matches_pre_m5_formula_exactly():
    r, v, gamma = R_EARTH + 50_000.0, 2000.0, math.radians(10.0)
    v_rel_new = dyn.relative_speed(v, gamma, r, LAT)  # default azimuth = 90 deg
    v_rel_new_explicit = dyn.relative_speed(v, gamma, r, LAT, azimuth_rad=math.pi / 2)
    assert v_rel_new == pytest.approx(v_rel_new_explicit, rel=1e-15)
    # Reconstruct the pre-M5 (in-plane-only) formula directly and confirm identity.
    v_atm = dyn.atmosphere_corotation_speed(r, LAT)
    v_radial = v * math.sin(gamma)
    v_tangential = v * math.cos(gamma)
    pre_m5 = math.hypot(v_radial, v_tangential - v_atm)
    assert v_rel_new == pytest.approx(pre_m5, rel=1e-12)


def test_polar_azimuth_adds_full_cross_track_term():
    # azimuth=0 (polar): in-plane atmosphere component vanishes, but the FULL v_atm
    # magnitude must still appear as a cross-track contribution to v_rel -- it must
    # NOT be silently dropped.
    r = R_EARTH
    v, gamma = 0.0, 0.0
    v_atm = dyn.atmosphere_corotation_speed(r, LAT)
    v_rel = dyn.relative_speed(v, gamma, r, LAT, azimuth_rad=0.0)
    assert v_rel == pytest.approx(v_atm, rel=1e-9)
    # And this must be strictly greater than what a (wrong) in-plane-only treatment
    # would give at this azimuth (which would incorrectly report v_rel=0 here).
    assert v_rel > 0.0


def test_relative_wind_increases_monotonically_toward_polar_when_matching_inplane_speed():
    # Give the vehicle exactly the IN-PLANE co-rotation speed for each azimuth (so the
    # in-plane relative-velocity term cancels to zero); what remains of v_rel is then
    # purely the cross-track term, which should grow monotonically from 0 (due east)
    # to the full v_atm magnitude (polar) as azimuth decreases from 90 to 0 deg.
    from ascent import inclination as inc

    r = R_EARTH
    azimuths_deg = [90.0, 70.0, 55.0, 45.0, 35.0, 0.0]
    v_rels = []
    for a_deg in azimuths_deg:
        az = math.radians(a_deg)
        v_inplane = inc.useful_rotational_boost(az, LAT, r=r)
        v_rels.append(dyn.relative_speed(v_inplane, 0.0, r, LAT, azimuth_rad=az))
    assert all(v_rels[i] < v_rels[i + 1] for i in range(len(v_rels) - 1))
    assert v_rels[0] == pytest.approx(0.0, abs=1e-6)                    # due east
    assert v_rels[-1] == pytest.approx(dyn.atmosphere_corotation_speed(r, LAT), rel=1e-9)  # polar


def test_cross_track_term_matches_inclination_module():
    # As above: match the vehicle's tangential speed to the in-plane co-rotation
    # component so only the cross-track term remains in v_rel, and confirm it equals
    # inclination.cross_track_rotational_component exactly.
    from ascent import inclination as inc

    r = R_EARTH + 100_000.0
    az = math.radians(40.743)  # ~55 deg inclination azimuth (see test_inclination.py)
    v_inplane = inc.useful_rotational_boost(az, LAT, r=r)
    v_rel = dyn.relative_speed(v_inplane, 0.0, r, LAT, azimuth_rad=az)
    cross = inc.cross_track_rotational_component(az, LAT, r=r)
    assert v_rel == pytest.approx(cross, rel=1e-6)
