"""M5 numerical verification checks A, F, G, I, J (DESIGN.md M5 S11).

Check D (polar limit) and E (due-east limit) live in test_inclination.py; check L
(atmospheric relative wind) lives in test_m5_relative_wind.py; check H (convergence)
lives in test_m5_convergence.py.
"""

import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, inclination as inc, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle
from ascent.insertion_search import find_drag_consistent_cutoff
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)

# The 28.5 deg guidance from DESIGN.md M4 S4 / M5 S1.
GUIDANCE_28P5 = (45.0, 20.0, 20.0)


def build_params(vehicle, azimuth_rad):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT, azimuth_rad=azimuth_rad,
    )


def initial_state(vehicle, azimuth_rad):
    v0 = inc.useful_rotational_boost(azimuth_rad, LAT)
    return np.array([R_EARTH, 0.0, v0, 0.0, vehicle.m0])


def evaluate(payload, azimuth_rad, guidance):
    vehicle = m4_vehicle(payload)
    params = build_params(vehicle, azimuth_rad)
    y0 = initial_state(vehicle, azimuth_rad)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    return find_drag_consistent_cutoff(vehicle, params, y0, *guidance, t_burn,
                                        TARGET_ALTITUDE)


# --- A. M4 baseline reproduction (28.5 deg via the generalized M5 pathway) ---------

def test_28p5_deg_via_generalized_pathway_reproduces_m4():
    az = inc.azimuth_from_inclination(math.radians(28.5), LAT)
    assert az == pytest.approx(math.pi / 2, abs=1e-9)
    res = evaluate(10_000.0, az, GUIDANCE_28P5)
    assert res.found_passing_cutoff
    assert res.cutoff_time == pytest.approx(250.87102188012523, rel=1e-6)
    assert res.post_atmosphere_state.apogee_altitude / 1000.0 == pytest.approx(
        400.0345022622496, rel=1e-6)
    assert res.insertion.circularization_delta_v == pytest.approx(96.52016790212474, rel=1e-4)

    # Refined boundary: 10,333 kg passes, 10,334 kg fails (DESIGN.md M4 S7).
    assert evaluate(10_333.0, az, GUIDANCE_28P5).found_passing_cutoff
    assert not evaluate(10_334.0, az, GUIDANCE_28P5).found_passing_cutoff


# --- F. Payload bookkeeping (burnout mass = dry + payload, at full depletion) -------

def test_burnout_mass_equals_dry_plus_payload_at_full_depletion():
    payload = 10_000.0
    az = inc.azimuth_from_inclination(math.radians(55.0), LAT)
    vehicle = m4_vehicle(payload)
    params = build_params(vehicle, az)
    y0 = initial_state(vehicle, az)
    control = controls.gravity_turn_control(vehicle.thrust, 43.0, math.radians(21.0), 20.0)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    result = run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5, rtol=1e-10,
                         atol=1e-10, terminal_depletion=True, include_apogee_event=False)
    assert result.m[-1] == pytest.approx(vehicle.m_dry + payload, abs=1e-2)


# --- G. Payload boundary (one below passes, one above fails) at a non-28.5 inclination --

def test_payload_boundary_at_mid_inclination():
    az = inc.azimuth_from_inclination(math.radians(45.0), LAT)
    guidance = (43.0, 21.0, 20.0)
    lo, hi = 8_000.0, 15_000.0
    assert evaluate(lo, az, guidance).found_passing_cutoff
    assert not evaluate(hi, az, guidance).found_passing_cutoff
    for _ in range(10):
        mid = 0.5 * (lo + hi)
        if evaluate(mid, az, guidance).found_passing_cutoff:
            lo = mid
        else:
            hi = mid
    assert evaluate(lo, az, guidance).found_passing_cutoff
    assert not evaluate(hi, az, guidance).found_passing_cutoff
    assert hi - lo < 20.0


# --- I. Guidance/cutoff-search repeatability at a non-28.5 inclination -------------

def test_repeatability_at_non_28p5_inclination():
    az = inc.azimuth_from_inclination(math.radians(70.0), LAT)
    guidance = (43.0, 21.0, 20.0)
    r1 = evaluate(10_000.0, az, guidance)
    r2 = evaluate(10_000.0, az, guidance)
    assert r1.found_passing_cutoff == r2.found_passing_cutoff
    if r1.found_passing_cutoff:
        assert r1.cutoff_time == r2.cutoff_time
        assert r1.insertion.circularization_delta_v == r2.insertion.circularization_delta_v
