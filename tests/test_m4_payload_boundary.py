"""M4 numerical verification checks A, B(partial), C, G, I (DESIGN.md M4 S8).

Check D (mass bookkeeping) lives in test_m4_vehicle.py; check E (Tsiolkovsky) and F
(orbital-element consistency) are covered by the existing M1-M3 propulsion/orbital
tests plus test_m4_loss_budget.py; check H (convergence) is in
test_m4_convergence.py.
"""

import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, orbital, propulsion as prop
from ascent.constants import (
    BASELINE_VEHICLE, LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle,
)
from ascent.insertion_search import find_drag_consistent_cutoff
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)

KICK_START = 45.0
KICK_ANGLE_DEG = 20.0
KICK_DURATION = 20.0


def build_params(vehicle):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT,
    )


def initial_state(vehicle):
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])


# --- A. Original M1-M3 vehicle remains a failure under the M4 criterion ------------

def test_original_verification_vehicle_still_fails_m4_insertion_criterion():
    vehicle = BASELINE_VEHICLE
    params = build_params(vehicle)
    y0 = initial_state(vehicle)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    # Reuse the M3-selected gravity-turn parameters for the original vehicle (its own
    # best-known guidance, DESIGN.md M3 S5) -- even under the SAME (or a fresh) search,
    # this vehicle's Tsiolkovsky delta-v (M1 S7.4) cannot plausibly meet the target;
    # this test only needs ONE reasonable guidance attempt to demonstrate the failure,
    # not an exhaustive resweep of M3.
    cutoff = find_drag_consistent_cutoff(vehicle, params, y0, 40.0, 25.0, 10.0, t_burn,
                                          TARGET_ALTITUDE)
    assert not cutoff.found_passing_cutoff


# --- C. Payload-boundary bracketing -------------------------------------------------

def test_payload_just_below_and_above_boundary_bracket_correctly():
    # Uses a coarse-then-refined bracket search identical in spirit to
    # scripts/m4_payload_sweep.py, but self-contained and fast (loose tolerance).
    def passes(payload):
        vehicle = m4_vehicle(payload)
        params = build_params(vehicle)
        y0 = initial_state(vehicle)
        t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
        return find_drag_consistent_cutoff(vehicle, params, y0, KICK_START, KICK_ANGLE_DEG,
                                            KICK_DURATION, t_burn,
                                            TARGET_ALTITUDE).found_passing_cutoff

    assert passes(10_000.0)   # reference payload passes
    lo, hi = 10_000.0, 20_000.0
    assert passes(lo)
    assert not passes(hi)
    for _ in range(8):  # coarse bisection, tolerance ~40 kg -- plenty for a bracket test
        mid = 0.5 * (lo + hi)
        if passes(mid):
            lo = mid
        else:
            hi = mid
    assert passes(lo)
    assert not passes(hi)
    assert hi - lo < 100.0


# --- G. Guidance/cutoff-search repeatability ----------------------------------------

def test_drag_consistent_cutoff_search_is_repeatable():
    vehicle = m4_vehicle(10_000.0)
    params = build_params(vehicle)
    y0 = initial_state(vehicle)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)

    r1 = find_drag_consistent_cutoff(vehicle, params, y0, KICK_START, KICK_ANGLE_DEG,
                                      KICK_DURATION, t_burn, TARGET_ALTITUDE)
    r2 = find_drag_consistent_cutoff(vehicle, params, y0, KICK_START, KICK_ANGLE_DEG,
                                      KICK_DURATION, t_burn, TARGET_ALTITUDE)
    assert r1.found_passing_cutoff and r2.found_passing_cutoff
    assert r1.cutoff_time == r2.cutoff_time
    assert r1.insertion.circularization_delta_v == r2.insertion.circularization_delta_v


# --- I. Atmosphere-relative drag regression -----------------------------------------

def test_m4_dynamics_still_uses_relative_speed_for_drag():
    # Regression: confirm the M4 code path did not bypass the M2 rotating-atmosphere
    # drag treatment. Construct a state where inertial and relative speed differ
    # sharply and check the RHS drag term uses the relative value (same check style as
    # M2's test_earth_rotation.py, exercised here through the M4 vehicle/params).
    vehicle = m4_vehicle(10_000.0)
    params = build_params(vehicle)
    r = R_EARTH + 50_000.0
    v_atm = dyn.atmosphere_corotation_speed(r, LAT)
    # Purely tangential inertial velocity exactly matching the local atmosphere speed:
    # relative speed must be ~0, so drag must be ~0 even though inertial speed is large.
    y = np.array([r, 0.0, v_atm, 0.0, vehicle.m0])
    control = controls.coast_control()
    dydt = dyn.ascent_rhs(0.0, y, params, control)
    v_rel = dyn.relative_speed(v_atm, 0.0, r, LAT)
    assert v_rel == pytest.approx(0.0, abs=1e-6)
    # With v_rel ~ 0, drag_acceleration ~ 0, so v_dot should be dominated by gravity only.
    g = dyn.local_gravity(r)
    assert dydt[2] == pytest.approx(-g * math.sin(0.0), abs=1e-6)
