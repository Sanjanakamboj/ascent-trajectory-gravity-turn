"""M3 guidance / gravity-turn verification (DESIGN.md M3 S6, checks A-H).

test_orbital.py covers check G (orbital-element consistency) in more depth; a repeat
of that check via a full-integration burnout state is also included here.
"""

import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, orbital, propulsion as prop
from ascent.constants import BASELINE_VEHICLE, LAUNCH_LATITUDE_DEG, R_EARTH
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE


def make_params():
    return dyn.AscentParams(
        isp=V.isp,
        reference_area=V.reference_area,
        drag_coefficient=V.drag_coefficient,
        m_min=V.m_min,
        latitude_rad=LAT,
    )


def initial_state():
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, V.m0])


# A guidance parameter set independently confirmed (scripts/m3_gravity_turn_sweep.py
# exploration) to survive powered flight all the way to burnout without a premature
# ground-impact "nose dive" -- see DESIGN.md M3 S2/S9 for why an early pitch-kick with
# this vehicle's thrust-to-weight often does NOT survive to burnout (a genuine, real
# gravity-turn sensitivity finding, not a bug). Tests that need a completed burn use
# this configuration; tests about crash-tolerant behavior (e.g. small-kick continuity)
# do not depend on reaching burnout at all.
GOOD_KICK_START = 40.0
GOOD_KICK_ANGLE = math.radians(10.0)
GOOD_KICK_DURATION = 10.0


# --- A. Zero pitch-kick limit -------------------------------------------------------

def test_zero_kick_duration_reduces_to_vertical_then_immediate_gravity_turn():
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=10.0,
                                             kick_angle_rad=math.radians(3.0),
                                             kick_duration=0.0)
    # Before kick_start_time: must be exactly vertical (chi = 90 deg), regardless of
    # kick_angle_rad, since the kick phase has zero width.
    thrust, chi = control(5.0, y0, params)
    assert chi == pytest.approx(math.pi / 2)
    # At/after kick_start_time, with zero duration, must already be zero-AoA (chi=gamma).
    y_mid = np.array([R_EARTH + 1000.0, 0.1, 500.0, math.radians(12.0), 400_000.0])
    thrust2, chi2 = control(15.0, y_mid, params)
    assert chi2 == pytest.approx(y_mid[3])  # chi == gamma exactly


def test_zero_kick_angle_never_deviates_from_vertical_or_gamma_tracking():
    params = make_params()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=10.0,
                                             kick_angle_rad=0.0, kick_duration=5.0)
    y_kick_phase = np.array([R_EARTH, 0.0, 400.0, math.radians(5.0), 450_000.0])
    thrust, chi = control(12.0, y_kick_phase, params)
    # kick_angle_rad = 0 -> the "kicked" angle equals vertical exactly.
    assert chi == pytest.approx(math.pi / 2)


# --- B. Small-kick continuity -------------------------------------------------------

def test_small_changes_in_kick_angle_produce_smooth_metric_changes():
    params = make_params()
    y0 = initial_state()
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)

    angles_deg = [1.0, 1.5, 2.0, 2.5, 3.0]
    energies = []
    for a_deg in angles_deg:
        control = controls.gravity_turn_control(V.thrust, kick_start_time=10.0,
                                                  kick_angle_rad=math.radians(a_deg),
                                                  kick_duration=5.0)
        result = run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5,
                             rtol=1e-9, atol=1e-9, terminal_depletion=True,
                             include_apogee_event=False)
        state = orbital.compute_orbital_state(result.r[-1], result.v[-1], result.gamma[-1])
        energies.append(state.specific_energy)

    energies = np.array(energies)
    diffs = np.diff(energies)
    # No wild jumps: consecutive differences should be comparable in magnitude (no
    # single step more than 5x any neighboring step), indicating a smooth response
    # away from any event/threshold boundary.
    for i in range(len(diffs) - 1):
        ratio = abs(diffs[i + 1]) / max(abs(diffs[i]), 1e-6)
        assert ratio < 5.0


# --- C. Zero-AoA gravity-turn condition ---------------------------------------------

def test_thrust_direction_aligns_with_velocity_after_kick_phase():
    params = make_params()
    y0 = initial_state()
    kick_start, kick_dur = 10.0, 5.0
    control = controls.gravity_turn_control(V.thrust, kick_start, math.radians(2.0), kick_dur)

    t_check = kick_start + kick_dur + 20.0
    y_sample = np.array([R_EARTH + 5000.0, 0.05, 600.0, math.radians(15.0), 420_000.0])
    thrust, chi = control(t_check, y_sample, params)
    assert chi == pytest.approx(y_sample[3], abs=1e-12)  # alpha = chi - gamma = 0 exactly


# --- D. Coast conservation (drag disabled) ------------------------------------------

def test_coast_after_burnout_conserves_energy_and_angular_momentum_no_drag():
    params = dyn.AscentParams(isp=V.isp, reference_area=V.reference_area,
                               drag_coefficient=0.0, m_min=V.m_min, latitude_rad=LAT)
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    burnout = run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5,
                          rtol=1e-10, atol=1e-10, terminal_depletion=True,
                          include_apogee_event=False)
    y_burnout = np.array([burnout.r[-1], burnout.theta[-1], burnout.v[-1],
                           burnout.gamma[-1], burnout.m[-1]])

    coast = controls.coast_control()
    result = run_ascent(y_burnout, (0.0, 600.0), params, coast, max_step=1.0,
                         rtol=1e-11, atol=1e-11, include_apogee_event=False)
    energy = 0.5 * result.v**2 - 3.986004418e14 / result.r
    ang_mom = result.r * result.v * np.cos(result.gamma)
    assert np.max(np.abs((energy - energy[0]) / energy[0])) < 1e-8
    assert np.max(np.abs((ang_mom - ang_mom[0]) / ang_mom[0])) < 1e-8


# --- E. Event consistency: burnout time vs analytic depletion time -----------------

def test_burnout_event_time_matches_analytic_depletion_time():
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn_analytic = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn_analytic * 1.5), params, control,
                         max_step=0.5, rtol=1e-10, atol=1e-10,
                         include_apogee_event=False)
    depletion_times = result.t_events[0]
    assert len(depletion_times) == 1
    assert depletion_times[0] == pytest.approx(t_burn_analytic, rel=1e-3)
    assert result.m[-1] == pytest.approx(V.m_min, abs=1e-2)
    assert np.all(result.m >= V.m_min - 1e-2)


# --- F. State reconstruction ----------------------------------------------------

def test_state_reconstruction_against_total_speed_and_gamma():
    v, gamma = 4321.0, math.radians(17.5)
    vr, vt = orbital.velocity_components(v, gamma)
    v_reconstructed = math.hypot(vr, vt)
    gamma_reconstructed = math.atan2(vr, vt)
    assert v_reconstructed == pytest.approx(v, rel=1e-12)
    assert gamma_reconstructed == pytest.approx(gamma, rel=1e-12)


# --- G. Orbital-element consistency (full-integration burnout state) ---------------

def test_burnout_orbital_elements_self_consistent():
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5,
                         rtol=1e-10, atol=1e-10, terminal_depletion=True,
                         include_apogee_event=False)
    state = orbital.compute_orbital_state(result.r[-1], result.v[-1], result.gamma[-1])
    if state.is_bound:
        energy_from_a = -3.986004418e14 / (2 * state.semi_major_axis)
        assert energy_from_a == pytest.approx(state.specific_energy, rel=1e-6)


# --- H. Guidance repeatability -------------------------------------------------------

def test_guidance_repeatability_identical_outputs():
    params = make_params()
    y0 = initial_state()
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)

    def run_once():
        control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                                  kick_angle_rad=GOOD_KICK_ANGLE,
                                                  kick_duration=GOOD_KICK_DURATION)
        return run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5,
                           rtol=1e-10, atol=1e-10, terminal_depletion=True,
                           include_apogee_event=False)

    r1 = run_once()
    r2 = run_once()
    np.testing.assert_array_equal(r1.t, r2.t)
    np.testing.assert_array_equal(r1.r, r2.r)
    np.testing.assert_array_equal(r1.v, r2.v)
    np.testing.assert_array_equal(r1.gamma, r2.gamma)
    np.testing.assert_array_equal(r1.m, r2.m)


# --- Velocity continuity across the pitch-kick (M3 S10 sanity check) ---------------

def test_velocity_continuous_across_pitch_kick_boundaries():
    # Check each phase boundary independently with a TIGHT symmetric window around it
    # (not a single t_eval array spanning the whole kick, whose middle two samples are
    # ~kick_duration apart and are *expected* to differ by the normal thrust
    # acceleration over that longer interval -- that would not be a discontinuity test).
    params = make_params()
    y0 = initial_state()
    kick_start, kick_dur = GOOD_KICK_START, GOOD_KICK_DURATION
    control = controls.gravity_turn_control(V.thrust, kick_start, math.radians(4.0), kick_dur)

    dt = 0.02
    for boundary in (kick_start, kick_start + kick_dur):
        t_eval = np.array([boundary - dt, boundary, boundary + dt])
        result = run_ascent(y0, (0.0, boundary + dt + 0.5), params, control,
                             max_step=0.005, rtol=1e-11, atol=1e-11, t_eval=t_eval,
                             include_apogee_event=False)
        dv = np.abs(np.diff(result.v))
        # Over a 0.02s window, even at max thrust acceleration (~19 m/s^2 near burnout)
        # the smooth speed change is at most ~0.4 m/s; 2 m/s is a generous bound that
        # would still catch a genuine step discontinuity (which would show up as O(10)
        # to O(100) m/s given this vehicle's speeds).
        assert np.all(dv < 2.0), f"boundary at t={boundary}: dv={dv}"
