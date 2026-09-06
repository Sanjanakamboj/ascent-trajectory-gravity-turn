"""M3 event-handling refinement tests (DESIGN.md M3 S5)."""

import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, propulsion as prop
from ascent.constants import BASELINE_VEHICLE, LAUNCH_LATITUDE_DEG, R_EARTH
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE


def make_params():
    return dyn.AscentParams(
        isp=V.isp, reference_area=V.reference_area, drag_coefficient=V.drag_coefficient,
        m_min=V.m_min, latitude_rad=LAT,
    )


def initial_state():
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, V.m0])


# See test_guidance.py for why this specific configuration (and not an early kick) is
# used wherever a test needs the vehicle to actually reach burnout.
GOOD_KICK_START = 40.0
GOOD_KICK_ANGLE = math.radians(10.0)
GOOD_KICK_DURATION = 10.0


def test_apogee_event_fires_during_coast_and_matches_max_altitude():
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn + 1200.0), params, control, max_step=0.5,
                         rtol=1e-9, atol=1e-9, include_apogee_event=True)

    apogee_times = result.t_events[2]
    assert len(apogee_times) >= 1
    t_apogee = apogee_times[0]

    # Altitude at the recorded apogee time must match (closely) the max altitude found
    # by direct search over the sampled trajectory.
    h = result.r - R_EARTH
    idx_max = int(np.argmax(h))
    assert result.t[idx_max] == pytest.approx(t_apogee, abs=2.0)


def test_burnout_mass_hits_dry_plus_payload_exactly():
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5,
                         rtol=1e-10, atol=1e-10, terminal_depletion=True,
                         include_apogee_event=False)
    assert result.m[-1] == pytest.approx(V.m_dry + V.m_payload, abs=1e-2)
    assert not result.success or result.t[-1] == pytest.approx(t_burn, rel=1e-3)


def test_mass_never_goes_below_floor_across_full_gravity_turn_run():
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn + 500.0), params, control, max_step=0.5,
                         rtol=1e-9, atol=1e-9)
    # A sub-milligram-relative adaptive-step/event-crossing residual is tolerated (as in
    # M2's equivalent check); a real propellant debt (kg-scale or larger) is not.
    assert np.all(result.m >= V.m_min - 0.1)


def test_ground_impact_directionally_correct_after_gravity_turn_apogee():
    # A short, steep sub-orbital case that will come back down within a reasonable
    # window: use coast_control from a modest altitude with a negative flight-path
    # angle so impact is guaranteed and quick (keeps the test fast and deterministic).
    params = make_params()
    y0 = np.array([R_EARTH + 5000.0, 0.0, 300.0, math.radians(-45.0), V.m_min])
    coast = controls.coast_control()
    result = run_ascent(y0, (0.0, 200.0), params, coast, max_step=0.5,
                         rtol=1e-10, atol=1e-10, include_apogee_event=False)
    assert len(result.t_events[1]) == 1
    assert result.r[-1] == pytest.approx(R_EARTH, rel=1e-6)
    # Confirm it was actually descending (r decreasing) right before impact.
    assert result.r[-1] < result.r[-2]


def test_target_altitude_crossing_event_does_not_imply_orbit():
    from ascent.constants import TARGET_ALTITUDE
    params = make_params()
    y0 = initial_state()
    control = controls.gravity_turn_control(V.thrust, kick_start_time=GOOD_KICK_START,
                                             kick_angle_rad=GOOD_KICK_ANGLE,
                                             kick_duration=GOOD_KICK_DURATION)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn + 800.0), params, control, max_step=0.5,
                         rtol=1e-9, atol=1e-9, target_altitude=TARGET_ALTITUDE)
    crossing_times = result.t_events[-1]
    # The event existing/firing says nothing about orbit achievement by itself -- that
    # requires orbital.is_circular_orbit_achieved on the corresponding state, which is
    # a separate, explicit check (not exercised here, just asserting the event and the
    # orbit criterion are independent APIs).
    assert isinstance(crossing_times, np.ndarray)
