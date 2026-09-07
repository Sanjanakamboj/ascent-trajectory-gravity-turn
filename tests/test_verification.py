"""Independent verification checks (M2 requirement 8: does not rely on 'solve_ivp
returned successfully' alone).

Checks implemented here:
  B. zero-thrust ballistic limit (energy / angular momentum conservation)
  C. constant radial-thrust short-duration sanity case
  D. mass-flow analytical check (integrated trajectory vs m(t) = m0 - |mdot| t)
  E. Tsiolkovsky consistency (propulsion-only reduction of the full ascent RHS)
  G. integrator convergence across tolerance/max-step settings

(Zero-drag limit and Earth-rotation checks (A, F) live in test_dynamics.py and
test_earth_rotation.py; dimensional/sign checks (H) also live in test_dynamics.py.)
"""

import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, propulsion as prop
from ascent.constants import (
    BASELINE_VEHICLE,
    LAUNCH_LATITUDE_DEG,
    MU_EARTH,
    R_EARTH,
)
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE


def make_params(cd=None, mu=None):
    return dyn.AscentParams(
        isp=V.isp,
        reference_area=V.reference_area,
        drag_coefficient=V.drag_coefficient if cd is None else cd,
        m_min=V.m_min,
        latitude_rad=LAT,
        mu=MU_EARTH if mu is None else mu,
    )


# --- B. Zero-thrust ballistic limit: energy / angular-momentum conservation ---------

def test_ballistic_coast_conserves_energy_and_angular_momentum():
    params = make_params(cd=0.0)  # isolate gravity: no drag in this check
    r0 = R_EARTH + 300_000.0
    v0 = 7500.0  # slightly below local circular speed -> bound elliptical orbit
    gamma0 = math.radians(5.0)
    y0 = np.array([r0, 0.0, v0, gamma0, V.m0])

    coast = controls.coast_control()
    result = run_ascent(y0, (0.0, 600.0), params, coast, max_step=1.0,
                         rtol=1e-11, atol=1e-11)
    assert result.success

    r, v, gamma = result.r, result.v, result.gamma
    energy = 0.5 * v**2 - MU_EARTH / r
    ang_mom = r * v * np.cos(gamma)

    assert np.max(np.abs((energy - energy[0]) / energy[0])) < 1e-8
    assert np.max(np.abs((ang_mom - ang_mom[0]) / ang_mom[0])) < 1e-8


# --- C. Constant radial-thrust short-duration sanity case ---------------------------

def test_short_duration_vertical_thrust_matches_manual_rhs_estimate():
    params = make_params()
    y0 = np.array([R_EARTH, 0.0, 0.0, math.radians(90.0), V.m0])
    control = controls.vertical_thrust_control(V.thrust)

    dydt0 = dyn.ascent_rhs(0.0, y0, params, control)  # analytic RHS at t=0

    dt = 1e-3
    result = run_ascent(y0, (0.0, dt), params, control, max_step=dt / 10,
                         rtol=1e-12, atol=1e-14)
    assert result.success

    # Over a very short dt, forward-Euler using the analytic RHS at t=0 should closely
    # match the integrator's result (both r, v change approx linearly over this dt).
    # A small residual (~second-order Euler truncation error, since T/m gives a large
    # acceleration) is expected; rel=1e-3 comfortably bounds it for r and mass.
    #
    # v's tolerance is deliberately looser (M5 finding): this dt straddles V_FLOOR
    # (v crosses it at ~6.6e-5 s into this 1e-3 s window). Below V_FLOOR, drag is
    # zeroed (dynamics.py's V_FLOOR guard, extended in M5 -- the inertial-velocity
    # direction needed to apply drag "anti-parallel to v" is itself undefined there),
    # so a_drag jumps from 0 to a nonzero value partway through this window. A
    # single-point (t=0) Euler estimate cannot track that step change to rel=1e-3; the
    # ~10% residual this produces is an expected consequence of that one legitimate
    # transition within the window, not a sign or order-of-magnitude error (which
    # rel=0.15 would still comfortably catch).
    y_euler = y0 + dydt0 * dt
    assert result.r[-1] == pytest.approx(y_euler[0], rel=1e-3)
    assert result.v[-1] == pytest.approx(y_euler[2], rel=0.15)
    assert result.m[-1] == pytest.approx(y_euler[4], rel=1e-9)


# --- D. Mass-flow analytical check ---------------------------------------------------

def test_integrated_mass_history_matches_analytic_linear_depletion():
    params = make_params()
    y0 = np.array([R_EARTH, 0.0, 0.0, math.radians(90.0), V.m0])
    control = controls.vertical_thrust_control(V.thrust)

    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    t_eval = np.linspace(0.0, t_burn * 0.999, 50)  # stay before depletion/cutoff
    result = run_ascent(y0, (0.0, t_burn * 0.999), params, control, max_step=0.1,
                         rtol=1e-10, atol=1e-10, t_eval=t_eval)
    assert result.success

    expected = np.array([prop.mass_at_time(V.m0, V.thrust, V.isp, t) for t in t_eval])
    np.testing.assert_allclose(result.m, expected, rtol=1e-6)


def test_propellant_depletion_event_fires_at_expected_time():
    params = make_params()
    y0 = np.array([R_EARTH, 0.0, 0.0, math.radians(90.0), V.m0])
    control = controls.vertical_thrust_control(V.thrust)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)

    result = run_ascent(y0, (0.0, t_burn * 1.5), params, control, max_step=0.5,
                         rtol=1e-10, atol=1e-10)
    assert result.success
    depletion_times = result.t_events[0]
    assert len(depletion_times) == 1
    assert depletion_times[0] == pytest.approx(t_burn, rel=1e-3)

    # mass must never go meaningfully below m_min anywhere in the trajectory -- a
    # sub-milligram-relative floating point/adaptive-step overshoot right at the event
    # crossing is tolerated (1e-2 kg out of a 500,000 kg vehicle), but a real propellant
    # debt (e.g. kg-scale or larger) is not.
    assert np.all(result.m >= V.m_min - 1e-2)


def test_ground_impact_event_terminates_descending_trajectory():
    params = make_params()
    # Ballistic (thrust-off) case launched at a shallow negative flight path angle from
    # a modest altitude so it comes back down within the integration window.
    y0 = np.array([R_EARTH + 2000.0, 0.0, 100.0, math.radians(-30.0), V.m_min])
    coast = controls.coast_control()
    result = run_ascent(y0, (0.0, 120.0), params, coast, max_step=0.5,
                         rtol=1e-10, atol=1e-10)
    assert result.success
    assert len(result.t_events[1]) == 1  # ground-impact event fired exactly once
    assert result.r[-1] == pytest.approx(R_EARTH, rel=1e-6)


# --- E. Tsiolkovsky consistency (propulsion-only reduction of the full ascent RHS) ---

def test_propulsion_only_reduction_matches_tsiolkovsky():
    # mu=0 removes gravity; cd=0 removes drag. What remains of v_dot is exactly
    # T/m (thrust aligned with velocity, alpha=0), i.e. the ideal-rocket-equation regime.
    params = make_params(cd=0.0, mu=0.0)
    y0 = np.array([R_EARTH, 0.0, 1.0, 0.0, V.m0])  # v>0, gamma=0 to avoid v-floor freeze

    def aligned_control(t, y, p):
        return V.thrust, y[3]  # chi = gamma -> alpha = 0

    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    result = run_ascent(y0, (0.0, t_burn * 0.999), params, aligned_control,
                         max_step=0.05, rtol=1e-11, atol=1e-11)
    assert result.success

    mf = result.m[-1]
    delta_v_numeric = result.v[-1] - result.v[0]
    delta_v_ideal = prop.tsiolkovsky_delta_v(V.isp, V.m0, mf)
    assert delta_v_numeric == pytest.approx(delta_v_ideal, rel=1e-4)


# --- G. Integrator convergence -------------------------------------------------------

def test_integrator_convergence_across_tolerances():
    params = make_params()
    y0 = np.array([R_EARTH, 0.0, 0.0, math.radians(88.0), V.m0])
    control = controls.vertical_thrust_control(V.thrust)
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    t_end = t_burn * 0.5

    settings = [
        dict(max_step=2.0, rtol=1e-6, atol=1e-6),
        dict(max_step=0.5, rtol=1e-9, atol=1e-9),
        dict(max_step=0.1, rtol=1e-12, atol=1e-12),
    ]
    finals = []
    for s in settings:
        result = run_ascent(y0, (0.0, t_end), params, control, **s)
        assert result.success
        finals.append((result.r[-1], result.v[-1], result.gamma[-1]))

    finals = np.array(finals)
    # differences between successive refinements must shrink (convergence), and the
    # tightest two settings must agree closely (near-converged solution).
    diffs = np.abs(np.diff(finals, axis=0))
    assert np.all(diffs[1] < diffs[0])
    np.testing.assert_allclose(finals[1], finals[2], rtol=1e-5)
