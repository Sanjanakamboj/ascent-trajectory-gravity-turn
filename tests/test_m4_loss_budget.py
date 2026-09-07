"""M4 delta-v/loss budget tests (DESIGN.md M4 S9)."""

import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, m4_vehicle
from ascent.loss_budget import compute_delta_v_budget
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)


def test_identity_holds_exactly_for_a_real_trajectory():
    vehicle = m4_vehicle(10_000.0)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=vehicle.drag_coefficient,
                               m_min=vehicle.m_min, latitude_rad=LAT)
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    y0 = np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])
    control = controls.gravity_turn_control(vehicle.thrust, 46.0, math.radians(22.0), 20.0)

    t_end = 250.0
    t_eval = np.linspace(0.0, t_end, 5000)
    result = run_ascent(y0, (0.0, t_end), params, control, max_step=0.2,
                         rtol=1e-10, atol=1e-10, terminal_depletion=False,
                         include_apogee_event=False, t_eval=t_eval)

    budget = compute_delta_v_budget(result.t, result.r, result.v, result.gamma, result.m,
                                     vehicle.isp, vehicle.drag_coefficient,
                                     vehicle.reference_area, LAT, vehicle.m0)

    assert budget.residual == pytest.approx(0.0, abs=1e-6)
    reconstructed = (budget.achieved_dv + budget.gravity_loss + budget.drag_loss
                      + budget.steering_loss)
    assert reconstructed == pytest.approx(budget.ideal_dv_to_cutoff, abs=1e-6)


def test_losses_are_nonnegative_for_a_climbing_trajectory():
    vehicle = m4_vehicle(10_000.0)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=vehicle.drag_coefficient,
                               m_min=vehicle.m_min, latitude_rad=LAT)
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    y0 = np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])
    control = controls.gravity_turn_control(vehicle.thrust, 46.0, math.radians(22.0), 20.0)
    t_eval = np.linspace(0.0, 250.0, 5000)
    result = run_ascent(y0, (0.0, 250.0), params, control, max_step=0.2, rtol=1e-10,
                         atol=1e-10, terminal_depletion=False, include_apogee_event=False,
                         t_eval=t_eval)
    budget = compute_delta_v_budget(result.t, result.r, result.v, result.gamma, result.m,
                                     vehicle.isp, vehicle.drag_coefficient,
                                     vehicle.reference_area, LAT, vehicle.m0)
    assert budget.gravity_loss > 0.0
    assert budget.drag_loss > 0.0
    assert budget.steering_loss > 0.0  # nonzero angle-of-attack during the pitch-kick


def test_zero_drag_zero_alpha_case_matches_pure_tsiolkovsky():
    # mu=0 (no gravity), Cd=0 (no drag), alpha=0 throughout: achieved_dv must equal the
    # ideal Tsiolkovsky value exactly (all three loss terms ~0).
    vehicle = m4_vehicle(10_000.0)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=0.0, m_min=vehicle.m_min,
                               latitude_rad=LAT, mu=0.0)
    y0 = np.array([R_EARTH, 0.0, 1.0, 0.0, vehicle.m0])

    def aligned_control(t, y, p):
        return vehicle.thrust, y[3]

    t_end = 60.0
    t_eval = np.linspace(0.0, t_end, 2000)
    result = run_ascent(y0, (0.0, t_end), params, aligned_control, max_step=0.05,
                         rtol=1e-11, atol=1e-11, terminal_depletion=False,
                         include_apogee_event=False, t_eval=t_eval)
    budget = compute_delta_v_budget(result.t, result.r, result.v, result.gamma, result.m,
                                     vehicle.isp, 0.0, vehicle.reference_area, LAT,
                                     vehicle.m0, mu=0.0)
    assert budget.gravity_loss == pytest.approx(0.0, abs=1e-9)
    assert budget.drag_loss == pytest.approx(0.0, abs=1e-9)
    assert budget.steering_loss == pytest.approx(0.0, abs=1e-6)
    assert budget.achieved_dv == pytest.approx(budget.ideal_dv_to_cutoff, rel=1e-6)


def test_pitch_kick_alone_produces_positive_steering_loss():
    # A pure vertical burn (chi always 90deg, gamma stays ~0 briefly at high v) has
    # alpha != 0 throughout if gamma != 90 -- construct a short case with a large,
    # sustained angle of attack and confirm steering_loss is clearly positive and
    # material (not a rounding artifact).
    # mu=0 isolates steering loss from gravity loss -- otherwise a near-90deg gamma
    # here would also rack up a large gravity loss (g*sin(gamma)~g), which is a
    # separate, already-tested term (test_losses_are_nonnegative_for_a_climbing_trajectory).
    # v0 is deliberately large (not the usual near-zero/floor value): with mu=0 (no
    # gravity turning) and gamma_dot = T*sin(alpha)/(m*v), a large v keeps gamma_dot
    # small, so gamma stays close to its initial 0 (i.e. alpha stays close to 90 deg,
    # chi - gamma) for the whole short test window instead of rapidly tracking chi.
    vehicle = m4_vehicle(10_000.0)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=0.0, m_min=vehicle.m_min,
                               latitude_rad=LAT, mu=0.0)
    y0 = np.array([R_EARTH, 0.0, 5000.0, 0.0, vehicle.m0])
    vertical = controls.vertical_thrust_control(vehicle.thrust)  # chi=90, gamma~0 -> alpha~90
    t_end = 2.0
    t_eval = np.linspace(0.0, t_end, 500)
    result = run_ascent(y0, (0.0, t_end), params, vertical, max_step=0.005, rtol=1e-12,
                         atol=1e-12, terminal_depletion=False, include_apogee_event=False,
                         t_eval=t_eval)
    budget = compute_delta_v_budget(result.t, result.r, result.v, result.gamma, result.m,
                                     vehicle.isp, 0.0, vehicle.reference_area, LAT,
                                     vehicle.m0, mu=0.0)
    # alpha~90deg means cos(alpha)~0, so nearly ALL of the ideal dv shows up as
    # steering loss and almost none as achieved dv.
    assert budget.steering_loss > 0.9 * budget.ideal_dv_to_cutoff
