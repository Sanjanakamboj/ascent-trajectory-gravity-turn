import math

import numpy as np
import pytest

from ascent import controls, dynamics as dyn
from ascent.constants import (
    BASELINE_VEHICLE,
    LAUNCH_LATITUDE_DEG,
    MU_EARTH,
    R_EARTH,
)

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE


def make_params(cd=None):
    return dyn.AscentParams(
        isp=V.isp,
        reference_area=V.reference_area,
        drag_coefficient=V.drag_coefficient if cd is None else cd,
        m_min=V.m_min,
        latitude_rad=LAT,
    )


# --- Dimensional / sign checks (M2 requirement H) -----------------------------------

def test_gravity_matches_mu_over_r_squared():
    r = R_EARTH + 400e3
    g = dyn.local_gravity(r)
    assert g == pytest.approx(MU_EARTH / r**2)


def test_gravity_decreases_with_altitude():
    g_surface = dyn.local_gravity(R_EARTH)
    g_orbit = dyn.local_gravity(R_EARTH + 400e3)
    assert 0 < g_orbit < g_surface


def test_altitude_radius_conversion():
    r = R_EARTH + 12_345.0
    h = r - R_EARTH
    assert h == pytest.approx(12_345.0)


def test_drag_acceleration_is_nonnegative_and_reduces_speed_gain():
    # Drag as implemented always enters v_dot as a subtraction (opposes motion along the
    # velocity direction, consistent with M1's dv/dt = (T-D)/m - g sin(gamma)).
    params = make_params()
    y = np.array([R_EARTH + 5000.0, 0.0, 300.0, math.radians(60.0), V.m0])
    control = controls.vertical_thrust_control(V.thrust)
    dydt_with_drag = dyn.ascent_rhs(0.0, y, params, control)

    params_no_drag = make_params(cd=0.0)
    dydt_no_drag = dyn.ascent_rhs(0.0, y, params_no_drag, control)

    # Removing drag must increase (or leave unchanged) v_dot -- drag can only subtract.
    assert dydt_no_drag[2] >= dydt_with_drag[2]
    assert dydt_no_drag[2] > dydt_with_drag[2]  # strictly, since v_rel > 0 here


def test_zero_drag_limit_removes_drag_exactly():
    params = make_params(cd=0.0)
    v_rel = 500.0
    a = dyn.drag_acceleration(v_rel, h=10_000.0, m=V.m0, cd=params.drag_coefficient,
                               area=params.reference_area)
    assert a == 0.0


def test_mass_decreases_only_while_thrusting():
    params = make_params()
    y = np.array([R_EARTH + 1000.0, 0.0, 200.0, math.radians(45.0), V.m0])
    coast = controls.coast_control()
    dydt = dyn.ascent_rhs(0.0, y, params, coast)
    assert dydt[4] == 0.0  # m_dot == 0 with thrust off

    thrusting = controls.vertical_thrust_control(V.thrust)
    dydt2 = dyn.ascent_rhs(0.0, y, params, thrusting)
    assert dydt2[4] < 0.0  # m_dot < 0 while thrusting


def test_mass_flow_clamped_at_m_min_even_if_thrust_commanded():
    params = make_params()
    y = np.array([R_EARTH + 1000.0, 0.0, 200.0, math.radians(45.0), V.m_min])
    thrusting = controls.vertical_thrust_control(V.thrust)
    dydt = dyn.ascent_rhs(0.0, y, params, thrusting)
    assert dydt[4] == 0.0


def test_zero_angle_of_attack_reduces_to_m1_equations():
    # chi == gamma (thrust aligned with velocity) must reduce exactly to the M1 S3
    # equations: v_dot = (T-D)/m - g sin(gamma); v*gamma_dot = (v^2/r - g) cos(gamma).
    params = make_params()
    r, v, gamma, m = R_EARTH + 20_000.0, 400.0, math.radians(30.0), 300_000.0
    y = np.array([r, 0.0, v, gamma, m])

    def aligned_control(t, y, p):
        return V.thrust, y[3]  # chi = gamma

    dydt = dyn.ascent_rhs(0.0, y, params, aligned_control)

    g = dyn.local_gravity(r)
    v_rel = dyn.relative_speed(v, gamma, r, LAT)
    a_drag = dyn.drag_acceleration(v_rel, r - R_EARTH, m, params.drag_coefficient,
                                    params.reference_area)
    expected_v_dot = (V.thrust) / m - a_drag - g * math.sin(gamma)
    expected_gamma_dot = (v**2 / r - g) * math.cos(gamma) / v

    assert dydt[2] == pytest.approx(expected_v_dot)
    assert dydt[3] == pytest.approx(expected_gamma_dot)


def test_gamma_frozen_below_velocity_floor():
    params = make_params()
    y = np.array([R_EARTH, 0.0, 0.0, math.radians(90.0), V.m0])
    control = controls.vertical_thrust_control(V.thrust)
    dydt = dyn.ascent_rhs(0.0, y, params, control)
    assert dydt[3] == 0.0  # gamma frozen at v=0
