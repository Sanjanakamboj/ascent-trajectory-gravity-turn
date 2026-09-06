import pytest

from ascent import propulsion as prop
from ascent.constants import BASELINE_VEHICLE, G0

V = BASELINE_VEHICLE


def test_m1_burn_duration_reproduced():
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    assert t_burn == pytest.approx(158.71288815789472, rel=1e-9)


def test_mass_flow_rate_matches_hand_calc():
    mdot = prop.mass_flow_rate(V.thrust, V.isp)
    assert mdot == pytest.approx(-2583.2810728774184, rel=1e-9)


def test_constant_thrust_mass_decreases_linearly():
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    times = [0.0, t_burn / 4, t_burn / 2, 3 * t_burn / 4, t_burn]
    mdot = prop.mass_flow_rate(V.thrust, V.isp)
    for t in times:
        expected = V.m0 + mdot * t
        assert prop.mass_at_time(V.m0, V.thrust, V.isp, t) == pytest.approx(expected)
    # linearity: equal time steps give equal mass decrements
    masses = [prop.mass_at_time(V.m0, V.thrust, V.isp, t) for t in times]
    diffs = [masses[i] - masses[i + 1] for i in range(len(masses) - 1)]
    for d in diffs:
        assert d == pytest.approx(diffs[0], rel=1e-9)


def test_integrated_mass_loss_equals_mdot_times_burn_time():
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    m_end = prop.mass_at_time(V.m0, V.thrust, V.isp, t_burn)
    assert (V.m0 - m_end) == pytest.approx(V.m_propellant, rel=1e-9)


def test_tsiolkovsky_matches_m1_hand_calc():
    mf = V.m0 - V.m_propellant
    dv = prop.tsiolkovsky_delta_v(V.isp, V.m0, mf)
    assert dv == pytest.approx(5044.928401454307, rel=1e-9)


def test_thrust_to_weight_at_liftoff_matches_m1():
    twr0 = prop.thrust_acceleration(V.thrust, V.m0) / G0
    assert twr0 == pytest.approx(1.5499686437264508, rel=1e-9)


def test_thrust_acceleration_increases_as_mass_falls():
    a1 = prop.thrust_acceleration(V.thrust, V.m0)
    a2 = prop.thrust_acceleration(V.thrust, V.m0 - 100_000)
    a3 = prop.thrust_acceleration(V.thrust, V.m_min)
    assert a1 < a2 < a3


def test_clamped_mass_flow_rate_stops_at_m_min():
    mdot_active = prop.clamped_mass_flow_rate(V.thrust, V.isp, V.m_min + 1.0, V.m_min)
    assert mdot_active < 0.0
    mdot_depleted = prop.clamped_mass_flow_rate(V.thrust, V.isp, V.m_min, V.m_min)
    assert mdot_depleted == 0.0
    mdot_below = prop.clamped_mass_flow_rate(V.thrust, V.isp, V.m_min - 1.0, V.m_min)
    assert mdot_below == 0.0
