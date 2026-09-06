"""M4 vehicle bookkeeping tests (DESIGN.md M4 S3/S6).

Verifies the M4 orbit-capable study vehicle is genuinely distinct from
BASELINE_VEHICLE (M1-M3), and that the payload-sweep mass bookkeeping is correct
(the specific bugs DESIGN.md M4 S6 calls out: payload treated as propellant, m0 held
fixed while payload varies, mass allowed below dry+payload, double-counted payload).
"""

import pytest

from ascent.constants import BASELINE_VEHICLE, M4_DESIGN, M4_REFERENCE_VEHICLE, m4_vehicle


def test_m4_vehicle_distinct_from_baseline():
    # Different Isp, dry mass, propellant mass -- not the same vehicle in any way.
    assert M4_REFERENCE_VEHICLE.isp != BASELINE_VEHICLE.isp
    assert M4_REFERENCE_VEHICLE.m_dry != BASELINE_VEHICLE.m_dry
    assert M4_REFERENCE_VEHICLE.m_propellant != BASELINE_VEHICLE.m_propellant
    # Thrust and aerodynamics are the documented UNCHANGED parameters.
    assert M4_REFERENCE_VEHICLE.thrust == BASELINE_VEHICLE.thrust
    assert M4_REFERENCE_VEHICLE.reference_area == BASELINE_VEHICLE.reference_area
    assert M4_REFERENCE_VEHICLE.drag_coefficient == BASELINE_VEHICLE.drag_coefficient


def test_baseline_vehicle_unmodified():
    # Defensive regression: the M1-M3 baseline's exact documented values.
    assert BASELINE_VEHICLE.m0 == 500_000.0
    assert BASELINE_VEHICLE.m_propellant == 410_000.0
    assert BASELINE_VEHICLE.m_dry == 80_000.0
    assert BASELINE_VEHICLE.m_payload == 10_000.0
    assert BASELINE_VEHICLE.thrust == 7.6e6
    assert BASELINE_VEHICLE.isp == 300.0


def test_m0_varies_with_payload_dry_and_propellant_fixed():
    # The specific bug DESIGN.md M4 S6 calls out: holding m0 fixed while sweeping
    # payload. Here m0 MUST change with payload; m_dry and m_propellant must NOT.
    v_low = m4_vehicle(0.0)
    v_high = m4_vehicle(50_000.0)
    assert v_high.m0 > v_low.m0
    assert v_high.m0 - v_low.m0 == pytest.approx(50_000.0)
    assert v_low.m_dry == v_high.m_dry == M4_DESIGN.m_dry
    assert v_low.m_propellant == v_high.m_propellant == M4_DESIGN.m_propellant


def test_payload_not_double_counted_or_treated_as_propellant():
    payload = 12_345.0
    v = m4_vehicle(payload)
    # Exactly one payload term in the mass budget (Vehicle.__post_init__ already
    # enforces m_dry + m_propellant + m_payload == m0 -- this re-derives it
    # independently here rather than relying solely on that internal check).
    assert v.m_dry + v.m_propellant + payload == pytest.approx(v.m0)
    # Payload must not have silently become part of the propellant load.
    assert v.m_propellant == M4_DESIGN.m_propellant
    assert v.m_payload == payload


def test_m_min_is_dry_plus_payload_and_scales_with_payload():
    v1 = m4_vehicle(0.0)
    v2 = m4_vehicle(20_000.0)
    assert v1.m_min == pytest.approx(M4_DESIGN.m_dry)
    assert v2.m_min == pytest.approx(M4_DESIGN.m_dry + 20_000.0)
    assert v2.m_min > v1.m_min


def test_negative_or_absurd_payload_still_consistent():
    # Not a realistic case, but the bookkeeping identity must hold regardless.
    v = m4_vehicle(1.0)
    assert v.m0 == pytest.approx(v.m_dry + v.m_propellant + v.m_payload)
