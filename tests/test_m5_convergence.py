"""M5 integrator convergence (check H, DESIGN.md M5 S13): 28.5 deg, one mid
inclination, and 90 deg, each at 3 solver settings.

Guidance/cutoff values below are taken directly from the committed authoritative
sweep CSV (scripts/m5_inclination_sweep_results.csv) -- NOT placeholders -- so this
test always exercises a real, previously-solved trajectory rather than guessing
parameters that might not converge or might hit an edge case (see DESIGN.md M5 S13:
an earlier draft of this test used an approximate/placeholder guidance+cutoff for the
90 deg case and it hung -- tracked down to exactly the near-polar v~0 edge case M5
had to fix in dynamics.py/controls.py; using the real solved values avoids the same
risk here and is more representative besides).
"""

import csv
import math
import os

import numpy as np
import pytest

from ascent import controls, dynamics as dyn, inclination as inc
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, m4_vehicle
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "scripts",
                         "m5_inclination_sweep_results.csv")


def _load_case(i_deg):
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            if float(row["inclination_deg"]) == i_deg:
                return row
    raise AssertionError(f"inclination {i_deg} not found in {CSV_PATH}")


def _run(i_deg, payload=None):
    row = _load_case(i_deg)
    payload = float(row["max_payload_kg"]) if payload is None else payload
    guidance = (float(row["kick_start_s"]), float(row["kick_angle_deg"]),
                float(row["kick_duration_s"]))
    cutoff_time = float(row["cutoff_time_s"])

    az = inc.azimuth_from_inclination(math.radians(i_deg), LAT)
    vehicle = m4_vehicle(payload)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=vehicle.drag_coefficient,
                               m_min=vehicle.m_min, latitude_rad=LAT, azimuth_rad=az)
    v0 = inc.useful_rotational_boost(az, LAT)
    y0 = np.array([R_EARTH, 0.0, v0, 0.0, vehicle.m0])
    control = controls.gravity_turn_control(vehicle.thrust, guidance[0],
                                             math.radians(guidance[1]), guidance[2])

    # NOTE: the tightest tier deliberately stops at 1e-9/1e-9, not 1e-12 (as M2-M4 used):
    # the near-polar case's early flight (small useful rotational boost -> small initial
    # v, near the dynamics.V_FLOOR transition -- DESIGN.md M5 S13) is genuinely stiff
    # there, and rtol=1e-10 or tighter was found to make a single integration take
    # minutes instead of seconds. 1e-9 is still a very tight tolerance and the
    # comparison below still verifies real convergence, just without going beyond the
    # point where cost becomes impractical for this specific regime.
    settings = [
        dict(max_step=2.0, rtol=1e-6, atol=1e-6),
        dict(max_step=0.5, rtol=1e-8, atol=1e-8),
        dict(max_step=0.1, rtol=1e-9, atol=1e-9),
    ]
    finals = []
    for s in settings:
        result = run_ascent(y0, (0.0, cutoff_time), params, control,
                             terminal_depletion=False, include_apogee_event=False, **s)
        finals.append((result.r[-1], result.v[-1], result.gamma[-1], result.m[-1]))
    return np.array(finals)


@pytest.mark.parametrize("i_deg", [28.5, 55.0, 90.0])
def test_convergence_across_inclinations(i_deg):
    finals = _run(i_deg)
    diffs = np.abs(np.diff(finals, axis=0))
    # Tightest two settings must agree closely regardless of inclination.
    assert diffs[1, 0] < 10.0       # r within 10 m
    assert diffs[1, 1] < 0.1        # v within 0.1 m/s
    assert diffs[1, 2] < 1e-5       # gamma within ~6e-4 deg
    assert diffs[1, 3] < 1e-2       # m within 10 g
