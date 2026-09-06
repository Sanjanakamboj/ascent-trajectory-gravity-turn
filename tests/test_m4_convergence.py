"""M4 integrator convergence (check H, DESIGN.md M4 S8/S11).

Repeats the near-boundary/reference payload powered-ascent-to-cutoff propagation at
three solver settings and confirms key outputs converge.
"""

import math

import numpy as np

from ascent import controls, dynamics as dyn, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, m4_vehicle
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
KICK_START = 45.0
KICK_ANGLE_DEG = 20.0
KICK_DURATION = 20.0
CUTOFF_TIME = 251.2047364808033  # from scripts/m4_payload_sweep_results.csv, payload=10,333 kg


def test_integrator_convergence_near_boundary_payload():
    vehicle = m4_vehicle(10_333.0)  # the maximum-passing payload (DESIGN.md M4 S7)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=vehicle.drag_coefficient,
                               m_min=vehicle.m_min, latitude_rad=LAT)
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    y0 = np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])
    control = controls.gravity_turn_control(vehicle.thrust, KICK_START,
                                             math.radians(KICK_ANGLE_DEG), KICK_DURATION)

    settings = [
        dict(max_step=2.0, rtol=1e-6, atol=1e-6),
        dict(max_step=0.3, rtol=1e-9, atol=1e-9),
        dict(max_step=0.05, rtol=1e-12, atol=1e-12),
    ]
    finals = []
    for s in settings:
        result = run_ascent(y0, (0.0, CUTOFF_TIME), params, control,
                             terminal_depletion=False, include_apogee_event=False, **s)
        finals.append((result.r[-1], result.v[-1], result.gamma[-1], result.m[-1]))

    finals = np.array(finals)
    diffs = np.abs(np.diff(finals, axis=0))
    # Tightest two settings must agree closely.
    assert diffs[1, 0] < 1.0        # r within 1 m
    assert diffs[1, 1] < 1e-2       # v within 1 cm/s
    assert diffs[1, 2] < 1e-6       # gamma within ~6e-5 deg
    assert diffs[1, 3] < 1e-3       # m within 1 g
