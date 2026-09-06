"""Prescribed thrust-direction controls for M2 verification.

These exist only to exercise the dynamics in ``dynamics.py`` -- none of them is a
guidance law, and none is tuned to reach orbit. M3 is where an actual (zero-AoA,
velocity-following) gravity-turn *guidance* law belongs; that law is deliberately
NOT implemented here.

Every control function has the signature required by ``dynamics.ascent_rhs``:

    control(t, y, params) -> (thrust_N, chi_rad)

where ``chi`` is the thrust-direction angle from local horizontal (same convention as
``gamma``, DESIGN.md S4).
"""

import math
from typing import Tuple

import numpy as np


def vertical_thrust_control(thrust: float):
    """Case A: constant thrust, held vertical (chi = 90 deg) for all time.

    Useful for isolating vertical-only dynamics (no curvature/turning terms driven by
    thrust misalignment).
    """

    def control(t, y, params) -> Tuple[float, float]:
        return thrust, math.pi / 2

    return control


def fixed_pitch_control(thrust: float, chi_rad: float):
    """Case B: constant thrust at a fixed, prescribed thrust-direction angle.

    ``chi_rad`` does not track gamma -- this is a genuinely prescribed (non-guided)
    pitch angle, held constant for the whole run.
    """

    def control(t, y, params) -> Tuple[float, float]:
        return thrust, chi_rad

    return control


def coast_control():
    """Case C: thrust off (T = 0). Direction is irrelevant when thrust is zero."""

    def control(t, y, params) -> Tuple[float, float]:
        return 0.0, 0.0

    return control


def vertical_rise_then_pitch_kick_control(thrust: float, kick_start_time: float,
                                           kick_angle_rad: float, chi_after_kick: float):
    """Diagnostic-only M1-style demonstration: vertical rise, then a fixed pitch-kick.

    NOT a gravity-turn guidance law and NOT optimized for orbit (DESIGN.md M1 S1.7
    describes the pitch-kick concept; this is a minimal prescribed realization of it
    for exercising the dynamics, per M2 scope).

    - For t < kick_start_time: thrust vertical (chi = 90 deg).
    - For t >= kick_start_time: thrust held at a single fixed angle
      ``chi_after_kick`` (e.g. 90 deg - kick_angle_rad), constant thereafter. This is a
      one-time kick to a new *fixed* pitch, not a continuously-updated guidance law.
    """

    def control(t, y, params) -> Tuple[float, float]:
        if t < kick_start_time:
            return thrust, math.pi / 2
        return thrust, chi_after_kick

    return control
