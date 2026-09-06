"""Thin integration wrapper around ``scipy.integrate.solve_ivp`` for the ascent ODE.

No guidance/optimization logic lives here -- this module only propagates whatever
``control`` function it is given (see ``controls.py``) and reports the result plus any
events.
"""

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy.integrate import solve_ivp

from .dynamics import (
    AscentParams,
    ControlFn,
    ascent_rhs,
    make_altitude_crossing_event,
    make_ground_impact_event,
    make_propellant_depletion_event,
)


@dataclass
class AscentResult:
    t: np.ndarray
    r: np.ndarray
    theta: np.ndarray
    v: np.ndarray
    gamma: np.ndarray
    m: np.ndarray
    t_events: list
    y_events: list
    success: bool
    message: str


def run_ascent(y0: np.ndarray, t_span, params: AscentParams, control: ControlFn,
                max_step: float = 0.5, rtol: float = 1e-9, atol: float = 1e-9,
                target_altitude: Optional[float] = None,
                dense_output: bool = False, t_eval: Optional[np.ndarray] = None):
    """Integrate the ascent ODE from y0 over t_span with the given control law.

    Events (DESIGN.md M2 "Events and termination"):
      - propellant depletion (non-terminal; logged for diagnostics)
      - ground impact, h = 0 while descending (terminal)
      - optional target-altitude crossing (non-terminal diagnostic), if
        ``target_altitude`` is given
    """
    events = [
        make_propellant_depletion_event(params),
        make_ground_impact_event(params),
    ]
    if target_altitude is not None:
        events.append(make_altitude_crossing_event(target_altitude, params))

    sol = solve_ivp(
        fun=lambda t, y: ascent_rhs(t, y, params, control),
        t_span=t_span,
        y0=y0,
        method="RK45",
        max_step=max_step,
        rtol=rtol,
        atol=atol,
        events=events,
        dense_output=dense_output,
        t_eval=t_eval,
    )

    r, theta, v, gamma, m = sol.y
    return AscentResult(
        t=sol.t, r=r, theta=theta, v=v, gamma=gamma, m=m,
        t_events=sol.t_events, y_events=sol.y_events,
        success=sol.success, message=sol.message,
    )
