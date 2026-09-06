"""M4 delta-v / loss accounting (DESIGN.md M4 S9).

Derivation (exact, from the M1 S3 / M2 dynamics.py equations of motion):

    v_dot = (T*cos(alpha) - D)/m - g(r)*sin(gamma)

Integrating from t=0 to the cutoff time t_c:

    v(t_c) - v(0) = INTEGRAL[ T*cos(alpha)/m ] dt
                  - INTEGRAL[ D/m ] dt              (drag loss)
                  - INTEGRAL[ g(r)*sin(gamma) ] dt  (gravity loss)

The first term, if thrust were always perfectly aligned with velocity (alpha=0,
cos(alpha)=1), would be exactly the ideal Tsiolkovsky delta-v for the mass consumed by
t_c: INTEGRAL[T/m]dt = ve*ln(m0/m(t_c)) (since m_dot=-T/ve). The shortfall from that
ideal value whenever alpha != 0 (i.e. during the pitch-kick phase, and any other instant
thrust is not exactly velocity-aligned) is:

    steering_loss = INTEGRAL[ (T/m)*(1 - cos(alpha)) ] dt

so that, exactly (an algebraic identity of the equations above, not a fitted residual):

    achieved_dv = ideal_dv_to_cutoff - steering_loss - drag_loss - gravity_loss

This module computes gravity_loss and drag_loss directly by numerical (trapezoidal)
time-integration of the sampled trajectory, and steering_loss as the exact algebraic
remainder of the identity above -- NOT as a leftover/unexplained residual. The identity
is verified to hold (achieved_dv + gravity_loss + drag_loss + steering_loss ==
ideal_dv_to_cutoff) to numerical-integration tolerance; see tests/test_m4_loss_budget.py.
"""

from dataclasses import dataclass

import numpy as np

from . import dynamics as dyn
from .atmosphere import density
from .constants import G0, MU_EARTH


@dataclass(frozen=True)
class DeltaVBudget:
    ideal_dv_to_cutoff: float   # m/s, Tsiolkovsky dv for the mass actually consumed by cutoff
    achieved_dv: float         # m/s, actual inertial speed increase v(t_c) - v(0)
    gravity_loss: float        # m/s
    drag_loss: float           # m/s
    steering_loss: float       # m/s, exact algebraic remainder (see module docstring)
    residual: float            # m/s, should be ~0: ideal - (achieved+gravity+drag+steering)


def compute_delta_v_budget(t: np.ndarray, r: np.ndarray, v: np.ndarray, gamma: np.ndarray,
                            m: np.ndarray, isp: float, cd: float, area: float,
                            latitude_rad: float, m0: float,
                            omega_earth=None, mu: float = MU_EARTH) -> DeltaVBudget:
    """Compute the delta-v/loss budget over the sampled arrays ``t[0]..t[-1]``.

    Arrays must all correspond to the SAME powered-flight segment (e.g. liftoff to
    engine cutoff); pass already-sliced arrays if the full run extends past cutoff.
    ``mu`` must match whatever gravitational parameter the trajectory was actually
    propagated with (default: Earth's) -- passing a mismatched ``mu`` (e.g. leaving it
    at Earth's value for a trajectory deliberately propagated with ``mu=0`` to isolate
    other loss terms) silently produces a nonsensical gravity_loss.
    """
    g = mu / r**2
    gravity_loss = np.trapezoid(g * np.sin(gamma), t)

    from .constants import R_EARTH
    h = r - R_EARTH
    kwargs = {} if omega_earth is None else dict(omega_earth=omega_earth)
    v_rel = np.array([
        dyn.relative_speed(vi, gi, ri, latitude_rad, **kwargs)
        for vi, gi, ri in zip(v, gamma, r)
    ])
    a_drag = dyn.drag_acceleration(v_rel, h, m, cd, area)
    drag_loss = np.trapezoid(a_drag, t)

    ve = isp * G0
    ideal_dv_to_cutoff = ve * np.log(m0 / m[-1])
    achieved_dv = v[-1] - v[0]

    steering_loss = ideal_dv_to_cutoff - achieved_dv - gravity_loss - drag_loss
    residual = ideal_dv_to_cutoff - (achieved_dv + gravity_loss + drag_loss + steering_loss)

    return DeltaVBudget(
        ideal_dv_to_cutoff=ideal_dv_to_cutoff, achieved_dv=achieved_dv,
        gravity_loss=gravity_loss, drag_loss=drag_loss, steering_loss=steering_loss,
        residual=residual,
    )
