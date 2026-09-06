"""Planar point-mass ascent dynamics (DESIGN.md S2-S6).

State vector (exact M1 convention, DESIGN.md S2):

    y = [r, theta, v, gamma, m]

    r      -- radial distance from Earth's center [m]
    theta  -- downrange central angle in the inertial orbital plane [rad]
    v      -- INERTIAL speed (ECI) [m/s]
    gamma  -- flight-path angle, angle between velocity and local horizontal,
              positive climbing [rad]
    m      -- instantaneous vehicle mass [kg]

Thrust direction / control interface
-------------------------------------
Physics is kept separate from guidance (per M2 scope: no gravity-turn law is hard-wired
here). A control input supplies, at each instant, the thrust magnitude and the thrust
*direction* expressed as an angle ``chi`` measured the same way as ``gamma`` (from local
horizontal, positive up). The angle between thrust and velocity directions is

    alpha = chi - gamma        (thrust angle of attack)

and the general (non-zero-alpha) planar equations of motion are:

    r_dot     = v * sin(gamma)
    theta_dot = v * cos(gamma) / r
    v_dot     = (T * cos(alpha) - D) / m  -  g(r) * sin(gamma)
    gamma_dot = [ T * sin(alpha) / m  +  (v**2 / r - g(r)) * cos(gamma) ] / v
    m_dot     = clamped mass-flow rate (propulsion.clamped_mass_flow_rate)

When the control law commands ``chi == gamma`` (zero angle of attack, thrust aligned
with velocity -- the M1 gravity-turn assumption), ``alpha = 0`` and this reduces exactly
to the M1 S3 equations:

    v_dot     = (T - D) / m - g(r) * sin(gamma)
    v*gamma_dot = (v**2 / r - g(r)) * cos(gamma)

Numerical handling of v -> 0 (liftoff singularity)
----------------------------------------------------
``gamma_dot`` above is obtained by dividing by ``v``, which is singular at liftoff
(v = 0, direction undefined). This is a well known feature of gravity-turn formulations,
not a bug: below a small velocity floor ``V_FLOOR`` the flight-path angle is *frozen*
(``gamma_dot = 0``) rather than integrated, matching the standard practice of holding the
vehicle vertical (or at its current attitude) until it has picked up enough speed for the
velocity direction to be well defined. This is an explicit, documented approximation.

Atmosphere-relative (relative-wind) velocity and its planar approximation
---------------------------------------------------------------------------
Drag must use the vehicle's speed relative to the *rotating atmosphere*, not the raw
inertial speed (DESIGN.md S3 "Note on relative wind"). The atmosphere is assumed to
co-rotate rigidly with the Earth. For the due-east baseline case (inclination = launch
latitude), the orbital plane's in-plane "horizontal" direction locally coincides with the
direction of Earth's rotation at the launch site, so -- within this planar model -- the
atmosphere's co-rotation velocity is approximated as purely tangential (horizontal, i.e.
along the same direction as the ``theta``/horizontal velocity component), with magnitude

    v_atm(r) = omega_earth * r * cos(lat)

which reduces to the M1 surface value ``omega_earth * R_earth * cos(lat)`` at ``r =
R_earth`` (h = 0). This is an explicit planar approximation: it is exact only for a
purely equatorial-tangent plane; for the due-east baseline it is the natural and
standard simplification, and is documented here rather than silently assumed.

Drag magnitude uses this relative speed, ``v_rel = |v_vec - v_atm_vec|`` (vector
subtraction of the two horizontal/radial components). Consistent with the M1 S3
equations -- which place the drag term only in the ``v_dot`` equation, not in the
``gamma_dot`` equation -- drag's *direction* is taken anti-parallel to the inertial
velocity vector in these two ODEs; only its *magnitude* uses the relative-wind speed.
A fully vector-resolved treatment (drag anti-parallel to v_rel, with components entering
both equations) is a possible future refinement, not implemented in M2 to avoid
introducing structure beyond the M1-documented equations.
"""

from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np

from .constants import G0, MU_EARTH, OMEGA_EARTH, R_EARTH
from . import atmosphere as atmo
from .propulsion import clamped_mass_flow_rate

V_FLOOR = 1e-3  # m/s, velocity below which gamma is frozen rather than integrated

# A control function takes (t, y, params) and returns (thrust_N, chi_rad).
ControlFn = Callable[[float, np.ndarray, "AscentParams"], Tuple[float, float]]


@dataclass(frozen=True)
class AscentParams:
    """Explicit, immutable simulation parameters (no hidden globals)."""

    isp: float
    reference_area: float
    drag_coefficient: float
    m_min: float
    latitude_rad: float
    mu: float = MU_EARTH
    r_earth: float = R_EARTH
    omega_earth: float = OMEGA_EARTH
    g0: float = G0


def local_gravity(r: float, mu: float = MU_EARTH) -> float:
    """Spherical gravitational acceleration g(r) = mu / r**2, magnitude [m/s^2]."""
    return mu / r**2


def atmosphere_corotation_speed(r: float, latitude_rad: float,
                                 omega_earth: float = OMEGA_EARTH) -> float:
    """Tangential (horizontal) speed of the rigidly-co-rotating atmosphere at radius r.

    v_atm(r) = omega_earth * r * cos(lat)

    See module docstring for the planar-approximation caveat. Reduces to the M1 surface
    rotational speed at r = R_earth.
    """
    return omega_earth * r * np.cos(latitude_rad)


def relative_speed(v: float, gamma: float, r: float, latitude_rad: float,
                    omega_earth: float = OMEGA_EARTH) -> float:
    """Speed relative to the rotating atmosphere, |v_vec - v_atm_vec|.

    Decomposes the inertial velocity into radial/tangential components, subtracts the
    (purely tangential) atmospheric co-rotation velocity, and returns the magnitude of
    the result.
    """
    v_radial = v * np.sin(gamma)
    v_tangential = v * np.cos(gamma)
    v_atm = atmosphere_corotation_speed(r, latitude_rad, omega_earth)
    v_rel_tangential = v_tangential - v_atm
    return np.sqrt(v_radial**2 + v_rel_tangential**2)


def drag_acceleration(v_rel: float, h: float, m: float, cd: float, area: float) -> float:
    """Drag deceleration magnitude D/m [m/s^2], using atmosphere-relative speed."""
    rho = atmo.density(h)
    drag_force = 0.5 * rho * v_rel**2 * cd * area
    return drag_force / m


def ascent_rhs(t: float, y: np.ndarray, params: AscentParams, control: ControlFn) -> np.ndarray:
    """Right-hand side dy/dt for y = [r, theta, v, gamma, m].

    ``control(t, y, params)`` must return ``(thrust_N, chi_rad)``.
    """
    r, theta, v, gamma, m = y
    h = r - params.r_earth

    thrust, chi = control(t, y, params)

    # Mass never drops below m_min; thrust/mass-flow are both zeroed once depleted.
    if m <= params.m_min:
        thrust = 0.0
    mdot = clamped_mass_flow_rate(thrust, params.isp, m, params.m_min, params.g0)

    g = local_gravity(r, params.mu)
    v_rel = relative_speed(v, gamma, r, params.latitude_rad, params.omega_earth)
    a_drag = drag_acceleration(v_rel, h, m, params.drag_coefficient, params.reference_area)

    alpha = chi - gamma

    r_dot = v * np.sin(gamma)
    theta_dot = v * np.cos(gamma) / r
    v_dot = (thrust * np.cos(alpha)) / m - a_drag - g * np.sin(gamma)

    if v < V_FLOOR:
        gamma_dot = 0.0
    else:
        gamma_dot = (thrust * np.sin(alpha) / m + (v**2 / r - g) * np.cos(gamma)) / v

    return np.array([r_dot, theta_dot, v_dot, gamma_dot, mdot])


# --------------------------------------------------------------------------------------
# Events (DESIGN.md M2 section: event handling)
# --------------------------------------------------------------------------------------

def make_propellant_depletion_event(params: AscentParams):
    """Non-terminal event: mass crosses m_min from above (propellant just depleted)."""

    def event(t, y):
        return y[4] - params.m_min

    event.terminal = False
    event.direction = -1  # crossing from above (mass decreasing) to m_min
    return event


def make_ground_impact_event(params: AscentParams):
    """Terminal event: altitude crosses zero while descending (r decreasing through R_earth)."""

    def event(t, y):
        return y[0] - params.r_earth

    event.terminal = True
    event.direction = -1
    return event


def make_altitude_crossing_event(target_altitude: float, params: AscentParams,
                                  terminal: bool = False):
    """Diagnostic (non-terminal by default) event: altitude crosses a target value."""

    def event(t, y):
        return (y[0] - params.r_earth) - target_altitude

    event.terminal = terminal
    event.direction = 0
    return event
