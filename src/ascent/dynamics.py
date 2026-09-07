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

V_FLOOR = 1.0  # m/s, velocity below which gamma is frozen rather than integrated. Raised
# from 1e-3 to 1.0 m/s during M5 (DESIGN.md M5 S13): gamma_dot's magnitude scales like
# 1/v near this threshold, and M1-M4 never actually integrated THROUGH the v=1e-3
# crossing in the interior of a real trajectory (v0 was always hundreds of m/s there).
# M5's near-polar direct-ascent case does cross it (v0 -> 0 as inclination -> 90 deg),
# and the resulting near-singular gamma_dot ~ 1/v transient right at a 1e-3 m/s
# threshold produced severe, genuine numerical stiffness (a real trajectory needed
# >170,000 adaptive steps for just 30 s of flight). A larger, still numerically tiny,
# threshold (1 m/s is still >2 orders of magnitude below any ascent speed of interest)
# makes that same transient far more benign (1/v ~ 15 rather than ~15,000) without
# changing the frozen-gamma PHILOSOPHY at all, and resolves the stiffness in practice.


# A control function takes (t, y, params) and returns (thrust_N, chi_rad).
ControlFn = Callable[[float, np.ndarray, "AscentParams"], Tuple[float, float]]


@dataclass(frozen=True)
class AscentParams:
    """Explicit, immutable simulation parameters (no hidden globals).

    ``azimuth_rad`` (M5): launch azimuth, clockwise from north, of THIS trajectory's
    launch plane. Defaults to 90 deg (due east), matching every M1-M4 case exactly.
    See ``relative_speed`` for how a non-90-deg azimuth is used.
    """

    isp: float
    reference_area: float
    drag_coefficient: float
    m_min: float
    latitude_rad: float
    mu: float = MU_EARTH
    r_earth: float = R_EARTH
    omega_earth: float = OMEGA_EARTH
    g0: float = G0
    azimuth_rad: float = np.pi / 2


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
                    omega_earth: float = OMEGA_EARTH, azimuth_rad: float = np.pi / 2) -> float:
    """Speed relative to the rotating atmosphere, |v_vec - v_atm_vec|.

    Decomposes the inertial velocity into radial/tangential (in-launch-plane)
    components and subtracts the atmosphere's IN-PLANE co-rotation component from the
    tangential one, exactly as in M1-M4 (default ``azimuth_rad = 90 deg``, due east).

    M5 (DESIGN.md M5 S3): for a launch azimuth other than due east, Earth's true
    (eastward) atmospheric rotation also has a component CROSS the launch plane,
    ``v_atm_total * cos(azimuth_rad)`` (``inclination.cross_track_rotational_component``).
    The vehicle's planar dynamics carry no cross-track velocity of their own (a
    standard direct-ascent simplification -- the vehicle is confined to its launch
    plane for all of M1-M5), so this cross-track atmosphere motion is entirely a
    relative-wind contribution and is included here as a third (perpendicular)
    component of the relative-velocity vector, rather than silently dropped:

        v_rel = sqrt(v_radial^2 + (v_tangential - v_atm*sin(Az))^2 + (v_atm*cos(Az))^2)

    At ``azimuth_rad = 90 deg`` (due east, all of M1-M4), ``cos(Az) = 0`` and this is
    IDENTICAL to the M1-M4 formula -- exact regression, not an approximation of it.
    """
    v_radial = v * np.sin(gamma)
    v_tangential = v * np.cos(gamma)
    v_atm_total = atmosphere_corotation_speed(r, latitude_rad, omega_earth)
    v_atm_inplane = v_atm_total * np.sin(azimuth_rad)
    v_atm_cross = v_atm_total * np.cos(azimuth_rad)
    v_rel_tangential = v_tangential - v_atm_inplane
    return np.sqrt(v_radial**2 + v_rel_tangential**2 + v_atm_cross**2)


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
    v_rel = relative_speed(v, gamma, r, params.latitude_rad, params.omega_earth,
                            params.azimuth_rad)
    # Below V_FLOOR the vehicle's own inertial-velocity DIRECTION is undefined (the
    # same reason gamma_dot is frozen below), so a drag magnitude computed from v_rel
    # cannot be meaningfully assigned a direction "anti-parallel to inertial v" either
    # -- applying it anyway would silently subtract an undefined-direction force from
    # v_dot. This was latent but invisible in M1-M4 (v0 was always well above V_FLOOR,
    # the due-east co-rotation speed); M5 exposed it directly: at a near-polar azimuth
    # (small useful in-plane boost, v0 near 0) the atmosphere's cross-track rotation
    # component alone gives a nonzero v_rel even at v=0, which produced a spurious
    # negative v_dot at t=0 (drag decelerating an already-stationary vehicle) before
    # this guard was added. Zeroing drag below V_FLOOR is a direct, minimal extension
    # of the existing V_FLOOR philosophy, not a new one.
    a_drag = (drag_acceleration(v_rel, h, m, params.drag_coefficient, params.reference_area)
              if v >= V_FLOOR else 0.0)

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


GROUND_EVENT_DEADBAND = 1.0  # m; see make_ground_impact_event docstring


def make_ground_impact_event(params: AscentParams):
    """Terminal event: altitude crosses zero while descending (r decreasing through R_earth).

    A small (1 m) deadband is subtracted from the trigger altitude. Every M1-M5
    trajectory starts at r = r_earth EXACTLY; whenever gamma is also exactly 0 there
    (e.g. the M5 near-polar case, where the vehicle's in-plane velocity is ~0 and
    dynamics.py freezes gamma_dot), r_dot is analytically exactly 0 for a stretch of
    time -- a genuine, real "flat start," not a descent. A found, real failure mode
    (M5): with the event armed at EXACTLY r = r_earth, floating-point noise in the
    integrator's internal stages during that flat interval could dip the event
    function a hair below zero, and solve_ivp's direction=-1 root-finding would then
    report a (spurious) ground impact at t=0. A 1 m deadband (utterly negligible
    against any real descent, which moves the vehicle by meters to kilometers) absorbs
    that noise without weakening real impact detection.
    """

    def event(t, y):
        return y[0] - (params.r_earth - GROUND_EVENT_DEADBAND)

    event.terminal = True
    event.direction = -1
    return event


def make_altitude_crossing_event(target_altitude: float, params: AscentParams,
                                  terminal: bool = False):
    """Diagnostic (non-terminal by default) event: altitude crosses a target value.

    NOTE: crossing a target altitude is NOT the same as achieving orbit there -- see
    ``orbital.is_circular_orbit_achieved``, which additionally requires the perigee to
    also be near the target (ruling out a suborbital trajectory that merely passes
    through the target altitude on its way up or down).
    """

    def event(t, y):
        return (y[0] - params.r_earth) - target_altitude

    event.terminal = terminal
    event.direction = 0
    return event


def make_apogee_event(terminal: bool = False):
    """Apogee event: radial velocity (r_dot = v*sin(gamma)) crosses zero from + to -.

    Fires at a local maximum of altitude -- during coast this is the true orbital
    apogee; it can in principle also fire during powered flight if the guidance law
    ever drives r_dot through zero while thrusting, which is a legitimate use of the
    same mathematical condition (a momentary altitude peak), not a bug.
    """

    def event(t, y):
        r, theta, v, gamma, m = y
        return v * np.sin(gamma)

    event.terminal = terminal
    event.direction = -1  # crossing from positive (climbing) to negative (descending)
    return event


def make_atmosphere_exit_event(params: AscentParams, terminal: bool = False):
    """Optional diagnostic event: altitude crosses the atmosphere model's H_MAX cutoff.

    Purely a diagnostic marker of when this simplified atmosphere model's density goes
    to exact zero (DESIGN.md atmosphere.py docstring) -- not a physical "edge of space"
    claim.
    """
    return make_altitude_crossing_event(atmo.H_MAX, params, terminal=terminal)
