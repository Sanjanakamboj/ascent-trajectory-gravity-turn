"""Orbital-state diagnostics derived from the planar point-mass ascent state.

Pure, independently-derived post-processing: nothing here feeds back into the dynamics
or guidance. Given an instantaneous inertial state ``(r, v, gamma)`` these functions
compute the standard two-body orbital elements the state implies *if propagated
unpowered/undragged from that instant* (osculating elements) -- they say nothing about
what actually happens next if thrust or drag continue to act.

Conventions match DESIGN.md S2/S4: ``v`` is inertial speed, ``gamma`` is flight-path
angle from local horizontal (positive climbing), ``r`` is radial distance from Earth's
center.
"""

import math
from dataclasses import dataclass
from typing import Optional

from .constants import MU_EARTH, R_EARTH


@dataclass(frozen=True)
class OrbitalState:
    r: float                       # m, radial distance from Earth's center
    altitude: float                # m, r - R_earth
    v: float                       # m/s, inertial speed
    gamma: float                   # rad, flight-path angle
    v_radial: float                # m/s
    v_tangential: float            # m/s
    specific_energy: float         # J/kg (m^2/s^2)
    specific_angular_momentum: float  # m^2/s
    eccentricity: float
    semi_major_axis: Optional[float]   # m; None if not a bound (elliptical) orbit
    perigee_radius: float          # m; defined for any conic
    perigee_altitude: float        # m
    apogee_radius: Optional[float]     # m; None if unbound (e >= 1)
    apogee_altitude: Optional[float]   # m; None if unbound
    circular_speed_here: float     # m/s, sqrt(mu/r) -- speed for a circular orbit at r
    is_bound: bool                 # specific_energy < 0
    intersects_earth: bool         # perigee_radius < r_earth


def velocity_components(v: float, gamma: float) -> tuple:
    """Decompose inertial speed/flight-path angle into (v_radial, v_tangential)."""
    return v * math.sin(gamma), v * math.cos(gamma)


def speed_and_gamma_from_components(v_radial: float, v_tangential: float) -> tuple:
    """Inverse of ``velocity_components``: reconstruct (v, gamma) from components."""
    v = math.hypot(v_radial, v_tangential)
    gamma = math.atan2(v_radial, v_tangential)
    return v, gamma


def compute_orbital_state(r: float, v: float, gamma: float,
                           mu: float = MU_EARTH, r_earth: float = R_EARTH) -> OrbitalState:
    """Compute the full set of osculating orbital diagnostics for one instantaneous state."""
    v_radial, v_tangential = velocity_components(v, gamma)

    specific_energy = 0.5 * v**2 - mu / r
    h_ang = r * v_tangential  # specific angular momentum, r * v * cos(gamma)

    p = h_ang**2 / mu  # semi-latus rectum, defined for any conic (h_ang could be 0)

    # e^2 = 1 + 2*energy*h^2/mu^2 ; clip tiny negative values from floating-point noise.
    e_sq = 1.0 + 2.0 * specific_energy * h_ang**2 / mu**2
    eccentricity = math.sqrt(max(e_sq, 0.0))

    is_bound = specific_energy < 0.0

    if abs(1.0 - eccentricity) < 1e-12:
        # Exactly circular relative to floating-point noise: p = r, use r directly to
        # avoid a 0/0-adjacent division below.
        perigee_radius = r
    else:
        perigee_radius = p / (1.0 + eccentricity)
    perigee_altitude = perigee_radius - r_earth

    if eccentricity < 1.0:
        semi_major_axis = -mu / (2.0 * specific_energy)
        apogee_radius = p / (1.0 - eccentricity) if eccentricity < 1.0 else None
        apogee_altitude = apogee_radius - r_earth if apogee_radius is not None else None
    else:
        # Parabolic/hyperbolic: unbound, no apogee, semi-major axis not a bound-orbit size.
        semi_major_axis = None
        apogee_radius = None
        apogee_altitude = None

    circular_speed_here = math.sqrt(mu / r)

    return OrbitalState(
        r=r, altitude=r - r_earth, v=v, gamma=gamma,
        v_radial=v_radial, v_tangential=v_tangential,
        specific_energy=specific_energy, specific_angular_momentum=h_ang,
        eccentricity=eccentricity, semi_major_axis=semi_major_axis,
        perigee_radius=perigee_radius, perigee_altitude=perigee_altitude,
        apogee_radius=apogee_radius, apogee_altitude=apogee_altitude,
        circular_speed_here=circular_speed_here, is_bound=is_bound,
        intersects_earth=perigee_radius < r_earth,
    )


def is_circular_orbit_achieved(state: OrbitalState, target_altitude: float,
                                altitude_tolerance: float = 10_000.0) -> bool:
    """Explicit orbit-achievement criterion (documented in DESIGN.md M3 section).

    Requires ALL of:
      1. the orbit is bound (specific_energy < 0),
      2. it does not intersect the Earth (perigee_radius > r_earth), and
      3. BOTH perigee and apogee altitude fall within ``altitude_tolerance`` of
         ``target_altitude`` -- i.e. a near-circular orbit AT the target altitude, not
         merely a trajectory that happened to cross that altitude on its way up.

    This deliberately rules out "max altitude > target" as a false positive: a suborbital
    lofted trajectory can have an apogee far above the target while its perigee is deep
    inside the Earth (or undefined/intersecting), which fails condition 2/3 here.
    """
    if not state.is_bound or state.intersects_earth:
        return False
    if state.apogee_altitude is None:
        return False
    return (
        abs(state.perigee_altitude - target_altitude) <= altitude_tolerance
        and abs(state.apogee_altitude - target_altitude) <= altitude_tolerance
    )
