"""M5: direct-ascent launch-azimuth / inclination geometry (DESIGN.md M5 S2).

Spherical-trigonometry relation (M1 S5, unchanged):

    cos(i) = cos(phi_launch) * sin(Az)

where ``i`` is the target orbital inclination, ``phi_launch`` is the launch-site
latitude, and ``Az`` is launch azimuth measured clockwise from true north (0=N,
90=E). This module provides small, independently-testable helpers around that
relation -- nothing here touches the ascent ODE or dynamics.

Direct-ascent accessibility: a site at latitude ``phi_launch`` can directly reach any
inclination ``i >= |phi_launch|`` (up to 90 deg here; retrograde is out of M5 scope).
Requesting ``i < |phi_launch|`` is not achievable by direct ascent from this site (it
would require a dogleg or an on-orbit plane change, neither modeled) and raises
``ValueError`` rather than silently clipping.

Rotational-velocity geometry: at the launch site, Earth's surface moves due east at
``v_rot = omega_earth * R_earth * cos(phi_launch)``. For a launch heading ``Az``
(measured clockwise from north), the component of that eastward velocity USEFUL along
the launch heading (i.e. the in-plane contribution an eastward-biased ascent keeps) is
the projection onto the heading unit vector:

    useful_boost(Az) = v_rot * sin(Az)

Substituting the inclination relation gives the compact closed form used throughout
M5:

    useful_boost(i) = omega_earth * R_earth * cos(i)

(independent of latitude once expressed in terms of i -- a standard, checkable
result). The REMAINING component of the site's rotational velocity, perpendicular to
the launch heading (i.e. out of the chosen launch plane), is:

    cross_track(Az) = v_rot * cos(Az) = omega_earth * R_earth * cos(phi_launch) * cos(Az)

This cross-track term is NOT captured by the planar ascent model's vehicle dynamics
(the vehicle is confined to its launch plane, by construction, for all M1-M5 direct-
ascent cases) but IS a real component of the atmosphere's true rotation -- see
DESIGN.md M5 S3 and ``dynamics.relative_speed``'s ``azimuth_rad`` parameter for how
this is retained (not silently dropped) in the drag calculation.
"""

import math

from .constants import OMEGA_EARTH, R_EARTH


def min_direct_ascent_inclination(launch_latitude_rad: float) -> float:
    """Minimum inclination directly reachable from this latitude (no dogleg)."""
    return abs(launch_latitude_rad)


def azimuth_from_inclination(inclination_rad: float, launch_latitude_rad: float) -> float:
    """Direct-ascent launch azimuth (rad, clockwise from north) for a target inclination.

    Returns the prograde-side (0 <= Az <= 90 deg) solution: Az=90 deg (due east) at
    i=latitude, decreasing toward Az=0 deg (due north) as i increases toward 90 deg.
    Raises ValueError for an inclination not directly reachable from this latitude
    (i < |latitude|) or i > 90 deg (retrograde/out of M5 scope -- see module docstring).
    """
    lat = launch_latitude_rad
    i_min = min_direct_ascent_inclination(lat)
    if inclination_rad < i_min - 1e-9:
        raise ValueError(
            f"Inclination {math.degrees(inclination_rad):.3f} deg is not directly "
            f"reachable from latitude {math.degrees(lat):.3f} deg (minimum direct-ascent "
            f"inclination is {math.degrees(i_min):.3f} deg). A dogleg or on-orbit plane "
            f"change would be required; neither is modeled in M1-M5."
        )
    if inclination_rad > math.pi / 2 + 1e-9:
        raise ValueError(
            f"Inclination {math.degrees(inclination_rad):.3f} deg is retrograde "
            f"(>90 deg), which is out of M5 scope by design (DESIGN.md M5 S15)."
        )
    ratio = math.cos(inclination_rad) / math.cos(lat)
    ratio = min(max(ratio, -1.0), 1.0)  # guard tiny floating-point overshoot at i=90
    return math.asin(ratio)


def inclination_from_azimuth(azimuth_rad: float, launch_latitude_rad: float) -> float:
    """Inverse of ``azimuth_from_inclination``: inclination reached by a given azimuth."""
    return math.acos(math.cos(launch_latitude_rad) * math.sin(azimuth_rad))


def site_rotational_speed(launch_latitude_rad: float, r: float = R_EARTH,
                           omega_earth: float = OMEGA_EARTH) -> float:
    """Total eastward co-rotation speed at radius r and this latitude, v_rot(r)."""
    return omega_earth * r * math.cos(launch_latitude_rad)


def useful_rotational_boost(azimuth_rad: float, launch_latitude_rad: float,
                             r: float = R_EARTH, omega_earth: float = OMEGA_EARTH) -> float:
    """In-launch-plane (useful) component of the site's rotational velocity.

    useful_boost = v_rot(r) * sin(Az) = omega_earth * r * cos(i)  (see module docstring
    for the equivalence, both forms given for cross-checking).
    """
    return site_rotational_speed(launch_latitude_rad, r, omega_earth) * math.sin(azimuth_rad)


def cross_track_rotational_component(azimuth_rad: float, launch_latitude_rad: float,
                                      r: float = R_EARTH,
                                      omega_earth: float = OMEGA_EARTH) -> float:
    """Out-of-launch-plane component of the site's rotational velocity (NOT modeled in
    the vehicle's planar dynamics, but retained in the atmosphere-relative-wind
    calculation -- see module docstring and DESIGN.md M5 S3).
    """
    return site_rotational_speed(launch_latitude_rad, r, omega_earth) * math.cos(azimuth_rad)
