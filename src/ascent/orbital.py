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


# ============================================================================
# Milestone 4: idealized apogee-circularization insertion criterion.
#
# Rationale (DESIGN.md M4 S2): requiring the RAW post-MECO state to already be a
# near-circular orbit (is_circular_orbit_achieved above) is the right check for "did
# the continuous main-engine burn alone insert into orbit," but a real ascent+
# circularization mission profile (ascent to a transfer-ellipse apogee near the target
# altitude, then a brief idealized "kick"/insertion impulse at apogee to raise perigee
# up to match) is the standard way payload-to-orbit capability is actually assessed.
# M4 uses the latter. The circularization burn is EXPLICITLY an idealized, impulsive,
# unconstrained-thrust insertion-stage burn -- it is never subtracted from the M4
# vehicle's own propellant/mass budget (finite-thrust upper-stage optimization is out
# of M4's scope), and it is reported as a diagnostic quantity, not silently treated as
# part of the continuous main-engine burn.
# ============================================================================

@dataclass(frozen=True)
class InsertionEvaluation:
    is_bound: bool
    apogee_altitude: Optional[float]        # m; None if unbound
    apogee_radius: Optional[float]          # m; None if unbound
    v_at_apogee: Optional[float]            # m/s; None if unbound (purely tangential there)
    circularization_delta_v: Optional[float]  # m/s; None if unbound
    reaches_bound_orbit: bool               # bound AND does not intersect Earth (M4 S2)
    meets_altitude_target: bool             # |apogee_altitude - target| <= altitude_tolerance
    circularization_within_cap: bool        # circularization_delta_v <= delta_v_cap
    meets_insertion_criterion: bool         # ALL of the above (M4's pass/fail)


def circularization_delta_v_at_apogee(state: OrbitalState) -> Optional[float]:
    """Idealized impulsive, tangential delta-v to circularize AT the orbit's own apogee.

    v_at_apogee = specific_angular_momentum / apogee_radius (purely tangential there,
    since radial velocity is exactly zero at apogee by definition). Returns None for an
    unbound (non-elliptical) state, where "apogee" is not defined.
    """
    if state.apogee_radius is None:
        return None
    v_apogee = state.specific_angular_momentum / state.apogee_radius
    v_circ_apogee = math.sqrt(MU_EARTH / state.apogee_radius)
    # For a valid ellipse, apogee speed is always <= local circular speed there.
    return max(v_circ_apogee - v_apogee, 0.0)


def evaluate_orbit_insertion(state: OrbitalState, target_altitude: float,
                              altitude_tolerance: float = 15_000.0,
                              delta_v_cap: float = 1500.0) -> InsertionEvaluation:
    """M4 insertion criterion (DESIGN.md M4 S2): ascent-to-apogee + idealized apogee
    circularization.

    PASS requires ALL of:
      1. the post-MECO osculating orbit is bound (state.is_bound),
      2. it does not itself already intersect the Earth as a standalone ellipse
         (``reaches_bound_orbit``) -- a sanity floor, not the headline criterion,
      3. the NATURAL apogee altitude (before circularizing) is within
         ``altitude_tolerance`` of ``target_altitude`` -- this is what an idealized,
         unconstrained circularization burn CANNOT fix: it can only change the
         orbit's shape at the altitude the ascent already reached, not move that
         altitude itself, and
      4. the idealized circularization delta-v needed at that apogee is no larger
         than ``delta_v_cap`` -- a documented, explicit engineering judgment call
         (not a silently smuggled-in constraint) bounding the insertion burn to a
         "plausible kick/insertion-stage impulse," not "a second full ascent burn."
         Because condition 3 already requires apogee near the target, condition 4
         mainly screens out cases with an extreme leftover eccentricity (e.g. a
         near-vertical lofted trajectory whose apogee happens to sit near 400 km but
         whose tangential speed there is nowhere close to orbital).

    ``altitude_tolerance`` defaults to 15 km (looser than M3's 10 km circular-orbit
    tolerance, since M4 explicitly circularizes rather than requiring the raw ascent
    to already be near-circular) and ``delta_v_cap`` defaults to 1500 m/s.
    """
    apogee_altitude = state.apogee_altitude
    apogee_radius = state.apogee_radius
    circ_dv = circularization_delta_v_at_apogee(state)
    v_apogee = (state.specific_angular_momentum / apogee_radius
                if apogee_radius is not None else None)

    reaches_bound_orbit = state.is_bound and not state.intersects_earth

    meets_altitude_target = (
        state.is_bound and apogee_altitude is not None
        and abs(apogee_altitude - target_altitude) <= altitude_tolerance
    )
    circularization_within_cap = circ_dv is not None and circ_dv <= delta_v_cap

    meets_insertion_criterion = (
        state.is_bound
        and reaches_bound_orbit
        and meets_altitude_target
        and circularization_within_cap
    )

    return InsertionEvaluation(
        is_bound=state.is_bound,
        apogee_altitude=apogee_altitude,
        apogee_radius=apogee_radius,
        v_at_apogee=v_apogee,
        circularization_delta_v=circ_dv,
        reaches_bound_orbit=reaches_bound_orbit,
        meets_altitude_target=meets_altitude_target,
        circularization_within_cap=circularization_within_cap,
        meets_insertion_criterion=meets_insertion_criterion,
    )


@dataclass(frozen=True)
class CutoffScanResult:
    found_passing_cutoff: bool
    cutoff_time: Optional[float]           # s; None if no passing cutoff was found
    cutoff_index: Optional[int]            # index into the scanned trajectory arrays
    state: Optional[OrbitalState]
    insertion: Optional[InsertionEvaluation]
    closest_apogee_miss_km: Optional[float]  # smallest |apogee-target| seen (m), even
                                              # among non-passing candidates -- for
                                              # honest "how close did it get" reporting
    closest_apogee_miss_index: Optional[int]  # index achieving that closest miss, even
                                               # when no candidate passed -- a much
                                               # better search starting-point than an
                                               # arbitrary fallback (e.g. t=0, where
                                               # gamma is trivially 0 for this vehicle's
                                               # initial condition)


def scan_best_insertion_cutoff(t, r, v, gamma, target_altitude: float,
                                altitude_tolerance: float = 15_000.0,
                                delta_v_cap: float = 1500.0) -> CutoffScanResult:
    """Scan a sampled powered-ascent trajectory for the best candidate engine-cutoff
    instant (DESIGN.md M4 S4): the sample time that MINIMIZES the idealized
    circularization delta-v among samples that (a) have gamma >= 0 (climbing, not
    already descending -- see DESIGN.md M4 S4 for why a descending cutoff state
    reliably fails the non-Earth-intersecting check even when its apogee number looks
    fine) and (b) meet the full M4 insertion criterion (``evaluate_orbit_insertion``).

    Returns ``found_passing_cutoff=False`` (with the closest apogee miss reported
    honestly) if no sample satisfies the criterion -- this function does not force a
    result to exist.
    """
    best = None
    best_apogee_miss = None
    best_apogee_miss_index = None

    for i in range(len(t)):
        if gamma[i] < 0.0:
            continue
        state = compute_orbital_state(r[i], v[i], gamma[i])
        if state.apogee_altitude is not None:
            miss = abs(state.apogee_altitude - target_altitude)
            if best_apogee_miss is None or miss < best_apogee_miss:
                best_apogee_miss = miss
                best_apogee_miss_index = i
        insertion = evaluate_orbit_insertion(state, target_altitude, altitude_tolerance,
                                              delta_v_cap)
        if insertion.meets_insertion_criterion:
            if best is None or insertion.circularization_delta_v < best[1].circularization_delta_v:
                best = (i, insertion, state)

    if best is None:
        return CutoffScanResult(
            found_passing_cutoff=False, cutoff_time=None, cutoff_index=None,
            state=None, insertion=None,
            closest_apogee_miss_km=(best_apogee_miss / 1000.0
                                     if best_apogee_miss is not None else None),
            closest_apogee_miss_index=best_apogee_miss_index,
        )

    idx, insertion, state = best
    return CutoffScanResult(
        found_passing_cutoff=True, cutoff_time=t[idx], cutoff_index=idx,
        state=state, insertion=insertion,
        closest_apogee_miss_km=(best_apogee_miss / 1000.0
                                 if best_apogee_miss is not None else None),
        closest_apogee_miss_index=best_apogee_miss_index,
    )
