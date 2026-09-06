"""M4: drag-consistent engine-cutoff / insertion search.

``orbital.scan_best_insertion_cutoff`` evaluates the osculating (idealized, two-body,
drag-free) orbital elements AT the candidate cutoff instant. That is exactly correct
for a cutoff altitude comfortably above the atmosphere model's ``H_MAX`` (100 km,
``atmosphere.py``), where density is already exactly zero and the true trajectory really
is Keplerian from that point on. But this vehicle's guidance reaches its best insertion
window at a fairly LOW cutoff altitude (order 70-110 km), still within the residual
(if small) modeled atmosphere -- so a real, drag-included coast from there up through
100 km loses some additional energy that the idealized instantaneous estimate does not
capture. On the reference-payload case this was a real, checked discrepancy: the
idealized cutoff-instant apogee was 386 km, but the drag-included coast up to 100 km
altitude gave 293 km -- large enough to matter, not a rounding effect.

This module fixes that with a two-stage search:
  1. A CHEAP coarse scan (no ODE re-integration -- polynomial evaluation of the
     powered-phase's dense/continuous solution) of the IDEALIZED apogee-miss across the
     whole burn, to bracket a promising neighborhood.
  2. A bounded scalar optimization (``scipy.optimize.minimize_scalar``) within that
     neighborhood using the actual DRAG-INCLUDED coast-to-vacuum objective -- this
     avoids the grid-resolution artifacts a fixed-step scan is prone to near a narrow
     optimum (found during development: two different, both "reasonable", scan
     resolutions disagreed on whether the same (kick_start, kick_angle, kick_duration)
     passed at all).
"""

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.optimize import minimize_scalar

from . import controls, orbital
from .atmosphere import H_MAX
from .constants import R_EARTH
from .dynamics import AscentParams
from .simulation import run_ascent


@dataclass(frozen=True)
class DragConsistentCutoff:
    found_passing_cutoff: bool
    cutoff_time: Optional[float]
    cutoff_state: Optional[np.ndarray]        # y at the engine-cutoff instant itself
    post_atmosphere_state: Optional[orbital.OrbitalState]  # state once h > H_MAX
    insertion: Optional[orbital.InsertionEvaluation]
    closest_apogee_miss_km: Optional[float]


def _coast_to_vacuum(y_cutoff: np.ndarray, params: AscentParams, t_cap: float = 1000.0):
    """Coast (thrust off, drag on) from y_cutoff until altitude exceeds H_MAX (or t_cap)."""
    coast = controls.coast_control()
    result = run_ascent(y_cutoff, (0.0, t_cap), params, coast, max_step=0.5,
                         rtol=1e-10, atol=1e-10, include_apogee_event=False,
                         target_altitude=H_MAX)
    h = result.r - R_EARTH
    idx = int(np.argmax(h >= H_MAX)) if np.any(h >= H_MAX) else len(h) - 1
    return result.r[idx], result.v[idx], result.gamma[idx]


def find_drag_consistent_cutoff(vehicle, params: AscentParams, y0: np.ndarray,
                                 kick_start: float, kick_angle_deg: float,
                                 kick_duration: float, t_burn_full: float,
                                 target_altitude: float, n_coarse: int = 800,
                                 bracket_half_width: float = 8.0) -> DragConsistentCutoff:
    """Find the best engine-cutoff time using the actual drag-included coast-to-vacuum
    state (see module docstring), via a cheap coarse bracket + expensive local refine.
    """
    control = controls.gravity_turn_control(
        vehicle.thrust, kick_start_time=kick_start,
        kick_angle_rad=math.radians(kick_angle_deg), kick_duration=kick_duration,
    )
    full_result = run_ascent(y0, (0.0, t_burn_full), params, control, max_step=0.5,
                              rtol=1e-10, atol=1e-10, terminal_depletion=False,
                              include_apogee_event=False, dense_output=True)
    sol = full_result.sol
    # The integration may have terminated before t_burn_full (e.g. a premature ground
    # impact for a bad guidance choice) -- the dense solution is only valid up to
    # whatever time was actually reached, NOT the nominal full-burn duration. Scanning
    # or evaluating past that produces extrapolation garbage (this was a real bug found
    # during development: a negative-r ValueError from evaluating sol(t) far outside
    # its valid range).
    t_reached = full_result.t[-1]
    if t_reached < 1.0:
        return DragConsistentCutoff(False, None, None, None, None, None)

    # Stage 1: cheap coarse scan of the IDEALIZED apogee-miss via the dense solution
    # (no re-integration -- just polynomial evaluation), to find a good neighborhood.
    coarse_ts = np.linspace(0.0, t_reached, n_coarse)
    coarse_states = sol(coarse_ts)  # shape (5, n_coarse)
    r_c, theta_c, v_c, gamma_c, m_c = coarse_states

    best_idx = None
    best_miss = None
    for i in range(n_coarse):
        if gamma_c[i] < 0.0:
            continue
        state = orbital.compute_orbital_state(r_c[i], v_c[i], gamma_c[i])
        if state.apogee_altitude is None:
            continue
        miss = abs(state.apogee_altitude - target_altitude)
        if best_miss is None or miss < best_miss:
            best_miss = miss
            best_idx = i

    if best_idx is None:
        return DragConsistentCutoff(False, None, None, None, None, None)

    t_guess = coarse_ts[best_idx]
    lo = max(t_guess - bracket_half_width, 0.0)
    hi = min(t_guess + bracket_half_width, t_reached)

    # Stage 2: expensive, drag-included objective, refined with a bounded scalar
    # optimizer (robust to the grid-resolution fragility a fixed-step scan showed).
    def objective(t_c):
        y_c = sol(t_c)
        if y_c[3] < 0.0:  # descending -- not a viable cutoff; penalize smoothly
            return 1e9 + abs(y_c[3]) * 1e6
        r_post, v_post, gamma_post = _coast_to_vacuum(y_c, params)
        state = orbital.compute_orbital_state(r_post, v_post, gamma_post)
        if state.apogee_altitude is None:
            return 1e9
        return abs(state.apogee_altitude - target_altitude)

    opt = minimize_scalar(objective, bounds=(lo, hi), method="bounded",
                           options=dict(xatol=1e-3))
    t_best = float(opt.x)

    y_best = sol(t_best)
    if y_best[3] < 0.0:
        return DragConsistentCutoff(False, None, None, None, None, best_miss / 1000.0)

    r_post, v_post, gamma_post = _coast_to_vacuum(y_best, params)
    state = orbital.compute_orbital_state(r_post, v_post, gamma_post)
    insertion = orbital.evaluate_orbit_insertion(state, target_altitude)

    miss_km = (abs(state.apogee_altitude - target_altitude) / 1000.0
               if state.apogee_altitude is not None else None)

    if not insertion.meets_insertion_criterion:
        return DragConsistentCutoff(False, None, None, None, None, miss_km)

    return DragConsistentCutoff(
        found_passing_cutoff=True, cutoff_time=t_best, cutoff_state=y_best,
        post_atmosphere_state=state, insertion=insertion,
        closest_apogee_miss_km=miss_km,
    )
