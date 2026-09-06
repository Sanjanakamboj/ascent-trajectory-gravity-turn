"""M4 payload-to-orbit solve: maximum payload capability at 400 km / 28.5 deg for the
M4 orbit-capable study vehicle, using the guidance found by m4_guidance_search.py.

Fixed guidance profile (DESIGN.md M4 S4). Engine cutoff time is NOT fixed -- for each
payload, ``insertion_search.find_drag_consistent_cutoff`` finds the best cutoff using
the actual drag-included coast up through the atmosphere model's vacuum altitude
(100 km), which is the physically consistent standard used throughout M4 (see
DESIGN.md M4 S4 for why the naive idealized-instantaneous-cutoff estimate is NOT used).

Method (DESIGN.md M4 S5):
  1. Coarse sweep from the reference payload (10,000 kg, known PASS) upward in fixed
     steps until a FAIL is found, to establish a bracket.
  2. Bisection refinement of that bracket (on the PASS/FAIL boolean -- the underlying
     transition is a near-step function, see DESIGN.md M4 S4/S9, so bisecting the
     boolean directly is more robust here than forcing a smooth-function root-finder
     like Brent onto a discontinuity).
  3. Required case set (DESIGN.md M4 S7): zero payload, a clearly passing payload, a
     near-boundary payload, the refined maximum payload, a slightly heavier failing
     payload -- plus extra low-payload points that surface a genuine, investigated
     non-monotonicity (DESIGN.md M4 S8 check B).

Usage:
    python scripts/m4_payload_sweep.py
Writes: scripts/m4_payload_sweep_results.csv
"""

import csv
import math
import os

import numpy as np

from ascent import dynamics as dyn, propulsion as prop
from ascent.atmosphere import density
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle
from ascent.insertion_search import find_drag_consistent_cutoff
from ascent.simulation import run_ascent
from ascent import controls

LAT = math.radians(LAUNCH_LATITUDE_DEG)

# Selected guidance (DESIGN.md M4 S4, from scripts/m4_guidance_search.py).
KICK_START = 45.0
KICK_ANGLE_DEG = 20.0
KICK_DURATION = 20.0
REFERENCE_PAYLOAD = 10_000.0

BRACKET_STEP = 1000.0       # kg, coarse upward step from the reference payload
BRACKET_MAX_PAYLOAD = 100_000.0  # kg, hard stop -- fail loudly if no FAIL found by here
BISECTION_TOLERANCE = 1.0   # kg


def build_params(vehicle):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT,
    )


def initial_state(vehicle):
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])


def evaluate_payload(payload_mass):
    """Find the best drag-consistent insertion cutoff for this payload, plus diagnostics
    (max-Q up to cutoff, burnout state) evaluated on the SAME trajectory used to find it.
    """
    vehicle = m4_vehicle(payload_mass)
    params = build_params(vehicle)
    y0 = initial_state(vehicle)
    t_burn_full = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)

    cutoff = find_drag_consistent_cutoff(vehicle, params, y0, KICK_START, KICK_ANGLE_DEG,
                                          KICK_DURATION, t_burn_full, TARGET_ALTITUDE)

    # Re-integrate the powered phase up to the identified cutoff time (or full burn, if
    # none found) with dense sampling, purely for max-Q / burnout-state reporting.
    control = controls.gravity_turn_control(vehicle.thrust, KICK_START,
                                             math.radians(KICK_ANGLE_DEG), KICK_DURATION)
    t_end = cutoff.cutoff_time if cutoff.found_passing_cutoff else t_burn_full
    t_eval = np.linspace(0.0, t_end, 3000)
    result = run_ascent(y0, (0.0, t_end), params, control, max_step=0.3, rtol=1e-10,
                         atol=1e-10, terminal_depletion=False, include_apogee_event=False,
                         t_eval=t_eval)

    h = result.r - R_EARTH
    v_rel = np.array([
        dyn.relative_speed(v, g, r, LAT) for v, g, r in zip(result.v, result.gamma, result.r)
    ])
    q = 0.5 * density(h) * v_rel**2

    extra = dict(
        t_burn_full=t_burn_full, max_q_kPa=float(np.max(q)) / 1000.0,
        burnout_time_s=result.t[-1], burnout_altitude_km=h[-1] / 1000.0,
        burnout_speed_ms=result.v[-1], burnout_gamma_deg=math.degrees(result.gamma[-1]),
        burnout_mass_kg=result.m[-1],
    )
    return vehicle, cutoff, extra


def find_bracket():
    payload = REFERENCE_PAYLOAD
    cutoff = evaluate_payload(payload)[1]
    if not cutoff.found_passing_cutoff:
        raise RuntimeError(
            f"Reference payload {REFERENCE_PAYLOAD} kg does not pass the insertion "
            f"criterion with the selected guidance -- cannot establish a bracket. "
            f"closest_apogee_miss_km={cutoff.closest_apogee_miss_km}"
        )
    lo = payload
    hi = None
    while payload < BRACKET_MAX_PAYLOAD:
        payload += BRACKET_STEP
        cutoff = evaluate_payload(payload)[1]
        if not cutoff.found_passing_cutoff:
            hi = payload
            break
        lo = payload
    if hi is None:
        raise RuntimeError(
            f"No failing payload found up to {BRACKET_MAX_PAYLOAD} kg -- cannot "
            f"bracket the payload boundary. Refusing to guess; widen BRACKET_MAX_PAYLOAD "
            f"deliberately if this is expected."
        )
    return lo, hi


def bisect_boundary(lo, hi):
    assert evaluate_payload(lo)[1].found_passing_cutoff
    assert not evaluate_payload(hi)[1].found_passing_cutoff
    while hi - lo > BISECTION_TOLERANCE:
        mid = 0.5 * (lo + hi)
        if evaluate_payload(mid)[1].found_passing_cutoff:
            lo = mid
        else:
            hi = mid
    return lo, hi


def build_row(payload_mass):
    vehicle, cutoff, extra = evaluate_payload(payload_mass)
    row = dict(
        payload_kg=payload_mass, m0_kg=vehicle.m0, m_dry_kg=vehicle.m_dry,
        m_propellant_kg=vehicle.m_propellant, m_min_kg=vehicle.m_min,
        t_burn_full_s=extra["t_burn_full"], max_q_kPa=extra["max_q_kPa"],
        burnout_time_s=extra["burnout_time_s"],
        burnout_altitude_km=extra["burnout_altitude_km"],
        burnout_speed_ms=extra["burnout_speed_ms"],
        burnout_gamma_deg=extra["burnout_gamma_deg"],
        burnout_mass_kg=extra["burnout_mass_kg"],
        unburned_propellant_kg=extra["burnout_mass_kg"] - vehicle.m_min,
        orbit_achieved=cutoff.found_passing_cutoff,
    )
    if cutoff.found_passing_cutoff:
        ins = cutoff.insertion
        row.update(
            apogee_altitude_km=ins.apogee_altitude / 1000.0,
            circularization_delta_v_ms=ins.circularization_delta_v,
            closest_apogee_miss_km=0.0,
        )
    else:
        row.update(
            apogee_altitude_km=float("nan"),
            circularization_delta_v_ms=float("nan"),
            closest_apogee_miss_km=cutoff.closest_apogee_miss_km,
        )
    return row


def main():
    print("=== M4 payload bracket search ===")
    lo, hi = find_bracket()
    print(f"bracket: PASS at {lo} kg, FAIL at {hi} kg")

    print("=== M4 payload boundary bisection ===")
    lo, hi = bisect_boundary(lo, hi)
    max_payload = math.floor(lo)
    just_failing = math.ceil(hi)
    print(f"refined boundary: PASS at {lo:.3f} kg, FAIL at {hi:.3f} kg "
          f"(tolerance {BISECTION_TOLERANCE} kg)")
    print(f"reported maximum payload (floor of PASS boundary): {max_payload} kg")
    assert evaluate_payload(max_payload)[1].found_passing_cutoff
    assert not evaluate_payload(just_failing)[1].found_passing_cutoff

    required_cases = sorted(set([
        0.0, 5000.0, 7500.0,               # low-payload monotonicity investigation
        REFERENCE_PAYLOAD,                  # clearly passing
        max(max_payload - 1000.0, 0.0),     # near-boundary (passing)
        float(max_payload),                 # refined maximum (passing, by construction)
        float(just_failing),                # just past the boundary (failing, by construction)
        max_payload + 500.0,                # slightly heavier, clearly failing
    ]))

    rows = [build_row(p) for p in required_cases]

    out_path = os.path.join(os.path.dirname(__file__), "m4_payload_sweep_results.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {len(rows)} cases to {out_path}")
    for row in rows:
        print(row)

    return dict(max_payload_kg=max_payload, boundary_lo=lo, boundary_hi=hi, rows=rows)


if __name__ == "__main__":
    main()
