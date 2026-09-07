"""M5: write the final authoritative inclination/payload table to CSV.

This script consolidates the results of the M5 guidance search and payload-boundary
solve into the committed authoritative CSV. It does NOT re-run the (expensive) search
itself for every inclination on every invocation -- see DESIGN.md M5 S6/S13 for why:
the near-polar case in particular required a targeted, hand-guided search (severe,
genuine numerical stiffness in the near-zero-rotational-boost regime, tracked down and
partially mitigated in dynamics.py/controls.py, but still far too slow for the
automatic coarse-to-fine search used at the other five inclinations -- a single
90 deg trajectory integration alone costs 15-30 s even after the V_FLOOR/deadlock
fixes, versus a fraction of a second at 28.5-70 deg).

Each row's (kick_start, kick_angle, kick_duration, cutoff_time, payload) was
independently found and verified (DESIGN.md M5 S6/S13): the five non-polar rows via
``scripts/m5_inclination_sweep.py``'s automatic coarse-to-fine search + drag-consistent
cutoff search (bracket+bisection at 0.1 kg tolerance), and the 90 deg row via a
targeted manual search over the SAME guidance-parameter space, verified with the SAME
drag-consistent (``insertion_search._coast_to_vacuum``) standard and the SAME M4
insertion criterion (``orbital.evaluate_orbit_insertion``, unmodified) -- confirmed
independently in this file's ``verify_row`` function, not merely asserted.

Usage:
    python scripts/m5_finalize_results.py
Writes: scripts/m5_inclination_sweep_results.csv
"""

import csv
import math
import os

import numpy as np

from ascent import controls, dynamics as dyn, inclination as inc, orbital, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle
from ascent.insertion_search import _coast_to_vacuum
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)

# (inclination_deg, payload_kg, kick_start_s, kick_angle_deg, kick_duration_s, cutoff_time_s)
CASES = [
    (28.5, 10334.0, 45.0, 20.0, 20.0, 251.2061383878732),
    (35.0, 10418.0, 45.0, 22.0, 20.0, 251.4820623959953),
    (45.0, 10327.0, 47.0, 30.0, 20.0, 251.78458053719987),
    (55.0, 10141.0, 46.0, 30.0, 18.0, 252.12050386916303),
    (70.0, 10258.0, 42.0, 30.0, 18.0, 252.7223685005077),
    (90.0, 10375.0, 39.7, 37.0, 20.0, 253.7002177608199),
]


def verify_row(i_deg, payload, ks, ka, kd, cutoff_time):
    """Independently re-integrate and re-verify one case against the M4 criterion."""
    az = inc.azimuth_from_inclination(math.radians(i_deg), LAT)
    vehicle = m4_vehicle(payload)
    params = dyn.AscentParams(isp=vehicle.isp, reference_area=vehicle.reference_area,
                               drag_coefficient=vehicle.drag_coefficient,
                               m_min=vehicle.m_min, latitude_rad=LAT, azimuth_rad=az)
    v0 = inc.useful_rotational_boost(az, LAT)
    y0 = np.array([R_EARTH, 0.0, v0, 0.0, vehicle.m0])
    t_burn_full = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    control = controls.gravity_turn_control(vehicle.thrust, ks, math.radians(ka), kd)

    t_eval = np.linspace(0.0, cutoff_time, 2000)
    result = run_ascent(y0, (0.0, cutoff_time), params, control, max_step=0.5,
                         rtol=1e-9, atol=1e-9, terminal_depletion=False,
                         include_apogee_event=False, t_eval=t_eval)

    y_cutoff = np.array([result.r[-1], result.theta[-1], result.v[-1],
                          result.gamma[-1], result.m[-1]])
    r_post, v_post, gamma_post = _coast_to_vacuum(y_cutoff, params)
    state = orbital.compute_orbital_state(r_post, v_post, gamma_post)
    insertion = orbital.evaluate_orbit_insertion(state, TARGET_ALTITUDE)
    assert insertion.meets_insertion_criterion, (
        f"Row for i={i_deg} deg, payload={payload} kg FAILED independent "
        f"re-verification against the M4 insertion criterion."
    )

    from ascent.atmosphere import density
    h = result.r - R_EARTH
    v_rel = np.array([
        dyn.relative_speed(v, g, r, LAT, azimuth_rad=az)
        for v, g, r in zip(result.v, result.gamma, result.r)
    ])
    q = 0.5 * density(h) * v_rel**2

    useful_boost = inc.useful_rotational_boost(az, LAT)
    cross_track = inc.cross_track_rotational_component(az, LAT)

    return dict(
        inclination_deg=i_deg, azimuth_deg=math.degrees(az),
        useful_rotational_boost_ms=useful_boost, cross_track_rotational_ms=cross_track,
        kick_start_s=ks, kick_angle_deg=ka, kick_duration_s=kd,
        max_payload_kg=payload, m0_kg=vehicle.m0, t_burn_full_s=t_burn_full,
        cutoff_time_s=result.t[-1], max_q_kPa=float(np.max(q)) / 1000.0,
        cutoff_altitude_km=h[-1] / 1000.0, cutoff_speed_ms=result.v[-1],
        cutoff_gamma_deg=math.degrees(result.gamma[-1]),
        apogee_km=state.apogee_altitude / 1000.0, perigee_km=state.perigee_altitude / 1000.0,
        circularization_dv_ms=insertion.circularization_delta_v,
    )


def main():
    rows = [verify_row(*case) for case in CASES]
    baseline = rows[0]["max_payload_kg"]
    for row in rows:
        row["payload_loss_vs_28p5_kg"] = baseline - row["max_payload_kg"]
        row["payload_loss_vs_28p5_pct"] = 100.0 * row["payload_loss_vs_28p5_kg"] / baseline
        print(row)

    out_path = os.path.join(os.path.dirname(__file__), "m5_inclination_sweep_results.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {len(rows)} verified rows to {out_path}")
    return rows


if __name__ == "__main__":
    main()
