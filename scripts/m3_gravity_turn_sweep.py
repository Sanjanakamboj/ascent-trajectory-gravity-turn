"""M3 gravity-turn guidance parameter sweep (bounded, documented, small).

Sweeps two guidance parameters over a small, interpretable grid, using the UNCHANGED
M1/M2 physics-verification baseline vehicle:

  - pitch-kick start time:  {20, 30, 40, 50, 60} s        (5 values)
  - pitch-kick angle:       {5, 10, 15, 20, 25} deg        (5 values)
  - pitch-kick duration:    fixed at 10 s for every case (not swept; a documented,
                             single choice, per DESIGN.md M3 S3)

=> 25 cases total.

Objective metric (DESIGN.md M3 S3): maximize burnout specific orbital energy
(``epsilon = v^2/2 - mu/r`` at the propellant-depletion state). This is a transparent,
single, physically meaningful metric -- NOT "force orbital insertion."

This script does NOT change vehicle mass, thrust, Isp, staging, payload, or target
altitude. It only explores guidance-law timing/angle for the existing vehicle.

Each case is integrated from t=0 to t_burn (propellant-depletion event set terminal, so
a case that survives stops exactly at burnout; a case whose guidance drives it into the
ground first stops at ground impact instead -- both are reported honestly).

Usage:
    python scripts/m3_gravity_turn_sweep.py
Writes: scripts/m3_sweep_results.csv
"""

import csv
import math
import os

import numpy as np

from ascent import controls, dynamics as dyn, orbital, propulsion as prop
from ascent.atmosphere import density
from ascent.constants import BASELINE_VEHICLE, LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE

KICK_START_TIMES = [20.0, 30.0, 40.0, 50.0, 60.0]     # s
KICK_ANGLES_DEG = [5.0, 10.0, 15.0, 20.0, 25.0]        # deg
KICK_DURATION = 10.0                                    # s, fixed for every case


def build_params():
    return dyn.AscentParams(
        isp=V.isp, reference_area=V.reference_area, drag_coefficient=V.drag_coefficient,
        m_min=V.m_min, latitude_rad=LAT,
    )


def initial_state():
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, V.m0])


def run_case(kick_start_time, kick_angle_deg, kick_duration=KICK_DURATION):
    params = build_params()
    y0 = initial_state()
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    control = controls.gravity_turn_control(
        V.thrust, kick_start_time=kick_start_time,
        kick_angle_rad=math.radians(kick_angle_deg), kick_duration=kick_duration,
    )
    result = run_ascent(y0, (0.0, t_burn), params, control, max_step=0.5,
                         rtol=1e-9, atol=1e-9, terminal_depletion=True,
                         include_apogee_event=False)

    h = result.r - R_EARTH
    v_rel = np.array([
        dyn.relative_speed(v, g, r, LAT) for v, g, r in zip(result.v, result.gamma, result.r)
    ])
    q = 0.5 * density(h) * v_rel**2
    max_q = float(np.max(q))
    max_alt_reached = float(np.max(h))

    reached_burnout = result.t[-1] >= t_burn * 0.999 and len(result.t_events[1]) == 0

    row = dict(
        kick_start_time_s=kick_start_time,
        kick_angle_deg=kick_angle_deg,
        kick_duration_s=kick_duration,
        reached_burnout=reached_burnout,
        max_dynamic_pressure_kPa=max_q / 1000.0,
        max_altitude_reached_km=max_alt_reached / 1000.0,
    )

    if reached_burnout:
        state = orbital.compute_orbital_state(result.r[-1], result.v[-1], result.gamma[-1])
        orbit_achieved = orbital.is_circular_orbit_achieved(state, TARGET_ALTITUDE)
        row.update(
            burnout_time_s=result.t[-1],
            burnout_altitude_km=h[-1] / 1000.0,
            burnout_speed_ms=result.v[-1],
            burnout_gamma_deg=math.degrees(result.gamma[-1]),
            specific_energy_J_per_kg=state.specific_energy,
            eccentricity=state.eccentricity,
            perigee_altitude_km=(state.perigee_altitude / 1000.0
                                  if not state.intersects_earth else float("nan")),
            apogee_altitude_km=(state.apogee_altitude / 1000.0
                                 if state.apogee_altitude is not None else float("nan")),
            intersects_earth=state.intersects_earth,
            orbit_achieved=orbit_achieved,
        )
    else:
        row.update(
            burnout_time_s=float("nan"),
            burnout_altitude_km=float("nan"),
            burnout_speed_ms=float("nan"),
            burnout_gamma_deg=float("nan"),
            specific_energy_J_per_kg=float("nan"),
            eccentricity=float("nan"),
            perigee_altitude_km=float("nan"),
            apogee_altitude_km=max_alt_reached / 1000.0,  # still meaningful: how high it got
            intersects_earth=True,
            orbit_achieved=False,
        )
    return row


def main():
    rows = []
    for kick_start in KICK_START_TIMES:
        for kick_angle in KICK_ANGLES_DEG:
            rows.append(run_case(kick_start, kick_angle))

    out_path = os.path.join(os.path.dirname(__file__), "m3_sweep_results.csv")
    fieldnames = list(rows[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} cases to {out_path}")

    survivors = [r for r in rows if r["reached_burnout"]]
    n_orbit = sum(1 for r in rows if r["orbit_achieved"])
    print(f"cases reaching burnout: {len(survivors)}/{len(rows)}")
    print(f"cases achieving 400 km circular orbit: {n_orbit}/{len(rows)}")

    if survivors:
        best = max(survivors, key=lambda r: r["specific_energy_J_per_kg"])
        print("best case (max burnout specific orbital energy):")
        for k, v in best.items():
            print(f"  {k}: {v}")
        return best
    else:
        print("No case reached burnout without ground impact -- reporting sweep as-is.")
        return None


if __name__ == "__main__":
    main()
