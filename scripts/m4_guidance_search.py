"""M4 guidance search: retune the M3 gravity-turn law for the M4 orbit-capable study
vehicle at its reference payload (10,000 kg).

Structured coarse-to-fine search over (kick_start_time, kick_angle, kick_duration).
For each combination, the engine cutoff time is found via
``insertion_search.find_drag_consistent_cutoff``, which evaluates candidates using the
ACTUAL drag-included coast up through the atmosphere model's vacuum cutoff (100 km),
not an idealized instantaneous-cutoff estimate (DESIGN.md M4 S4 explains why this
distinction turned out to matter: an early idealized-only design point that looked like
it inserted into a 386 km orbit actually settled to 293 km once residual drag between
70-100 km altitude was correctly included).

Objective: minimize the idealized circularization delta-v (computed on the
drag-consistent post-atmosphere state) among (kick_start, kick_angle, kick_duration)
combinations that meet the M4 insertion criterion.

Usage:
    python scripts/m4_guidance_search.py
"""

import math

import numpy as np

from ascent import dynamics as dyn, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, M4_REFERENCE_VEHICLE, R_EARTH, TARGET_ALTITUDE
from ascent.insertion_search import find_drag_consistent_cutoff

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = M4_REFERENCE_VEHICLE


def build_params(vehicle=V):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT,
    )


def initial_state(vehicle=V):
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])


def evaluate(kick_start, kick_angle_deg, kick_duration, vehicle=V):
    params = build_params(vehicle)
    y0 = initial_state(vehicle)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    return find_drag_consistent_cutoff(vehicle, params, y0, kick_start, kick_angle_deg,
                                        kick_duration, t_burn, TARGET_ALTITUDE)


def coarse_search():
    kick_starts = [35.0, 45.0, 55.0]
    kick_angles = [15.0, 20.0, 25.0, 30.0]
    kick_durations = [15.0, 20.0, 25.0]

    results = []
    for ks in kick_starts:
        for ka in kick_angles:
            for kd in kick_durations:
                res = evaluate(ks, ka, kd)
                if res.found_passing_cutoff:
                    results.append((ks, ka, kd, res))
    return results


def fine_search(center_ks, center_ka, center_kd):
    kick_starts = np.arange(center_ks - 4.0, center_ks + 4.01, 2.0)
    kick_angles = np.arange(max(center_ka - 6.0, 1.0), center_ka + 6.01, 3.0)
    kick_durations = np.arange(max(center_kd - 8.0, 5.0), center_kd + 8.01, 4.0)

    results = []
    for ks in kick_starts:
        for ka in kick_angles:
            for kd in kick_durations:
                res = evaluate(ks, ka, kd)
                if res.found_passing_cutoff:
                    results.append((ks, ka, kd, res))
    return results


def main():
    print("=== M4 guidance search: coarse stage ===")
    coarse = coarse_search()
    print(f"{len(coarse)} coarse combinations meet the insertion criterion")
    for ks, ka, kd, res in coarse:
        ins = res.insertion
        print(f"  ks={ks} ka={ka} kd={kd} -> cutoff_t={res.cutoff_time:.2f} "
              f"circ_dv={ins.circularization_delta_v:.1f} "
              f"apogee_km={ins.apogee_altitude/1000:.2f}")

    if not coarse:
        raise RuntimeError("No coarse combination met the insertion criterion; widen "
                            "the search range before proceeding.")

    best_coarse = min(coarse, key=lambda r: r[3].insertion.circularization_delta_v)
    print(f"\nBest coarse combination: ks={best_coarse[0]}, ka={best_coarse[1]}, "
          f"kd={best_coarse[2]}, circ_dv={best_coarse[3].insertion.circularization_delta_v:.1f}")

    print("\n=== M4 guidance search: fine stage ===")
    fine = fine_search(best_coarse[0], best_coarse[1], best_coarse[2])
    print(f"{len(fine)} fine combinations meet the insertion criterion")
    candidates = fine if fine else [best_coarse]
    best = min(candidates, key=lambda r: r[3].insertion.circularization_delta_v)
    ks, ka, kd, res = best
    ins = res.insertion
    print(f"\nSelected M4 guidance: kick_start={ks} s, kick_angle={ka} deg, "
          f"kick_duration={kd} s, cutoff_time={res.cutoff_time:.2f} s")
    print(f"  apogee = {ins.apogee_altitude/1000:.2f} km")
    print(f"  circularization delta-v = {ins.circularization_delta_v:.1f} m/s")
    print(f"  meets_insertion_criterion = {ins.meets_insertion_criterion}")

    return dict(kick_start=ks, kick_angle_deg=ka, kick_duration=kd,
                cutoff_time=res.cutoff_time)


if __name__ == "__main__":
    main()
