"""M4 headline figure: payload mass vs. insertion capability.

Runs its own finer payload sweep (denser than the required-case CSV from
scripts/m4_payload_sweep.py) purely for plotting, using the same selected guidance and
the same drag-consistent insertion search as the rest of M4.

Usage:
    python scripts/m4_payload_figure.py
Writes: figures/m4_payload_capability.png
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ascent import dynamics as dyn, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle
from ascent.insertion_search import find_drag_consistent_cutoff

LAT = math.radians(LAUNCH_LATITUDE_DEG)
KICK_START = 45.0
KICK_ANGLE_DEG = 20.0
KICK_DURATION = 20.0
MAX_PAYLOAD_KG = 10_333.0  # from scripts/m4_payload_sweep.py


def build_params(vehicle):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT,
    )


def initial_state(vehicle):
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])


def evaluate(payload):
    vehicle = m4_vehicle(payload)
    params = build_params(vehicle)
    y0 = initial_state(vehicle)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    return find_drag_consistent_cutoff(vehicle, params, y0, KICK_START, KICK_ANGLE_DEG,
                                        KICK_DURATION, t_burn, TARGET_ALTITUDE)


def main():
    low_payloads = [0.0, 2500.0, 5000.0, 6000.0, 6500.0]
    fine_payloads = list(np.arange(7000.0, 11700.0, 250.0))
    payloads = sorted(set(low_payloads + fine_payloads))

    rows = []
    for p in payloads:
        res = evaluate(p)
        if res.found_passing_cutoff:
            rows.append(dict(payload=p, passed=True,
                              apogee_km=res.post_atmosphere_state.apogee_altitude / 1000.0,
                              circ_dv=res.insertion.circularization_delta_v))
        else:
            rows.append(dict(payload=p, passed=False,
                              apogee_km=(res.closest_apogee_miss_km + TARGET_ALTITUDE / 1000.0
                                         if res.closest_apogee_miss_km is not None else np.nan),
                              circ_dv=np.nan))
        print(rows[-1])

    payload_arr = np.array([r["payload"] for r in rows])
    passed_arr = np.array([r["passed"] for r in rows])
    apogee_arr = np.array([r["apogee_km"] for r in rows])
    circ_dv_arr = np.array([r["circ_dv"] for r in rows])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("M4 payload capability -- M4 orbit-capable study vehicle, 400 km / 28.5 deg",
                 fontsize=11)

    ax = axes[0]
    pass_mask = passed_arr
    ax.scatter(payload_arr[pass_mask] / 1000.0, apogee_arr[pass_mask], color="tab:green",
               label="PASS (meets insertion criterion)", zorder=3)
    ax.scatter(payload_arr[~pass_mask] / 1000.0, apogee_arr[~pass_mask], color="tab:red",
               marker="x", label="FAIL", zorder=3)
    ax.axhline(TARGET_ALTITUDE / 1000.0, color="black", linestyle="--", linewidth=1,
               label="target altitude (400 km)")
    ax.axhspan((TARGET_ALTITUDE - 15_000) / 1000.0, (TARGET_ALTITUDE + 15_000) / 1000.0,
               color="black", alpha=0.06, label="+/- 15 km tolerance")
    ax.axvline(MAX_PAYLOAD_KG / 1000.0, color="tab:blue", linestyle=":", linewidth=1.5,
               label=f"max payload ({MAX_PAYLOAD_KG/1000:.2f} t)")
    ax.set_xlabel("payload [t]")
    ax.set_ylabel("apogee altitude [km]\n(FAIL points: closest achievable, for reference)")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.scatter(payload_arr[pass_mask] / 1000.0, circ_dv_arr[pass_mask], color="tab:green",
               zorder=3)
    ax.axvline(MAX_PAYLOAD_KG / 1000.0, color="tab:blue", linestyle=":", linewidth=1.5,
               label=f"max payload ({MAX_PAYLOAD_KG/1000:.2f} t)")
    ax.set_xlabel("payload [t]")
    ax.set_ylabel("idealized circularization delta-v [m/s]\n(PASSing cases only)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    out_path = os.path.join(os.path.dirname(__file__), "..", "figures",
                             "m4_payload_capability.png")
    fig.savefig(out_path, dpi=150)
    print(f"figure saved to {os.path.abspath(out_path)}")


if __name__ == "__main__":
    main()
