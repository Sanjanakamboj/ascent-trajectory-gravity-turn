"""M4 representative orbit-capable ascent trajectory (reference payload = 10,000 kg).

Produces the M4 headline ascent-trajectory figure and the delta-v/loss budget for the
selected guidance (kick_start=45 deg... see below) on the M4 orbit-capable study
vehicle -- explicitly NOT the M1-M3 verification vehicle. Uses the drag-consistent
engine-cutoff search (``insertion_search.find_drag_consistent_cutoff``), the same
standard used throughout M4 (DESIGN.md M4 S4).

Usage:
    python scripts/m4_trajectory.py
Writes: figures/m4_orbit_capable_trajectory.png, figures/m4_delta_v_budget.png
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ascent import atmosphere as atmo
from ascent import controls, dynamics as dyn, propulsion as prop
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle
from ascent.insertion_search import find_drag_consistent_cutoff
from ascent.loss_budget import compute_delta_v_budget
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
REFERENCE_PAYLOAD = 10_000.0
KICK_START = 45.0
KICK_ANGLE_DEG = 20.0
KICK_DURATION = 20.0


def build_params(vehicle):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT,
    )


def initial_state(vehicle):
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, vehicle.m0])


def run_representative(payload=REFERENCE_PAYLOAD):
    vehicle = m4_vehicle(payload)
    params = build_params(vehicle)
    y0 = initial_state(vehicle)
    t_burn_full = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)

    cutoff = find_drag_consistent_cutoff(vehicle, params, y0, KICK_START, KICK_ANGLE_DEG,
                                          KICK_DURATION, t_burn_full, TARGET_ALTITUDE)
    assert cutoff.found_passing_cutoff, "reference payload must pass -- check guidance/vehicle"

    control = controls.gravity_turn_control(vehicle.thrust, KICK_START,
                                             math.radians(KICK_ANGLE_DEG), KICK_DURATION)
    t_eval = np.linspace(0.0, cutoff.cutoff_time, 4000)
    powered_result = run_ascent(y0, (0.0, cutoff.cutoff_time), params, control, max_step=0.2,
                                 rtol=1e-10, atol=1e-10, terminal_depletion=False,
                                 include_apogee_event=False, t_eval=t_eval)

    # Continue past cutoff as an unpowered (but drag-included) coast, for the figure --
    # this is the SAME physical coast used to find the cutoff, just continued further
    # (through vacuum and on to apogee) for plotting.
    y_cutoff = cutoff.cutoff_state
    coast = controls.coast_control()
    coast_result = run_ascent(y_cutoff, (0.0, 3000.0), params, coast, max_step=1.0,
                               rtol=1e-10, atol=1e-10, include_apogee_event=True)

    return vehicle, powered_result, cutoff, coast_result


def combine_for_plot(powered_result, coast_result):
    t = np.concatenate([powered_result.t, powered_result.t[-1] + coast_result.t[1:]])
    r = np.concatenate([powered_result.r, coast_result.r[1:]])
    v = np.concatenate([powered_result.v, coast_result.v[1:]])
    gamma = np.concatenate([powered_result.gamma, coast_result.gamma[1:]])
    m = np.concatenate([powered_result.m, coast_result.m[1:]])
    return t, r, v, gamma, m


def make_trajectory_figure(t, r, v, gamma, m, cutoff_t, apogee_t, path, vehicle, cutoff):
    h = r - R_EARTH
    v_rel = np.array([dyn.relative_speed(vi, gi, ri, LAT) for vi, gi, ri in zip(v, gamma, r)])
    q = 0.5 * atmo.density(h) * v_rel**2

    fig, axes = plt.subplots(3, 2, figsize=(11, 9))
    fig.suptitle(
        "M4 orbit-capable ascent trajectory -- M4 STUDY VEHICLE (NOT the M1-M3 "
        "verification vehicle)\n"
        f"payload={vehicle.m_payload/1000:.1f} t, kick_start={KICK_START:.0f}s, "
        f"kick_angle={KICK_ANGLE_DEG:.0f}deg, kick_duration={KICK_DURATION:.0f}s, "
        f"cutoff at {cutoff_t:.1f}s (before full depletion)",
        fontsize=9.5,
    )

    idx_kick_start = int(np.searchsorted(t, KICK_START))
    idx_kick_end = int(np.searchsorted(t, KICK_START + KICK_DURATION))
    idx_cutoff = int(np.searchsorted(t, cutoff_t))
    idx_max_q = int(np.argmax(q[:idx_cutoff + 1]))
    idx_apogee = int(np.searchsorted(t, apogee_t)) if apogee_t is not None else None

    panels = [
        (axes[0, 0], h / 1000.0, "altitude [km]"),
        (axes[0, 1], v, "inertial speed [m/s]"),
        (axes[1, 0], np.degrees(gamma), "flight-path angle [deg]"),
        (axes[1, 1], q / 1000.0, "dynamic pressure [kPa]"),
        (axes[2, 0], m / 1000.0, "mass [t]"),
    ]
    for ax, y, ylabel in panels:
        ax.plot(t, y)
        ax.set_xlabel("t [s]")
        ax.set_ylabel(ylabel)
        for idx, color in [(idx_kick_start, "tab:orange"), (idx_kick_end, "tab:orange"),
                           (idx_cutoff, "tab:red"), (idx_max_q, "tab:green")]:
            ax.axvline(t[idx], color=color, linestyle="--", linewidth=1, alpha=0.7)
        if idx_apogee is not None:
            ax.axvline(t[idx_apogee], color="tab:purple", linestyle=":", linewidth=1, alpha=0.7)

    axes[2, 1].axis("off")
    ins = cutoff.insertion
    text = (
        "orange dashed: pitch-kick start/end\n"
        "red dashed: commanded engine cutoff\n"
        "green dashed: max-Q (up to cutoff)\n"
        "purple dotted: apogee (drag-included coast)\n\n"
        f"apogee (post-atmosphere, drag-consistent) =\n"
        f"  {ins.apogee_altitude/1000:.1f} km (target 400 km +/- 15 km)\n"
        f"idealized circularization dv = {ins.circularization_delta_v:.1f} m/s\n"
        f"(a separate, unconstrained insertion-stage\n"
        f"impulse -- NOT drawn from this vehicle's own\n"
        f"propellant budget; see DESIGN.md M4 S2)\n\n"
        f"M4 insertion criterion met: {ins.meets_insertion_criterion}"
    )
    axes[2, 1].text(0.0, 0.5, text, fontsize=8.5, va="center")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"figure saved to {os.path.abspath(path)}")


def make_loss_budget_figure(budget, path):
    fig, ax = plt.subplots(figsize=(7, 5))
    labels = ["achieved\ndv", "gravity\nloss", "drag\nloss", "steering\nloss"]
    values = [budget.achieved_dv, budget.gravity_loss, budget.drag_loss, budget.steering_loss]
    colors = ["tab:blue", "tab:red", "tab:orange", "tab:purple"]
    bottoms = np.cumsum([0] + values[:-1])
    for label, value, bottom, color in zip(labels, values, bottoms, colors):
        ax.bar(0, value, bottom=bottom, color=color, label=label, width=0.5)
    ax.axhline(budget.ideal_dv_to_cutoff, color="black", linestyle="--", linewidth=1,
               label=f"ideal dv to cutoff ({budget.ideal_dv_to_cutoff:.0f} m/s)")
    ax.set_xticks([])
    ax.set_ylabel("delta-v [m/s]")
    ax.set_title("M4 delta-v / loss budget (reference payload, to engine cutoff)", pad=14)
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(path, dpi=150)
    print(f"figure saved to {os.path.abspath(path)}")


def main():
    vehicle, powered_result, cutoff, coast_result = run_representative()
    cutoff_t = powered_result.t[-1]

    t, r, v, gamma, m = combine_for_plot(powered_result, coast_result)
    apogee_events = coast_result.t_events[2] if len(coast_result.t_events) > 2 else np.array([])
    apogee_t = cutoff_t + apogee_events[0] if len(apogee_events) > 0 else None

    fig_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
    make_trajectory_figure(t, r, v, gamma, m, cutoff_t, apogee_t,
                            os.path.join(fig_dir, "m4_orbit_capable_trajectory.png"),
                            vehicle, cutoff)

    budget = compute_delta_v_budget(
        powered_result.t, powered_result.r, powered_result.v, powered_result.gamma,
        powered_result.m, vehicle.isp, vehicle.drag_coefficient, vehicle.reference_area,
        LAT, vehicle.m0,
    )
    make_loss_budget_figure(budget, os.path.join(fig_dir, "m4_delta_v_budget.png"))

    print("=== M4 representative trajectory (reference payload) ===")
    print(f"vehicle m0={vehicle.m0} kg, payload={vehicle.m_payload} kg")
    print(f"cutoff time = {cutoff_t:.2f} s, altitude = {(powered_result.r[-1]-R_EARTH)/1000:.2f} km")
    print(f"burnout speed = {powered_result.v[-1]:.1f} m/s, gamma = "
          f"{math.degrees(powered_result.gamma[-1]):.3f} deg")
    print(f"burnout mass = {powered_result.m[-1]:.1f} kg (m_min={vehicle.m_min} kg, "
          f"unburned propellant = {powered_result.m[-1]-vehicle.m_min:.1f} kg)")
    ins = cutoff.insertion
    print(f"apogee (drag-consistent, post-atmosphere) = {ins.apogee_altitude/1000:.2f} km, "
          f"circularization dv = {ins.circularization_delta_v:.1f} m/s")
    print(f"insertion criterion met: {ins.meets_insertion_criterion}")
    print("\n--- delta-v / loss budget (to cutoff) ---")
    print(f"ideal dv to cutoff (Tsiolkovsky, mass consumed by cutoff) = "
          f"{budget.ideal_dv_to_cutoff:.1f} m/s")
    print(f"achieved inertial dv                                     = "
          f"{budget.achieved_dv:.1f} m/s")
    print(f"gravity loss                                             = "
          f"{budget.gravity_loss:.1f} m/s")
    print(f"drag loss                                                = "
          f"{budget.drag_loss:.1f} m/s")
    print(f"steering loss (exact remainder, see loss_budget.py)      = "
          f"{budget.steering_loss:.1f} m/s")
    print(f"identity residual (should be ~0)                         = "
          f"{budget.residual:.6e} m/s")


if __name__ == "__main__":
    main()
