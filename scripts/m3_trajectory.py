"""M3 selected gravity-turn trajectory — full profile, figures, and convergence study.

Uses the guidance parameters selected by scripts/m3_gravity_turn_sweep.py (maximum
burnout specific orbital energy among cases that survive to burnout):

    kick_start_time = 40 s, kick_angle = 25 deg, kick_duration = 10 s

on the UNCHANGED M1/M2 physics-verification baseline vehicle. This script does not
change the vehicle in any way; it only propagates further (through coast) than the
sweep did, to show apogee/impact behavior and produce the M3 trajectory figures.

Usage:
    python scripts/m3_trajectory.py
Writes: figures/m3_gravity_turn_trajectory.png, figures/m3_velocity_components.png
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ascent import atmosphere as atmo
from ascent import controls, dynamics as dyn, orbital, propulsion as prop
from ascent.constants import BASELINE_VEHICLE, LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE

KICK_START = 40.0
KICK_ANGLE_DEG = 25.0
KICK_DURATION = 10.0


def build_params(cd=None):
    return dyn.AscentParams(
        isp=V.isp, reference_area=V.reference_area,
        drag_coefficient=V.drag_coefficient if cd is None else cd,
        m_min=V.m_min, latitude_rad=LAT,
    )


def initial_state():
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    return np.array([R_EARTH, 0.0, v_rot, 0.0, V.m0])


def build_control():
    return controls.gravity_turn_control(V.thrust, KICK_START, math.radians(KICK_ANGLE_DEG),
                                          KICK_DURATION)


def run_full_trajectory(rtol=1e-9, atol=1e-9, max_step=0.5):
    params = build_params()
    y0 = initial_state()
    control = build_control()
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)
    t_end = t_burn + 1500.0
    t_eval = np.linspace(0.0, t_end, 6000)
    return run_ascent(y0, (0.0, t_end), params, control, max_step=max_step,
                       rtol=rtol, atol=atol, target_altitude=TARGET_ALTITUDE,
                       t_eval=t_eval), t_burn


def summarize(result, t_burn):
    h = result.r - R_EARTH
    v_rel = np.array([
        dyn.relative_speed(v, g, r, LAT) for v, g, r in zip(result.v, result.gamma, result.r)
    ])
    q = 0.5 * atmo.density(h) * v_rel**2

    idx_burnout = int(np.searchsorted(result.t, t_burn))
    idx_burnout = min(idx_burnout, len(result.t) - 1)
    # Ascent max-Q (aerodynamic-loads-relevant) vs. the global max, which -- for a
    # suborbital trajectory whose perigee is deep inside the Earth -- occurs instead
    # during the high-speed terminal re-entry/impact dive, a physically real but
    # distinct phenomenon from ascent max-Q. Both are reported; only the ascent value
    # is marked as "max-Q" on the figure to avoid a misleading label.
    idx_max_q_ascent = int(np.argmax(q[:idx_burnout + 1]))
    idx_max_q = int(np.argmax(q))
    idx_max_alt = int(np.argmax(h))

    burnout_state = orbital.compute_orbital_state(result.r[idx_burnout], result.v[idx_burnout],
                                                   result.gamma[idx_burnout])
    orbit_achieved = orbital.is_circular_orbit_achieved(burnout_state, TARGET_ALTITUDE)

    apogee_times = result.t_events[2] if len(result.t_events) > 2 else np.array([])
    impact_times = result.t_events[1]

    return dict(
        h=h, q=q, idx_burnout=idx_burnout, idx_max_q=idx_max_q,
        idx_max_q_ascent=idx_max_q_ascent, idx_max_alt=idx_max_alt,
        burnout_state=burnout_state, orbit_achieved=orbit_achieved,
        apogee_times=apogee_times, impact_times=impact_times,
    )


def make_main_figure(result, summary, path):
    t = result.t
    fig, axes = plt.subplots(3, 2, figsize=(11, 9))
    fig.suptitle(
        "M3 gravity-turn trajectory — unchanged M1/M2 verification vehicle\n"
        f"kick_start={KICK_START:.0f}s, kick_angle={KICK_ANGLE_DEG:.0f}deg, "
        f"kick_duration={KICK_DURATION:.0f}s — orbit achieved: {summary['orbit_achieved']}",
        fontsize=10,
    )

    def mark(ax, idx, label, color):
        ax.axvline(t[idx], color=color, linestyle="--", linewidth=1, alpha=0.7)

    idx_b = summary["idx_burnout"]
    idx_q = summary["idx_max_q_ascent"]
    idx_kick_start = int(np.searchsorted(t, KICK_START))
    idx_kick_end = int(np.searchsorted(t, KICK_START + KICK_DURATION))

    panels = [
        (axes[0, 0], summary["h"] / 1000.0, "altitude [km]"),
        (axes[0, 1], result.v, "inertial speed [m/s]"),
        (axes[1, 0], np.degrees(result.gamma), "flight-path angle [deg]"),
        (axes[1, 1], summary["q"] / 1000.0, "dynamic pressure [kPa]"),
        (axes[2, 0], result.m / 1000.0, "mass [t]"),
    ]
    for ax, y, ylabel in panels:
        ax.plot(t, y)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("t [s]")
        for idx, color in [(idx_kick_start, "tab:orange"), (idx_kick_end, "tab:orange"),
                           (idx_b, "tab:red"), (idx_q, "tab:green")]:
            mark(ax, idx, "", color)

    apogee_t = summary["apogee_times"]
    if len(apogee_t) > 0:
        for ax, _, _ in panels:
            ax.axvline(apogee_t[0], color="tab:purple", linestyle=":", linewidth=1, alpha=0.7)

    axes[2, 1].axis("off")
    q_ascent_kPa = summary["q"][summary["idx_max_q_ascent"]] / 1000.0
    q_global_kPa = summary["q"][summary["idx_max_q"]] / 1000.0
    legend_text = (
        "orange dashed: pitch-kick start/end\n"
        "red dashed: burnout (propellant depletion)\n"
        "green dashed: max-Q DURING ASCENT "
        f"({q_ascent_kPa:.1f} kPa)\n"
        "purple dotted: apogee (if reached in window)\n\n"
        f"NOTE: this suborbital trajectory's perigee is\n"
        f"deep inside the Earth, so it re-enters and hits\n"
        f"the ground; the terminal re-entry dive reaches a\n"
        f"much higher dynamic pressure ({q_global_kPa:.0f} kPa, not\n"
        f"marked above) than ascent max-Q -- a real but\n"
        f"separate phenomenon from ascent aerodynamic loads.\n\n"
        f"Baseline vehicle does NOT reach 400 km circular\n"
        f"LEO with this (or any swept) gravity-turn\n"
        f"parameter set -- see DESIGN.md M3 S9. Diagnostic/\n"
        f"supporting only, not a validated flight result."
    )
    axes[2, 1].text(0.0, 0.5, legend_text, fontsize=8.5, va="center")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"figure saved to {os.path.abspath(path)}")


def make_velocity_figure(result, summary, path):
    v_radial = result.v * np.sin(result.gamma)
    v_tangential = result.v * np.cos(result.gamma)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    fig.suptitle(
        "M3 gravity-turn trajectory — velocity components and ground track (diagnostic)",
        fontsize=10,
    )

    axes[0].plot(result.t, v_radial, label="radial")
    axes[0].plot(result.t, v_tangential, label="tangential")
    axes[0].set_xlabel("t [s]")
    axes[0].set_ylabel("velocity component [m/s]")
    axes[0].legend()

    axes[1].plot(np.degrees(result.theta), summary["h"] / 1000.0)
    axes[1].set_xlabel("downrange angle theta [deg]")
    axes[1].set_ylabel("altitude [km]")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"figure saved to {os.path.abspath(path)}")


def convergence_study():
    settings = [
        dict(max_step=2.0, rtol=1e-6, atol=1e-6),
        dict(max_step=0.5, rtol=1e-9, atol=1e-9),
        dict(max_step=0.1, rtol=1e-12, atol=1e-12),
    ]
    params = build_params()
    y0 = initial_state()
    control = build_control()
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)

    rows = []
    for s in settings:
        result = run_ascent(y0, (0.0, t_burn), params, control, terminal_depletion=True,
                             include_apogee_event=False, **s)
        h = result.r - R_EARTH
        v_rel = np.array([
            dyn.relative_speed(v, g, r, LAT)
            for v, g, r in zip(result.v, result.gamma, result.r)
        ])
        q = 0.5 * atmo.density(h) * v_rel**2
        state = orbital.compute_orbital_state(result.r[-1], result.v[-1], result.gamma[-1])
        rows.append(dict(
            settings=s, burnout_time=result.t[-1], burnout_altitude_km=h[-1] / 1000.0,
            burnout_speed=result.v[-1], burnout_gamma_deg=math.degrees(result.gamma[-1]),
            max_q_kPa=float(np.max(q)) / 1000.0, apogee_altitude_km=(
                state.apogee_altitude / 1000.0 if state.apogee_altitude is not None
                else float("nan")),
        ))
    return rows


def main():
    result, t_burn = run_full_trajectory()
    summary = summarize(result, t_burn)

    fig_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
    make_main_figure(result, summary, os.path.join(fig_dir, "m3_gravity_turn_trajectory.png"))
    make_velocity_figure(result, summary, os.path.join(fig_dir, "m3_velocity_components.png"))

    bs = summary["burnout_state"]
    print("=== M3 selected gravity-turn trajectory ===")
    print(f"kick_start={KICK_START}, kick_angle={KICK_ANGLE_DEG} deg, kick_duration={KICK_DURATION} s")
    print(f"burnout time = {t_burn:.2f} s, altitude = {bs.altitude/1000:.2f} km, "
          f"speed = {bs.v:.1f} m/s, gamma = {math.degrees(bs.gamma):.2f} deg")
    print(f"specific orbital energy = {bs.specific_energy:.3e} J/kg, eccentricity = {bs.eccentricity:.4f}")
    print(f"perigee altitude = {bs.perigee_altitude/1000:.2f} km "
          f"(intersects Earth: {bs.intersects_earth})")
    print(f"apogee altitude = {bs.apogee_altitude/1000 if bs.apogee_altitude else float('nan'):.2f} km")
    print(f"orbit achieved (400 km circular criterion): {summary['orbit_achieved']}")
    print(f"apogee event times: {summary['apogee_times']}")
    print(f"ground impact event times: {summary['impact_times']}")
    print(f"max dynamic pressure DURING ASCENT (t<=burnout) = "
          f"{summary['q'][summary['idx_max_q_ascent']]/1000:.2f} kPa")
    print(f"max dynamic pressure OVERALL (terminal re-entry dive) = "
          f"{summary['q'][summary['idx_max_q']]/1000:.2f} kPa at "
          f"t={result.t[summary['idx_max_q']]:.1f} s")
    print(f"max altitude reached in window = {np.max(summary['h'])/1000:.2f} km")

    print("\n=== Integrator convergence study ===")
    for row in convergence_study():
        print(row)


if __name__ == "__main__":
    main()
