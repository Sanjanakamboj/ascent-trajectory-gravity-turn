"""M2 diagnostic trajectory — prescribed control, NOT optimized for orbit.

Runs the unchanged M1 physics-verification baseline vehicle through a simple prescribed
vertical-rise + one-time pitch-kick control (controls.vertical_rise_then_pitch_kick_control),
using the M2 dynamics/atmosphere/propulsion modules, and reports the resulting altitude,
speed, flight-path angle, mass, and dynamic pressure history.

This is explicitly a physics-verification diagnostic, not a guided/optimized ascent:
no gravity-turn guidance law is applied (that begins in M3), and the pitch-kick angle
used here was not tuned to reach orbit or to maximize any performance metric.

Usage:
    python scripts/m2_diagnostic_trajectory.py
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ascent import atmosphere as atmo
from ascent import controls, dynamics as dyn, propulsion as prop
from ascent.constants import BASELINE_VEHICLE, LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
V = BASELINE_VEHICLE


def build_params():
    return dyn.AscentParams(
        isp=V.isp,
        reference_area=V.reference_area,
        drag_coefficient=V.drag_coefficient,
        m_min=V.m_min,
        latitude_rad=LAT,
    )


def main():
    params = build_params()

    # Diagnostic-only prescribed control (DESIGN.md M1 S1.7 concept; NOT a guidance
    # law): vertical for the first 10 s, then a fixed pitch-kick to 88 deg from
    # horizontal (a small, untuned 2 deg kick), held constant thereafter.
    control = controls.vertical_rise_then_pitch_kick_control(
        thrust=V.thrust,
        kick_start_time=10.0,
        kick_angle_rad=math.radians(2.0),
        chi_after_kick=math.radians(88.0),
    )

    # Initial condition per DESIGN.md S2: at liftoff the vehicle is already co-rotating
    # with the Earth, so its INERTIAL initial velocity is v_rot (purely horizontal,
    # gamma0 = 0), not zero. Using v0 = 0 here (as an earlier draft of this script did)
    # is inconsistent with M1 and was found, during M2 verification, to produce a
    # spurious "relative wind" of ~v_rot at t=0 (the co-rotating atmosphere then appears
    # to be moving relative to a vehicle that isn't actually co-rotating) -- an artifact,
    # not a real liftoff dynamic pressure spike. See DESIGN.md M2 section for the full
    # discussion, including the reason a purely radial ("vertical") thrust applied to
    # this nonzero-horizontal-velocity initial state curves gamma gradually rather than
    # producing an immediate steep climb (a real, correctly-modeled consequence of
    # describing liftoff in the inertial frame, analogous to a radial burn on a circular
    # orbit) -- not a bug, and not a case selected to make the trajectory look better.
    v_rot = dyn.atmosphere_corotation_speed(R_EARTH, LAT)
    y0 = np.array([R_EARTH, 0.0, v_rot, 0.0, V.m0])
    t_burn = prop.burn_duration(V.m_propellant, V.thrust, V.isp)

    # Integrate through burnout and well beyond, to see whether the vehicle is still
    # ascending, coasting, or has come back down within a generous diagnostic window.
    t_end = t_burn + 600.0
    t_eval = np.linspace(0.0, t_end, 4000)

    result = run_ascent(
        y0, (0.0, t_end), params, control,
        max_step=0.2, rtol=1e-9, atol=1e-9,
        target_altitude=TARGET_ALTITUDE, t_eval=t_eval,
    )

    h = result.r - R_EARTH
    q = 0.5 * atmo.density(h) * np.array([
        dyn.relative_speed(v, g, r, LAT)
        for v, g, r in zip(result.v, result.gamma, result.r)
    ]) ** 2

    idx_burnout = int(np.searchsorted(result.t, t_burn))
    idx_burnout = min(idx_burnout, len(result.t) - 1)

    max_alt_idx = int(np.argmax(h))
    max_speed_idx = int(np.argmax(result.v))
    max_q_idx = int(np.argmax(q))

    impact_times = result.t_events[1]
    depletion_times = result.t_events[0]

    v_circ = math.sqrt(3.986004418e14 / (R_EARTH + TARGET_ALTITUDE))
    reached_orbit = (
        abs(h[-1] - TARGET_ALTITUDE) < 1000.0
        and abs(result.v[-1] - v_circ) < 50.0
        and abs(result.gamma[-1]) < math.radians(1.0)
    )

    print("=== M2 diagnostic trajectory: prescribed control, NOT optimized for orbit ===")
    print(f"integration success: {result.success} ({result.message})")
    print(f"burn duration (M1 hand calc): {t_burn:.2f} s")
    print(f"burnout time index t = {result.t[idx_burnout]:.2f} s")
    print(f"burnout mass = {result.m[idx_burnout]:.1f} kg (m_min = {V.m_min:.1f} kg)")
    print(f"max altitude = {h[max_alt_idx]/1000:.2f} km at t = {result.t[max_alt_idx]:.1f} s")
    print(f"max speed (inertial) = {result.v[max_speed_idx]:.1f} m/s at t = {result.t[max_speed_idx]:.1f} s")
    print(f"max dynamic pressure = {q[max_q_idx]/1000:.2f} kPa at t = {result.t[max_q_idx]:.1f} s")
    print(f"propellant depletion event fired at t = {depletion_times}")
    print(f"ground impact event fired at t = {impact_times}")
    print(f"final altitude = {h[-1]/1000:.2f} km, final speed = {result.v[-1]:.1f} m/s, "
          f"final gamma = {math.degrees(result.gamma[-1]):.2f} deg")
    print(f"target circular speed at {TARGET_ALTITUDE/1000:.0f} km = {v_circ:.1f} m/s")
    print(f"400 km circular-orbit conditions achieved: {reached_orbit}")

    make_figure(result.t, h, result.v, result.gamma, result.m, q)

    return dict(
        t=result.t, h=h, v=result.v, gamma=result.gamma, m=result.m, q=q,
        t_burn=t_burn, max_alt=h[max_alt_idx], max_speed=result.v[max_speed_idx],
        max_q=q[max_q_idx], burnout_mass=result.m[idx_burnout],
        impact_times=impact_times, reached_orbit=reached_orbit,
    )


def make_figure(t, h, v, gamma, m, q):
    fig, axes = plt.subplots(3, 2, figsize=(11, 9))
    fig.suptitle(
        "M2 physics-verification trajectory — prescribed control, NOT optimized for orbit\n"
        "(diagnostic/supporting only; see DESIGN.md M2 section)",
        fontsize=10,
    )

    axes[0, 0].plot(t, h / 1000.0)
    axes[0, 0].set_ylabel("altitude [km]")
    axes[0, 0].set_xlabel("t [s]")

    axes[0, 1].plot(t, v)
    axes[0, 1].set_ylabel("inertial speed [m/s]")
    axes[0, 1].set_xlabel("t [s]")

    axes[1, 0].plot(t, np.degrees(gamma))
    axes[1, 0].set_ylabel("flight-path angle [deg]")
    axes[1, 0].set_xlabel("t [s]")

    axes[1, 1].plot(t, m / 1000.0)
    axes[1, 1].set_ylabel("mass [t]")
    axes[1, 1].set_xlabel("t [s]")

    axes[2, 0].plot(t, q / 1000.0)
    axes[2, 0].set_ylabel("dynamic pressure [kPa]")
    axes[2, 0].set_xlabel("t [s]")

    axes[2, 1].axis("off")
    axes[2, 1].text(
        0.0, 0.5,
        "Baseline vehicle does NOT reach\n400 km circular LEO (see DESIGN.md S7.4:\n"
        "M1 ideal delta-v ~5.05 km/s vs a rough\n~8.9-9.5 km/s LEO budget).\n\n"
        "No gravity-turn guidance law is applied\n(that begins in M3); pitch profile here\n"
        "is a fixed, untuned diagnostic kick.",
        fontsize=9, va="center",
    )

    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
    out_path = os.path.join(out_dir, "m2_diagnostic_trajectory.png")
    fig.savefig(out_path, dpi=150)
    print(f"figure saved to {os.path.abspath(out_path)}")


if __name__ == "__main__":
    main()
