"""M5: direct-ascent payload-to-orbit vs. target inclination for the M4 orbit-capable
study vehicle -- the AUTOMATIC coarse-to-fine guidance search used to find 5 of the 6
authoritative points (28.5-70 deg).

NOTE (DESIGN.md M5 S6/S13): the 90 deg (polar) point is NOT reliably solvable by this
automatic search in practical time -- the near-zero rotational boost there triggers
genuine, severe numerical stiffness (documented and partially mitigated in
dynamics.py/controls.py's V_FLOOR handling), making each trial trajectory 15-30x more
expensive than at other inclinations. 90 deg was instead solved with a targeted manual
search (see DESIGN.md M5 S6). The single, independently-reverified source of truth for
ALL SIX points is ``scripts/m5_finalize_results.py`` and its output CSV -- run THAT
script to (re)produce the committed ``m5_inclination_sweep_results.csv``. This script
remains as a runnable record of the automatic-search METHOD (and reproduces the same
28.5-70 deg numbers), not as the final data pipeline.

For each target inclination:
  1. compute the direct-ascent launch azimuth (inclination.azimuth_from_inclination),
  2. re-optimize guidance (kick_start, kick_angle, kick_duration) in a small local
     search seeded from the nearest already-solved inclination (coarse-to-fine per
     inclination, per DESIGN.md M5 S6 -- NOT a frozen 28.5-deg guidance reused
     everywhere),
  3. solve the maximum payload via the SAME bracket+bisection method and the SAME M4
     insertion criterion as scripts/m4_payload_sweep.py (unchanged: bound,
     non-Earth-intersecting, apogee within 15 km of 400 km, circularization dv <=
     1500 m/s).

Usage:
    python scripts/m5_inclination_sweep.py
Writes: scripts/m5_inclination_sweep_results.csv (superseded by
        scripts/m5_finalize_results.py's output of the same filename -- re-run THAT
        script, not this one, to restore the committed authoritative CSV)
"""

import csv
import math
import os

import numpy as np

from ascent import dynamics as dyn, inclination as inc, propulsion as prop
from ascent.atmosphere import density
from ascent.constants import LAUNCH_LATITUDE_DEG, R_EARTH, TARGET_ALTITUDE, m4_vehicle
from ascent.insertion_search import find_drag_consistent_cutoff
from ascent import controls
from ascent.simulation import run_ascent

LAT = math.radians(LAUNCH_LATITUDE_DEG)
REFERENCE_PAYLOAD = 10_000.0
BISECTION_TOLERANCE = 0.1
BRACKET_STEP = 1000.0
BRACKET_MAX_PAYLOAD = 100_000.0

AUTHORITATIVE_INCLINATIONS_DEG = [28.5, 35.0, 45.0, 55.0, 70.0, 90.0]
FINE_INCLINATIONS_DEG = list(np.arange(28.5, 90.001, 5.0))
if 90.0 not in FINE_INCLINATIONS_DEG:
    FINE_INCLINATIONS_DEG.append(90.0)

# Known-good starting guidance for 28.5 deg (DESIGN.md M4 S4).
SEED_GUIDANCE = {28.5: (45.0, 20.0, 20.0)}


def build_params(vehicle, azimuth_rad):
    return dyn.AscentParams(
        isp=vehicle.isp, reference_area=vehicle.reference_area,
        drag_coefficient=vehicle.drag_coefficient, m_min=vehicle.m_min,
        latitude_rad=LAT, azimuth_rad=azimuth_rad,
    )


def initial_state(vehicle, azimuth_rad):
    v0 = inc.useful_rotational_boost(azimuth_rad, LAT)
    return np.array([R_EARTH, 0.0, v0, 0.0, vehicle.m0])


def evaluate_guidance(payload, azimuth_rad, ks, ka, kd, n_coarse=400):
    vehicle = m4_vehicle(payload)
    params = build_params(vehicle, azimuth_rad)
    y0 = initial_state(vehicle, azimuth_rad)
    t_burn = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    return find_drag_consistent_cutoff(vehicle, params, y0, ks, ka, kd, t_burn,
                                        TARGET_ALTITUDE, n_coarse=n_coarse)


def local_guidance_search(payload, azimuth_rad, seed, spans=((3, 1), (3, 2), (3, 2))):
    """Small local grid around ``seed`` = (kick_start, kick_angle, kick_duration).

    ``spans`` gives (n_points, step) per dimension, e.g. (3,1) -> seed-1,seed,seed+1.
    Returns (best_guidance, best_circ_dv), or (None, None) if nothing in the grid passes.
    """
    ks0, ka0, kd0 = seed
    (n_ks, s_ks), (n_ka, s_ka), (n_kd, s_kd) = spans
    ks_vals = [ks0 + s_ks * (j - n_ks // 2) for j in range(n_ks)]
    ka_vals = [max(ka0 + s_ka * (j - n_ka // 2), 1.0) for j in range(n_ka)]
    kd_vals = [max(kd0 + s_kd * (j - n_kd // 2), 5.0) for j in range(n_kd)]

    best = None
    for ks in ks_vals:
        for ka in ka_vals:
            for kd in kd_vals:
                res = evaluate_guidance(payload, azimuth_rad, ks, ka, kd)
                if res.found_passing_cutoff:
                    if best is None or res.insertion.circularization_delta_v < best[3]:
                        best = (ks, ka, kd, res.insertion.circularization_delta_v)
    return (None, None) if best is None else (best[:3], best[3])


# Search levels applied CONSISTENTLY to every inclination (DESIGN.md M5 S6/S13): a
# genuine bug was found and fixed here during development -- an earlier version of
# this search returned as soon as ANY span found a passing point, which gave
# inclinations whose seed already passed in the narrowest span (e.g. 28.5 deg, seeded
# from its own known-good M4 guidance) a much SHALLOWER search than inclinations whose
# narrow span found nothing and were forced to widen (e.g. 45 deg). That depth
# asymmetry, not the underlying orbital mechanics, was producing an apparent payload
# INCREASE with inclination (28.5 deg under-optimized relative to 45 deg). Fixed by
# always evaluating LEVELS 1-2 for every inclination and keeping the best result
# across both, widening further only if nothing passes at all.
SEARCH_LEVELS = [
    ((3, 1), (3, 2), (3, 2)),
    ((5, 2), (5, 4), (5, 4)),
]
WIDER_SEARCH_LEVELS = [
    ((7, 1), (7, 2), (7, 2)),
    ((7, 4), (7, 6), (7, 6)),
]


def solve_guidance_for_inclination(i_deg, seed):
    """Consistent-depth local search (levels 1-2 always; wider only if needed)."""
    az = inc.azimuth_from_inclination(math.radians(i_deg), LAT)

    overall_best = None
    overall_best_dv = None
    for spans in SEARCH_LEVELS:
        found, dv = local_guidance_search(REFERENCE_PAYLOAD, az, seed, spans)
        if found is not None and (overall_best_dv is None or dv < overall_best_dv):
            overall_best, overall_best_dv = found, dv

    if overall_best is not None:
        return az, overall_best

    for spans in WIDER_SEARCH_LEVELS:
        found, dv = local_guidance_search(REFERENCE_PAYLOAD, az, seed, spans)
        if found is not None and (overall_best_dv is None or dv < overall_best_dv):
            overall_best, overall_best_dv = found, dv
        if overall_best is not None:
            return az, overall_best

    raise RuntimeError(
        f"No passing guidance found for inclination {i_deg} deg within the widened "
        f"local search around seed {seed}. Refusing to silently report a FAIL here "
        f"without a wider, deliberate search."
    )


def find_bracket(azimuth_rad, guidance):
    ks, ka, kd = guidance
    payload = REFERENCE_PAYLOAD
    res = evaluate_guidance(payload, azimuth_rad, ks, ka, kd)
    if not res.found_passing_cutoff:
        raise RuntimeError(f"Reference payload does not pass at azimuth={math.degrees(azimuth_rad)}")
    lo, hi = payload, None
    while payload < BRACKET_MAX_PAYLOAD:
        payload += BRACKET_STEP
        res = evaluate_guidance(payload, azimuth_rad, ks, ka, kd)
        if not res.found_passing_cutoff:
            hi = payload
            break
        lo = payload
    if hi is None:
        raise RuntimeError("No failing payload found up to BRACKET_MAX_PAYLOAD")
    return lo, hi


def bisect_boundary(azimuth_rad, guidance, lo, hi):
    ks, ka, kd = guidance
    while hi - lo > BISECTION_TOLERANCE:
        mid = 0.5 * (lo + hi)
        res = evaluate_guidance(mid, azimuth_rad, ks, ka, kd)
        if res.found_passing_cutoff:
            lo = mid
        else:
            hi = mid
    return lo, hi


def full_diagnostics(payload, azimuth_rad, guidance, n_coarse=400):
    # IMPORTANT: n_coarse must match whatever resolution was used to determine
    # PASS/FAIL for this payload (evaluate_guidance's default) -- a genuine bug found
    # during development: recomputing the cutoff here with the library's own default
    # n_coarse=800 (different from the n_coarse=400 used for bisection) could, right at
    # a razor-edge payload boundary, disagree with the bisection's own PASS verdict.
    ks, ka, kd = guidance
    vehicle = m4_vehicle(payload)
    params = build_params(vehicle, azimuth_rad)
    y0 = initial_state(vehicle, azimuth_rad)
    t_burn_full = prop.burn_duration(vehicle.m_propellant, vehicle.thrust, vehicle.isp)
    cutoff = find_drag_consistent_cutoff(vehicle, params, y0, ks, ka, kd, t_burn_full,
                                          TARGET_ALTITUDE, n_coarse=n_coarse)
    t_end = cutoff.cutoff_time if cutoff.found_passing_cutoff else t_burn_full
    control = controls.gravity_turn_control(vehicle.thrust, ks, math.radians(ka), kd)
    t_eval = np.linspace(0.0, t_end, 3000)
    result = run_ascent(y0, (0.0, t_end), params, control, max_step=0.3, rtol=1e-10,
                         atol=1e-10, terminal_depletion=False, include_apogee_event=False,
                         t_eval=t_eval)
    h = result.r - R_EARTH
    v_rel = np.array([
        dyn.relative_speed(v, g, r, LAT, azimuth_rad=azimuth_rad)
        for v, g, r in zip(result.v, result.gamma, result.r)
    ])
    q = 0.5 * density(h) * v_rel**2
    return vehicle, cutoff, dict(
        t_burn_full=t_burn_full, max_q_kPa=float(np.max(q)) / 1000.0,
        burnout_time_s=result.t[-1], burnout_altitude_km=h[-1] / 1000.0,
        burnout_speed_ms=result.v[-1], burnout_gamma_deg=math.degrees(result.gamma[-1]),
        burnout_mass_kg=result.m[-1],
    )


def solve_inclination(i_deg, seed):
    az, guidance = solve_guidance_for_inclination(i_deg, seed)
    lo, hi = find_bracket(az, guidance)
    lo, hi = bisect_boundary(az, guidance, lo, hi)
    max_payload = math.floor(lo)
    just_failing = math.ceil(hi)
    assert evaluate_guidance(max_payload, az, *guidance).found_passing_cutoff
    assert not evaluate_guidance(just_failing, az, *guidance).found_passing_cutoff

    vehicle, cutoff, extra = full_diagnostics(max_payload, az, guidance)
    if not cutoff.found_passing_cutoff:
        raise RuntimeError(
            f"full_diagnostics disagreed with the bisection's PASS verdict at "
            f"payload={max_payload} kg, inclination={i_deg} deg -- resolution "
            f"mismatch between bracket/bisection and diagnostics evaluation."
        )
    useful_boost = inc.useful_rotational_boost(az, LAT)
    cross_track = inc.cross_track_rotational_component(az, LAT)
    ins = cutoff.insertion

    row = dict(
        inclination_deg=i_deg, azimuth_deg=math.degrees(az),
        useful_rotational_boost_ms=useful_boost, cross_track_rotational_ms=cross_track,
        kick_start_s=guidance[0], kick_angle_deg=guidance[1], kick_duration_s=guidance[2],
        max_payload_kg=max_payload, just_failing_payload_kg=just_failing,
        m0_kg=vehicle.m0, t_burn_full_s=extra["t_burn_full"],
        cutoff_time_s=extra["burnout_time_s"], max_q_kPa=extra["max_q_kPa"],
        cutoff_altitude_km=extra["burnout_altitude_km"],
        cutoff_speed_ms=extra["burnout_speed_ms"],
        cutoff_gamma_deg=extra["burnout_gamma_deg"],
        apogee_km=ins.apogee_altitude / 1000.0 if ins.apogee_altitude else float("nan"),
        circularization_dv_ms=ins.circularization_delta_v,
    )
    return row, guidance


def main():
    rows = []
    guidance_map = {}
    seed = SEED_GUIDANCE[28.5]
    for i_deg in AUTHORITATIVE_INCLINATIONS_DEG:
        print(f"=== solving inclination {i_deg} deg (seed {seed}) ===", flush=True)
        row, guidance = solve_inclination(i_deg, seed)
        rows.append(row)
        guidance_map[i_deg] = guidance
        seed = guidance  # seed the next inclination from this one
        print(row, flush=True)

    baseline_payload = rows[0]["max_payload_kg"]
    for row in rows:
        row["payload_loss_vs_28p5_kg"] = baseline_payload - row["max_payload_kg"]
        row["payload_loss_vs_28p5_pct"] = 100.0 * row["payload_loss_vs_28p5_kg"] / baseline_payload

    # NOTE: written under a DISTINCT filename from the committed authoritative CSV
    # (which is produced by scripts/m5_finalize_results.py and independently
    # re-verifies every row, including the 90 deg point this automatic search cannot
    # reliably reach) -- this avoids any risk of silently overwriting the committed
    # results with a re-run of the (slower, automatic-search-only, non-polar) method.
    out_path = os.path.join(os.path.dirname(__file__),
                             "m5_inclination_sweep_automatic_search_28p5_to_70.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {len(rows)} authoritative rows to {out_path}")

    # Fine sweep for the headline curve: seed from the nearest authoritative
    # inclination's solved guidance. Uses a CHEAP single-level (27-combination) local
    # search only (not the full consistent-depth search used for the authoritative
    # table) to keep runtime practical -- this is exactly the "coarse re-optimization
    # + seed/interpolate neighboring inclinations" policy the task allows when full
    # re-optimization at every point is prohibitive. If even that finds nothing, the
    # nearest authoritative guidance is reused AS-IS (still bracket+bisected fresh at
    # this inclination, so the reported payload is not fabricated) rather than
    # widening further, since the authoritative table already carries the
    # fully-searched result at that neighboring point.
    #
    # NOTE: this exploratory fine sweep was dropped from the final M5 deliverable
    # (DESIGN.md M5 S13) once the near-polar search demonstrated the same severe cost
    # would apply to any inclination near 90 deg; the code is left here as a runnable
    # record of the approach, not as part of the final data pipeline.
    fine_rows = []
    auth_sorted = sorted(guidance_map.keys())
    for i_deg in FINE_INCLINATIONS_DEG:
        nearest = min(auth_sorted, key=lambda a: abs(a - i_deg))
        seed = guidance_map[nearest]
        az = inc.azimuth_from_inclination(math.radians(i_deg), LAT)
        found, _ = local_guidance_search(REFERENCE_PAYLOAD, az, seed, SEARCH_LEVELS[0])
        guidance = found if found is not None else seed
        try:
            lo, hi = find_bracket(az, guidance)
            lo, hi = bisect_boundary(az, guidance, lo, hi)
            max_payload = math.floor(lo)
            useful_boost = inc.useful_rotational_boost(az, LAT)
            fine_rows.append(dict(inclination_deg=i_deg, azimuth_deg=math.degrees(az),
                                   useful_rotational_boost_ms=useful_boost,
                                   max_payload_kg=max_payload,
                                   kick_start_s=guidance[0], kick_angle_deg=guidance[1],
                                   kick_duration_s=guidance[2]))
            print(f"fine i={i_deg}: max_payload={max_payload} kg (guidance {guidance})", flush=True)
        except RuntimeError as e:
            print(f"fine i={i_deg}: FAILED to solve -- {e}", flush=True)

    fine_path = os.path.join(os.path.dirname(__file__), "m5_inclination_fine_results.csv")
    with open(fine_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fine_rows[0].keys()))
        writer.writeheader()
        writer.writerows(fine_rows)
    print(f"wrote {len(fine_rows)} fine rows to {fine_path}")

    return rows, fine_rows


if __name__ == "__main__":
    main()
