"""M6 canonical final summary: prints a single consolidated report of the project's
headline results, reading ONLY from already-committed, authoritative CSV outputs and
already-existing package APIs. Duplicates no trajectory equations and re-solves no
guidance search -- this is a reporting script, not another physics computation.

Usage:
    python scripts/final_portfolio_summary.py
"""

import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")


def _read_csv(name):
    path = os.path.join(SCRIPTS, name)
    with open(path) as f:
        return list(csv.DictReader(f))


def section(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def m3_summary():
    section("M3 — verification vehicle, gravity-turn guidance sweep (400 km / 28.5 deg target)")
    rows = _read_csv("m3_sweep_results.csv")
    achieved = [r for r in rows if r["orbit_achieved"] == "True"]
    print(f"cases swept: {len(rows)}")
    print(f"cases achieving 400 km circular insertion: {len(achieved)} / {len(rows)}")
    print("Conclusion: the M1-M3 baseline vehicle (Isp=300s) is USEFUL for dynamics/physics")
    print("verification but does NOT reach the 400 km target under any tested gravity-turn")
    print("guidance -- expected, given the Delta-v deficit documented in DESIGN.md S7.4. It")
    print("is not, and must not be presented as, an orbit-capable launch vehicle.")


def m4_summary():
    section("M4 — orbit-capable study vehicle, payload-to-orbit @ 400 km / 28.5 deg")
    rows = _read_csv("m4_payload_sweep_results.csv")
    passing = [r for r in rows if r["orbit_achieved"] == "True"]
    best_pass = max(passing, key=lambda r: float(r["payload_kg"]))
    # "First failing" = the lightest failing case that is heavier than the best pass
    # (i.e. the FAIL side of the boundary bracket), not just any low-payload failure
    # (the M4 sweep also includes documented low-payload failures, e.g. 0/5000 kg,
    # from a separate non-monotonicity check -- DESIGN.md M4 S8 check B).
    heavier_failing = [r for r in rows if r["orbit_achieved"] == "False"
                        and float(r["payload_kg"]) > float(best_pass["payload_kg"])]
    first_fail = min(heavier_failing, key=lambda r: float(r["payload_kg"]))
    print(f"cases in the CSV: {len(rows)}")
    print(f"maximum passing payload:   {float(best_pass['payload_kg']):.0f} kg "
          f"(apogee {float(best_pass['apogee_altitude_km']):.3f} km, "
          f"circularization dv {float(best_pass['circularization_delta_v_ms']):.2f} m/s)")
    print(f"first failing payload:     {float(first_fail['payload_kg']):.0f} kg "
          f"(apogee miss {float(first_fail['closest_apogee_miss_km']):.1f} km)")
    print("bisection tolerance: 1.0 kg (script-level); reported headline figure rounds to")
    print("the nearest 10-100 kg -- see DESIGN.md M4 S9 / M6 audit for resolution statement.")
    print("guidance used: kick_start=45s, kick_angle=20deg, kick_duration=20s (fixed, M4 S4)")
    print("insertion criterion: bound + non-Earth-intersecting + apogee within 15 km of")
    print("400 km target + idealized circularization dv <= 1500 m/s (DESIGN.md M4 S2)")


def m4_loss_budget():
    section("M4 — delta-v / loss budget (reference payload = 10,000 kg case)")
    print("Definitions (DESIGN.md M4 S6): achieved_dv = ideal_dv_to_cutoff - steering_loss")
    print("  - drag_loss - gravity_loss, all evaluated over the SAME powered-flight interval")
    print("  to the drag-consistent engine cutoff. This is an internal-consistency budget for")
    print("  this vehicle/guidance/model, not a universal launch-vehicle delta-v accounting.")
    print("See figures/m4_delta_v_budget.png and DESIGN.md M4 S6 for the numerical breakdown.")


def m5_summary():
    section("M5 — direct-ascent payload vs. target inclination (28.5-90 deg)")
    rows = _read_csv("m5_inclination_sweep_results.csv")
    rows.sort(key=lambda r: float(r["inclination_deg"]))
    print(f"{'i(deg)':>7} {'Az(deg)':>8} {'v_boost(m/s)':>13} {'v_cross(m/s)':>13} "
          f"{'payload(kg)':>12} {'circ dv(m/s)':>13}")
    for r in rows:
        print(f"{float(r['inclination_deg']):7.1f} {float(r['azimuth_deg']):8.2f} "
              f"{float(r['useful_rotational_boost_ms']):13.1f} "
              f"{float(r['cross_track_rotational_ms']):13.1f} "
              f"{float(r['max_payload_kg']):12.0f} "
              f"{float(r['circularization_dv_ms']):13.2f}")

    payloads = [float(r["max_payload_kg"]) for r in rows]
    baseline = next(float(r["max_payload_kg"]) for r in rows if float(r["inclination_deg"]) == 28.5)
    polar = next(float(r["max_payload_kg"]) for r in rows if float(r["inclination_deg"]) == 90.0)
    boosts = [float(r["useful_rotational_boost_ms"]) for r in rows]

    print(f"\npayload range across the 6 authoritative inclinations: "
          f"{min(payloads):.0f}-{max(payloads):.0f} kg (spread "
          f"{100*(max(payloads)-min(payloads))/baseline:.1f}% of baseline)")
    print(f"28.5 deg baseline payload: {baseline:.0f} kg")
    print(f"90 deg (polar) payload:    {polar:.0f} kg "
          f"({'+' if polar >= baseline else ''}{polar - baseline:.0f} kg, "
          f"{'+' if polar >= baseline else ''}{100*(polar-baseline)/baseline:.2f}%)")
    monotonic_boost = all(boosts[i] > boosts[i + 1] for i in range(len(boosts) - 1))
    monotonic_payload = all(payloads[i] >= payloads[i + 1] for i in range(len(payloads) - 1))
    print(f"useful rotational boost strictly monotonic decreasing with inclination: "
          f"{monotonic_boost}")
    print(f"maximum payload monotonic (non-increasing) with inclination: {monotonic_payload}")
    print("Interpretation (DESIGN.md M5 S11, M6 audit): the geometric rotational-assistance")
    print("effect IS monotonic; the SEARCH-OPTIMIZED payload is NOT. For this vehicle, the")
    print("payload differences across inclination are small enough that the guidance-search")
    print("landscape's ruggedness competes with the (smaller) pure rotational-assistance")
    print("effect. The 90 deg case is a direct-ascent inclination effect, NOT a plane-change")
    print("delta-v penalty -- no plane change is modeled anywhere in this project.")


def solver_test_summary():
    section("Verification / test summary")
    print("Run `pytest -W error` for the authoritative, current count.")
    print("As of this M6 audit: 119 tests passing (81 from M1-M4 unchanged, 38 from M5),")
    print("covering: unit/analytical checks (M1), integration/event-handling checks (M2-M3),")
    print("mass-bookkeeping and payload-boundary checks (M4-M5), inclination geometry and")
    print("atmosphere-relative-wind regression tests (M5), and integrator convergence checks")
    print("at 28.5/55/90 deg (M5).")


def main():
    m3_summary()
    m4_summary()
    m4_loss_budget()
    m5_summary()
    solver_test_summary()
    print()


if __name__ == "__main__":
    main()
