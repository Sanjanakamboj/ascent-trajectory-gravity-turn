# Ascent Trajectory / Gravity Turn

Launch-vehicle ascent simulation portfolio project. Goal: produce (1) an ascent
trajectory profile for a single-stage-baseline launch vehicle flying a gravity-turn
ascent to a target circular LEO, and (2) a payload-to-orbit trade curve versus target
orbital inclination.

This project is built **milestone-by-milestone**. Each milestone is committed and pushed
individually; no milestone is started before the previous one is checkpointed and
approved.

## Status

- [x] **M1 — Mission definition, equations, conventions, hand calculations, verification
      plan, limitations.** See [`DESIGN.md`](DESIGN.md). No integrator yet.
- [x] **M2 — Atmosphere + propulsion + point-mass ascent ODE + basic trajectory
      verification.** See [`DESIGN.md` §12](DESIGN.md#12-milestone-2--physics-implementation-and-verification).
      Physics verification only, no guidance law; baseline vehicle does **not** reach
      400 km circular LEO with a fixed pitch profile (expected — see §12.9).
- [x] **M3 — Gravity-turn guidance, full trajectory profile, event handling, numerical
      convergence.** See [`DESIGN.md` §13](DESIGN.md#13-milestone-3--gravity-turn-guidance-orbital-diagnostics-and-one-verified-trajectory).
      Verified guidance law + 25-case parameter sweep on the unchanged vehicle; 0/25
      cases achieve 400 km circular orbit (expected, given the M1 §7.4 Δv deficit).
- [x] **M4 — Orbit-capable study vehicle + payload-to-orbit solve at 400 km / 28.5°.**
      See [`DESIGN.md` §14](DESIGN.md#14-milestone-4--orbit-capable-study-vehicle-and-payload-to-orbit-solve).
      A new, explicitly separate **M4 study vehicle** (higher Isp, higher propellant
      fraction — NOT the M1–M3 verification vehicle, which remains a documented FAILURE
      under the M4 criterion) achieves 400 km circular insertion (ascent + idealized
      circularization) up to **10,333 kg payload**.
- [ ] M5 — Inclination sweep, launch-azimuth/Earth-rotation coupling, payload-vs-
      inclination curve.
- [ ] M6 — Independent validation, sensitivity analysis, figure audit, packaging.

## Important note on figures

Any figure produced before the trajectory model is validated (M3's convergence study and
the M9 checks in `DESIGN.md`) is **diagnostic/supporting only**, not a validated result.
Figures will be explicitly labeled as such until the verification plan in `DESIGN.md`
§9 has been executed.

## Repository layout

```
DESIGN.md       mission definition, equations, conventions, hand calcs, verification
                plan, limitations, M2-M4 implementation notes (start here)
src/ascent/     simulation package
    constants.py       Earth constants, M1-M3 baseline vehicle, M4 study vehicle
    atmosphere.py       simplified exponential density model
    propulsion.py       mass-flow / thrust / Tsiolkovsky helpers
    dynamics.py         planar point-mass ascent ODE + Earth-rotation handling
    controls.py         gravity-turn guidance law + prescribed thrust-direction profiles
    simulation.py       solve_ivp wrapper + event handling
    orbital.py          orbital-element diagnostics + orbit-insertion criteria
    insertion_search.py M4 drag-consistent engine-cutoff search
    loss_budget.py       M4 delta-v / loss accounting
tests/          pytest test suite
scripts/        m2_diagnostic_trajectory.py, m3_gravity_turn_sweep.py, m3_trajectory.py,
                m4_guidance_search.py, m4_payload_sweep.py, m4_trajectory.py,
                m4_payload_figure.py
figures/        generated figures (diagnostic/supporting unless noted as validated)
```

### M4 orbit-capable study vehicle — payload capability (current headline result)

![M4 payload capability](figures/m4_payload_capability.png)

**This is a different vehicle from the M1–M3 verification vehicle** (higher Isp, 450 s
vs. 300 s; higher propellant fraction — see `DESIGN.md` §14.3), explicitly labeled the
"M4 orbit-capable study vehicle." Using a retuned gravity-turn guidance law and a
standard ascent + idealized-circularization insertion (`DESIGN.md` §14.2), it achieves
400 km circular orbit for payloads up to **10,333 kg**; the original M1–M3 vehicle
remains a confirmed failure under the same criterion (`DESIGN.md` §14.9 check A).

![M4 orbit-capable trajectory](figures/m4_orbit_capable_trajectory.png)

### M3 gravity-turn trajectory (unchanged verification vehicle — does not reach orbit)

![M3 gravity-turn trajectory](figures/m3_gravity_turn_trajectory.png)

Verified gravity-turn guidance law (vertical rise -> pitch-kick -> zero-AoA turn) on the
unchanged M1/M2 baseline vehicle, selected from a 25-case parameter sweep by maximum
burnout specific orbital energy — see `DESIGN.md` §13. **Diagnostic/supporting only,
not a validated result.** No case in the sweep reaches 400 km circular orbit, which is
the expected outcome given the Δv deficit documented in `DESIGN.md` §7.4.

### M2 diagnostic trajectory (earlier, prescribed-control-only baseline)

![M2 diagnostic trajectory](figures/m2_diagnostic_trajectory.png)

Fixed near-vertical pitch profile (no guidance law) on the unchanged M1 baseline
vehicle — see `DESIGN.md` §12.9. Diagnostic/supporting only, kept for comparison
against the M3 result (`DESIGN.md` §13.7).

## Development

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -W error
```

## Scope and limitations

See `DESIGN.md` §10 for the full, explicit list of modeling limitations (point-mass only,
no 6-DOF, no winds, simplified atmosphere, constant Cd, no structural/load model, no
throttling, no staging unless later added, no range constraints, no operational
flight-safety claim).
