# DESIGN.md — Ascent Trajectory / Gravity Turn

**Milestone 1: mission definition, governing equations, conventions, hand calculations,
verification plan, and limitations.**

> Status: design/derivation only. No integrator, no ODE solver, no trajectory code exists
> yet. Everything numerical in this document is a *hand calculation* used to sanity-check
> the scenario, not a simulation result. Any plot produced before the M3/M4 solver is
> validated must be labeled **diagnostic/supporting**, not proof of vehicle performance.

---

## 1. Mission scenario

A single-stage-baseline expendable launch vehicle ("Ascent-1") lifts off from a fixed
launch site, executes a vertical rise + pitch-kick + gravity-turn ascent through a simple
atmosphere, and targets insertion into a circular low Earth orbit (LEO). The scenario is
deliberately sized to exercise every physical effect this portfolio project needs:
gravity-turn steering, aerodynamic drag, propellant mass depletion, a single well-defined
propulsion stage (extensible to staging later), a concrete circular-orbit target, and a
payload-vs-inclination trade driven by launch-site latitude and Earth rotation.

### 1.1 Launch site

| Quantity | Value | Notes |
|---|---|---|
| Site | Cape Canaveral–like | representative, not an operational range model |
| Latitude, `lat` | 28.5° N | drives minimum direct-ascent inclination |
| Longitude | not modeled | planar (in-plane) ascent only in M1–M4; azimuth/geometry handled analytically for the inclination trade |

### 1.2 Earth constants

| Symbol | Value | Units | Notes |
|---|---|---|---|
| `mu` (GM) | 3.986004418e14 | m³/s² | standard gravitational parameter |
| `R_earth` | 6,378,137 | m | WGS84 equatorial radius, used as spherical Earth radius (no oblateness/J2 in M1–M5) |
| `omega_earth` | 7.2921159e-5 | rad/s | Earth sidereal rotation rate |
| `g0` | 9.80665 | m/s² | standard gravity, used only to convert Isp → exhaust velocity |

### 1.3 Target orbit

| Quantity | Value |
|---|---|
| Target altitude (circular) | 400 km |
| Target radius `r_target` | 6,778,137 m |
| Baseline target inclination | 28.5° (= launch latitude; due-east launch, minimum-inclination direct ascent) |
| Later milestone (M5) | inclination swept over a range including values > latitude (dogleg/plane-change cost not modeled — see Limitations) and retrograde |

### 1.4 Launch azimuth convention

Launch azimuth `Az` is measured **clockwise from true north** (0° = north, 90° = east,
180° = south, 270° = west) — standard aeronautical/range convention. Baseline case:
`Az = 90°` (due east), which gives the minimum reachable direct-ascent inclination
(`i_min = lat`) and the maximum possible assist from Earth's rotation.

### 1.5 Vehicle: "Ascent-1" (single-stage baseline)

Sized as a **numerical-methods testbed**, not a claim that this exact mass budget reaches
orbit (see §7.4 and Limitations — single-stage-to-orbit with chemical propulsion at
Isp = 300 s is known to be extremely mass-fraction-constrained; that tension is expected
and is exactly what M4's payload-solve milestone will quantify).

| Quantity | Symbol | Value | Units |
|---|---|---|---|
| Initial (gross liftoff) mass | `m0` | 500,000 | kg |
| Propellant mass | `m_p` | 410,000 | kg |
| Dry mass (structure + engines, excl. payload) | `m_dry` | 80,000 | kg |
| Payload mass (baseline; trade variable in M4/M5) | `m_pay` | 10,000 | kg |
| Consistency check | `m_dry + m_p + m_pay = m0` | 80,000 + 410,000 + 10,000 = 500,000 ✓ | kg |
| Thrust (constant, single stage) | `T` | 7.6e6 | N |
| Specific impulse (constant, altitude-averaged) | `Isp` | 300 | s |
| Reference (frontal) area | `A_ref` | 10.75 (d = 3.7 m) | m² |
| Drag coefficient (constant) | `Cd` | 0.3 | – |

Constant `Isp` and constant `Cd` are explicit simplifications (see Limitations); M2 may
later introduce altitude-dependent Isp (sea-level ↔ vacuum) and Mach-dependent Cd without
changing the M1 equations' structure.

### 1.6 Atmosphere model (planned for M2, not implemented yet)

US Standard Atmosphere 1976 — tabulated/exponential density and pressure vs. altitude,
`rho(h)`, `P(h)`. M1 only reserves its role in the drag term; no atmosphere code exists
yet.

### 1.7 Pitch-over / gravity-turn control law (concept only)

1. **Vertical rise**: from liftoff to a small altitude/time (e.g. until clear of the
   tower, ~ a few seconds), thrust held vertical, `γ = 90°`.
2. **Pitch-kick**: a brief, small commanded tilt of the thrust vector (a few degrees) away
   from vertical, in the plane defined by the launch azimuth, to seed a non-zero
   horizontal velocity component. This is the only phase where thrust is *not* aligned
   with the velocity vector.
3. **Gravity turn**: thrust held aligned with the instantaneous velocity vector
   (zero angle-of-attack) for the remainder of powered flight; the trajectory's curvature
   from vertical to horizontal is then driven passively by gravity, not by continued
   active steering. This is the standard "gravity turn" approximation used throughout the
   ascent literature.

### 1.8 Engine cutoff (MECO) condition

Baseline cutoff is **propellant depletion**: burn ends when `m_p` has been consumed
(`t_burn = m_p / mdot`, §7). A later milestone may add an alternate cutoff
(target speed/altitude reached) as a solver stopping condition, but the two must not be
silently conflated — M1 defines depletion as the single authoritative baseline condition.

---

## 2. Reference frame

Ascent dynamics are formulated as a **planar, two-body point-mass model in the
Earth-Centered Inertial (ECI) frame**, restricted to the orbital plane defined by the
launch site's position and inertial velocity vector at liftoff (i.e., the plane fixed by
the chosen launch azimuth). Polar coordinates `(r, θ)` are used within that plane, origin
at Earth's center.

- `r` — radial distance from Earth's center (`h = r - R_earth` is altitude).
- `θ` — downrange central angle, measured in the inertial orbital plane from the launch
  point's initial position vector, increasing in the direction of flight.
- `v` — **inertial** speed (ECI), magnitude of the velocity vector.
- `γ` — flight-path angle, the angle between the velocity vector and the local horizontal
  (perpendicular to the local radius vector), **positive climbing (up)**.
- `m` — instantaneous vehicle mass.

Earth's rotation enters this inertial-frame formulation entirely through the **initial
condition**: at liftoff the vehicle is already co-rotating with the launch site, so its
inertial velocity is not zero but `v0 = omega_earth * R_earth * cos(lat)`, directed due
east (tangent to the latitude circle), before the engines light. This is what allows an
easterly launch to "keep" part of Earth's rotational speed toward orbital insertion, and
is why the equations of motion below are correct for insertion checks (which are inertial
by definition) while the drag term below needs an explicit correction (§3, Note on
relative wind).

---

## 3. Governing equations (planar point-mass ascent)

**Kinematics:**
```
dr/dt     = v sin(γ)
dθ/dt     = v cos(γ) / r
```

**Dynamics** (thrust aligned with velocity vector, i.e. zero angle-of-attack gravity turn;
during the brief pitch-kick phase only, thrust direction is offset from velocity by the
small commanded kick angle — see §1.7):
```
dv/dt     = (T - D) / m  -  g(r) sin(γ)
v dγ/dt   = [ v²/r  -  g(r) ] cos(γ)
```

**Mass depletion** (constant thrust / constant Isp, while propellant remains):
```
dm/dt     = - T / (Isp * g0)          (0 once m_p is exhausted; T → 0 at MECO)
```

**Local gravity** (spherical Earth, no J2 in M1–M5):
```
g(r)      = mu / r²                    (directed toward Earth's center, i.e. -r̂)
```

**Aerodynamic drag:**
```
D         = 1/2 * rho(h) * v_rel² * Cd * A_ref
```

**Note on relative wind (important, and a known source of ascent-code bugs — see
Verification Plan):** the atmosphere co-rotates with the Earth. Drag must be computed
using the vehicle's speed **relative to the rotating atmosphere**,
`v_rel = v - v_atm(r, lat)` (vector subtraction of the local atmospheric co-rotation
velocity from the inertial velocity), **not** the raw inertial speed `v` used everywhere
else in the equations above. `v` (inertial) is what matters for orbital energy and the
circular-orbit insertion check; `v_rel` (rotating-frame/relative) is what matters for
drag. M2 must implement both explicitly and keep them distinct.

**Altitude:**
```
h = r - R_earth
```

---

## 4. Sign / angle conventions

| Convention | Definition |
|---|---|
| Flight-path angle `γ` | angle between velocity vector and local horizontal; **positive = climbing**. `γ ≈ 90°` at liftoff (nominally vertical, though `v ≈ 0` there so direction is set by the vertical-rise phase, not by `γ` itself); `γ → 0°` at circular-orbit insertion (horizontal flight). |
| Pitch / thrust direction | equal to the velocity-vector direction (zero angle-of-attack) during gravity turn; offset from velocity by the small commanded kick angle only during the pitch-kick phase (§1.7). |
| Launch azimuth `Az` | measured **clockwise from true north**: 0° = N, 90° = E, 180° = S, 270° = W. |
| Heading | measured the same way as azimuth — **clockwise from north** (not counterclockwise from east/math convention). |
| Orbital inclination `i` | standard range 0°–180°, angle between the orbital plane and Earth's equatorial plane. `i = 0°`: equatorial prograde. `i = 90°`: polar. `i = 180°`: equatorial retrograde. `i < 90°` = prograde (same rotational sense as Earth); `i > 90°` = retrograde. |
| Inertial vs. rotating-Earth velocity | `v` (state variable, used in all dynamics/energy/insertion equations) is **inertial (ECI)**. `v_rel` (drag only) is speed **relative to the rotating atmosphere**. These are explicitly different quantities and must never be interchanged (§3 Note). |
| Downrange angle `θ` | increases in the direction of flight within the inertial orbital plane, measured from the launch point's initial position vector. |

---

## 5. Launch-site inclination constraint

Spherical trigonometry relates launch azimuth, launch latitude, and the resulting orbital
inclination (flat/point-mass geometry, no dogleg or plane-change maneuver):

```
cos(i) = cos(lat) * sin(Az)
```

**Consequences (baseline latitude = 28.5°):**

- **Minimum directly reachable prograde inclination**: `i_min = lat = 28.5°`, achieved at
  `Az = 90°` (due east). No smaller inclination is directly reachable from this site
  without an explicit plane-change/dogleg maneuver (not modeled — see Limitations).
- **Retrograde is reachable**: e.g. `Az = 270°` (due west) gives
  `cos(i) = -cos(lat)` → `i = 180° - lat = 151.5°`. Physically, a westward launch must
  first cancel the eastward co-rotation velocity `v0` and then build speed in the opposite
  sense — this is the most rotation-penalized case.
- **Polar orbit** (`i = 90°`): requires `sin(Az) = 0`, i.e. `Az = 0°` or `180°` (due
  north/south). Here the co-rotation velocity `v0` is perpendicular to the direction of
  travel — it neither helps nor directly opposes the burn, but it is *not* zero, and a
  correct polar-orbit insertion must still account for it as a nonzero orthogonal velocity
  component in the final inertial-velocity vector sum (otherwise the achieved inclination
  drifts slightly off 90°).
- **Earth-rotation boost vs. azimuth/inclination**: the usable rotational assist along the
  direction of flight scales like `v_rot * sin(Az) = v_rot * cos(i)/cos(lat) `... more
  directly, the assist component along the flight direction is maximal (`+v_rot`) at
  `Az = 90°` (due east, `i = lat`), zero at `Az = 0°/180°` (polar), and works maximally
  against the vehicle (`-v_rot`, requiring the largest extra Δv) at `Az = 270°` (due west,
  `i = 180° - lat`). This azimuth/inclination-dependent Δv coupling is exactly what M5's
  inclination sweep and payload-vs-inclination curve will quantify; **M1 establishes the
  physics only and does not optimize or sweep it yet.**

---

## 6. Angle/velocity conventions used for the calculations in §7

- Baseline case uses `Az = 90°` (due east) so that `i = lat` and the rotational velocity
  `v0` is exactly aligned (parallel) with the required orbital velocity direction — this
  is the simplest, most favorable case and is why a plain scalar subtraction
  (`v_circ - v_rot`) is valid for the §7 hand estimate. For any other azimuth this becomes
  a vector sum, deferred to M5.

---

## 7. Hand calculations (baseline case, `Az = 90°`, `lat = 28.5°`, `h_target = 400 km`)

All values computed directly from the constants and vehicle parameters in §1
(script used: ad-hoc, not committed — these are one-time sanity numbers, not test
fixtures).

| Quantity | Formula | Result |
|---|---|---|
| Surface rotational velocity | `v_rot = omega_earth * R_earth * cos(lat)` | **408.7 m/s** |
| Target circular-orbit speed (400 km) | `v_circ = sqrt(mu / r_target)` | **7668.6 m/s** |
| Exhaust velocity | `ve = Isp * g0` | 2942.0 m/s |
| Mass flow rate | `mdot = T / ve` | 2583.3 kg/s |
| Burn duration | `t_burn = m_p / mdot` | **158.7 s** |
| Mass ratio | `m0 / (m0 - m_p)` = `500,000 / 90,000` | 5.556 |
| Ideal Δv (Tsiolkovsky) | `ve * ln(mass_ratio)` | **5044.9 m/s** |
| Thrust-to-weight at liftoff | `T / (m0 * g0)` | **1.550** (> 1, liftoff is possible ✓) |
| Rough vehicle Δv needed, no losses, due-east | `v_circ - v_rot` | 7259.8 m/s |
| Rough total Δv budget incl. losses | `(v_circ - v_rot) + gravity-loss(~1500) + drag-loss(~150)` | **≈ 8910 m/s** (typical order for real LEO ascents is ~9.3–9.5 km/s including a proper gravity-turn loss profile; ~8.9 km/s is a reasonable rough floor) |
| Dynamic-pressure sanity estimate (representative max-Q-like point, not simulated) | `q = 1/2 * rho * v_rel²` at `rho ≈ 0.4 kg/m³` (~11–13 km alt, US Std Atm), `v_rel ≈ 450 m/s` | **≈ 40.5 kPa** (≈ 846 psf — same order of magnitude as real vehicles' max-Q, e.g. Falcon 9 ~30 kPa / 620 psf) |

### 7.1 Interpretation — thrust-to-weight and burn

`T/W0 = 1.55` confirms the vehicle can lift off (`T > W0`) with reasonable margin for
gravity losses during the initial low-speed phase. `t_burn ≈ 159 s` is a plausible
single-stage burn duration for this thrust class.

### 7.2 Interpretation — dynamic pressure

The order-of-magnitude match (~40 kPa hand estimate vs. ~30 kPa for a real vehicle at
max-Q) is a reasonable sanity check on `Cd`, `A_ref`, and the assumed density/speed at a
representative mid-ascent point. This is **not** a simulated max-Q — no trajectory has
been integrated yet — it only checks that the drag-term inputs are of a physically
sensible size before any ODE code is written.

### 7.3 Consistency — Tsiolkovsky vs. mass budget

`m0/(m0-m_p) = 500,000/90,000 = 5.556` uses exactly the masses declared in §1.5
(`m0 = 500,000`, `m_p = 410,000`, `m_dry+m_pay = 90,000`), so the Δv figure is internally
consistent with the declared mass budget — no separate/undocumented mass assumption was
used.

### 7.4 Honest gap: ideal Δv vs. required Δv budget

The vehicle's ideal (vacuum, single-burn) Tsiolkovsky Δv is **5045 m/s**, while a rough
total Δv budget to close 400 km circular LEO from this site (due east, including
approximate gravity and drag losses) is **≈ 8900–9500 m/s**. The baseline single-stage
vehicle as sized in §1.5 does **not** close this budget on paper. This is expected and is
**not** a scenario-definition error: single-stage-to-orbit with `Isp = 300 s` chemical
propulsion is a well-known hard problem (mass-ratio requirements are extreme), which is
precisely why real orbital vehicles stage. Per the task requirement, §1.5 defines "at
minimum a clearly defined single-stage baseline that can later be extended" — this
baseline is intentionally a **numerical-methods testbed** for the gravity-turn ODE, drag
model, and mass depletion, not a claim that it reaches orbit unmodified. M4's
payload-to-orbit solve will determine, from the actual integrated trajectory (not this
rough scalar estimate), what payload mass (if any, at this propellant/thrust sizing)
closes the achieved orbit — and may instead motivate an explicit staged extension before
M4, which will be decided at that checkpoint, not silently assumed now.

---

## 8. Future milestone architecture

| Milestone | Scope |
|---|---|
| **M1** (this document) | Mission definition, equations, conventions, hand calculations, verification plan, assumptions. No integrator. |
| **M2** | Atmosphere model (US Standard Atmosphere 1976), propulsion model, point-mass ascent ODE implementation, basic vertical/pitched-trajectory verification (no full gravity-turn guidance logic yet). |
| **M3** | Gravity-turn guidance (pitch-kick + zero-AoA turn), full trajectory profile, MECO/event handling, numerical convergence study (timestep/tolerance). |
| **M4** | Payload-to-orbit solve at the single baseline inclination (28.5°) — determine achievable payload mass (or vehicle sizing implications per §7.4) via the validated integrator. |
| **M5** | Inclination sweep, launch-azimuth/Earth-rotation Δv coupling (§5), payload-vs-inclination trade curve. |
| **M6** | Independent validation, sensitivity analysis, figure audit (diagnostic vs. validated), portfolio packaging. |

No milestone beyond M1 is started in this commit.

---

## 9. Verification plan (to be executed starting M2/M3, defined now)

1. **Zero-drag limit** — set `Cd = 0` (or `rho ≡ 0`): trajectory must match a drag-free
   powered-ascent case; comparison against an independent zero-drag integration or
   analytic powered-flight segment.
2. **Zero-thrust ballistic limit** — set `T = 0` after some initial velocity: the point
   mass must follow a Keplerian conic (ellipse) exactly; check against the vis-viva
   equation and conservation of specific orbital energy/angular momentum.
3. **Constant-thrust mass-flow check** — verify `m(t) = m0 - mdot*t` exactly matches
   direct integration of `dm/dt = -T/(Isp*g0)` for the duration `T` is nonzero, and that
   mass strictly stops decreasing at `t_burn` (no propellant "debt").
4. **Tsiolkovsky consistency** — with drag and gravity artificially disabled and thrust
   held exactly antiparallel... i.e. in the pure 1-D/no-gravity-loss idealization, the
   integrator's Δv over the burn must reproduce `ve * ln(m0/mf)` to within integrator
   tolerance.
5. **Circular-orbit speed check** — at the target radius `r_target`, confirm
   `v_circ = sqrt(mu/r_target)` is reproduced by the code's own orbit-mechanics helper,
   independent of the ascent integration, and that at MECO (if orbit is achieved) `v ≈
   v_circ` and `γ ≈ 0` simultaneously (both conditions, not just speed).
6. **Energy / angular-momentum sanity checks** — after MECO (thrust off), specific orbital
   energy `epsilon = v²/2 - mu/r` and angular momentum `h_ang = r*v*cos(γ)` must be
   constant along the coast (no thrust, no drag) to within numerical tolerance.
7. **Integrator timestep/tolerance convergence** — repeat the same case at halving
   timesteps (or tightening adaptive tolerance) and confirm state-vector convergence
   (e.g. Richardson extrapolation or simple convergence-order check) before trusting any
   trajectory as more than diagnostic.
8. **Atmosphere limiting behavior** — confirm `rho(h) → 0` smoothly/monotonically as
   `h → infinity` (or beyond the model's tabulated range) and that drag vanishes
   accordingly; confirm `rho(0)` matches the standard sea-level value used as the model's
   reference.
9. **Earth-rotation contribution cross-check** — independently verify
   `v0 = omega_earth * R_earth * cos(lat)` against the code's initial-condition setup, and
   confirm that a due-east baseline run's inertial initial velocity is nonzero and equal
   to this value (a common bug is starting the inertial-frame integration at `v=0`).
10. **Independent simplified ascent comparison** — cross-check the full integrator against
    a simplified closed-form or 1-D vertical/flat-Earth powered-ascent approximation over
    a short early-flight window where the simplification is expected to hold, as an
    independent sanity check outside the point-mass polar-coordinate formulation itself.

---

## 10. Model limitations (stated honestly, up front)

- **Point-mass approximation** — no vehicle geometry, no center-of-mass/center-of-pressure
  offset, no bending modes.
- **No 6-DOF attitude dynamics** — attitude is not simulated; thrust direction is
  prescribed analytically (vertical rise → pitch-kick → velocity-aligned), not the result
  of a flight-control-system/attitude-dynamics model.
- **No winds** — atmosphere is static relative to the rotating Earth; no gusts, shear, or
  weather.
- **Simplified atmosphere** — US Standard Atmosphere 1976 (planned, M2) is a 1-D,
  time-invariant, latitude/season-independent model.
- **Constant/simplified `Cd`** — a single scalar drag coefficient (`Cd = 0.3`) is used
  across the whole ascent; no true Mach-number-dependent drag rise (e.g. transonic drag
  bump) is modeled unless a later milestone explicitly adds it.
- **No structural/bending/load model** — no max-Q structural constraint enforcement, no
  loads analysis beyond the dynamic-pressure sanity number in §7.
- **No engine throttling** — constant thrust and constant Isp unless a later milestone
  explicitly implements throttling.
- **No staging** — single-stage baseline only, unless and until a later milestone
  explicitly adds a staging event (see §7.4 for why this may become necessary).
- **No launch-range constraints** — no downrange safety corridor, no abort modes, no
  range-safety destruct logic.
- **No operational flight-safety claim** — this is a portfolio/engineering-methods
  project; it must not be read as validated for, or representative of, any real flight
  safety, licensing, or operational decision.

---

## 11. Verification/derivation cross-reference

Every equation in §3 is used by at least one check in §9 (drag term → checks 1 & 8;
thrust/mass term → checks 3 & 4; gravity/orbital-energy term → checks 2, 5 & 6; Earth
rotation initial condition → check 9). No equation is introduced in §3 without a
corresponding verification path, and no verification check depends on an equation not
derived in §3.

---

## 12. Milestone 2 — physics implementation and verification

> **M2 is physics verification, not an optimized ascent.** No gravity-turn guidance law
> (velocity-following, zero angle-of-attack after a pitch-kick) is implemented here —
> that is explicitly M3 scope. All prescribed controls below are fixed/open-loop and
> were not tuned to reach orbit or to improve any performance metric.

### 12.1 M1 re-audit (before any code was written)

Before implementing M2, the M1 hand calculations were independently recomputed from the
constants and vehicle parameters in §1. All values matched exactly
(`v_circ = 7668.6 m/s`, `v_rot = 408.7 m/s`, `mdot = 2583.3 kg/s`, `t_burn = 158.7 s`,
`mass_ratio = 5.556`, `dv_ideal = 5044.9 m/s`, `T/W0 = 1.550`). **No arithmetic or
definition error was found**; the M1 vehicle is retained unchanged as the
**physics-verification baseline** for M2 (distinct from an "orbit-capable vehicle,"
which remains an open later-milestone design question — see §7.4, unchanged).

### 12.2 Atmosphere model implemented (`src/ascent/atmosphere.py`)

A **single-exponential** density model — a deliberate, documented scope decision, and a
discrepancy from the M1 §1.6 note that a layered US Standard Atmosphere 1976 table was
"planned": a single exponential is deterministic, trivial to verify analytically, and
sufficient for M2's goal of verifying the dynamics/drag plumbing. It is explicitly not a
high-fidelity launch atmosphere.

```
rho(h) = RHO0 * exp(-h / SCALE_HEIGHT)   for 0 <= h <= H_MAX
rho(h) = RHO0                            for h < 0        (clamped to sea level)
rho(h) = 0                               for h > H_MAX    (vacuum approximation)
```

`RHO0 = 1.225 kg/m^3` (sea-level reference), `SCALE_HEIGHT = 8500 m`,
`H_MAX = 100,000 m`. Tested for: sea-level reference value, positivity, monotonic
decrease, exact exponential scaling (`rho(8500)/rho(0) = e^-1` etc.), exact-zero
high-altitude limit, and negative-altitude clamping (`tests/test_atmosphere.py`, 7
tests).

### 12.3 Propulsion model implemented (`src/ascent/propulsion.py`)

```
mdot           = -T / (Isp * g0)
a_T            = T / m
m(t)           = m0 + mdot * t                 (mdot negative; unclamped form)
t_burn         = m_propellant / |mdot|
dv_ideal       = Isp * g0 * ln(m0 / mf)
```

plus `clamped_mass_flow_rate`, which forces `mdot = 0` once `mass <= m_min`
(`m_min = m_dry + m_payload`) — mass can never be driven below dry+payload regardless of
commanded thrust. Verified against the exact M1 hand-calc values (`t_burn`, `mdot`,
`dv_ideal`, `T/W0`), linear mass depletion, integrated-mass-loss consistency, and
monotonically increasing thrust acceleration as mass falls (`tests/test_propulsion.py`,
8 tests).

### 12.4 Point-mass ascent ODE implemented (`src/ascent/dynamics.py`)

State `y = [r, theta, v, gamma, m]`, exactly the M1 §2 convention. General form (thrust
direction `chi`, angle of attack `alpha = chi - gamma`, physics separated from any
guidance law via a `control(t, y, params) -> (thrust, chi)` callback):

```
r_dot     = v * sin(gamma)
theta_dot = v * cos(gamma) / r
v_dot     = (T * cos(alpha) - D) / m  -  g(r) * sin(gamma)
gamma_dot = [ T * sin(alpha) / m + (v**2/r - g(r)) * cos(gamma) ] / v      (v >= V_FLOOR)
gamma_dot = 0                                                              (v <  V_FLOOR)
m_dot     = clamped_mass_flow_rate(...)
g(r)      = mu / r**2
```

When `chi == gamma` (zero angle of attack, thrust aligned with velocity — the M1
gravity-turn assumption), `alpha = 0` and this reduces **exactly** to the M1 §3
equations; verified algebraically in
`test_dynamics.py::test_zero_angle_of_attack_reduces_to_m1_equations`.

**`V_FLOOR = 1e-3 m/s`, gamma-freeze at liftoff.** `gamma_dot` divides by `v` and is
singular as `v -> 0`. Below `V_FLOOR`, `gamma_dot` is held at exactly 0 rather than
integrated — a standard, explicitly documented handling of the well-known gravity-turn
liftoff singularity, not a bug.

### 12.5 Earth rotation / atmosphere-relative velocity (`dynamics.py`)

Preserved exactly as the distinction M1 §3 required. The rotating atmosphere is
approximated, within this planar model, as purely tangential (horizontal) with radius-
dependent magnitude `v_atm(r) = omega_earth * r * cos(lat)`, which reduces to the M1
surface value at `r = R_earth` (`tests/test_earth_rotation.py::
test_surface_rotational_speed_matches_m1_hand_calc`, exact match to `408.7 m/s`). This
is an explicit planar approximation: exact for a purely equatorial-tangent plane, and
the natural simplification for the due-east baseline case (`inclination = latitude`)
where the orbital plane's local horizontal locally coincides with Earth's rotation
direction.

Relative (air-relative) speed is computed from the full radial+tangential velocity
vector, not just the tangential component:

```
v_radial     = v * sin(gamma)
v_tangential = v * cos(gamma)
v_rel        = sqrt(v_radial**2 + (v_tangential - v_atm(r))**2)
```

Drag magnitude uses `v_rel`; consistent with the M1 §3 equations (which place drag only
in the `v_dot` equation, not `gamma_dot`), drag's *direction* is taken anti-parallel to
the inertial velocity in these two ODEs — only its magnitude uses the relative-wind
speed. A fully vector-resolved drag treatment is a possible future refinement, not
implemented in M2.

Four dedicated tests (`test_earth_rotation.py`) were written specifically to fail if
inertial speed were substituted for atmosphere-relative speed: a stationary/co-rotating
pad vehicle has exactly zero relative speed; a due-east inertial velocity gets the
correct rotation subtraction (and is shown to differ materially from raw inertial
speed); the zero-Earth-rotation limit reduces to plain inertial speed; and a purely
radial (vertical) inertial velocity still picks up the correct nonzero relative
tangential component from atmospheric rotation. All pass.

**Genuine issue found and corrected during M2 (not an M1 error, an M2 implementation
choice):** an early draft of the diagnostic script (§12.9) initialized the ascent state
with `v0 = 0` (idealized "at rest," ignoring the pad's co-rotation speed). This is
inconsistent with M1 §2's own statement that "at liftoff the vehicle is already
co-rotating... its inertial velocity is not zero." With `v0 = 0`, the relative-speed
formula above produces a spurious ~409 m/s "relative wind" at `t = 0` (since the
vehicle's encoded tangential inertial speed is 0 while the atmosphere's is `v_rot`),
which showed up as an unphysical maximum dynamic pressure exactly at liftoff. This was
corrected by initializing the diagnostic trajectory with the physically correct
`v0 = v_rot`, `gamma0 = 0` (purely horizontal, co-rotating) initial condition, which
resolves the artifact — max dynamic pressure now occurs at `t ~= 54.5 s`,
`h ~= 10.5 km`, matching the expected real-vehicle max-Q altitude range and the M1 §7
hand-estimate order of magnitude closely. See §12.9. This correction only affects the
diagnostic script's initial condition, not any dynamics equation, unit test, or
documented M1 result.

One further, real (not a bug) consequence of this corrected initial condition is worth
recording: because the vehicle's inertial velocity at liftoff is already substantially
horizontal (`v_rot ~= 409 m/s`) rather than zero, a purely radial ("vertical") thrust
command is nearly perpendicular to the velocity vector at that instant, so it initially
curves `gamma` upward relatively gradually rather than producing an immediate steep
climb — mathematically the same effect as a radial burn on a near-circular orbit
changing orbit shape faster than it changes speed. This is a correctly-modeled
consequence of describing ascent in the true inertial frame with Earth rotation
included, not a defect, and is exactly the kind of subtlety M2's verification-first
approach is meant to surface honestly rather than paper over.

### 12.6 Prescribed controls (`src/ascent/controls.py`)

Four control functions, none of them a guidance law:

- `vertical_thrust_control` — Case A, constant thrust held vertical (`chi = 90 deg`).
- `fixed_pitch_control` — Case B, constant thrust at a fixed, prescribed `chi`.
- `coast_control` — Case C, thrust off.
- `vertical_rise_then_pitch_kick_control` — a minimal, explicitly diagnostic-only
  realization of the M1 §1.7 pitch-kick concept (vertical, then one fixed pitch angle
  held constant); not tuned, not a continuously-updated guidance law.

### 12.7 Events and termination (`dynamics.py`, `simulation.py`)

| Event | Direction | Terminal? |
|---|---|---|
| Propellant depletion (`m` crosses `m_min` from above) | `-1` | No (logged for diagnostics) |
| Ground impact (`r` crosses `R_earth` from above) | `-1` | Yes |
| Target-altitude crossing (optional, diagnostic) | `0` | No (default) |

Mass is independently guaranteed never to drop below `m_min` by the RHS's own clamp
(`propulsion.clamped_mass_flow_rate`), verified directly in
`test_verification.py::test_propellant_depletion_event_fires_at_expected_time` (mass
stays within 1e-2 kg of `m_min`, i.e. a ~2e-8 relative floating-point/adaptive-step
residual at the event crossing, not a real propellant debt) and in the same test's
depletion-time check (fires within 0.1% of the analytic `t_burn`). Ground-impact
termination is verified on a descending ballistic case
(`test_ground_impact_event_terminates_descending_trajectory`).

### 12.8 Independent verification results

All required checks (M2 objective §8 A-H) pass:

| Check | Result |
|---|---|
| A. Zero-drag limit | `Cd=0` removes drag exactly (`drag_acceleration(...) == 0.0`); zero-drag RHS gives strictly larger `v_dot` than the drag-on case |
| B. Zero-thrust ballistic limit | Specific orbital energy and angular momentum constant to within `< 1e-8` relative over a 600 s coast (drag disabled to isolate gravity) |
| C. Constant radial-thrust sanity case | Forward-Euler estimate from the analytic RHS at `t=0` matches a `dt=1e-3 s` integration to `rel=1e-3` (residual consistent with 2nd-order Euler truncation) |
| D. Mass-flow analytical check | Integrated mass history matches `m(t) = m0 - |mdot| t` to `rtol=1e-6` over the whole burn |
| E. Tsiolkovsky consistency | With `mu=0, Cd=0` (propulsion-only reduction of the full ascent RHS, not a separately hand-coded formula), integrated `delta_v` matches `Isp*g0*ln(m0/mf)` to `rel=1e-4` |
| F. Earth-rotation check | `atmosphere_corotation_speed(R_earth, lat)` reproduces the M1 hand value `408.7 m/s` exactly |
| G. Integrator convergence | 3 tolerance/step settings (`rtol` 1e-6 -> 1e-9 -> 1e-12); successive differences shrink, tightest two agree to `rtol=1e-5` |
| H. Dimensional/sign checks | Gravity decreases with altitude and matches `mu/r**2` exactly; drag only ever subtracts from `v_dot`; mass strictly decreases only while thrusting (`m_dot < 0` thrusting, `== 0` coasting or depleted); `h = r - R_earth` exact |

45 tests total (`pytest -W error`), all passing.

### 12.9 M2 diagnostic trajectory (`scripts/m2_diagnostic_trajectory.py`)

Unchanged M1 physics-verification baseline vehicle, integrated with
`vertical_rise_then_pitch_kick_control` (vertical for the first 10 s, then a single
fixed, untuned pitch-kick to `chi = 88 deg`, held constant for the rest of the burn — no
gravity-turn guidance law). Initial condition: `r0 = R_earth`, `v0 = v_rot = 408.7 m/s`,
`gamma0 = 0` (co-rotating with the launch site, per §2 and the correction in §12.5).

**Figure:** [`figures/m2_diagnostic_trajectory.png`](figures/m2_diagnostic_trajectory.png)
— explicitly titled *"M2 physics-verification trajectory — prescribed control, NOT
optimized for orbit"*. This is diagnostic/supporting only, not a validated result (per
the repo-wide rule in `README.md`).

**Results:**

| Quantity | Value |
|---|---|
| Burn duration (matches M1 hand calc) | 158.7 s |
| Burnout mass | 90,000.0 kg (== `m_min` exactly, propellant fully depleted) |
| Max altitude | 903.8 km, at t = 594.8 s |
| Max speed (inertial) | 3543.6 m/s, at t = 158.8 s (burnout) |
| Max dynamic pressure | 34.1 kPa, at t = 54.5 s, h = 10.5 km |
| Ground impact within 758.7 s window | No (still descending, h = 802.7 km, gamma = -67.8 deg, at end of window) |
| 400 km circular-orbit conditions achieved | **No** |

**Interpretation:** the vehicle flies well past 400 km altitude (apogee ~904 km) because
the fixed, near-vertical pitch profile used here does not turn the trajectory toward
horizontal the way a real gravity-turn guidance law would — it is not a guided ascent.
It does **not** achieve circular-orbit conditions at 400 km (`v = 7668.6 m/s`,
`gamma ~= 0`) at any point; by burnout it has `v = 3543.6 m/s` at `gamma ~= 76 deg`,
nowhere near the required speed or horizontal flight-path angle, consistent with the
Δv deficit already identified in M1 §7.4 (ideal Δv ~5.05 km/s vs a rough ~8.9-9.5 km/s
LEO budget) compounded by this pitch profile not being an efficient (velocity-aligned)
gravity turn. **This "no" is the expected, honest result for this milestone** — no
assumption was retuned to force a different answer, per the M2 task instructions.

The max-dynamic-pressure result (34.1 kPa at ~10.5 km altitude) is a strong,
independent cross-check on the atmosphere/drag/relative-velocity implementation: it
lands within the same order of magnitude as, and at very close to the same altitude as
assumed by, the M1 §7 hand estimate (~40.5 kPa, assumed at 11-13 km) and real vehicles'
typical max-Q (e.g. ~30 kPa for Falcon 9).

### 12.10 Discrepancies from M1 (summary)

1. Atmosphere: single-exponential model implemented instead of the layered US Standard
   Atmosphere 1976 that M1 §1.6 listed as "planned" — a documented M2 scope decision
   (§12.2), not an error; a higher-fidelity table remains a possible future refinement.
2. Diagnostic-script initial condition corrected mid-M2 from an M1-inconsistent
   `v0 = 0` to the M1-documented `v0 = v_rot, gamma0 = 0` (§12.5) — this is a correction
   toward, not away from, M1's stated physics, and did not change any M1 equation,
   constant, or hand-calculation result.

No other M1 equation, constant, or hand-calculation result was altered.

---

## 13. Milestone 3 — gravity-turn guidance, orbital diagnostics, and one verified trajectory

> **M3 is physics/guidance verification on the unchanged vehicle, not vehicle tuning.**
> No mass, thrust, Isp, staging, payload, or target-altitude change was made anywhere in
> this milestone.

### 13.1 M2 reconfirmation (before any M3 code was written)

`scripts/m2_diagnostic_trajectory.py` was re-run unchanged and reproduced its M2
values exactly: burn duration 158.71 s, burnout mass 90,000.0 kg, max dynamic pressure
34.06 kPa (at ~10.5 km altitude, t ~= 54.5 s), max altitude 903.80 km, orbit not
achieved. No discrepancy found; M3 proceeded on the unchanged M2 baseline.

### 13.2 Gravity-turn guidance law (`src/ascent/controls.gravity_turn_control`)

Deterministic, parameterized, and separate from `dynamics.py` (a `control(t, y, params)`
callback, same interface as every other M2 control). Three phases, triggered by **time**
(not altitude — see the note below):

1. **Vertical rise** (`t < kick_start_time`): `chi = 90 deg` (thrust radial/vertical).
2. **Pitch-kick** (`kick_start_time <= t < kick_start_time + kick_duration`, only if
   `kick_duration > 0`): thrust held at a **fixed** offset from vertical,
   `chi = 90 deg - kick_angle`, for the whole kick duration. This is a thrust-*direction*
   change only; it does not modify the velocity state (§13.10 confirms continuity). If
   `kick_duration == 0`, this phase has zero width and `kick_angle` is unused (check A).
3. **Zero-angle-of-attack gravity turn** (`t >= kick_start_time + kick_duration`):
   `chi = gamma(t)`, re-evaluated every RHS call — thrust tracks the *current* velocity
   direction, so gamma is bent by gravity/dynamics, not prescribed directly (§13.10).

Thrust cutoff at propellant depletion is enforced by `dynamics.ascent_rhs`'s existing
mass-flow clamp, not duplicated in the guidance function.

**Why time, not altitude, triggers the kick:** an altitude-triggered switch would need
its own root-finding event mid-integration to restart the control law; a time trigger is
simpler, fully deterministic, and sufficient for a bounded, interpretable M3 sweep. Each
sweep row still reports the resulting altitude at which the kick occurs, so the
altitude information is not lost, only derived rather than used as the trigger.

**A real gravity-turn sensitivity, found and documented (not a bug):** because the
correct M2 initial condition is `v0 = v_rot ~= 409 m/s` (horizontal, `gamma0 = 0`, not
zero — DESIGN.md §12.5), and because the zero-AoA `gamma_dot` equation
(`[T sin(alpha)/m + (v^2/r - g) cos(gamma)] / v`) has NO thrust term once `alpha = 0`,
`gamma` decreases continuously after the kick ends whenever `v` is below local circular
speed (true for essentially all of ascent) — the well-known mechanism that makes a
gravity turn work. If the kick happens too early (before enough speed has built up),
`gamma` decays through zero and the vehicle noses over into the ground **before
burnout**. Exploration confirmed this cleanly: with this vehicle, kick-start times of 10,
20, and 30 s crash into the ground well before the 158.7 s burn completes, **regardless
of kick angle (5-30 deg tested) or kick duration (5 or 10 s tested)** — only
`kick_start_time >= ~32 s` (for a 15 deg/10 s kick; the exact boundary depends on angle)
survives to burnout. This is a genuine, physically-grounded finding about this vehicle's
thrust-to-weight and rotation-driven initial condition, not an implementation defect,
and it directly shaped the sweep range chosen below.

### 13.3 Orbital diagnostics (`src/ascent/orbital.py`)

Independently derived from `(r, v, gamma)`:

```
v_radial, v_tangential = v*sin(gamma), v*cos(gamma)
specific_energy        = v^2/2 - mu/r
h_ang                  = r * v_tangential                      (specific angular momentum)
p                      = h_ang^2 / mu                          (semi-latus rectum)
eccentricity           = sqrt(1 + 2*specific_energy*h_ang^2/mu^2)
perigee_radius         = p / (1 + eccentricity)                (defined for any conic)
semi_major_axis        = -mu / (2*specific_energy)             (bound orbits only, e<1)
apogee_radius          = p / (1 - eccentricity)                (bound orbits only, e<1)
circular_speed_here    = sqrt(mu/r)
```

Unbound (`eccentricity >= 1`, parabolic/hyperbolic) states correctly report
`semi_major_axis = None` and `apogee_altitude = None` — no bound-orbit size is fabricated
for them. A `perigee_radius < R_earth` sets `intersects_earth = True` (the current
osculating conic, if propagated unpowered/undragged, would hit the ground before
completing an orbit — i.e. this is a suborbital trajectory even if its instantaneous
altitude is high). Tested in `tests/test_orbital.py`: circular case (energy = `-mu/2r`
exactly, perigee = apogee = altitude), elliptical case cross-checked against vis-viva
(`v^2 = mu(2/r - 1/a)`), hyperbolic case (no apogee, e >= 1), a suborbital case
correctly flagged `intersects_earth = True`, and an explicit "high apogee alone is not
orbit" case.

### 13.4 Orbit-achievement criterion (explicit, before the sweep)

```
is_circular_orbit_achieved(state, target_altitude, altitude_tolerance=10 km) :=
    state.is_bound
    AND NOT state.intersects_earth
    AND |perigee_altitude - target_altitude| <= altitude_tolerance
    AND |apogee_altitude  - target_altitude| <= altitude_tolerance
```

This requires a near-circular orbit **at** the target altitude — both perigee and
apogee within 10 km of 400 km — and explicitly does **not** accept "max altitude
> 400 km" as sufficient (a lofted suborbital trajectory can have an apogee far above
400 km while its perigee is deep inside the Earth; `tests/test_orbital.py::
test_max_altitude_alone_does_not_imply_orbit_achieved` checks exactly this).

### 13.5 Guidance parameter sweep (`scripts/m3_gravity_turn_sweep.py`)

Bounded 2D grid, **5 x 5 = 25 cases**, unchanged vehicle:

| Parameter | Values swept |
|---|---|
| Pitch-kick start time | 20, 30, 40, 50, 60 s |
| Pitch-kick angle | 5, 10, 15, 20, 25 deg |
| Pitch-kick duration | fixed at 10 s for every case (not swept) |

**Objective metric:** maximize burnout specific orbital energy
(`epsilon = v^2/2 - mu/r` at the propellant-depletion state) among cases that survive
powered flight to burnout. This is transparent and single-valued; it does not attempt to
force orbital insertion.

**Results (full CSV: [`scripts/m3_sweep_results.csv`](scripts/m3_sweep_results.csv)):**
17/25 cases reached burnout (8 crashed into the ground before burnout — all at
`kick_start_time = 20 s`, plus 3 of 5 at `kick_start_time = 30 s`, consistent with
§13.2's finding). **0/25 cases achieved the 400 km circular orbit criterion (§13.4).**
This is reported honestly, not hidden — it is the expected result given M1 §7.4's
already-documented Δv shortfall.

**Selected case (max burnout specific energy among survivors):**
`kick_start_time = 40 s`, `kick_angle = 25 deg`, `kick_duration = 10 s`
(`specific_energy = -5.217e7 J/kg`). Energy at `kick_start_time = 40 s` is fairly flat
across kick angle (`-5.218e7` to `-5.224e7 J/kg`, a ~0.2% spread) — the timing of the
kick matters far more than its magnitude for this objective, another genuine,
interpretable finding from the sweep, not a cherry-picked result.

### 13.6 Selected trajectory (`scripts/m3_trajectory.py`)

Figures: [`figures/m3_gravity_turn_trajectory.png`](figures/m3_gravity_turn_trajectory.png)
(altitude, speed, flight-path angle, dynamic pressure, mass vs. time, with pitch-kick,
burnout, ascent max-Q, and apogee marked) and
[`figures/m3_velocity_components.png`](figures/m3_velocity_components.png) (radial vs.
tangential velocity; altitude vs. downrange angle). Both titled
*"M3 gravity-turn trajectory — unchanged M1/M2 verification vehicle"* and explicitly
labeled diagnostic/supporting, not a validated flight result.

**Burnout state:** t = 158.71 s, altitude = 42.19 km, speed = 4451.3 m/s,
gamma = 1.43 deg, specific orbital energy = -5.218e7 J/kg, eccentricity = 0.681.
**Perigee altitude = -5160 km (`intersects_earth = True`), apogee altitude = 43.1 km.**
**400 km circular-orbit criterion: NOT achieved.**

The vehicle follows this trajectory to apogee at t ~= 175.7 s (h ~= 43.1 km, matching
the analytic burnout-state apogee to within 0.03 km — see §13.9 check G/apogee
cross-check) and then, because its perigee is deep inside the Earth, re-enters and
impacts the ground at t ~= 297.2 s.

**Max dynamic pressure during ascent (t <= burnout) = 70.1 kPa**, occurring essentially
at burnout itself (this trajectory's altitude stays low enough, 42 km, that dynamic
pressure is still rising when propellant runs out, rather than peaking mid-flight the
way the M2 diagnostic's steeper trajectory did — a genuine feature of this particular
low-loft guidance solution, not an error). **Max dynamic pressure overall = 1803 kPa**,
reached during the high-speed terminal re-entry dive (t ~= 291 s) — a real but distinct
phenomenon from ascent aerodynamic loads, and is called out separately on the figure so
it is not mistaken for ascent max-Q.

### 13.7 Comparison against the M2 prescribed-control diagnostic

| Quantity | M2 diagnostic (vertical + one fixed 88 deg pitch) | M3 selected gravity turn |
|---|---|---|
| Burnout time | 158.7 s | 158.7 s (unchanged vehicle/burn) |
| Burnout altitude | (not applicable; M2 tracked full flight) | 42.2 km |
| Burnout speed | 3543.6 m/s | 4451.3 m/s |
| Burnout gamma | ~76 deg (steep) | 1.4 deg (nearly horizontal) |
| Max altitude reached | 903.8 km | 43.1 km (apogee) |
| Max-Q (ascent) | 34.1 kPa @ ~10.5 km | 70.1 kPa @ ~burnout (42 km) |
| Orbit achieved | No | No |

The M3 guidance law achieves a much higher burnout speed and a much shallower burnout
flight-path angle than M2's fixed near-vertical pitch — i.e. it IS a more efficient use
of the same propellant toward orbital insertion (higher specific energy) — but trades
away altitude entirely (42 km vs. 904 km apogee) to do it, and still falls well short of
orbital energy. Both are honest, unglossed results for this vehicle.

### 13.8 Event handling refinement

| Event | Direction | Terminal? | M3 status |
|---|---|---|---|
| Propellant depletion | -1 | configurable (`terminal_depletion`, default False, per M2) | Verified: fires within 0.1% of analytic `t_burn`; mass never drops more than ~0.1 kg below `m_min` (~2e-7 relative) across a full gravity-turn run |
| Ground impact | -1 | Yes (always) | Verified directionally correct (`r` strictly decreasing into the event) on both a simple ballistic case and the gravity-turn sweep's crash cases |
| **Apogee (new)** | -1 (`r_dot` crossing + to -) | No (default; can be set terminal) | Verified: fires during coast, and the apogee time matches the sampled-trajectory altitude maximum to within 2 s |
| Target-altitude crossing | 0 | No (default) | Explicitly verified NOT to imply orbit achievement (separate API, `orbital.is_circular_orbit_achieved`) |
| Atmosphere-exit (H_MAX, diagnostic) | 0 | No (default) | Available (`make_atmosphere_exit_event`); a thin wrapper over the altitude-crossing event at the atmosphere model's cutoff |

"Liftoff / positive vertical motion" was considered and deliberately **not** added as a
separate event: this model has no hold-down phase — thrust is on and the state already
has nonzero (co-rotating) velocity from `t = 0` — so there is no ambiguous liftoff
instant for an event to detect.

### 13.9 Independent verification results (checks A-H, DESIGN.md M3 objective)

| Check | Result |
|---|---|
| A. Zero pitch-kick limit | `kick_duration=0`: kick phase has zero width, `chi=gamma` immediately after `kick_start_time` regardless of `kick_angle`; `kick_angle=0`: kicked angle equals vertical exactly. Both verified directly on the control function. |
| B. Small-kick continuity | 5 angles (1-3 deg) at a fixed kick timing: no burnout-energy jump more than 5x any neighboring step — smooth response, away from the crash/survive boundary |
| C. Zero-AoA condition | After the kick phase, `chi == gamma` exactly (`abs=1e-12`) at an arbitrary sampled state |
| D. Coast conservation | Burnout state propagated with `Cd=0`: specific energy and angular momentum constant to `< 1e-8` relative over a 600 s coast |
| E. Event consistency | Depletion event time matches analytic `t_burn` to `rel=1e-3`; burnout mass matches `m_dry+m_payload` to `abs=1e-2` kg |
| F. State reconstruction | `(v, gamma) -> (v_radial, v_tangential) -> (v, gamma)` round-trips to `rel=1e-12` |
| G. Orbital-element consistency | `-mu/(2a)` reconstructed from the burnout state's own semi-major axis matches its `specific_energy` to `rel=1e-6`; independently, `sqrt(mu*a*(1-e^2))` reconstructs `specific_angular_momentum` to `rel=1e-6` in `test_orbital.py` |
| H. Guidance repeatability | Identical `(t, r, v, gamma, m)` arrays (exact equality) across two runs of the same parameters |

Plus a dedicated velocity-continuity sanity check (§13.10) and the full M2 test suite
(37 tests, unchanged) still passing. **60 tests total, all passing under
`pytest -W error`.**

### 13.10 Physical sanity checks (M2 objective §10)

- **No discontinuous velocity change at the pitch-kick.** The kick changes only the
  commanded thrust *direction*; it never modifies the state vector `y`. Verified with a
  tight (±0.02 s) window around each phase boundary: `|Δv| < 2 m/s` (a generous bound —
  a genuine discontinuity would show as O(10-100) m/s given this vehicle's speeds).
  (An earlier draft of this test used a single non-uniform `t_eval` array whose middle
  gap spanned the entire kick duration — a ~4.9 s window, not ~0.1 s — and flagged a
  smooth ~11 m/s change as if it were a discontinuity; that was a test-construction bug,
  not a dynamics bug, and was corrected to bracket each boundary tightly and
  independently.)
- **Drag still uses atmosphere-relative velocity**, unchanged from M2 (`dynamics.
  relative_speed`, not touched in M3).
- **Atmosphere still rotates consistently with M2** (`atmosphere_corotation_speed`,
  unchanged).
- **Thrust stops exactly at propellant depletion** — enforced by the unchanged M2 mass-
  flow clamp; §13.9 check E confirms the cutoff mass and timing.
- **Coast phase consumes no mass** — `coast_control` commands `T=0`, and
  `clamped_mass_flow_rate` returns exactly 0 whenever thrust is 0 (unchanged M2 logic,
  re-verified in `test_events_m3.py`).
- **The flight path is bent by gravity, not prescribed.** `gamma(t)` is never set
  directly by the guidance law; only the thrust-direction angle `chi` is prescribed, and
  `gamma` is always obtained by integrating the dynamics ODE (§13.2).

### 13.11 Integrator convergence (selected trajectory, powered phase)

Three settings, `(max_step, rtol, atol)`:

| Setting | Burnout altitude [km] | Burnout speed [m/s] | Burnout gamma [deg] | Max-Q [kPa] | Apogee altitude [km] |
|---|---|---|---|---|---|
| `(2.0, 1e-6, 1e-6)` | 42.127 | 4452.115 | 1.4466 | 70.413 | 43.086 |
| `(0.5, 1e-9, 1e-9)` | 42.158 | 4452.072 | 1.4552 | 70.153 | 43.129 |
| `(0.1, 1e-12, 1e-12)` | 42.158 | 4452.072 | 1.4552 | 70.153 | 43.129 |

Burnout time is identical (event-detected, `t = 158.71288815789472 s`) at all three
settings, as expected. The tightest two settings agree to 5-6 significant figures on
every reported quantity (altitude within 0.5 m, speed within 0.001 m/s, gamma within
1e-5 deg, max-Q within 0.0005 kPa, apogee within 0.5 m) — convergence is not claimed
merely because `solve_ivp` reported success.

### 13.12 Was 400 km circular orbit achieved?

**No**, for any of the 25 swept cases, including the selected (best-by-objective) one.
This is the expected, honest outcome given M1 §7.4's already-documented ideal-Δv
shortfall (~5.05 km/s available vs. a rough ~8.9-9.5 km/s LEO budget), now confirmed
against an actual integrated, verified gravity-turn guidance law rather than a rough
scalar estimate. No vehicle parameter was changed to try to close this gap — that
remains explicitly out of scope until a later milestone (M1 §7.4, restated in the M3
task instructions).

### 13.13 M3 limitations

- The pitch-kick trigger is time-based, not altitude-based (§13.2); an altitude-based
  trigger would need its own guidance-restart event and was judged unnecessary
  complexity for a bounded M3 sweep.
- The guidance law's "kicked" angle is held fixed (not linearly ramped) during the kick
  phase — a simpler, equally defensible choice, documented in `controls.py`.
- The sweep is a 2D grid (kick timing x kick angle) with kick duration fixed at 10 s;
  duration itself was not swept, per the instruction to keep the search small.
- The orbit-achievement altitude tolerance (10 km) is a documented, somewhat arbitrary
  diagnostic choice, not derived from a mission requirement.
- As in M2, drag's magnitude uses atmosphere-relative speed but its direction is taken
  anti-parallel to the inertial velocity (not to the relative-velocity vector) in the
  `v_dot`/`gamma_dot` equations — unchanged from M2, not revisited in M3.
- No structural/thermal load limit is checked against the max-Q or terminal re-entry
  dynamic-pressure values reported in §13.6 — they are reported as diagnostics only.

---

## 14. Milestone 4 — orbit-capable study vehicle and payload-to-orbit solve

> **M4 asks a narrower question than "reach orbit": at 400 km / 28.5°, what
> vehicle/payload combinations actually achieve the defined insertion condition?**
> The M1–M3 verification vehicle is retained, unmodified, as a documented FAILURE case
> under the M4 criterion (§14.10 check A) — it is never silently replaced. A new,
> explicitly separate **M4 orbit-capable study vehicle** is defined instead (§14.3).

### 14.1 M3 reconfirmation (before any M4 code was written)

The M3-selected case was independently reproduced exactly: kick_start=40 s,
kick_angle=25°, kick_duration=10 s → burnout t=158.71 s, altitude=42.19 km,
speed=4451.3 m/s, gamma=1.43°, specific energy=-5.218e7 J/kg, e=0.6811, perigee
intersects Earth, 400 km circular orbit = **FAIL**. This remains the authoritative M3
result and is not reinterpreted here.

### 14.2 M4 orbit-insertion criterion (defined before any vehicle/guidance search)

DESIGN.md M3 §13.4 already defined a strict "raw ascent alone must already be
near-circular" criterion. M4 instead evaluates a standard **ascent-to-apogee +
idealized apogee-circularization** mission profile (`orbital.evaluate_orbit_insertion`,
`orbital.circularization_delta_v_at_apogee`): propagate the powered ascent, identify a
cutoff/insertion instant, and compute the impulsive tangential delta-v needed to
circularize at that instant's (post-atmosphere) apogee. **The circularization burn is
explicitly an idealized, unconstrained-thrust insertion-stage impulse — it is never
subtracted from the M4 vehicle's own propellant/mass budget** (finite-thrust
upper-stage optimization is out of M4 scope, per the task's explicit scope boundary);
it is reported as a diagnostic quantity.

**PASS requires ALL of:**
1. the post-cutoff osculating orbit is bound (`is_bound`),
2. it does not itself intersect the Earth as a standalone ellipse (`reaches_bound_orbit`
   — a sanity floor, not the headline criterion),
3. the **natural apogee altitude** (before circularizing) is within **15 km** of
   400 km — an idealized circularization burn cannot move WHERE the apogee is, only
   the orbit's shape at that altitude, and
4. the idealized circularization delta-v needed there is **≤ 1500 m/s** — a documented,
   explicit engineering judgment call (not a silently smuggled-in constraint) bounding
   the insertion burn to a plausible kick/insertion-stage impulse, not a second full
   ascent burn.

This explicitly distinguishes **"orbit-capable / reaches a bound LEO"**
(`reaches_bound_orbit` alone) from **"meets the 400 km insertion criterion"**
(`meets_insertion_criterion`, all four conditions) — `tests/test_m4_insertion.py`
includes a dedicated case with high apogee and low (in-cap) circularization delta-v
that nonetheless fails because the natural apogee is far from the target, so "reaches
orbit" is never conflated with "meets this mission's target."

**A genuine discrepancy found and corrected during M4 (not an M1–M3 error):** the
guidance search initially evaluated candidate cutoffs using the **idealized**
(instantaneous, drag-free) osculating apogee at the cutoff sample. Because this
vehicle's best insertion window occurs at a fairly low cutoff altitude (order
70–110 km) — still within the (small but nonzero) modeled atmosphere below
`atmosphere.H_MAX` (100 km) — a real drag-included coast from there up through 100 km
loses measurable additional energy that the idealized instantaneous estimate misses
entirely. On the reference-payload case this was checked directly: the idealized
cutoff-instant apogee was 386 km, but the actual drag-included coast up to 100 km
altitude settled to 293 km — large enough to invalidate the naive criterion, not a
rounding effect. **Fix:** `insertion_search.find_drag_consistent_cutoff` evaluates every
candidate cutoff by actually coasting (thrust off, drag on) from that instant up to
100 km altitude before computing orbital elements, and this is the standard used
throughout M4 (guidance search, payload sweep, final trajectory) — the idealized
`orbital.scan_best_insertion_cutoff` from M3 is retained only as a cheap first-pass
bracketing step, never as the final answer.

### 14.3 M4 orbit-capable study vehicle (`constants.m4_vehicle`, `M4_DESIGN`)

Explicitly separate from `BASELINE_VEHICLE` (M1–M3), which is **never modified**.
Smallest reasonable change set, changing only what M1 §7.4 already identified as the
blocking constraint (propellant mass fraction and Isp); aerodynamics, launch site,
target altitude, and the dynamics/atmosphere framework are all unchanged:

| Quantity | BASELINE_VEHICLE (M1–M3) | M4 study vehicle | Changed? |
|---|---|---|---|
| Dry mass | 80,000 kg | 50,000 kg | **yes** |
| Propellant mass | 410,000 kg | 440,000 kg | **yes** |
| Isp | 300 s | 450 s | **yes** |
| Thrust | 7.6 MN | 7.6 MN | no |
| Reference area | 10.75 m² | 10.75 m² | no |
| Drag coefficient | 0.3 | 0.3 | no |
| Payload (reference) | 10,000 kg | 10,000 kg (swept in M4) | reference unchanged |

**Analytical justification (before any propagation), reference payload = 10,000 kg:**

| Quantity | Value |
|---|---|
| m0 | 500,000 kg |
| mf (dry+payload) | 60,000 kg |
| Mass ratio | 8.333 |
| Exhaust velocity (ve = Isp·g0) | 4413.0 m/s |
| **Ideal Δv (Tsiolkovsky)** | **9356.7 m/s** |
| Thrust-to-weight at liftoff | 1.550 (unchanged from M1 — same thrust, same m0) |
| Mass flow rate | 1722.19 kg/s |
| Full-depletion burn duration | 255.49 s |

An Isp of 450 s (vs. M1's 300 s) is representative of a higher-performance (LH2/LOX-
class) propulsion assumption — it is the dominant, explicitly justified change, because
at Isp=300 s no physically plausible structural mass fraction closes the LEO delta-v
budget (M1 §7.4: a mass ratio of ~20+ would be required). The propellant fraction here
implies a ~10.2% dry-mass-fraction-of-stack (50,000 / (50,000+440,000)), aggressive but
comparable to serious SSTO study-vehicle proposals (e.g. VentureStar/X-33-class
targets) — **not** an existing operational vehicle, hence "study vehicle" throughout.
No staging is used (single-stage, per the task's preference to isolate payload
capability cleanly).

### 14.4 Guidance retuning (`scripts/m4_guidance_search.py`)

Same three-phase `controls.gravity_turn_control` law as M3 (vertical rise → fixed-angle
pitch-kick → zero-AoA gravity turn), retuned for the new vehicle via a structured
coarse-to-fine search:

- **Coarse**: kick_start ∈ {35,45,55} s × kick_angle ∈ {15,20,25,30}° × kick_duration ∈
  {15,20,25} s = 36 combinations (a wider 180-combination exploratory sweep was run
  during development and located the same optimum; the committed script uses this
  smaller, still-representative grid so it re-runs in a practical amount of time).
- **Fine**: a local refinement grid around the best coarse point.
- For every combination, the engine-cutoff time is found via
  `find_drag_consistent_cutoff` (§14.2), NOT full propellant depletion — DESIGN.md
  §14.5 explains why an early, commanded cutoff is used.
- **Objective**: minimize the idealized circularization delta-v (§14.2) among
  combinations that meet the full insertion criterion.

Running the committed search: only 1/36 coarse combinations meet the insertion
criterion at all — kick_start=45 s, kick_angle=20°, kick_duration=20 s (circ_dv=96.5
m/s, apogee=400.03 km) — underscoring the same guidance sensitivity already
documented in M3 §13.2 and rediscovered here (§14.2). The fine local refinement around
it finds 13 passing neighbors, including an essentially tied alternative (kick_start=47
s, kick_angle=17°, kick_duration=16 s; circ_dv=96.5 m/s, apogee=400.06 km) — statistically
indistinguishable from the coarse optimum.

**Selected guidance**: kick_start = 45 s, kick_angle = 20°, kick_duration = 20 s (the
coarse-stage optimum, confirmed not meaningfully improved upon by the fine stage).
At the reference payload this gives cutoff_time ≈ 250.87 s, natural (drag-consistent,
post-100 km) apogee ≈ 400.03 km, and circularization Δv ≈ 96.5 m/s.

### 14.5 Commanded engine cutoff (not full propellant depletion)

Unlike M1–M3, the M4 insertion is achieved with a **commanded/guided engine cutoff**
before full propellant depletion (`controls.with_cutoff`, and directly via
`find_drag_consistent_cutoff`'s identified cutoff time) — this is standard real-vehicle
practice ("identify a suitable cutoff/insertion point" per the task instructions), and
was found to be necessary here: the vehicle's very high thrust-to-weight makes the
zero-AoA gravity-turn phase extremely sensitive to timing (a well-known consequence of
the `gamma_dot` equation losing its thrust term once alpha=0 — DESIGN.md M3 §13.2's
finding, which recurs here). At the reference payload, cutoff occurs at t≈250.87 s
against a full-depletion time of 255.49 s, leaving ≈8,800 kg of unburned propellant
onboard as ordinary dead mass (not jettisoned, not idealized away — DESIGN.md §14.6).

An earlier design iteration attempted to instead **resize propellant** so that full
depletion coincided with the insertion window, removing the need for a cutoff layer.
This was tried and rejected: reducing propellant mass changes the vehicle's ENTIRE
mass-vs-time history (a lighter vehicle accelerates faster from t=0), so the state at
the same nominal burn time is a materially different (much more energetic) trajectory
than the original heavier vehicle's state at that same instant — it does not simply
"stop the original trajectory earlier." Only an explicit, independent commanded-cutoff
control (leaving the mass history otherwise unchanged) reproduces the state actually
found by the guidance search. This is documented rather than silently discovered and
discarded.

### 14.6 Mass bookkeeping (DESIGN.md M4 S6)

`m0 = m_dry + m_propellant + payload`; `m_propellant` and `m_dry` are FIXED
(`M4VehicleDesign`) across the payload sweep — **only `payload` (and therefore `m0`)
varies**. `m_min = m_dry + payload` (propellant fully depleted with payload still
attached). Because M4 uses a commanded cutoff before full depletion, the mass AT
CUTOFF is `m0 - mdot·t_cutoff > m_min` (some propellant remains unburned) — this is
distinct from, and does not weaken, the underlying floor: `tests/test_dynamics.py`
(M2, unchanged) and `tests/test_m4_vehicle.py` both confirm mass can never be driven
below `m_min` regardless of commanded thrust, and if a case is propagated to full
depletion (no cutoff), burnout mass equals `m_min` exactly (M1–M3 behavior, unchanged).
`tests/test_m4_vehicle.py` specifically targets the bugs named in the task instructions:
payload treated as propellant, m0 held fixed while payload varies, mass allowed below
`m_min`, and double-counted payload — all four are explicitly tested and pass.

### 14.7 Payload-to-orbit solve (`scripts/m4_payload_sweep.py`)

Method: (1) coarse upward sweep from the reference payload (10,000 kg, confirmed PASS)
in 1000 kg steps until a FAIL is found, bracketing the boundary; (2) bisection of that
bracket directly on the PASS/FAIL boolean, to a 1 kg tolerance — plain bisection was
used rather than forcing a smooth-function root-finder (e.g. Brent) onto what turned
out to be a near-step transition (§14.9); (3) the required case set (below) plus extra
low-payload points investigating a genuine non-monotonicity (§14.9 check B). The
bracket search fails loudly (raises) if the reference payload itself fails, or if no
failing payload is found by 100,000 kg — it does not silently guess a bracket.

**Headline figure**: [`figures/m4_payload_capability.png`](figures/m4_payload_capability.png)
(`scripts/m4_payload_figure.py`) — payload vs. apogee altitude (PASS/FAIL, target band,
max-payload line) and payload vs. idealized circularization Δv, from a denser
(non-required-case) sweep across the transition region.

**Result: maximum payload = 10,333 kg** (10,334 kg confirmed failing; the transition is
sharp — see §14.9).

### 14.8 Required payload cases

Full CSV: [`scripts/m4_payload_sweep_results.csv`](scripts/m4_payload_sweep_results.csv).

| Payload [kg] | m0 [kg] | Burnout alt [km] | Burnout v [m/s] | Burnout γ [deg] | Max-Q [kPa] | Apogee [km] | Circ. Δv [m/s] | Orbit achieved |
|---|---|---|---|---|---|---|---|---|
| 0 | 490,000 | 108.9 | 9194.8 | 3.63 | 29.4 | — | — | **FAIL** (closest apogee miss 0.015 km — extremely close, but fails other conditions, §14.9) |
| 5,000 | 495,000 | 90.3 | 8861.1 | 1.88 | 28.8 | — | — | **FAIL** (closest miss 0.033 km) |
| 7,500 | 497,500 | 80.4 | 7946.1 | 0.96 | 28.5 | 400.01 | 104.8 | **PASS** |
| **10,000 (reference)** | 500,000 | 72.4 | 8010.5 | 0.114 | 28.2 | 400.03 | 96.5 | **PASS** |
| 9,333 (near-boundary) | 499,333 | 74.6 | 7983.1 | 0.339 | 28.3 | 399.95 | 96.8 | **PASS** |
| **10,333 (maximum)** | 500,333 | 71.3 | 8031.1 | 0.0017 | 28.2 | 399.89 | 96.8 | **PASS** |
| 10,334 (just failing) | 500,334 | 71.3 | 8538.4 | 0.032 | 28.2 | closest miss ≈393 km | — | **FAIL** |
| 10,833 (clearly failing) | 500,833 | 69.6 | 8509.1 | -0.139 | 28.1 | closest miss ≈388 km | — | **FAIL** |

(Exact floating-point values are in `scripts/m4_payload_sweep_results.csv`. The two
"failing" rows' burnout state is reported at full propellant depletion, since no
passing engine-cutoff instant was found for those payloads.)

### 14.9 Numerical verification (checks A–J)

| Check | Result |
|---|---|
| A. Original vehicle remains a failure | `BASELINE_VEHICLE` re-run through `find_drag_consistent_cutoff` with a representative M3-style guidance attempt: **no passing cutoff found** — consistent with M1 §7.4 and M3 §13.12 |
| B. Payload monotonicity | **Investigated, real, explained non-monotonicity found**: 0 kg and 5,000 kg payload FAIL — remarkably, within 15-33 METERS of the target apogee — while 7,500-10,333 kg PASS. Root cause (not a bug): the fixed guidance profile was tuned for the 10,000 kg reference; a LIGHTER vehicle has more excess energy and, at the same kick timing, "overshoots" past the narrow viable insertion window before its perigee clears the ground (traced directly: at payload=0, even climbing states (γ>0) extremely near the target apogee still have `intersects_earth=True`). This is a property of using one FIXED guidance profile across the whole payload range, not a defect in the dynamics/orbital-element code, and it underscores just how narrow this vehicle's viable insertion window is. It does not affect the maximum-payload result, which sits in the monotonically-behaved region (7,500 kg → 10,333 kg PASS → 10,334 kg FAIL) |
| C. Payload-boundary bracketing | 10,333 kg confirmed PASS and 10,334 kg confirmed FAIL by direct construction (`assert`s in `m4_payload_sweep.py::main` and `tests/test_m4_payload_boundary.py`) |
| D. Mass bookkeeping | `tests/test_m4_vehicle.py`: m0 varies correctly with payload (dry/propellant fixed), `m_min = m_dry+payload`, payload never double-counted or treated as propellant |
| E. Tsiolkovsky consistency | `tests/test_m4_loss_budget.py::test_zero_drag_zero_alpha_case_matches_pure_tsiolkovsky` — achieved Δv matches `Isp·g0·ln(m0/mf)` to `rel=1e-6` for the M4 vehicle at the reference payload, mu=0/Cd=0/alpha=0 reduction |
| F. Energy/orbital-element consistency | `orbital.py`'s existing self-consistency checks (M3, unchanged) plus the exact identity check in the loss budget (§14.11) |
| G. Guidance repeatability | `tests/test_m4_payload_boundary.py::test_drag_consistent_cutoff_search_is_repeatable` — identical cutoff time and circularization Δv across two runs |
| H. Integrator convergence | §14.12 |
| I. Atmosphere-relative drag | `tests/test_m4_payload_boundary.py::test_m4_dynamics_still_uses_relative_speed_for_drag` — regression check that M4 code paths still route drag through `dynamics.relative_speed` (M2, unchanged) |
| J. M1–M3 regression | All 60 pre-existing tests still pass unchanged (verified before every commit) |

### 14.10 Δv / loss accounting (reference payload, to engine cutoff)

Exact derivation from the M1 §3 / M2 `dynamics.py` equation of motion
(`v_dot = (T·cos(alpha) - D)/m - g(r)·sin(gamma)`), implemented in
`src/ascent/loss_budget.py`:

```
ideal_dv_to_cutoff = ve * ln(m0 / m(cutoff))              (mass ACTUALLY consumed by cutoff)
achieved_dv        = v(cutoff) - v(0)
gravity_loss       = INTEGRAL[ g(r) * sin(gamma) ] dt      (trapezoidal, over the sampled run)
drag_loss          = INTEGRAL[ D/m ] dt                    (trapezoidal)
steering_loss      = ideal_dv_to_cutoff - achieved_dv - gravity_loss - drag_loss   (exact
                                                            algebraic remainder, NOT a
                                                            fitted/unexplained residual --
                                                            see loss_budget.py docstring)
```

| Term | Value |
|---|---|
| Ideal Δv to cutoff (Tsiolkovsky, mass consumed by t=250.9 s) | 8807.4 m/s |
| Achieved inertial Δv | 7601.8 m/s |
| Gravity loss | 617.2 m/s |
| Drag loss | 28.0 m/s |
| **Steering loss** (exact remainder — nonzero angle-of-attack during the pitch-kick) | 560.4 m/s |
| Identity residual (`ideal - (achieved+gravity+drag+steering)`) | 0.0 m/s (exact, to numerical-integration precision) |

Figure: [`figures/m4_delta_v_budget.png`](figures/m4_delta_v_budget.png). Drag loss is
small (this vehicle spends most of its early, high-density-atmosphere time at fairly
low speed); gravity and steering losses are comparable in size, both driven by the
same early, still-mostly-vertical phase of flight.

### 14.11 Integrator convergence

Three solver settings (`max_step`, `rtol`, `atol`) applied to the maximum-passing
payload (10,333 kg) powered-ascent-to-cutoff propagation (cutoff t=251.20 s):

| Setting | Final r [m] | Final v [m/s] | Final γ [deg] | Final m [kg] |
|---|---|---|---|---|
| (2.0, 1e-6, 1e-6) | 6,449,192.81 | 8031.860 | -0.02586 | 67,711.3726 |
| (0.3, 1e-9, 1e-9) | 6,449,474.82 | 8031.142 | 0.00168 | 67,711.3726 |
| (0.05, 1e-12, 1e-12) | 6,449,475.14 | 8031.141 | 0.00171 | 67,711.3726 |

Tightest two settings agree to within 0.3 m in r, 8e-4 m/s in v, ~3e-5 deg in gamma,
and effectively exact in mass — convergence is not claimed merely from solver success.

(`tests/test_m4_convergence.py` asserts the tightest two settings agree to within 1 m
in r, 1 cm/s in v, ~6e-5 deg in gamma, and 1 g in mass — all satisfied.)

### 14.12 Comparison: M3 vs. M4

| | M3 | M4 |
|---|---|---|
| Vehicle | Unchanged M1 verification vehicle | New, explicitly separate "M4 orbit-capable study vehicle" |
| Isp | 300 s | 450 s |
| Dry / propellant mass | 80,000 / 410,000 kg | 50,000 / 440,000 kg |
| Guidance | Retuned gravity-turn, full-depletion burnout | Retuned gravity-turn, **commanded early cutoff** |
| Insertion method | Raw ascent state must already be near-circular | Ascent-to-apogee + idealized circularization impulse |
| 400 km circular orbit | **FAIL** (0/25 sweep cases) | **PASS**, up to 10,333 kg payload |

M4 does not "fix" M3's vehicle — it defines and analyzes a genuinely different,
explicitly labeled vehicle, and answers a genuinely different (narrower, standard
ascent+circularization) mission question at the same target orbit and inclination.

### 14.13 M4 limitations

- The payload sweep uses ONE fixed guidance profile (kick_start/angle/duration) across
  all payloads; §14.9 check B's investigated non-monotonicity at very low payload is a
  direct consequence of this simplification, not re-tuned away.
- The commanded engine cutoff is found by a numerical search (coarse bracket + bounded
  scalar refinement), not a closed-form guidance law; a real vehicle would need an
  onboard-computable cutoff trigger, which is out of scope here.
- The circularization burn is fully idealized (impulsive, unconstrained thrust, not
  drawn from the vehicle's own propellant) — finite-thrust upper-stage design is
  explicitly excluded from M4.
- The 15 km altitude tolerance and 1500 m/s circularization-delta-v cap are documented,
  reasoned engineering choices, not derived from a specific mission requirement.
- Inclination remains fixed at 28.5° throughout M4 (no inclination trade — that is M5).
- As in M2/M3, drag's direction (not just magnitude) is taken anti-parallel to inertial
  velocity rather than to the relative-wind vector; unchanged, not revisited in M4.
