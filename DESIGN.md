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
