"""Physical constants and the M1 baseline scenario/vehicle parameters.

All values are taken verbatim from ``DESIGN.md`` (M1). Nothing here is tuned or altered
for M2 — this module exists so that M2 code and tests share a single, documented source
of truth instead of re-typing numbers.
"""

from dataclasses import dataclass

# --- Earth constants (DESIGN.md S1.2) -----------------------------------------------
MU_EARTH = 3.986004418e14      # m^3/s^2, standard gravitational parameter
R_EARTH = 6_378_137.0          # m, WGS84 equatorial radius (spherical-Earth approximation)
OMEGA_EARTH = 7.2921159e-5     # rad/s, Earth sidereal rotation rate
G0 = 9.80665                   # m/s^2, standard gravity (Isp -> exhaust velocity only)

# --- Launch site (DESIGN.md S1.1) ---------------------------------------------------
LAUNCH_LATITUDE_DEG = 28.5

# --- Target orbit (DESIGN.md S1.3) --------------------------------------------------
TARGET_ALTITUDE = 400e3        # m
BASELINE_INCLINATION_DEG = 28.5

# --- Atmosphere reference values (DESIGN.md S1.6 / M2 atmosphere.py) ----------------
# See atmosphere.py for the actual model and full documentation of these parameters.
ATMOSPHERE_RHO0 = 1.225        # kg/m^3, sea-level reference density
ATMOSPHERE_SCALE_HEIGHT = 8500.0  # m
ATMOSPHERE_H_MAX = 100_000.0   # m, treated as vacuum above this altitude


@dataclass(frozen=True)
class Vehicle:
    """M1 physics-verification baseline vehicle ("Ascent-1").

    This is the *physics-verification baseline* used to exercise and verify the
    dynamics, atmosphere, and propulsion models. It is deliberately **not** claimed to
    be an "orbit-capable vehicle" — DESIGN.md S7.4 documents that its ideal Tsiolkovsky
    delta-v (~5.045 km/s) falls well short of a realistic LEO delta-v budget
    (~8.9-9.5 km/s). Whether/how to close that gap (staging, resizing, etc.) is an
    explicit later design decision, not something silently changed here.
    """

    m0: float = 500_000.0       # kg, initial (gross liftoff) mass
    m_propellant: float = 410_000.0  # kg
    m_dry: float = 80_000.0     # kg, structure + engines, excludes payload
    m_payload: float = 10_000.0  # kg
    thrust: float = 7.6e6       # N, constant
    isp: float = 300.0          # s, constant (altitude-averaged simplification)
    reference_area: float = 10.75  # m^2 (d = 3.7 m)
    drag_coefficient: float = 0.3  # constant, simplified

    @property
    def m_min(self) -> float:
        """Minimum allowed mass: dry mass + payload (propellant fully depleted)."""
        return self.m_dry + self.m_payload

    def __post_init__(self):
        total = self.m_dry + self.m_propellant + self.m_payload
        if abs(total - self.m0) > 1e-6:
            raise ValueError(
                f"Vehicle mass budget inconsistent: m_dry+m_propellant+m_payload="
                f"{total} != m0={self.m0}"
            )


BASELINE_VEHICLE = Vehicle()
