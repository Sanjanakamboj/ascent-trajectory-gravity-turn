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


# ============================================================================
# Milestone 4: orbit-capable STUDY vehicle -- explicitly distinct from
# BASELINE_VEHICLE above. See DESIGN.md M4 S3 for the full sizing rationale.
#
# BASELINE_VEHICLE (M1-M3) remains completely unchanged and is never modified or
# replaced by anything below. Nothing in M4 overwrites it.
# ============================================================================

@dataclass(frozen=True)
class M4VehicleDesign:
    """Fixed structural/propulsive parameters of the M4 "orbit-capable study vehicle".

    Only ``payload`` varies across the M4 payload sweep (see ``m4_vehicle``); m0 is
    DERIVED from payload, not held fixed (DESIGN.md M4 S6 -- holding m0 fixed while
    varying payload is an explicitly documented bookkeeping bug to avoid).

    Differences from BASELINE_VEHICLE (M1-M3), and why (DESIGN.md M4 S3):
      - Isp raised from 300 s to 450 s: representative of a higher-performance
        (LH2/LOX-class) propulsion assumption; needed because at Isp=300s no
        physically plausible structural mass fraction closes the LEO delta-v budget
        (M1 S7.4). This is the dominant, explicitly justified change.
      - Propellant mass fraction raised sharply (dry mass cut from 80,000 kg to
        50,000 kg, propellant raised from 410,000 kg to 440,000 kg, at the same
        500,000 kg gross liftoff mass reference) -- an aggressive but documented
        ~10% dry-mass-fraction-of-stack assumption, comparable to serious SSTO
        study-vehicle proposals (e.g. VentureStar/X-33-class targets), NOT an
        existing operational vehicle. Labeled "study vehicle" throughout for this
        reason.
      - Thrust (7.6 MN) and aerodynamic parameters (reference area 10.75 m^2,
        Cd=0.3) are UNCHANGED from BASELINE_VEHICLE, per the instruction to change
        the smallest reasonable set of parameters and keep the aerodynamic/dynamics
        framework identical.
      - Launch site, target altitude, atmosphere, and Earth-rotation treatment are
        all unchanged.
    """

    m_dry: float = 50_000.0     # kg, structural/engine dry mass (fixed across payload sweep)
    m_propellant: float = 440_000.0  # kg (fixed across payload sweep)
    thrust: float = 7.6e6       # N, unchanged from BASELINE_VEHICLE
    isp: float = 450.0          # s, raised from BASELINE_VEHICLE's 300 s (see above)
    reference_area: float = 10.75  # m^2, unchanged from BASELINE_VEHICLE
    drag_coefficient: float = 0.3  # unchanged from BASELINE_VEHICLE


M4_DESIGN = M4VehicleDesign()
M4_REFERENCE_PAYLOAD = 10_000.0  # kg, same reference payload value as BASELINE_VEHICLE


def m4_vehicle(payload_mass: float) -> Vehicle:
    """Construct the M4 orbit-capable study vehicle for a given payload mass.

    m0 = m_dry + m_propellant + payload_mass -- m_dry and m_propellant are FIXED
    (M4_DESIGN); only payload (and therefore m0) varies. This is the payload-sweep
    convention documented in DESIGN.md M4 S6: adding payload to an otherwise-fixed
    vehicle increases gross liftoff mass, which is what should degrade performance
    as payload increases (the expected monotonic trade, M4 S8 check B) -- silently
    holding m0 fixed while sweeping payload would hide that degradation entirely.
    """
    d = M4_DESIGN
    m0 = d.m_dry + d.m_propellant + payload_mass
    return Vehicle(
        m0=m0, m_propellant=d.m_propellant, m_dry=d.m_dry, m_payload=payload_mass,
        thrust=d.thrust, isp=d.isp, reference_area=d.reference_area,
        drag_coefficient=d.drag_coefficient,
    )


M4_REFERENCE_VEHICLE = m4_vehicle(M4_REFERENCE_PAYLOAD)
