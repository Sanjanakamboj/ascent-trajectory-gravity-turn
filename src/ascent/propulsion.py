"""Propulsion / mass-flow model (DESIGN.md S3, S7).

Pure functions only — no hidden state, no reference to a specific vehicle object beyond
what is passed in explicitly.
"""

from .constants import G0


def mass_flow_rate(thrust: float, isp: float, g0: float = G0) -> float:
    """Propellant mass-flow rate [kg/s], returned as a negative number (mass loss).

    mdot = -T / (Isp * g0)
    """
    return -thrust / (isp * g0)


def exhaust_velocity(isp: float, g0: float = G0) -> float:
    """Effective exhaust velocity ve = Isp * g0 [m/s]."""
    return isp * g0


def thrust_acceleration(thrust: float, mass: float) -> float:
    """Thrust acceleration magnitude a_T = T / m [m/s^2]."""
    return thrust / mass


def mass_at_time(m0: float, thrust: float, isp: float, t: float, g0: float = G0) -> float:
    """Mass after burning at constant (thrust, isp) for duration t (t >= 0), unclamped.

    m(t) = m0 - |mdot| * t
    """
    mdot = mass_flow_rate(thrust, isp, g0)
    return m0 + mdot * t


def burn_duration(propellant_mass: float, thrust: float, isp: float, g0: float = G0) -> float:
    """Time to deplete ``propellant_mass`` at constant (thrust, isp)."""
    mdot = mass_flow_rate(thrust, isp, g0)
    return propellant_mass / (-mdot)


def tsiolkovsky_delta_v(isp: float, m0: float, mf: float, g0: float = G0) -> float:
    """Ideal rocket equation: delta-v = Isp * g0 * ln(m0 / mf)."""
    import math

    return isp * g0 * math.log(m0 / mf)


def clamped_mass_flow_rate(thrust: float, isp: float, mass: float, m_min: float,
                            g0: float = G0) -> float:
    """Mass-flow rate, forced to zero once ``mass`` has reached ``m_min``.

    This is the safety clamp referenced in DESIGN.md's M2 section: integration must
    never be allowed to consume propellant below the vehicle's dry+payload mass. Once
    ``mass <= m_min`` the returned mdot is exactly 0.0, regardless of the requested
    thrust — the caller (dynamics.py) is expected to also drop the thrust term itself
    (a mass-depleted vehicle produces no thrust), but this function independently
    guarantees mass cannot go negative-propellant even if that is forgotten.
    """
    if mass <= m_min:
        return 0.0
    return mass_flow_rate(thrust, isp, g0)
