"""Ascent trajectory / gravity-turn simulation package.

Milestone 1 documented the mission definition, governing equations, and conventions in
``DESIGN.md`` (no code). Milestone 2 (current) implements and verifies the reusable
physics foundation: ``atmosphere``, ``propulsion``, ``dynamics`` (planar point-mass
ascent ODE), and ``controls`` (prescribed, non-guided thrust-direction profiles used only
to exercise the dynamics). No gravity-turn guidance law, payload optimization, or
inclination sweep is implemented yet -- those begin in later milestones.
"""

__version__ = "0.2.0"
