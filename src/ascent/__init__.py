"""Ascent trajectory / gravity-turn simulation package.

Milestone 1 documented the mission definition, governing equations, and conventions in
``DESIGN.md`` (no code). Milestone 2 implemented and verified the reusable physics
foundation (``atmosphere``, ``propulsion``, ``dynamics``, ``controls``). Milestone 3
added a verified gravity-turn guidance law and orbital diagnostics (``orbital``) on the
unchanged M1 vehicle, which does not reach orbit. Milestone 4 (current) defines a
separate, explicitly labeled M4 orbit-capable study vehicle
(``constants.m4_vehicle``), a drag-consistent engine-cutoff search
(``insertion_search``), and delta-v/loss accounting (``loss_budget``), and solves for
maximum payload capability at 400 km / 28.5 deg. No inclination sweep, launch-azimuth
trade, or staging optimization is implemented yet -- those begin in later milestones.
"""

__version__ = "0.2.0"
