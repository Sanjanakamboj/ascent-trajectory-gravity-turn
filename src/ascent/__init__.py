"""Ascent trajectory / gravity-turn simulation package.

Milestone 1 documented the mission definition, governing equations, and conventions in
``DESIGN.md`` (no code). Milestone 2 implemented and verified the reusable physics
foundation (``atmosphere``, ``propulsion``, ``dynamics``, ``controls``). Milestone 3
added a verified gravity-turn guidance law and orbital diagnostics (``orbital``) on the
unchanged M1 vehicle, which does not reach orbit. Milestone 4 defined a
separate, explicitly labeled M4 orbit-capable study vehicle
(``constants.m4_vehicle``), a drag-consistent engine-cutoff search
(``insertion_search``), and delta-v/loss accounting (``loss_budget``), and solved for
maximum payload capability at 400 km / 28.5 deg. Milestone 5 added direct-ascent
launch-azimuth/inclination geometry (``inclination``) and a payload-vs-inclination
trade across the launch site's full direct-ascent range (28.5-90 deg). Milestone 6
(final) is a documentation/audit milestone: independent CSV reproduction, a
guidance-search sensitivity study on the M5 non-monotonicity, and README/DESIGN.md
consolidation -- no new trajectory physics. No staging, dogleg, on-orbit plane change,
or retrograde-launch trade is implemented.
"""

__version__ = "0.2.0"
