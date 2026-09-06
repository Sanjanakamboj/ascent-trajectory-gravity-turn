"""Package-level smoke test.

Solver-behavior tests live in the dedicated test_*.py files added in Milestone 2
(test_atmosphere.py, test_propulsion.py, test_dynamics.py, test_verification.py).
"""

import ascent


def test_package_imports():
    assert ascent.__version__ == "0.2.0"
