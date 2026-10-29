"""Physical constants and the benchmark protocol shared by every stage."""

from __future__ import annotations

AVOGADRO = 6.02214076e23  # 1/mol
KJ_PER_MOL_PER_EV = 96.48533212
GAS_CONSTANT_KJ_PER_MOL_K = 8.314462618e-3

# The ml-peg protocol for every box: NPT equilibration, then production, both
# discarded and kept as in the benchmark analysis. See DECISIONS.md (protocol).
EQUILIBRATION_PS = 200.0
PRODUCTION_PS = 2000.0
