"""Physical constants and the benchmark protocol shared by every stage."""

from __future__ import annotations

AVOGADRO = 6.02214076e23  # 1/mol
KJ_PER_MOL_PER_EV = 96.48533212
GAS_CONSTANT_KJ_PER_MOL_K = 8.314462618e-3

# The ml-peg protocol for every box: NPT equilibration (discarded), then
# production. 1 ns of production rather than 2 ns, for cost. See DECISIONS.md
# (boxes and protocol).
EQUILIBRATION_PS = 200.0
PRODUCTION_PS = 1000.0
# The Sage reference runs are longer than the protocol, so the power analysis can
# compare it with 2 ns of production. Every Sage summary uses only the protocol's
# window (EQUILIBRATION_PS to EQUILIBRATION_PS + PRODUCTION_PS).
SAGE_RUN_PS = 2200.0
