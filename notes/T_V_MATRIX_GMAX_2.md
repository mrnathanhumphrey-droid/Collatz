# T_V_MATRIX_GMAX_2 — not executed

**Date:** 2026-05-12. Phase 2 deliverable, conditional on Phase 1 closing.

---

## Status: NOT EXECUTED

Phase 1 (T_V_RECURSION.md) landed at **H_M_RECURSION_UNDERSPECIFIED**. The recursion M_{n+1}^{ab}(g, c) → Σ M_n^{a'b'}(g', c') does not close on V_M^{(g_max)} as defined in cross_freq materials. Two distinct structural obstructions:

1. **Phase offsets θ_{v,g} = 2^v·ẽ_g/3** produced by the unit-shuffle step are generically not of the form ẽ_{G''} for any G''; hence the moments that emerge are outside V_M.

2. **Shift index G = v' + g - v** takes ODD values for g = 2 even (under the surviving parity constraints), generating odd-G moments outside V_M = span{even-g}.

Together, these mean no finite-rank Q matrix on V_M^{(g_max=2)} can represent T_V faithfully. The block structure is:

  - **g = 0 row**: T_diag = (1/5)·[[1,1],[4,4]] on (P_+, P_−) — diagonal block, RIGOROUS (CROSS_FREQ_PHASE1_EXPANSION §6, result_77_T_diagonal.py).
  - **g = 0 → g = 2 coupling**: Off-diagonal feed via the off-diagonal correction Off_{n+1} = (3·Σ_g W_±(g)·X̄_n(c; g)) — RIGOROUS shape (CROSS_FREQ_PHASE1_EXPANSION §7).
  - **g = 2 → g = ? feedback**: BLOCKED. The recursion produces phase-twisted moments at odd G, outside V_M.

The 2x2 block on the (g=0) slice ALONE is just T_diag (rate 1, not 1/2). The full T_V at g_max = 2 would require knowing the g = 2 → g = 2 self-coupling and the g = 2 → g = 0 backward coupling — neither of which is derivable in V_M alone.

**Matrix not constructed.** See T_V_RECURSION.md §4–§9 for the obstruction; T_V_DISPOSITION.md for routing.
