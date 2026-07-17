# T_V_SPECTRUM — not executed

**Date:** 2026-05-12. Phase 3 deliverable, conditional on Phase 2 completing.

## Status: NOT EXECUTED

Phase 2 (T_V_MATRIX_GMAX_*.md) not executed because Phase 1 (T_V_RECURSION.md) landed H_M_RECURSION_UNDERSPECIFIED.

The only block that IS rigorously known: T_diag = (1/5)·[[1, 1], [4, 4]] on the g = 0 slice, spectrum {0, 1} (eigenvectors (1, -1) and (1, 4)). This is rigorously derived in `result_77_T_diagonal.py` and rehearsed in CROSS_FREQ_PHASE1_EXPANSION §6.

T_diag's spectrum {0, 1} does NOT have an eigenvalue at 1/2. The rate-1/2 sits at the off-diagonal correction Off_n, which (per Phase 1) doesn't reduce to a closed operator on V_M.
