# T_V_MATRIX_GMAX_4 — not executed

**Date:** 2026-05-12. Phase 2 deliverable, conditional on Phase 1 closing.

## Status: NOT EXECUTED

Same obstruction as T_V_MATRIX_GMAX_2. Adding g = 4 to the truncation doesn't help; in fact:

- g = 2: survival pattern (v even, v' odd) → odd-G moments (G = v'+2-v odd).
- g = 4: survival pattern (v odd, v' even) → odd-G moments (G = v'+4-v odd).
- g = 6: survival pattern more complex (ẽ_6 ≡ 0 mod 3 requires mod-9 refinement; expect mixed odd/even-G).

Adding g = 4 introduces more rows but each row points to moments OUTSIDE V_M^{(g_max=4)} = span{g ∈ {0, 2, 4}}. The truncation cannot be closed in any obvious way.

See T_V_RECURSION.md §4–§9 for obstruction.
