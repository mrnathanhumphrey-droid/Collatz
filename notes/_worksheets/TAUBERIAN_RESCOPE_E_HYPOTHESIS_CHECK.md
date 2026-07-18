# TAUBERIAN_RESCOPE_E_HYPOTHESIS_CHECK (Alberts 2508.20814 × inputs)

**Date:** 2026-05-13.

## h × I matrix

Try construction: N(X) = Σ_{k ≤ X} ε_k, N̂(X) = Σ_{k ≤ X} |ε_k|. L(s, N) = Σ ε_k k^{-s} (Dirichlet series).

| Hypothesis | (1) ε_k | (2) C1 | (3) C2 | (4) BT |
|---|---|---|---|---|
| h_1: N, N̂ : ℝ≥0 → ℂ with N̂ nondecreasing | SATISFIED (Σ |ε_k| of nonnegative terms). | | | |
| h_2: |N(X)| ≤ N̂(X) | SATISFIED by triangle inequality. | | | |
| h_3: L(s, F) absolutely convergent on Re s > σ_a | UNVERIFIABLE — requires knowing decay rate of ε_k. Empirically ε_k ≈ O(2^{-k}) for k≤6 (then jump at k=7). If sequence is *eventually* O(k^{-c}) for some c > 0, then σ_a = 1 − c is the abscissa. PLAUSIBLE but UNVERIFIABLE without full sequence. | | | |
| h_4: L(s, F) meromorphic continuation to Re s ≥ σ_a − δ with finitely many poles | UNVERIFIABLE — **Mode H circular**. The meromorphic continuation property is essentially the closure target. There is no proven functional equation for the Dirichlet series Σ ε_k k^{-s}. | | | |
| h_5: Twisted-moment bound ∫_T^{2T} |L(σ+it, F) Z^{it}| dt ≪ T^η (log T)^β on left edge σ = σ_a − δ | UNVERIFIABLE — this is the load-bearing analytic input. Without (h_4), there is no L on the left edge to bound, period. | | | |
| h_6: constants σ_a, δ, η̃, η, β exist | UNVERIFIABLE | | | |
| h_7: b = order of pole of L(s, N̂) at s = σ_a | UNVERIFIABLE | | | |

For SELECTED, we'd need:
- A proven analytic continuation of L(s, N) = Σ ε_k k^{-s} (Mode H trap).
- An explicit twisted-moment bound at the left edge.

Neither is supplied by inputs (1)-(4). Input (4) BT actively suggests that the global/adelic structure required for L-function analytic continuation is NOT available for a Syracuse-derived single-place Dirichlet series.

**Theorem E disposition: BLOCKER (h_4 and h_5 Mode H circular).**

---

## Comparison with C and D

Theorem E is an *enhancement* of C and D: it replaces pointwise vertical bounds with twisted-moment bounds, giving better error terms. But the underlying analytic-continuation hypothesis (h_4) is the SAME as in C and D. So E inherits the same Mode H circularity.

**Aggregate disposition for E: BLOCKER.**
