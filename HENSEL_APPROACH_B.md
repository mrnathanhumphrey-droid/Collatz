# HENSEL Phase 2 Approach B — Recursive Hensel series with cancellation identity

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

## Disposition: NOT_TRIGGERED (Approach A succeeded provisionally)

Approach B was queued as fallback if Approach A failed to close the closed form. Approach A produced a structurally exact Hensel-lifted closed form by direct digit-wise reduction of the inner Gauss sum (see `HENSEL_APPROACH_A.md`). Approach B is therefore not strictly needed for closure at family level.

However, this document records the **alternative framing** in case Approach A's numerical verification (Phase 3) fails or reveals an unexpected residual, and as a check on Approach A's structural soundness.

## Reformulation of the "recursive Hensel series"

The brief described `δ_k(C_a), ε_k(C_a)` as recursive Hensel-correction series. The Phase 1 articulation collapsed these:

> **δ_k = (next base-p digit of C_a-1, divided appropriately)** — they are just digits of the canonical representative, not abstract recursive series.

Specifically, with `C_a = 1 + p·c_1 + p²·c_2 + p³·c_3 + ... + p^{r-1}·c_{r-1} mod p^r`:
- `s_0 = c_1, s_1 = c_2, s_2 = c_3, ..., s_{r-2} = c_{r-1}`
- `s*(r)(C_a) = s_0 + p·s_1 + p²·s_2 + ... + p^{r-2}·s_{r-2} mod p^{r-1}`

The "δ_k" recursion is trivial digit extraction. There is NO non-trivial recursive series.

## Telescoping/cancellation identity

The original Approach B asked: does the recursive series δ_{k+1} = f(δ_k, ε_k, C_a) telescope?

Under the Phase 1 collapse, this question reformulates as: does the expansion `P_a(s*(r)) = Σ_{j=2}^r (-1)^{j-1} (p·s*(r))^j / (j·(j-1)) mod p^{r+1}` simplify further?

The generating identity `(1+y)·log(1+y) = Σ_{j=2}^∞ (-1)^j · y^j / (j·(j-1))` is itself a "telescoping" of the Taylor series for log via the (1+y) shift. **No further simplification known at the level of formal series.**

In particular:
- The coefficients `1/(j·(j-1))` admit partial-fraction `1/(j-1) − 1/j`. So `Σ_{j=2}^r (-1)^{j-1} y^j / (j(j-1)) = Σ_{j=2}^r (-1)^{j-1} y^j · (1/(j-1) − 1/j)`. Partial sums in this form are alternating telescoping, but the result is just the same `(1+y)·log(1+y)` (a known identity).
- **No hypergeometric identity** mechanically reduces this further. The series `(1+y)·log(1+y)` is non-elementary.

## Cochrane Theorem 2 reapplication

Cochrane Theorem 2 (which proves T78.1 in result_78_FINAL.md) gives complete-sum vanishing under degp H+ = 0. This was a separate construction at q=3.

Re-applying at r ≥ 4 family level: the "polynomial identification" g(u) = c·(1+p)^u − p·a·u in the Cochrane framework has degp H+ structure determined by the principal-unit subgroup. The H+ polynomial mod p has the same "constant non-zero" property as at q=3 (per T78.2_p analog). **Cochrane Theorem 2 does NOT directly produce the Hensel-lifted closed form** — it gives complete-sum vanishing and Fourier sparsity, but the inner Gauss sum G_p(a) still needs the saddle-point evaluation we did in Approach A.

## VMV (deferred to Approach C)

Vinogradov mean-value approach is Approach C's territory.

## Verdict

> **APPROACH_B_NOT_TRIGGERED.** Approach A's direct saddle-point reduction succeeded provisionally; no telescoping or hypergeometric identity beyond the `(1+y)·log(1+y)` series identity is needed at family level. If Approach A's Phase 3 numerical verification fails, B should be revisited with focus on:
> 1. Why the closed-form polynomial isn't capturing G_p(a) exactly — possible: small-prime denominator handling, or hidden higher-order terms.
> 2. Whether a "telescoping over digits" reformulation reveals a cleaner sum-structure.

For now: **NOT_TRIGGERED.** Approach A is the route.
