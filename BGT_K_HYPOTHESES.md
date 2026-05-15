# BGT_K — Second-order Karamata representation (Hawkes Thm 3.2)

## Phase 0 — verbatim statement

**Source:** Hawkes 2RV paper Theorem 3.2 + Theorem 3.3. Lines 404-693 of `arxiv_2311.02655_Second_Order_Regular_Variation_Hawkes.txt`.

> **Theorem 3.2 (Second-order Karamata representation theorem).** For α ∈ R, ρ ≤ 0 and A ∈ A_∞^ρ, [a function F is in 2RV_∞^{α,ρ}(A) if and only if it admits a representation involving the integrated auxiliary function A].

[Full statement requires reading more pages; extract captures the headline.]

> **Theorem 3.3 (Extended second-order Karamata theorem).** For α ≠ 0, ρ ≤ 0 and A ∈ A_∞^ρ, [extension to integral transforms ∫_{t_0}^t s^θ F(s) ds].

## Hypothesis types

- h_1: F is 2RV (first-order RV + second-order limit condition with auxiliary A).
- h_2: ρ ≤ 0, A ∈ A_∞^ρ (auxiliary slowly varying with index ρ ≤ 0, decays to 0).

## Phase 1

Same as candidate D: |ε_k| is rapidly (not regularly) varying; L=|ε_k|·2^k fails slow-variation at k=7 jump.

**Phase 1 verdict: NO_FIT** — fails at h_1 (first-order RV prerequisite for 2RV).

## Phase 2

If h_1 had held, the representation would refine the prefactor of L(k) at second order, identifying a *slow secular drift* in L within plateau. The drift would be encoded in A(k), itself slowly varying with index ρ ≤ 0.

This would be useful for a structural picture but does NOT deliver polynomial-in-A bound. The (1/2)^n rate is first-order; 2RV refines prefactor, not rate.

## Phase 3

Same as D — single-regime, k=7 jump structurally blocked.

## Disposition: NO_FIT (inherits D's failure at first-order RV)
