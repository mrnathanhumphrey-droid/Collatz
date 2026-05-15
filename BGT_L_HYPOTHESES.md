# BGT_L — Extended second-order Karamata theorem (Hawkes Thm 3.3)

## Phase 0 — verbatim statement

**Source:** Hawkes 2RV paper Theorem 3.3 (extended second-order Karamata). Line 618 of `arxiv_2311.02655_Second_Order_Regular_Variation_Hawkes.txt`.

> **Theorem 3.3 (Extended second-order Karamata theorem).** For α ≠ 0, ρ ≤ 0 and A ∈ A_∞^ρ, [the integral transform of F has explicit 2RV structure with auxiliary function modified by an explicit factor].

## Hypothesis types

- h_1: F ∈ 2RV (inherits from D / K).
- h_2: α ≠ 0 (non-zero first-order index).

## Phase 1

| hyp | check | verdict |
|---|---|---|
| h_1 | inherits from D, K: fails | **FAILED.** |
| h_2 | requires first-order RV with non-zero index | L is slowly varying (index 0, if anywhere) | FAILED (excludes our regime). |

**Phase 1 verdict: NO_FIT.**

## Disposition: NO_FIT (inherits D, plus h_2 excludes index 0 regime)
