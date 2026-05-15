# BGT_M — Potter bounds

## Phase 0 — verbatim statement

**Source:** Kevei 2019 Theorem 11. Lines 332-358.

> **Theorem 11 (Potter bounds).** (i) Let ℓ be a slowly varying function. Then for each A > 1, δ > 0 there exists x_0 such that for each x, y ≥ x_0
> ℓ(x) / ℓ(y) ≤ A · max{ (x/y)^δ, (y/x)^δ }.
>
> (iii) If f ∈ RV_ρ then for each A > 1, δ > 0 there exist x_0 > 0 such that for x, y ≥ x_0
> f(x) / f(y) ≤ A · max{ (x/y)^{ρ+δ}, (x/y)^{ρ−δ} }.

## Hypothesis types

- h_1: ℓ slowly varying (or f RV with index ρ).

## Phase 1

| hyp | check | verdict |
|---|---|---|
| h_1 (ℓ = L = |ε_k|·2^k slowly varying) | fails at k=7 jump (same as A) | **FAILED.** |
| h_1 (test L empirically for Potter-bound-style ratio control) | within plateau k=2..6, L(x)/L(y) ratios are in [0.86, 1.07] — clean Potter behavior with A ≈ 1.1, δ ≈ 0.1. Across the k=6→7 boundary: L(7)/L(6) = 4.72, vastly exceeds Potter for any small δ. | **FAILED at jump.** |

**Phase 1 verdict: NO_FIT.**

## Disposition: NO_FIT (downstream of slow-variation failure at k=7)

Potter bounds are derived consequence of slow variation; same failure mode as A.
