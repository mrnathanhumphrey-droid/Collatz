# BGT_F — Karamata's monotone density theorem

## Phase 0 — verbatim statement

**Source:** Kevei 2019 Theorem 14. Lines 541-558.

> **Theorem 14.** Let U(x) = ∫_0^x u(t) dt ~ c · x^ρ · ℓ(x) as x → ∞ for c ≥ 0, ρ ≥ 0, ℓ slowly varying, and assume that u is ultimately monotone. Then
> u(x) ~ c · ρ · x^{ρ−1} · ℓ(x).

## Hypothesis types

- h_1: U is the integral of u (sequence analog: cumulative sum).
- h_2: U is regularly varying with index ρ ≥ 0.
- h_3 (load-bearing): u is ULTIMATELY MONOTONE.

## Phase 1 — hypothesis × ε_k matrix

ε_k has sign pattern (+, +, −, −, −, −, −, −) and across the jump |ε_k| values are non-monotone (decreasing to k=6 at 2.49×10^{-4}, then jumping to k=7 at 1.18×10^{-3}, k=8 at 7.46×10^{-4}). So |ε_k| is not monotone either.

| hyp | check | verdict |
|---|---|---|
| h_3 | u = ε_k is sign-mixed, non-monotone | **FAILED.** |
| h_3 | u = |ε_k| has decline + jump + decline, non-monotone | **FAILED.** |

**Phase 1 verdict: NO_FIT.**

(Per pre-registered priors: NO_FIT 75% for F was the highest NO_FIT prior. Confirmed.)

## Phase 2 / Phase 3

N/A — Phase 1 fails categorically.

## Disposition: NO_FIT (categorical)

Syracuse ε_k is sign-mixed (k=1,2 positive, k≥3 negative) and |ε_k| is non-monotone at k=7. Monotone density theorem fails at the monotonicity hypothesis directly.
