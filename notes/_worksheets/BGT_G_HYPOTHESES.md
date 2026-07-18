# BGT_G — Karamata's integration theorem (Karamata's Lemma)

## Phase 0 — verbatim statement

**Source:** Kevei 2019 Theorem 12 (Karamata's theorem, direct part). Lines 467-475.

> **Theorem 12 (Karamata's theorem, direct part).** Let f ∈ RV_ρ be locally bounded on [a, ∞). Then
> (i) for σ ≥ −(ρ + 1)
> x^{σ+1} f(x) / ∫_a^x t^σ f(t) dt → σ + ρ + 1;
> (ii) for σ < −(ρ + 1)
> x^{σ+1} f(x) / ∫_x^∞ t^σ f(t) dt → −(σ + ρ + 1).

(Companion converse: Theorem 13.)

## Hypothesis types

- h_1: f ∈ RV_ρ — first-order regularly varying with some index ρ.
- h_2: f locally bounded on [a, ∞).

## Phase 1 — hypothesis × ε_k matrix

Same h_1 failure as candidate D: |ε_k| is not first-order RV (it's *exponentially* decaying, hence rapidly varying); and L(k)=|ε_k|·2^k fails slow variation at k=7 jump.

| hyp | check | verdict |
|---|---|---|
| h_1 (f = |ε_k|) | rapidly varying, not RV | **FAILED categorically.** |
| h_1 (f = L) | slow variation fails at k=7 | **FAILED at jump.** |
| h_2 | locally bounded — yes (only 8 points) | SATISFIED. |

**Phase 1 verdict: NO_FIT.**

## Phase 2

If h_1 had held with f = L slowly varying (ρ = 0), the conclusion at σ = 0 would give: ∫_0^k L(t) dt ~ k · L(k). Equivalently the partial sum Σ_{j≤k} L(j) ~ k · L(k). Empirically Σ_{j≤6} L(j) ≈ 0.187, 6·L(6) ≈ 0.191 — yes, within plateau this holds well. But it fails at k=7 (Σ_{j≤7} L(j) ≈ 0.337 vs 7·L(7) ≈ 1.053).

This is informative: Karamata's lemma "almost" holds within plateau but fails at jump. Useful for chain-side consistency check, but doesn't deliver closure.

## Phase 3

Same as A and D — single-regime, blocked at k=7.

## Disposition: NO_FIT

Same structural reason as D — first-order RV hypothesis fails for ε_k.
