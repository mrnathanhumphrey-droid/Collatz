# BGT_H — Goldie subexponential / regularly varying tail

## Phase 0 — verbatim statement

**Source:** Jessen-Mikosch 2007 §3 (Karamata's Tauberian for tail distributions). Sample-style: "F is subexponential if F^{*2}(x) / F(x) → 2 as x → ∞." Operates on distribution function tails.

Standing definition (per Bingham 1988 + Jessen-Mikosch): X has regularly-varying tail with index −α if P(|X| > x) = x^{−α} · ℓ(x) for slowly-varying ℓ.

## Hypothesis types

- h_1: F is a probability distribution function (P(X ≤ x)) with X ≥ 0.
- h_2: F̄(x) := 1 − F(x) is regularly varying with index −α < 0.

## Phase 1

| hyp | check | verdict |
|---|---|---|
| h_1 | ε_k is not a probability distribution; |ε_k|·2^k is bounded by ~0.2, not normalizable to a distribution in the standard sense | **FAILED categorically.** |
| h_2 | |ε_k| decays exponentially (geometric ~2^{-k}), not power-law (subexponential or heavy-tailed) | **FAILED categorically.** |

**Phase 1 verdict: NO_FIT (categorical).**

(Per pre-registered priors: NO_FIT 75% — matches.)

## Disposition: NO_FIT (categorical)

Subexponential / heavy-tail framework operates on power-law-decaying tails of distributions; ε_k is light-tailed (geometric decay) and not a distribution. Wrong category.
