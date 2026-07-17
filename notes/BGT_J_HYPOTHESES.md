# BGT_J — Hazard-rate framework

## Phase 0 — verbatim statement

**Source:** arxiv_2504.11655 Theorem 2.4 (Karamata's UCT), Theorem 2.7 (Karamata representation via hazard rate), Theorem 2.11 (Karamata's theorem), Theorem 2.15 (von Mises). Hazard rate framework around generalized hazard rate h_F(x) = F'(x) / (1 − F(x)).

> **Theorem 2.7 (Karamata's Representation, hazard-rate form).** [Standard Karamata representation with the slowly-varying part encoded via the integrated hazard rate.]

> **Theorem 2.15 (von Mises condition).** Let g : R+ → R+ be differentiable. If the von Mises condition lim_{t→∞} t · g'(t) / g(t) = ρ holds, then g is regularly varying with index ρ.

## Hypothesis types

- h_1: F is a CDF (probability distribution function).
- h_2: F differentiable, hazard rate h_F well-defined.
- h_3 (von Mises): t · g'(t) / g(t) → ρ.

## Phase 1

| hyp | check | verdict |
|---|---|---|
| h_1 | ε_k not a CDF | **FAILED categorically** (same as H). |
| h_3 (apply von Mises to L=|ε_k|·2^k as a continuous interpolation) | k · L'(k) / L(k): take L'(k) ≈ L(k+1) − L(k). Values: k=2: 2·(+0.0026)/0.038 ≈ 0.14; k=3: 3·(−0.0015)/0.041 ≈ −0.11; k=4: 4·(−0.0024)/0.039 ≈ −0.25; k=5: 5·(−0.005)/0.037 ≈ −0.67; k=6: 6·(+0.119)/0.032 ≈ +22.3 (jump). | **FAILED at k=6 (jump term explodes).** Within plateau (k=2..5) the von Mises quotient drifts negative, doesn't converge to a fixed ρ. |

**Phase 1 verdict: NO_FIT.**

## Disposition: NO_FIT (categorical at h_1, jump at h_3)

Hazard-rate framework operates on CDFs of random variables. ε_k is not a CDF, and the von Mises condition fails at k=6 jump on any reasonable interpolation.
