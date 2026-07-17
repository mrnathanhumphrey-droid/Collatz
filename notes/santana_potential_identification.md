# Santana-framework potential ϕ identification — outcome (c) sparse data

**Date:** 2026-05-03. Brief: identify the specific potential ϕ within Santana 2026's bridge theorem that produces our empirical Gibbs form P(q | j) ∝ exp(α(j)·q). Numerical: `santana_potential_identification.py`. CSV: `santana_potential_data.csv`.

## 1. Empirical inputs

From Result 34:
- α(2) ≈ 0
- α(4) = −3.02
- α(5) = −2.30
- P(j=2) = 0.9379, P(j=4) = 0.0237, P(j=5) = 0.0379
- m_j = (4^j − 1)/3: m_2=5, m_4=85, m_5=341
- ⟨v|j⟩: 2.0 (j=2), 2.146 (j=4), 2.05 (j=5)

j=3 and j=6 forbidden by number-theoretic constraint (Result 7); next allowed j values are {7, 8, 10, 11, ...}.

## 2. Key diagnostic: α(j) is NON-MONOTONE

α(j) trajectory: 0 → −3.02 → −2.30 from j=2 → 4 → 5. The direction reverses between j=4 and j=5.

**This non-monotonicity rules out every monotone-in-j candidate:**

| Candidate | Fit | R² (3 pts) | Verdict |
|---|---|---|---|
| Linear: a·j + b | best fit | trivially fits 2 pts | not monotone-consistent |
| Logarithmic: a·log(j) + b | — | — | monotone, fails non-monotonicity |
| Inverse: a/j + b | — | — | monotone |
| log(m_j): a·log(m_j) + b | — | — | monotone (m_j strictly increasing in j) |
| log(P(j)): a·log(P(j)) + b | — | — | P(j) generally decreasing in j |

All monotone forms FAIL because the empirical α(j) is non-monotone.

## 3. Non-monotone candidates tested

### j mod 3 cycle hypothesis: FALSIFIED

Since j ≡ 0 mod 3 is forbidden, allowed j cycles through mod 3 = {2, 1, 2, 1, 2, ...}. If α(j) depended only on j mod 3:
- j=2 (mod 3 = 2): α = 0
- j=5 (mod 3 = 2): α should also = 0; empirical = −2.30 ❌

Falsified. j mod 3 alone doesn't determine α(j).

### Single-parameter forms with α(2) ≡ 0 baseline

Fit one parameter c on j=4, predict α(5):

| Form | c | α(5) pred | gap |
|---|---|---|---|
| c·log(m_j/m_2) | −1.066 | −4.50 | −2.20 |
| **c·log(P(2)/P(j))** | **−0.821** | **−2.63** | **−0.33** ← best |
| c·(j−2) | −1.510 | −4.53 | −2.23 |
| c·log(j/2) | −4.357 | −3.99 | −1.69 |
| c·j(j−2) | −0.378 | −5.66 | −3.36 |
| c·(j%3 − 2) | +3.020 | 0.00 | +2.30 |

**Best single-parameter candidate: α(j) ≈ −0.82·log(P(2)/P(j)).** Predicts α(5) = −2.63 vs empirical −2.30 (gap 0.33). Within ~15% but not within sampling precision.

### Three-parameter fits

α = a·log(P(j)) + b·log(m_j) + c fits 3 data points exactly by construction (3 params, 3 points) — no validation possible. Not informative.

## 4. Reverse-engineered ψ(j) from α(j) (Step 2)

ψ(j) = −log Z(j) where Z(j) = Σ_q exp(α(j)·q) over band midpoints {0.125, 0.375, 0.625, 0.875, 0.975}:

| j | α(j) | Z(j) | ψ(j) = −log Z(j) |
|---|---|---|---|
| 2 | 0.00 | 5.000 | −1.6094 |
| 4 | −3.02 | 1.283 | −0.2492 |
| 5 | −2.30 | 1.650 | −0.5005 |

ψ(j) is fully determined by α(j); has no independent information. Also non-monotone (since α(j) is).

## 5. Tao K_h = 3/log(4/3) connection (Step 5)

Test: α(j) / (log(m_j)/log(4/3)) should be constant if ϕ is K_h-related:

| j | log(m_j)/log(4/3) | α(j) | ratio |
|---|---|---|---|
| 2 | 5.59 | 0.00 | 0.000 |
| 4 | 15.44 | −3.02 | −0.196 |
| 5 | 20.27 | −2.30 | −0.114 |

No constant ratio; no clean K_h-related closed form visible.

## 6. ⟨v|j⟩ cross-check (Step 6)

Esscher tilt of Geom(1/2) on v: solve E_w[v] = ⟨v|j⟩_emp for w_j:

| j | ⟨v|j⟩ | implied w_j |
|---|---|---|
| 2 | 2.000 | +0.000 |
| 4 | 2.146 | +0.095 |
| 5 | 2.050 | +0.035 |

w_j range is tiny (~0 to +0.095), while α(j) range is large (~−3 to 0). Different scales, no obvious linear relation.

**This is consistent with Result 36 follow-up 3** (v_t conditionally independent within σ-band). v|j and σ-band|j are nearly orthogonal observables — the j-class predicts σ-band membership strongly (via α(j)) but barely affects the per-step v-distribution within band (via w_j ≈ 0).

The framework MATCH at the q-marginal (Gibbs form P(q|j)) does not extend to the v-marginal (P(v|j) is essentially Geom(1/2) regardless of j).

## 7. Hölder continuity (Step 3) — UNDETERMINED

ϕ(j, q) = α(j)·q + ψ(j) is bounded across j ∈ {2, 4, 5} (range of α: [−3.02, 0]). For higher j:
- IF α(j) bounded across all j: ϕ bounded → standard Pollicott-Urbański framework applies.
- IF α(j) → −∞ as j → ∞: ϕ unbounded → needs Bilbao-Lucena discontinuous-map extensions.

Cannot determine from current data resolution. **NEEDS MORE j-STRATIFIED DATA.**

## 8. Verdict: outcome (c)

**Three data points (j=2, 4, 5) is INSUFFICIENT to identify a non-monotone functional form for α(j) cleanly.**

What's established:
- Non-monotonicity of α(j) is a real structural fact (not noise)
- α(j) ≈ −0.82·log(P(2)/P(j)) is the closest single-parameter fit (within ~15%)
- ψ(j) is determined fully once α(j) is known (via normalization)
- The framework match holds at the σ-band marginal level but not at the v-marginal level
- Tao K_h has no clean ratio with α(j)

What's NOT established:
- Closed-form α(j) (3-param fits are exact-by-construction, not validating)
- Hölder/boundedness of ϕ across all j
- Equilibrium-state variational verification (requires Z(j) values across more j)

### To resolve to outcome (a) or (b)

Generate σ-band-stratified P(q|j) statistics for j ∈ {7, 8, 10, 11}. Each adds one data point. With 7 total:
- Non-monotone candidates with one cycle-parameter (e.g., α(j) = f(j mod 3) + smooth correction) become identifiable
- Asymptotic behavior of α(j) as j → ∞ becomes visible
- Hölder/boundedness check becomes possible
- ψ(j) and Z(j) closed-form candidates can be cross-validated

Until this data exists: framework match is at the **vocabulary level** ("our findings ARE Gibbs equilibrium states for some ϕ within Santana's framework") but NOT at the **technical level** ("the specific ϕ is [closed form]").

## 9. Honest scope statement

The brief estimated 3-4 hours for full identification + Hölder + variational + K_h + ⟨v|j⟩ verification. With only 3 data points, Steps 1-6 collectively yield outcome (c) within ~30 minutes; the remaining 2.5+ hours of variational/Hölder verification is moot until more j-data exists. Honest delivery: outcome (c) classification with diagnostic of non-monotonicity, best-fit single-parameter candidate (within 15%), and explicit recommendation for j-stratified data generation.

**For Lagarias/Tao-level discussion:** the Gibbs equilibrium-state interpretation is REAL (P(q|j) ∝ exp(α(j)·q) at R² ≥ 0.994 per Result 34), and identifies our work as residing within Santana's framework formally. The specific ϕ remains under-determined empirically. The best statement to convey: "Within Santana's bridge theorem, our empirical Gibbs form has tilt parameters α(2) ≈ 0, α(4) = −3.02, α(5) = −2.30 — non-monotone in j, with no clean closed-form identifiable from 3 data points. Generating j ∈ {7, 8, 10, 11} would lock the functional form."

## 10. Files

- `santana_potential_identification.py` — fitting and diagnostic code
- `santana_potential_data.csv` — j-stratified α, ψ, m_j, P_j, ⟨v|j⟩, log(m_j)/log(4/3)
- `santana_potential_log.txt` — full diagnostic log
- `santana_potential_identification.md` — this document
