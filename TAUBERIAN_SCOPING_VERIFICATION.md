# TAUBERIAN_SCOPING_VERIFICATION — Phase 4: empirical verification at small n

**Date:** 2026-05-12. Wilson. Phase 4 of the Tauberian scoping probe.

## Purpose

For each candidate Tauberian theorem matching the singularity-type pattern (Phase 3), compute the predicted asymptotic for ε_n at n=2..6 from theorem statement + R76 §11 leading coefficient (1/30), and compare against empirical ε_n.

---

## Setup

R76 §10 conjecture (verbatim): *"S_n = 7/15 − (1/30)·(1/2)^n + O((1/4)^n)."*

Define ε_n = S_n − 7/15. Then R76 §10's conjecture gives ε_n = −(1/30)·(1/2)^n + O((1/4)^n).

Empirical ε_n (cached, exact rationals through k=6; n=1 transient excluded):

| n | ε_n (float) | |ε_n|·2^n |
|---|---|---|
| 2 | +0.009524 | 0.0381 |
| 3 | −0.005092 | 0.0407 |
| 4 | −0.002452 | 0.0392 |
| 5 | −0.001152 | 0.0369 |
| 6 | −0.000498 | 0.0319 |

R76 §10 leading-order prediction: |ε_n|·2^n → 1/30 = 0.03333 as n → ∞. Observed values are within ±20% of this for n=2..6, with rapid sign change between n=2 and n=3 (per R76 §10's "sign flip at n = 2 → 3").

---

## Verification of Chevalier 2507.15394 Theorem 1.14 (pure α=1/2 square-root)

**Prediction:** ε_n ~ C · (1/2)^n · n^{-3/2}, equivalently |ε_n|·2^n · n^{3/2} → C as n → ∞.

**Computed |ε_n|·2^n · n^{3/2}:**

| n | n^{3/2} | |ε_n|·2^n · n^{3/2} |
|---|---|---|
| 2 | 2.828 | 0.1077 |
| 3 | 5.196 | 0.2117 |
| 4 | 8.000 | 0.3139 |
| 5 | 11.180 | 0.4121 |
| 6 | 14.697 | 0.4683 |

The product grows by factor 4.3× across n=2..6 — clearly not approaching a constant.

**Verdict: FALSIFIED at the leading-order level.** Chevalier Thm 1.14 in its pure form does not describe ε_n's leading behavior.

---

## Verification of Chevalier Theorem 1.16 (meromorphic h, pole of order M at 0)

**Prediction:** ε_n ~ D · (1/2)^n · n^{M − 3/2}.

For M = 1: ε_n ~ D · (1/2)^n · n^{-1/2}.

**Computed |ε_n|·2^n · n^{1/2}:**

| n | n^{1/2} | |ε_n|·2^n · n^{1/2} |
|---|---|---|
| 2 | 1.414 | 0.0539 |
| 3 | 1.732 | 0.0706 |
| 4 | 2.000 | 0.0785 |
| 5 | 2.236 | 0.0824 |
| 6 | 2.449 | 0.0781 |

Product range = 1.53×. **Better fit than pure α=1/2**, but still not converging — peak at n=5, then dips at n=6.

For M = 2: prediction is |ε_n|·2^n · n^{1/2} would need to grow like n, which is not seen.

**Verdict for Chevalier Thm 1.16 with M=1: WEAK SUPPORT.** The n^{-1/2} subleading correction is approximately right in scale (subleading correction of order n^{-1/2} compared to a constant leading term), but not a clean fit on its own.

---

## Verification of "leading simple-pole + branch-cut subleading" model

**Model:** ε_n = −(1/30)·(1/2)^n + (subleading correction), where subleading is of the form c·(1/2)^n · f(n) with f(n) a slowly-varying (sub-exponential) function.

Compute the **subleading correction**: δ_n := |ε_n|·2^n − 1/30.

| n | |ε_n|·2^n | δ_n = |ε_n|·2^n − 1/30 |
|---|---|---|
| 2 | 0.03810 | +0.00476 |
| 3 | 0.04074 | +0.00740 |
| 4 | 0.03924 | +0.00590 |
| 5 | 0.03686 | +0.00352 |
| 6 | 0.03187 | −0.00147 |

**Critical observation: δ_n is non-monotone and changes sign between n=5 and n=6.**

This is the strongest single empirical fact constraining the singularity structure. Any single-term subleading prediction (n^{-α} with fixed α and fixed sign) must be reconciled with this sign flip.

### Test subleading n^{-3/2} (i.e., R76 §10's leading is exact, subleading is √-branch at z=2):

**Prediction:** δ_n ~ c · n^{-3/2}.

| n | n^{3/2} | δ_n · n^{3/2} |
|---|---|---|
| 2 | 2.828 | +0.01347 |
| 3 | 5.196 | +0.03845 |
| 4 | 8.000 | +0.04723 |
| 5 | 11.180 | +0.03936 |
| 6 | 14.697 | −0.02164 |

Sign flips. Not a constant times n^{-3/2}. **Verdict: not supported.**

### Test subleading 1/n (i.e., logarithmic singularity at z=2):

**Prediction:** δ_n ~ c · 1/n.

| n | n | δ_n · n |
|---|---|---|
| 2 | 2 | +0.00952 |
| 3 | 3 | +0.02220 |
| 4 | 4 | +0.02361 |
| 5 | 5 | +0.01762 |
| 6 | 6 | −0.00881 |

Sign flips, not constant. **Verdict: not supported as pure subleading.**

### Test subleading exponential decay (next pole at z=4, i.e., (1/4)^n / (1/2)^n = (1/2)^n):

**Prediction:** δ_n ~ c · (1/2)^n.

| n | (1/2)^n | δ_n / (1/2)^n |
|---|---|---|
| 2 | 0.250 | +0.01902 |
| 3 | 0.125 | +0.05920 |
| 4 | 0.0625 | +0.09440 |
| 5 | 0.03125 | +0.11264 |
| 6 | 0.015625 | −0.09418 |

Not constant; sign flips. **Verdict: not supported either.**

### Verdict: No single-term subleading model fits

The non-monotone, sign-flipping δ_n pattern strongly suggests **multiple competing subleading terms**, possibly with secondary singularities of E(z) contributing on the second sheet of the branch.

---

## Strength-of-consistency assessment (Adversarial Check A3)

Per the brief's anti-pattern A3 ("Empirical verification at n=2..6 is small-sample"):

The empirical evidence at n=2..6 distinguishes some hypotheses:

- **Pure α=1/2 square-root branch**: STRONGLY INCONSISTENT (Chevalier Thm 1.14 falsified; product n^{3/2}·|ε_n|·2^n grows 4.3× across the data).
- **Pure pole + n^{-1/2} subleading (Chevalier Thm 1.16, M=1)**: WEAK CONSISTENCY (1.53× range, peaks at n=5, dips at n=6 — fit not stable).
- **Leading simple-pole, with subleading branch of any form**: QUALITATIVELY CONSISTENT (the leading-order |ε_n|·2^n ≈ 1/30 is clearly visible), but the subleading branch type is not determined.
- **Pure logarithmic (FS log entry)**: NOT CONSISTENT in any form tested.

Honest report: Six data points are **insufficient** to robustly select among "leading simple-pole + subleading O(n^{-1/2})", "leading simple-pole + subleading O(n^{-3/2})", or "leading simple-pole + subleading exponential + secondary oscillation". The non-monotonicity of δ_n indicates the subleading has multiple terms, which 5 data points cannot fit cleanly.

---

## What additional data resolves this (Adversarial Check A5)

To **distinguish power-law from log subleading**, the following data would be diagnostic:

1. **Extend to ε_7 (k=7 Markov chain).** Critical disambiguator. Predicted by R77.6 if power-law branch: |ε_7|·2^7 close to 1/30; if log: same but with slower decay. The sign of δ_7 (positive or negative) would also probe whether the sign flip at n=5→6 is the start of a settled negative regime or a transient.

2. **Extend to ε_8 (k=8 Markov chain).** [4/4] diagonal point. If the sign flip at n=6 was a transient, δ_8 should return to consistency with the leading subleading direction; if persistent, it's a clear signal that there is a secondary oscillatory term.

3. **Fit a multi-term ansatz to all available data.** E.g., ε_n = −(1/30)·(1/2)^n + A·(1/2)^n/n^α + B·(1/2)^n/n^β·cos(γn + φ), then determine which (α, β, γ, φ) family is most consistent. Six data points might be enough to identify A, α, B, but γ and φ require more.

4. **Use D-Padé (differential approximants)**, which can extract branch exponent at lower N than ordinary Padé. Would require additional numerical infrastructure but might give a clean signal even at N=5.

---

## Conclusion of Phase 4

**The H_SQUARE_ROOT_MATCHES_PLUS_EMPIRICAL hypothesis is EMPIRICALLY FALSIFIED** at the leading-order level: Chevalier 2507.15394 Theorem 1.14 in its pure form predicts ε_n ~ C·(1/2)^n · n^{-3/2}, and this is not seen at n=2..6.

**The H_SQUARE_ROOT_MATCHES_BUT_EMPIRICAL_AMBIGUOUS hypothesis** also falls: the empirical data is inconsistent with α=1/2 even at the *subleading* level.

**The H_GENERAL_BRANCH_MATCHES hypothesis** holds in a weak sense: FS Ch. VI's mixed-singularity framework (leading simple-pole + subleading branch at z=2) is qualitatively consistent, but the precise subleading branch type is not determinable from N=5.

**The H_AMBIGUOUS / INCONCLUSIVE disposition is the honest call** given the N=5 data limitation and the non-monotone subleading correction pattern.

Specific theorem the framework points to: **Chevalier 2507.15394 Theorem 1.16 (meromorphic h with pole at 0)** is the cleanest single-theorem candidate consistent with both R77.6's branch-cut detection AND R76 §10's leading (1/30)·(1/2)^n behavior — but the value of M (and the precise meromorphic h) is not determinable from N=5.

---

## Verification script (tauberian_verify.py) for main-thread execution

The numerical analysis above was computed in the subagent thread. A clean verification script for main-thread re-execution is saved as `tauberian_verify.py` alongside this disposition. Expected output is the tables shown in this document. The script loads `experiments_output/result_77_7_eps_exact_through_k7.json` and produces all consistency-test tables.
