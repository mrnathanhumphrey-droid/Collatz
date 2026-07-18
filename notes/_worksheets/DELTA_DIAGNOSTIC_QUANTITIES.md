# DELTA_DIAGNOSTIC_QUANTITIES — five diagnostic quantities at n=2..6

**Date:** 2026-05-12. Wilson (analyst) — leading-vs-subleading diagnostic, Phase 1 of 3.

## Purpose

Tabulate, at each n ∈ {2,3,4,5,6}, the five quantities specified in the brief:

1. **|ε_n|·2^n** — leading-order-removed coefficient (should → 1/30 if simple-pole exact).
2. **δ_n := |ε_n|·2^n − 1/30** — subleading correction.
3. **sign(δ_n)** — verify sign-flip claim between n=5 and n=6.
4. **δ_n / δ_{n-1}** — consecutive ratios (n=3..6).
5. **δ_n / ε_n** — relative size of correction.

## Source

ε_n from `experiments_output/result_77_7_eps_exact_through_k7.json` (exact rationals, cached). All quantities below are computed with exact Fractions and reported as floats for readability.

For reference:
- ε_2 = 1/105
- ε_3 = −5191/1019445
- ε_4, ε_5, ε_6 = large-numerator rationals (see JSON)

## Quantities

| n | ε_n (float) | \|ε_n\|·2^n | δ_n = \|ε_n\|·2^n − 1/30 | sign(δ_n) | δ_n / δ_{n−1} | δ_n / ε_n |
|---|---:|---:|---:|:---:|---:|---:|
| 2 | +9.5238×10⁻³ | 0.0380952 | +0.00476190 | + | —      | +0.5000 |
| 3 | −5.0920×10⁻³ | 0.0407358 | +0.00740250 | + | +1.5546 | −1.4538 |
| 4 | −2.4523×10⁻³ | 0.0392364 | +0.00590310 | + | +0.7974 | −2.4072 |
| 5 | −1.1518×10⁻³ | 0.0368563 | +0.00352297 | + | +0.5969 | −3.0586 |
| 6 | −4.9791×10⁻⁴ | 0.0318659 | −0.00146744 | **−** | **−0.4165** | +2.9472 |

Exact δ_2 = 4/105 − 1/30 = (24 − 7)/210 ... wait: 4/105 = 24/630; 1/30 = 21/630; δ_2 = 3/630 = **1/210** ≈ 0.00476190.

## Observations

### (a) Leading-order |ε_n|·2^n is approximately 1/30 ≈ 0.03333

Range: [0.0319, 0.0407]. All within ±22% of 1/30. Consistent with R76 §10's prediction that leading coefficient is exactly 1/30. The leading-order claim is *not* contradicted by the data; the deviation is the δ_n subleading correction we are analyzing.

### (b) δ_n confirmed non-monotone with sign-flip between n=5 and n=6

Signs across n=2..6: **+, +, +, +, −**. Sign-flip at n=5→6 confirmed.

δ_n magnitude trajectory: 0.00476 ↗ 0.00740 (peak) ↘ 0.00590 ↘ 0.00352 ↘ |−0.00147|. Single peak at n=3, then monotone decay in magnitude through n=6.

### (c) Consecutive ratios δ_n/δ_{n−1} vary wildly

- δ_3/δ_2 = +1.5546
- δ_4/δ_3 = +0.7974
- δ_5/δ_4 = +0.5969
- δ_6/δ_5 = **−0.4165**

The ratios are NOT converging to any specific value. They trend downward (1.55 → 0.80 → 0.60 → −0.42) and the last step crosses zero. This is incompatible with a single geometric form (c·ρ^n with real ρ) and incompatible with a power-law correction n^{-α}·(1/2)^n (which would give ratios approaching 1/2 from above monotonically). The downward-then-flip pattern is the signature of either (i) a two-mode superposition where modes interfere, (ii) genuine irregularity from non-asymptotic-yet small n, or (iii) a richer singularity structure.

### (d) δ_n / ε_n is large and growing

δ_n / ε_n ranges from 0.5 (n=2) to ~−3 (n=4,5) and back to ~+3 (n=6). The subleading correction is *not* small compared to ε_n at observed n; it is comparable in magnitude to the leading −(1/30)·(1/2)^n contribution. This means R76 §10's −(1/30)·(1/2)^n is not yet the dominant signal at n=2..6 — the data is still in a pre-asymptotic regime where the subleading "correction" δ_n·(1/2)^n competes with the −(1/30)·(1/2)^n leading term.

This is a critical finding: **at n=2..6, the subleading and leading are of similar magnitude, NOT a small perturbation on top of a large leading.** This complicates any ansatz fit, because asymptotic forms assume the leading dominates.

## What this tells us before any ansatz fit

The non-monotone δ_n with sign-flip at n=6 already constrains the candidate forms:

- Single-term **c·ρ^n** with real ρ: cannot match sign flip after monotone-positive run — **incompatible without further structure**.
- Single-term **c·n^{-α}·(1/2)^n** with real α: same problem — **incompatible**.
- Single-term **c·log(n)·(1/2)^n**: same — **incompatible**.
- **Two-term c_1·ρ_1^n + c_2·ρ_2^n**: 4 parameters fit 4 points exactly; held-out n=6 is the test.
- **c·cos(ωn+φ)·(1/2)^n**: 3 parameters; can in principle accommodate sign-flip via the cosine; held-out n=6 is the test.

Phase 2 (ANSATZE.md) runs the formal fits.
