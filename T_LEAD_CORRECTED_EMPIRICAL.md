# T_LEAD_CORRECTED_EMPIRICAL — empirical cross-check at ε_1..ε_13

**Date:** 2026-05-12. Phase 4 of T_lead corrected-rate probe.

---

## 1. ε_n trajectory (sources)

| k | ε_k | source |
|---|---|---|
| 1 | +2.000000e−01 = +1/5 | exact Fraction (R77.7 v2 cache) |
| 2 | +9.523810e−03 = +1/105 | exact Fraction |
| 3 | −5.091986e−03 | exact Fraction |
| 4 | −2.452258e−03 | exact Fraction |
| 5 | −1.151747e−03 | exact Fraction |
| 6 | −4.979057e−04 | exact Fraction |
| 7 | −1.175237e−03 | **exact Fraction (R77.7 v2, 2026-05-12)** |
| 8 | −7.455464e−04 | float64 power-iter on K_8 (4374 states), 1e-15 cross-check |
| 9 | −7.520257e−06 | float64 power-iter on K_9 (13122 states), 1e-15 cross-check |
| 10 | +7.207509e−04 | float64 matrix-free + eigs + FFT (39366), 1e-15 |
| 11 | +1.501967e−03 | float64 + FFT (118098), 1e-15 |
| 12 | +2.274714e−03 | float64 + FFT (354294), 1e-15 |
| 13 | +2.948247e−03 | FFT on cached pi_13_truncated, 1e-19 |

## 2. Sign pattern + envelope

Sign: `+ + − − − − − − − + + + +` (single zero-crossing at k=9→10).

Envelope |ε_n|·2^n:
- k=2..6: 0.0381, 0.0407, 0.0392, 0.0368, 0.0319 (nearly stable near 0.04 — the OLD rate-1/2 envelope reading)
- k=7..9: 0.1504, 0.1908, 0.00385 (the 4.7× jump at k=7; near-zero at k=9)
- k=10..13: 0.7380, 3.0762, 9.3199, 24.16 (steadily growing — NOT bounded)

The envelope is **GROWING** at k ≥ 10. This means |ε_n| at large n does NOT decay faster than (1/2)^n — it grows in absolute terms relative to (1/2)^n envelope.

## 3. Ratios |ε_n / ε_{n-1}|

| k | ratio | value |
|---|---|---|
| 2 | ε_2/ε_1 | 0.04762 |
| 3 | ε_3/ε_2 | 0.5347 |
| 4 | ε_4/ε_3 | 0.4816 |
| 5 | ε_5/ε_4 | 0.4697 |
| 6 | ε_6/ε_5 | 0.4322 |
| 7 | ε_7/ε_6 | **2.3604** |
| 8 | ε_8/ε_7 | 0.6344 |
| 9 | ε_9/ε_8 | 0.01009 |
| 10 | ε_10/ε_9 | 95.85 |
| 11 | ε_11/ε_10 | 2.0840 |
| 12 | ε_12/ε_11 | 1.5145 |
| 13 | ε_13/ε_12 | 1.2961 |

Decelerating ratios at k=11..13 (2.084 → 1.515 → 1.296) — approaching some limit from above. Geometric mean of k=11..13 ratios:

  geomean(2.084, 1.515, 1.296) = exp((ln 2.084 + ln 1.515 + ln 1.296)/3)
                                ≈ exp((0.7345 + 0.4153 + 0.2594)/3)
                                ≈ exp(0.4697)
                                ≈ **1.599**

So |ε_n/ε_{n-1}| in late trajectory ≈ 1.6, much larger than 43/45 ≈ 0.956.

## 4. Hadamard radius (|ε_n|^(1/n))

| k | |ε_k|^(1/k) |
|---|---|
| 7 | 0.3753 |
| 8 | 0.4360 |
| 9 | 0.2752 |
| 10 | 0.4847 |
| 11 | 0.5543 |
| 12 | 0.6045 |
| 13 | 0.6388 |

|ε_n|^(1/n) is GROWING at k=10..13 (matches PADE_NUMERICAL_DISPOSITION's Hadamard reading 2.06→1.81→1.66→1.57 for the inverse, the radius of convergence of E(z)).

**At k=13, |ε_n|^(1/n) ≈ 0.64**, NOT 43/45 ≈ 0.956 and NOT 0.984.

## 5. Comparison to T_lead prediction λ = 43/45 ≈ 0.9556

If T_lead's eigenvalue 43/45 were the asymptotic rate, we'd expect:
- |ε_n/ε_{n-1}| → 43/45 ≈ 0.9556 (or its negative, with sign alternation)
- |ε_n|^(1/n) → 43/45 ≈ 0.9556

**Neither holds in the data through k=13.**

Late |ε_n/ε_{n-1}| ≈ 1.6 (74% gap from 0.956)
Late |ε_n|^(1/n) ≈ 0.64 (33% gap from 0.956)

Both diagnostics MISS the T_lead prediction of 43/45 by large margins.

## 6. Comparison to candidate ρ ≈ 0.984 (STATE.md slow-mode)

If ρ ≈ 0.984:
- |ε_n/ε_{n-1}| → 0.984 (geomean at k=11..13 is 1.6, gap 0.6 = 62%)
- |ε_n|^(1/n) → 0.984 (at k=13 we have 0.64, gap 0.35 = 35%)

ρ = 0.984 ALSO misses the empirical data substantially through k=13. The empirical data through k=13 is NOT yet in any single-rate asymptotic regime — neither 43/45, nor 0.984, nor 1/2.

## 7. The Padé/Hadamard reading (PADE_NUMERICAL_DISPOSITION)

PADE_NUMERICAL_DISPOSITION lands H_TWO_SINGULARITIES_VISIBLE with leading singularity radius at |z| ≈ 1.57 at n=13 (trending inward). Translated to rate-language:
- Hadamard at n=13 says rate(ε) ≈ 1/1.57 ≈ 0.637 (matches |ε_13|^(1/13) ≈ 0.639 ✓)
- z=2 reading (rate 1/2) is REFUTED as the leading singularity by n=10
- z=1.016 reading (rate ≈ 0.984) is NOT YET supported (Hadamard at n=13 hasn't reached 1/0.984 ≈ 1.02 — it's still at 1.57)

So the EMPIRICAL data through n=13 says the leading-singularity radius is at |z| ≈ 1.5..1.7 INSIDE z=2 but OUTSIDE the predicted slow-mode z≈1.02. The data is in a TRANSIENT regime.

## 8. So where does T_lead's 43/45 fit?

T_lead's eigenvalue 43/45 ≈ 0.956 corresponds to singularity at |z| ≈ 1/0.956 ≈ 1.046.

This is CLOSE to but DIFFERENT FROM the predicted slow-mode at z ≈ 1.016 (rate 0.984). Both are in the [1.0, 1.1] zone of z-radius.

**Critical observation:** the Hadamard at n=13 is 1.57 and TRENDING INWARD (n=10..13 → 2.06, 1.81, 1.66, 1.57). At this rate of inward trend, when does Hadamard reach 1.046 (T_lead's prediction)?

Roughly: log(radius) decreases linearly in (some n-related quantity). Fitting the four points:
- n=10: r=2.06
- n=11: r=1.81
- n=12: r=1.66
- n=13: r=1.57

log(r) values: 0.723, 0.594, 0.508, 0.452. Differences: −0.129, −0.086, −0.056 — DECELERATING decay. So Hadamard is decelerating, suggesting it's approaching an asymptote.

Linear fit to log(r) vs 1/n: at large n, predict r → some limit. With only 4 points, projections are uncertain, but the decelerating decay suggests the asymptotic r is NOT 1.016 (would require slowing decay to that value) — it's possibly larger.

**Tentative reading:** the empirical late-trajectory radius might actually be heading toward 1.046 (T_lead's prediction) rather than 1.016. The 43/45 from cross-freq machinery may be a BETTER predictor of the long-term asymptote than the prior-session 0.984 fit.

But this is interpretation under heavy uncertainty — the data through n=13 doesn't distinguish 1.016 from 1.046 reliably (both are within the noise of where the Hadamard trend is heading).

## 9. Honest verdict for Phase 4

**T_lead's prediction of eigenvalue 43/45 ≈ 0.956 (radius 1.046) is in the same ballpark as the empirical asymptote direction, but the data through n=13 is too transient to confirm any specific value.**

What we CAN say:
- 1/2 (R77.3 falsified): empirically refuted (Hadamard already inside z=2 by n=10).
- ρ ≈ 0.984 (STATE prior-session two-mode fit): consistent with overall direction, but n=13 data hasn't reached this asymptote.
- 43/45 ≈ 0.956 (T_lead from cross-freq machinery): consistent with overall direction, slightly closer to z=1 from the empirical trend than 0.984.
- ~1.6 (instantaneous late-trajectory geomean ratio): NOT an asymptotic rate — it's a transient.

The empirical data is currently in a TWO-SINGULARITY transient (PADE_NUMERICAL_DISPOSITION) and cannot CLEANLY decide between 43/45 and 0.984. Both are within reasonable distance of where the data is heading.

## 10. The exact ε_7 cross-check

R77.7 v2 produced ε_7 as exact Fraction. The float reproduces |ε_7| = 1.175e−3 exactly at machine precision (cross-checked against the prior power-iteration). |ε_7|·2^7 = 0.1504, matching the |ε_n|·2^n envelope's 4.7× jump from k=6.

This confirms the empirical reading is precise — no rounding ambiguity at k=7. The transient regime extends through at least k=8..13.

## 11. Bottom-line empirical comparison

|Predicted rate | Source | Predicted |ε_n/ε_{n-1}| | Late empirical (k=11..13 geomean) | Gap |
|---|---|---|---|---|
|1/2 = 0.5 | R77.3 (FALSIFIED) | 0.5 | 1.599 | 220% |
|43/45 = 0.9556 | T_lead corrected | 0.9556 | 1.599 | 67% |
|0.984 | STATE prior-session | 0.984 | 1.599 | 63% |
|0.827 | order-3 recurrence (STATE Tier 3) | 0.827 | 1.599 | 93% |

NONE of the candidate rates fit the LATE empirical ratios (because the data is transient, not asymptotic). The Hadamard radius is the right diagnostic, and at n=13 it predicts rate ≈ 0.64 — strictly transient.

The decelerating-inward trend of Hadamard is qualitatively consistent with both 43/45 and 0.984 as long-term asymptotes, but n=13 data does NOT distinguish them.
