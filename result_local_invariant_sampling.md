# Local-Invariant Soil Sampling on Syracuse Trajectories

Date: 2026-05-05.  Ensemble: all odd integers coprime to 3 in [1, 10^7] (excluding the cycle vertex n=1) → **3,333,332 starting points**, iterated under accelerated Syracuse map T(n) = (3n+1)/2^v_2(3n+1) until return to 1 or 10000 step cap.

Bonferroni-corrected significance threshold: p < 0.05/5 = **0.0100** (5 tests).

Total Syracuse steps recorded: 179,787,888.
Trajectories that hit the step cap (didn't return to 1): 0.
Iteration runtime: 28.1 s.

## Verdict summary

| test | invariant | null model | statistic | p-value | verdict |
|---|---|---|---:|---:|:---:|
| T1 | v_2(3n+1) marginal | geometric P(v_2=k) = 2^(-k) | χ²=1504114.5, df=24 | 0.00e+00 | **STRUCTURE** |
| T2 | v_2 \| residue (mod 27) | independence | χ²=27307163.5 | 0.00e+00 | **STRUCTURE** |
| T2' | v_2 \| residue (mod 81) | independence | χ²=69807125.4 | 0.00e+00 | (same verdict) |
| T3.5 | (mod 5, mod 27) joint | CRT product | χ²=17715329.2, df=68 | 0.00e+00 | **STRUCTURE** |
| T3.7 | (mod 7, mod 27) joint | CRT product | χ²=28680259.4, df=102 | 0.00e+00 | **STRUCTURE** |
| T3.11 | (mod 11, mod 27) joint | CRT product | χ²=46656114.8, df=170 | 0.00e+00 | **STRUCTURE** |
| T4 | return-time τ \| start mod 27 | starting residue is process variable | F=1.00, ANOVA | 4.58e-01 | **NULL** |
| T5 | visit-freq vs π_3 stationary | empirical → π_3 | Fisher χ²=24867.9, df=36 | 0.00e+00 | **STRUCTURE** |

## T1 — v_2 marginal vs geometric

| k | observed | expected (geom) | rel dev |
|---:|---:|---:|---:|
| 1 | 89,904,599 | 89,893,944 | +0.012% |
| 2 | 43,326,350 | 44,946,972 | -3.606% |
| 3 | 22,344,788 | 22,473,486 | -0.573% |
| 4 | 14,615,862 | 11,236,743 | +30.072% |
| 5 | 5,368,604 | 5,618,372 | -4.446% |
| 6 | 2,163,990 | 2,809,186 | -22.967% |
| 7 | 1,041,648 | 1,404,593 | -25.840% |
| 8 | 478,823 | 702,296 | -31.820% |
| 9 | 224,451 | 351,148 | -36.081% |
| 10 | 228,176 | 175,574 | +29.960% |

(See `result_local_invariant_sampling_v2_dist.csv` for k=1..29.)

## T2 — v_2 conditional on residue

Test of independence between v_2 levels (1..29, dropping empty bins) and residue class (coprime mod 27 / mod 81). Under the geometric null with CRT independence, v_2 should factor from any mod-3^k residue.

- mod 27: χ² = 27307163.54, df = 357, p = 0.0000e+00
- mod 81: χ² = 69807125.41, df = 1113, p = 0.0000e+00

## T3 — sibling primes

| prime p | χ² | df | p | verdict-this-pair |
|---:|---:|---:|---:|:---:|
| 5 | 17715329.23 | 68 | 0.0000e+00 | STRUCTURE |
| 7 | 28680259.44 | 102 | 0.0000e+00 | STRUCTURE |
| 11 | 46656114.80 | 170 | 0.0000e+00 | STRUCTURE |

## T4 — return-time conditional on starting residue

τ(n) ≈ 3.4713·log(n) + 1.4568  (Pearson R = 0.1438)
Lagarias prediction: 2/log(4/3) = 6.9521

ANOVA on residuals by starting mod 27 class: F = 0.9966, p = 4.5835e-01.

| start mod 27 | N_starts | mean(τ - τ_pred) | std |
|---:|---:|---:|---:|
| 1 | 185,185 | +0.0078 | 23.7518 |
| 2 | 185,185 | +0.0352 | 23.8425 |
| 4 | 185,185 | +0.0670 | 23.8562 |
| 5 | 185,186 | +0.0180 | 23.8976 |
| 7 | 185,186 | +0.0061 | 23.8995 |
| 8 | 185,185 | -0.0041 | 23.9273 |
| 10 | 185,185 | +0.0921 | 23.9439 |
| 11 | 185,185 | +0.0165 | 23.9702 |
| 13 | 185,185 | -0.0571 | 23.8372 |
| 14 | 185,185 | +0.0422 | 23.8866 |
| 16 | 185,185 | +0.0398 | 23.9405 |
| 17 | 185,185 | -0.0277 | 23.8415 |
| 19 | 185,185 | +0.0403 | 23.9184 |
| 20 | 185,185 | -0.1452 | 23.8058 |
| 22 | 185,185 | -0.0741 | 23.9152 |
| 23 | 185,185 | +0.0006 | 23.8596 |
| 25 | 185,185 | -0.0258 | 23.9089 |
| 26 | 185,185 | -0.0314 | 23.8359 |

## T5 — visit frequencies vs π_3

Fisher combined p over 18 starting residues: 0.0000e+00.

| start mod 27 | visits total | χ² vs π_3 | p |
|---:|---:|---:|---:|
| 1 | 9,989,667 | 409031.84 | 0.0000e+00 |
| 2 | 9,994,725 | 206306.03 | 0.0000e+00 |
| 4 | 10,000,614 | 312821.95 | 0.0000e+00 |
| 5 | 9,991,566 | 518064.30 | 0.0000e+00 |
| 7 | 9,989,370 | 919911.40 | 0.0000e+00 |
| 8 | 9,987,459 | 387946.71 | 0.0000e+00 |
| 10 | 10,005,280 | 251287.77 | 0.0000e+00 |
| 11 | 9,991,245 | 271634.75 | 0.0000e+00 |
| 13 | 9,977,618 | 281277.64 | 0.0000e+00 |
| 14 | 9,996,036 | 361270.96 | 0.0000e+00 |
| 16 | 9,995,583 | 353570.78 | 0.0000e+00 |
| 17 | 9,983,066 | 284175.35 | 0.0000e+00 |
| 19 | 9,995,671 | 309037.31 | 0.0000e+00 |
| 20 | 9,961,333 | 244917.90 | 0.0000e+00 |
| 22 | 9,974,508 | 253390.28 | 0.0000e+00 |
| 23 | 9,988,307 | 680843.68 | 0.0000e+00 |
| 25 | 9,983,426 | 760008.83 | 0.0000e+00 |
| 26 | 9,982,414 | 337821.51 | 0.0000e+00 |

## Interpretation

Of 5 tests: **1 NULL, 4 STRUCTURE** (Bonferroni α = 0.0100).

STRUCTURE was detected in some tests. See per-test sections above for effect sizes and per-residue/per-prime patterns. Note: with N ~ 10^7 starts and ~10^8 total steps, even small relative deviations from the null can register at very small p — interpret effect sizes (relative deviations and per-residue magnitudes), not just p-values, before claiming structural significance.

**Caveat on power.** This test has very high power to detect small deviations (large N). A STRUCTURE flag at p < 0.01 may represent a relative effect of 0.1–1% that doesn't carry useful arithmetic information. The effect-size tables above should be inspected before drawing structural conclusions.
