# Result: eps_8 = S_8 - 7/15

**Date:** 2026-05-05.  Float64 power iteration on K_8 (4374 states, M=4374), with scipy.sparse.linalg.eigs Arnoldi cross-check.

## Verdict

AMBIGUOUS — 0.5 < |eps_8/eps_7| < 0.7. Neither clean decay nor sustained bouncing. eps_9 (k=9 chain, ~13k states, ~3-5 min compute) would disambiguate.

## Headline numbers

- eps_8 = `-7.4554636729e-04` (power iteration)
- eps_8 = `-7.4554636729e-04` (scipy.eigs cross-check; agreement 1.78e-15)
- eps_8 = `-7.4554636729e-04` (FFT cross-check; agreement 3.33e-16)
- |eps_8/eps_7| = **0.634380**
- |eps_7/eps_6| = 2.360360 (prior, for context)

## Ratio trajectory (k=1..8)

| k → k+1 | |eps_{k+1}/eps_k| |
|---|---|
| 1 → 2 | 0.047619 |
| 2 → 3 | 0.534659 |
| 3 → 4 | 0.481592 |
| 4 → 5 | 0.469668 |
| 5 → 6 | 0.432305 |
| 6 → 7 | 2.360360 |
| 7 → 8 | 0.634380 |

## eps_k table (k=1..8)

| k | eps_k | source |
|---|---|---|
| 1 | +2.0000000000e-01 | exact rational (cached) |
| 2 | +9.5238095238e-03 | exact rational (cached) |
| 3 | -5.0919863259e-03 | exact rational (cached) |
| 4 | -2.4522582483e-03 | exact rational (cached) |
| 5 | -1.1517469151e-03 | exact rational (cached) |
| 6 | -4.9790566522e-04 | float64 power iter, K_6 (729 states) |
| 7 | -1.1752368304e-03 | float64 power iter, K_7 (1458 states) |
| 8 | -7.4554636729e-04 | float64 power iter, K_8 (4374 states) + eigs cross-check |

## Computation diagnostics

| step | wall time | iterations | residual |
|---|---|---|---|
| K_6 build | 0.06s | — | — |
| K_6 power iter | 0.00s | 7 | 2.59e-16 |
| K_7 build | 0.58s | — | — |
| K_7 power iter | 0.00s | 8 | 2.81e-16 |
| K_8 build | 7.08s | — | — |
| K_8 power iter | 0.02s | 9 | 2.84e-16 |
| K_8 eigs cross-check | 0.05s | — | — |

**Cross-checks at k=8:**
- power iter pi vs eigs vector: |·|_1 diff = 5.66e-16, max-abs = 3.90e-18
- S_8 via X_8 - X_7 vs S_8 via FFT: diff = 3.33e-16
- pi_8 sum = 1.000000000000000
- K_8 row-sum max deviation from 1 = 2.22e-16

## Files

- `result_epsilon_8.py` — script
- `result_epsilon_8.csv` — eps_k and ratios for k=1..8
- `result_epsilon_8.md` — this writeup