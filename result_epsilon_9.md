# Result: eps_9 = S_9 - 7/15

**Date:** 2026-05-05.  Float64 power iteration on K_9 (13122 states, M=13122), with scipy.sparse.linalg.eigs Arnoldi cross-check and FFT verification.

## Verdict

ACCELERATED DECAY — |eps_9/eps_8| < 0.40, faster than pre-spike. The k=7 spike was a one-off and post-spike trajectory is faster than the k=2..6 trend. Asymptotic rate may be > 1/2.

## Headline numbers

- eps_9 = `-7.5202571564e-06` (power iteration)
- eps_9 = `-7.5202571564e-06` (eigs cross-check; agreement 0.00e+00)
- eps_9 = `-7.5202571555e-06` (FFT cross-check; agreement 9.44e-16)
- |eps_9/eps_8| = **0.010087**
- |eps_8/eps_7| = 0.634380 (prior, for context)

## Ratio trajectory (k=1..9)

| k → k+1 | |eps_{k+1}/eps_k| |
|---|---|
| 1 → 2 | 0.047619 |
| 2 → 3 | 0.534659 |
| 3 → 4 | 0.481592 |
| 4 → 5 | 0.469668 |
| 5 → 6 | 0.432305 |
| 6 → 7 | 2.360360 |
| 7 → 8 | 0.634380 |
| 8 → 9 | 0.010087 |

## eps_k table (k=1..9)

| k | eps_k | source |
|---|---|---|
| 1 | +2.0000000000e-01 | exact rational (cached) |
| 2 | +9.5238095238e-03 | exact rational (cached) |
| 3 | -5.0919863259e-03 | exact rational (cached) |
| 4 | -2.4522582483e-03 | exact rational (cached) |
| 5 | -1.1517469151e-03 | exact rational (cached) |
| 6 | -4.9790566522e-04 | float64 power iter (486 states) |
| 7 | -1.1752368304e-03 | float64 power iter (1458 states) |
| 8 | -7.4554636729e-04 | float64 power iter (4374 states) |
| 9 | -7.5202571564e-06 | float64 power iter (13122 states) + eigs + FFT cross-checks |

## Computation diagnostics (k=9)

- K_9 build time (vectorized bincount): 0.75s
- Power iter: 10 iterations, residual = 3.11e-16, 0.31s
- eigs Arnoldi cross-check: 0.67s, leading eval = 1.000000000000
- |pi9_power - pi9_eigs|_1 = 4.68e-16, max-abs = 3.90e-18
- S_9 via X_9 - X_8 vs S_9 via FFT: diff = 9.44e-16
- pi_9 sum = 1.000000000000000
- K_9 row-sum max deviation from 1 = 2.22e-16

## Files

- `result_epsilon_9.py` — script
- `result_epsilon_9.csv` — eps_k and ratios for k=1..9
- `result_epsilon_9.md` — this writeup