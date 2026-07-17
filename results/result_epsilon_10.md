# Result: eps_10 = S_10 - 7/15

**Date:** 2026-05-05.  Float64 matrix-free power iteration on K_10 (39366 states, M=39366). Dense storage (12.4 GB) avoided via chunked bincount per orbit-power. scipy.eigs and FFT cross-checks.

## Verdict

OSCILLATION CONFIRMED — |eps_10/eps_9| = 95.84 >> 1. The k=9 collapse was a near-zero node, and eps_10 rebounds to a magnitude consistent with the next half-cycle peak. The k-space trajectory of eps_k carries a non-trivial oscillating component (complex eigenpair in K).

## Headline numbers

- eps_10 = `+7.2075091711e-04` (power iter)
- eps_10 = `+7.2075091712e-04` (eigs cross-check; agreement 3.55e-15)
- eps_10 = `+7.2075091711e-04` (FFT cross-check; agreement 6.11e-16)
- |eps_10/eps_9| = **95.841259**
- |eps_9/eps_8| = 0.010087 (prior)

## Ratio trajectory (k=1..10)

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
| 9 → 10 | 95.841259 |

## eps_k table (k=1..10)

| k | eps_k | source |
|---|---|---|
| 1 | +2.0000000000e-01 | exact rational (cached) |
| 2 | +9.5238095238e-03 | exact rational (cached) |
| 3 | -5.0919863259e-03 | exact rational (cached) |
| 4 | -2.4522582483e-03 | exact rational (cached) |
| 5 | -1.1517469151e-03 | exact rational (cached) |
| 6 | -4.9790566522e-04 | float64 dense power iter |
| 7 | -1.1752368304e-03 | float64 dense power iter |
| 8 | -7.4554636729e-04 | float64 dense power iter |
| 9 | -7.5202571564e-06 | float64 dense power iter |
| 10 | +7.2075091711e-04 | float64 matrix-free + eigs + FFT |

## Computation diagnostics (k=10)

- States: 39,366, M = 39,366
- Init MatVecK: 0.01s
- Power iter: 11 iterations, residual = 7.65e-16, 154.31s total
- eigs Arnoldi cross-check: 290.28s, leading eval = 1.000000000000
- |pi10_power - pi10_eigs|_1 = 8.62e-16, max-abs = 1.73e-18
- S_10 via X_10 - X_9 vs S_10 via FFT: diff = 6.11e-16
- pi_10 sum = 1.000000000000000

## Files

- `result_epsilon_10.py` — script (with matrix-free K_10)
- `result_epsilon_10.csv` — eps_k and ratios for k=1..10
- `result_epsilon_10.md` — this writeup