# Result: eps_11 = S_11 - 7/15

**Date:** 2026-05-05.  Float64 matrix-free power iteration on K_11 (118,098 states, M=118,098). FFT cross-check.

## Verdict

EXPLORE — |eps_11/eps_10| = 2.0839, sign flip = False. Doesn't fit known model regimes cleanly.

## Headline numbers

- eps_11 = `+1.5019670121e-03` (power iter)
- eps_11 = `+1.5019670121e-03` (FFT cross-check; agreement 2.00e-15)
- |eps_11/eps_10| = **2.083892**
- sign(eps_11) = +

## Two-mode model check

Model: eps_k ≈ scale · ρ^k · cos(k·θ + φ),  ρ = 0.984, θ = 0.68 rad (period ≈ 9.2)

- eps_11 predicted (scale-calibrated to eps_10): `+5.5147e-04`
- eps_11 actual: `+1.5020e-03`
- actual / predicted: `2.7236`

## Ratio trajectory (k=1..11)

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
| 10 → 11 | 2.083892 |

## eps_k table (k=1..11)

| k | eps_k | sign | source |
|---|---|:---:|---|
| 1 | +2.0000000000e-01 | + | exact rational |
| 2 | +9.5238095238e-03 | + | exact rational |
| 3 | -5.0919863259e-03 | - | exact rational |
| 4 | -2.4522582483e-03 | - | exact rational |
| 5 | -1.1517469151e-03 | - | exact rational |
| 6 | -4.9790566522e-04 | - | float64 power iter (cached) |
| 7 | -1.1752368304e-03 | - | float64 power iter (cached) |
| 8 | -7.4554636729e-04 | - | float64 power iter (cached) |
| 9 | -7.5202571564e-06 | - | float64 power iter (cached) |
| 10 | +7.2075091711e-04 | + | float64 power iter (cached) |
| 11 | +1.5019670121e-03 | + | float64 matrix-free + FFT |

## Computation diagnostics (k=11)

- States: 118,098, M = 118,098
- Init MatVecK: 0.03s
- Power iter: 12 iters, residual = 6.99e-16, 1443.04s total (120.3s per matvec)
- FFT cross-check vs X-formula: diff = 2.00e-15
- pi_11 sum = 1.000000000000000

## Files

- `result_epsilon_11.py` — script
- `result_epsilon_11.csv` — eps_k and ratios for k=1..11
- `result_epsilon_11.md` — this writeup