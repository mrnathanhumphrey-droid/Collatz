# Result: eps_7 = S_7 - 7/15 (asymptotic rate)

**Date:** 2026-05-05.  Float64 power iteration on K_7 (1458 states).

## Verdict: **RATIO REVERSED UP — non-monotone trajectory; structurally informative**

- eps_7 = `-1.1752368304e-03`
- |eps_7/eps_6| = **2.360360**

**Cross-validated** (`result_epsilon_7_verify.py`): float64 power iteration and scipy.eigs left-dominant eigenvector both give the same eps_7 to 1e-15. At k=5 the float64 path matches the cached exact rational to 1e-15. The non-monotone result is real, not numerical noise.

## Structural reading: S_k is non-monotone

|ε_k| has a local minimum at k=6 and bounces back at k=7:

| k | S_k | eps_k | \|eps_k\| |
|---|---|---|---|
| 5 | 0.46551 | -0.00115 | 1.15e-3 |
| 6 | **0.46617** | -0.00050 | **4.98e-4** ← local min |
| 7 | 0.46549 | -0.00118 | 1.18e-3 ← bounced back |

This **kills the "|eps_n|·2^n stable near 0.04" rate-1/2 envelope claim** that was load-bearing through k=5. Updated envelope:

| n | \|eps_n\|·2^n |
|---|---|
| 2 | 0.0381 |
| 3 | 0.0407 |
| 4 | 0.0392 |
| 5 | 0.0368 |
| 6 | 0.0319 |
| 7 | **0.1504** ← 4× envelope |

The "stable near 0.04" was the visible portion of a longer-period oscillation. Asymptotic-rate question is now wide open. Plausible candidate forms (no refit per brief): damped oscillation with complex-conjugate rates ρ·e^±iθ, longer-period structure, or non-elementary asymptotic shape.

## Ratio trajectory (k=1..7)

| transition | |eps_{k+1}/eps_k| |
|---|---|
| 1 -> 2 | 0.047619 |
| 2 -> 3 | 0.534659 |
| 3 -> 4 | 0.481592 |
| 4 -> 5 | 0.469668 |
| 5 -> 6 | 0.432305 |
| 6 -> 7 | 2.360360 |

## eps_k table (k=1..7)

| k | eps_k | source |
|---|---|---|
| 1 | +2.0000000000e-01 | exact rational (cached) |
| 2 | +9.5238095238e-03 | exact rational (cached) |
| 3 | -5.0919863259e-03 | exact rational (cached) |
| 4 | -2.4522582483e-03 | exact rational (cached) |
| 5 | -1.1517469151e-03 | exact rational (cached) |
| 6 | -4.9790566522e-04 | float64 power iter (487 states) |
| 7 | -1.1752368304e-03 | float64 power iter (1458 states) |

## Computation diagnostics

| step | wall time | iterations | residual |
|---|---|---|---|
| K_6 build | 0.06s | — | — |
| K_6 power iter | 0.00s | 7 | 2.59e-16 |
| K_7 build | 0.55s | — | — |
| K_7 power iter | 0.00s | 8 | 2.81e-16 |

Total wall: 0.61s (brief estimated 1-4 hours; actual sub-minute on this hardware).

## Files

- `result_epsilon_7.py` — script
- `result_epsilon_7.csv` — eps_k and ratios for k=1..7
- `result_epsilon_7.md` — this writeup
- `result_epsilon_7_diagnostic.md` — strategy + cost