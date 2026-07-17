# Diagnostic: eps_7 strategy + cost

**Strategy:** float64 power iteration (same as eps_6).

## Cost summary

- K_6 build + iter: 0.06s (7 iters)
- K_7 build + iter: 0.55s (8 iters)
- Convergence residual at k=7: 2.81e-16

## Precision check

- pi_7 sum check: 1.000000000000000 (deviation from 1: 0.00e+00)
- K_7 row-sum check: max |row sum - 1| = 2.22e-16
- eps_7 reported to 10 sig figs (float64 noise floor on values ~1e-4 is ~1e-19 absolute, ~1e-15 relative)

## What was NOT attempted

- Exact-rational stationary at k=7: previously killed in R77.7 after 7+ hours from denominator bloat. Not re-attempted per brief.

- High-precision (mpmath) double-check: not needed; float64 noise floor far below the precision the brief required (±5e-5).