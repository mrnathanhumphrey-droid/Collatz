# Diagnostic: eps_6 strategy + cost

**Strategy used:** A (float64 power iteration on K_6).

## Why Strategy A

- K_6 has **486 states**, not 1458 as the brief stated. (The 1458-state chain is k=7, where R77.7 died.) k=6 is well within tractable range.
- Float64 power iteration converges at rate equal to the second eigenvalue of K_6 (empirically near 1/2 per R66), reaches 1e-14 residual in ~50 iterations, sub-second total.
- Required precision (eps_6 to ±0.0005) is far below float64's ~15 sig digits.
- Strategy B (R75 recursion) would still need pi_5 and the off-diagonal mass at k=5->6, which scales similarly and offers no speedup.

## Computational cost

- K_6 build: ~milliseconds
- Power iteration: 7 iters to residual 2.59e-16, 0.00s wall time
- Total: <1 second

## Precision achieved

- pi_6 sum check: 1.000000000000000 (deviation from 1: 0.00e+00)
- K_6 row-sum check: max |row sum - 1| = 1.11e-16
- eps_6 reported to 10 significant digits (float64 ~15 digits available)

## What was NOT done

- Strategy B (R75 recursion via Theorem 75.2): not necessary; Strategy A produced sufficient precision in <1 second. Strategy B would only be advantageous if it yielded exact rationals; for the two-mode prediction test 10-digit float is already definitive.
- Exact Fraction stationary at k=6 (would give exact S_6 for algebraic-form check): NOT FIRED in this run; could be added as a separate pass if the structural check (step 5 of brief) becomes load-bearing. Rough cost estimate: ~2-30 minutes given k=5 took 5 sec at 162 states (cubic scaling + denominator bloat).