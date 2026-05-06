# Result: eps_12 = S_12 - 7/15

**Date:** 2026-05-05.  Float64 matrix-free power iteration on K_12 (354,294 states, M=354,294). FFT cross-check.

## Verdict

RECURRENCE CONFIRMED: eps_12 = +2.2747e-03 lies within pre-registered band [+1.50e-3, +2.5e-3]. Order-3 recurrence model strengthened.

## Headline numbers

- eps_12 = `+2.2747137206e-03` (power iter, X-formula)
- eps_12 = `+2.2747137206e-03` (FFT cross-check; agreement 1.33e-15)

- eps_11 = `+1.5019670121e-03`
- |eps_12 / eps_11| = 1.514490
- sign(eps_12) = +  (eps_11 was +)

## Pre-registered prediction

From order-3 linear recurrence fit on eps_2..eps_11:
- Predicted band: `[+1.50e-3, +2.5e-3]`
- Measured: `+2.2747137206e-03`
- Outcome: in band

## eps_k table (k=6..k)

| k | eps_k | sign | source |
|---|---|:---:|---|
| 6 | -4.9790566522e-04 | - | float64 power iter (cached) |
| 7 | -1.1752368304e-03 | - | float64 power iter (cached) |
| 8 | -7.4554636729e-04 | - | float64 power iter (cached) |
| 9 | -7.5202571564e-06 | - | float64 power iter (cached) |
| 10 | +7.2075091711e-04 | + | float64 power iter (cached) |
| 11 | +1.5019670121e-03 | + | float64 power iter (cached) |
| 12 | +2.2747137206e-03 | + | this run |

## Compute diagnostics

- States n = 354,294
- M = 354,294
- Chunk size = 256
- Power iter: 13 iters, final residual = 7.40e-16
- Total iter time: 19003.2s (5.28 hours)
- Per-matvec: 1461.8s
- FFT cross-check vs X-formula: diff = 1.33e-15

## Method notes

Method: matrix-free power iteration. K_k is *dense* per row (M nonzeros / row, M ≈ n), so a sparse csr representation would not save memory; instead we exploit the multiplicative structure: each row j of K is the histogram of (2^{-(v+1)} / Z_v) summed over v=0..M-1 onto state idx((q*r_j+1)*2^{-v-1} mod N). bincount over chunked v keeps peak memory bounded by chunk*n*8 bytes (here 0.73 GB).

Aitken acceleration was wired in (componentwise Δ² applied when residual ratio is in [0.7, 0.99]) but typically not triggered — power iter on K_k converges sharply once a subspace alignment is reached, so the linear-convergence regime is brief.