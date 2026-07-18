# Result: mode amplitude decomposition of pi_k onto K_k eigenvectors

**Date:** 2026-05-05.

## Setup

| k | n_states | pi power-iter steps | top-20 eigvalue range |
|---|---|---|---|
| 5 | 162 | 6 | 1.0000e+00 ↓ 2.4317e-06 |
| 6 | 486 | 7 | 1.0000e+00 ↓ 2.4727e-04 |
| 7 | 1458 | 8 | 1.0000e+00 ↓ 1.0464e-03 |

## Reconstruction quality (top-20 spans pi_k how well)

| k | ||pi - Σ a_i v_i||_inf | rel to ||pi||_inf |
|---|---|---|
| 5 | 3.9003e-02 | 0.8634 |
| 6 | 2.0613e-02 | 0.9092 |
| 7 | 1.0676e-02 | 0.9396 |

## Per-k top-5 amplitudes |a_i|² and variance fractions

### k = 5

| rank | |λ| | arg(λ) | |a_i| | |a_i|² | variance frac |
|---|---|---|---|---|---|
| 1 | 1.0000e+00 | +0.0000 | 7.8567e-02 | 6.1728e-03 | 1.0000 |
| 2 | 3.4282e-04 | +0.0000 | 3.5996e-17 | 1.2957e-33 | 0.0000 |
| 3 | 3.4278e-04 | +1.2568 | 1.9950e-17 | 3.9799e-34 | 0.0000 |
| 4 | 3.4278e-04 | -1.2568 | 1.9950e-17 | 3.9799e-34 | 0.0000 |
| 5 | 3.4272e-04 | +2.5134 | 3.5997e-17 | 1.2958e-33 | 0.0000 |

### k = 6

| rank | |λ| | arg(λ) | |a_i| | |a_i|² | variance frac |
|---|---|---|---|---|---|
| 1 | 1.0000e+00 | +0.0000 | 4.5361e-02 | 2.0576e-03 | 1.0000 |
| 2 | 1.2267e-03 | +0.5241 | 1.0083e-17 | 1.0167e-34 | 0.0000 |
| 3 | 1.2267e-03 | -0.5241 | 1.0083e-17 | 1.0167e-34 | 0.0000 |
| 4 | 1.2256e-03 | +1.5718 | 2.2009e-17 | 4.8441e-34 | 0.0000 |
| 5 | 1.2256e-03 | -1.5718 | 2.2009e-17 | 4.8441e-34 | 0.0000 |

### k = 7

| rank | |λ| | arg(λ) | |a_i| | |a_i|² | variance frac |
|---|---|---|---|---|---|
| 1 | 1.0000e+00 | +0.0000 | 2.6189e-02 | 6.8587e-04 | 1.0000 |
| 2 | 2.9589e-03 | +0.0000 | 5.0415e-18 | 2.5417e-35 | 0.0000 |
| 3 | 2.9577e-03 | +0.8985 | 1.2036e-17 | 1.4487e-34 | 0.0000 |
| 4 | 2.9577e-03 | -0.8985 | 1.2036e-17 | 1.4487e-34 | 0.0000 |
| 5 | 2.9549e-03 | +1.7963 | 8.7826e-18 | 7.7133e-35 | 0.0000 |

## Pre-registered questions

### Q1: Mode-crossing

- k=5 dominant non-trivial: rank 5, |a|² = 1.2958e-33
- k=6 dominant non-trivial: rank 4, |a|² = 4.8441e-34
- k=7 dominant non-trivial: rank 3, |a|² = 1.4487e-34
- Mode-crossing observed: **YES**

### Q2: Amplitude growth/decay with k

Of 19 non-trivial lift-tracked modes:
- growing (a_5 < a_6 < a_7): **0**
- decaying (a_5 > a_6 > a_7): **0**
- non-monotone: **19**

### Q3: Connection to ε_k

Reference: |ε_5|, |ε_6|, |ε_7| = 1.1517e-03, 4.9791e-04, 1.1752e-03.
|ε_k| has a local minimum at k=6 and bounces back at k=7 (ratio 0.43, 2.36).

Modes with same qualitative non-monotone shape (a_6 < a_5, a_7 > a_6): 17 found:
- rank 2 at k=5: |a|² = 1.2957e-33, 1.0167e-34, 1.4487e-34
- rank 5 at k=5: |a|² = 1.2958e-33, 2.0764e-35, 8.4463e-35
- rank 6 at k=5: |a|² = 1.2958e-33, 2.0764e-35, 8.4463e-35
- rank 7 at k=5: |a|² = 1.9866e-34, 2.0764e-35, 8.4463e-35
- rank 8 at k=5: |a|² = 1.9866e-34, 2.0764e-35, 8.4463e-35
- rank 9 at k=5: |a|² = 7.5301e-35, 2.0764e-35, 8.4463e-35
- rank 10 at k=5: |a|² = 7.5301e-35, 2.0764e-35, 8.4463e-35
- rank 11 at k=5: |a|² = 3.9798e-34, 2.7860e-35, 2.9904e-35
- rank 12 at k=5: |a|² = 9.7639e-34, 2.7860e-35, 2.9904e-35
- rank 13 at k=5: |a|² = 9.7639e-34, 2.7860e-35, 2.9904e-35

## Files

- `mode_amplitudes_probe.py` — script
- `pi_k{5,6,7}.npy` — cached stationary distributions
- `mode_amplitudes_k{5,6,7}.csv` — per-k amplitude tables
- `mode_amplitudes_comparison.csv` — cross-k tracked via lift
- `mode_amplitudes_findings.md` — this writeup