# Result: eps_6 = S_6 - 7/15 (two-mode prediction test)

**Date:** 2026-05-05.

## Verdict: **TWO-MODE INCOMPLETE / amplitudes wrong (need refit)**

- Computed eps_6 (float64 power iteration on K_6, residual 3e-16 after 7 iterations) = `-4.9790566522e-04`.
- Two-mode prediction (A,B = -0.4918, +1.3268): -5.864348e-03.
- Pure (1/2)^k prediction: -7.684375e-03.
- Pure (1/3)^k prediction: +1.820027e-03.

## eps_k table (k=1..6)

| k | eps_k | predicted (k=1..5 fit) | residual |
|---|---|---|---|
| 1 | +2.000000e-01 | +1.963667e-01 | +3.633333e-03 |
| 2 | +9.523810e-03 | +2.447222e-02 | -1.494841e-02 |
| 3 | -5.091986e-03 | -1.233426e-02 | +7.242273e-03 |
| 4 | -2.452258e-03 | -1.435725e-02 | +1.190499e-02 |
| 5 | -1.151747e-03 | -9.908668e-03 | +8.756921e-03 |
| 6 | -4.979057e-04 | -5.864348e-03 | +5.366442e-03 |

## Two-mode refit comparison

| fit window | A | B | SS_res |
|---|---|---|---|
| k=1..5 (recomputed) | -0.491778 | +1.326776 | 5.0752e-04 |
| k=1..6              | -0.486345 | +1.318142 | 5.3592e-04 |

Amplitude drift: |ΔA|/|A| = 0.0110, |ΔB|/|B| = 0.0065

Both amplitudes stable to <5% — two-mode characterization is robust under inclusion of k=6.

## |eps_{k+1}/eps_k| ratios

| transition | |ratio| |
|---|---|
| 1 -> 2 | 0.047619 |
| 2 -> 3 | 0.534659 |
| 3 -> 4 | 0.481592 |
| 4 -> 5 | 0.469668 |
| 5 -> 6 | 0.432305 |

|ratio| at k=4->5 was 0.4697; new k=5->6 ratio = 0.432305.

## Files

- `result_epsilon_6.py` — script (Strategy A, float64 power iteration)
- `result_epsilon_6.csv` — eps_k values for k=1..6, predicted vs actual
- `result_epsilon_6.md` — this writeup
- `result_epsilon_6_diagnostic.md` — strategy + cost diagnostics