# Low-v_2 Residue Admissibility — empirical horizon test

**Date:** 2026-05-05.  N = 3,333,332 odd-coprime-to-3 starts in [1, 10^7]; iterated up to 1000 Syracuse steps (most collapse to 1 in <100 steps).

**Admissibility horizon definition.**  k_max(m) = first step k where cumulative Σ_{i≤k} v_2(3n_i+1) ≥ m. Beyond k_max, the residue r_k mod 2^m is no longer determined by the initial residue mod 2^m alone — bits beyond 2^m are required to continue the trace deterministically.

Bit-budget arithmetic predicts: low-v_2 trajectories take MORE steps to consume m bits → LARGER k_max. The brief's hypothesis is the opposite (smaller k_max for low-v).

## Verdict

BIT-BUDGET CONFIRMED. Low-v_2 trajectories have SYSTEMATICALLY LARGER admissibility horizons than typical trajectories. The brief's hypothesis (smaller k_max for low-v) is REFUTED — the math predicts and the data confirms that lower per-step v_2 means more steps to exhaust m bits.

## Ensemble stats

| metric | value |
|---|---:|
| total v_2 mean | 107.52 |
| total v_2 median | 104 |
| trajectory length mean | 53.94 |
| trajectory length median | 51 |
| avg v_2 / step (mean) | 2.1006 |
| avg v_2 / step (median) | 2.0185 |

## Horizon reached fraction by m (full ensemble)

| m | frac reached | mean horizon | median horizon |
|---:|---:|---:|---:|
| 6 | 1.0000 | 3.50 | 3.0 |
| 8 | 1.0000 | 4.50 | 4.0 |
| 10 | 1.0000 | 5.50 | 6.0 |
| 12 | 1.0000 | 6.50 | 6.0 |
| 16 | 1.0000 | 8.50 | 8.0 |
| 20 | 1.0000 | 10.50 | 11.0 |

## Horizons stratified by bottom-q on (avg v_2 per step)

Definition: 'low-v_2' = trajectories with avg v_2 per step in the bottom-q quantile. Comparing mean k_max(m) for (low) vs (rest).

Ratio = mean horizon (low) / mean horizon (rest).  Ratio > 1 ↔ low-v has LARGER horizons (bit-budget prediction).  Ratio < 1 ↔ low-v has SMALLER horizons (brief's hypothesis).

| q | m | n_low | mean h (low) | mean h (rest) | median h (low) | median h (rest) | ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.001 | 6 | 3,347 | 4.47 | 3.50 | 5.0 | 3.0 | 1.2763 |
| 0.001 | 8 | 3,347 | 5.83 | 4.50 | 6.0 | 4.0 | 1.2963 |
| 0.001 | 10 | 3,347 | 7.19 | 5.50 | 7.0 | 5.0 | 1.3080 |
| 0.001 | 12 | 3,347 | 8.54 | 6.50 | 9.0 | 6.0 | 1.3139 |
| 0.001 | 16 | 3,347 | 11.21 | 8.50 | 11.0 | 8.0 | 1.3198 |
| 0.001 | 20 | 3,347 | 13.87 | 10.50 | 14.0 | 10.0 | 1.3217 |
| 0.010 | 6 | 33,719 | 4.21 | 3.49 | 4.0 | 3.0 | 1.2066 |
| 0.010 | 8 | 33,719 | 5.49 | 4.49 | 6.0 | 4.0 | 1.2237 |
| 0.010 | 10 | 33,719 | 6.77 | 5.49 | 7.0 | 5.0 | 1.2340 |
| 0.010 | 12 | 33,719 | 8.04 | 6.48 | 8.0 | 6.0 | 1.2399 |
| 0.010 | 16 | 33,719 | 10.57 | 8.48 | 11.0 | 8.0 | 1.2470 |
| 0.010 | 20 | 33,719 | 13.09 | 10.47 | 13.0 | 10.0 | 1.2494 |
| 0.050 | 6 | 170,875 | 4.05 | 3.47 | 4.0 | 3.0 | 1.1665 |
| 0.050 | 8 | 170,875 | 5.27 | 4.46 | 5.0 | 4.0 | 1.1811 |
| 0.050 | 10 | 170,875 | 6.48 | 5.45 | 6.0 | 5.0 | 1.1903 |
| 0.050 | 12 | 170,875 | 7.70 | 6.44 | 8.0 | 6.0 | 1.1965 |
| 0.050 | 16 | 170,875 | 10.13 | 8.41 | 10.0 | 8.0 | 1.2041 |
| 0.050 | 20 | 170,875 | 12.55 | 10.39 | 13.0 | 10.0 | 1.2082 |

**Cells with ratio > 1 (low has larger horizon):** 18/18
**Cells with ratio < 1 (low has smaller horizon):** 0/18

## Files

- `result_low_v2_residue_admissibility.py` — script
- `result_low_v2_residue_admissibility.csv` — per-(q, m) horizon stats
- `result_low_v2_residue_admissibility.md` — this writeup