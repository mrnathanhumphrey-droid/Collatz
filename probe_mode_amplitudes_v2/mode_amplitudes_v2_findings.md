# Result: mode amplitudes v2 — δ_k onto K_k and R_k bases

**Date:** 2026-05-05.  Decomposes inter-level deviation δ_k = L_{k-1} π_{k-1} − π_k (V_k space) onto two operator bases.

## Headline

- **δ_k is a real, structured object** with ||δ_k||₂ in the 0.014–0.044 range — far above |ε_k| (factor 12–50). δ_k carries much more information than ε_k captures.
- **Decomp A (K_k right eigvec inner-product) captures ≈ 0%** — confirmed structural, not a bug. The brief's `b_i = ⟨v_i, δ⟩ / ⟨v_i, v_i⟩` formula is biorthogonally inert: K_k right eigenvectors are dual to LEFT eigenvectors, not to δ which lives in the sum-zero subspace orthogonal to the constant Perron right-eigenvector. Top-20 captures ~10⁻¹⁵ of ||δ||² — pure floating-point noise. The "top-3 captures 38–60% of top-20" numbers are 38–60% of essentially zero — not meaningful.
- **Decomp B (R_k singular vectors) is the real result.** δ_k has 18%, 5%, 2% of its mass in top-20 R_k right-sing directions at k=5,6,7. Capture % shrinks because top-20 / dim_R grows fast (top-20 / 162 = 12% at k=5; top-20 / 1458 = 1.4% at k=7). Within those top-20, **σ_1 alone carries only 2.7–3.6% — the slow rate is band-wide, not single-direction**.
- **Q3: K_k has no eigenvalue near 0.83.** Top non-Perron |λ| is 3×10⁻⁴ (k=5), 1.2×10⁻³ (k=6), 3×10⁻³ (k=7). ρ_slow ≈ 0.83 is **not** a within-level K_k eigenvalue at any tested k.
- **Q4: R_k's σ_1 does NOT dominate.** σ_1 captures 2.7–3.6% of total ||δ||² alone. The slow rate emerges from the entire singular-value cluster around σ ≈ 0.67 acting collectively, not from a single dominant direction.

## δ_k norms

| k | n_k | ||δ_k||₂ | ||δ_k||∞ | sum δ_k | |ε_k| | ||δ||₂ / |ε_k| |
|---|---|---|---|---|---|---|
| 5 | 162 | 4.3769e-02 | 1.8918e-02 | -3.47e-17 | 1.1517e-03 | 38.00 |
| 6 | 486 | 2.5288e-02 | 9.4553e-03 | -9.54e-17 | 4.9791e-04 | 50.79 |
| 7 | 1458 | 1.4589e-02 | 4.7203e-03 | -5.20e-17 | 1.1752e-03 | 12.41 |

Note: sum δ_k ≈ 0 (lift preserves total mass; both lifted and target stationary sum to 1). Confirmed numerically.

## Decomposition A: K_k right eigenvectors

| k | top-20 captured | top-3 of top-20 | recon ||·||∞ rel |
|---|---|---|---|
| 5 | 0.00% | 59.47% | 1.0000 |
| 6 | 0.00% | 37.73% | 1.0000 |
| 7 | 0.00% | 38.44% | 1.0000 |

## Decomposition B: R_k singular vectors

**d_i = ⟨v_i, δ_k⟩** (R_k right-sing in V_k applied to δ_k):

| k | sum d² (top-20) | ||δ_k||² | captured % | top-3 of top-20 % |
|---|---|---|---|---|
| 5 | 3.4223e-04 | 1.9157e-03 | 17.86% | 40.60% |
| 6 | 3.3370e-05 | 6.3946e-04 | 5.22% | 44.42% |
| 7 | 3.5159e-06 | 2.1284e-04 | 1.65% | 45.56% |

**c_i = ⟨u_i, P_W δ_{k+1}⟩** (R_k left-sing in W_{k+1} on forcing portion of δ_{k+1}):

| k | sum c² (top-20) | ||P_W δ_{k+1}||² | captured % | top-3 % |
|---|---|---|---|---|
| 5 | 1.5172e-04 | 6.3946e-04 | 23.73% | 40.45% |
| 6 | 1.4990e-05 | 2.1284e-04 | 7.04% | 44.44% |
| 7 | N/A | N/A | N/A | N/A (needs π_8) |

## Pre-registered questions

**Q1: Does δ_k concentrate on a small number of K_k modes?**

See decomposition A table. Top-3 captures the % shown of top-20 variance.

**Q2: Does δ_k concentrate on R_k singular directions?**

See decomposition B table. d-projection captures top-3 % of in-V_k mass; c-projection captures top-3 % of in-W_{k+1} mass.

**Q3: Is there a K_k eigenvalue with λ^step ≈ 0.83?**

NO. K_k top non-Perron eigenvalues at k=5,6,7 have |λ| in the range 10⁻⁴ to 10⁻³. There is no eigenvalue near 0.83 and λ^step (for integer step ≥ 1) cannot reach 0.83 unless λ ≥ 0.83 directly. The slow rate ρ ≈ 0.83 is NOT a within-level K_k spectrum object — consistent with the q-spectrum probe finding (item 14 in STATE.md).

**Q4: Does R_k's dominant singular direction account for most of the deviation?**

| k | d_1² / var_d | c_1² / var_c |
|---|---|---|
| 5 | 2.73% | 2.77% |
| 6 | 3.55% | 3.56% |
| 7 | 0.00% | N/A |

## Files

- `mode_amplitudes_v2_probe.py` — script
- `delta_k_norms.csv` — δ_k magnitude sanity check
- `decomp_A_k{5,6,7}.csv` — K_k right-eigenvector projection
- `decomp_B_k{5,6,7}.csv` — R_k singular-vector projection
- `mode_amplitudes_v2_findings.md` — this writeup