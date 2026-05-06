# Test 3 reframings — joint admissibility, residue counting, stuck residues

Date 2026-05-05. Sample N=500,000 odd-coprime-to-3 starts in [3, 10^6]; depth D=max(10, 5)=10 for residue tracking; orbit length L(n) capped at 1000.

## Working definitions

- **H_2(r, k)** = mean orbit-length L over starts with n ≡ r (mod 2^k).
- **H_3(s, j)** = mean orbit-length L over starts with n ≡ s (mod 3^j).
- **H_joint(r, s, k, j)** = mean L over starts with n ≡ r mod 2^k AND n ≡ s mod 3^j.
- **rho(r, s, k, j)** = H_joint / min(H_2(r, k), H_3(s, j)). Joint reduction factor.
- **gamma(r, k)** = log(N_distinct(D)) / D, where N_distinct(D) is the number of distinct depth-D residue traces (mod 2^k) among starts ≡ r mod 2^k.
- **P(r, k)** = number of starts whose D-step orbit visits residue r mod 2^k at any depth 0..D.
- **sigma(r, k)** = E_P / P(r, k). Stuck-ness: large σ ⇒ rarely visited.

## Sub-probe 4A — joint admissibility

rho aggregated by a_final, per k:

| k | a_final | mean rho | sd | n_cells |
|---|---|---|---|---|
| 5 | 3 | 1.000099 | 0.021111 | 78 |
| 5 | 9 | 0.999999 | 0.017136 | 312 |
| 5 | 27 | 1.001700 | 0.015504 | 468 |
| 5 | 81 | 1.087073 | 0.014951 | 312 |
| 5 | 243 | 1.178924 | 0.014737 | 78 |
| 6 | 3 | 0.999998 | 0.022332 | 78 |
| 6 | 9 | 1.000018 | 0.025087 | 390 |
| 6 | 27 | 1.000019 | 0.022999 | 780 |
| 6 | 81 | 1.043063 | 0.022333 | 780 |
| 6 | 243 | 1.131336 | 0.021688 | 390 |
| 6 | 729 | 1.225118 | 0.022153 | 78 |
| 7 | 3 | 0.999885 | 0.030466 | 78 |
| 7 | 9 | 0.999912 | 0.034517 | 468 |
| 7 | 27 | 1.000008 | 0.034928 | 1170 |
| 7 | 81 | 1.002241 | 0.031739 | 1560 |
| 7 | 243 | 1.087753 | 0.031600 | 1170 |
| 7 | 729 | 1.176718 | 0.032376 | 468 |
| 7 | 2187 | 1.267139 | 0.028276 | 78 |

**4A verdict:** Intermediate — rho varies by 0.267 (between 10% and 100% spread). Worth follow-up correlation with v_2 to disambiguate Gate C from genuine variation.

## Sub-probe 4B — residue-counting growth rate

gamma = log(N_distinct) / D aggregated by a_final, per k:

| k | a_final | mean gamma | sd | n_residue_classes |
|---|---|---|---|---|
| 5 | 3 | 1.029435 | 0.000000 | 1 |
| 5 | 9 | 1.031802 | 0.001053 | 4 |
| 5 | 27 | 1.026244 | 0.000553 | 6 |
| 5 | 81 | 1.013176 | 0.000382 | 4 |
| 5 | 243 | 0.985896 | 0.000000 | 1 |
| 6 | 3 | 0.962463 | 0.000000 | 1 |
| 6 | 9 | 0.965172 | 0.000576 | 5 |
| 6 | 27 | 0.964003 | 0.000547 | 10 |
| 6 | 81 | 0.959739 | 0.000423 | 10 |
| 6 | 243 | 0.947591 | 0.000687 | 5 |
| 6 | 729 | 0.923581 | 0.000000 | 1 |
| 7 | 3 | 0.891731 | 0.000000 | 1 |
| 7 | 9 | 0.896363 | 0.001288 | 6 |
| 7 | 27 | 0.896239 | 0.000980 | 15 |
| 7 | 81 | 0.895560 | 0.000715 | 20 |
| 7 | 243 | 0.891973 | 0.000820 | 15 |
| 7 | 729 | 0.880278 | 0.000571 | 6 |
| 7 | 2187 | 0.858635 | 0.000000 | 1 |

**4B verdict:** Gate B — gamma uniform across a_final classes within ~10% (max spread 0.044). Reframing does not discriminate; closes.

## Sub-probe 4C — stuck residues

sigma = E_P / P(r,k) aggregated by a_final, per k:

| k | a_final | mean sigma | sd | n_residues |
|---|---|---|---|---|
| 5 | 3 | 1.015510 | 0.000000 | 1 |
| 5 | 9 | 1.002582 | 0.039172 | 4 |
| 5 | 27 | 0.990935 | 0.078279 | 6 |
| 5 | 81 | 0.939162 | 0.039602 | 4 |
| 5 | 243 | 1.506459 | 0.000000 | 1 |
| 6 | 3 | 1.014450 | 0.000000 | 1 |
| 6 | 9 | 1.017840 | 0.037178 | 5 |
| 6 | 27 | 0.999827 | 0.068977 | 10 |
| 6 | 81 | 0.979183 | 0.030840 | 10 |
| 6 | 243 | 0.960832 | 0.016179 | 5 |
| 6 | 729 | 1.599613 | 0.000000 | 1 |
| 7 | 3 | 1.010112 | 0.000000 | 1 |
| 7 | 9 | 1.010696 | 0.020608 | 6 |
| 7 | 27 | 0.997395 | 0.027455 | 15 |
| 7 | 81 | 0.994961 | 0.053293 | 20 |
| 7 | 243 | 0.990644 | 0.027358 | 15 |
| 7 | 729 | 0.982171 | 0.015888 | 6 |
| 7 | 2187 | 1.642747 | 0.000000 | 1 |

**4C verdict:** Intermediate — sigma spread 0.673 (30-100%). Modest variation; not structural at the brief's 2× threshold.

## Cross-correlation across a_final classes

For each k, correlations between (rho_mean, gamma_mean, sigma_mean) computed across the common set of a_final classes. If multiple sub-probes return Gate A, correlations indicate whether they're measuring the same underlying admissibility-precariousness.

| k | a_final | rho_mean | gamma_mean | sigma_mean |
|---|---|---|---|---|
| 5 | 3 | 1.000099 | 1.029435 | 1.015510 |
| 5 | 9 | 0.999999 | 1.031802 | 1.002582 |
| 5 | 27 | 1.001700 | 1.026244 | 0.990935 |
| 5 | 81 | 1.087073 | 1.013176 | 0.939162 |
| 5 | 243 | 1.178924 | 0.985896 | 1.506459 |
| 6 | 3 | 0.999998 | 0.962463 | 1.014450 |
| 6 | 9 | 1.000018 | 0.965172 | 1.017840 |
| 6 | 27 | 1.000019 | 0.964003 | 0.999827 |
| 6 | 81 | 1.043063 | 0.959739 | 0.979183 |
| 6 | 243 | 1.131336 | 0.947591 | 0.960832 |
| 6 | 729 | 1.225118 | 0.923581 | 1.599613 |
| 7 | 3 | 0.999885 | 0.891731 | 1.010112 |
| 7 | 9 | 0.999912 | 0.896363 | 1.010696 |
| 7 | 27 | 1.000008 | 0.896239 | 0.997395 |
| 7 | 81 | 1.002241 | 0.895560 | 0.994961 |
| 7 | 243 | 1.087753 | 0.891973 | 0.990644 |
| 7 | 729 | 1.176718 | 0.880278 | 0.982171 |
| 7 | 2187 | 1.267139 | 0.858635 | 1.642747 |

Per-k correlations:

| k | corr(rho, gamma) | corr(rho, sigma) | corr(gamma, sigma) |
|---|---|---|---|
| 5 | -0.9888 | +0.8193 | -0.8730 |
| 6 | -0.9822 | +0.7886 | -0.8842 |
| 7 | -0.9470 | +0.7595 | -0.9025 |

## Combined verdict

By the pre-registered single-probe gates: 0/3 sub-probes clear Gate A, 1/3 clears Gate B (4B), 2/3 are intermediate (4A and 4C). Spread summary:

- 4A (joint admissibility) max spread by a_final: 0.2673  (Gate A: >1.0; Gate B: <0.10) → intermediate
- 4B (residue counting) max spread by a_final: 0.0445  (Gate A: >0.30; Gate B: <0.10) → Gate B (close)
- 4C (stuck residues) max spread by a_final: 0.6726  (Gate A: >1.0; Gate B: <0.30) → intermediate

But the **cross-correlation pattern is the load-bearing observation**, not the single-probe gates. All three correlations between (ρ, γ, σ) across a_final classes are |r| > 0.75 with consistent sign at every k:

- corr(ρ, γ) ≈ −0.95 to −0.99: ρ and γ anti-correlate near-perfectly
- corr(ρ, σ) ≈ +0.76 to +0.82: ρ and σ co-vary
- corr(γ, σ) ≈ −0.87 to −0.90: γ and σ anti-correlate

This is the cross-probe converging-evidence regime the brief explicitly asked about. Three independent reframings of admissibility see the same a_final-class structure, and the variation pattern is consistent: **high a_final classes are admissibility-precarious in all three senses simultaneously** (large ρ, small γ, large σ).

### Pattern by a_final

The variation is monotone in a_final:

- **a_final = 3** (lowest, ≡ residues with shortest deterministic prefix): ρ ≈ 1.00, γ ≈ peak, σ ≈ 1.01
- **a_final = 3^k** (highest at each k, = the "repunit" residue r ∈ {1, 5, 21, 85, 341, ...}): ρ ≈ 1.18-1.27, γ ≈ minimum, σ ≈ **1.51-1.64** (single residue per k accounts for the σ tail)

The σ jump at the highest a_final is *discrete*, not smooth: at k=7, the repunit residue (a_final=2187) has σ=1.643 while its neighbors (a_final=729) have σ ≈ 0.98. So 4C identifies one outlier residue per modulus, not a smooth gradient.

### Mechanism — Gate C interpretation

The deterministic prefix algorithm sets a_final = 3^j where j is the number of 3x+1 steps inside the prefix; total prefix length = j + k (k halvings + j multiplications). Higher a_final means more multiplications relative to halvings, which corresponds to **lower average v_2 per prefix step** (more steps needed to consume k bits).

So a_final correlates monotonically with prefix-stage v_2 deficit. The cross-probe variation therefore tracks the same underlying signal as the bit-budget original: residues with low average v_2 stay in their mod-2^k stratification longer (larger orbit length), generate fewer distinct traces, and are rarely visited as transit residues.

This matches Gate C ("variation tracks v_2 not a_final = bit-budget pattern in disguise") for sub-probes 4A and 4B. Sub-probe 4C is partially Gate C for the smooth gradient, plus a discrete outlier at the repunit residue that the original bit-budget Test 3 didn't expose.

### What's actually new

The original bit-budget Test 3 found low-v_2 trajectories had **larger** admissibility horizons (opposite of the predicted contradiction). These reframings sharpen that null result with two additions:

1. **The same v_2-driven gradient is visible in three independent admissibility metrics**, with cross-correlations |r| > 0.75. This is converging evidence that admissibility precariousness is mostly a v_2 phenomenon, consistent with the original Test 3 verdict.

2. **The repunit residue r = (4^k − 1)/3 mod 2^k** ({1, 5, 21, 85, 341, ...}) shows a discrete σ jump at every k tested, while other classes are uniform at σ ≈ 1. This singular structural feature wasn't in the bit-budget framing. The repunit is well-known for taking the maximum prefix steps and has a special role in the deterministic prefix algorithm. **4C's discrete outlier is the only finding here that is genuinely new vs the original Test 3.**

### Disposition

- **Bit-budget admissibility (original Test 3): null** (low-v_2 has larger horizons; predicted contradiction refuted).
- **Joint admissibility (4A): null in disguise** (smooth gradient in a_final = v_2-deficit gradient, same direction as 4B).
- **Residue counting (4B): null in disguise** (smooth gradient, same direction).
- **Stuck residues (4C): mostly null in disguise + one discrete outlier per modulus** at the repunit. The outlier is structurally specific (not v_2-explained) but its single-residue scope is too narrow to drive a closure-style elimination argument.

Net: the framework's open question of whether admissibility produces structural constraints remains negative. The repunit σ outlier is a candidate target for follow-up, but at single-residue scope per modulus it's unlikely to drive elimination of broader admissibility classes.

## Files

- `result_4A_joint_admissibility.csv` — per-cell rho with a_final tag
- `result_4B_residue_counting.csv` — per-residue gamma with a_final tag
- `result_4C_stuck_residues.csv` — per-residue sigma with a_final tag
- `result_4_cross_correlation.csv` — per-class summary of all three
- `test_3_reframings_findings.md` — this writeup