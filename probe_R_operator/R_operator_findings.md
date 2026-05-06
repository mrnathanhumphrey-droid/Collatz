# R_k forcing operator — singular value spectrum

**Date:** 2026-05-05 (post-compact). Probe target: forcing operator R_k = P_{W_{k+1}} ∘ K_{k+1} ∘ L_k (V_k -> W_{k+1}), where W_{k+1} = V_{k+1} \ominus L_k(V_k). Tests whether the rate-determining structure missing from the prior R-operator probe (which used the W -> W self-map and returned zero) lives in the V -> W forcing block.

## Convention note

Brief specified `V_k = span(pi_k)` (1-dim) but also requested top-20 singular values — those two are inconsistent (a 1-dim input gives at most one nonzero singular value). Adopted the interpretation V_k = full level-k coprime function space (n_k = 2·3^{k-1}-dim), consistent with the prior R-operator probe's V_k notation. R_k is therefore an n_k-dim → 2·n_k-dim linear map and yields n_k singular values.

## Construction

- Row-vector convention: `u @ L_k @ K_{k+1} @ P_{W_{k+1}}^T` transports a level-k function through the lift, the next-level Markov dynamics, and the orthogonal projection onto W_{k+1}.
- L_k is the uniform fiber lift L_k[i, j] = 1/3 if coprime_{k+1}[j] mod 3^k = coprime_k[i] else 0; row-stochastic.
- W_{k+1} basis P_W obtained from QR of L_k^T (orthonormal complement of L_k's row-space inside V_{k+1}). Sanity check max|P_W @ L_k^T| reported per k below — should be ~ machine eps.
- R_k as a column-vector linear map V_k -> W_{k+1} is `R_k = P_W @ K_{k+1}^T @ L_k^T`, shape (dim(W_{k+1}), n_k). SVD via `numpy.linalg.svd(R_k, full_matrices=False)`.

## Pre-registered reference rates

| label | value | source |
|---|---|---|
| rho_slow (order-3 real root) | 0.826934 | `result_renormalization_recurrence_fits.csv` |
| order-3 complex pair magnitude | 0.192080 | same |
| order-2 top root | 0.312245 | same |
| rate-1/2 (legacy walked-back claim) | 0.5000 | R75/R76 |
| 1.0 | broken projection diagnostic | walk-back gate |
| 0.0 | trivial forcing diagnostic | walk-back gate |

## Per-k summary

| k | dim(V_k) | dim(W_{k+1}) | sigma_1 | sigma_2 | sigma_3 | sigma_20 | sigma_min | orth_resid |
|---|---|---|---|---|---|---|---|---|
| 5 | 162 | 324 | 6.705541e-01 | 6.705541e-01 | 6.703471e-01 | 6.581581e-01 | 1.640694e-17 | 2.78e-17 |
| 6 | 486 | 972 | 6.706155e-01 | 6.706155e-01 | 6.705925e-01 | 6.691314e-01 | 2.199363e-18 | 2.78e-17 |
| 7 | 1458 | 2916 | 6.706223e-01 | 6.706223e-01 | 6.706198e-01 | 6.704561e-01 | 5.091373e-19 | 2.78e-17 |

## Top 20 singular values per k

### k = 5  (V_k dim = 162, W_6 dim = 324)

| rank | sigma | log10(sigma) |
|---|---|---|
| 1 | 6.7055407475e-01 | -0.1736 |
| 2 | 6.7055407475e-01 | -0.1736 |
| 3 | 6.7034711010e-01 | -0.1737 |
| 4 | 6.7034711010e-01 | -0.1737 |
| 5 | 6.6952450702e-01 | -0.1742 |
| 6 | 6.6952450702e-01 | -0.1742 |
| 7 | 6.6891304674e-01 | -0.1746 |
| 8 | 6.6891304674e-01 | -0.1746 |
| 9 | 6.6730529763e-01 | -0.1757 |
| 10 | 6.6730529763e-01 | -0.1757 |
| 11 | 6.6631714822e-01 | -0.1763 |
| 12 | 6.6631714822e-01 | -0.1763 |
| 13 | 6.6399724492e-01 | -0.1778 |
| 14 | 6.6399724492e-01 | -0.1778 |
| 15 | 6.6267720167e-01 | -0.1787 |
| 16 | 6.6267720167e-01 | -0.1787 |
| 17 | 6.5975021192e-01 | -0.1806 |
| 18 | 6.5975021192e-01 | -0.1806 |
| 19 | 6.5815806049e-01 | -0.1817 |
| 20 | 6.5815806049e-01 | -0.1817 |

### k = 6  (V_k dim = 486, W_7 dim = 972)

| rank | sigma | log10(sigma) |
|---|---|---|
| 1 | 6.7061549886e-01 | -0.1735 |
| 2 | 6.7061549886e-01 | -0.1735 |
| 3 | 6.7059245939e-01 | -0.1735 |
| 4 | 6.7059245939e-01 | -0.1735 |
| 5 | 6.7050036664e-01 | -0.1736 |
| 6 | 6.7050036664e-01 | -0.1736 |
| 7 | 6.7043136542e-01 | -0.1736 |
| 8 | 6.7043136542e-01 | -0.1736 |
| 9 | 6.7024764828e-01 | -0.1738 |
| 10 | 6.7024764828e-01 | -0.1738 |
| 11 | 6.7013303618e-01 | -0.1738 |
| 12 | 6.7013303618e-01 | -0.1738 |
| 13 | 6.6985862859e-01 | -0.1740 |
| 14 | 6.6985862859e-01 | -0.1740 |
| 15 | 6.6969898805e-01 | -0.1741 |
| 16 | 6.6969898805e-01 | -0.1741 |
| 17 | 6.6933528394e-01 | -0.1744 |
| 18 | 6.6933528394e-01 | -0.1744 |
| 19 | 6.6913142554e-01 | -0.1745 |
| 20 | 6.6913142554e-01 | -0.1745 |

### k = 7  (V_k dim = 1458, W_8 dim = 2916)

| rank | sigma | log10(sigma) |
|---|---|---|
| 1 | 6.7062232663e-01 | -0.1735 |
| 2 | 6.7062232663e-01 | -0.1735 |
| 3 | 6.7061976615e-01 | -0.1735 |
| 4 | 6.7061976615e-01 | -0.1735 |
| 5 | 6.7060952504e-01 | -0.1735 |
| 6 | 6.7060952504e-01 | -0.1735 |
| 7 | 6.7060184506e-01 | -0.1735 |
| 8 | 6.7060184506e-01 | -0.1735 |
| 9 | 6.7058136864e-01 | -0.1735 |
| 10 | 6.7058136864e-01 | -0.1735 |
| 11 | 6.7056857349e-01 | -0.1736 |
| 12 | 6.7056857349e-01 | -0.1736 |
| 13 | 6.7053787333e-01 | -0.1736 |
| 14 | 6.7053787333e-01 | -0.1736 |
| 15 | 6.7051997026e-01 | -0.1736 |
| 16 | 6.7051997026e-01 | -0.1736 |
| 17 | 6.7047906372e-01 | -0.1736 |
| 18 | 6.7047906372e-01 | -0.1736 |
| 19 | 6.7045606283e-01 | -0.1736 |
| 20 | 6.7045606283e-01 | -0.1736 |

## sigma_1 vs reference rates

| k | sigma_1 | closest reference | distance |
|---|---|---|---|
| 5 | 0.670554 | rho_slow (order-3 real root) (0.8269) | 0.1564 |
| 6 | 0.670615 | rho_slow (order-3 real root) (0.8269) | 0.1563 |
| 7 | 0.670622 | rho_slow (order-3 real root) (0.8269) | 0.1563 |

## Cross-k decay/growth of sigma_1

| k -> k+1 | sigma_1(k) | sigma_1(k+1) | ratio | log_3 ratio |
|---|---|---|---|---|
| 5 -> 6 | 0.670554 | 0.670615 | 1.000092 | 0.0001 |
| 6 -> 7 | 0.670615 | 0.670622 | 1.000010 | 0.0000 |

## Structural findings (not pre-registered, but striking)

The data carries clean structure beyond the pre-registered checks. Worth flagging
explicitly before any interpretation of sigma_1 alone:

### 1. Exact rank = 2·n_k/3 at every k

| k | n_k | rank(R_k) at tol 1e-10 | n_k/2 (for ref) | 2·n_k/3 (matches!) |
|---|---|---|---|---|
| 5 | 162  | 108  | 81  | 108 |
| 6 | 486  | 324  | 243 | 324 |
| 7 | 1458 | 972  | 729 | 972 |

R_k has an exact (n_k/3)-dimensional null space. The non-rank singular values
drop precipitously (sigma_{2n_k/3} ~ 0.49, sigma_{2n_k/3 + 1} ~ 1e-15: 14
orders of magnitude). This is structural, not numerical.

### 2. Singular value pairing: each sigma appears (at least) twice

Every consecutive pair sigma_{2i-1} = sigma_{2i} agrees to ~16 digits.
This holds across all three k values and through the full nontrivial
spectrum. The forcing operator has a 2-fold algebraic degeneracy at every
singular value — likely from the K_- = σ K_+ σ chain symmetry (proved
elsewhere; sibling 3x±1 study) propagating into R_k via the lift's mod-3
reflection structure. Worth a follow-up algebraic verification.

### 3. Spectrum supported on a closed band [~0.488, ~0.671]

Not a discrete spectrum dominated by sigma_1. The nontrivial singular values
fill a band:

| k | sigma_1 (top of band) | sigma_{rank} (bottom of band) | spread |
|---|---|---|---|
| 5 | 0.6705541 | 0.4879767 | 0.1826 |
| 6 | 0.6706155 | 0.4879530 | 0.1827 |
| 7 | 0.6706223 | 0.4879504 | 0.1827 |

Both endpoints are k-stable to 4-5 significant figures; the band width is
constant. This is a **continuous-band-like spectrum** at the scales tested,
not an isolated eigenvalue picture. The "rate" interpretation in terms of
sigma_1 alone is therefore the wrong frame — what governs the inter-level
forcing is the full band's distribution.

### 4. sigma_1 and sigma_min (nontrivial) both k-stable

| k | sigma_1 | sigma_{rank} |
|---|---|---|
| 5 | 0.6705541 | 0.4879767 |
| 6 | 0.6706155 | 0.4879530 |
| 7 | 0.6706223 | 0.4879504 |

Differences: sigma_1 increases by 6.1e-5 then 6.8e-6 (ratio 0.11);
sigma_min decreases by 2.4e-5 then 2.6e-6 (ratio 0.11). Both are
geometrically converging to fixed limits — strong evidence for a
well-defined k → ∞ asymptotic spectrum. Geometric extrapolation of
sigma_1 → 0.67063 ± 1e-5; sigma_min → 0.48795 ± 1e-5. Neither matches
3/(2√5) = 0.67082 (off 2e-4) or any other obvious closed form I checked.

## Verdict

**Brief's third walk-back gate fires:** the singular values cluster but
do not match rho_slow ≈ 0.83, rate 1/2, or any order-2/3 recurrence root
within 0.05. Per the brief, reporting and pausing for analysis rather
than over-interpreting.

What the probe DOES establish:

1. **R_k as a V → W forcing operator is non-trivial** (sigma_1 ≈ 0.67,
   rank = 2n_k/3, rank-stable across k). Walk-back gate "sigma_1 ≈ 0
   means trivial forcing" does NOT fire.
2. **R_k is a contraction with sigma_1 < 1** (no broken-projection
   warning).
3. **The spectrum is band-supported, not single-eigenvalue.** Trying to
   match a *single* rate (rho_slow, 1/2, etc.) to sigma_1 misses the
   structure — what's converging is a closed band on roughly [0.488, 0.671].
4. **The rate-determining structure for ε_k convergence is not localized
   in sigma_1.** If ε_k convergence rate ≈ 0.83 is real, it would have
   to come from a non-singular-vector property of R_k (e.g., a specific
   linear functional applied to specific input directions, or interaction
   with subsequent renormalization steps).

What's open:

- **What is 0.6706 algebraically?** The k-stability is striking; the
  number doesn't obviously match Plancherel quantities, simple rationals,
  or roots of small-degree polynomials I checked. Worth investigating
  via the dominant singular vector u_1 — its support pattern in W_{k+1}
  may reveal what subspace it's selecting. (Vectors saved in npz files.)
- **Does R_k iterated produce ε_k?** The natural next probe: build the
  composite map V_k → W_{k+1} → ... up several levels. If ε_k arises
  from a specific functional applied to (R_{k-1} ∘ R_{k-2} ∘ ... )(u_0),
  the slow rate 0.83 would emerge from the composition, not single R_k.
- **Why exact rank 2n_k/3?** The kernel of R_k has dimension exactly
  n_k/3. What level-k functions, when uniformly fiber-lifted and evolved
  one K step, stay inside the lifted subspace L_k(V_k)? Likely tied
  to the Z/3 algebra of the q=3 case.

sigma_1 across k = [5, 6, 7]: 0.670554, 0.670615, 0.670622  (extrap. → ~0.67063)
sigma_min(nontrivial) across k: 0.487977, 0.487953, 0.487950  (extrap. → ~0.48795)
rank(R_k)/n_k across k: 0.667, 0.667, 0.667  (= 2/3 exactly)

## Timings

| k | t_K | t_lift+P_W | t_build_R | t_SVD |
|---|---|---|---|---|
| 5 | 0.0s | 0.0s | 0.0s | 0.0s |
| 6 | 0.0s | 0.1s | 0.0s | 0.1s |
| 7 | 0.1s | 0.8s | 0.1s | 0.7s |

## Files

- `R_k5_singular_values.csv` — full singular spectrum
- `R_k5_dominant_vectors.npz` — top-20 u/v vectors
- `R_k6_singular_values.csv` — full singular spectrum
- `R_k6_dominant_vectors.npz` — top-20 u/v vectors
- `R_k7_singular_values.csv` — full singular spectrum
- `R_k7_dominant_vectors.npz` — top-20 u/v vectors
- `R_operator_findings.md` — this writeup
