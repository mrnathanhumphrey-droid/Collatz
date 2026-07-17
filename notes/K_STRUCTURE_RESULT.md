# K_k structure: spectrum is {1, 0, ..., 0}; rate-1/2 of ε_n is OBSERVABLE-driven, not Markov-spectral

**Date:** 2026-05-15. Direct corroboration via exact-Q computation. Strengthens R77.4 erratum from "empirical |λ_2| ≈ 10⁻³ (numerical noise on ill-conditioned matrix)" to an exact structural statement.

## Statement

> **K_k has spectrum {1, 0, 0, ..., 0} — exactly one nonzero eigenvalue at λ = 1 (the stationary), all others zero with Jordan structure at 0 of chain length k. Equivalently, K_k^k has rank 1; K_k mixes to stationary in EXACTLY k Markov steps from any starting distribution.**

The W_{k-1} subspace (3-fiber-zero-mean functions inside V_k) is mapped to 0 by K_k exactly. The image of K_k lives in the lift subspace T(V_{k-1}) ⊂ V_k, and K_k|_{T(V_{k-1})} ≅ K_{k-1} under the lift isomorphism.

## Proof

Structural row-equality. The level-k Markov chain transitions from r ∈ (Z/3^k)* are:

  r ↦ ((3r + 1) · 2^{-v}) mod 3^k, weighted 2^{-v}.

For r' = r + 3^{k-1} (the next lift in the same 3-fiber at level k-1):

  3r' + 1 = 3r + 1 + 3^k ≡ 3r + 1 (mod 3^k).

So 3r' + 1 ≡ 3r + 1 mod 3^k, and ((3r' + 1) · 2^{-v}) mod 3^k = ((3r + 1) · 2^{-v}) mod 3^k. The transition target depends only on r mod 3^{k-1}, not on which lift.

Therefore K_k(r, s) = K_k(r + 3^{k-1}, s) = K_k(r + 2·3^{k-1}, s) for all s. The 3 rows in each 3-fiber are EXACTLY equal.

This means K_k factors as
  K_k = ψ_k ∘ T_sum,
where T_sum: V_k → V_{k-1} is the 3-fiber summing operator (kernel = W_{k-1}, the 3-fiber-zero-mean subspace) and ψ_k: V_{k-1} → V_k inserts the level-(k-1) data into level-k rows.

By marginal consistency (Syracuse Markov is coherent under mod-3^{k-1} reduction):

  T_sum ∘ ψ_k = K_{k-1}.

By the rank-trace-trick (nonzero spec(AB) = nonzero spec(BA) for any compatible A, B):

  **nonzero spec(K_k) = nonzero spec(T_sum ∘ ψ_k) = nonzero spec(K_{k-1}).**

By induction down to K_1: K_1 is 2×2 with both rows equal (since 3·1+1 ≡ 1 and 3·2+1 ≡ 1, both mod 3), hence rank(K_1) = 1, spec(K_1) = {1, 0}.

Therefore **nonzero spec(K_k) = {1}** for all k ≥ 1. ∎

## Numerical verification (exact-Q)

[`K_structure_verify.py`](K_structure_verify.py) computed at V_MAX=16, exact rationals:

| k | N_k | rank(K_k) | rank(K_k²) | rank(K_k³) | trace(K) | trace(K²) | trace(K³) |
|---|-----|-----------|------------|------------|----------|-----------|-----------|
| 2 | 6   | 2         | 1          | —          | 1        | 1         | 1         |
| 3 | 18  | 6         | 2          | 1          | 1        | 1         | 1         |

Trace pattern: trace(K^m) = sum of eigenvalues to power m. With spec = {1, 0, 0, ..., 0}, trace(K^m) = 1 for all m. **Verified exactly to all checked m.**

Rank pattern: rank(K_k^m) = N_{k-m} = 2·3^{k-m-1} for m = 1..k-1, and rank(K_k^k) = 1. So **K_k converges to stationary in EXACTLY k steps, not via spectral decay**.

Jordan structure at 0: chain length k (Jordan block of nilpotency index k accounts for the rank pattern).

## What this means for the c=7/45 subdominant rate

R77.4 erratum's "|λ_2| ≈ 10⁻³ to 10⁻⁵" finding for K_k was NUMERICAL NOISE from float arithmetic on an ill-conditioned matrix. The true spectrum of K_k contains no information about decay rates beyond "converges in k steps."

Therefore **the rate-1/2 envelope of ε_n = S_n − 7/15 is NOT a Markov-chain spectral-gap phenomenon at any finite level.** Combined with R77.5 (multi-resolution decomposition V_{k+1} = T(V_k) ⊕ W_k, with R_k structurally in W_k by marginal consistency) and R77.6 (E(z) := Σ ε_n z^n has a BRANCH CUT singularity at z=2 from Padé diagonal monotone convergence, not a simple pole), the picture is:

> **Rate-1/2 lives in the structured overlap of the moment observable φ_n with the multi-resolution residual sequence {R_k}, not in any operator's spectrum at finite truncation.**

Specifically (per R77.5 §5):
  ε_n = Σ_k ⟨φ_n, lift_n(R_k)⟩.

The L² magnitudes alone give ‖lift_n(R_k)‖² · 3^k ≈ 0.155 (R77.5 Stage 1 exact), with ‖lift_n(R_k)‖² = ‖R_k‖² · 3^{n-k-1} = 0.155·3^{n-2k-1}. A general Cauchy-Schwarz bound gives |⟨φ_n, lift_n(R_k)⟩| ≤ ‖φ_n‖·‖lift_n(R_k)‖, which produces a 3^{-1/2}·... bound, not 2^{-n}.

The 2^{-n} rate is a SPECIFIC INNER-PRODUCT CANCELLATION between φ_n (the bilinear pair-form moment functional from R76) and the residual R_k pattern, not a magnitude bound. The branch cut at z=2 (R77.6) is the analytic encoding of this cancellation pattern.

## Why this is useful

1. **Confirms R77.5's displacement** of R77.2's "find λ_2 = 1/2 at finite truncation" framing. The finite-truncation operator simply does not have such structure.
2. **Confirms R77.4 erratum's spirit** while correcting its quantitative reading — the |λ_2| values reported there are numerical noise, the true spectrum is exactly {1, 0, ..., 0}.
3. **Localizes the open question:** the c=7/45 subdominant rate (1/2)^n is fully encoded in the φ_n × R_k inner-product structure, NOT in any operator spectrum. Three concrete attack paths remain:
   - **Inverse-limit transfer operator**: build the Tao recursion on L²(Ẑ_3^×) (the projective completion of ⊕_k V_k) and locate the continuous-spectrum branch cut at λ=1/2.
   - **Generating-function Tauberian transfer**: prove E(z) has a singularity at z=2 of definite type (power-law α vs logarithmic), then apply Flajolet-Sedgewick Ch. VI singularity analysis to get rigorous |ε_n| ≤ C·n^α·2^{-n}.
   - **Direct bilinear-pair inner-product computation**: extend M_n(η) = ⟨φ_n, π_n⟩-like structure analytically and prove the (1/2)^n cancellation pattern explicitly.

## Mode-of-action note for future work

R77.4 erratum's reported eigenvalue trend |λ_2(K_k)| ≈ 10⁻⁶, 10⁻⁵, 10⁻⁴, 10⁻³ at k=3,4,5,6 (growing with k) **was numerical noise on an ill-conditioned matrix**. The true rank/eigenvalue structure (exact-Q) shows K_k is rank N_{k-1}, not full rank, and the only nonzero eigenvalue is 1.

This is a useful reminder: **for ill-conditioned matrices (rank-deficient, large dim), trust exact-Q rank/trace checks over float eigenvalue computation.** numpy.linalg.eigvals on a singular matrix returns small-magnitude artifacts that look like real eigenvalues but aren't.

## Files

- [`K_W_restricted_spectrum.py`](K_W_restricted_spectrum.py) — preliminary float probe: showed K_W ≈ 0 at ~10⁻¹⁷ across k=2..5
- [`K_structure_verify.py`](K_structure_verify.py) — exact-Q verification at k=2, 3: rank, trace(K^m), Jordan chain length
- [`experiments_output/K_W_restricted_spectrum.json`](experiments_output/K_W_restricted_spectrum.json) — float spectrum data
- [`experiments_output/K_structure_verify.json`](experiments_output/K_structure_verify.json) — exact-Q rank/trace data
- [`result_77_4_K_spectrum_erratum.md`](result_77_4_K_spectrum_erratum.md) — earlier empirical erratum (this writeup is the stronger structural follow-up)
- [`result_77_5_inter_level_residual.md`](result_77_5_inter_level_residual.md) — multi-resolution decomposition (companion result)
- [`result_77_6_generating_function.md`](result_77_6_generating_function.md) — branch cut at z=2 (companion result)
