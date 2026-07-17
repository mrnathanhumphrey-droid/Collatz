# Inter-level Tao transfer U_n + twisted endomorphism Phi_omega: spectrum is CONTINUOUS-ON-CIRCLE, not discrete CC pair

**Date:** 2026-05-15. Follow-up to the K_k structural lemma (K_k has spectrum {1, 0, ..., 0} via 3-fiber row equality + marginal consistency). Built the alternative inter-level operator Phi_omega = T^omega ∘ U_n: V_n → V_n via twisted-fiber projection, hunting for a discrete CC pair carrying empirical period-9.2 oscillation.

## Headline

> **Two structural findings + one negative.**
>
> (1) **The Fourier-side Tao transfer operator U_n: V_n^Fourier → V_{n+1}^Fourier maps V_n entirely into W_n.** Exact (proven by 3rd-root-of-unity phase cancellation): for any f ∈ V_n, T_sum(U_n(f)) = 0 because Σ_{a=0,1,2} e^{-2πi a·2^{-v}/3} = 0 for any v. Numerically verified at n=2,3,4,5 (residual ~10^{-13}). So U_n's entire image lives in the 3-fiber-zero-mean complement.
>
> (2) **U_n's singular spectrum is non-trivial with pair structure.** Singular values come in pairs (multiplicity 2 each) — signature of complex-conjugate symmetry under ξ → -ξ. At n=4, sigma(U_n) spans 0.77 to 1.29 with at least 8 distinct values. This is a richer structure than K_k's trivial {1, 0} spectrum.
>
> (3) **Phi_omega = T^omega ∘ U_n: V_n → V_n has complex spectrum, but it's CONTINUOUS-ON-CIRCLE, not a discrete CC pair.** At n=5, the top 10 eigenvalues of Phi_{omega_3} all have |z| ≈ 0.319 (tight cluster) with arguments spanning ~[-1.28, -0.81] rad continuously distributed — NOT concentrated at any θ ≈ 2π/9.2 ≈ 0.683 rad. Similarly Phi_{omega_3^2} clusters on a circle of radius ≈ 0.587.
>
> **Period-9.2 CC pair is NOT in this operator family at finite truncation.** Consistent with R77.6 + PADE_NUMERICAL_DISPOSITION's branch-cut / continuous-spectrum reading. The CC pair, if asymptotically real, must be either (a) a continuous-spectrum resonance (saddle point in resolvent), (b) a feature of a DIFFERENT operator family (likely the bilinear pair-form T_M acting on M_n, not U_n on μ̂_n), or (c) a finite-n transient that doesn't survive to the limit.

## Structural lemma: U_n's image is in W_n

**Claim:** for any f: (Z/3^n)* → C, the function (U_n f): (Z/3^{n+1})* → C satisfies Σ_{a=0,1,2} (U_n f)(ξ + a·3^n) = 0 for every ξ ∈ (Z/3^n)*.

**Proof.** U_n is the Fourier-side Tao recursion:

  (U_n f)(ξ') := (1/Z) Σ_{v≥1} 2^{-v} · e^{-2πi ξ' 2^{-v}/3^{n+1}} · f(ξ' · 2^{-v} mod 3^n).

For ξ' = ξ + a·3^n with a ∈ {0, 1, 2}:
- The argument of f: (ξ + a·3^n) · 2^{-v} mod 3^n = ξ · 2^{-v} mod 3^n (same for all a, since a·3^n · 2^{-v} ≡ 0 mod 3^n).
- The phase: e^{-2πi (ξ + a·3^n) 2^{-v}/3^{n+1}} = e^{-2πi ξ 2^{-v}/3^{n+1}} · e^{-2πi a · 2^{-v}/3} = e^{-2πi ξ 2^{-v}/3^{n+1}} · ω_3^{-a · (2^{-v} mod 3)}.

Since 2^{-v} mod 3 ∈ {1, 2} (coprime to 3), Σ_{a=0}^{2} ω_3^{-a·c} = 0 for any c ∈ {1, 2}.

So the sum over a vanishes term by term in v. ∎

**Numerical confirmation:** ||T_sum ∘ U_n|| ≈ 10⁻¹³ across n=2,3,4,5 (machine epsilon).

## Phi_omega spectrum across n

Defined Phi_omega := T^omega ∘ U_n where T^omega: V_{n+1} → V_n is the **twisted fiber selection**:

  T^omega(g)(ξ) := (1/3) [g(ξ) + ω · g(ξ + 3^n) + ω² · g(ξ + 2·3^n)].

For ω = 1: T^omega = (1/3) · T_sum → Phi_1 = 0 on Image(U_n) by structural lemma.
For ω = ω_3, ω_3²: non-trivial projections onto Z3-eigenspaces of the 3-fiber translation.

Spectrum of Phi_omega at n=2..5 (top singular value of each):

| n   | dim V_n | top |λ| of Phi_{ω_3}  | top |λ| of Phi_{ω_3²} |
|-----|---------|-----------------------|------------------------|
| 2   | 6       | 0.294 (3 distinct ×2) | 0.526 (rank=2 ×2 +)    |
| 3   | 18      | 0.316 (9 distinct ×2) | 0.578                  |
| 4   | 54      | 0.318 (continuous)    | 0.588                  |
| 5   | 162     | 0.319 (continuous)    | 0.587                  |

**Across-n trend: top moduli converge to 0.319 (ω_3 case) and 0.587 (ω_3² case)**. These are STABLE asymptotic radii.

**Arguments / phases:** at n=5, top 10 Phi_{ω_3} eigenvalues have args ∈ {-1.01, -1.08, -1.11, -0.99, -1.18, -0.91, -0.89, -1.21, -0.81, -1.28} rad. Continuously distributed in a band around -1.05 rad; no concentration at θ = 2π/9.2 ≈ 0.683.

**Conclusion: the spectrum is a CONTINUOUS DISTRIBUTION ON A CIRCLE.** At finite n, this manifests as N_n discrete points clustering toward the circle; in the inverse limit n→∞, this becomes continuous spectrum on a circle of radius 0.319 (for Phi_{ω_3}) or 0.587 (for Phi_{ω_3²}).

This is the signature of CONTINUOUS spectrum, matching R77.6's branch-cut reading.

## Why no discrete CC pair

The Phi_omega operator family captures the "twisted residual after one Tao transfer step at level n+1, projected back to level n via 3-fiber-translation twist." Its finite-n spectrum has dim(V_n) = 2·3^{n-1} eigenvalues clustering on a circle.

A DISCRETE CC pair would correspond to a SINGLE complex eigenvalue (with conjugate) that is **isolated** in the spectrum. For Phi_omega's spectrum to have an isolated eigenvalue at, say, 0.319 · e^{±i·0.683} (the period-9.2 candidate), there would need to be a discrete eigenvalue separated from the other points. There isn't.

The probe DEFINITIVELY rules out the discrete-CC-pair picture for the family {Phi_omega: ω ∈ {ω_3, ω_3²}}.

## The radii 0.319 and 0.587

Neither matches:
- c=7/45 = 0.156 (subdominant amplitude)
- T_lead's eigenvalue 43/45 ≈ 0.956
- Empirical Hadamard ρ ≈ 0.984 at n=10..13

But ratio: 0.587 / 0.319 ≈ 1.840 ≈ √(2π/π)? Not obvious.

0.319² + 0.587² = 0.102 + 0.345 = 0.446. Not clean.

0.319 ≈ 1/π = 0.318 (very close). 0.587 ≈ 7/12 = 0.583 (very close, 0.7% gap). These might be coincidences.

If 0.319 and 0.587 are exact rationals at the inverse limit, they're not immediately interpretable in terms of the c=7/45 structure. Future probe: compute at higher n (n=6, 7) and check rate of convergence + exact-rational identification.

## What this probe ruled out

1. **Period-9.2 CC pair as discrete eigenvalue of Phi_omega family.** ✗
2. **Top |eigenvalue| of Phi_omega matching empirical 0.984 or 43/45.** ✗
3. **A simple square endomorphism on V_n with the c=7/45 asymptotic rate.** ✗

## What this probe established

1. **U_n maps V_n → W_n exactly** (3-fiber cancellation). Structural fact.
2. **U_n has non-trivial pair-structured spectrum** (richer than K_k).
3. **Phi_omega family has continuous-spectrum-on-circle structure** at n→∞.
4. **The asymptotic circle radii 0.319 (ω_3) and 0.587 (ω_3²) are stable across n** — suggests these are well-defined limits of finite-n operators.

## Next probe candidates

The c=7/45 CC pair, if it exists, may live in:

**Option D1: Bilinear pair-form operator T_M directly.** M_n(η) = Σ_ξ μ̂_n(ξ) μ̂_n*(ξη) is bilinear in μ̂_n. The recursion M_{n+1} ← T_M(M_n) is non-square. Build T_M as a matrix; compute its spectrum at small n. T_lead's 43/45 emerges as a 2x2 projection; the full T_M on V_n might have CC pair eigenvalues.

**Option D2: Iterate Phi_omega and track angles.** Compute Phi_omega^k for k = 2, 3, ..., 9, look at how eigenvalue arguments rotate. If a stable rotation rate emerges across n at iterate k, that's the CC structure.

**Option D3: Phi_omega's INVARIANT SUBSPACES.** Decompose V_n via Phi_omega's eigenspaces. Project onto the dominant eigenspace; iterate; see if the projected dynamics has discrete CC structure.

**Option D4: Direct attack on T_lead extension.** R77 §3 derived T_lead from leading Off_lin at v + v' = 2. Extend to include v + v' = 3, 4, ... (sub-leading off-diagonal). The "extended T_lead" might have CC pair eigenvalues. This is the project-internal cross_freq machinery.

## Files

- [`interlevel_U_spectrum.py`](interlevel_U_spectrum.py) — U_n singular value probe (showed pair structure + T_sum @ U_n = 0)
- [`interlevel_twisted_endomorphism.py`](interlevel_twisted_endomorphism.py) — Phi_omega complex spectrum probe (this writeup)
- [`experiments_output/interlevel_U_spectrum.json`](experiments_output/interlevel_U_spectrum.json)
- [`experiments_output/interlevel_twisted_endomorphism.json`](experiments_output/interlevel_twisted_endomorphism.json)
- [`K_STRUCTURE_RESULT.md`](K_STRUCTURE_RESULT.md) — K_k structural lemma (companion)
- [`T_LEAD_CORRECTED_DISPOSITION.md`](T_LEAD_CORRECTED_DISPOSITION.md) — T_lead 43/45 result (foundation)
- [`PADE_NUMERICAL_DISPOSITION.md`](PADE_NUMERICAL_DISPOSITION.md) — empirical CC pair evidence (the target structure)
