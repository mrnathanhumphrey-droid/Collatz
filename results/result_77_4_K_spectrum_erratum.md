# R77.4 Erratum / Follow-up — K_k spectrum is NOT the rate operator

**Date:** 2026-05-04. Companion to [result_77_4_operator_shape.md](result_77_4_operator_shape.md). Documents an empirical finding that invalidates the framing (not the verdict) of R77.x.

## Headline

> **The level-k Markov transition matrix K_k mixes in essentially one step at every level k = 3..6.** Its second-largest eigenvalue has |λ_2| ~ 10⁻⁵ to 10⁻³, and it has **zero eigenvalues anywhere near 1/2** (no eigenvalue in [0.3, 0.7] at any k). Therefore the (1/2)^n decay of ε_n cannot be a spectral-gap phenomenon at the level of K_k.

R77.x (R77.1–R77.4) implicitly framed ε_n's envelope shape as a question about the spectrum of "the rate operator T". This follow-up shows that whatever T is, **it is not K_k** — the within-level transition operator on coprime residues mod 3^k.

## Empirical evidence

Direct computation via [result_77_4_K_spectrum.py](result_77_4_K_spectrum.py): build K_k as a float matrix, compute full spectrum via numpy, count eigenvalues in concentric bands around λ = 1/2.

| k | N_k = dim K_k | \|λ_2\| (second-largest) | closest eigenvalue to 1/2 | count in \|λ − 1/2\| < 0.1 |
|---|---|---|---|---|
| 3 | 18  | 3.19×10⁻⁶ | distance 0.4999968 | 0 |
| 4 | 54  | 2.95×10⁻⁵ | distance 0.4999706 | 0 |
| 5 | 162 | 3.38×10⁻⁴ | distance 0.4996618 | 0 |
| 6 | 486 | 1.23×10⁻³ | distance 0.4989341 | 0 |

The closest eigenvalue to 1/2 at every k is approximately λ_2 itself (real, very near 0), giving |λ_2 − 0.5| ≈ 0.5 by triangle inequality. **No K_k eigenvalue lives near 1/2 at any tested level.**

K_k has eigenvalue 1 (simple, the stationary) and a tight cluster near 0 (with |λ| growing slowly with k). This is a strong-mixing operator, not a slow-decay one.

## Why this invalidates R77.x's framing

R77.2 conjectured a "T_3 companion matrix with spectrum {1/2, 1/4, 1/8}". R77.3 falsified that 3-mode geometric ansatz over Q. R77.4 then asked which operator-shape (Jordan / log / power-law) governs the envelope ε_n.

In all of these, the implicit assumption was that **some operator on a finite-dimensional state space** has a spectral feature at λ = 1/2 that reproduces the (1/2)^n envelope. The companion matrix from the empirical recursion was DESCRIPTIVE — it codes the recursion that fits ε_n — but R77.x carried an unstated PRESUMPTION that this companion would relate to K_k or a natural operator derived from it.

**This empirical result rules that out.** K_k has no spectral feature near 1/2 at any level, so:
1. The companion-matrix spectrum from R77.3 (≈ {0.534, 0.144, −0.084}) is **not** a numerical approximation to a sub-spectrum of K_k.
2. The envelope (1/2)^n is **not** a spectral-gap phenomenon at level k.
3. R77.4's hypothesis classes (Jordan / log / power) were curve-fitting against the *envelope*, not against any operator-theoretic prediction grounded in K_k.

## What this does NOT change

- **R77.4 verdict (M) holds** for what it actually said: "of the five hypotheses fitted to e_n on n=2..6, no single shape dominates after small-sample correction; Jordan ruled out by direction." The verdict was about the envelope's curve-fitting, not about K_k.
- **R77.3 verdict (β) holds**: the 3-mode geometric ansatz over Q is approximate, not exact. Confirmed independently.
- **c = 7/45 rate-1/2 rigor status is UNCHANGED**: still gated on Tao Prop 1.17 effective C_A. The Nisoli-bypass route (R77.2) remains closed (R77.3) and the operator-shape route (R77.4) remains parked. This finding does not open or close any closure path.

## What this DOES change — reframing

The right operator-theoretic question is no longer "what spectrum does K_k have near 1/2?" but rather:

> **What is the operator governing the inter-level refinement π_k → π_{k+1}, and what does its spectrum look like near 1/2?**

The state space at level k+1 is a 3-fold refinement of level k's: each coprime residue r mod 3^k corresponds to {r, r+3^k, r+2·3^k} mod 3^{k+1}. The "lift" operator L_{k→k+1}: R^{N_k} → R^{N_{k+1}} sends π_k to the uniform extension over the 3 fibers; the actual π_{k+1} differs from L_{k→k+1}·π_k by a residual whose magnitude empirically scales like (1/2)^n.

The spectrum of this **residual operator** (or its companion on the projective limit) is what R77.x was actually trying to characterize. Computing it directly is the natural next probe; cf. R77.4 follow-up option (2) parked pending direction.

## Numerical aside: where does the 1/N decay in |λ_2| come from?

The trend |λ_2(K_k)| ≈ 0.04/N_k (rough fit: 3.19e-6, 2.95e-5, 3.38e-4, 1.23e-3 against N_k = 18, 54, 162, 486; ratio |λ_2|·N_k = 5.7e-5, 1.6e-3, 5.5e-2, 0.60) actually **grows superlinearly** in N_k — the chain becomes *less* sharply mixing as k grows, but stays well-mixing by absolute standards. At k=6 the spectral gap is 1 − 1.23×10⁻³ ≈ 0.9988, still essentially 1.

This monotone growth toward... something (maybe a finite asymptote, maybe slowly diverging) is itself an interesting object. But it is not 1/2 and shows no signs of approaching 1/2.

## Files

- [result_77_4_K_spectrum.py](result_77_4_K_spectrum.py) — script
- [experiments_output/result_77_4_K_spectrum_data.csv](experiments_output/result_77_4_K_spectrum_data.csv) — full eigenvalue listings for k = 3..6 (720 rows)
- [result_77_4_K_spectrum_erratum.md](result_77_4_K_spectrum_erratum.md) — this writeup

## Recommended next moves (parked pending direction)

1. **Build the inter-level residual operator** R_k: π_k ↦ π_{k+1} − L_{k→k+1}·π_k (after suitable embedding to a common space) and compute its spectrum across k. Spectral feature near 1/2 here would be the actual empirical signature R77.x was reaching for.
2. **Generating function on ε_n** — compute Σ ε_n z^n's analytic structure via Padé / partial-sum diagnostics to identify singularity type (pole vs branch cut vs essential). Sub-day; cheap.
3. **Re-examine T78.1–78.6 for inter-level vs within-level confusion** — verify none of the structural anchors silently assumed K_k's spectrum has the rate-1/2 feature. If any did, those theorems need re-statement (likely cosmetic, but should be checked).
