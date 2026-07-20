# Probe R4 — the edge density (thread 3 closing derivation) — **scaling laws NOT confirmed at finite L**

**Date:** 2026-07-20  CPU, dense L=2,3. Probe `probes/probe_thread3_R4.py`. Gates the edge-density law's three
pre-registered SHAPES. **Method validated; the pre-registered finite-L scalings are refuted — an honest negative,
with a clear diagnosis (the 1/θ pole is an L→∞ continuum feature, not present in the discrete finite-L band).**

Object: each band mode's exact shell contribution product_j = A_j·(λ_j − 1/3), A_j = ⟨1|r_j⟩⟨ℓ_j|v₀⟩ the
spectral amplitude (= ψ_kin·φ_j·g_j, per R3-C). Extracted robustly per k-pair via a grouped 2×2 biorthogonality
inversion (individual near-degenerate matching is singular; grouped cond ≈ 1.0–1.7 — clean).

## Method validation (the cross-check that says the extraction is sound)
The **partner** (real, near 1/3) product = **+0.144** → 3×0.144 = **0.432 ≈ the full-chain plateau 0.439**. So the
dense per-mode amplitude method is correct (and *more* accurate than R3's 2×2, which undershot to 0.114 by
dropping subdominant mixing). The band-mode numbers below are therefore trustworthy, not artifacts.

## R4-A — the 1/θ ratio: **REFUTED** (ratio 0.013, pre-reg 2.0)
At L=3, the two members of the k=1 doublet:
| mode | λ | amplitude A_j | product_j |
|---|---|---|---|
| k=1 member 1 | 0.235+0.183j | **2e-15 (DARK)** | ~5e-16 |
| k=1 member 2 | 0.238+0.183j | −1.0e-4+3.0e-4j | 6.5e-5 |
| k=2 member 1 | 0.020+0.184j | +1.6e-3+1.7e-2j | 6.3e-3 |
| k=2 member 2 | 0.028+0.172j | +2.5e-3−3.1e-2j | 1.1e-2 |

**k=1 total product 6.5e-5 ≪ k=2 total 5.0e-3** (ratio **0.013**, not 2). One doublet member is a **symmetry
dark state** (A≈0 — the symmetric agreement functional/init doesn't excite it), the other couples only weakly.
The top band mode barely couples to the agreement observable — the **opposite** of the 1/θ prediction (which
wanted k=1 the *larger*). L=2 is worse (super-critical + the seat formula θ_k=2πk/3^{L−1} mis-identifies the
top pair there).

## R4-B — L-invariant normalized density: **REFUTED**
(3^{L−1}/2π)·θ₁·|product|: L=2 → 9.2e-2, L=3 → 6.5e-5. Ratio L3/L2 = 7e-4 — nowhere near invariant (dominated by
the k=1 product collapsing from L=2's mis-identified mode to L=3's dark-doublet).

## R4-C — arg → ±π/2: **REFUTED**
arg(product_{k=1}) = −2.32 (L=3), −1.39 (L=2); π/2 = 1.571. Not at the pole's sine-transform phase.

## Diagnosis (constructive — this is where the honesty points)
The edge-density law derives 7/15 as the L→∞ **continuum** residue of a band whose mode density diverges as 1/θ
at the edge. **At finite L=2,3 the band is not a continuum — it is 1–2 discrete modes per k.** The k=1 "edge"
mode is the doublet, whose near-degeneracy splits it into one **dark** member and one weakly-coupled member, so
the discrete finite-L spectrum **cannot** exhibit the continuum 1/θ scaling the pre-registration assumed. The
laws (R4-A/B/C) are L→∞ statements tested against finite-L discrete data — a category mismatch, not (necessarily)
a false derivation. The continuum limit is exactly where the EP wall (L=4 unresolvable) blocks the numerics, so
**this gate cannot confirm or refute the edge-density law as stated** — it can only report that the finite-L
discrete band does not display the predicted shapes.

**What would test it:** either (i) the closed-form edge amplitude Ĝ from the D-FORM + Real-T1 normalizations
directly (pure symbolic, no eigen-extraction — the pen's route), or (ii) a per-mode coupling g_j·φ_j extracted
via per-band 2×2 blocks (kin, mode-j) rather than the spectral amplitude, in case the "coupling-overlap product"
convention differs from A_j(λ_j−1/3). The R4-D conventions are frozen below for either.

## R4-D — convention dump (frozen verbatim for the pen)
- basis: M = build_M_gen(3, L, 2, [2^−d]) dense; states (a, b, γ), a,b ∈ ⟨2⟩ mod 3^L, γ ∈ ℤ/3^L.
- right eigvecs r_j = columns of `eig(M)`; left ℓ_j = columns of `eig(Mᵀ)`, matched by eigenvalue; **bilinear
  (un-conjugated) inner product** ℓ·r; grouped 2×2 inversion S = LᵀR for near-degenerate pairs.
- readout ⟨1| = all-ones (a_m = 1ᵀMᵐv₀ = P(pair agrees to depth m)); init v₀ = δ(1,1,0) = index 0.
- amplitude A_j = ⟨1|r_j⟩·[(LᵀR)⁻¹ Lᵀv₀]_j; product_j = A_j(λ_j − 1/3); θ_k = 2πk/3^{L−1}; pair = 2 upper-half
  modes nearest the seat σ(θ_k); no ℓ²-renormalization of eigvecs.

## Status
Extraction **validated** (partner 0.144 → 3× = plateau). Edge-density pre-registered shapes (1/θ ratio 2:1,
L-invariance, arg ±π/2) **not confirmed at finite L** — the k=1 doublet has a dark member and weak coupling, so
the discrete finite-L band does not show the continuum 1/θ scaling. **This is a limit-vs-finite-L category gap,
not a clean refutation of the derivation**: the continuum test is blocked by the EP wall, so the edge-density law
must be closed symbolically (the pen's Ĝ). Honest negative; magnitudes 0-for-27 preserved; conventions frozen.
