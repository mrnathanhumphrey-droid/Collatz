# C2_SYRACUSE_WINDOW — Window identification for Syracuse μ_n

**Date:** 2026-05-12. Cluster 2 cut-and-project probe, Phase 3.

Building on Candidate B from Phase 2: (G, H, L) = (ℝ, ∏'_p ℚ_p, ℚ).

---

## (a) Syracuse support → window

**Level-n support:** μ_n is supported on (ℤ/3^n)\* ⊂ ℤ/3^n. Lifted to ℤ_3 (level n=∞), this is ℤ_3\* = ℤ_3 \ 3ℤ_3, the 3-adic units.

**Window in BMP coordinates** (H = ∏'_p ℚ_p):

W_Syracuse = (ℤ_3\*) × ∏_{p≠3} ℤ_p ⊂ H

Then:
- ⋏(W_Syracuse) := π_G(L ∩ (G × W_Syracuse)) = {q ∈ ℚ : q ∈ ℤ_3\* AND q ∈ ℤ_p for all p ≠ 3} = {n ∈ ℤ : gcd(n, 3) = 1} = **3-coprime integers in ℤ**.

This is the asymptotic SUPPORT of Syracuse μ_n (as n → ∞ and after lifting to ℝ).

## Density / regularity of W_Syracuse

- **Relatively compact in H:** ℤ_3\* is closed in ℤ_3 (it's the complement of the open subgroup 3ℤ_3 — wait, 3ℤ_3 is open AND closed in ℤ_3, so ℤ_3\* = ℤ_3 \ 3ℤ_3 is open AND closed, i.e., clopen). Each ∏_{p≠3} ℤ_p factor is compact. So W_Syracuse is **compact** (hence relatively compact). ✓
- **Measurable:** Yes, clopen sets are Borel. ✓
- **Non-empty interior:** ℤ_3\* is open in ℤ_3 (clopen ⊃ open), and ∏_{p≠3} ℤ_p is open in ∏'_p ℚ_p in restricted-product topology (basic open set). So W_Syracuse has non-empty interior. ✓
- **Boundary θ_H(∂W) = 0:** Since W_Syracuse is clopen, ∂W_Syracuse = ∅, and θ_H(∅) = 0. ✓

**So W_Syracuse is a REGULAR WINDOW for the support.** The unweighted comb δ_{3-coprime integers} = ω_{1_{W_Syracuse}} is a regular model set, and the strict regular-model-set Theorem 5.9 (1606.08831 form) applies.

This is structurally just BMP's F_1 (squarefree at p=3 only, i.e., 3-coprime numbers). The diffraction is computable explicitly: γ̂ supported on rationals with denominator a power of 3 (cubefree-at-3 since the window is 1 mod 3), with specific Bragg intensities.

**For the unweighted SUPPORT of Syracuse, Cluster 2 fires cleanly.** This isn't the closure target object, but it establishes that Step 1 of the cut-and-project encoding succeeds for the support layer.

---

## (b) Syracuse weights → window-supported function

The actual closure target is **|μ̂_n(ξ)|² for Syracuse μ_n with Markov weights**.

The weights are NOT uniform over (ℤ/3^n)\*. From `result_78_FINAL.md` and the R76/R77 P_+/P_− decomposition:
- μ_n splits along the principal-unit subgroup structure of (ℤ/3^n)\*.
- The Fourier mass of f(u) = e_q(c·4^u) concentrates on {a ≡ 1 mod 3} ⊂ ℤ/3^r — same "+" coset structure.
- Numerically, μ_n is non-uniform with a documented 7/15 ↔ 8/15 split between (P_+, P_−) at level n=1.

To encode the weights via the cut-and-project scheme, we'd need a function h: H → ℂ such that ω_h ≡ Σ_n μ_∞(n) δ_n on G = ℝ.

**Natural candidate for h:**

h(y) = ρ(y_3) · ∏_{p≠3} 1_{ℤ_p}(y_p)

where ρ: ℤ_3 → ℝ_{≥0} is the **Radon–Nikodym density of μ w.r.t. Haar measure on ℤ_3**.

Open analytic questions about ρ:
- **Does ρ exist as an L^1 function?** I.e., is μ absolutely continuous w.r.t. Haar on ℤ_3? Markov stationary distributions on ℤ_3 can be absolutely continuous, singular continuous, or atomic. For Syracuse, this is not established.
- **Is ρ continuous?** This is the regularity question — for the weighted-comb Theorem 5.4 to fire, we need h ∈ P_K(H), which requires continuity of ρ on ℤ_3.
- **Is ρ positive-definite as a function on the LCA group ℤ_3?** Positive-definiteness of ρ on the LCA group ℤ_3 means ρ̂ is a positive measure on Ẑ_3 = ℚ_3/ℤ_3 = ℤ[1/3]/ℤ.

## What positive-definiteness of ρ means for closure

ρ ∈ P_K(ℤ_3) means:
- (a) ρ continuous on ℤ_3
- (b) ρ compactly supported (automatic; ℤ_3 compact)
- (c) Σ_i Σ_j c_i c̄_j ρ(x_i − x_j) ≥ 0 for all finite (x_i, c_i) — i.e., ρ is a positive-definite function on the LCA group ℤ_3

(c) is the load-bearing condition. By Bochner B.3, this is equivalent to ρ̂ being a positive measure on the dual ℤ[1/3]/ℤ.

ρ̂(ξ) for ξ ∈ ℤ[1/3]/ℤ is essentially the Fourier transform Σ_n μ_∞(n) e_{3^∞}(−ξn), which at level n becomes precisely the μ̂_n(ξ) we're trying to bound.

So **ρ ∈ P_K(ℤ_3)** is equivalent to saying μ_n's Fourier transforms admit a uniform positive-definite extension — which is a strong structural property that essentially encodes the polynomial-in-A target indirectly.

**This is the Mode H circularity made explicit:** the cut-and-project encoding of Syracuse μ_n requires a property (h ∈ P_K) that is essentially the target Fourier-decay property.

---

## (c) Weak version: drop the weights to the unweighted indicator

If we restrict attention to the **support indicator** 1_{W_Syracuse}, Theorem 5.9 fires cleanly:
- γ̂_support = dens(L)² · ω_{|1̂_{W_Syracuse}|²}
- dens(L) = 1 (BMP normalization, ℚ in ℝ × ∏'_p ℚ_p has cocompact fundamental domain [0,1] × ∏_p ℤ_p of volume 1)
- 1̂_{W_Syracuse}(ξ) factors over the BMP product: for ξ = (ℓ_∞, (ℓ_p)) ∈ ℝ × ∏̂'_p ℚ_p,

|1̂_{W_Syracuse}(ξ)| = [δ_{0,(ℓ_p)_{p≠3}}] · |1̂_{ℤ_3\*}(ℓ_3)|

where the bracket is 1 iff all (ℓ_p)_{p≠3} are zero (i.e., we're on the rationals-with-denominator-a-power-of-3 axis), else 0.

1̂_{ℤ_3\*}(ℓ_3) on ℚ_3/ℤ_3:
- ℓ_3 = 0: |1̂| = θ_{ℤ_3}(ℤ_3\*) = 2/3
- ℓ_3 = 1/3 (or any character of order 3 on ℤ_3 — i.e., trivial on 3ℤ_3): |1̂| = ζ-character computation. Specifically 1̂_{ℤ_3\*} on a character χ of order 3 is ∫_{ℤ_3\*} χ(x) dx = (∫_{ℤ_3} − ∫_{3ℤ_3}) χ = 1·[χ = 1] − (1/3)·[χ|_{3ℤ_3} = 1] = (for χ non-trivial on ℤ_3 but trivial on 3ℤ_3) 0 − 1/3 = −1/3.
- ℓ_3 = 1/9 (character of order 9, non-trivial on 3ℤ_3): |1̂| = 0 (since χ averaged over any 3-adic coset of a subgroup it's non-trivial on, integrates to 0). Actually need careful calculation here: more generally |1̂_{ℤ_3\*}(ξ)| for ξ ∈ ℤ[1/3]/ℤ has specific Gauss-sum form.

**This is the cleanly-computable diffraction of the unweighted support.** It does NOT bound |μ̂_n(ξ)|, only |Σ_{n ∈ 3-coprime, |n| ≤ N} χ(n)| / N as N → ∞.

---

## Disposition for Phase 3

The window W_Syracuse = (ℤ_3\*) × ∏_{p≠3} ℤ_p is:
- A **regular window** in H = ∏'_p ℚ_p (clopen, hence rel. compact, measurable, open interior, empty boundary)
- The natural support window for the asymptotic Syracuse support set (3-coprime integers in ℤ)

But the Syracuse measure μ_n has **non-uniform weights** that are NOT given by an indicator 1_W. To encode the weights, we need a function h on H, which in turn requires:
- The Radon–Nikodym density ρ of μ w.r.t. Haar on ℤ_3 to exist
- ρ continuous (or at least Riemann integrable) on ℤ_3
- For Theorem 5.4: ρ positive-definite, i.e., ρ ∈ P_K(ℤ_3)

These conditions on ρ are NOT established for Syracuse μ and are essentially equivalent to the polynomial-in-A target itself. **The window-on-support story works (regular weak model set, Theorem 5.9 fires for the unweighted indicator). The weight encoding hits Mode H target-object circularity.**
