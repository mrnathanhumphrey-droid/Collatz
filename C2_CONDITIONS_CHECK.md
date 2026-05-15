# C2_CONDITIONS_CHECK — Cut-and-project conditions for Syracuse encoding

**Date:** 2026-05-12. Cluster 2 cut-and-project probe, Phase 4.

Verification for the best candidate (Candidate B from Phase 2):
**G = ℝ, H = ∏'_p ℚ_p (rational adelic finite part), L = ℚ embedded diagonally.**
**W_Syracuse = (ℤ_3\*) × ∏_{p≠3} ℤ_p.**

This is the **support-only** scheme (unweighted indicator). For the weighted-comb scheme (Candidate B-weighted), see §(d) below.

---

## (a) π_G injective on L

L = ℚ ⊂ ℝ × ∏'_p ℚ_p, diagonal embedding q ↦ (q, q, q, …).

π_G : ℝ × ∏'_p ℚ_p → ℝ, π_G(q, (q_p)) = q (the Archimedean component).

For (q_1, …) ∈ L and (q_2, …) ∈ L with π_G(q_1, …) = π_G(q_2, …): both Archimedean parts equal q_1 = q_2 ∈ ℝ. Since the embedding ℚ ↪ ℝ × ∏'_p ℚ_p is diagonal, q_1 = q_2 forces equality on all components. **π_G is injective on L.** ✓

This is the BMP setup verbatim (page 35).

---

## (b) π_H(L) dense in H

π_H(L) = ℚ (image of the diagonal in the H-coordinate is ℚ embedded in ∏'_p ℚ_p).

**By the Strong Approximation Theorem for ℚ in the adeles**, ℚ is dense in ∏'_p ℚ_p (BMP, page 35: "the denseness of π_int(Q) in ∏ Q_p is equivalent to the Strong Approximation Theorem"). ✓

---

## (c) W_Syracuse regular

From Phase 3 (a):
- relatively compact: ✓ (compact)
- measurable: ✓ (clopen ⊂ Borel)
- non-empty interior: ✓ (clopen → interior = set)
- θ_H(∂W_Syracuse) = 0: ✓ (clopen → boundary = ∅, Haar measure = 0)

**W_Syracuse is a regular window.** ✓

---

## (d) Maximal density (for weak-model-set version, 1512.00912 Thm 5.9)

dens(L) · θ_H(W_Syracuse) = 1 · θ_{ℤ_3}(ℤ_3\*) · ∏_{p≠3} θ_{ℤ_p}(ℤ_p) = 1 · (2/3) · 1 = **2/3**.

Counting density of 3-coprime integers: |{n ∈ [−N, N] : gcd(n, 3) = 1}| / (2N) → 2/3 as N → ∞. ✓ Maximal density holds.

---

## (e) Aggregate: support scheme verification

All three regular-model-set cut-and-project conditions hold for the SUPPORT encoding. Theorem 5.9 (1606.08831 regular form) fires:

> γ_support = dens(L) · ω_{1_{W_Syracuse} * 1̃_{W_Syracuse}} = ω_{1_{W_Syracuse} * 1̃_{W_Syracuse}}
> γ̂_support = dens(L)² · ω_{|1̂_{W_Syracuse}|²} = ω_{|1̂_{W_Syracuse}|²}

This is the diffraction of δ_{3-coprime integers}, supported on ℚ ∩ (ℤ[1/3]/ℤ)-axis i.e. rationals with denominator a power of 3.

**The unweighted support encoding cleanly satisfies Theorem 5.9.**

---

## (f) Weighted-comb scheme verification (Candidate B-weighted)

Now the load-bearing question: does the weighted version with h built from Syracuse μ work?

Required: h ∈ P_K(H) for Theorem 5.4, i.e., positive-definite continuous compactly supported function on H.

The natural h(y) = ρ(y_3) · ∏_{p≠3} 1_{ℤ_p}(y_p) where ρ is the Radon–Nikodym density of μ_∞ on ℤ_3.

(f1) **Compact support of h:** Yes — h vanishes outside (ℤ_3) × ∏_{p≠3} ℤ_p which is compact. ✓

(f2) **Continuity of ρ on ℤ_3:** UNKNOWN. Markov stationary distributions on totally disconnected compact groups can be absolutely continuous, singular, or atomic. For Syracuse the existence of a continuous (or even L^1) density is not established. The closest available results (R76/R77) show μ_n's mass concentrates on the principal-unit "+" coset, which is consistent with a *singular* measure on ℤ_3 supported on a sub-lattice structure — NOT consistent with a continuous density.

(f3) **Positive-definite ρ on ℤ_3:** By Bochner B.3, ρ ∈ P_K(ℤ_3) means ρ̂ is a positive measure on Ẑ_3 = ℤ[1/3]/ℤ. This means the Fourier coefficients of μ_∞ form the moments of a positive measure on the dual. This is a strong analytic statement — equivalent to a global polynomial-in-A character-sum control over all 3-power moduli.

(f4) **Verdict on h ∈ P_K(H):** UNVERIFIED. Conditions (f2) and (f3) are essentially the polynomial-in-A Fourier-decay statement on μ_n itself, dressed up in different language. The hypothesis Theorem 5.4 needs IS the conclusion the framework would deliver. **This is a Mode H target-object circularity: the encoding hypothesis is the closure target.**

---

## (g) Adversarial checks

### (A1) Set vs measure
Confirmed: Theorem 5.9 in its regular form (1606.08831) is for ω_{1_W} (unweighted comb). The weighted extension is Theorem 5.4 (1606.08831) / Theorem 5.4 (1512.00912), which requires h ∈ P_K(H) — the load-bearing hypothesis we cannot verify for Syracuse.

### (A2) Density condition
**Maximal density holds** (computed in (d) above: dens = 2/3). So even if other strict conditions failed, the weak-model-set Theorem 5.9 of 1512.00912 would fire for the support comb. The density condition is NOT the failure point — it's the weight encoding.

### (A3) Bochner foundation
For the unweighted support comb δ_Λ: γ_support is positive-definite by construction (autocorrelation of a measure). Bochner applies directly.

For the weighted comb ω_h: positive-definiteness of γ_μ = μ_∞ * μ̃_∞ on ℝ is by construction. Bochner applies. **But** to verify that Theorem 5.4's hypothesis h ∈ P_K(H) on the internal-space function holds requires positive-definiteness of ρ on ℤ_3, which is non-trivial — that's exactly the polynomial-in-A statement we're trying to establish.

### (A4) Pisot/Salem warning
The diffraction the framework delivers for the unweighted support is **pure-point on rationals with cubefree-at-3 denominators** — this is algebraic-resonance Fourier structure, not Rajchman decay (consistent with Feng/Kahane's Pisot warning since the dilation by 3 is integer-Pisot). The polynomial-in-A target needs the WEIGHTS' decay, not the support's diffraction.

### (A5) Target-object trap (Mode H)
**Active and load-bearing.** The probe asks "is μ_n cut-and-project encodable?" but the framework's natural object is either:
- δ_Λ (unweighted comb) — Λ = 3-coprime integers, regular weak model set; diffraction known
- ω_h (weighted comb) — requires h ∈ P_K(H) which is the target

The autocorrelation γ_μ = μ_∞ * μ̃_∞ on ℝ is a valid object of analysis, but the cut-and-project framework's mechanical delivery requires going through ω_h. The encoding pathway μ_n → γ_μ → γ̂_μ via Theorem 5.9 only works in the unweighted regime.

### (A6) Inherited-claim audit (Mode E)
**BMP's "super-singular model set" framing** is structurally correct and rigorously stated (BMP §"Further connections", pages 34–36; verified verbatim in C2_THEOREM_5_9_HYPOTHESES.md (h)). The intuition transplants TECHNICALLY for the support layer — 3-coprime integers really are a (weak) model set with maximal density in the BMP super-singular sense, and Theorem 5.9 delivers the diffraction of their unweighted Dirac comb.

BMP's framing does NOT, however, derive a polynomial Fourier bound for weighted Markov-chain-stationary measures on ℤ_3\*. The framing extends to a SET on ℤ via the BMP scheme; extension to a MEASURE on the SET requires Theorem 5.4 weighted-comb form, whose hypothesis is the target property.

So the inherited claim is technically rigorous as stated (super-singular model SET for Syracuse SUPPORT), but does NOT automatically carry over to Syracuse measure μ_n. The framing's transplant is **partial** in exactly the way Mode E warns about — the supporting framework is correct, but its scope is the support not the weights.

---

## Aggregate verdict

| Layer | Condition | Status |
|---|---|---|
| Support, π_G injective | yes | ✓ |
| Support, π_H(L) dense | yes (Strong Approx) | ✓ |
| Support, W regular | yes (W is clopen) | ✓ |
| Support, max density | 2/3 = dens(L)·θ_H(W) | ✓ |
| Support, Thm 5.9 fires | YES | ✓ (regular model set scheme) |
| Weights, h ∈ P_K(H) for Thm 5.4 | requires continuity + positive-definiteness of ρ on ℤ_3 | UNVERIFIED — circular with closure target |
| Weights, Thm 5.4 fires | depends on h ∈ P_K | NO (mode H) |

**The support encoding works cleanly. The weighted-comb encoding hits Mode H target-object circularity.**

Since the polynomial-in-A closure target IS the weighted |μ̂_n(ξ)|² object — not the unweighted support diffraction — the framework's mechanical machinery does not fire automatically for the closure target.

This is **H_C2_ENCODING_PARTIAL** disposition: the support layer satisfies the cut-and-project conditions and Theorem 5.9 applies cleanly there. The weight layer fails because Theorem 5.4's hypothesis on h is essentially the closure target itself.
