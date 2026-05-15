# C2_DISPOSITION — Cluster 2 cut-and-project probe, top-level summary

**Date:** 2026-05-12. Cluster 2 cut-and-project encoding probe.

---

## DISPOSITION: **H_C2_ENCODING_PARTIAL** (with Mode H target-object trap)

The support layer of Syracuse μ_n cleanly satisfies the cut-and-project conditions in the BMP super-singular scheme. The weighted layer, which is what the polynomial-in-A closure target needs, fails because the natural encoding requires the weight function h to be in P_K(H) — positive-definite continuous compactly supported — which is essentially the target Fourier-decay property dressed in different language.

**Cluster 2 does NOT mechanically deliver the polynomial-in-A bound for Syracuse μ_n.** It does deliver a known result for the unweighted support set (3-coprime integers), but that's not the closure target.

Pre-registered favored outcome H_C2_ENCODING_PARTIAL at ~35% — confirmed.

---

## Phase summary

### Phase 0 — PDF access ✓
All load-bearing transcripts at C:/tmp/crystal/. Verbatim quotation of Theorem 5.9 (both forms), cut-and-project scheme definition, Definition 5.2, Theorem 5.4, Theorem 5.7, Remark 5.11, and BMP super-singular framing achieved. See C2_THEOREM_5_9_HYPOTHESES.md.

### Phase 1 — Two Theorem 5.9 forms identified

1. **Regular form (1606.08831 §5):** G = ℝ^d, H compactly generated LCA, L lattice in G × H, W regular window (rel. compact, measurable, non-empty interior, ∂W has Haar-zero), π_G injective on L, π_H(L) dense in H. Conclusion: γ̂ = dens(L)² · ω_{|1̂_W|²}.

2. **Weak-model-set form (1512.00912 §5.5, Theorem 5.9):** G σ-compact LCA, L lattice, W rel. compact measurable (no regularity), maximal density along a van Hove sequence. Conclusion: same diffraction formula.

3. **Weighted-comb form (Theorem 5.4, 1606.08831 / 1512.00912):** (G, H, L) cut-and-project, h ∈ P_K(H) (positive-definite continuous compactly supported). Conclusion: ω̂_h = dens(L) · ω_{ĥ}. Holds **without projection assumptions** (1606.08831 verbatim).

### Phase 2 — Candidate (G, H, L) verdicts

| Candidate | Verdict | Why |
|---|---|---|
| A: G = ℤ_3, H = ∏_{p≠3} ℤ_p, L = ℤ | FAILS | ℤ is dense in Ẑ; L is not discrete in ℤ_3 × ∏_{p≠3} ℤ_p |
| B: G = ℝ, H = ∏'_p ℚ_p, L = ℚ (BMP) | VALID scheme | Strong Approximation; matches BMP F_k construction |
| B-weighted: B + h = ρ for ρ = dμ/dHaar on ℤ_3 | PARTIAL (circular) | h ∈ P_K(H) ≡ polynomial-in-A target |
| C: G = H = ℤ_3 | FAILS | No non-trivial lattice in profinite × profinite |
| D: B with finite-level windows | Works at finite level | Window encodes set not weights |

Candidate B is the natural BMP-canonical embedding and IS a valid cut-and-project scheme. Its window for the Syracuse support is W_Syracuse = ℤ_3\* × ∏_{p≠3} ℤ_p.

### Phase 3 — Window status

W_Syracuse = ℤ_3\* × ∏_{p≠3} ℤ_p:
- Relatively compact: ✓ (compact, clopen)
- Measurable: ✓ (Borel)
- Non-empty interior: ✓ (clopen)
- ∂W has Haar measure 0: ✓ (boundary is empty)
- **REGULAR WINDOW.** ✓

⋏(W_Syracuse) = {n ∈ ℤ : gcd(n, 3) = 1} = 3-coprime integers in ℤ. Density 2/3 = dens(L) · θ_H(W), maximal density confirmed.

### Phase 4 — Conditions verification

For the SUPPORT layer (unweighted indicator 1_{W_Syracuse}):
- (a) π_G injective on L: ✓ (BMP diagonal embedding of ℚ)
- (b) π_H(L) dense in H: ✓ (Strong Approximation Theorem)
- (c) W regular: ✓ (W is clopen)
- (d) Maximal density: ✓ (2/3 verified)

**Theorem 5.9 fires for the unweighted Dirac comb of the 3-coprime integer support. Disposition for SUPPORT: confirmed cut-and-project encoding (regular weak model set in BMP's super-singular sense).**

For the WEIGHTED layer (ω_h with h built from Syracuse μ_∞):
- Requires h ∈ P_K(H) (Theorem 5.4)
- Reduces to: ρ = dμ_∞/dHaar exists, is continuous, and is positive-definite on ℤ_3
- These are NOT independently established for Syracuse μ_n
- The closest we know is that μ_n is mass-concentrated on the "+" coset (R76/R77), consistent with μ_∞ being SINGULAR — not absolutely continuous
- Even if ρ existed, positive-definiteness of ρ on ℤ_3 ≡ μ_n character sums are moments of a positive measure on the dual ≡ polynomial-in-A target

**Theorem 5.4 does NOT fire for Syracuse μ_n. The encoding's load-bearing hypothesis is the closure target itself (Mode H circularity).**

### Phase 5 — Adversarial checks

(A1) Set vs measure: Set encoding (Theorem 5.9) ≠ measure encoding (Theorem 5.4). Syracuse is a measure.
(A2) Density: NOT the failure point — maximal density holds (2/3).
(A3) Bochner: applies, but to apply Theorem 5.4 we need ρ ∈ P_K(ℤ_3) which is the target.
(A4) Pisot/Salem: Algebraic-resonance pure-point diffraction at rationals with cubefree-at-3 denominators is what Theorem 5.9 delivers for the unweighted support — consistent with Pisot scaling-by-3. NOT the same as a polynomial-in-A Rajchman-type decay for Syracuse μ_n weights.
(A5) Target-object trap (Mode H): **active and load-bearing.** Framework's natural objects are δ_Λ (unweighted set) and ω_h (positive-definite weight function on H), neither of which is Syracuse μ_n unconditionally.
(A6) Inherited-claim audit (Mode E): BMP super-singular framing transplants technically rigorously for the SUPPORT, partially for the measure. Documented in C2_THEOREM_5_9_HYPOTHESES.md (h).

---

## What Cluster 2 actually delivers for Syracuse

**The cut-and-project framework delivers the diffraction of the 3-coprime-integer SET on ℝ:** γ̂ supported on rationals with cubefree-at-3 denominator, intensities computable via BMP's product formula. This is a known and well-understood result analogous to BMP's F_k for k = 1.

**The cut-and-project framework does NOT deliver a polynomial-in-A bound on |μ̂_n(ξ)|² for the Syracuse Markov stationary measure.** The natural extension (Theorem 5.4 weighted-comb form) requires the weight function ρ = dμ_∞/dHaar on ℤ_3 to be in P_K(ℤ_3) — positive-definite continuous compactly supported. This is essentially equivalent to the polynomial-in-A target itself, hence circular.

---

## Comparison with the five-probe map (POLYNOMIAL_IN_A_LANDSCAPE.md)

| Probe | Framework | Object | Verdict |
|---|---|---|---|
| 1 | L²-flattening | smooth IFS | strategic collapse in discrete |
| 2 | SL_2 / Furstenberg | P¹-stationary | T1 transfer fails |
| 3 | Cocycle Dolgopyat | C^{1+α} cocycle | UNI fails (linear exclusion) |
| R1 | ARHW + smoothing | self-conformal | structural pincer |
| R2 | Drift conditions | transient mixing | wrong target |
| **C2** | **Cut-and-project / BMP** | **Markov stationary** | **support: fires; weights: Mode H circular** |

Cluster 2 doesn't add a sixth structural failure — it adds a **partial result** with a clean mechanical delivery for the support (consistent with BMP F_k) and a circularity at the weight layer. The closure target (weights) is not delivered, but the encoding analysis is informative.

---

## Recommendation for next move

**Routing forward:**

1. **Cluster 2 closes for closure purposes** — no polynomial-in-A bound on μ_n weights via cut-and-project. The framework's mechanical delivery stops at the support layer.

2. The Cluster 1 direct attack (Tauberian arc on Σ_k |π̂_k(ξ)|² generating series + Bourgain-Konyagin discrete sum-product) remains the sole forward route, as identified in POLYNOMIAL_IN_A_LANDSCAPE.md.

3. **Possible salvage of Cluster 2 (note for the Tauberian arc):** the support encoding gives a Plancherel decomposition: Σ_x μ_∞(x) χ(x) = Σ_{x ∈ 3-coprime} ρ(x) χ(x), where the **sum** factors are controlled by the BMP super-singular Bragg measure ω_{|1̂_W|²} (computable explicitly). The weights ρ are then a multiplicative factor over the support. This Plancherel decomposition is structurally informative for the Tauberian-arc generating-series approach but does not itself unblock the polynomial-in-A bound.

4. **No new structural failure beyond what the five-probe map already documented.** Cluster 2's role is to formally close the BMP/PSF/cut-and-project family as the sixth structural framework family that has been honestly probed.

---

## Mode lessons learned

- **Mode H (target-object trap) closes at the framework boundary:** Cluster 2 frames Syracuse μ_n as a *cut-and-project measure*, but the framework's mechanical theorem requires h ∈ P_K(H), which is essentially the closure target. The circularity becomes the disposition.

- **Mode E (inherited-claim audit) PASSES:** BMP's super-singular model set framing transplants rigorously for the SUPPORT layer (not analogy, actual theorem). It does not extend to the weight layer except as an asymptotic-Plancherel decomposition.

- **The unweighted-support diffraction is a real result, just not the closure target.** Worth recording in the Tauberian-arc input toolkit if needed.

---

## Files

- C2_THEOREM_5_9_HYPOTHESES.md — verbatim hypothesis statements (Phase 1)
- C2_SYRACUSE_CANDIDATES.md — five candidate (G, H, L) schemes evaluated (Phase 2)
- C2_SYRACUSE_WINDOW.md — window identification + weighted-comb status (Phase 3)
- C2_CONDITIONS_CHECK.md — three conditions + adversarial checks (Phase 4)
- C2_DISPOSITION.md — this document (top-level summary)
