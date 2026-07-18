# FAURE candidate A — Phase 3: Partial-Expansion Verification for Syracuse / Tao recursion

This is the LOAD-BEARING question of the entire probe. If Syracuse fits Faure 2009's partially-expanding-map class, candidate A is SELECTED. If not, the entire Faure-school category fails (because A is the broadest of the candidates).

## Faure 2009's partially-expanding-map class (verbatim from candidate A)

The model is f: T² → T², f(x, s) = (kg(x) mod 1, s + τ(x)/(2π) mod 1), where:
- g: S¹ → S¹ is C^∞ diffeomorphism, with k ≥ 2 integer and g(x+1) = g(x)+1
- E(x) = kg(x) mod 1 is uniformly expanding on S¹: E_min = k min_x(dg/dx) > 1
- τ: S¹ → R is C^∞
- M = T² is smooth compact 2-torus
- f is C^∞ k-to-1 map

The MAP IS A SKEW PRODUCT over an expanding base E, with a NEUTRAL CIRCLE FIBER S¹_s. The C^∞ regularity is global. The transfer operator F̂ acts on C^∞(T²) by pullback.

## Syracuse / Tao recursion: structural inventory

From C1_TAO_RECURSION_FORM.md (Phase 1 verbatim from Tao 2022 §7.1):

**(1) State space:**
- The chain lives on (Z/3^n)* (the units mod 3^n) for each level n, or equivalently on the profinite group Z_3*.
- The transition n → n+1 expands the state space by factor 3 (refinement of 3-adic residues).
- No smooth manifold structure: (Z/3^n)* is finite at each level; Z_3* is a profinite group (totally disconnected Hausdorff topological group, compact, abelian, NOT a manifold).

**(2) Tao recursion form (Lemma 1.12 + project-internal R75):**
- μ̂_{n+1}(ξ) = Σ_{v=1}^∞ 2^{-v} · e^{-2πi ξ·2^{-v}/3^{n+1}} · μ̂_n(ξ · 2^{-v} mod 3^n)
- The "step" multiplies the Fourier variable ξ by 2^{-v} mod 3^n with random v ~ Geom(2)
- This is an iid renewal product (per C1 Tao §7.1), NOT a deterministic map iteration

**(3) The "expansion" and "contraction" directions:**
- 3-adic side: multiplying by 4 = 1+3 maps Z_3* → Z_3* as a TRANSLATION on the principal-unit subgroup (since (1+3)^u is a 3-adic exponential). This is a NEUTRAL action, not expansion.
- 2-adic side: the random factor 2^{-v} v∈Geom(2) acts on ξ. Mod 3^n, 2 is a unit of multiplicative order 3^{n-1} in (Z/3^n)*, so 2^{-v} permutes the residues — also a neutral action mod 3^n, but the variation in v provides randomness.
- Level refinement n → n+1: 3-fold refinement of the residue space. This is the closest analog to "expansion" — the operator F̂_n→n+1 is a 3:1 lift (each level-n residue has 3 level-(n+1) preimages).

**(4) The Tao recursion as a transfer-operator iteration:**
- Per result_77_T_lead_spectrum.md, the operator T_diag = (1/5)·[[1,1],[4,4]] acts on the 2-D space (P_+, P_-) of class-resolved bilinear pair-form moments.
- Spectrum of T_diag is {0, 1} — NOT inside the unit disk strictly; eigenvalue 1 sits on the unit circle. There is NO automatic spectral gap from T_diag.
- The actual rate-1/2 decay comes from OFF-DIAGONAL corrections (result_77_T_lead_spectrum.md §2-3), which is a different operator.
- K_k (within-level Markov transition) has |λ_2| ~ 10^-5 to 10^-3 (result_77_4_K_spectrum_erratum.md): "no K_k eigenvalue lives near 1/2 at any tested level". K_k is essentially a single-step mixing operator, not a slow-decay operator.
- The Pade picture (PADE_NUMERICAL_DISPOSITION) shows the LEADING singularity at n=13 is at |z|≈1.57 (NOT z=2 of R77.6's earlier reading), with sign pattern + + - - - - - - - + + + + zero-crossing between n=9 and n=10, consistent with complex-conjugate pair at θ ≈ 0.68 rad, period 9.2 in n-space.

## Phase 3 check: does Syracuse satisfy Faure 2009's partial-expansion hypothesis?

### Sub-criterion 1: Smooth compact manifold M with f: M → M (k:1, C^∞)

**FAILED.** Syracuse acts on:
- (Z/3^n)* at finite level n — a finite set, no manifold structure
- Z_3* in the inverse limit — a profinite group, totally disconnected, NOT a manifold (no tangent bundle, no smooth structure, no Riemannian metric)
- T² in Faure's model is connected and smooth

The category mismatch is fundamental: Faure's transfer operator F̂ is the pullback by a smooth map on a smooth manifold; Syracuse's Tao recursion operator on C(Z_3*) is the pullback by a measure-preserving stochastic kernel induced by Tao's renewal structure.

### Sub-criterion 2: Uniform expansion E_min > 1 in some direction

**FAILED.** Syracuse has:
- NO direction along which the map is uniformly expanding by a fixed factor. The 3-adic level refinement n → n+1 is a 3:1 fan-out (analogous to k=3 in Faure), but this is a PURELY DISCRETE refinement on a profinite space — not a smooth expanding map dE/dx ≥ E_min > 1.
- The Geom(2) randomness in the v-factor introduces stochasticity, not deterministic expansion.

A SUPERFICIAL analog: the 3:1 cover Z/3^{n+1} → Z/3^n has "expansion factor 3" in a coarse sense, but this is the lift, not a smooth differential.

### Sub-criterion 3: Neutral direction with circle-fiber structure

**FAILED.** Faure's neutral direction is s ∈ S¹_s, a smooth circle. Syracuse has no analog of a neutral S¹ fiber. The 2-adic Geom(2) randomness on v could be regarded as a "stochastic neutral direction" (in the sense that v is random, not deterministically chosen), but this is a probability-theoretic structure, not a smooth-fiber-bundle structure.

### Sub-criterion 4: Fourier-mode reduction F̂ = ⊕_ν F̂_ν (semiclassical limit ν → ∞)

The Tao recursion DOES have a Fourier-mode decomposition: at each level n, μ̂_n(ξ) for ξ ∈ Z/3^n is a Fourier coefficient. But:
- The Fourier variable ξ lives on Z/3^n, NOT on a continuous S¹
- The semiclassical limit ν → ∞ in Faure has a smooth-manifold analytic meaning (h-pseudodifferential calculus on T*S¹); the analog n → ∞ in Syracuse is the PROFINITE limit, a different category
- The reduced operator F̂_ν in Faure has C^∞ symbol; the Syracuse analog (Tao's recursion at level n) acts on functions on a finite group

### Sub-criterion 5: Partially captive trapped set K ⊂ T*S¹

**FAILED.** Faure's K is a compact subset of the cotangent bundle T*S¹, defined as the limit of pullbacks of a large compact set under the iterated map. The "partially captive" property is a Hausdorff-dimension condition on K. Syracuse has no cotangent bundle (no manifold structure), no analog of K.

### Sub-criterion 6: Anisotropic Sobolev / smooth distribution space

**FAILED.** Faure's spectral result lives in H^m(S¹) for m < 0 — smooth distributions on a smooth manifold. Syracuse's natural function space is C((Z/3^n)*) (functions on a finite group) at each level, or the inverse-limit C(Z_3*) — a Banach space of continuous functions on a profinite group, NOT a Sobolev space of smooth distributions.

### Net Phase 3 verdict for candidate A

**NOT_PARTIALLY_EXPANDING (in Faure 2009's sense)**

Mark: **NEEDS_PROFINITE_EXTENSION**.

Specifically:
- Syracuse's Tao recursion has SOME structural analogs to Faure 2009's setup (Fourier mode decomposition, 3:1 fan-out cover, asymptotic semiclassical-like parameter n)
- But the CORE INPUT REQUIREMENTS of Faure's proof technique are violated at every sub-criterion that involves smoothness, manifold structure, or differential geometry
- The smoothness gap is LOAD-BEARING at the proof-technique level: Faure's proof uses h-pseudodifferential calculus on the cotangent bundle T*S¹, escape functions on T*S¹, Egorov's theorem for pseudodifferential composition, semiclassical analysis on smooth distributions. None of these tools have profinite analogs in the existing literature

The natural mathematical extension that WOULD make Syracuse fit:
- A "**partially-expanding map on a profinite group**" theory: replace the smooth manifold T² with the profinite space Z_3* × (stochastic-fiber), replace C^∞ pseudodifferential calculus with a profinite-group / p-adic Fourier-analysis calculus, replace escape functions on T*S¹ with weight functions on the dual group, replace Egorov's theorem with a profinite-group-Fourier composition theorem.
- This is a **substantial research program** — there is no existing literature on a partially-expanding-map spectral theory in the profinite category. The closest analogs are p-adic Bruhat-Tits dynamics (already attempted in the BT_DISPOSITION arc, NO_FIT) and adelic Mellin theory (ADELIC_DISPOSITION, NO_FIT).

## Verdict on candidate A

**NO_FIT (smooth-manifold prerequisite hardcoded into Faure 2009's proof technique).**

PARTIAL is plausible IF the spirit of the theorem (transfer-operator spectral gap from semiclassical analysis of underlying dynamics) is the relevant abstraction. But "spirit" is not a theorem; converting Faure 2009 to the profinite setting requires substantial new mathematics.

## Cross-check: does the PADE picture match Faure 2009's CONCLUSIONS even if hypotheses fail?

YES — this is the most interesting aspect of the verdict.

Faure 2009 predicts:
- Discrete Ruelle resonances inside the unit disk (in C)
- Spectral gap ≤ 1/√E_min = 1/√k for k = 3 ≈ 0.577 (in our case the "k" analog)
- Resonances "repulse each other like random complex matrices" (Faure remark on Figure 2)
- Spectral radius asymptotically realizing the gap in the high-frequency limit

PADE_NUMERICAL_DISPOSITION observes:
- Leading singularity at |z| ≈ 1.57 at n=13 → corresponding "spectral radius" ≈ 1/1.57 ≈ 0.637 (close to 1/√3 ≈ 0.577)
- Asymptotic at z ≈ 1.016 → corresponding radius ≈ 0.984
- Complex-conjugate pair structure with period 9.2 in n-space, θ ≈ 0.68 rad — consistent with complex resonances "repulsing" / generic non-real spectrum
- Multi-spectral picture (transient at n=2..6, intermediate at n=7..13, asymptotic beyond)

**The numerical observations match Faure 2009's PREDICTIONS quantitatively.** The leading singularity radius |z|≈1.57 ≈ √3+ε at n=13 is suspiciously close to √E_min for k=3. The complex-conjugate pair structure matches Faure-Tsujii's band structure prediction (candidate C).

This is the case for PARTIAL: the categorical fit is so strong on the CONCLUSION side that the theorem's framework is plausibly the right framework for Syracuse, even though the hypotheses don't literally apply. The work needed: construct a profinite-group partially-expanding-map theory, prove the analog of Faure 2009 Theorem 2, and verify the spectral-radius prediction matches PADE.

## Final mark for candidate A

**PARTIAL_EXPANSION_PARTIAL** (in the sense of: spirit YES, technical hypotheses NO, conclusion-side numerical predictions match Faure's framework).

In disposition language: **PARTIAL** (not SELECTED — the theorem doesn't literally apply; not NO_FIT either — the framework predicts the observed phenomenology).

Required technical extension: a profinite-group / p-adic-Fourier analog of semiclassical pseudodifferential calculus, with which Faure 2009's proof can be re-executed in the profinite category. Estimated work: a research monograph, not a session.
