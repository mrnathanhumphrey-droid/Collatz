# C2_SYRACUSE_CANDIDATES — candidate (G, H, L) schemes for Syracuse μ_n

**Date:** 2026-05-12. Cluster 2 cut-and-project probe, Phase 2.

## Setup: Syracuse μ_n's structural origin

Per `result_78_FINAL.md`, `T_LEAD_CORRECTED_DISPOSITION.md`, and the seven-probe trail:

- μ_n is the stationary distribution of the **Syracuse Markov chain on (ℤ/3^n)\***, the **units** mod 3^n.
- |supp(μ_n)| = φ(3^n) = 2·3^{n−1}.
- As n → ∞, the supports nest into ℤ_3\* = ℤ_3 \ 3ℤ_3 (the 3-adic units), with limit measure μ on ℤ_3\*.
- Weights μ_n(x) for x ∈ (ℤ/3^n)\* are NOT uniform — they encode the Markov chain's stationary structure. Numerically they are mass-concentrated on the "+" class structure (from R76/R77; the F̂ support {a ≡ 1 mod 3} match in R78 §"Crucial observation" is at the **dual-side** of the same 3-adic class decomposition).
- Closure target object is |μ̂_n(ξ)|² for ξ ∈ ℤ/3^n; i.e., G_target = ℤ/3^n with dual ℤ/3^n, lifting to G_∞ = ℤ_3 with dual Ĝ_∞ = ℚ_3/ℤ_3 = ℤ[1/3]/ℤ.

So Syracuse μ is a measure on ℤ_3 (not a *point set* in the BMP sense), with weights given by a Markov-chain stationary distribution.

This already places Syracuse in the **weighted-comb** regime (Theorem 5.4 / Theorem 5.9 of 1512.00912) rather than the strict regular-model-set regime (Theorem 5.9 of 1606.08831 for 1_W).

---

## Candidate A — Direct ℤ_3 with adelic complement

**G = ℤ_3** (compact, hence σ-compact, abelian LCA — qualifies for 1512.00912's σ-compact-G version).

**H = ∏_{p ≠ 3} ℤ_p** (Tychonoff product of compact abelian groups; compactly generated since each ℤ_p is compact and the product over infinitely many is compactly generated as an LCA group).

**L = ?** Candidates:
- L_1 = ℤ embedded diagonally: n ↦ (n mod 3^∞, (n mod p^∞)_{p≠3}). But ℤ_3 has no "n mod 3^∞" in a finitary sense — rather n embeds into ℤ_3 via n ↦ n (its 3-adic value). So L_1 = {(n, (n)_{p≠3}) : n ∈ ℤ} as a subset of ℤ_3 × ∏_{p≠3} ℤ_p.
- Is L_1 a lattice (discrete + cocompact) in ℤ_3 × ∏_{p≠3} ℤ_p?
  - **DISCRETENESS FAILS.** ℤ is **dense** in ℤ_3 (as ℤ_3 is the 3-adic closure of ℤ). So π_{ℤ_3}(L_1) = ℤ is dense in ℤ_3. This means L_1 itself is dense in ℤ_3 × ∏_{p≠3} ℤ_p, not discrete.
  - More fundamentally: ℤ_3 × ∏_{p≠3} ℤ_p ≅ Ẑ (profinite integers) ≅ ∏_p ℤ_p, and ℤ embeds densely in Ẑ. No discrete subgroup of Ẑ is non-trivial (Ẑ has no non-trivial discrete subgroups since it's profinite/compact and totally disconnected with no isolated points outside 0).

**Verdict A: FAILS at "L is a lattice." π_G is injective on L_1 vacuously (since L_1 maps injectively to ℤ_3), but L_1 is not discrete in ℤ_3 × ∏_{p≠3} ℤ_p. Candidate A is structurally NOT a cut-and-project scheme.**

The intuition that fails: in BMP, the lattice is **ℤ embedded in ℝ × ∏'_p ℚ_p**, where the Archimedean ℝ factor is what makes ℤ discrete (since ℤ is discrete in ℝ). Strip ℝ off and replace by ℤ_3, and discreteness evaporates because ℤ is dense in ℤ_3.

---

## Candidate B — Full BMP-style adelic, G = ℝ (super-singular framing)

**G = ℝ** (real line, σ-compact LCA — matches BMP).

**H = ∏'_p ℚ_p** (rational adeles' finite part: restricted direct product, with all-but-finitely-many components in ℤ_p).

**L = ℚ** embedded diagonally in ℝ × ∏'_p ℚ_p (BMP's setup).

π_G(L) = ℚ ⊂ ℝ dense, π_H(L) = ℚ ⊂ ∏'_p ℚ_p dense (Strong Approximation Theorem). L is a discrete cocompact subgroup of ℝ × ∏'_p ℚ_p (fundamental domain [0,1] × ∏_p ℤ_p has volume 1). So **this IS a cut-and-project scheme.**

**Window for Syracuse:** The Syracuse Markov chain lives on ℤ_3\* asymptotically. In BMP coordinates this means the window in H = ∏'_p ℚ_p must:
- project to ℤ_3\* in the p=3 component
- project to ℤ_p in the p ≠ 3 components (to recover integer constraint)
- have Markov-chain-stationary-weighted mass within the 3-adic component (not uniform indicator)

For an unweighted Syracuse-support set Λ ⊂ ℤ:
**W_Syracuse = (ℤ_3\*) × ∏_{p≠3} ℤ_p ⊂ ∏'_p ℚ_p.**

Then Λ(W_Syracuse) = π_G(L ∩ (G × W_Syracuse)) = {n ∈ ℚ : n ∈ ℤ_3\* AND n ∈ ℤ_p for all p ≠ 3} = {n ∈ ℤ : gcd(n, 3) = 1} = the **3-coprime integers in ℤ**.

**Note: π_G(L ∩ (G × W_Syracuse)) is NOT (ℤ/3^n)\* for fixed n. It is the full set of 3-coprime integers — i.e., the asymptotic n=∞ support.**

This is a genuine cut-and-project encoding for the **set of 3-coprime integers in ℤ**, which is dens(L)·θ_H(W_Syracuse) = 1·(2/3) = 2/3 fraction of integers.

**But this is the unweighted SET. Syracuse μ_n has non-uniform weights given by the Markov stationary distribution.** That requires Candidate B-weighted (next).

---

## Candidate B-weighted — Candidate B with weighted comb

Take (G, H, L) as in Candidate B, but instead of 1_W for the indicator of the 3-coprime support, use a weight function:

**h: ∏'_p ℚ_p → ℂ**, h(y) = μ(y_3) · ∏_{p≠3} 1_{ℤ_p}(y_p),

where μ on ℤ_3\* is the limiting Syracuse stationary distribution and y_3 is the p=3 component of y.

Then ω_h = Σ_{(x,y)∈L} h(y) δ_x is the weighted Dirac comb on ℝ.

**Theorem 5.4** would deliver ω̂_h = dens(L) · ω_{ĥ}.

But Theorem 5.4 requires **h ∈ P_K(H)**: positive-definite, continuous, COMPACTLY SUPPORTED. h built from μ on ℤ_3\* — is this in P_K?
- Continuous on ∏'_p ℚ_p (with restricted product topology): need μ continuous on ℤ_3\*. Syracuse μ is a measure, not a function; the **density of μ** w.r.t. Haar on ℤ_3 might not be continuous (Markov-chain stationary distributions on ℤ_3 are typically singular continuous or atomic at level n, and continuity on ℤ_3 is a deep question — this is precisely the polynomial-in-A Fourier-decay question!).
- Compactly supported: yes, on (ℤ_3\*) × ∏_{p≠3} ℤ_p which is compact.
- Positive-definite: positive-definite on the LCA group ∏'_p ℚ_p means ĥ is a positive measure on Ĥ. This is a STRONG condition on μ — not guaranteed.

**So Theorem 5.4 doesn't fire directly for h = μ. The hypothesis P_K(H) is exactly the condition we're trying to verify for Syracuse μ.**

This is a circularity: the polynomial-in-A bound on |μ̂_n| would imply h's continuity/regularity to apply Theorem 5.4, and Theorem 5.4 then gives back |1̂_W| structure. Without independent verification of h ∈ P_K, Theorem 5.4 can't be applied.

**Verdict B-weighted: PARTIAL — encoding exists structurally but the load-bearing hypothesis h ∈ P_K(H) on the weight function is itself unverified for Syracuse μ. The target-object trap (Mode H) is closing here: the property we'd derive (μ̂_n decay) is the property we'd need to assume (h positive-definite hence Fourier-transformable as a measure).**

---

## Candidate B-unweighted (3-coprime set in ℤ)

Returning to the unweighted version of B: the set of 3-coprime integers as a subset of ℤ.

This is a **weak model set with maximal density**, exactly analogous to BMP's F_k for k=1 in the p=3 component only. Theorem 5.9 of 1512.00912 fires:
- maximal density: lim |{n ∈ [−N, N] : gcd(n, 3) = 1}|/(2N) = 2/3 = dens(L)·θ_H(W_Syracuse). YES, by uniform-distribution / density of 3-coprime integers.
- W = (ℤ_3\*) × ∏_{p≠3} ℤ_p: relatively compact (compact, in fact), measurable.

So the **set of 3-coprime integers IS a weak model set with maximal density**, and its diffraction is computable as γ̂ = dens(L)² · ω_{|1̂_W|²}.

**But this is NOT the Syracuse stationary measure μ_n.** It's the unweighted support set. The Syracuse μ_n is a *probability measure on the support* with Markov-chain-derived weights. The diffraction of the unweighted comb δ_Λ gives no information about |μ̂_n(ξ)| for the weighted measure.

The closure question concerns Σ_x μ_n(x) χ(x), a character sum **weighted** by μ_n. The unweighted diffraction tells us how the SET of 3-coprime integers diffracts (gives a known computation: rationals with cubefree denominators, factoring through Riemann ζ), not how the weighted Markov measure diffracts.

**Verdict B-unweighted: VALID as a cut-and-project encoding, but the OBJECT is wrong (Mode H, target-object selection trap). The unweighted set is not the polynomial-in-A target.**

---

## Candidate C — Local-only G = H = ℤ_3

**G = ℤ_3, H = ℤ_3.**

For G × H = ℤ_3 × ℤ_3, the natural diagonal L = ℤ_3 (image: {(x, x) : x ∈ ℤ_3}) is **not discrete** in ℤ_3 × ℤ_3 (it's a closed subgroup of full Haar measure, not discrete). To get a discrete L, we'd need L = (3^k ℤ_3) × (3^j ℤ_3) for some shifts, but these are open subgroups, not lattices (not relatively dense as discrete sets in ℤ_3 × ℤ_3 — they're either full or empty in measure terms).

The fundamental issue: **ℤ_3 has no non-trivial lattice (cocompact discrete subgroup)** because ℤ_3 is profinite, totally disconnected, and any discrete subgroup is trivial (consists of isolated points; in a profinite group, isolated points have no accumulation, so discrete subgroups are finite; finite subgroups of ℤ_3 are trivial since ℤ_3 is torsion-free).

**Verdict C: FAILS. No non-trivial lattice in ℤ_3 × ℤ_3. The local-only candidate has no valid (G, H, L).**

---

## Candidate D — Sidestep: 3-adic-only window in BMP scheme (refinement of B)

Same (G, H, L) = (ℝ, ∏'_p ℚ_p, ℚ) as Candidate B, but choose the window to encode finer 3-adic structure:

**W_n = (W_3,n) × ∏_{p≠3} ℤ_p**, where W_3,n ⊂ ℤ_3 is the set of x ∈ ℤ_3 with x mod 3^n ∈ supp(μ_n) (i.e., the lift to ℤ_3 of the level-n Syracuse support).

For each n, W_3,n is open and compact (clopen) in ℤ_3, hence W_n is regular in ∏'_p ℚ_p. Λ(W_n) is the set of integers whose mod-3^n residue is in the Syracuse support of level n.

**Verdict D: works as a finite-level approximation but DOESN'T deliver μ_n's weighted measure — it delivers the unweighted level-n support set, which is just a coset union, not the stationary distribution.**

Stacking the W_n into a tower as n → ∞ converges to the full support (3-coprimes) and again loses the weights.

---

## Summary of candidates

| Candidate | G | H | L | Verdict | Failure mode |
|---|---|---|---|---|---|
| A | ℤ_3 | ∏_{p≠3} ℤ_p | ℤ diagonal | FAILS | L not discrete (ℤ dense in Ẑ) |
| B | ℝ | ∏'_p ℚ_p | ℚ diagonal | VALID scheme but wrong object | Gives unweighted 3-coprime set, not μ_n |
| B-weighted | same | same | same + h = μ | PARTIAL (circular) | h ∈ P_K(H) is exactly what we'd need to derive |
| C | ℤ_3 | ℤ_3 | — | FAILS | No non-trivial lattice in profinite × profinite |
| D | ℝ | ∏'_p ℚ_p | ℚ | VALID at finite level | Window encodes set, not weights |

---

## Disposition for Phase 2

**The natural cut-and-project encoding of Syracuse μ_n has two distinguishable layers:**

1. **Support layer.** The set of 3-coprime integers ⋏(W_Syracuse) IS a weak model set with maximal density in the BMP scheme (G=ℝ, H=∏'_p ℚ_p, L=ℚ). Theorem 5.9 fires for the **unweighted support** but tells us about the diffraction of δ_Λ = Σ_{x ∈ 3-coprime} δ_x, which has known form (ω_{|1̂_W|²} on rationals with cubefree denominator).

2. **Weight layer.** The Markov-chain weights μ_n on the support require a **weighted Dirac comb** ω_h with h built from μ on ℤ_3\*. The applicable theorem is 5.4 (weighted PSF analog), but it requires h ∈ P_K(H): positive-definite continuous compactly supported. This is a strong analytic condition on μ that is NOT a free consequence of μ being a Markov stationary measure.

The polynomial-in-A bound on |μ̂_n(ξ)| is a STATEMENT about μ_n's smoothness/decay at the Fourier-analytic level, which is exactly what's needed to verify h ∈ P_K. So **the framework requires the conclusion to set up its hypothesis**: a Mode H target-object circularity.

The best candidate, B, encodes the SUPPORT correctly but not the WEIGHTS. The weighted refinement (B-weighted) is structurally PARTIAL — the encoding pieces (G, H, L, candidate h) all exist, but the load-bearing hypothesis h ∈ P_K is the target-object property itself.

**Verdict: H_C2_ENCODING_PARTIAL (with Mode H target-object trap).** Candidate B is the BMP-canonical encoding; it captures Syracuse's support as a weak model set with maximal density, but the polynomial-in-A bound concerns the weights, which require Theorem 5.4 weighted-comb form whose hypothesis P_K(H) on h is circular with what we're trying to prove.
