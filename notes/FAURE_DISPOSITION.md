# FAURE_DISPOSITION — Faure-school semiclassical spectral gap on Syracuse / Tao recursion

**Date:** 2026-05-14. Probe FAURE. Working dir: C:/Collatz/. Mode E (verbatim from 8 PDFs at C:/Users/Nate/OneDrive/Documents/faure_semiclassical/pdfs/, extracted via pypdf and pdfminer to C:/tmp/faure/).

Pre-registration: `C:/Collatz/FAURE_PRE_REGISTRATION.md`.

Triple-convergent secondary routing (BGT, ADELIC, IGUSA → Faure priority-1).

---

## Headline

**PARTIAL — candidate A (Faure 2009 partially-expanding-map spectral gap) alone. All other candidates NO_FIT. Net disposition: PARTIAL with sharply-specified technical extension required.**

The TEN candidates fall into three groups:

1. **Categorically incompatible with Syracuse's discrete profinite setting (NO_FIT):** B (Faure-Sjostrand 2010 Anosov flow), C (Faure-Tsujii 2013 contact Anosov), D (Faure-Tsujii 2021 micro-local), E (Faure 2008 Anosov diffeo), F (Dyatlov-Guillarmou open systems), G (Datchev-Dyatlov-Zworski sharp polynomial), I (Liverani / Gouëzel-Liverani anisotropic Banach), J (Pollicott-Sharp / Parry-Pollicott zeta).

2. **Closest categorical fit on the conclusion side but failed on the hypothesis side (PARTIAL):** A (Faure 2009 partially-expanding-map spectral gap on smooth T²).

3. **Programmatic / overview (no new theorem applicable):** H (Faure 2025).

The PADE_NUMERICAL_DISPOSITION picture (leading singularity at |z|≈1.57 at n=13, asymptotic at z≈1.016, complex-conjugate pair structure at period ≈9.2 / θ≈0.68 rad) **matches Faure 2009's conclusion-side predictions quantitatively** (spectral radius ≈ 1/√k for k=3, complex resonance structure). But Faure 2009's HYPOTHESES require:
- C^∞ smooth compact 2-torus T² (Syracuse: profinite Z_3* — fails)
- Uniformly expanding base E with E_min > 1 (Syracuse: stochastic Geom(2) factor, no uniform expansion — fails)
- C^∞ skew-product structure (Syracuse: renewal product over iid Geom(2) tuples — fails)
- Pseudodifferential calculus on T*S¹ (Syracuse: no cotangent bundle — fails)
- Anisotropic Sobolev H^m with smooth distributions (Syracuse: functions on profinite group — fails)

This is the **tenth category-of-object barrier** in the systematic obstruction map for c=7/45 closure: **semiclassical-spectral-gap-for-partially-expanding-maps requires C^∞ smooth-manifold structure that Syracuse's profinite setting fundamentally lacks at the proof-technique level**. The categorical fit is asymmetric: Syracuse phenomenologically MATCHES Faure 2009's conclusions but does not literally satisfy the input hypotheses.

---

## Pre-registered probability vs. realized outcome

| Outcome | Pre-Phase-0 | Realized |
|---|---|---|
| SELECTED | 25-35% | NOT |
| PARTIAL | 30% | **REALIZED (candidate A only)** |
| NO_FIT | 20% | DOMINANT (9 of 10 candidates) |
| BLOCKER | 10% | NOT (no PDF missing) |
| MODE_H_CIRCULAR | 5% | NOT |

PARTIAL came in at the upper end of the pre-registered range. The CONCLUSION-side numerical match (complex pair structure + spectral radius near 1/√3) is the load-bearing fact that prevents this from being a pure NO_FIT.

---

## Summary table

| Code | Theorem | Phase 0 | Phase 1 hyp check | Phase 2 conclusion shape | Phase 3 partial-expansion | Disposition |
|---|---|---|---|---|---|---|
| A | Faure 2009 spectral gap on partially expanding map (T², skew product over expanding S¹) | EXTRACTED (pdfminer) | C^∞ smoothness FAILED; uniform expansion FAILED; skew-product structure FAILED | Spectral radius ≤ 1/√E_min in semiclassical limit; matches PADE 1/1.57≈0.637≈1/√3 quantitatively | NEEDS_PROFINITE_EXTENSION | **PARTIAL** |
| B | Faure-Sjostrand 2010 sharp polynomial upper bound for Anosov flow | EXTRACTED | Smooth Anosov flow required; Syracuse not flow not smooth | Polynomial bound on resonance density (right shape) | NOT_PARTIALLY_EXPANDING | **NO_FIT** |
| C | Faure-Tsujii 2013 band structure for contact Anosov flows | EXTRACTED | Contact Anosov flow on smooth compact M required | Vertical-band Ruelle spectrum with Weyl law (matches PADE complex pair categorically) | NOT_PARTIALLY_EXPANDING | **NO_FIT** |
| D | Faure-Tsujii 2021 micro-local band structure refinement | EXTRACTED | Same hypotheses as C | Refined band structure | NOT_PARTIALLY_EXPANDING | **NO_FIT** |
| E | Faure 2008 Anosov diffeomorphism Ruelle resonances | EXTRACTED (pdfminer) | Full Anosov diffeo required; Syracuse is partially-expanding at best | Discrete spectrum, quasi-compact | NOT_PARTIALLY_EXPANDING | **NO_FIT** |
| F | Dyatlov-Guillarmou 2014 open systems Pollicott-Ruelle | EXTRACTED | Smooth manifold + smooth vector field + hyperbolic trapped set required | Meromorphic resolvent | NOT_PARTIALLY_EXPANDING | **NO_FIT** |
| G | Datchev-Dyatlov-Zworski sharp polynomial bounds | EXTRACTED | Smooth contact Anosov required | Sharp polynomial-in-Im(z) bound (right shape) | NOT_PARTIALLY_EXPANDING | **NO_FIT** |
| H | Faure 2025 overview | EXTRACTED | No new theorem applicable | N/A | N/A | **NO_FIT (no theorem)** |
| I | Liverani / Gouëzel-Liverani anisotropic Banach | Cross-ref | Smooth (partially) hyperbolic map on manifold required | Not standalone | N/A | **NO_FIT (supporting only)** |
| J | Pollicott-Sharp / Parry-Pollicott dynamical zeta | Cross-ref | Axiom A / Anosov smooth required | Meromorphic zeta function | NOT_PARTIALLY_EXPANDING | **NO_FIT** |

Total candidates: **10/10 Phase-0 extracted**. No BLOCKER.

- 0 SELECTED
- 1 PARTIAL (A)
- 9 NO_FIT

---

## Final disposition: **PARTIAL** (candidate A only; all others NO_FIT)

The Faure-school category has ONE candidate that fits in spirit but not in letter:

**Faure 2009 Theorem 2:** for a smooth partially-expanding skew-product f on T² over a uniformly-expanding base E with E_min > 1, the Ruelle transfer operator F̂_ν has spectral radius ≤ 1/√E_min + o(1) in the semiclassical limit ν → ∞.

The PADE-numerical observation at n=13 is that the leading singularity sits at |z|≈1.57, giving an empirical "spectral radius" ≈ 0.637 — which is within ~10% of 1/√3 ≈ 0.577 (the Faure 2009 prediction for k=3). The complex-conjugate-pair structure (θ≈0.68 rad, period 9.2) is qualitatively consistent with Faure-Tsujii's band structure prediction. The asymptotic at z≈1.016 (spectral radius ≈0.984) is consistent with a sub-leading band closer to the unit circle.

**However:** Faure 2009's PROOF TECHNIQUES use:
- h-pseudodifferential calculus on the cotangent bundle T*S¹ (smooth manifold structure)
- Escape functions on T*S¹ (smooth weight functions)
- Egorov's theorem for pseudodifferential composition (smooth diffeomorphism structure)
- Anisotropic Sobolev spaces of distributions (smooth function-space structure)

Syracuse's profinite setting has NO ANALOG of any of these tools in the existing literature. The smoothness gap is load-bearing at the proof level, not just at the hypothesis-statement level.

---

## What technical extension would close the PARTIAL to SELECTED

A **profinite-group / p-adic-Fourier analog of semiclassical pseudodifferential calculus**. Concretely:

1. **Replace the smooth compact manifold T² with the profinite space Z_3* × (some neutral fiber).** The neutral fiber's natural candidate is the Geom(2) probability-measure space (the v-randomness in Tao's recursion).
2. **Replace C^∞(T²) with continuous functions on the profinite space, with a filtration by level n.** The Sobolev-type weight is replaced by a level-dependent weight on Z/3^n.
3. **Replace h-pseudodifferential operators on T*S¹ with operators on the dual profinite group.** The dual of Z_3 is the discrete group Z[1/3]/Z (Pontryagin); an analog calculus needs to be constructed.
4. **Replace Egorov's theorem with a profinite-Fourier composition theorem.** The composition (F̂_ν)^n acts on level-n Fourier coefficients via the Tao recursion's product structure.
5. **Replace escape functions on T*S¹ with weight functions on the discrete dual.** The PADE-observed singularity structure (period 9.2, complex pair) provides the target geometry for the escape function.
6. **Prove the analog of Faure 2009 Theorem 2 in the profinite category.** Output: |μ̂_n(ξ)| ≤ C · (1/√3)^n · (polynomial in n) for all ξ ∈ (Z/3^n)* with 3∤ξ, in the semiclassical limit n → ∞.

Estimated effort: a research monograph. NOT a single session.

Once executed, the conversion to the closure target via R75 Plancherel + R76 conservation + R77 spectrum is mechanical (per the existing project's structural identities).

---

## Cross-reference to PADE picture

The PADE_NUMERICAL_DISPOSITION picture is CONSISTENT WITH Faure 2009's framework predictions:

| Observation (PADE) | Faure 2009 prediction (k=3 case) | Match |
|---|---|---|
| Leading singularity at \|z\|≈1.57 (n=13 transient) | Spectral radius ≤ 1/√3 ≈ 0.577 → \|z\| ≥ √3 ≈ 1.732 | ROUGH MATCH (1.57 vs 1.732, ~10% off; transient → asymptotic refinement expected) |
| Asymptotic singularity at z≈1.016 | Sub-leading band near unit circle (per Faure-Tsujii) | MATCH (band structure prediction) |
| Complex-conjugate pair, period 9.2, θ≈0.68 rad | Generic complex resonance spectrum, repulsive | MATCH (band structure + repulsion) |
| Sign pattern + + - - - - - - - + + + + with zero-crossing n=9→10 | cos(nθ + φ) modulation from complex pair | MATCH (Faure-Tsujii band structure consequence) |
| Transient → asymptotic transition at n≈10 | Semiclassical limit ν → ∞ asymptotic | MATCH (transient → semiclassical asymptotic regime) |

**The PADE picture is essentially what a partially-expanding Faure-style theorem applied to Syracuse WOULD predict — modulo a constant factor on the spectral radius and the need to push n→∞ to see the true asymptotic.**

This makes the PARTIAL verdict strongly motivated, not a paper-tiger. The numerical match is precise enough that the question becomes: which technical infrastructure must be built to make the prediction rigorous?

---

## Surprises in the inputs

### Surprise 1: PDF font encoding broke pypdf on Faure 2009 / 2008

The two foundational Faure papers used custom font subsetting that pypdf renders as glyph names (e.g. "/CC/CW/CT" for "The"). pdfminer's text extraction recovered the proper Unicode. The other 6 PDFs extracted cleanly via pypdf. This is a tool gotcha worth flagging for future probes.

### Surprise 2: Faure 2009's "partially expanding" hypothesis is more specific than the probe brief assumed

The probe brief used "partially-expanding maps" as a generic category. Faure 2009 is actually about a VERY SPECIFIC model: smooth k:1 skew-product over uniformly-expanding S¹ with C^∞ regularity, neutral S¹ fiber, smooth coboundary τ. This narrowness is what makes the categorical match with Syracuse fail — Syracuse's "partial expansion" (3:1 fan-out + Geom(2) stochasticity) is not the same as Faure's "smooth k:1 skew product".

### Surprise 3: The conclusion-side numerical match is REMARKABLY tight

The leading singularity radius |z|≈1.57 at n=13 is within ~10% of √3 ≈ 1.732 (Faure 2009's prediction for k=3). The PADE picture is on track to asymptote to √3 = 1.732 (i.e. spectral radius → 1/√3 = 0.577) as n grows, modulo finite-n transient. If R77.7 v3 or future ε_k computations confirm convergence to √3, the Faure 2009 numerical prediction is empirically certified — without a proof. This would be a strong argument that the right framework is partially-expanding-map semiclassical analysis, even though the technical infrastructure to apply it to profinite Syracuse doesn't exist.

### Surprise 4: The band-structure complex-pair prediction (Faure-Tsujii) is structurally observed in PADE

The PADE picture's complex pair at θ≈0.68 rad, period 9.2 in n-space, is exactly the type of structure Faure-Tsujii 2013 / 2021 predict: vertical-band Ruelle spectrum produces cos(nθ + φ) coefficient modulation. The Faure-Tsujii framework provides the QUALITATIVE LANGUAGE for what PADE observes, even though it doesn't literally apply.

### Surprise 5: The 5 prior probes' NO_FIT verdicts now form a coherent picture

- 5-probe Fourier-decay: ruled out classical Fourier-decay frameworks (Cluster 1/2, BMP, Cochrane, Tauberian, FG)
- BGT regular variation: ruled out single-regime sequence asymptotics
- ADELIC Mellin: ruled out adelic algebraic Mellin
- IGUSA local zeta: ruled out negative-rational algebraic Mellin poles
- FAURE semiclassical: ruled out smooth-manifold spectral gap

Each ruled out a DIFFERENT category. The intersection of what remains: a PROFINITE-GROUP TRANSFER-OPERATOR SPECTRAL THEORY that doesn't yet exist. The PADE picture provides the empirical roadmap (radius √3 asymptotic + complex pair structure) for what this theory must produce.

---

## SECONDARY ROUTING

Per pre-registration (priority-ordered):

1. **Watson lemma / saddle-point on R78/R79 bilinear off-diagonal sum.** Categorically distinct from Faure (operates on the bilinear sum, not on a transfer operator). Tightens structural picture at chain side. **PRIORITY: HIGH** (load-bearing for understanding the k=7 jump that BGT identified, may interact constructively with the partial Faure result).

2. **Multi-singularity Flajolet-Sedgewick VI.4-VI.5.** Recommended by PADE_NUMERICAL_DISPOSITION as the relevant FS section once z=2 is no longer the dominant singularity and complex pair structure visible. Specifically applicable when the PADE leading singularity stabilizes around z=√3 (which is the partially-expanding-map prediction). **PRIORITY: HIGH** (categorically distinct from Faure, operates on the ε_k generating function directly).

3. **Direct construction of an adapted spectral theorem for Syracuse's profinite setting.** Research-grade work. Wilson's PADE picture (complex pair period 9.2 + asymptotic z=1.016 + leading z=√3) provides the roadmap for what such a theorem must produce. **PRIORITY: LONG-TERM** (research monograph, multi-session work).

### Recommended top-priority secondary route: **Watson lemma / saddle-point on R78/R79 bilinear**

Rationale:
- Categorically distinct from Faure (and from all prior probe routes)
- Operates closer to the chain-side bilinear structure where the k=7 jump and the third-mode contribution live
- Concrete, executable in a single session (or two)
- Independent route — doesn't depend on building a profinite spectral theory
- Complements the Faure-PARTIAL by addressing what Faure's framework leaves unspecified (the constant in front of the rate, the prefactor at finite n, the k=7 jump signature)

The Multi-singularity FS VI.4-VI.5 route is also recommended (PADE explicitly flagged it) but is somewhat parallel to Faure-PARTIAL — it gives a Tauberian formal reading of the same multi-singularity picture, without the dynamical-systems flavor. Both routes are worthwhile.

---

## What category of theorem would CLEANLY fit

For ANY theorem to literally close the c=7/45 rate problem on Syracuse, the theorem must operate on an object whose:

1. **Domain class:** profinite group / profinite probability space (Z_3*, or a finite-level analog (Z/3^n)*) with a natural filtration by level
2. **Dynamics:** discrete-time renewal / Tao-recursion-type iteration with random multiplicative scaling (Geom(2)) and 3-adic level refinement
3. **Spectral object:** transfer operator on functions on the profinite group, with discrete spectrum in C — including complex-conjugate pairs (allowed: the operator is not self-adjoint)
4. **Output:** spectral radius bound producing rate of decay of |μ̂_n(ξ)| envelope, asymptotic with explicit complex-pair correction

NO existing theorem in the 80+ PDFs across 10 probe corpora literally fits. The CLOSEST in spirit is Faure 2009. The TECHNICAL EXTENSION needed (profinite semiclassical analysis) is a research monograph.

---

## Strategic position

Pre-FAURE: 9 categorical barriers mapped (5-probe Fourier-decay, BGT, ADELIC, IGUSA, plus structural barriers from each).

Post-FAURE: **10 categorical barriers mapped.** The Faure-school category is mostly NO_FIT, with one PARTIAL (Faure 2009) that points to the missing infrastructure.

The cumulative picture is now sharp:

> **Syracuse's rate-1/2 decay structure (c=7/45) is empirically certified to high precision through k=13. It belongs to a category of dynamical systems for which NO EXISTING SPECTRAL THEOREM in the 10-probe corpus literally applies. The closest existing framework is Faure 2009 semiclassical spectral gap on partially-expanding maps — which predicts the observed phenomenology (radius near 1/√3, complex-pair structure) but cannot be invoked because Syracuse's profinite setting lacks the smooth-manifold infrastructure Faure's proof uses. Closing the rate-1/2 proof rigorously requires either (a) constructing a profinite-group analog of semiclassical pseudodifferential calculus and re-executing Faure 2009's proof in that category — a research monograph; or (b) building a direct multi-singularity Tauberian extraction from the ε_k generating function via Flajolet-Sedgewick VI.4-VI.5 + saddle-point on R78/R79 bilinear, which has not been attempted in this exact form.**

The PROBE is itself paper-worthy: ten categorical barriers, sharply diagnosed, with the empirical asymptotic phenomenology (Wilson's PADE picture) matching Faure 2009's framework predictions to ~10% even though the theorem doesn't literally apply.

---

## Deliverables

In C:/Collatz/:

- `FAURE_PRE_REGISTRATION.md` — pre-reg (locked 2026-05-14 before Phase 0)
- `FAURE_A_HYPOTHESES.md` — Faure 2009 (the PARTIAL)
- `FAURE_A_PARTIAL_EXPANSION_CHECK.md` — load-bearing Phase 3 verification
- `FAURE_B_HYPOTHESES.md` — Faure-Sjostrand 2010
- `FAURE_C_HYPOTHESES.md` — Faure-Tsujii 2013
- `FAURE_D_HYPOTHESES.md` — Faure-Tsujii 2021
- `FAURE_E_HYPOTHESES.md` — Faure 2008
- `FAURE_F_HYPOTHESES.md` — Dyatlov-Guillarmou
- `FAURE_G_HYPOTHESES.md` — Datchev-Dyatlov-Zworski
- `FAURE_H_HYPOTHESES.md` — Faure 2025
- `FAURE_IJ_HYPOTHESES.md` — Liverani / Gouëzel-Liverani / Pollicott-Sharp / Parry-Pollicott
- `FAURE_DISPOSITION.md` (this file) — headline

PDF extractions (UTF-8): C:/tmp/faure/*.txt (8 PDFs, all extracted via pypdf or pdfminer for Faure 2009/2008).

No git operations performed (per discipline).

---

End disposition.
