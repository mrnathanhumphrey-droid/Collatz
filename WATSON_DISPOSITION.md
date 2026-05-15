# WATSON_DISPOSITION — Watson lemma / saddle-point on R78/R79 bilinear off-diagonal sum

**Date:** 2026-05-14. Probe WATSON. Working dir: C:/Collatz. Mode E (technique statements verbatim from C:/Users/Nate/OneDrive/Documents/watson_saddle_point/pdfs/ via pypdf; extracted to C:/tmp/watson/).

Pre-registration: `C:/Collatz/WATSON_PRE_REGISTRATION.md`. Locked before Phase 2 execution.

Secondary routing source: FAURE_DISPOSITION priority-1.

---

## Headline

**PARTIAL — saddle-point applies exactly at r=3 (R78.6 confirmed) and Darboux/multi-saddle on f(z) matches PADE singularity FORM qualitatively, but the technique does NOT supply the rigorous asymptotic RATE (the empirical κ=0.522 inter-a cancellation), nor does it close the conversion to |μ̂_n(ξ)| bound. Specific gap identified: rigorous transfer-operator analyticity in the profinite setting — IDENTICAL gap to FAURE_DISPOSITION.**

The probe was the strongest-motivated single-session route in the entire arc (per FAURE_DISPOSITION secondary-routing priority). Outcome lands at the upper bound of the pre-registered PARTIAL probability range — concrete match against PADE on FORM, no match on RATE.

---

## Pre-registered probability vs. realized outcome

| Outcome | Pre-Phase-2 | Realized |
|---|---|---|
| SELECTED | 25-35% | NOT |
| PARTIAL | 35% | **REALIZED** |
| NO_FIT | 15% | NOT |
| BLOCKER | 10% | NOT (sum form extracted in WATSON_R78_SUM_FORM.md) |
| MODE_H_CIRCULAR | 5% | NOT |

---

## Summary table per technique

| Code | Technique | Hypothesis match | Conclusion shape | Conversion to closure | Disposition |
|---|---|---|---|---|---|
| A | Watson's lemma (real Laplace) | FAILS (no Laplace structure in bilinear sum) | Standard Σ Γ a_n x^{-(α+βn+1)} not applicable | — | **NO_FIT** |
| B | Laplace's method | FAILS (no exponentially-weighted integrand at level r) | — | — | **NO_FIT** |
| C | Saddle-point / steepest descent | APPLICABLE (R78.6 already uses Cochrane Prop 4 saddle at r=3) | Gaussian neighborhood of saddle s*∈Z/p; exact at r=3 | r ≥ 4 needs Hensel | **PARTIAL** (used at r=3 RIGOROUS; r ≥ 4 OPEN) |
| D | Stationary phase | APPLICABLE on the discrete sum (= Poisson summation analog, already used in PATH2_BILINEAR Inner-Plancherel) | Cancels c_2 variable cleanly | Gives |T_p| ≤ 2N at r=3 | **PARTIAL** (used; bound is κ=1, doesn't reach κ=0.522) |
| E | Mellin-Barnes contour | APPLICABLE in principle for ε_k ↔ f(z) singularity | Same as Darboux | Requires analyticity of f(z) (NOT supplied) | **PARTIAL** (form match; rigor depends on transfer-op theory) |
| F | Multi-saddle / complex pair | FORM MATCH (PADE complex pair at θ≈0.68 rad period 9.2) | cos(kθ+φ) modulation in ε_k | Singularity locations require transfer-op input | **PARTIAL** (form match; rate not rigorous from saddle alone) |
| G | p-adic saddle-point (Cochrane Prop 4) | EXACT at r=3 (R78.6) | ψ_lead = e_q(P_a(s*)) machine-precision | Hensel-lifted r ≥ 4 OPEN | **PARTIAL** (exact at r=3; family extension open) |
| H | Watson + Borel resummation | NOT NEEDED (series doesn't diverge) | — | — | **NO_FIT** |
| I | Uniform/coalescing (Airy-type) | NOT APPLICABLE (PADE singularities not coalescing) | — | — | **NO_FIT** |
| J | Darboux's method | APPLICABLE for ε_k coefficient asymptotic | a_n ~ ρ^{-n} · n^{-α} · cos(nθ+φ) FORM matches PADE | Singularity location requires rigor | **PARTIAL** (form match; rigor depends on transfer-op theory) |

**Net by category:** C+D+E+F+G+J all PARTIAL (form match, rate gap). A+B+H+I NO_FIT.

---

## Final disposition: **PARTIAL**

The PARTIAL is sharp:

> **What saddle-point/Darboux supplies:** The functional form of the ε_k asymptotic (real exponential decay + complex-conjugate-pair cosine modulation) matches PADE qualitative structure. The saddle-point at the level-r side (R78.6) is exact at r=3 and verified. The Inner-Plancherel reduction (PATH2_BILINEAR Attempt G+) gives a rigorous κ ≤ 1 bound at r=3.
>
> **What saddle-point/Darboux does NOT supply:**
> 1. Rigorous identification of the singularity ρ_∞ of f(z) = Σ ε_k z^k. PADE estimates ρ ≈ 1.016 (asymptotic slow-mode) but this is finite-data inference; the saddle-point technique operates on the integrand, not on the analyticity properties of f.
> 2. Inter-a phase cancellation that produces R79b's empirical κ = 0.522 at the level-r side. Classical saddle-point + Plancherel + Cauchy-Schwarz give only κ = 1.
> 3. Hensel-lifted closed form of ψ_true(a) at r ≥ 4. Phase deviation up to 160° at r=5 (verified in Phase 2).
> 4. Rigorous conversion to |μ̂_n(ξ)| bound. Requires transfer-operator analyticity — the same FAURE_DISPOSITION-identified missing infrastructure.

**The gap matches FAURE_DISPOSITION exactly:** rigorous spectral / analytic theory of the transfer operator T_M (or its dual on f(z)) in the profinite setting. Saddle-point CONFIRMS the form predicted by Faure 2009 (cos-modulated exponential decay), but cannot itself establish the rate rigorously without the smooth-manifold infrastructure that profinite Syracuse lacks.

---

## What technical step would close PARTIAL to SELECTED

The same step as in FAURE_DISPOSITION: a **profinite analytic transfer-operator theory** that:

1. Establishes f(z) = Σ ε_k z^k is meromorphic in some annulus around |z| = ρ_∞
2. Identifies the rightmost singularity (real or complex pair) with explicit ρ_∞
3. Provides rigorous error bars on the Darboux coefficient asymptotic

OR an alternative route through:

4. Direct band-l¹ analysis of ĥ_{r,ℓ} (R79's identified open path)
5. Bourgain-Konyagin sum-product on the multiplicative subgroup ⟨4⟩

None of these are single-session work. Saddle-point + Darboux SUPPLY THE TARGET FORM; they don't close it.

---

## Cross-reference to PADE picture

The Phase 2 numerical computation reproduces (and partially adjudicates) the PADE picture:

| Element | PADE | WATSON Phase 2 | Status |
|---|---|---|---|
| Hadamard radius at k=13 | 1.57 | 1.565 | EXACT MATCH (same data) |
| Transient → asymptotic transition | n ≈ 10 | r ≈ 6..8 (T_p side); k ≈ 10..13 (ε_k side) | CONSISTENT |
| Leading singularity at k=13 | \|z\| ≈ 1.57 (complex pair?) | Hadamard trend says decreasing; multi-saddle fit ambiguous | CONSISTENT (within data) |
| Period 9.2 in n-space | from sign-pattern analysis | Free fits prefer period ~30 or ~4; period 9.2 RSS 4.6× worse | INCONSISTENT (one sign change insufficient to constrain) |
| Slow-mode ρ ≈ 1.016 | from STATE.md k-space analysis | NOT VISIBLE in k=2..13 data; transient regime dominates | NEEDS k ≥ 20 DATA |
| Faure prediction √3 | Faure 2009 Thm 2, k=3 | Held-√3 fit RSS = 1.75e-5, comparable to held-1.57 | CONSISTENT within margin |

The PADE picture is not strongly contradicted by WATSON Phase 2; it's also not strongly confirmed
by independent saddle-point inference. The two pictures are AGNOSTICALLY CONSISTENT.

---

## Cross-reference to Faure 2009 prediction

Faure 2009 predicts spectral radius 1/√3 ≈ 0.577 in semiclassical limit ν → ∞. Watson Phase 2:

- Hadamard tail at k=13 gives 1/ρ_13 = 0.639 (10% above 1/√3)
- Held-√3 multi-saddle fit RSS = 1.75e-5
- Held-1.57 multi-saddle fit RSS = 1.71e-5

**The two predictions (Faure asymptotic 1/√3, PADE transient 1.57) are within 10% of each other
and the data is consistent with both at k=10..13.** Faure's asymptotic prediction (semiclassical
limit) would manifest as ρ_k → √3 as k → ∞; PADE's transient (n=13) at 1.57 may be the same
phenomenon at finite k.

**Decisive test:** compute ε_k=14, 15, 16, ...; if ρ_k stabilizes at √3, Faure's prediction is
empirically confirmed.

---

## Triple-PARTIAL synthesis

Combined with BGT PARTIAL and FAURE PARTIAL:

- **BGT PARTIAL** (sequential RV / Kendall K2 fires in plateau k=2..6, fails at k=7 jump):
  identifies plateau-to-transient transition at k=7 as the structural boundary.
- **FAURE PARTIAL** (Faure 2009 framework predicts cos-modulated exponential decay; smooth-manifold
  proof-machinery doesn't transfer to profinite): identifies the spectral-gap framework as the
  right CATEGORY.
- **WATSON PARTIAL (this probe)**: saddle-point supplies form, doesn't supply rate. Identifies
  transfer-operator analyticity as the missing rigor.

All three PARTIALs converge on the SAME structural gap: **need a profinite analog of the smooth
semiclassical transfer-operator theory**. The PADE phenomenology + Faure framework + saddle/Darboux
shape all match QUALITATIVELY; the RIGOROUS execution requires the missing infrastructure.

This is the **eleventh barrier** in the obstruction map for c=7/45 closure:

> The saddle-point / Darboux machinery is universal enough to apply, and predicts the right
> functional form, but does not deliver rigorous asymptotic rates without transfer-operator
> analytic input. The combined BGT + FAURE + WATSON PARTIAL findings reinforce that the missing
> infrastructure (profinite transfer-operator analyticity / spectral gap) is the load-bearing
> rigor gap. No SECONDARY ROUTE in the existing 11-probe corpus closes this gap; each ruled-out
> route demonstrates a different categorical barrier.

---

## SECONDARY ROUTING

Per pre-registration (priority-ordered):

1. **Multi-singularity Flajolet-Sedgewick VI.4-VI.5** (PADE's flagged next step). Specifically:
   - §VI.4 handles complex-conjugate-pair singularities and cos(nθ+φ) coefficient asymptotics
   - §VI.5 handles two real singularities at distinct |z|
   - These are the Darboux machinery in clean form; their RIGOROUS application to f(z) still
     requires KNOWING f(z) is analytic with the claimed singularity structure (the gap).
   - **PRIORITY: MODERATE** (refines form, doesn't bridge the rigor gap)

2. **Direct construction of profinite-analog semiclassical theory.** Research monograph (per
   FAURE_DISPOSITION priority-3). **PRIORITY: LONG-TERM**.

3. **NEW (emergent from WATSON Phase 2): Compute ε_k=14..20+ and test whether ρ_k stabilizes
   at √3 (Faure) or 1.016 (STATE.md slow-mode).** This is a concrete numerical task. If ρ_k
   tightens around one of the two predictions, the structural picture sharpens; if neither, a
   third structural feature is present. **PRIORITY: HIGH** for empirical clarification.

### Recommended top-priority secondary route: **Numerical extension of ε_k to k=20+**

Rationale:
- Adjudicates Faure (√3) vs PADE slow-mode (1.016) directly from data
- Concrete, executable (per existing K_p computation infrastructure used for k=1..13)
- Doesn't depend on building new analytic infrastructure
- Resolves the agnostic CONSISTENCY between Faure and PADE found in WATSON Phase 2
- If ρ_k stabilizes at √3: empirical certification of Faure 2009 framework (without proof)
- If ρ_k continues to 1.016: PADE slow-mode confirmed; framework search continues for the
  proper transfer-operator theory at that radius

The Flajolet-Sedgewick VI.4-VI.5 secondary route is also worthwhile but is somewhat parallel —
it gives the formal Tauberian reading of the same multi-singularity picture without independent
data.

---

## What category of theorem would CLEANLY close the saddle-point / Darboux route

For saddle-point or Darboux to RIGOROUSLY apply to f(z) for the closure problem, we need:

1. **Analyticity of f(z) in {|z| < ρ_∞}**: from a transfer-operator theorem stating T_M has
   spectral radius < 1/ρ_∞ on the relevant invariant subspace (Faure-style result, in profinite
   setting).
2. **Identification of the rightmost singularity** ρ_∞ with algebraic-type properties (pole order,
   branch cut character): from explicit spectral decomposition of T_M.
3. **Bounded extension across the singularity** for Darboux/contour-deformation to work: from
   resolvent / Riesz-projector estimates.

NO existing theorem in the 11-probe corpus literally supplies (1) + (2) + (3) for Syracuse. The
closest in spirit is Faure 2009 (FAURE_DISPOSITION); the next-closest is Pollicott-Sharp dynamical
zeta functions (FAURE candidate J, NO_FIT). Both require smooth-manifold infrastructure.

---

## Surprises in Phase 2 execution

### Surprise 1: r=2..6 direct T_p computation gives κ ≈ 1.17 (super-trivial)

R79b reports κ = 0.522 at r=8..20. Our direct computation at r=2..6 finds κ ≈ 1.17 (log-linear
fit), meaning |T_p| GROWS slightly faster than N at small r. The transient-to-asymptotic
transition for T_p(r) is at r ≈ 6..8, matching the PADE transient-to-asymptotic transition at
k ≈ 10. **Two independent transitions at similar location is a structural signature.**

### Surprise 2: Multi-saddle fits are underdetermined

The ε_k k=2..13 data has only 1 sign change; multi-mode fits with 4-6 parameters can fit the
single change at multiple competing parameter combinations. PADE's specific (ρ_1=1.016,
ρ_2=1.57, θ_2=2π/9.2) is one of MANY consistent models. The data alone doesn't strongly select.

### Surprise 3: Free 5-parameter fit prefers ρ_1 = 0.54 (a growing mode!)

If we let the 5-parameter multi-saddle model freely fit, it picks ρ_1 = 0.54 (modulus < 1)
meaning the "asymptotic" mode is GROWING. This is a fitting artifact: the data at k=2..13 is
not yet in the asymptotic regime, so the optimizer fits transient growth as the "asymptotic"
mode. Reinforces that **the asymptotic ρ_∞ cannot be inferred from k ≤ 13 alone**.

### Surprise 4: Phase deviation of ψ_lead at r=4..5 is LARGE

R78.6 says saddle-exact at r=3. At r=4: phase deviation up to 88°. At r=5: up to 160°. These are
not small Hensel corrections — they're substantial structural deviations. R79b reports
|Σ 1̂·ψ_lead| / |Σ 1̂·ψ_true| ≈ 0.4-0.6 (factor-2 gap); our phase-level findings are even more
dramatic at the individual-a level.

### Surprise 5: Held-√3 and held-1.57 fits have similar RSS

Multi-saddle fits with ρ_2 fixed at √3 (Faure prediction) vs 1.57 (PADE Hadamard at n=13) give
RSS = 1.75e-5 vs 1.71e-5 — within 2% of each other. **The data at k=2..13 cannot
distinguish the two predictions.** This is the agnostic-consistency finding; it's also the
strongest argument for extending ε_k computation.

---

## Strategic position

Pre-WATSON: 10 categorical barriers mapped (5-probe Fourier, BGT, ADELIC, IGUSA, FAURE).

Post-WATSON: **11 categorical barriers mapped.** Saddle-point / Watson lemma / Darboux are
universally applicable tools that supply the EXPECTED FUNCTIONAL FORM for ε_k decay and for
T_p(r) bilinear cancellation. They do NOT deliver rigorous asymptotic rates without the
transfer-operator analytic infrastructure that profinite Syracuse lacks.

The 11-barrier picture is now even sharper:

> Syracuse's c=7/45 decay structure is empirically certified to high precision through k=13.
> Eleven categorical analysis routes — Fourier-decay (5-probe), regular-variation, adelic,
> Igusa local zeta, Faure semiclassical, saddle-point/Watson — have each been mapped and
> each either NO_FIT or PARTIAL with a sharply-identified rigor gap. The intersection of
> what's needed across PARTIALs is: a profinite analytic transfer-operator theory with
> spectral / Riesz-projector estimates for the Tao recursion's T_M. The Watson/saddle-point
> route confirms the expected functional form (complex-conjugate-pair singularity → cos-
> modulated exponential decay) but cannot close the rigor gap on its own.

The PROBE adds two concrete deliverables beyond pre-existing FAURE_DISPOSITION:

1. **Confirmation of triple-PARTIAL convergence**: BGT (plateau), FAURE (form), WATSON (rate).
   All three identify the same missing infrastructure.
2. **Recommended high-priority numerical follow-up**: extend ε_k to k=20+ to adjudicate the
   Faure-vs-PADE ρ_∞ predictions agnostic consistency.

---

## Deliverables

In C:/Collatz/:

- `WATSON_PRE_REGISTRATION.md` — pre-reg (locked 2026-05-14 before Phase 2)
- `WATSON_R78_SUM_FORM.md` — load-bearing-sum extraction + technique routing
- `WATSON_J_DARBOUX_HYPOTHESES.md` — Darboux's method (Temme §2.4) verbatim
- `WATSON_C_SADDLEPOINT_HYPOTHESES.md` — saddle-point / steepest descent (Manton) verbatim
- `WATSON_ASYMPTOTIC.md` — Phase 2 numerical results
- `WATSON_COMPARISON.md` — Phase 3 comparison to PADE + Faure + ε_k k=8..13
- `WATSON_DISPOSITION.md` (this file) — headline

- `watson_phase2_asymptotic.py` — multi-saddle fits to ε_k
- `watson_phase2b_period.py` — period scan
- `watson_phase2c_padeguided.py` — PADE-constrained fits + Hadamard + ratio test
- `watson_phase2d_saddle_extension.py` — direct T_p(r) computation r=2..6 + ψ_lead vs ψ_true

PDF extractions (UTF-8): C:/tmp/watson/*.txt (5 PDFs, all via pypdf).

No git operations performed (per discipline).

---

## Honest scope

This probe was a SINGLE-SESSION executive probe (per brief). The asymptotic computation is
numerically tight (multi-start fits, Hadamard analysis, direct T_p computation). The disposition
PARTIAL is robust: form-match + rate-gap.

What this probe did NOT attempt:
- Numerical extension of ε_k beyond k=13 (deferred to top secondary-route)
- Hensel-lifted explicit closed form of ψ_true at r ≥ 4 (open problem flagged in R78.6)
- Construction of profinite transfer-operator theory (research monograph, FAURE-aligned)
- Flajolet-Sedgewick VI.4-VI.5 formal Tauberian extraction (secondary-route candidate)

The probe's contribution is the SHARP DIAGNOSIS that saddle-point/Darboux supply form but not
rate, and the TRIPLE-PARTIAL CONVERGENCE that BGT + FAURE + WATSON identify the SAME missing
infrastructure.

---

End disposition.
