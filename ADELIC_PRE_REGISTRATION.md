# ADELIC_PRE_REGISTRATION — adelic Mellin / Tate-style probe

**Date:** 2026-05-13. Probe ADELIC. Working dir: C:/Collatz/. Locked before Phase 0 PDF extraction.

---

## Probe scope

Single-theorem selection for the c=7/45 closure question via the adelic Mellin / Tate-style framework. Tests whether any theorem establishes the singularity structure of F_∞(s) for Syracuse μ_n's adelic Mellin from chain-side data without circularly assuming closure.

Load-bearing prior: BT_DISPOSITION's archimedean-place finding (the 1-attractor visibility requires adelic / global / multi-place substrate; 3-adic and 2-adic local data alone cannot see it).

---

## Inputs (verbatim, no inheritance)

1. **Syracuse μ̂_n form + R75 Plancherel decomposition.** Per C1_TAO_RECURSION_FORM.md and c_seven_forty_fifth.md.
   - μ̂_n(ξ) = E[χ(2^{-a_1} + 3 · 2^{-a_{[1,2]}} + … + 3^{n-1} · 2^{-a_{[1,n]}})] with a_i iid Geom(2) on ℕ+1, χ(x) = exp(-2πi ξ (x mod 3^n) / 3^n).
   - State space: (Z/3^n)*.
   - R75: c = (1/3) · lim_{n→∞} S_∞ with S_n = Σ_{ξ ∈ Z/3^n, 3∤ξ} |μ̂_n(ξ)|² → 7/15.
   - Plancherel split S_{k+1} = S_k + Off-diag(k) (proved).

2. **BT archimedean-place finding (verbatim).**
   - From BT_DISPOSITION.md §"Negative case Q3": "The 1-attractor is an *archimedean* phenomenon — r_n stops being > 1 archimedean-ly. The 3-adic place can't see archimedean convergence; the 2-adic place is consumed driving the recursion forward. So the attractor lives at the *archimedean* place, which the p-adic Bruhat-Tits machinery is built to ignore. Any substrate that can see the attractor must include the archimedean place — i.e., must work over **adelic** or **global** geometry, not p-adic-only."
   - Salvageable picture: M_v = [[3,1],[0,2^v]] is a v-indexed family of hyperbolic isometries of T_3 with eigenvalues {3, 2^v}, all sharing fixed point ∞ ∈ ∂T_3, forming a pencil of axes in the Borel subgroup B ⊂ PGL_2(ℚ_3) = stabilizer of ∞.
   - Structural reason for adelic visibility: the attractor is a global integer-1, not a 3-adic limit — requires real place ∞ + 2-adic place + 3-adic place simultaneously, since recursion is "multiply by 3 (3-adic-driven), divide by power of 2 (2-adic-driven), measure size at ∞ (archimedean)".

3. **R77 T_diag two-mode + R77.6 branch-cut.** Per result_77_T_lead_spectrum.md and result_77_6_generating_function.md.
   - T_diag = (1/5)·[[1,1],[4,4]], eigenvalues {0, 1}, eigenvectors (1,−1) and (1,4).
   - Off-diagonal rate ½ empirical at k=2..6.
   - R77.6: f̃(z) := (E(z) − ε_1·z)/z² where E(z) = Σ ε_n z^n has branch-cut singularity at z = 2 (verdict G-branch-cut, type indeterminate at N = 5; BGT_DISPOSITION confirms via k=7 jump in ε_k post-2026-05-12 data).

4. **R78 (1+3)^u algebraic structure.** Per result_78.md and GENERALIZATION_R78_SPECIFIC_FEATURES.md.
   - F̂_p^full magnitude formula has algebraic substrate (1+3)^u; principal-unit-coset 1 + 3·Z_3 mod 3^{r+1} parametrization, with truncated 3-adic log L_3(1+3) as phase.
   - Two load-bearing R78-specific features: cubic phase degree (feature 3) and one-parameter principal-unit setup (feature 5). Both ARE present here.
   - Cochrane Theorem 2 attack on Kalafatelis eq 190: FAILED (D=0 obstruction).
   - van der Corput differencing (R79): FAILED (B=1 gives rate 0.73, B=2 worsens; band-l¹ cancellation required, not pointwise).

5. **ε_k multi-regime + BGT PARTIAL (per eps_exact_through_k8 JSON + BGT_DISPOSITION.md).**
   - ε_k exact rationals through k=8:
     ε_1 = 1/5 = +0.200, ε_2 = 1/105 = +9.52e-3, ε_3 ≈ −5.09e-3, ε_4 ≈ −2.45e-3, ε_5 ≈ −1.15e-3, ε_6 ≈ −4.98e-4, ε_7 ≈ −1.176e-3 (wait — re-check; actually from JSON these are the deviations S_k − 7/15, signs match table 4 in c_seven_forty_fifth.md).
   - |ε_k|·2^k pattern: plateau k=2..6 at ~0.04; k=7 jump to 0.150; k=8 = 0.191 (per memory line).
   - Sign pattern (+,+,−,−,−,−,−,−) confirmed in JSON: ε_1>0, ε_2>0, ε_3..ε_8 all negative.
   - BGT PARTIAL: sequential RV fires in plateau k=2..6; fails across jump.
   - k=7 jump independently confirmed by Chevalier 1.16 (Tauberian) and BGT slow-variation as STRUCTURAL, not finite-N artifact.

---

## Candidate theorems

Phase 0 must locate each verbatim in the pulled corpus (C:/Users/Nate/OneDrive/Documents/adelic_mellin/pdfs/):

- **A.** Tate's adelic Mellin functional equation (Ramakrishnan-Valenza Ch. 7 review; Binder; Poonen). L(s, χ) for idele class characters; functional equation via adelic Poisson summation.
- **B.** Local archimedean Mellin factor (Tate gamma factor structure). Singularity structure of F_∞(s) for a measure on ℚ_∞^* = ℝ^*.
- **C.** Local non-archimedean Mellin factor (Tate epsilon / local L-factor). Singularity structure of F_p(s) for a measure on ℚ_p^*.
- **D.** Adelic Poisson summation (Tate Lemma 2.2.1 in Binder; Lemma 4 in Ramakrishnan-Valenza review).
- **E.** Cartwright-Kaimanovich-Woess 1994 (random walks on affine group Aff(F) of local field F; limit theorems, harmonic functions).
- **F.** Anker-Schapira-Trojan 2007 & 2013 (heat-kernel sharp asymptotics on affine buildings Ã_r).
- **G.** Chambert-Loir-Tschinkel 2009 (Igusa integrals + adelic volume asymptotics; Mellin transforms of height functions).
- **H.** Bourgain-Gamburd-Sarnak (generalized Selberg 3/16 + affine sieve via SL_2(Z) infinite-index subgroups).
- **I.** Kontorovich 2014 (levels of distribution for affine sieve; polynomial-in-A bounds).
- **J.** Saloff-Coste Notices survey (random walks + invariant diffusions on locally compact groups).
- **K.** Any candidate surfaced during Phase 0 reading.

---

## Pre-registered priors

| Code | Theorem | Prior P(SELECTED) | Rationale |
|---|---|---|---|
| A | Tate adelic Mellin functional equation | 18% | Mode H escape candidate; closure-direct IF Syracuse admits Tate-style FE |
| B | Local archimedean Mellin factor | 8% | Depends on what natural archimedean object Syracuse pushes to |
| C | Local non-archimedean Mellin factor | 25% | High — the 3-adic factor F_3(s) should be computable from R77 T_diag + R78 (1+3)^u |
| D | Adelic Poisson summation | 12% | Proof tool inside Tate FE, not standalone closure |
| E | Cartwright-Kaimanovich-Woess | 22% | Highest among non-Tate candidates; directly addresses random walks on Aff(local field) and (1+3)^u lives in Aff(ℚ_3) Borel B |
| F | Anker-Schapira-Trojan heat kernel | 15% | Heat-kernel framework match noted earlier; asymmetry problem is real (BT and FG flagged) |
| G | Chambert-Loir-Tschinkel Igusa | 10% | Adelic factorization inherent but Mellin-of-height is different object than Markov-stationary Mellin |
| H | BGS affine sieve | 8% | Different object (infinite-index SL_2(Z) subgroup vs profinite (Z/3^n)*) |
| I | Kontorovich levels | 7% | Tracks BGS prior |
| J | Saloff-Coste | 3% | Already failed in FG probe via reversibility |
| K | Surfaced candidate | reserved | TBD |

Sum: ~128% — priors are non-disjoint; multiple candidates could simultaneously fire (e.g., A+D both via functional equation; C and E both via 3-adic factor).

P(at least one SELECTED) ≈ 12% (honest, given seven prior NO_FIT/PARTIAL outcomes).

---

## Pre-registered honest priors (final-disposition outcomes, exhaustive, sum to 100%)

- SELECTED: 12% — at least one candidate fires cleanly, deriving F_∞ singularity from F_p data.
- PARTIAL: 30% — framework fires for one place factor (likely F_3 via R77 T_diag) but archimedean F_∞ requires separate analysis.
- NO_FIT: 30% — adelic framework requires structural inputs Syracuse doesn't have (eighth category-of-object barrier).
- MODE_H_CIRCULAR: 15% — F_∞ analytic continuation hypothesis IS the closure target.
- BLOCKER: 13% — lit pull missing key Tate-extension paper.

---

## Decision rules

For each candidate K:

- **SELECTED:** hypotheses SATISFIED by inputs (1)+(2)+(3)+(4)+(5) verbatim AND theorem produces explicit polynomial-in-A bound on |μ̂_n(ξ)| OR closes Mode H gap (derives F_∞ singularity from chain-side data).
- **PARTIAL:** hypotheses SATISFIED but conclusion qualitative not quantitative, or quantitative with non-explicit constants. Report gap.
- **NO_FIT:** at least one hypothesis FAILED.
- **BLOCKER:** theorem statement UNVERIFIABLE in available corpus.
- **MODE_H_CIRCULAR:** theorem requires F_∞ analytic continuation = closure target.

---

## Adelic factorization tagging (Phase 3)

Each candidate gets ONE of:

- **ARCHIMEDEAN_VISIBLE:** delivers results at the archimedean place specifically.
- **NON_ARCH_ONLY:** only at non-archimedean factor (p-adic only).
- **GLOBAL_BUT_PLACE_BLIND:** global statement but doesn't distinguish places.
- **ADELIC_FACTORIZATION_INHERENT:** factorization across all places is built in.

Pre-classification (subject to verbatim verification):
- A: ADELIC_FACTORIZATION_INHERENT
- B: ARCHIMEDEAN_VISIBLE (by definition; question is whether hypotheses fit Syracuse)
- C: NON_ARCH_ONLY (by definition)
- D: ADELIC_FACTORIZATION_INHERENT
- E: NON_ARCH_ONLY (local field, doesn't see archimedean place of ℚ)
- F: NON_ARCH_ONLY (affine building over p-adic field) unless Riemannian-symmetric-space analog invoked
- G: ADELIC_FACTORIZATION_INHERENT
- H: GLOBAL_BUT_PLACE_BLIND (uses archimedean ordering but doesn't produce a factor F_∞)
- I: as H
- J: depends on group; mostly GLOBAL_BUT_PLACE_BLIND

---

## Locked thresholds

- "Hypothesis SATISFIED" requires verbatim match of hypothesis to one of inputs (1)-(5), with no ambiguity.
- "Conclusion delivered" requires explicit statement about F_∞(s) (singularity location, type, prefactor) OR a polynomial-in-A bound |μ̂_n(ξ)| ≤ n^{-A} with explicit constants.
- "Mode H circular" if a hypothesis includes "F_∞(s) has meromorphic continuation past Re(s) = σ_0" with σ_0 = the relevant critical line.

---

## Pre-registered SECONDARY ROUTING (if NO_FIT/PARTIAL)

In order of preference:

1. **Igusa local zeta** — operates on R78 (1+3)^u directly. Categorically distinct from prior probes. Bypasses multi-regime obstruction. (Recommended top priority by BGT_DISPOSITION.)
2. **Faure 2009 semiclassical spectral gap** — partially-expanding maps spectral gap, directly addresses R77 off-diagonal rate-½. (Just added to FG corpus per BGT_DISPOSITION secondary routing.)
3. **Heat-kernel narrowing** if heat-kernel candidate F is PARTIAL — specific to Anker-Schapira-Trojan on Ã_r tree.
4. **Watson lemma / saddle-point on R78/R79 bilinear** — operates closer to chain-side; complementary to BGT.

---

## Mode E discipline

- Theorem hypotheses VERBATIM from PDF (no inheritance from project files).
- Mode H awareness ELEVATED for ADELIC: candidates that require "analytic continuation of F_∞ past a vertical line" or "meromorphic continuation of L-function" as INPUT are Mode H circular. Candidates that DERIVE F_∞ from F_p via functional equation are Mode H escapes — these are the load-bearing candidates.
- Functional equation specifically: check whether Syracuse admits a TATE-STYLE functional equation. If yes (rare), this is closure route. If no, find which natural completion of Syracuse to the adeles does admit one.
- pypdf for extraction (NOT pdftotext). UTF-8 file write before read.
- Don't propose options at the end. Report data + structure.
- Don't git commit.

---

## Anticipated structural questions

(1) **Does Syracuse μ_n admit a natural lift to Q_3^*, then to A_Q^*?** The state space (Z/3^n)* projects to Z_3^* via inverse limit. Z_3^* lifts to Q_3^* trivially (just include 0 — but μ_n is supported away from 0, so the lift is well-defined). The adelic lift requires placing the measure at the 3-adic component and a TRIVIAL/DELTA structure at all other places. This is NOT a natural adelic measure in the Tate sense (which requires non-trivial archimedean component). The question is whether Tate's FE applies to such a "purely 3-adic" measure.

(2) **R77.6 branch-cut at z = 2: which adelic local factor singularity does this correspond to?** If F_3(s) is the 3-adic Mellin factor of μ_n's lift, its singularities live on Re(s) = 0 (the natural Tate pole structure for 3-adic Mellin) or on Re(s) = 1 (functional-equation reflection). The translation z → s requires identifying z with q^s for some prime power q; for q = 3, z = 2 corresponds to s = log_3(2) ≈ 0.631, which is NOT a standard Tate pole location for ℚ_3. For q = ? ≠ 3, this might correspond to a different prime. Pre-register: if the branch-cut at z=2 maps cleanly to a Tate singularity at some s_*, this is structural evidence; if not, the generating function and adelic Mellin are different objects.

(3) **R77 T_diag spectrum {0, 1} on (1,-1) and (1, 4):** Compatible with classical Cartier-Hashizume p-adic spherical function spectrum? The Cartier 1973 formula gives γ(z) = (3^z + 3^{1-z})/4 with branch points at z = 0 and z = 1 (matching the {0, 1} spectrum). This is suggestive but the Tao recursion is not a K-bi-invariant walk on T_3 (per BT_PHASE0_THEOREMS T7), so the spherical-function machinery doesn't directly apply.

---

End pre-registration. Phase 0 proceeds.
