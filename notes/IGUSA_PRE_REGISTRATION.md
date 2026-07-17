# IGUSA_PRE_REGISTRATION

**Date locked:** 2026-05-14 (before any Phase 1 hypothesis check).
**Probe:** IGUSA — single-theorem selection on the R78 (1+3)^u algebraic substrate; target the pole at s = log_3(2) ≈ 0.631 corresponding to R77.6's branch-cut at z=2 via z = 3^s.
**Mode:** E (verbatim hypotheses extracted from PDFs at C:/Users/Nate/OneDrive/Documents/igusa_local_zeta/pdfs/, pypdf to C:/tmp/igusa/).
**Working dir:** C:/Collatz/.

---

## Substrate inputs (verbatim, no inheritance from elsewhere)

(1) R78 (1+p)^u polynomial substrate (`result_78.md`, `C1_R78_DISAMBIGUATION.md`, `GENERALIZATION_R78_SPECIFIC_FEATURES.md`):

  - The Kalafatelis partial-sum object is `S = Σ_{u=0}^{p^{r-1}−1} e_{p^{r+1}}( c · (1+p)^u − p² m u )` with p=3.
  - The natural polynomial-in-u identification (Cochrane T2 substrate) is:
    g(u) = c · Σ_{k=0}^r C(u, k) · p^k − p² m u   mod p^{r+1}.
  - This polynomial has **D = deg_p H+ = 0** in Cochrane's notation (R78 finding): after factoring τ=1, H(u) mod p has only its constant term non-zero. Higher-derivative polynomial behavior is trivial mod p.
  - The "polynomial substrate" available to Igusa is therefore **g(u) ∈ Z_p[u]** (univariate of degree r), not "(1+3)^u" itself (which is exponential, not polynomial in u over Z).
  - Alternative polynomialization (Postnikov / R78-feature-5): write u = u_0 + p·s and view phase function in s. Phase becomes P_a(s) = p s − C_a · L_p(1+ps), the cubic-phase substrate for R78's chain.

(2) F̂_p magnitude theorem (`FHAT_THEOREM_VERIFICATION_RESULTS.md`):

  THEOREM (verified at 33 cells, p ∈ {3..31}, r ∈ {1..6}): for f_p(u) := e_M(c·(1+p)^u), M = p^{r+1}, u ∈ Z/p^r:
  - supp(F̂_p^full) = {p·a (mod M) : a ∈ Z/p^r, a ≡ c mod p}, cardinality p^{r-1}.
  - |F̂_p^full(ξ)| = p^{(r+3)/2} on support, uniformly.
  - Vanishing off support.

  This is the *Fourier-magnitude saturation*; the Igusa zeta of g(u) (or of the dual sub-support polynomial) is the candidate object whose pole structure should reproduce this magnitude formula via Mellin-of-magnitude.

(3) R77.6 branch-cut numerical anchor (`PADE_NUMERICAL_DISPOSITION.md`, `T_LEAD_CORRECTED_DISPOSITION.md`):

  - R77.6 read a branch-cut at z = 2 in f̃(z) = generating function of ε_k.
  - PADE_NUMERICAL_DISPOSITION (n ≤ 13) **REFUTES** z=2 as the leading singularity: Hadamard radius at n=13 is 1.57, trending inward toward 1.046 (= 1/(43/45), per T_lead) or 1.016 (= 1/0.984, slow-mode fit).
  - Original numerical anchor for this probe: z = 2 ↔ s = log_3(2) ≈ 0.631.
  - **Updated anchors (sub-leading / candidate locations):**
    - s = log_3(2) ≈ 0.6309 (z=2; R77.6 reading, sub-leading per PADE_NUMERICAL_DISPOSITION)
    - s = log_3(45/43) ≈ 0.0413 (z=1.046; T_lead corrected 43/45 reading)
    - s = log_3(1/0.984) ≈ 0.0147 (z=1.016; STATE.md slow-mode reading)
  - The probe targets ALL THREE as candidate pole locations; the load-bearing question is whether *any* Igusa-style theorem of (1+3)^u substrate produces ANY of these.

(4) Multi-regime ε_k pattern (`result_77_7_eps_exact_through_k8_v2_vec_pool.json`):

  |ε_k|·2^k values k=1..8 (decimal): 0.4, 0.0381, ≈0.0102, ≈0.0392, ≈0.0337, ≈0.0316, ≈0.302 (jump), ≈0.220.

  Plateau k=2..6 at ~0.03-0.04, jump at k=7 to ~0.3, post-jump k=8 ≈0.22. Multi-regime structure load-bearing for Phase 2 conclusion-shape check.

(5) Tao recursion form (`C1_TAO_RECURSION_FORM.md`):

  μ̂_n(ξ) = E[χ_n(2^{-a_1} + 3·2^{-a_{[1,2]}} + … + 3^{n-1} 2^{-a_{[1,n]}})] — additive character on a 2D renewal walk over (a_1,…,a_n) ∈ (ℕ+1)^n.
  Project-internal one-step recursion: μ̂_{n+1}(ξ) = Σ_v 2^{-v} A_v(ξ) μ̂_n(ξ·2^{-v} mod 3^n) with A_v a phase coefficient.

---

## Candidate list (verbatim PDF location pre-Phase-0)

| Code | Theorem | Source |
|---|---|---|
| A | Igusa rationality theorem | A_Denef_Bourbaki.txt §1.3.2 (Theorem (Igusa)) |
| B | Denef-Hoornaert Newton-polyhedron explicit formula | E_Bories_Veys_Surfaces.txt §0.4 Thm 0.27; J_Newton_Polyhedra.txt Thm 2.9 |
| C | Bories-Veys non-degenerated surface singularities | E_Bories_Veys_Surfaces.txt Thm 0.12 |
| D | Monodromy conjecture (Igusa-Denef-Loeser) | C_Veys_Monodromy.txt Conj 2.12; J_Newton_Polyhedra.txt Conj 3.1 |
| E | Motivic Igusa (Denef-Loeser 1998) | B_Denef_Loeser_Motivic.txt; J_Newton_Polyhedra.txt Thm 5.5 |
| F | Cluckers / Bories-Cluckers polynomial mapping | I_Bories_Cluckers_Polymap.txt Thm 2.5 |
| G | Stationary phase / Igusa local singular series (Watson on Igusa) | A_Denef_Bourbaki.txt §1.4 + 4.1 (asymptotic of F(y) as y→0) |
| H | Bernstein-Sato b-function ↔ Igusa pole relation | C_Veys_Monodromy.txt Thm 2.6, Cor 2.7 |
| I | Igusa for monomials / monomial-like polynomials (Denef Bourbaki §1.3, Veys Thm 3.4) | A_Denef_Bourbaki.txt §1.3, C_Veys_Monodromy.txt Thm 3.4 |
| J | Igusa for c·(1+p)^u-type substrates **specifically** | Open search across corpus |
| K | Any candidate surfaced during reading | (TBD) |

---

## Pre-registered priors per candidate

| Code | Prior | Reasoning |
|---|---|---|
| A (Rationality) | LOW for closure (foundational) | Gives Z(s,f) rational; candidate-pole list non-explicit without resolution data. Not enough to predict pole at log_3(2). |
| B (Denef-Hoornaert) | **HIGHEST** | Explicit formula via Newton polyhedron Γ(f); pole list reads off facets of Γ(g). If g(u) = c · Σ C(u,k) p^k − p² m u has facets producing s = −σ(v)/m(v) ∈ {log_3(2), log_3(45/43), …}, then SELECTED. |
| C (Bories-Veys surfaces) | LOW | Restricted to n=2 (surface singularities, polynomial in 3 vars per Thm 0.12). Syracuse g(u) is **univariate** — wrong dimension. |
| D (Monodromy) | MODERATE | Predicts exp(2πi Re(s_0)) is a monodromy eigenvalue. Phase 2 check: if Igusa pole comes out, monodromy of g(u) should have eigenvalue exp(2πi · log_3(2)). |
| E (Motivic) | MODERATE | Motivic level cleaner; specializes to p-adic. Same input poly. |
| F (Bories-Cluckers polymap) | LOW | Polynomial mapping (multi-component f); R78 substrate is single polynomial. |
| G (Stationary phase / Watson on Igusa) | HIGH for asymptotic | If a pole exists, the asymptotic of F(y) = ∫ e_M(c·(1+p)^u) ··· involves the pole's residue. Gives the |F̂_p^full| magnitude reconstruction. |
| H (Bernstein-Sato) | MODERATE | b-function of g(u) — if known, gives pole list. For univariate g of degree r, b_g(s) is computable. |
| I (Monomial / linear) | **HIGHEST if reduces** | g(u) mod p is **constant** (R78's D=0 finding) — degenerates. Need to check whether g(u) reduces to a monomial form after lift. |
| J (Specific c·(1+p)^u literature) | UNKNOWN — open search | If exists, conclusive. |

---

## Decision rules (locked)

For each candidate K:

- **SELECTED:** Phase 0 verbatim ✓ + Phase 1 every hypothesis SATISFIED for at least one polynomialization of R78 substrate + Phase 2 produces Z(s, f, 3) with pole at s ∈ {log_3(2), log_3(45/43), log_3(1/0.984)} + Phase 3 conversion to closure via R75/R76/R77 explicit (or constructive sketch).
- **PARTIAL:** Candidate fires (Phase 1 ✓) but Phase 2 produces pole list NOT containing any of the targets, OR contains a target but conversion to closure fails / requires additional structural input.
- **NO_FIT:** ≥1 hypothesis FAILED at Phase 1 across all polynomializations.
- **BLOCKER:** Theorem statement UNVERIFIABLE in corpus.
- **MODE_H_CIRCULAR:** Theorem requires |μ̂_n(ξ)| analytic continuation or similar closure-target object as INPUT.

**Locked thresholds:**
- "Pole at log_3(2)" means |s_actual − log_3(2)| < 0.01 (i.e., within 1% of the target).
- "Pole at log_3(45/43)" means |s_actual − 0.0413| < 0.005.
- "Pole at log_3(1/0.984)" means |s_actual − 0.0147| < 0.005.
- Real-part match suffices (poles can have imaginary shifts 2πik/(N_i log p) by Igusa Thm).

**Mode H risk specifically LOW for Igusa:** Igusa hypotheses are polynomial-data (degree, Newton polyhedron, non-degeneracy). The closure target s = log_3(2) is derivable FROM polynomial data, not the other way. Mode H circular fingerprint should be absent unless the polynomial f itself requires closure-target data to define.

---

## Pre-registered SECONDARY ROUTING if NO_FIT/PARTIAL

1. **Faure 2009 semiclassical spectral gap** — partially expanding map transfer operator on Tao recursion (categorically distinct from Igusa; addresses chain-side operator).
2. **Watson lemma / saddle-point on R78/R79 bilinear off-diagonal sum** — direct asymptotic in B=1 vdC framework; bypasses Igusa polynomial-extraction concern.
3. **Adelic Mellin functional-equation extension with new substrate** — return to ADELIC route with explicit (1+3)^u algebra fed in; cross-frequency M-saturation.
4. **Direct computation of Z(s, g, 3) for the specific Syracuse polynomial via Mainfile thesis rationality proof** — last resort. Compute the embedded resolution of g(u) ∈ Z_3[u] explicitly and read off candidate poles.

---

## Honest probability priors (pre-Phase-1)

- SELECTED: 20-25%
- PARTIAL: 30%
- NO_FIT: 30%
- BLOCKER: 10%
- MODE_H_CIRCULAR: 5-10%

Adjusted (down) from the brief's 20-30% SELECTED estimate after noting:
- **R78's polynomial g(u) is UNIVARIATE.** All Newton-polyhedron / explicit-formula machinery (B, C, F) is designed for multi-variable polynomials. Univariate degenerate case is handled by Igusa rationality (A) but the candidate-pole list reduces to roots of g — and g's roots in Z_3 control the pole list.
- **R78's D=0 finding** says g(u) mod 3 is a CONSTANT (≡ c, a unit). This means g has NO ROOTS in Z_3 (it's a unit everywhere — `|g(u)|_3 = 1` uniformly!). The Igusa zeta `Z(s, g, 3) = ∫_{Z_3} |g(u)|^s du = ∫_{Z_3} 1 du = 1` — TRIVIAL. No poles at all.

This is a structural pre-warning: the natural univariate polynomialization of the Syracuse substrate gives a TRIVIAL Igusa zeta (no poles). The non-trivial Igusa zeta must operate on a different polynomial — likely the Cochrane-Pinner cubic phase P_a(s) (R78 Feature 4 / 5), which is a polynomial in s with parameters a. This shifts the probe to a *parametrized* Igusa zeta (Z(s, f_a, 3) as a function of a), or to a *higher-dimensional* substrate where both s and a are integration variables.

---

End pre-registration. Locked 2026-05-14 before any hypothesis check.
