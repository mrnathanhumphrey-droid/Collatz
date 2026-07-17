# TAUBERIAN_RESCOPE_PRE_REGISTRATION

**Date:** 2026-05-13.
**Probe:** Tauberian arc re-scope — single-theorem selection for c=7/45 closure.
**Mode:** E (verbatim theorem hypotheses from PDF, no inheritance from prior project files).

This file is committed BEFORE any selection or matrix evaluation.

---

## Context — prior probe arcs closed NO-GO

1. 5-probe modern Fourier-decay arc (`POLYNOMIAL_IN_A_LANDSCAPE.md`) — continuous-smooth-dynamical mismatch.
2. Cluster C1 (Cochrane / BC / HB / BGK discrete exp sums) — single-complete-sum mismatch (`C1_DISPOSITION.md`).
3. Cluster C2 (BMP / PSF / cut-and-project) — encoding partial; weight layer Mode-H circular (`C2_DISPOSITION.md`).
4. BT (Bruhat-Tits + BKL billiards) — p-adic-only substrate cannot see archimedean 1-attractor (`BT_DISPOSITION.md`).

Tauberian arc is the SOLE remaining live route. Categorically aligned: operates on generating-series / Dirichlet-series / Fourier-coefficient asymptotics → "different object" than the four closed routes.

---

## Inputs (verbatim, no inheritance)

(1) ε_k EXACT-RATIONAL SEQUENCE, k = 1..8. Pulled directly from `experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json` (exact Fraction):

| k | ε_k | |ε_k| · 2^k |
|---|-----|-----------|
| 1 | 1/5 | 0.400000 |
| 2 | 1/105 | 0.038095 |
| 3 | -5191/1019445 | 0.040736 |
| 4 | -1.135e16 / 4.627e18 | 0.039236 |
| 5 | -9.243e58 / 8.026e61 | 0.036856 |
| 6 | -4.924e188 / 9.890e190 | 0.031866 |
| 7 | huge negative rational | 0.150430 |
| 8 | huge negative rational | 0.190860 |

**Empirical observations on the input (1):**
- k=2..6: normalized magnitudes nearly flat in the 0.03-0.04 band — consistent with |ε_k| ~ C · 2^{-k}, i.e. geometric decay at rate 1/2.
- k=7, 8: sharp jump back up to 0.15, 0.19. **The sequence is NOT monotonically decreasing in normalized norm.**
- All ε_k for k ≥ 3 are negative (sign pattern +, +, -, -, -, -, -, -).
- ε_1 is the only "boundary" term (the only one with a known closed simple form: 1/5).

This sign-then-jump pattern is itself informational and may bear on which Tauberian theorem (if any) applies.

(2) C1 RENEWAL-WALK FORM. Tao's working form (Tao 1909.03562 eq 7.5):
S_χ(n) is a **product over a 2D renewal walk** {(j, b_{[1,j]}) : j ∈ [n/2]} of conditional expectations of an additive character, indexed by a Geom(2)^n tuple-space:
  |S_χ(n)| ≤ E ∏_{j∈[n/2]} | f(3^{2j-2} 2^{-b_{[1,j]}}, b_j) |
The phase function within each f-factor is 2-adic exponential 2^{-l+1} viewed mod 3^n. Domain: tuple-space (ℕ+1)^n with Geom(2) weights, NOT Z/3^n.

(3) C2 BMP F_1 SUPPORT DIFFRACTION. The Syracuse SUPPORT (3-coprime integers) is a regular weak model set in the BMP super-singular scheme (G=ℝ, H=∏'_p ℚ_p, L=ℚ, W=ℤ_3* × ∏_{p≠3} ℤ_p), maximal density 2/3, ∂W=∅. Diffraction: γ̂_support = ω_{|1̂_W|²} supported on rationals with denominator a power of 3 — pure point. The WEIGHTS layer (Syracuse μ_n Markov stationary distribution) requires h ∈ P_K(H), which is essentially the polynomial-in-A target itself (Mode H circular).

(4) BT ARCHIMEDEAN-PLACE FINDING. The 1-attractor of Collatz lives at the **archimedean place**. The 3-adic Bruhat-Tits machinery is built to ignore the archimedean place; the 2-adic place is consumed driving the recursion forward. Any substrate that can see the 1-attractor must include the archimedean place — i.e., must be adelic/global, not p-adic-only. Concretely: ℚ in ℝ × ∏'_p ℚ_p (full BMP super-singular scheme).

---

## Candidate list (read VERBATIM from PDF; no paraphrase)

**Primary (per brief):**
- A. Flajolet-Sedgewick "Analytic Combinatorics" Ch. VI singularity analysis (Theorems VI.1, VI.3, VI.4 and VI.5).
- B. Chevalier 2507.15394 Theorem 1.16 — Tauberian square-root + pole of order M, explicit M parameter. **PRIOR FAVORITE per memory.**
- C. Korevaar 2002 Bull. AMS — Wiener-Ikehara (Theorem 4.2) and the Newman-style Analytic Theorem (Theorem 6.1 + Section 5 finite-form).
- D. Newman 1980 / Zagier 1997 — Analytic Theorem (Section 5 of Zagier; the "if g is holomorphic on Re(z) ≥ 0 then ∫f(t) dt converges and equals g(0)" form).
- E. Alberts 2508.20814 — Explicit Tauberian for Dirichlet series with twisted-moment vertical bound input.
- F. Singha Roy 2511.15928 — Landau-Selberg-Delange for L-functions (Theorem 1.1).
- G. Tao 2020 Notorious Collatz (slide deck — likely a literature pointer rather than an applicable theorem; check for Lagarias references).

**Secondary (fall-back if A-F all fail):**
- H1. Borwein survey "Century of Tauberian."
- H2. Holland Abel limit theorem.
- H3. Ingham Euler-Maclaurin arxiv 1910.03036.
- H4. Selberg-Delange remarks.
- H6. Häggström basics of Tauberian.
- H9. Guide to Tauberian arithmetic applications (arxiv 2504.16233).
- H10-H12. Lagarias 3x+1 overview / generalizations / stochastic models (likely literature pointers).

Reading order: B, A, C, D first (these are the most likely to fit a generating-series object with explicit asymptotic-expansion deliverable). Then E, F. Then G (literature pointer). Then H if needed.

If any reading turns up a theorem in the candidate corpus not in this list, add it.

---

## Disposition categories (locked thresholds)

For each candidate K we will determine one of:

- **SELECTED.** All hypotheses SATISFIED by inputs (1)+(2)+(3)+(4) verbatim. Where the theorem has an explicit numeric parameter (e.g. Chevalier M, FS α, Wiener-Ikehara A) — that parameter is realizable from inputs (1)-(4) without invoking unproven properties of μ_n or ε_k.

- **NO_FIT.** At least one hypothesis FAILED for every reading-direction of every candidate object construction (full μ̂_n, generating series of ε_k, generating series of |ε_k|·2^k, partial-sum sequences, etc.). Concretely: there is no way to put inputs (1)-(4) into the theorem's hypothesis slots without violating a verbatim condition.

- **PARTIAL.** Hypotheses SATISFIED but the theorem's load-bearing parameter (M, α, T, η, β, …) is infeasible from inputs OR is the target itself (Mode H circular). Report what value of the parameter would unblock and what would have to be true of the inputs for that value to be realized.

- **BLOCKER.** Hypotheses UNVERIFIABLE from the inputs as given. Specify which input is missing (e.g. need ε_9..ε_K for K-term asymptotic; need ρ ∈ P_K(ℤ_3) which is the C2 Mode-H trap; need explicit Dirichlet-series functional equation; need vertical-moment bound on twisted L-integral).

## Locked thresholds

- **For B (Chevalier 1.16):** SELECTED requires (i) there exists a generating series g(z) = Σ b_n z^n with b_n related to ε_n by an explicit rational-functional translation (signs and the |·|·2^k normalization are tolerable), (ii) g(z) admits a continuous extension to D (closed unit disc) holomorphic on D̊, (iii) ∃ meromorphic h_p in a neighborhood of D(1,1)^{1/2} with a single pole at 0 of multiplicity M ≥ 1 such that g(z) = h_p(√(1-z)) for all z ∈ D̊, (iv) M is realizable from inputs (i.e. computable from ε_k pattern, not arbitrary). If any of (i)-(iv) fails for every reasonable g construction → NO_FIT or PARTIAL (PARTIAL if h_p exists but M infinite/unbounded).

- **For A (FS VI.1-VI.5):** SELECTED requires (i) a generating series f(z) = Σ a_n z^n in the input that is Δ-analytic in some Δ(φ, R) (analytic continuation past unit circle except for boundary singularities at finite set of points on |z|=1), (ii) singular expansion of the form f(z) = σ(z) + O(τ(z)) as z → ζ with σ in the standard scale S = {(1-z)^{-α} λ(z)^β : α, β ∈ C}. SELECTED gives [z^n]f(z) ~ ζ^{-n} · n^{α-1} (log n)^β / Γ(α). If standard-scale expansion is not provable from inputs → BLOCKER or NO_FIT.

- **For C (Wiener-Ikehara, Korevaar 4.2):** SELECTED requires (i) an arithmetic monotone-nondecreasing-after-shift "counting function" S(t) constructible from ε_k or |ε_k|·2^k or |μ̂_n(ξ)|² partial sums, (ii) Laplace-Stieltjes transform f(z) = ∫ e^{-zt} dS(t) analytic on Re(z) > 1, (iii) f(z) - A/(z-1) has continuous extension to Re(z) = 1. Output: e^{-t}S(t) → A. If S(t) cannot be made monotone nondecreasing → NO_FIT.

- **For D (Newman-Zagier Analytic Theorem):** SELECTED requires (i) f: [0, ∞) → ℝ bounded and locally integrable, (ii) g(z) := ∫_0^∞ f(t) e^{-zt} dt extends holomorphically to Re(z) ≥ 0. Concludes ∫_0^∞ f(t) dt converges. Need to construct f from ε_k or μ̂_n. If f cannot be made bounded → NO_FIT.

- **For E (Alberts 2508.20814):** SELECTED requires (i) N, N̂ : ℝ≥0 → ℂ with N̂ nondecreasing and |N(X)| ≤ N̂(X), (ii) Mellin transform L(s, F) absolutely convergent on Re(s) > σ_a, (iii) L(s, F) has meromorphic continuation to Re(s) ≥ σ_a - δ with finitely many poles, (iv) twisted-moment bound ∫_T^{2T} |L(σ+it, F) Z^{it}| dt ≪ T^η (log T)^β at the left edge σ = σ_a - δ.

- **For F (Singha Roy LSD 1.1):** SELECTED requires (i) {a_n} ⊂ ℂ has property P(ν, {α_χ}_χ; c_0, Ω) — Dirichlet-character decomposition of the Dirichlet series F(s) := Σ a_n n^{-s} with explicit log-derivative description in terms of Dirichlet L(s, χ) characters mod q, (ii) growth bound Σ_{x<n≤2x} |a_n| ≤ κ x^{1/ν}.

- **For G (Tao 2020):** SELECTED only if Tao identifies a SPECIFIC Tauberian theorem that fits the c=7/45 form. If only a literature pointer to Lagarias bibliography → BLOCKER (insufficient: a pointer is not a theorem statement).

- **For H1-H12 (secondaries):** Only invoked if all of A-F + G land NO_FIT or BLOCKER. Threshold for SELECTED: same as primary — all hypotheses verbatim-satisfied with realizable parameters.

---

## Mode discipline reminders

- **Mode E (verbatim hypothesis extraction).** Every hypothesis sourced from PDF text, not memory or prior project files.
- **Mode H (target-object trap).** If a candidate's hypothesis is *equivalent to the closure target* (e.g. requiring μ_n Fourier-decay-on-Dirichlet-side as a hypothesis), it gets PARTIAL not SELECTED.
- **Lit-scan-narrow-first** (per feedback memory): no scanning beyond the 20-PDF corpus on first pass. If all 20 fail, the disposition reports what category of theorem is missing, not "go find more."
- **Hunt-don't-decide.** This file does not recommend; it reports.

---

## Pre-registered outcomes (probabilistic — committed before evaluation)

| Candidate | SELECTED | PARTIAL | NO_FIT | BLOCKER |
|---|---|---|---|---|
| A (FS VI.1-VI.5) | 20% | 30% | 35% | 15% |
| B (Chevalier 1.16) | 15% | 35% | 30% | 20% |
| C (Wiener-Ikehara) | 10% | 25% | 50% | 15% |
| D (Newman-Zagier) | 10% | 20% | 55% | 15% |
| E (Alberts explicit) | 5% | 20% | 60% | 15% |
| F (LSD Singha Roy) | 5% | 20% | 60% | 15% |
| G (Tao 2020) | 0% | 5% | 30% | 65% (pointer-only) |
| H secondaries | aggregate 5% SELECTED | 30% PARTIAL | 50% NO_FIT | 15% BLOCKER |

Overall "≥1 SELECTED" prior: ~40%. Overall "all PARTIAL or BLOCKER" prior: ~45%. Overall "all NO_FIT" prior: ~15%.

---

## Files to produce (committed before reading PDFs)

- `TAUBERIAN_RESCOPE_PRE_REGISTRATION.md` (this file)
- `TAUBERIAN_RESCOPE_A_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_A_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_B_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_B_HYPOTHESIS_CHECK.md`, `TAUBERIAN_RESCOPE_B_M_PARAMETER.md`
- `TAUBERIAN_RESCOPE_C_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_C_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_D_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_D_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_E_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_E_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_F_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_F_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_G_HYPOTHESES.md` (optional, may collapse to "literature pointer" finding)
- `TAUBERIAN_RESCOPE_H_HYPOTHESES.md` (only if needed; combined for secondaries)
- `TAUBERIAN_RESCOPE_DISPOSITION.md` — top-level summary

No git operations. Nathan commits manually.
