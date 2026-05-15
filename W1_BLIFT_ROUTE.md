# W1 Deliverable B — Route choice

**Date:** 2026-05-14
**Decision:** **Route 2 — verbatim citation of Hasebe-Saigo 2014 (arXiv:1306.0137,
Nagoya Math. J. 215).**

---

## 1. Comparison

### Route 1 (direct construction in abelian-B)

- Build the lift by fiberwise reduction: B = L^∞(Ω_B, ν) abelian, conditional
  expectation E_B(·)(ω) integrates over residual randomness at accumulator
  history ω.
- Apply scalar HS 2011 / Hasebe monograph Thm 3.26 ω-by-ω, then check
  ω-measurability + integrability to lift.
- Estimated effort: 1-3 days.
- Risk: low (B is abelian; the fiberwise scalar theorem is verbatim, and the
  ω-measurability is clean because Syracuse atoms are ω-measurable by
  construction).
- Output: a project-internal proof.

### Route 2 (verbatim citation)

- Cite Hasebe-Saigo 2014 Theorem 3.4 (moment-cumulant formula) + Proposition
  3.5 (B-extensivity), specialized to abelian B and to the Syracuse
  conditional-expectation setup.
- Estimated effort: 4-8 hours (writing the specialization + checking
  hypotheses match).
- Risk: lower (the theorem is in a published peer-reviewed journal, Nagoya
  Math. J. 215 (2014), with the operator-valued moment-cumulant
  correspondence established at full generality — no commutativity assumption
  on B is invoked).
- Output: pointer to literature + Mode-E hypothesis match.

---

## 2. Why Route 2 wins

1. **It exists verbatim.** Hasebe-Saigo 2014 is a published paper whose
   abstract states explicitly: "we introduce operator-valued monotone
   cumulants ... and show the moment-cumulant formula." Theorem 3.4 is the
   B-amalgamated formula
   `ϕ(X_1 ⋯ X_n) = Σ_{π ∈ M(n)} (1/|π|!) K_π(X_1, ..., X_n).`
   Proposition 3.5 is B-extensivity `K_n^{N.X} = N K_n^X.`
   No project-internal construction is needed.

2. **Hypotheses match Syracuse setup.** HS 2014 §2.1 states the hypotheses:
   B is "a unital algebra" (no commutativity assumption), A is a unital
   algebra containing B, and ϕ : A → B is a B-linear conditional expectation
   with ϕ(b) = b for b ∈ B. Our Syracuse `(A, E_B, B)` from
   AMALG_FREENESS_SETUP.md satisfies all three:
   - A = W*({T_j, M_{b_{[1,k]}}}) is a unital von Neumann algebra
     (hence a unital ∗-algebra)
   - B = W*({M_{b_{[1,j]}}}) is a unital subalgebra (in fact abelian, but
     that's a strict specialization, not a restriction on the theorem)
   - E_B : A → B is the conditional expectation onto B, B-bilinear with
     E_B(b) = b for b ∈ B (standard fact for vN-algebra conditional
     expectations).
   The abelian property of B in the Syracuse case is a strict
   specialization. The theorem applies without modification.

3. **Skeide's monotone-independence-over-B definition (HS 2014 Def 2.2) is
   the operative definition.** Our Syracuse monotone-independence diagnostic
   (Task 1 numerical probe at 10⁶ separation, M_3_alt = 0.1078, plus
   Hasebe-monograph Defn 1.21 peak-rule factorization) matches exactly the
   condition of HS 2014 Def 2.2:
   > "A family of subalgebras (A_λ)_{λ∈Λ} over B is said to be monotone
   > independent over B if `ϕ(X_1 ⋯ X_n) = ϕ(X_1 ⋯ X_{i-1} ϕ(X_i)
   > X_{i+1} ⋯ X_n)` holds for any X_i ∈ A_{λ_i} whenever i satisfies
   > λ_{i-1} < λ_i and λ_i > λ_{i+1}"
   The X̃_j operators (centered off-diagonal corrections at step j) generate
   subalgebras `A_j := B⟨X̃_j⟩_0` over B in the sense of HS 2014 §2.1,
   indexed by j ∈ {1, 2, 3, ...} with the natural total order.

4. **Multivariate version is the version we need.** Our Syracuse application
   uses three operators (X̃_{j_1}, X̃_{j_2}, X̃_{j_1}) at two distinct
   levels j_1, j_2. HS 2014 develops the multivariate theory explicitly:
   X = (X_1, ..., X_r) random vectors, multilinear moments
   `μ_{i_1,...,i_n}^X(b_1,...,b_n) = ϕ(b_1 X_{i_1} b_2 ⋯ b_n X_{i_n})`.
   The repetition of indices (j_1 appearing twice) is built into the
   multi-index `(i_1, ..., i_n)`.

5. **Avoids Mode-H risk.** Route 1 would require constructing the
   ω-measurability + integrability argument from scratch; even with
   abelian-B simplifying things, there is a risk of assuming the
   conclusion at some step. Route 2 has zero Mode-H exposure: the
   theorem is proved in HS 2014; we just match hypotheses.

6. **Effort saving.** 4-8 hours vs 1-3 days, with strictly lower risk.

---

## 3. Disposition

Route 2 chosen. Proceed to W1.C with HS 2014 Theorem 3.4 + Proposition
3.5 as the verbatim source, specialized to Syracuse (A, E_B, B) with
B abelian and the X̃_j family indexed by j ∈ N with the natural total
order.

Route 1 (fiberwise direct construction) is filed as a backup — it would
serve as an independent project-internal proof if for some reason the
HS 2014 citation became unavailable. It is not needed for closure of W1.
