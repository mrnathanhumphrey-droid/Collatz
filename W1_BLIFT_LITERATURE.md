# W1 Deliverable A — B-amalgamated monotone literature pull

**Date:** 2026-05-14
**Mode:** E (verbatim citation discipline; cached PDF read via Read tool;
magic-byte verification documented below)

---

## 1. Newly pulled paper (load-bearing for W1)

**Hasebe, T. and Saigo, H. (2014). "On operator-valued monotone independence."
Nagoya Mathematical Journal **215**, 151-167. arXiv:1306.0137v2 [math.OA],
19 March 2014. 13 pages.**

- arXiv abstract: "We investigate operator-valued monotone independence, a
  noncommutative version of independence for conditional expectation. First
  we introduce operator-valued monotone cumulants to clarify the whole theory
  and show the moment-cumulant formula. As an application, one can obtain an
  easy proof of Central Limit Theorem for operator-valued case. Moreover, we
  prove a generalization of Muraki's formula for the sum of independent random
  variables and a relation between generating functions of moments and
  cumulants."
- PDF source: `https://arxiv.org/pdf/1306.0137`
- Retrieved 2026-05-14 via WebFetch; cached at
  `C:/Users/Nate/.claude/projects/c--As-Above-So-Below-Master/<session-id>/tool-results/webfetch-1778781330594-42hs85.pdf`
  (188.7 KB).
- **Magic-byte verification:** the cached blob was successfully read via Read
  tool as a 13-page PDF document (cat-n rendered page text inline), confirming
  valid `%PDF` header semantics through structural integrity of all 13 pages.
- **Target stable location:** `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`
  (file-move blocked by current sandbox; user can drop the cached copy into the
  closure-hunt folder. The cached path remains valid for re-reading until the
  WebFetch cache is recycled.)

This paper is the principal source for Route 2 (verbatim citation). It is
self-contained for the operator-valued (B-amalgamated) monotone cumulant
theory: B denotes any unital algebra (no commutativity assumption used in the
core theorems), ϕ : A → B is a B-linear conditional expectation with
ϕ(b) = b for b ∈ B, and the moment-cumulant correspondence is proved at full
generality.

**Section-by-section summary:**

- §2.1 — Operator-valued probability space `(A, B, ϕ)` (Definition 2.1
  defines the multilinear functional `μ_{i1,...,in}^X(b_1,...,b_n) =
  ϕ(b_1 X_{i1} b_2 ⋯ b_n X_{in})`).
- §2.2 — Skeide's definition of monotone independence over B
  (Definition 2.2 — verbatim below in W1.C).
- §2.3 — Dot operation `N.X := X^{(1)} + ⋯ + X^{(N)}` for monotone iid
  copies (Definition 2.3), with associativity-up-to-state
  (Proposition 2.4).
- §3 — Monotone partitions M(n) (same combinatorial object as in HS 2011
  / Hasebe monograph; verbatim with Fig. 1). The B-amalgamated multilinear
  functional `A_π` is defined recursively via Speicher's [13] interval-block
  contraction rule. **Lemma 3.1** (polynomiality of `ϕ((N.X_1)⋯(N.X_n))`
  in N), **Remark 3.2** (universality of coefficients), **Definition 3.3**
  (n-th joint cumulant = coefficient of N), **Theorem 3.4 (moment-cumulant
  formula, verbatim):**
  > `ϕ(X_1 ⋯ X_n) = Σ_{π ∈ M(n)} (1/|π|!) K_π(X_1, ..., X_n).`
  **Proposition 3.5 (B-extensivity, verbatim):**
  > `K_n^{N.X} = N K_n^X.`
  **Theorem 3.6** — operator-valued CLT (limit = monotone pair-partition
  formula; arcsine limit recovered in scalar case).
- §4 — Generating functions on `Mul^r⟦B⟧`; modified-composition `F .` G`
  operation (Definition 4.4); **Theorem 4.9 (extended Muraki formula
  for operator-valued monotone independence):**
  > `μ^{X+Y} = μ^X .` μ^Y` for X, Y monotone independent over B.
  **Corollary 4.10** — differential equations
  `(d/dt) μ^{t.X} = κ^X .` μ^{t.X} = μ^{t.X} ⋆ κ^X.`

**What is proved verbatim, with explicit relevance to Syracuse:**
1. (T3.4) The B-amalgamated moment-cumulant formula via M(n).
2. (P3.5) B-amalgamated additivity / extensivity of cumulants.
3. (T3.6) Operator-valued CLT (not needed for W1 directly but anchors
   the framework).
4. (T4.9) B-amalgamated Muraki formula.

**No commutativity assumption on B** is invoked anywhere in §2-§3. The
abelian-B case of Syracuse is a strict specialization in which the
multilinear functionals collapse to commutative pointwise products of
B-valued (= L^∞(Ω_B)-valued) functions; all theorems apply verbatim.

---

## 2. Companion paper (cited but not load-bearing for W1)

**Popa, M. (2008). "A combinatorial approach to monotonic independence over
a C*-algebra." Pacific Journal of Mathematics **237**(2), 299-325.
arXiv:math/0612570v3 [math.OA].**

- arXiv source: `https://arxiv.org/abs/math/0612570`
- Status: cited as reference [10] in Hasebe-Saigo 2014.
- Role: single-random-variable case of the operator-valued theory. HS 2014
  explicitly generalizes Popa's results to multivariate random vectors and
  introduces cumulants; for our Syracuse application, multivariate is
  required (X̃_{j_1}, X̃_{j_2}, X̃_{j_1} at three distinct levels j_1, j_2).
- Pull deferred: not strictly needed for W1 closure, since HS 2014 contains
  the full multivariate B-amalgamated machinery. Recorded here for the lit
  ledger.

---

## 3. What was already in the closure-hunt corpus

- `hasebe_saigo_2011_monotone_cumulants.pdf` (arXiv:0907.4896v3) — scalar
  monotone cumulants. Foundation paper, but state-valued (`φ: A → C`).
- `hasebe_monotone_probability_theory_monograph.pdf` — 131pp monograph.
  §1 monotone independence; §3 scalar monotone cumulants; §8 random-matrix
  applications (asymptotic operator-valued, not closed-form B-amalgamated).
- `muraki_2003_five_independences_kyoto_precursor.pdf` — Muraki's
  five-independence theorem (axiomatic foundations).
- `memoirs.pdf` = Speicher 1998 Memoirs Vol. 132 No. 627 — B-amalgamated
  free probability, the framework that HS 2014 generalizes from "free"
  to "monotone."

The decisive new pull is Hasebe-Saigo 2014.

---

## 4. Disposition

**Route 2 is viable as a verbatim citation route.** Hasebe-Saigo 2014
Theorem 3.4 + Proposition 3.5 are exactly the lifted statements needed
for the Syracuse application (abelian-B is a strict specialization,
no extra adaptation required at the theorem level).

Proceed to W1.B (route choice) and W1.C (rigorous lift statement) using
Hasebe-Saigo 2014 as the verbatim source.
