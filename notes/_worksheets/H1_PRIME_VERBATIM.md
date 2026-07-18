# H1' — Verbatim statement of HS 2014 Defn 2.2 and the operative reading for Syracuse

**Date:** 2026-05-14
**Mode:** E — verbatim PDF extraction via `pypdf 6.10.2`, no paraphrase.
**Source:** Hasebe, T. & Saigo, H. (2014). "On operator-valued monotone independence." Nagoya Math. J. 215, 151–167. arXiv:1306.0137v2.
**PDF location:** `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf` (193 KB, 13 pages, magic-byte verified upstream).

---

## 1. Verbatim extraction — page 3, the definition

The operative definition is HS 2014 Defn 2.2 (page 3). The page-2 preamble of §2 ("Operator-valued monotone independence", §2.1 "Preliminary concepts") is required to know the meaning of the symbols. I quote the full preamble verbatim, then the definition.

### 1.1 Preamble — verbatim, p. 2

> "Involutions on algebras are not essential in the scope of this paper, so we do not consider them below. In this paper, B denotes a unital algebra and A a unital algebra containing B as a subalgebra. We assume that the unit of B coincides with that of A. In this paper, algebras can be considered over any commutative field such as R and C. We say that C is a subalgebra of A over B if C is a subalgebra of A and bc ∈ C for all b ∈ B and c ∈ C. A subalgebra of A over B may not contain the unit of A.
>
> For X_1, … , X_r ∈ A, let B⟨X_1, … , X_r⟩_0 denote the subalgebra of A over B consisting of finite sums of elements of {b_1 X_{i_1} b_2 … X_{i_n} b_{n+1} : b_i ∈ B, n ≥ 1, i_1, … , i_n ∈ {1, … , r}}. Note that, in general, B is not contained in B⟨X_1, … , X_r⟩_0.
>
> Let D be another unital algebra containing B as a subalgebra. A map f from A to D is called B-linear if f(b_1 x b_2 + y) = b_1 f(x) b_2 + f(y) for all b_1, b_2 ∈ B and x, y ∈ A. A B-linear map h is called a B-homomorphism if h(xy) = h(x) h(y) for any x, y ∈ A. A B-linear map ϕ with values in B is called a conditional expectation if ϕ(b) = b for b ∈ B. From now on we assume that ϕ is a conditional expectation in the above sense. A triple (A, B, ϕ) is called an algebraic probability space or a noncommutative probability space, as in the case B = C."

### 1.2 Defn 2.2 — verbatim, p. 3

> "The monotone independence over B was introduced by Skeide [12].
>
> **Definition 2.2.** Let Λ be an index set equipped with a linear order <. A family of subalgebras (A_λ)_{λ∈Λ} over B is said to be monotone independent over B if
>
>   ϕ(X_1 … X_n) = ϕ(X_1 … X_{i−1} ϕ(X_i) X_{i+1} … X_n)
>
> holds for any X_i ∈ A_{λ_i} whenever i satisfies λ_{i−1} < λ_i and λ_i > λ_{i+1} (one of the inequalities is eliminated when i = 1 or i = n). Independence for random vectors X_λ = (X_{λ,1}, … , X_{λ,k_λ}), λ ∈ Λ is defined by considering the subalgebras A_λ := B⟨X_{λ,1}, … , X_{λ,k_λ}⟩_0."

### 1.3 Reading

The definition is a **peak-rule factorization through ϕ**. Concretely: in any monomial X_1 … X_n with X_i ∈ A_{λ_i}, identify positions i where the index sequence has a **strict local maximum** (λ_{i−1} < λ_i and λ_i > λ_{i+1}; with the obvious convention at the endpoints). At every such peak, the algebra element X_i may be replaced by its conditional expectation ϕ(X_i) ∈ B without changing ϕ of the product.

Crucially, Defn 2.2 is a closure condition: the substitution can be **iterated**. After one substitution, the resulting word has X_i replaced by ϕ(X_i) ∈ B, which by the B-linearity of ϕ slides through to combine with neighboring B-elements. The reduced word has at most n−1 algebra factors, and the new index sequence may have new peaks, on which the rule applies again. Iterating reduces every word to one whose ϕ can be evaluated by repeatedly pulling out peaks.

This is exactly the form of monotone independence one needs for HS 2014 Thm 3.4 (the moment–cumulant formula via monotone partitions M(n)). Thm 3.4 is established in [HS 2014 §3] **given** Defn 2.2 as the structural hypothesis on the family (A_λ). Verifying Defn 2.2 closes H1'.

### 1.4 What Defn 2.2 does NOT require (audit-relevant)

Reading the quoted text adversarially against the W1 audit caveats:

- **No normality, faithfulness, traciality, or modularity assumption.** ϕ need only be a B-linear, B-fixing conditional expectation; B and A need only be unital algebras over a commutative field. There is no positivity, ∗-structure, or modular operator anywhere in Defn 2.2.
- **No commutativity of B.** The definition is stated over a unital algebra B (potentially noncommutative). The Syracuse setting is the strict specialization B-abelian, well within scope.
- **No size restriction on Λ.** Λ is "an index set equipped with a linear order"; ℕ with the natural order (the Syracuse index set j = 1, 2, …) is a valid choice.
- **No 2-sided B-bimodule structure beyond what is given.** The subalgebras A_λ are "over B" in the sense of HS 2014 p. 2 (closed under left B-multiplication; "B is in general not contained in A_λ"). This is the **left-only** definition. The Syracuse Off_j operators inherit this from the way the prior-accumulator multiplications act on them.

---

## 2. What H1' asserts for the Syracuse X̃_j family

### 2.1 The algebraic data

Fix level n. Take:

- **A:** the unital algebra of bounded operators on ℓ²((ℤ/3^n)\* ) acting on the Tao bilinear pair form (so A contains all multiplication operators M_φ for B-measurable functions φ, all shift operators σ_v: f(ξ) ↦ f(ξ·2^{−v} mod 3^n) for v ∈ ℤ, and all twist operators e^{iθ(b)} for B-measurable phase functions θ(b)).
- **B = B_marginal:** the abelian unital subalgebra of A generated by **multiplication operators in the prior accumulator filtration**:

   `B_marginal = vN({M_{b_{[1, k]}} : k ≥ 0})`

  with the convention `b_{[1, 0]} = 0` (i.e. B_marginal contains the unit). Concretely, every element of B_marginal is a multiplication operator by a measurable function of the running accumulator sequence.
- **ϕ = E_{B_marginal}:** the canonical conditional expectation onto B_marginal — average over all (v_{2k−1}, v_{2k}) ~ Geom(2) iid pairs at the indices that are not pinned by B_marginal's σ-algebra at the point of evaluation. (This is a bona fide conditional expectation in the algebraic-probability sense of HS 2014 p. 2: it is B-linear and fixes B.)

- **X̃_j = Off_j − ϕ(Off_j)** for j = 1, 2, … where Off_j is the off-diagonal correction operator at step j defined in `AMALG_FREENESS_MOMENT_CALCULATION.md §7`:

   `Off_j(f)(ξ) = Σ_{v ≠ v', v, v' ≥ 1} 2^{−v} 2^{−v'} · exp(−2πi · ξ · 3^{2j−2} · 2^{−b_{[1, j−1]}} · (2^{−v} − 2^{−v'}) / 3^n) · f(ξ · 2^{−(v + v')} mod 3^n)`

  X̃_j is the marginal-centered off-diagonal at step j. By construction `ϕ(X̃_j) = 0`.

- **A_j:** the subalgebra of A over B_marginal generated by X̃_j in the HS 2014 sense:

   `A_j := B_marginal⟨X̃_j⟩_0`
   `   = { Σ_finite b_0 X̃_j b_1 X̃_j … b_{m−1} X̃_j b_m : b_i ∈ B_marginal, m ≥ 1 }`

  Per HS 2014 p. 2, A_j does **not** in general contain the unit (it is the *non-unital* subalgebra over B_marginal). This matches Defn 2.2's framing precisely.

### 2.2 The index set and linear order

Take Λ = ℕ with the **natural order** j < j+1 < j+2 < … . The Syracuse iteration index. There is no ambiguity in choice of order: the iteration is forward-in-time, j = 1 is the first step, j = 2 the second, and so on.

### 2.3 H1' — the hypothesis

**H1' (level-graded monotone independence over B_marginal).** The family of subalgebras (A_j)_{j ≥ 1} defined in §2.1, with the index set Λ = ℕ and the natural order, satisfies HS 2014 Definition 2.2: for any sequence of indices (j_1, … , j_n) and any X_i ∈ A_{j_i}, whenever i is a peak in the sense j_{i−1} < j_i > j_{i+1} (with endpoint convention),

   `ϕ(X_1 … X_n) = ϕ(X_1 … X_{i−1} ϕ(X_i) X_{i+1} … X_n)`.

### 2.4 Why this is the right H1'

The W1 audit (caveat 1.2 in `TRACK_A_INTEGRATION.md §1.3`) named the residual hypothesis as monotone independence over `B_marginal`, not over the strict full B = vN(b_{[1, j]}) of `AMALG_FREENESS_SETUP.md §5`. The choice matters: under strict full B, the X̃_j collapse to algebraic zero (the reading-A control in the diagnostic JSON yields M_3_alt ~ 10^{−18}, confirmed). Under B_marginal, the diagnostic is non-trivial (M_3_alt = 0.10783, reading B, sum_entries) and lands in the κ_3^B slot allowed by HS 2014 Thm 3.4 after centering kills κ_1^B (see `W1_BLIFT_VERIFICATION.md §3`).

The marginal B_marginal is the σ-algebra of **prior** accumulator information `b_{[1, k]}` for k ≤ j−1, NOT including b_{[1, j]}. This is the standard level-graded filtration used in noncommutative stochastic-process settings: at step j the operator X̃_j is evaluated against a "what was known before step j" σ-algebra.

### 2.5 Status

H1' is the **single residual hypothesis** between the Track A leading c = 7/45 derivation and a publication-grade unconditional theorem. Per `TRACK_A_INTEGRATION.md §1.4`, closing H1' makes the leading derivation rigorous.

---

## 3. Files

- PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`
- Track A integration: `C:/Collatz/TRACK_A_INTEGRATION.md` §1.4 (residual H1')
- Operator definition Off_j: `C:/Collatz/AMALG_FREENESS_MOMENT_CALCULATION.md` §7
- Numerical anchor: `C:/Collatz/experiments_output/monotone_diagnostic_n3.json`
- W1 sanity check at n=3: `C:/Collatz/W1_BLIFT_VERIFICATION.md`
- Companion low-order checks: `C:/Collatz/H1_PRIME_LOW_ORDER_CHECKS.md`
- Companion structural argument: `C:/Collatz/H1_PRIME_STRUCTURAL_ARGUMENT.md`
- Disposition: `C:/Collatz/H1_PRIME_DISPOSITION.md`
