# Deliverable A — Hasebe-Saigo monotone cumulants, verbatim definitions

**Date:** 2026-05-14
**Status:** Mode E — verbatim extracts. All quotations directly from PDFs at
`C:/Users/Nate/OneDrive/Documents/closure hunt/`. Source-page cited inline.

## Mode-E declaration on operator-valued (B-amalgamated) status

Hasebe-Saigo 2011 (arXiv:0907.4896v3, 13pp) develops monotone cumulants in the
**scalar** algebraic-probability-space setting (`A` a unital ∗-algebra over C,
`φ: A → C` a state). The Hasebe 2024 monograph (131pp) collects the same
material plus extensions; its Section 1 defines monotone independence for
subalgebras of an nc-probability space and its Section 3 redoes the cumulant
theory carefully. **Neither source proves a full operator-valued (B-amalgamated)
moment-cumulant correspondence with B-valued extensivity at every order.**
What is available verbatim is:

1. The scalar moment-cumulant formula via monotone partitions (HS 2011 Thm 6.1;
   Hasebe monograph Thm 3.26).
2. The scalar additivity / extensivity axiom (M3) for monotonically iid sums.
3. The monotone-product factorization of mixed moments for two subalgebras
   (Hasebe monograph Defn 1.15) and several subalgebras (Defn 1.21).
4. The reciprocal-Cauchy-transform composition law `H_{X+Y} = H_X ∘ H_Y` for
   monotonically independent X, Y (Muraki 2003 Thm 4; HS 2011 Defn 5.3).

The B-amalgamated lift is **not yet a verbatim theorem** in the closure-hunt
corpus. The natural lift — replacing `φ: A → C` with `E_B: A → B` and reading
the moment formulas with B-valued products — is consistent with Hasebe
monograph §8 (random-matrix applications of monotone independence, where the
state is a normalized trace and the analog of B is the diagonal subalgebra),
but the monograph there works with **asymptotic** monotone independence in
the large-N limit (Hasebe monograph §8.2), not a closed-form B-amalgamated
cumulant theory at finite scale.

**Operational reading used in Deliverables B-D.** We apply the scalar
Hasebe-Saigo formula "fiberwise" over the B-spectrum: conditioning on the
B-measurable accumulator `b_{[1,j]}`, the residual operator at step j is a
random variable whose distribution is parameterized by b_{[1,j-1]} (the
"prior" accumulator) and which contributes a scalar-valued monotone moment
at each fiber. Marginalizing the scalar formula over the B-distribution
recovers a B-valued (=function-of-accumulator) expectation. This is the
**marginal-centering reading** of `AMALG_FREENESS_MOMENT_CALCULATION.md §8`,
validated by the 2026-05-14 numerical probe (reading B,
`monotone_diagnostic_n3.json`, `M_3_alt = 0.1078` vs reading-A control
`3.3 × 10⁻¹⁸`).

The remaining Mode-E gap is the **rigorous lift of the scalar HS theorem
to the operator-valued setting**. We document where this gap matters in
Deliverable C.

---

## Section 1 — Scalar monotone independence (verbatim)

### HS 2011 §2 Definition 2.5 — Monotone independence (p. 2-3)

> "Let `Λ` be the index set equipped with a linear order `<`. `{A_λ}` is
> **monotone independent** if
> `φ(a_1...a_i...a_n) = φ(a_i) φ(a_1...a_{i-1} a_{i+1}...a_n)`
> holds when i satisfies `λ_{i-1} < λ_i` and `λ_i > λ_{i+1}` (one of the
> inequalities is eliminated when i = 1 or i = n)."

### Hasebe monograph Defn 1.15 — Two-subalgebra version (p. 12)

> "Let (A, φ) be a nc-probability space. Subalgebras B, C ⊆ A are called
> **monotonically independent** if
> `φ(c_0 b_1 c_1 b_2 c_2 ⋯ b_n c_n) = φ(b_1 b_2 ⋯ b_n) φ(c_0) φ(c_1) φ(c_2) ⋯ φ(c_n)` ... (1.5)
> for all n ≥ 1, b_1, ..., b_n ∈ B and c_1, ..., c_{n-1} ∈ C and c_0, c_n ∈ C ∪ {1_Ã}."

### Hasebe monograph Defn 1.21 — Several subalgebras (p. 15)

> "A family of subalgebras (A_i)_{i∈I} of A is called **monotonically
> independent** if for any i_1, i_2, ..., i_n ∈ I and a_1 ∈ A_{i_1}, ...,
> a_n ∈ A_{i_n}, we have
> ```
> φ(a_1 a_2 ⋯ a_n) =
>   { φ(a_k) φ(a_1 ⋯ a_{k-1} a_{k+1} ⋯ a_n)    if 2 ≤ k ≤ n−1, i_{k−1} < i_k > i_{k+1},
>     φ(a_1) φ(a_2 a_3 ⋯ a_n)                  if i_1 > i_2,
>     φ(a_n) φ(a_1 a_2 ⋯ a_{n−1})              if i_{n−1} < i_n. }
> ```

**Remark 1.16(a) (p. 13):** "Monotone independence has an 'asymmetric' nature:
B and C being monotonically independent does not imply that C and B are
monotonically independent."

---

## Section 2 — Cumulant axioms (verbatim)

### Hasebe monograph Defn 3.3 — Monotone cumulants (p. 30)

> "A rule that associates with each nc-probability space (A, φ) and random
> variable x ∈ A a sequence of complex numbers (κ_n(x))_{n≥1} is called
> **monotone cumulants** if the following conditions hold.
>
> **(M1) Polynomiality.** There are universal polynomials P_n^▷(t_1, ..., t_{n−1}),
> n ≥ 2, with convention P_1^▷ = 0, such that for any nc-probability space (A, φ)
> and x ∈ A,
> `κ_n(x) = φ(x^n) + P_n^▷(φ(x), φ(x^2), ..., φ(x^{n−1})),  n ∈ N.`
>
> **(M2) Homogeneity.** For any nc-probability space (A, φ) and any x ∈ A,
> `κ_n(λ x) = λ^n κ_n(x),  n ∈ N, λ ∈ C.`
>
> **(M3) Extensivity.** If N ∈ N and x_1, ..., x_N are monotonically independent
> and identically distributed (monotonically iid) random variables in a
> nc-probability space (A, φ), then
> `κ_n(x_1 + x_2 + ⋯ + x_N) = N κ_n(x_1),  n ∈ N.`"

### HS 2011 Definition 4.5 (p. 5)

> "Let `r_n = r_n(X)` be the coefficient of N in `M_n(N.X)` (or the coefficient
> of t in `m_n(t)`). We call `r_n` the **n-th monotone cumulant** of X."

(Here `M_n(N.X) := φ((X^(1) + ⋯ + X^(N))^n)` for monotonically iid copies. The
notation `r_n` in HS 2011 is the same as `κ_n` in the monograph; we use κ_n
below.)

### HS 2011 Definition 5.3 — Monotone convolution (p. 8)

> "The monotone convolution `μ ⊳ ν` of probability measures μ and ν is defined
> by the relation `H_{μ ⊳ ν}(z) = H_μ ∘ H_ν(z)`, Im z ≠ 0, where H_μ is defined by
> `H_μ(z) = [∫_R μ(dx)/(z − x)]^{−1}`."

---

## Section 3 — Moment-cumulant formula via monotone partitions (verbatim)

### Hasebe monograph Defn 3.23 (p. 39)

> "Let T be a totally ordered finite set. An ordered set partition `π = (ρ, ≤)`
> of T is called a **monotone set partition** if
> - ρ is a noncrossing set partition,
> - if B, B' ∈ ρ satisfies B ⪯ B' then B ≤ B'.
>
> The set of monotone set partitions of T is denoted by `M(T)`."

(Here `B ⪯ B'` is the covering relation: B covers B' if min B ≤ i ≤ max B
for all i ∈ B'. So the linear order on blocks refines the covering partial
order: inner blocks come **higher** in the linear order than outer blocks.)

### HS 2011 Theorem 6.1 / Hasebe monograph Theorem 3.26 (p. 40)

> "On any nc-probability space (A, φ) and for any x ∈ A, we have
> `φ(x^n) = Σ_{π ∈ M(n)} (1 / |π|!) κ_π(x),  n ∈ N.`              (3.13)"

Here `κ_π := κ_{|B_1|} κ_{|B_2|} ⋯ κ_{|B_k|}` for π = (B_1, ..., B_k) (the
product depends only on the underlying set partition ρ, not the linear order).

### Cardinality (Hasebe monograph Prop 3.25, p. 39)

> "The cardinality of M(T) is `(|T| + 1)! / 2`."

So `|M(1)| = 1`, `|M(2)| = 3`, `|M(3)| = 12`, `|M(4)| = 60`.

### Differential recursion (Hasebe monograph Prop 3.15, p. 34)

Let `m_n^▷(t) := φ((x_1 + ... + x_t)^n)` extended from N to R by polynomial
continuation (justified by Prop 3.8). Then
> "`m_0^▷(t) ≡ 1` and for n ≥ 1,
> `d m_n^▷(t) / dt = Σ_{ℓ=0}^{n−1} (ℓ + 1) κ_{n−ℓ}(x) m_ℓ^▷(t)`,
> with initial conditions m_0^▷(0) = 1 and m_n^▷(0) = 0 for n ≥ 1."

### Explicit low-order cumulants (Hasebe monograph Ex. 3.12, eqs 3.3-3.4, p. 34)

> "`κ_1(x) = φ(x);`
> `κ_2(x) = φ(x^2) − φ(x)^2;`
> `κ_3(x) = φ(x^3) − (5/2) φ(x^2) φ(x) + (3/2) φ(x)^3.`"

(κ_1 and κ_2 coincide with classical cumulants; κ_3 differs from the classical
`E[X^3] − 3 E[X^2] E[X] + 2 E[X]^3` by the coefficients `5/2` vs `3` and `3/2`
vs `2` — the gap encodes the broken symmetry of monotone vs commutative
convolution.)

---

## Section 4 — Composition of monotone convolution (verbatim)

### Muraki 2003 Theorem 4 (p. 4)

> "Let X_1, X_2, ..., X_n ∈ A be monotonically independent self-adjoint random
> variables on a C*-probability space (A, φ), in the natural order of
> {1, 2, ..., n}. Then
> `H_{X_1 + X_2 + ⋯ + X_n}(z) = H_{X_1}(H_{X_2}(⋯ (H_{X_n}(z))⋯))`."

### Hasebe monograph Theorem 1.27 (p. 16) — Two-variable shifted version

> "Let (A, φ) be a unital C*-probability space. Suppose that x, y ∈ A are
> monotonically independent. Then for all z ∈ C with |z| < 1/(‖x‖ + ‖y‖) we have
> |M_y(z)| < 1/‖x‖ and
> `M_{x+y}(z) = M_x(M_y(z)),`                                    (1.9)
> where `M_x(z) := z φ((1_A − z x)^{−1}) = Σ_{n≥0} φ(x^n) z^{n+1}`."

### Hasebe monograph eq. (1.10) — Expansion of `φ((x+y)^n)` (p. 17)

> "`φ((x + y)^n) = Σ_{ℓ=0}^n Σ_{k_0+...+k_ℓ = n−ℓ, k_i ≥ 0}
>                  φ(x^ℓ) φ(y^{k_0}) φ(y^{k_1}) ⋯ φ(y^{k_ℓ})`."

This is the **scalar monotone-convolution moment formula** that feeds the
differential recursion. The pattern (monotone-independent y on both sides
of x, with x appearing ℓ times and the y-blocks each factoring) is the
combinatorial content of monotone Wick contraction.

---

## Section 5 — Central limit theorem (verbatim, for orientation)

### HS 2011 Theorem 5.1 (p. 7) / Hasebe monograph Theorem 3.18 (p. 36)

> "Let (A, φ) be a C*-algebraic probability space. Let X^(1), ..., X^(N), ...
> be identically distributed, monotone independent self-adjoint random
> variables with φ(X^(1)) = 0 and φ((X^(1))^2) = 1. Then the probability
> distribution of `X̄_N := (X^(1) + ⋯ + X^(N))/√N` converges weakly to the
> **arcsine law** with mean 0 and variance 1."

(This is the foil to the free CLT, which gives the semicircle; and the
classical CLT, which gives the Gaussian. The arcsine signature is the
diagnostic distribution of monotone convolution.)

### Hasebe monograph eq. (3.7) — Limit moments

> "`lim_{N→∞} φ((X̄_N)^n) = { (n−1)!! / (n/2)!  if n is even, 0 if n is odd. }`"

(For comparison: classical Gaussian gives `(n−1)!!`; free semicircle gives
the Catalan numbers `C_{n/2}`; monotone arcsine gives `(n−1)!! / (n/2)!`.)

---

## Section 6 — What is NOT in the closure-hunt PDFs (Mode-E gap log)

1. **Operator-valued (B-amalgamated) monotone cumulants.** Hasebe-Saigo 2011
   and Hasebe monograph both work with `φ: A → C` (scalar-valued state). The
   monograph §8 uses monotone independence for random-matrix outliers but
   passes to the large-N asymptotic regime, where the auxiliary algebra is
   the diagonal subalgebra of M_N(C) and the analysis is asymptotic, not a
   closed-form `E_B`-valued cumulant theory.

2. **B-extensivity at every order.** The scalar (M3) axiom of HS 2011 gives
   `κ_n(N.x) = N κ_n(x)`. The B-valued analog `κ_n^B(N.X) = N κ_n^B(X)` —
   with `N` a B-valued scalar and `κ_n^B` taking values in B — is the
   conjectural lift used informally in Deliverable C. **No verbatim theorem
   stating this lift appears in the closure-hunt corpus.**

3. **Mixed monotone cumulants.** HS 2011 §6 last paragraph: "this approach
   to cumulants is also applicable to the theory of mixed cumulants for
   classical, free, Boolean and other notions of independence, which will
   be presented in another paper [6]." The reference [6] is "T. Hasebe and
   H. Saigo, in preparation" (status as of 2010). Whether a published
   follow-up exists with explicit mixed monotone cumulants is **not
   resolved** in the corpus.

4. **Operator-valued moment-cumulant formula on monotone partitions over
   B.** Conjecturally,
   `E_B(X^n) = Σ_{π ∈ M(n)} (1/|π|!) κ_π^B(X)`
   should hold with `κ_π^B(X) := κ_{|B_1|}^B(X) ⋯ κ_{|B_k|}^B(X)` where
   the products are taken in B (the abelian accumulator subalgebra in our
   Syracuse setting, so the products are commutative scalar-functions of
   the accumulators). **This formula is not verbatim in the corpus.**

**Operational decision.** Deliverables B and C use the conjectural B-valued
lift but flag it explicitly as Mode-E "conjectural extension of HS 2011
scalar Thm 3.26". The Syracuse case has B abelian (commutative scalar
functions of accumulators), which simplifies the lift considerably: the
"products in B" are commutative pointwise products, and the lift reduces
fiberwise to the scalar HS theorem at each fixed accumulator realization.
The conjectural step is the **uniformity** of the formula across the
B-spectrum and the polynomial structure of P_n^▷ in the B-valued moments.

---

## Files

- HS 2011 PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2011_monotone_cumulants.pdf` (13 pages, arXiv:0907.4896v3)
- Hasebe monograph PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_monotone_probability_theory_monograph.pdf` (131 pages)
- Muraki 2003 PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/muraki_2003_five_independences_kyoto_precursor.pdf` (8 pages)
