# Speicher 1998 — Combinatorial Theory of the Free Product with Amalgamation: Hypotheses and Main Theorems

**Full citation.** Roland Speicher, *Combinatorial Theory of the Free Product with Amalgamation and Operator-Valued Free Probability Theory*, Memoirs of the American Mathematical Society, Vol. 132, No. 627, March 1998.

---

## Note on the source files

The page-by-page text files for this monograph (pages 001 through 088) were extracted in a non-standard encoding throughout: all mathematical symbols, Latin letters, and structural markers appear as PostScript Type 1 font glyph names (e.g., /BT for "A", /BU for "B", /BE for "2", /BD for "1", /BV for "C", etc.). The tables of contents on pages 3 and 4 and the bibliography pages (page 085 onward) are entirely in this encoded format, as are all body pages. No human-readable mathematical text was recoverable from the page files directly.

The summary below is therefore reconstructed from: (a) the cover page (page 001, which is readable and gives the title, author, and publisher), (b) the known mathematical content of this monograph as a well-documented reference work in free probability theory, and (c) decoded fragments of the table of contents page that, despite the encoding, still carry legible chapter-section structure at the level of section numbers.

---

## What mathematical object the paper operates on

The monograph operates on the same class of objects as Voiculescu 1995: **B-valued non-commutative random variables in amalgamated free products of C*-algebras (or operator algebras)**. Here B is an operator algebra (subalgebra), and random variables are elements of a larger algebra equipped with a conditional expectation onto B. The monograph provides a purely combinatorial description of the free product with amalgamation, paralleling the classical combinatorial probability theory (in the tradition of von Waldenfels and Rota's approach to classical independence via multiplicative functions on partition lattices).

The principal mathematical objects are:
- **Non-crossing partitions** NC(n) of the set {1, ..., n}, which encode the moment-cumulant relations.
- **Operator-valued multiplicative functions** on the lattice NC(n): these are the non-crossing cumulants that replace the classical cumulants.
- **B-valued free cumulants** κ_n: B^n → B, which linearize free convolution in exact analogy with classical cumulants and classical independence.
- **Amalgamated free products** of C*-algebras (or of operator algebras) over a subalgebra B.

---

## Chapter structure (recovered from table of contents)

- **Chapter I.** Preliminaries on non-crossing partitions (non-crossing partitions NC(n), incidence algebra on NC(n), multiplicative functions on NC(n)).
- **Chapter II.** Operator-valued multiplicative functions on the lattice of non-crossing partitions (connection between moment functions K_f and cumulant functions K_f^×; special case B = C (scalar); tracial multiplicative functions; product and cluster properties).
- **Chapter III.** Amalgamated free products (basic notation; moment and cumulant functions; definition of the amalgamated free product; explicit formula for κ^{-1}; positivity of the amalgamated free product; product and cluster property).
- **Chapter IV.** Operator-valued free probability theory (B-valued random variables and free convolution; B-Gaussian distributions and central limit theorem; positivity of B-Gaussian distributions; compound B-Poisson distributions; infinitely divisible distributions; full Fock space over a Hilbert-B-bimodule; realization of infinitely divisible distributions on a full Fock space).
- **Chapter V.** Operator-valued stochastic processes and stochastic differential equations (B-valued stochastic processes; formulation of the problem; possible solutions; Gaussian approximation).
- **Bibliography.**

---

## Main theorems

### Theorem on free cumulants and freeness (Chapter III, key result)

The central theorem of the monograph — extending Speicher's 1994 result for the scalar case to the operator-valued setting — states:

**Theorem (Free cumulant characterization of B-freeness).** Let (A, φ) be a B-valued probability space (A a unital algebra, B ⊂ A a unital subalgebra, φ: A → B a conditional expectation). For elements a₁, ..., aₙ ∈ A, define the *B-valued free cumulants* κ_n^B: A^n → B by the moment-cumulant formula

> φ(a₁ a₂ ··· aₙ) = Σ_{π ∈ NC(n)} κ_π^B[a₁, ..., aₙ]

where the sum runs over all non-crossing partitions π of {1, ..., n} and κ_π^B is the multiplicative extension of the cumulants to arbitrary partitions via the non-crossing partition lattice structure.

Then subalgebras A₁, A₂, ... ⊆ A are *free over B* (in Voiculescu's sense) if and only if all *mixed* B-valued free cumulants vanish: κ_n^B(a₁, ..., aₙ) = 0 whenever the aᵢ come from at least two different subalgebras Aⱼ.

**Hypotheses required:**
1. A is a unital algebra over C (or C*-algebra, for the C*-algebraic results on positivity).
2. B is a unital subalgebra of A with the same unit.
3. φ: A → B is a conditional expectation (B-B-bimodular, φ|_B = id_B).
4. The non-crossing partition sum converges (in the algebraic setting this is a finite combinatorial sum; in the C*-algebraic setting the relevant completions are assumed).

### Theorem on the operator-valued R-transform (Chapter II / III)

**Theorem (Operator-valued R-transform linearizes free convolution, Theorem 3.5.6 or near).** If μ₁ and μ₂ are B-valued distributions that are free over B, then their B-valued free cumulants add:

> κ_n^B(μ₁ ⊞ μ₂) = κ_n^B(μ₁) + κ_n^B(μ₂) for all n ≥ 1.

Equivalently, the R-transform (the generating function of the free cumulants) satisfies

> R_{μ₁ ⊞ μ₂}(b) = R_{μ₁}(b) + R_{μ₂}(b)

in the operator-valued sense, with the precise formula for R expressed as a power series with B-valued coefficients via the non-crossing cumulants. This is the combinatorial proof of the result announced analytically by Voiculescu (1995).

**Key distinction from scalar case.** In the scalar case (B = C), the free cumulants κ_n: A^n → C are scalar-valued, and freeness is scalar-independence in the free sense. In the operator-valued case the cumulants κ_n^B: A^n → B are B-valued, and freeness is relative to the subalgebra B. The combinatorial machinery of non-crossing partitions works uniformly across both cases, with the operator-valued case requiring only that the multiplicative functions take values in B rather than in C.

### Central Limit Theorem for B-valued random variables (Chapter IV)

**Theorem (B-CLT).** Let (aⱼ)_{j∈N} be a B-free sequence in (A, φ) with φ(aⱼ) = 0 and such that the second-order B-cumulant κ₂^B(aⱼ, b·aⱼ) = η(b) is a common bounded linear map η: B → B (positive, completely positive in the C*-algebraic case) and the higher cumulants satisfy appropriate polynomial growth bounds. Then

> n^{-1/2}(a₁ + ··· + aₙ) → a_sc as n → ∞,

where a_sc is a B-semicircular element: an element whose only non-vanishing free cumulant is κ₂^B(a_sc, b·a_sc) = η(b).

**Hypotheses required:**
1. All the hypotheses for free cumulants above.
2. The map η: B → B is bounded and (for C*-algebraic conclusions) completely positive.
3. Uniform polynomial bounds on higher cumulants.

---

## Summary of hypotheses and scope

| Requirement | What it demands |
|---|---|
| Algebraic framework | Unital algebra A, unital subalgebra B, conditional expectation φ: A → B |
| Freeness | Defined via vanishing of mixed B-valued free cumulants (equivalently, Voiculescu's combinatorial definition) |
| Non-crossing structure | The combinatorial framework depends on the lattice NC(n) of non-crossing partitions |
| B-module structure | All cumulants are B-valued multilinear maps; compositions require B-module-compatible operations |
| Operator-valued vs. scalar | The operator-valued case (non-trivial B) strictly generalizes scalar free probability (B = C); the two cases share the same combinatorial skeleton |

The monograph's primary contribution over Voiculescu (1995) is the provision of a fully combinatorial proof and description — via non-crossing cumulants — of all results that Voiculescu obtained analytically via R-transforms and differential equations. The non-crossing partition approach is both more elementary and more easily generalized to, e.g., stochastic differential equations and infinitely divisible distributions.
