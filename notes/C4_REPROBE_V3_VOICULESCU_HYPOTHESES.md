# Voiculescu 1995 — Operator-Valued Free Convolution: Hypotheses and Main Theorems

**Full citation.** Dan Voiculescu, "Operations on certain non-commutative operator-valued random variables," *Astérisque* 232 (1995), pp. 243–275.

---

## What mathematical object the paper operates on

The paper operates on **B-valued non-commutative random variables**, where B is a fixed unital algebra over C (in the analytic sections, a unital Banach algebra). A B-valued random variable is an element a of a larger unital algebra A that contains B as a subalgebra (with the same unit), together with a conditional expectation φ: A → B satisfying φ(b₁ a b₂) = b₁ φ(a) b₂ and φ(b) = b for all b ∈ B. The scalar case B = C recovers the ordinary (scalar-valued) free probability framework developed by Voiculescu in his earlier work.

The paper generalizes the addition and multiplication of free pairs of scalar non-commutative random variables to this operator-valued setting, replacing scalar states with conditional expectations and scalar free products with free products with amalgamation over B.

---

## Core definition: B-freeness

**Definition 1.2 (B-freeness).** Let (A, φ) be as above with B ⊂ Aᵢ ⊂ A for i ∈ I. The family (Aᵢ)ᵢ∈I is called *free* (over B, with respect to φ) if

> φ(a₁ a₂ ··· aₙ) = 0

whenever aⱼ ∈ Aᵢⱼ with i₁ ≠ i₂ ≠ ··· ≠ iₙ and φ(aⱼ) = 0 for 1 ≤ j ≤ n.

A family of subsets or elements is called free if the corresponding generated subalgebras (over B ∪ {element}) are free in the above sense.

---

## Main theorem: the B-valued R-transform and additive free convolution (Theorem 4.9)

This is the central computational result. To state it, one needs the Cauchy transform Gμ and the R-transform R of a B-valued distribution.

**Setup.** For a B-valued distribution μ ∈ SE_B (the set of symmetric B-valued distributions, defined via the symmetric moments μ(S_n(b,...,b))), define the generating series

> Gμ(b) = Σ_{n≥0} μ(b(Tb)ⁿ) (formally in b ∈ B near 0),

and write Gμ(b) = b + b·T_μ(b)·b. Let K_μ be the compositional inverse of Gμ (i.e., K_μ(G_μ(b)) = G_μ(K_μ(b)) = b near 0), and write

> K_μ(b)⁻¹ = b⁻¹ + R_μ(b),

defining the B-valued R-transform R_μ as a germ of a B-valued analytic function at 0 ∈ B.

**Theorem 4.9.** *Assume B is a Banach algebra and μ ∈ SE_B is such that G_μ(b) is analytic in some neighborhood of 0 ∈ B. Let K and R be germs of B-valued analytic functions at 0 ∈ B such that*

> *K(G_μ(b)) = G_μ(K(b)) = b*

*and*

> *K(b)⁻¹ = b⁻¹ + R(b)*

*for b ∈ GL(B) in some neighborhood of 0. Then*

> *R(b) = Σ_{n≥0} SR_{n+1}(μ)(b^{⊗n})*

*where the SR_{n+1}(μ) are given by the canonical element with distribution μ.*

**Consequence (Section 4.11).** If μ₁, μ₂ ∈ SE_B are freely independent B-valued distributions, then

> R_{μ₁ ⊞ μ₂}(b) = R_{μ₁}(b) + R_{μ₂}(b).

That is, the B-valued R-transform linearizes additive free convolution over B, exactly as in the scalar case.

**Precise hypotheses required:**
1. B is a unital Banach algebra over C.
2. μ is a symmetric B-valued distribution for which G_μ is analytic in a neighborhood of 0 ∈ B (equivalently, the symmetric moments define continuous multilinear maps Bⁿ → B with appropriate growth bounds).
3. The two distributions μ₁ and μ₂ are free over B (in the sense of Definition 1.2), not merely scalar-free; their freeness is relative to the subalgebra B.

---

## Multiplication theorem and free exponential (Sections 5.7–5.10)

The paper also treats multiplicative free convolution for symmetric B-valued distributions via a differential equation (the "semigroup for ⊠" equation, analogous to Section 4 for ⊞). The **free exponential map**

> fexp : ∏_{n≥1} SX_n(B) → SE_{B,1}

is shown (Proposition 5.8) to be a bijection, and

> S(μ₁ ⊠ μ₂) is completely determined by Sμ₁ and Sμ₂

(Proposition 5.9), parallel to Theorem 4.9 for addition.

---

## Free Central Limit Theorem (Theorem 8.4)

**Theorem 8.4.** *Assume B is a Banach algebra and (aⱼ)_{j∈N} is a B-free sequence of random variables such that:*

> *1. φ(aⱼ) = 0 for all j ∈ N;*
> *2. there is a bounded linear map η: B → B such that*
>    *lim_{n→∞} n⁻¹ Σ_{1≤j≤n} φ(aⱼ b aⱼ) = η(b);*
> *3. there are constants C_k (k ≥ 1) such that*
>    *sup_{j} ‖φ(aⱼ b₁ aⱼ ··· b_k aⱼ)‖ ≤ C_k ‖b₁‖ ··· ‖b_k‖.*

*Let S_n = n^{-1/2}(a₁ + ··· + aₙ). Then the distribution of S_n converges pointwise to the B-semicircular distribution with canonical form A*(1) + A(η).*

---

## Summary of hypotheses

| Requirement | What it demands |
|---|---|
| Algebraic ambient | A is a unital algebra over C containing B as a subalgebra |
| Conditional expectation | φ: A → B is B-B-bimodular, φ(b) = b |
| Freeness | Subalgebras are free over B (not scalar free): alternating products of centered elements centered on each factor vanish under φ |
| Banach setting (Theorem 4.9) | B is a Banach algebra, G_μ is analytic near 0 |
| Independence | Freeness is over the subalgebra B, permitting dependence via B (this is the key generalization over scalar free probability) |

The central fact is that freeness over B is a strictly weaker condition than scalar freeness when B is non-trivial. Two B-valued random variables can be B-free even though they are not scalar-free; their joint distribution is controlled entirely by their individual B-valued distributions and their B-freeness relation. The R-transform formula then computes the sum.
