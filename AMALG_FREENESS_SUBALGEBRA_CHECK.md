# Phase 2 — Validity of B as an Amalgamation Subalgebra

**Date:** 2026-05-14  
**Depends on:** AMALG_FREENESS_SETUP.md

---

## 1. Is B abelian?

**Yes.** The elements of B are multiplication operators M_{b_{[1,j]}} on L²(Ω), where b_{[1,j]}(ω) = Σ_{i=1}^j v_i(ω) are real-valued measurable functions on Ω = (ℕ+1)^∞. Multiplication operators by real-valued functions commute:

> M_{b_{[1,j]}} · M_{b_{[1,k]}} = M_{b_{[1,j]} · b_{[1,k]}} = M_{b_{[1,k]}} · M_{b_{[1,j]}}

Hence B is commutative. This is the simplest possible amalgamation setting in Voiculescu's framework (Voiculescu 1995, §1.1): abelian B is already covered by the scalar free probability framework in the regime where B acts on A from outside, but here B ⊂ A, which is the proper operator-valued setting.

---

## 2. Is B ⊂ A in a natural way?

**Yes.** The T_j operators depend on b_{[1,j-1]} as a parameter in their phase factors. Specifically, from C1_TAO_RECURSION_FORM.md and result_77_T_lead_spectrum.md:

> (T_j f)(ξ) = Σ_{v ≥ 1} 2^{-v} · exp(-2πi ξ · 3^{2j-2} · 2^{-b_{[1,j]}} · [phase] / 3^n) · f(ξ · 2^{-v} mod 3^n)

The accumulator b_{[1,j]} enters T_j as a multiplication-by-function factor in the phase. Therefore M_{b_{[1,j]}} is a natural sub-operation of T_j: it appears as a coefficient (phase multiplier) inside the definition of T_j. This means:

> B = W*({M_{b_{[1,j]}}}) ⊂ W*({T_j, M_{b_{[1,j]}}}) = A

**The embedding is natural**, arising from the factored structure of the Tao recursion: each T_j factors as (phase-multiplication by M_{b_{[1,j-1]}}) composed with a shift-and-weight operator.

---

## 3. Is φ well-defined?

**Yes.** The conditional expectation φ = E_B requires:

**(a) φ(b) = b for all b ∈ B.**  
Since B is the σ-algebra of {b_{[1,j]}}, conditioning onto B is the identity on B-measurable operators: E_B(M_{b_{[1,j]}}) = M_{b_{[1,j]}}. ✓

**(b) φ is B-B-bilinear: φ(b₁ · a · b₂) = b₁ · φ(a) · b₂ for all b₁, b₂ ∈ B, a ∈ A.**  
In the measure-theoretic setting, the conditional expectation E[· | G] for a sub-σ-algebra G satisfies E[g₁ · X · g₂ | G] = g₁ · E[X | G] · g₂ whenever g₁, g₂ are G-measurable. Here g_i = M_{b_{[1,j_i]}} ∈ B, so B-bilinearity follows from the standard tower property of conditional expectation. ✓

**(c) φ maps A to B.**  
φ(T_j) is the conditional expectation of T_j onto the σ-algebra of running sums. Since T_j depends on b_{[1,j]} (via phase) and on v_j (via the geometric weight and shift), φ(T_j) = E_B(T_j) is a B-valued operator: it is a function of b_{[1,j]} only, obtained by averaging out the remaining randomness in v_j given b_{[1,j]}. This is well-defined. ✓

---

## 4. Is the unit shared?

**Yes.** Both A and B contain the identity operator I. The conditional expectation φ satisfies φ(I) = I ∈ B. ✓

---

## 5. Does B admit a trace making φ a trace-preserving map?

The natural trace on A_n (at finite level n) is the operator trace on B(H_n) normalized by 1/dim(H_n). The restriction of this trace to B is the integral with respect to the stationary measure on the Geom(2) paths, which is a well-defined probability measure. φ is trace-preserving in this sense. ✓

---

## 6. Structural remark: B is abelian but not central

Although B is abelian, B is **not** in the center of A. Specifically, T_j and M_{b_{[1,j]}} do not commute in general, because T_j involves a shift ξ ↦ ξ · 2^{-v} that is not in the domain of M_{b_{[1,j]}} (the accumulator acts on the probability space, the shift acts on the frequency space Z/3^n). This non-centrality is essential: if B were central, the B-amalgamated freeness condition would reduce to scalar freeness (trivially failing or holding by different reasons).

---

## Summary

| Check | Result |
|---|---|
| B is abelian | Yes — all M_{b_{[1,j]}} commute |
| B ⊂ A naturally | Yes — M_{b_{[1,j]}} appear as phase-multiplication factors inside T_j |
| φ(b) = b for b ∈ B | Yes — standard conditional expectation identity |
| φ is B-B-bilinear | Yes — tower property of conditional expectation |
| Unit is shared | Yes |
| B is central in A | No — T_j and M_{b_{[1,j]}} do not commute |

**Conclusion:** B is a valid amalgamation subalgebra for the Voiculescu (A, φ, B) framework. All three components of the operator-valued probability space are well-defined and satisfy the hypotheses of Voiculescu (1995) Definition 1.2 and Speicher (1998) Chapter III.
