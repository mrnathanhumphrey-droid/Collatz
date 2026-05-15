# C1 Phase 1 — Tao's exp-sum form for Syracuse μ̂_n(ξ), verbatim

**Date:** 2026-05-12
**Source:** Tao 1909.03562 §7.1 (pp. 32–34 of v2 arXiv PDF). Extracted via pypdf to `C:/tmp/crystal/tao_1909.txt`.

---

## 1. The Fourier coefficient is a renewal-process expectation, NOT a single complete sum

Tao defines χ on Z[1/2] as the additive character (eq 7.1, p. 33):

> "Let χ = χ_{n,ξ} : Z[1/2] → C denote the character
> χ(x) := e^{-2πi ξ (x mod 3^n) / 3^n}                             (7.1)"

The Fourier coefficient of Syrac(Z/3^nZ) at ξ ≠ 0 mod 3 is then **Proposition 7.1** (p. 33):

> "S_χ(n) := E χ( 2^{-a_1} + 3^1 2^{-a_{[1,2]}} + ⋯ + 3^{n-1} 2^{-a_{[1,n]}} )    (7.2)
> obeys the estimate S_χ(n) ≪_A n^{-A} for any A > 0, where … a_{[i,j]} := a_i + ⋯ + a_j."

Here (a_1, …, a_n) ≡ Geom(2)^n. **The "exponential sum" here is an expectation over an n-fold iid Geom(2) random vector**, not a single sum over x ∈ Z/3^n. There are 4^n + lower-order terms when expanded as a sum over (a_1,…,a_n) ∈ ℕ^n weighted by 2^{-(a_1+…+a_n)} — an *exponentially weighted sum over a tuple-index set*.

## 2. Tao's pair-grouping → product over a 2D renewal process

Tao groups adjacent terms (pp. 33–34): for j ∈ [n/2] set b_j := a_{2j-1} + a_{2j} (Pascal(2,1/2) on ℕ+2). Then:

> "S_χ(n) = E [ ∏_{j∈[n/2]} f( 3^{2j-2} 2^{-b_{[1,j]}}, b_j ) ] · g( 3^{n-1} 2^{-b_{[1,⌊n/2⌋]}} )"
> (n odd; g factor dropped for n even)

where (eq 7.4, p. 34):

> "f(x,b) := E ( χ( x(2^{a_2} + 3) ) | a_1 + a_2 = b )     (with (a_1, a_2) ≡ Geom(2)^2)"

and g(x) := E χ(x · 2^{-Geom(2)}).

Then via triangle inequality (eq 7.5):

> "|S_χ(n)| ≤ E ∏_{j∈[n/2]} | f( 3^{2j-2} 2^{-b_{[1,j]}}, b_j ) |"

**Crucial structural feature:** Tao's bound is a **product over a 2-dimensional renewal walk** {(j, b_{[1,j]}) : j ∈ [n/2]}, not a single character sum. The phase function within a single f(x,3) is:

> "f(x,3) = (1/2) χ(5x) + (1/2) χ(7x) = (χ(5x)/2)(1 + χ(2x))     [proof of Lemma 7.2]
> |f(x,3)| = |1 + χ(2x)| / 2 = cos(π θ(j, l))"

with θ(j, l) := { ξ · 3^{2j-2} (2^{-l+1} mod 3^n) / 3^n } (eq 7.8).

## 3. The "inner factor A_v(ξ)" the brief references

The brief's stated form **μ̂_{n+1}(ξ) = Σ_v 2^{-v} A_v(ξ) μ̂_n(2^{-v} ξ mod 3^n)** is the **project-internal one-step recursion** (R75/R77 in `result_77_T_lead_spectrum.md` line 14), derived from peeling off the *first* Geom(2) factor in eq (7.2) and re-expressing the residual as a Fourier coefficient on Z/3^n. **It is NOT directly Tao's working form.** Tao's §7 works with the full product structure (eq 7.5), not the one-step peel.

In the one-step-peel form, A_v(ξ) is the contribution at a_1 = v:

  μ̂_{n+1}(ξ) = E_{a_1} [ χ_{n+1}( 2^{-a_1} ) · μ̂_n^{(a_1)}(ξ) ]
            = Σ_{v ≥ 1} 2^{-v} · χ_{n+1}( 2^{-v} ) · [Fourier coefficient of shifted-rescaled measure at ξ on Z/3^{n+1}]

After the level-(n+1) → level-n reduction (using ξ' = ξ · 2^{-v}? · 3^{-1} mod 3^n? — the precise bookkeeping is in Tao's §6 derivation of Prop 1.14 from Prop 1.17, and in R75's "Tao recursion → diagonal/off-diagonal split"), the inner factor in the project-internal form is:

  A_v(ξ) = phase factor × (level-shift bookkeeping) — at minimum **A_v(ξ) is a single complex number, NOT itself a non-trivial exp-sum**.

In other words: **the project-internal A_v(ξ) is a phase factor / level-bookkeeping coefficient at a fixed v, and the "exp-sum nature" is carried by the outer sum over v (countable, geometrically weighted), not by A_v itself.**

## 4. Answers to brief's Phase 1 (a)-(d)

**(a) Exact form.** Tao's working form is the product (eq 7.5):
  |S_χ(n)| ≤ E ∏_{j∈[n/2]} |f( 3^{2j-2} 2^{-b_{[1,j]}}, b_j )|
The project-internal one-step peel (μ̂_{n+1} = Σ_v 2^{-v} A_v μ̂_n) is a derived recursion, where A_v is a phase coefficient.

**(b) Phase function of A_v(ξ).** In Tao's eq 7.4, the relevant phase is contained in χ( x(2^{a_2}+3) ) with x = 3^{2j-2} 2^{-b_{[1,j]}}. Substituting χ(y) = e^{-2πi ξ (y mod 3^n)/3^n}:

  phase(j, a_2, b) = ξ · 3^{2j-2} · 2^{-b_{[1,j]}} · (2^{a_2} + 3)   mod 3^n

This is **NOT a polynomial in a_2 — it's a 2-adic exponential function (2^{a_2}) viewed mod 3^n.** It is also **NOT a rational function of bounded degree in any of the "summation variables".** It is a *p-adic exponential function of a tuple-index*.

**(c) Summation domain.** The "x" being summed over is **not in Z/3^n at all**. The summation is over tuples (a_1, …, a_n) ∈ (ℕ+1)^n with iid Geom(2) weights 2^{-a_i}. The argument that lands in Z/3^n is the *value* of the linear combination 2^{-a_1} + 3·2^{-a_{[1,2]}} + … mod 3^n, but this is the *image* of the summation, not its domain. There is no multiplicative subgroup structure on the domain.

**(d) Multiplicative character.** χ in Tao's eq 7.1 is **purely additive** (it's e^{-2πi ξ ·/ 3^n}, a group homomorphism Z[1/2] → C^* in the additive variable). **No multiplicative character χ appears anywhere in §7.** The weight structure does not factor as ρ(x) χ(g(x)) for any multiplicative χ.

---

## Summary of category

The Syracuse Fourier coefficient μ̂_n(ξ) is **an expectation of an additive character composed with a 2-adic-linear-combination evaluated mod 3^n**. After pair-grouping it becomes **a product (along a 2D renewal walk) of conditional expectations of an additive character**. The structural ingredients are:

1. **Additive character only.** No multiplicative χ.
2. **Phase function is 2-adic exponential in tuple variables.** Not polynomial / rational of bounded degree in (Z/3^n)^*.
3. **No multiplicative subgroup domain.** The tuple-domain (ℕ+1)^n is unrelated to (Z/3^n)^*.
4. **Product-of-renewal-walk structure**, not single-sum.

This category mismatch is the load-bearing finding for Phase 3.
