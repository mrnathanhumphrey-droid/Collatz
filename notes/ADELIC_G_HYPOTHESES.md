# ADELIC_G — Chambert-Loir / Tschinkel 2009 (Igusa integrals + adelic volume asymptotics)

**Source:** C:/tmp/adelic/ChambertLoir_Tschinkel_2009_Igusa_Integrals_Adelic_Geometry.txt. Main results Theorems 1.1.1, 1.2.1, 1.3.1.

## Verbatim Theorem 1.2.1 (p. 5-6)

> "Theorem 1.2.1. Assume that v is archimedean.
> If a_v(L,D) > 0, then b_v(L,D) ≥ 1 and there exists a positive real number c such that
>    V(B) ∼ c B^{a_v(L,D)} (log B)^{b_v(L,D) − 1}.
> If a_v(L,D) = 0, then there exists a positive real number c such that
>    V(B) ∼ c B^{a_v(L,D)} (log B)^{b_v(L,D)}."

## Verbatim Theorem 1.3.1 (p. 7)

> "Theorem 1.3.1. When B → ∞, one has the following asymptotic expansion
>    V(B) ∼ [1 / (a(L,D) (b(L,D) − 1)!)] · B^{a(L,D)} (log B)^{b(L,D) − 1} ∫_{X(𝔸_F)} H_E(x)^{-1} dτ_X(x)."

## Verbatim Mellin definition (p. 6)

> "To prove this theorem, we introduce the Mellin transform Z(s) = ∫_{U(F_v)} ∥f_L(x)∥_v^s dτ_{(X,D),v}(x) and establish its analytic properties. We regard Z(s) as an integral over the compact analytic manifold X(F_v) of the function ∥f_L∥_v^s with respect to a singular measure, connecting the study of such zeta functions with the theory of Igusa local zeta functions, see [28, 29]. In particular, we show that Z(s) is holomorphic for Re(s) > a_v(L,D) and admits a meromorphic continuation to some half-plane {Re(s) > a_v(L,D) − ε}, with a pole at s = a_v(L,D) of order b_v(L,D). This part of the proof works over any local field. If v is archimedean (and ε > 0 is small enough), then Z(s) has no other pole in this half-plane. Our volume estimate then follows from a standard Tauberian theorem."

## Hypotheses isolated

- **h1 (GEOMETRIC SETUP):** X smooth projective algebraic variety over a number field F (or local field). D, L effective divisors. ω_{X(D)} = canonical line bundle twisted by D.
- **h2 (HEIGHT FUNCTION):** H_L : U(𝔸_F) → ℝ_+ continuous proper, defined by H_L((x_v)) = ∏_v ∥f_L(x_v)∥_v^{-1}.
- **h3 (MEASURE):** τ_X = Tamagawa-like measure on adelic space X(𝔸_F), built from local metrics on ω_X(D) via the construction of §2.4.
- **h4 (MELLIN TRANSFORM Z(s)):** Z(s) := ∫_{X(𝔸_F)} H_L(x)^{-s} dτ_X(x), or local version Z_v(s) := ∫_{U(F_v)} ∥f_L(x)∥_v^s dτ_{(X,D),v}(x).
- **h5 (ANALYTIC HYPOTHESIS):** Z(s) admits meromorphic continuation past the first pole at s = a(L,D); this is proved by the paper for archimedean v.
- **CONCLUSION:** Volume V(B) = τ_X({H_L ≤ B}) has asymptotic c B^{a(L,D)} (log B)^{b(L,D) − 1}.

## Hypothesis × input check

| Hyp | (1) μ_n | (4) R78 (1+3)^u |
|---|---|---|
| h1 (X smooth projective variety) | FAILED. Syracuse μ_n lives on a *finite group* (Z/3^n)*; not on an algebraic variety. There is no canonical line bundle, no divisor D, no height function. | RELATED — R78 (1+3)^u sits inside ℚ_3*, which is the F_v-points of 𝔾_m as an algebraic variety. So if we view μ_n as living on 𝔾_m(ℚ_3) restricted to the principal-unit coset, we have a *variety* structurally. But it's the ℚ_3-points of a *group* (𝔾_m) — and CLT operates on smooth projective varieties + divisors, not on group ℚ_3-points directly. |
| h2 (height function) | NOT_PROVIDED — no natural height on (Z/3^n)* or even on the principal-unit coset of ℚ_3*. The Mellin transform Z(s) requires H_L : variety → ℝ_+. | If we used the 3-adic norm \|·\|_3 as a "height", then on the principal-unit coset 1 + 3 ℤ_3, the height is constant (everything has norm 1). Doesn't probe interesting structure. |
| h3 (Tamagawa measure) | NOT_PROVIDED | — |
| h4 (Mellin transform Z(s)) | The CLT Mellin Z_v(s) = ∫ ∥f_L∥_v^s dτ is an INTEGRAL OVER THE SPACE; our R75 Mellin is a SUM ON THE PONTRYAGIN DUAL of the finite group. Different objects. | — |
| h5 (meromorphic continuation past first pole) | This IS the analytic statement that, if it held, would close the rate problem. But the hypotheses don't apply, so the statement is just an abstract template. | — |

## Disposition for G

**NO_FIT (categorical).**

CLT's Igusa-integral framework is for *Mellin transforms of height functions on adelic spaces of algebraic varieties*. Syracuse μ_n is a *finite-group measure*, not a height function on a variety. The Mellin transforms involved (Z_v(s) vs μ̂_n(ξ)) are different objects.

The "Mellin transform of height function" is structurally close to the "Igusa local zeta" Z(s; f) = ∫_{F_v} |f(x)|^s dx for f a polynomial, which IS a candidate secondary route from BGT_DISPOSITION. CLT generalizes Igusa to height functions on arbitrary varieties, but Syracuse doesn't supply a variety or a height. The closer relative is **plain Igusa local zeta on (1+3)^u**, which is the BGT-flagged secondary route.

**Structural near-miss:** R78's (1+3)^u algebraic structure IS the kind of polynomial whose Igusa local zeta one could compute. But that would be candidate K (Igusa local zeta, NOT in CLT's adelic-Mellin form), and it would close only F_3 (3-adic factor), not F_∞.

**Adelic factorization tag:** **ADELIC_FACTORIZATION_INHERENT** (Z(s) = ∏_v Z_v(s) factors over places by construction). But this factorization doesn't help: the variety/height hypothesis class is empty for Syracuse.

## Mode H circular fingerprint

Mild: h5 (meromorphic continuation past first pole) is structurally similar to the closure target. CLT's theorem PROVES h5 under their geometric/height hypotheses — that's the substantive content. If Syracuse fit the hypotheses, this would be a closure escape. But it doesn't.
