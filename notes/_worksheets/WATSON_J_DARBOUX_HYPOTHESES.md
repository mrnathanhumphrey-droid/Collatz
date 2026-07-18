# WATSON_J_DARBOUX_HYPOTHESES — Darboux's method, verbatim from Temme 2013

**Date:** 2026-05-14. Mode E verbatim extraction.

## Source

Temme, N.M. (2013). "Uniform Asymptotic Methods for Integrals." arXiv:1308.1547.

C:/Users/Nate/OneDrive/Documents/watson_saddle_point/pdfs/Temme_2013_Uniform_Asymptotic_Methods_Integrals.pdf, §2.4 "Generating functions; Darboux's method", pages 8-9.

## Verbatim theorem statement and method

> ### 2.4 Generating functions; Darboux's method
>
> The classical orthogonal polynomials, and many other special functions, have generating functions
> of the form
>
> > G(z, w) = Σ_{n=0}^∞ F_n(z) w^n.    (2.14)
>
> The radius of convergence may be finite or infinite, and may depend on the variable z. ...
>
> From the generating function a representation in the form a Cauchy integral follows:
>
> > F_n(z) = (1/(2πi)) ∮_C G(z, w) dw / w^{n+1},    (2.16)
>
> where C is a circle around the origin inside the domain where G(z, w) is analytic as a function of w.
>
> When the function G(z, w) has simple algebraic singularities, an asymptotic expansion of F_n(z)
> can usually be obtained by deforming the contour C around the branch points or other singularities
> of G(z, w) in the w-plane.

## Worked example (Legendre polynomials, verbatim)

> 1 / √(1 − 2 x w + w²) = Σ_{n=0}^∞ P_n(x) w^n,  −1 ≤ x ≤ 1, |w| < 1.    (2.17)
>
> [Singular points: w_± = e^{±iθ} where x = cos θ. Contour deformed to two loops C_± around branch cuts.]
> [Substitute w = w_+ e^s, get contour C_+ in s-plane. Contribution P_n^+(cos θ) is computed via Watson's lemma:]
>
> P_n^+(cos θ) = (e^{-(n+1/2) iθ + iπ/4} / (π √(2 sin θ))) · ∫_0^∞ e^{-ns} f_+(s) ds/√s,    (2.19)
>
> where f_+(s) = √s / (e^s − 1) · (1 − e^{-2iθ}) / (e^s − e^{-2iθ}).    (2.20)
>
> Expanding f_+ in powers of s, we can use Watson's lemma to obtain the large n expansion.
>
> [...] The contribution from the singularity at w_− can be obtained in the same way. It is the
> complex conjugate of the contribution from w_+, and we have
>
> P_n(cos θ) = 2 Re P_n^+(cos θ).

## The key formula for complex-conjugate-pair singularities

From the Legendre example, the COEFFICIENT asymptotic from a generating function with conjugate-pair
algebraic singularities at w_± = ρ e^{±iθ} takes the form:

> a_n ~ 2 Re[ c · ρ^{-n} · n^{-α} · e^{-i n θ} · (correction series in 1/n) ]
>      = 2 |c| ρ^{-n} n^{-α} cos(n θ + arg c) · (1 + O(1/n))

where α is the algebraic-singularity exponent (1/2 for square-root branch points).

## Darboux's method (verbatim continuation)

> The way of handling coefficients of power series is related to Darboux's method, in which again
> the asymptotic behavior is considered of the coefficients of a power series f(z) = Σ a_n z^n.
> A comparison function g, say, is needed with the same relevant singular point(s) as f. When g has
> an expansion g(z) = Σ b_n z^n, in which the coefficients b_n have known asymptotic behavior, then,
> under certain conditions on f(z) − g(z) near the singularity, it is possible to find asymptotic
> forms for the coefficients a_n.

The Darboux method's "comparison function" g is chosen to capture the singularity exactly; the
remainder f − g is regular (or has a milder singularity) and contributes lower-order terms.

## Multi-singularity / coalescing case (verbatim)

Temme §2.4 (end):
> When relevant singularities are in close proximity, or even coalescing, we need uniform methods
> and for a uniform treatment we refer to [44]. In our example of the Legendre polynomials, uniform
> methods are needed to deal with small values of θ; in that case J-Bessel functions are needed.

This is the regime where complex-conjugate singularities approach the real axis. For PADE's
θ≈0.68 rad (~39°), the singularities are NOT close to the real axis, so non-uniform
(standard Darboux) applies.

## Applicability check for the ε_k generating function

Hypotheses required for Darboux on f(z) = Σ ε_k z^k:

1. **f(z) is analytic in some disk around z=0** — YES, ε_k bounded gives radius ≥ some R > 0.
2. **Singularities of f(z) on |z|=ρ are algebraic** (branch points / poles / branch-cuts).
   — REQUIRES VERIFICATION. PADE detection identifies POLE locations; whether the singularity is
   algebraic vs essential is a structural question. If the underlying transfer-operator spectrum
   has discrete eigenvalues, the generating function has SIMPLE POLES at z = 1/λ for each
   eigenvalue λ. That's algebraic (pole order 1).
3. **The complex-conjugate-pair structure produces cos(nθ+φ) modulation** — DERIVED above.

## What Darboux predicts for ε_k

If the generating function f(z) = Σ ε_k z^k has:
- A real singularity at z = ρ_1 (real, asymptotic slow-mode)
- A complex-conjugate-pair at z = ρ_2 e^{±iθ_2}

Then ε_k for large k satisfies:

```
ε_k ~ C_1 · ρ_1^{-k} + 2 |C_2| · ρ_2^{-k} · cos(k θ_2 + φ_2) + O(ρ_3^{-k})
```

The DOMINANT term is the one with smallest ρ — the closest singularity to the origin.

If 1/ρ_1 > 1/ρ_2 (i.e., ρ_1 < ρ_2): real singularity dominates, no oscillation visible.
If 1/ρ_1 < 1/ρ_2 (ρ_1 > ρ_2): complex pair dominates, cos(kθ + φ) modulation is visible.

## PADE picture as Darboux input

From PADE_NUMERICAL_DISPOSITION:
- Real (asymptotic) singularity predicted at z ≈ 1.016 (slow-mode from STATE.md ρ ≈ 0.984)
- Complex pair predicted at |z| ≈ 1.57, θ ≈ 0.68 rad (n=13 transient leading)

These have ρ_1 = 1.016 (asymptotic dominant) and ρ_2 = 1.57 (transient leading, becoming
sub-dominant as k grows).

At k ~ 13: ρ_2 dominates (smaller |z|), oscillation visible.
At k → ∞: ρ_1 dominates (smallest |z|), oscillation faded.

This is exactly Temme's "comparison function g + correction f − g" structure with g picking up
the leading singularity at ρ_2 and the remainder picking up ρ_1.

## Connection to Faure 2009 prediction

Faure 2009 predicts spectral radius ≤ 1/√3 = 0.577, hence the dominant singularity of f(z) at
|z| ≥ √3 ≈ 1.732.

PADE leading at |z|≈1.57 vs Faure prediction √3≈1.732: 10% gap. The Faure prediction is for the
SEMICLASSICAL LIMIT, i.e., asymptotic in r → ∞. PADE at n=13 is in the TRANSIENT regime. Both
predictions converge as n → ∞ to the same "leading singularity at √3."

**Faure prediction's z=√3 is the DARBOUX LEADING SINGULARITY** — the closest pole / branch point
of the transfer operator's generating function.

Wait — closer look: ρ_1 (slow-mode) at 1.016 is CLOSER to origin than √3 ≈ 1.732. So if the
slow-mode is real, it dominates asymptotic over the Faure complex pair. The Faure prediction would
then describe the SUBDOMINANT singularity (the complex pair), not the leading.

Unless the slow-mode at z=1.016 is itself a complex pair very close to the real axis, in which case
Faure's prediction at z=√3 would describe a different (sub-leading) feature.

This needs Phase 2 numerical verification.
