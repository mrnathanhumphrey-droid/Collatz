# WATSON_C_SADDLEPOINT_HYPOTHESES — saddle-point / steepest descent, verbatim from Manton DAMTP

**Date:** 2026-05-14. Mode E verbatim extraction.

## Source

Manton, J. (2012). "Asymptotic Methods." Cambridge Mathematical Tripos Part II, Lent 2012.

C:/Users/Nate/OneDrive/Documents/watson_saddle_point/pdfs/Manton_Cambridge_DAMTP_Asymptotic_Methods.pdf, "The Method of Steepest Descent", pp. 10-12.

## Verbatim theorem

> ### The Method of Steepest Descent
>
> In this section, we take the idea of Laplace-type integrals and the stationary phase method one
> step further, by extending the results to integrals in the complex plane of the form
>
> > I(x) = ∫_C f(z) e^{x φ(z)} dz   as  x → ∞
>
> where x is a real parameter. Here, C is some curve in the complex plane, and the functions f(z)
> and φ(z) are analytic in a domain containing C.
>
> Write z = p + iq, with p = Re z, q = Im z, and let φ(z) = u(p, q) + iv(p, q), where the functions
> u, v are real-valued. Since φ is analytic, these functions satisfy the Cauchy–Riemann equations
>
> > ∂u/∂p = ∂v/∂q and ∂u/∂q = −∂v/∂p
>
> which are saying that the gradients of u and v have the same magnitude, but are orthogonal, i.e.
> the curves of constant u and v are at right angles. Note that the gradient of u is perpendicular
> to the curves of constant v, and vice versa.
>
> The 'optimal' contour is one of 'steepest descent' of u: in that case, the integral becomes
> Laplace-like, dominated by the part where u has its maximum. Such a path of steepest descent is
> parallel, but in opposite direction, to ∇u, i.e. parallel to the lines of constant v. Thus the
> phase is constant on such a path.
>
> By Cauchy's theorem, we can deform the contour of integration. The idea of the method of steepest
> descent is to deform the original contour to some other contour on which the phase is constant
> as far as possible.
>
> Paths of steepest descent typically start at saddle points of u, or at infinity. Note that, since
> ∇²u = 0, the function u has no local maxima or minima.

## Saddle-point asymptotic formula (verbatim)

> Let us now determine the basic contribution from a simple saddle. Let C be a steepest descent
> path from a saddle at z_0. Consider the integral
>
> > I_0 = ∫_C e^{x φ(z)} dz   where  φ(z) = α e^{iβ} (z − z_0)²
>
> where α > 0 and −π < β ≤ π. Changing variables to y = −i e^{iβ/2} (z − z_0) moves the contour
> onto the positive real axis. Then
>
> > I_0 = i e^{-iβ/2} ∫_{C'} e^{-x α y²} dy = (i/2) √(π / (xα)) e^{-iβ/2} = i √(π / (2 φ''(z_0) x))
>
> More generally, by Taylor expanding around a general saddle point, we obtain the leading
> asymptotics for a path of steepest descent starting at a saddle point,
>
> > I(x) = ∫_C f(z) e^{x φ(z)} dz ~ i f(z_0) e^{φ(z_0) x} √(π / (2 φ''(z_0) x))   as  x → ∞

## Multi-saddle / endpoint phase mismatch (verbatim)

> If the phases at the endpoints of the integral differ, or if they are different from the phase
> at the saddle point, we need two or three paths of steepest descent, joined up at infinity.

## Stirling formula application (verbatim)

> Γ(x + 1) = ∫_0^∞ t^x e^{-t} dt = ∫_0^∞ e^{x log t − t} dt = x^{x+1} ∫_0^∞ e^{x (log t − t)} dt
>
> where we have changed variables to bring the integral into the general form of integrals that
> we have been considering above. The meromorphic function φ(t) = log t − t has a saddle on the
> real axis, at t = 1, where φ(t) viewed as a real function has its maximum.
>
> The method of steepest descent and the Laplace method are the same in this case. Noting that
> there are two contributions to the integral, one from each side of the saddle, we obtain the
> expansion
>
> > Γ(x + 1) ~ √(2π) x^{x + 1/2} e^{-x}    as  x → ∞

## Hypothesis structure (codified)

Standard saddle-point requires:
1. **An asymptotic parameter x → ∞**
2. **An integrand of the form f(z) e^{x φ(z)} along contour C**
3. **φ analytic with isolated saddle(s) z_0 with φ'(z_0) = 0**
4. **f analytic in a domain containing the saddle**

## Why this is applicable to R78/R79 bilinear sum (with re-parametrization)

The Cochrane Prop 4 p-adic saddle-point at q = p^{r+1} can be viewed as a steepest-descent on a
**p-adic contour**. The "asymptotic parameter" is r (the level), and the saddle is in Z/p (selecting
s* ∈ Z/p from the s ∈ Z/p^r summation range).

For r = 3 (R78.6), the saddle is EXACT — all higher-order corrections vanish mod q. For r ≥ 4,
Hensel lifting introduces sub-leading saddle corrections analogous to the second-derivative
correction in classical steepest descent.

**The hypothesis "analytic φ" maps to "polynomial P_a(s) mod q with controlled p-adic derivative
structure" — analytic-over-Z_p saddle.**

## Multi-saddle for the complex-conjugate-pair PADE picture

When two saddles z_0, z_0* are present (complex-conjugate by reality of f, φ), Manton's prescription
is two paths of steepest descent joined at infinity. Each contributes one Gaussian piece:

```
I(x) ~ i f(z_0) e^{φ(z_0) x} √(π / (2 φ''(z_0) x)) + complex conjugate
     = 2 Re[ i f(z_0) e^{φ(z_0) x} √(π / (2 φ''(z_0) x)) ]
     = 2 |f(z_0)| · |e^{φ(z_0) x}| · √(π / |φ''(z_0) x|) · cos(x Im φ(z_0) + arg f(z_0) + π/4 − (1/2) arg φ''(z_0))
```

For a Padé generating function with complex-pair singularities at z = ρ e^{±iθ}, parameterize
w = log z, so the singularities become w = log ρ ± iθ. The coefficient extraction
ε_k = (1/(2πi)) ∮ f(z) z^{-k-1} dz becomes a saddle-point problem with parameter k, saddle at
z = ρ e^{±iθ}, and Gaussian neighborhood.

The result:
```
ε_k ~ A · ρ^{-k} · k^{-α} · cos(k θ + φ)
```

where:
- ρ = modulus of the conjugate-pair singularity
- θ = argument
- α = singularity exponent (1/2 for square-root, 0 for simple pole, etc.)
- A, φ = constants from the saddle-point integration

## PADE Substitution Check

PADE prediction:
- |ρ_2| ≈ 1.57 at n=13 (transient leading)
- θ_2 ≈ 0.68 rad
- Period 9.2 in n-space ↔ 2π/θ_2 = 2π/0.68 ≈ 9.24 — **matches**.

So multi-saddle Darboux gives:
```
ε_k ~ A · (1.57)^{-k} · k^{-α} · cos(k · 0.68 + φ) + B · (1.016)^{-k} · k^{-β}
```

The **second term** (slow-mode at z=1.016) becomes asymptotically dominant when:
```
(1.016)^{-k} >> (1.57)^{-k}    ⟹    k log(1.57/1.016) >> 0   ⟹    k log(1.55) >> 0
```

This is satisfied for all k → ∞ but the prefactor `(1.57)^{-k} / (1.016)^{-k} = (1.016/1.57)^k`
decays slowly: at k=13, ratio = (0.647)^13 ≈ 0.0024. So the complex pair contributes 0.24% of the
slow-mode at k=13, dominating only if A/B >> 400.

Empirically PADE finds the complex pair dominating at k=13. This means A/B is ~ 400-1000 in
absolute terms, i.e. the transient prefactor is LARGE.

## What does this give for the closure target?

The closure target is a rigorous polynomial-in-A bound on |μ̂_n(ξ)|. The asymptotic for ε_k says
ε_k decays like (1.016)^{-k} · k^{-α} asymptotically — that's a polynomial-in-k bound MULTIPLIED
by an exponential factor.

**If 1.016 > 1 (which it is): ε_k → 0 exponentially.** This means S_n → S_∞ exponentially, which
means |μ̂_n(ξ)| has a positive spectral gap, which IS the closure target.

**Watson/saddle-point applied to f(z) of ε_k, combined with PADE singularity data, predicts an
exponential decay rate of ε_k governed by the closest singularity of f(z) — at z=1.016 asymptotically
or z=1.57 transiently.**

This is the same fact as "the transfer operator's leading subdominant eigenvalue has modulus
< 1 / 1.016 ≈ 0.984."

## What's still missing

A RIGOROUS proof of:
1. f(z) is analytic in |z| < 1 (or some explicit positive radius)
2. The closest singularity is at z = 1.016 (or wherever PADE actually converges to)
3. The singularity is algebraic of known type (so coefficient asymptotic is rigorous, not just heuristic)

For (1) and (3): these depend on transfer-operator theory for the Tao recursion — exactly what
FAURE_DISPOSITION identifies as the missing infrastructure.

For (2): PADE estimates ρ but with finite-n noise; rigorous bound requires either an analytic
argument (transfer operator) or a more refined numerical Hadamard-radius result with explicit
error bars.

This is exactly the **PARTIAL** boundary: technique applies, predicts the right structure, but
rigorous prerequisites for executing the saddle-point on f(z) are the missing transfer-operator
infrastructure.
