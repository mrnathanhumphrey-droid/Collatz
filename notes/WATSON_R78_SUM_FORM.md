# WATSON_R78_SUM_FORM — extraction of the load-bearing sum and technique routing

**Date:** 2026-05-14. Mode E. Pre-registration locked.

## Two distinct asymptotic objects in the project

The probe targets BOTH of these (related but distinct):

### Object 1 — R78/R79 bilinear sum (LEVEL-r object; the per-r quantity)

```
T_p(r) := Σ_{a ∈ supp_r} 1̂_r(p·a) · ψ_r(a),   |supp_r| = p^{r-1}
```

with `1̂_r(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u)`, `ψ_r(a) = G_r(a)/√q = e_q(P_a(s*(C_a)))` at r=3 (R78.6).

Empirical (R79b): |T_3(r)| ∝ N^{0.522 ± 0.008} at r=8..20.

Asymptotic regime: **r → ∞ at fixed p** (here p=3). Parameter: N = p^{r-1}.

### Object 2 — ε_k generating function (LEVEL-k object; the c=7/45 deviation sequence)

```
f(z) := Σ_{k=2}^{∞} ε_k z^k,   ε_k := S_k − S_∞ − corrections (deviations from c=7/45 fit)
```

Numerical (PADE): ε_k k=1..13 with sign pattern `+ + − − − − − − − + + + +`, leading singularity
at |z|≈1.57 at n=13, slow-mode predicted at z≈1.016 (asymptotic), complex-conjugate pair at
θ≈0.68 rad period 9.2 in n-space.

Asymptotic regime: **k → ∞**. Parameter: k.

### Critical observation — these are different asymptotic regimes

T_p(r) is a per-r bilinear quantity at p=3; ε_k is a per-k coefficient of a generating function
whose singularity structure governs ε_k's large-k decay. They are linked through:

```
S_k = M_k(1) = Σ_ξ |μ̂_k(ξ)|²    (R75 Plancherel)
ε_k ≈ S_k − S_∞    (subdominant deviation; R76/R77 reduce to off-diagonal contribution)
```

The off-diagonal contribution to ε_k at level k is bounded by polynomial-in-q saving on T_p(r) at
level r ~ k. So a polylog-strength bound on T_p ↔ a polylog-strength bound on ε_k.

## Technique routing decision

Two-pronged. Apply different techniques to the two objects:

### Object 1 (T_p sum at level r) → **p-adic saddle-point (technique G) + multi-saddle (F)**

The saddle-point is the natural fit. Cochrane Prop 4 + truncated p-adic log + saddle is already
the existing machinery at r=3. The question is asymptotic-in-r behavior.

The "asymptotic parameter" for saddle-point is NOT the standard `x→∞ in ∫ e^{x φ(z)} dz`. It is the
LEVEL r, and the saddle is a finite-set saddle (s ∈ Z/p^r) selecting a discrete saddle s*(a) ∈ Z/p.
The saddle is degenerate in the classical-analysis sense (no second-derivative-controlled Gaussian
neighborhood). It is regular in the p-adic sense (`v_p` of dP/ds determines the saddle's "weight").

### Object 2 (f(z) generating function for ε_k) → **Darboux's method (technique J)**

Temme §2.4, verbatim:

> When the function G(z, w) has simple algebraic singularities, an asymptotic expansion of F_n(z)
> can usually be obtained by deforming the contour C around the branch points or other singularities
> of G(z, w) in the w-plane.

Plus (cross-reference Flajolet-Sedgewick Ch. VI as recommended in FAURE secondary routing):

> The way of handling coefficients of power series is related to Darboux's method, in which again
> the asymptotic behavior is considered of the coefficients of a power series f(z) = Σ a_n z^n. A
> comparison function g, say, is needed with the same relevant singular point(s) as f. When g has
> an expansion g(z) = Σ b_n z^n, in which the coefficients b_n have known asymptotic behavior, then,
> under certain conditions on f(z) − g(z) near the singularity, it is possible to find asymptotic
> forms for the coefficients a_n.

And Manton on multi-saddle: when there is a complex-conjugate pair of saddles, each contributes
i √(π/2φ''(z0)x) · f(z0) e^{φ(z0)x}, and the total coefficient gets a `2 Re(...)` factor producing
cos(nθ + φ) modulation. This is the **exact predicted structure** from PADE (θ≈0.68 rad period 9.2).

## Phase 2 will execute these two threads in parallel

Thread α — Object 1: apply p-adic saddle-point to T_p(r) at general r (Hensel-lifted) and derive
the rate of decay of T_p(r) / N as r → ∞.

Thread β — Object 2: apply Darboux + multi-saddle to f(z) using PADE-fitted singularity structure
(real z≈1.016 + complex pair at z=|z|·e^{±iθ}), derive asymptotic for ε_k at large k.

Cross-check: does the Object-1 rate match what Object-2 predicts on the same ε_k?

## Why the bilinear bound IS the saddle-point question

The R78.4-78.6 derivation already used saddle-point at r=3. The "what's missing for general r"
(per R78.6 Status: RIGOROUS at r=3 only) is precisely the Hensel-lifted saddle correction at r≥4.

This means: **the saddle-point technique has ALREADY BEEN PARTIALLY APPLIED** in the project. The
new question is the asymptotic-in-r behavior, which is what Phase 2 Thread α addresses.

## Reference: standard saddle-point asymptotic

Manton verbatim, p.11:
> By Taylor expanding around a general saddle point, we obtain the leading asymptotics for a
> path of steepest descent starting at a saddle point,
>
> I(x) = ∫_C f(z) e^{x φ(z)} dz ~ i f(z0) e^{φ(z0)x} √(π / 2 φ''(z0) x)  as  x → ∞

For multi-saddle (Manton p.11):
> If the phases at the endpoints of the integral differ, or if they are different from the phase
> at the saddle point, we need two or three paths of steepest descent, joined up at infinity.

For Darboux's method on coefficient asymptotics, the singularity at |z|=ρ contributes ~ ρ^{-k}
times a polynomial in k. A complex-conjugate pair at z = ρ e^{±iθ} contributes 2 Re(c · ρ^{-k} e^{-ikθ})
= 2 |c| ρ^{-k} cos(kθ + arg c), the exact PADE-observed modulation pattern.

## Watson lemma — assess applicability

Cohn Nebraska notes, verbatim:
> if f(t) ~ Σ a_n t^{α + βn} as t → 0+, ∫_0^b f(t) e^{-x t} dt ~ Σ Γ(α + βn + 1) a_n x^{-(α+βn+1)}
> as x → ∞.

To apply Watson to a SUM over Z/p^r, we would need to either:
(i) Re-express the sum as a Laplace-transform integral, OR
(ii) Apply Watson INDIRECTLY via a generating-function f(z) for the coefficients.

Option (ii) is exactly Darboux's method with the Gamma function-converting step — and that's the
right tool. Direct Watson on the bilinear sum has no natural Laplace structure (no positive parameter
multiplying a continuous variable).

**Watson NO_FIT to Object 1. Watson SUPPORT (via Gamma functions in Darboux output) to Object 2.**

## Stationary phase — assess applicability

Manton p.8 verbatim:
> Generalised Fourier Integrals: The Stationary Phase Method ... case, the integral is said to have
> a stationary phase at c.

For T_p, the sum is `Σ_{a in supp} 1̂(p·a) · e_q(P_a(s*(C_a)))`. Two pieces of phase: the 1̂
geometric kernel and the saddle phase P_a(s*). The "stationary phase" technique applies if there
is a slowly-varying amplitude + rapidly-oscillating phase. Here:
- e_q(P_a(s*)) varies as a function of a (q=p^{r+1} is large; phase is moderate)
- 1̂(p·a) has its own oscillation structure

This IS stationary phase, but on a discrete sum, not a continuous integral. The discrete-stationary-
phase analog is **Poisson summation** combined with stationary phase. The bilinear inner-Plancherel
in PATH2_BILINEAR uses exactly this (Poisson Σ → frequency-space sum).

**Stationary phase ≡ Plancherel inner step (already applied in PATH2_BILINEAR).**

## What about Mellin-Barnes?

Temme §2.5 verbatim:
> Mellin convolution integrals are of the form F_λ(x) = ∫_0^∞ t^{λ−1} h(xt) f(t) dt, and they
> reduce to the standard form of Watson's lemma when h(t) = e^{-t}. For the general case ...
> A main step in the method to obtain asymptotic expansions of the integral in (2.21) is the use
> of Mellin transforms and their inverses. These inverses can be viewed as Mellin-Barnes integrals.

Mellin-Barnes converts coefficient ↔ singularity. The PADE picture is precisely a singularity-of-f(z)
picture. Mellin-Barnes IS the natural tool for extracting ε_k asymptotic from f(z) singularity data.

**Mellin-Barnes APPLICABLE to Object 2. Synonym with Darboux for coefficient asymptotic via
contour over generating function in this discrete-coefficient setting.**

## Conclusion of routing

- **Object 1** (T_p bilinear sum, level-r): saddle-point already used at r=3, generalization is
  Hensel-lifted saddle. This is what R78 path 4 has flagged as open. Phase 2 Thread α will set up
  the Hensel-lifted saddle and derive its asymptotic order in r.

- **Object 2** (ε_k coefficient sequence, level-k): Darboux/Mellin-Barnes on f(z) using PADE
  singularity structure. Phase 2 Thread β will extract the leading + sub-leading + complex-pair
  asymptotic for ε_k.

Both threads will be executed and compared in Phase 3 against:
- ε_k numerical data k=8..13
- PADE singularity structure (z≈1.016, z=1.57 transient, complex pair)
- Faure prediction (radius √3)
