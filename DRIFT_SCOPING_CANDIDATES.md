# DRIFT_SCOPING_CANDIDATES — Phase 2, candidate Foster–Lyapunov functions for Syracuse

## State space: which Markov chain?

Two natural Markov-chain formulations of "Syracuse":

**Chain A (lift, on Z⁺):** The Collatz/Syracuse trajectory X_n ∈ Z⁺ under T(x) = (3x+1)/2^{v_2(3x+1)} for odd x. (For even x, X_{n+1} = x/2.) Irreducibility presumed (Collatz conjecture). Stationary distribution: open question — almost certainly **no proper stationary distribution exists on Z⁺** because typical orbits descend to {1, 2, 4} (Tao's "almost bounded values"), so the chain is either transient on the absorbing cycle or absorbing.

**Chain B (truncation, on Z/3^n Z \ {0}):** The Markov chain with transition kernel
> K_n(r → s) = P[(3r + 1) · 2^{−v} ≡ s (mod 3^n)], v ~ Geom(1/2)
This is the chain whose stationary measure is the Syracuse μ_n that appears in Tao Prop 1.17. **Finite state space** (|supp π_n| = 2 · 3^{n−1}). **Irreducible, aperiodic, finite ⇒ unique stationary π_n.** This is the relevant chain for the c = 7/45 question.

**The polynomial-in-A Fourier bound needs π̂_n(ξ) decay for π_n on Chain B.**

## Reversibility

**Chain B is NOT reversible** in general. The transition r → s involves multiplication by 2^{−v} mod 3^n with v ~ Geom(1/2), then translation by 2^{−v} (3 r + 1)/3 = (r + 2^{−v}/3) mod 3^n. The forward kernel K(r, s) ≠ K(s, r) generically because:
1. The 3-adic valuation v_2 is not preserved under reversal (going s → r requires inverting (3r+1) · 2^{−v} mod 3^n, picking up a different v distribution).
2. The class structure r ∈ {1, 2} mod 3 evolves under T in a one-directional way (v even ↦ class +, v odd ↦ class − per R66 chain rule); the reverse transition requires a different stochastic structure.

**Consequence.** Taghvaei–Mehta 2005.08145 Theorem 1 (Poincaré inequality from drift) requires reversibility (Assumption 1). For Chain B, must use Proposition 2 (non-reversible extension): need BOTH K_n and its L²(π_n) adjoint K_n† to satisfy (4)–(5) with the SAME V, K, λ, b, α. This is a **strictly stronger** hypothesis.

## Candidate Lyapunov functions

### (C1) V(x) = 1 + dist(x, 0) (linear-growth on Z⁺ lift)

**Only applicable to Chain A (Z⁺). Not applicable to Chain B (finite state space).**

On Chain A: K V(x) = E[V(T(x))] for x odd large. (3x+1)/2 ≈ 3x/2 with prob 1/2 (v=1), (3x+1)/4 ≈ 3x/4 with prob 1/4 (v=2), etc. Geometric mean of multiplier: E[log(3/2^v)] = log 3 − log 2 · E[v] = log 3 − 2 log 2 ≈ 1.099 − 1.386 < 0, so contraction in **log**-norm.

For V(x) = 1 + x: K V(x) ≈ 1 + (3/2)x ≥ V(x), so NO drift in linear V on Chain A. Even at v=1 (prob 1/2), (3x+1)/2 ≈ 1.5 x, expansion. Drift inequality fails for linear V.

For V(x) = 1 + log(1 + x) (logarithmic): K V ≈ V + log(3) − 2 log(2) = V − 0.288 + lower order, satisfying drift in log-V (modulo small-set details). **Logarithmic drift exists on Chain A** under standard Collatz heuristic.

**Obstruction for c = 7/45 question.** π for Chain A doesn't exist as a probability measure (orbit absorbs to {1,2,4}). The "stationary distribution" is the trivial absorbing measure on the 1-cycle. This is the **wrong chain** for the Fourier bound.

### (C2) V(x) = 1 + v_3(x) on (Z/3^n Z)\* (3-adic valuation)

**Applicable to Chain B.** But on Chain B's state space (Z/3^n Z)\ {0}, v_3 is bounded above by n. So V: state-space → [1, n+1] is bounded. Bounded V on a finite state space is trivially Foster-Lyapunov: any constant or bounded function satisfies a drift inequality (since K V ≤ ‖V‖_∞ ≤ c for some c). **Drift inequality holds but is uninformative.**

The minorization condition (5) on Chain B with finite state space: for any small set K ⊂ Z/3^n Z, minorization with α > 0 holds automatically by Doeblin's condition on a finite state space (positive density Markov chain on finite state).

**Conclusion.** On Chain B, every bounded V satisfies drift; Doeblin's condition gives a spectral gap with constant **depending on n** (the state-space cardinality). The constants λ, b, α, ν(K) in Taghvaei–Mehta Theorem 1 will scale with n, giving β = λ/(1 + 2b/α) → 0 as n → ∞ at some rate **inherited from the chain's mixing time on (Z/3^n Z)**.

This is structurally the same problem as the L²-flattening Translation A failure (Probe 1, §2.2): the "uniform-in-n" spectral gap is what's needed for a polynomial-in-A bound, and the framework supplies a gap **with n-dependent constant** rather than a uniform gap.

### (C3) V(x) = 1 + log(1 + dist(x, 0)) on (Z/3^n Z) (logarithmic)

Same finite-state-space issue as (C2): bounded V on a finite state space, trivially drift-satisfying with **n-dependent constants**. Same obstruction.

### (C4) V(x) = −log π_n(x) (entropy / Csiszár-style)

V(x) = −log π_n(x) where π_n is the stationary distribution itself. Self-referential — would require knowing π_n to construct V. Even setting aside circularity:
- K V(x) = −Σ_y K(x, y) log π_n(y). The drift K V ≤ (1 − λ) V + b · 1_K is exactly the **modified log-Sobolev inequality** on Chain B with constant λ controlling entropy decay rate.
- Modified-log-Sobolev gives stronger mixing than spectral gap, but **still doesn't give Fourier decay of π_n** — it gives KL-divergence decay of P^t μ_0 toward π_n.

Same fundamental obstruction.

## Pre-conclusion (Phase 2)

On Chain B (the relevant chain for the c = 7/45 question), **every reasonable Foster–Lyapunov candidate is bounded** (finite state space). Drift inequality holds trivially. The framework's output is:
- Spectral gap β_n > 0 for each fixed n
- β_n's **dependence on n is NOT bounded below by the framework's constants** without additional structure

The constants (λ, b, α, ν(K)) inherent to (4)–(5) for any bounded V on (Z/3^n Z)\ {0} have n-dependence inherited from chain irreducibility/aperiodicity, not from a uniform geometric property. Concretely:
- α (minorization on K via Doeblin) decays at least as 3^{−n} (since π_n(K) ≤ 1/|state space| · |K| ≤ |K|/(2·3^{n−1}) for natural choices of K).
- ν(K) is the minorization measure mass on K, similarly 3^{−n}-scaled.
- The resulting β = λ/(1 + 2b/α) decays at least exponentially in n.

**The framework gives an n-dependent (exponentially small) spectral gap, not a uniform-in-n polynomial-in-A bound.**

## Reversibility check (adversarial A2)

Even granting a candidate V with usable drift, Taghvaei–Mehta Theorem 1 (the cleanest "drift → Poincaré → spectral gap" theorem) requires reversibility (Assumption 1). Chain B is not reversible. Must use Proposition 2 (non-reversible), which requires K_n† to satisfy (4)–(5) with the SAME V. This is **strictly harder** and the proposition's β_+ bound (eq. 15) is in terms of ‖K_n‖_{L²₀(π_n)} for the operator K_n on functions, not for the kernel reversal.

The non-reversibility is a **technical obstruction** stacked on top of the structural one (n-dependent constants).

## State-space caveat (adversarial A3)

The natural Syracuse Markov chain on Z (Chain A) has NO proper stationary measure (orbits absorb). The relevant chain for the Fourier question is Chain B on Z/3^n Z, which is **a different Markov chain at each level n**. The drift framework applies separately to each Chain B at level n; uniform-in-n constants do not follow from any single application.

To get uniform-in-n constants, one would need a "limit chain" on Z_3 (3-adic integers). On Z_3, Chain A's projection T̄ does have a natural Markov-chain interpretation under Geom(1/2)-distributed v_2, and a stationary measure μ_∞ on Z_3 exists (the limit of μ_n under appropriate normalization). But μ_∞ on Z_3 is **3-adic**: its "Fourier coefficients" are characters of Z_3, which are precisely the level-n characters of Z/3^n Z. The framework's output for the Z_3 chain, even if uniform constants were obtainable, is still TV/L²-mixing, NOT Fourier decay of μ_∞.

## Summary of Phase 2

- Chain A on Z⁺: drift functions exist (log-V works), but no proper stationary measure for the c = 7/45 question.
- Chain B on Z/3^n Z: bounded V's trivially satisfy drift, but the framework's spectral-gap constants depend on n with no uniform-in-n control.
- Chain on Z_3 (limit): Foster–Lyapunov could in principle give uniform constants, but the conclusion (L²/TV mixing) is the wrong-flavor bound — not Fourier decay of μ_∞.

In each case, EITHER the chain has a usable stationary measure but trivial uniform-in-n spectral gap, OR the chain has the right stationary measure but the framework's output is TV/L²-mixing, NOT Fourier decay.
