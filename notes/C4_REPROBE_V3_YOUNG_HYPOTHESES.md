# Young 1986 — Random Perturbations of Matrix Cocycles: Hypotheses and Main Theorems

**Full citation.** Lai-Sang Young, "Random perturbations of matrix cocycles," *Ergodic Theory and Dynamical Systems* 6 (1986), pp. 627–637.

---

## What mathematical object the paper operates on

The paper operates on **GL(2, R)-valued cocycles over a measure-preserving dynamical system**. Specifically, the base system is a homeomorphism f: X → X of a compact metric space X preserving an ergodic Borel probability measure m. The cocycle is a continuous map A: X × Ω → GL(2, R) (or, in the unperturbed case, A₀: X → GL(2, R)), where Ω is a compact metric space of perturbation parameters. The paper studies what happens to the Lyapunov exponents of the cocycle when the matrix map A₀ is composed at each step with a small random element drawn from a probability measure on the perturbation space.

The main objects are:
- Lyapunov exponents λ₁ ≥ λ₂ of a GL(2, R) cocycle over (X, m, f).
- Oseledec splittings E¹(x) ⊕ E²(x) of R² over X.
- Perturbed cocycles obtained by composing A₀(f^n x) with small random matrices at each step.
- Invariant measures on X × P¹ (X times the projective line) for the associated Markov chain.

---

## Main theorems

### Setup and Lyapunov exponents of the perturbed system

Fix ε > 0. Let Ω be a compact metric space with a distinguished point ω̄, and let A: X × Ω → GL(2, R) be continuous with A₀(x) = A(x, ω̄). For each ε > 0 let ν_ε be a Borel probability measure on Ω. Define the perturbed cocycle by (x, ω) ↦ A(x, ω₀) over the system (f × σ): (X × Ω^Z, m × ν_ε^∞) where σ is the shift. By Oseledec's theorem, this system has Lyapunov exponents λ₁^ε ≥ λ₂^ε.

### Uniform Density Condition (UDC)

The key hypothesis on the perturbation family, stated precisely:

**Uniform Density Condition (UDC).** There exists K > 0 such that for every ε > 0, for m-a.e. x and every u ∈ P¹, the image of ν_ε under the map ω ↦ A(x, ω)u is absolutely continuous with respect to Lebesgue measure on P¹_fx, with density p_ε(x, u) satisfying:

> (1) supp p_ε(x, u) is contained in the ε-neighbourhood of A₀(x)u.
> (2) p_ε(x, u) ≤ K / ε.

In words: the stochastic perturbation acts on P¹ with a density that is supported near the deterministic image of u, and whose height is bounded uniformly in x and u (up to the factor 1/ε).

---

### Theorem 1 (stability of Lyapunov exponents)

**Theorem 1.** *Let A: X × Ω → GL(2, R) and {ν_ε} satisfy (UDC). Then*

> *λ_i^ε → λ_i (i = 1, 2) as ε → 0.*

That is, the Lyapunov exponents of the randomly perturbed cocycle converge to those of the unperturbed cocycle A₀ as the perturbation size tends to zero.

**Precise hypotheses:**
1. X is a compact metric space; m is an ergodic Borel probability measure on X.
2. f: X → X is a homeomorphism preserving m.
3. A₀: X → GL(2, R) is continuous.
4. (A, {ν_ε}) satisfies the Uniform Density Condition.
5. The perturbation space Ω is compact metric with a base point ω̄, and A is jointly continuous.

Note: the paper explicitly treats the 2-dimensional case (GL(2, R)) because the proof method — controlling invariant measures on X × P¹ via density estimates — is specific to dimension 2. The author notes that much of the framework (Lemma 2.2, the argument with invariant measures on X × P¹) applies to GL(n, R), but the decisive estimates are carried out only for n = 2.

---

### Theorem 2 (stability of Oseledec splittings)

**Theorem 2.** *Under the same hypotheses as Theorem 1, assuming additionally that λ₁ > λ₂ and ε sufficiently small (so λ₁^ε > λ₂^ε also), let R(x, ω) = E¹(x, ω) ⊕ E²(x, ω) denote the Oseledec decomposition of the perturbed cocycle over (x, ω) ∈ X × Ω^Z. Let δ > 0 be given. Then for all sufficiently small ε > 0, for (m × ν_ε^∞)-a.e. (x, ω),*

> *lim_{n→∞} (1/(2n+1)) · #{|k| < n : ∠(Eⁱ(f^k x, σ^k ω), Eⁱ(f^k x)) < δ, i = 1,2} > 1 − δ.*

In words: for small ε, the perturbed Oseledec subspaces are within angle δ of the unperturbed Oseledec subspaces for a fraction > 1 − δ of the time steps.

**Additional hypothesis (beyond those of Theorem 1):** λ₁ > λ₂ (the unperturbed exponents are distinct, so the Oseledec splitting exists and is non-trivial).

---

### Theorem 3' (continuity of exponents under perturbation of f, for fixed ε > 0)

**Theorem 3'.** *Let M be a compact oriented surface with Riemannian metric and Riemannian measure m. For f ∈ Diff^1_m(M) (C¹ diffeomorphisms of M preserving m), consider the Lyapunov exponents λ_i^ε(f) of Df composed with ε-rotations. For fixed ε > 0, the map f ↦ λ_i^ε(f) is continuous (i = 1, 2).*

This contrasts with the ε = 0 case, where Mañé's theorem states that λ₁(f) = λ₂(f) = 0 on a residual set outside the Anosov components.

**Additional hypotheses:** M is a compact oriented surface (dim 2); f is C¹ and measure-preserving; ε > 0 is held fixed.

---

### The special case: rotation perturbations (Theorems 1', 2')

In the simplest concrete setting (stated as Theorems 1' and 2' before the general versions), the perturbation consists of composing A₀(f^n x) with an independent rotation R_{ω_n}, where ω_n are i.i.d. uniform on [−ε, ε]. The UDC is easily verified in this case. Theorem 1' gives λ_i^ε → λ_i, and Theorem 2' gives the stability of the Oseledec splitting in the same sense.

---

## Summary of hypotheses

| Requirement | Content |
|---|---|
| Dynamical base | (X, m, f): compact metric space, ergodic Borel probability measure, homeomorphism |
| Cocycle | A₀: X → GL(2, R) continuous; perturbations A: X × Ω → GL(2, R) continuous |
| Perturbation class | Satisfies UDC: density of the perturbed action on P¹ is supported near A₀u, bounded by K/ε |
| Dimension | Strictly 2 (the density estimates are specific to GL(2, R) / the projective line P¹) |
| For splitting stability | λ₁ > λ₂ (non-degenerate Oseledec splitting) |
| Independence | The perturbations ω_n are i.i.d. (or, in the general version, form a product measure on Ω^Z) |

---

## Remarks on near-independence

Young's UDC explicitly allows the perturbed distributions ν_ε to depend on both x (the base point) and the current direction u ∈ P¹; the constant K in the density bound is uniform across all x and u. The key feature that the proof exploits is not exact independence between different time steps, but rather the absolute continuity of the conditional distributions on P¹ with a uniform density bound. The argument works step by step through the Markov chain on X × P¹, not through any joint independence assumption at multiple time steps simultaneously.

In particular, the theorem does not assume that the matrices A(x, ω) and A(f^m x, ω') at different times are independent of each other as random variables (they can share the same underlying randomness Ω); what is assumed is that at each step the *conditional* distribution of A(x, ω)u on P¹_fx is absolutely continuous with bounded density. This is a substantially weaker requirement than scalar independence.
