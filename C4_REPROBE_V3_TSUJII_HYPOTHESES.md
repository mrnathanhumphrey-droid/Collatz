# Tsujii 2010 — Quasi-Compactness of Transfer Operators for Contact Anosov Flows: Hypotheses and Main Theorems

**Full citation.** Masato Tsujii, "Quasi-compactness of transfer operators for contact Anosov flows," arXiv:0806.0732v3 (8 April 2010). Published version: *Ergodic Theory and Dynamical Systems*.

---

## What mathematical object the paper operates on

The paper operates on **transfer operators (also called Koopman operators or Ruelle–Perron–Frobenius operators) associated to contact Anosov flows on closed smooth manifolds**. Specifically:

- The phase space M is an orientable closed (compact, no boundary) C^r manifold of dimension 2d+1, equipped with a C^r contact form α.
- The flow F^t: M → M is a C^r contact Anosov flow: it preserves α, and the tangent bundle splits as TM = E^c ⊕ E^s ⊕ E^u (center = flow direction, stable, unstable), with exponential contraction/expansion on E^s and E^u.
- The **transfer operator** is L^t(u)(z) = u ∘ F^t(z), acting initially on C^r(M).

The paper constructs a scale of **Hilbert spaces B^β** (for 0 < β < (r−1)/2) embedded between Sobolev spaces:

> W^s(M) ⊂ B^β ⊂ W^{-s}(M) for s > β, and W^s(M) ⊃ B^β for s < -β.

These are spaces of distributions on M, constructed via a modified Littlewood–Paley decomposition adapted to the Darboux chart system for the contact structure. The transfer operator extends to a bounded operator on B^β, and on this Hilbert space it is quasi-compact.

---

## Setup and definitions

**Contact Anosov flow.** Let d ≥ 1 and r ≥ 3 be integers. Let M be an orientable closed C^r manifold of dimension 2d+1 and α a C^r contact form on M (so ω = α ∧ (dα)^d is a volume form). The flow F^t: M → M preserves α and is Anosov: there exist λ₀ > 0 and C > 0 with

> ‖DF^t_z|_{E^s}‖ ≤ C · 2^{−λ₀t} and ‖DF^{-t}_z|_{E^u}‖ ≤ C · 2^{−λ₀t} for all t ≥ 0, all z ∈ M.

Since F^t preserves α, the subbundles E^s and E^u lie in the null space of α, giving dim E^s = dim E^u = d. The Reeb vector field v (defined by α(v) = 1 and v ∈ null dα) generates the flow: E^c = span{v}.

Let Λ₀ > 0 be a constant (at least dλ₀) such that

> |det(DF^{-t}_z|_{E^u})| ≤ C · 2^{−Λ₀t} for all t ≥ 0, all z ∈ M.

**Transfer operator.** L^t: C^r(M) → C^r(M) is defined by L^t(u)(z) = u ∘ F^t(z).

---

## Main theorem: Theorem 1.1 (quasi-compactness)

**Theorem 1.1.** *For each 0 < β < (r − 1)/2, there exists a Hilbert space B^β, which is contained in W^s(M) for s < −β and contains W^s(M) for s > β, such that the transfer operator L^t for sufficiently large t extends to a bounded operator on B^β and the essential spectral radius of the extension L^t: B^β → B^β is bounded by*

> *max{2^{−Λ₀t/2}, 2^{−βλ₀t}} < 1.*

**Precise hypotheses of Theorem 1.1:**
1. M is an orientable closed C^r manifold, dimension 2d+1, with r ≥ 3.
2. α is a C^r contact form on M.
3. F^t: M → M is a C^r contact Anosov flow preserving α.
4. λ₀ > 0 and Λ₀ > 0 are the hyperbolicity exponents (Λ₀ ≥ dλ₀).
5. 0 < β < (r − 1)/2.
6. t is sufficiently large (the operator L^t: B^β → B^β is bounded for all sufficiently large t, with the essential spectral radius bound holding uniformly for large t).

The contact hypothesis is load-bearing: the estimate 2^{−Λ₀t/2} on the central part (the part of the operator associated with frequencies aligned along the flow direction) uses the non-integrability of the contact form α essentially. For non-contact hyperbolic flows this estimate fails and the corresponding argument breaks down.

---

## Local chart version: Theorem 3.2

The proof reduces Theorem 1.1 to the following local claim:

**Theorem 3.2.** *There exist positive constants λ∗ and Λ∗ such that the operator L(G, g) for any G: V' → V in H(λ∗, Λ∗) and g ∈ C^r(V') extends to a bounded operator L(G, g): B^β_ν → B^β_{ν'} for any 0 < β < (r−1)/2 and ν, ν' ≥ 2β + 2d + 2. Further, for any ε > 0 and 0 < β < (r−1)/2, there exist ν∗ ≥ 2β + 2d + 2, C∗ > 0, and a family of norms ‖·‖_{(λ)} on B^β_{ν∗} (all equivalent to the standard norm), such that if G: V' → V belongs to H(λ, Λ) with λ ≥ λ∗ and Λ ≥ Λ∗ with Λ ≥ dλ, and if g ∈ C^r(V'), then there exists a compact operator K(G, g): B^β_{ν∗} → B^β_{ν∗} such that*

> *‖L(G, g) − K(G, g)‖_{(λ)} ≤ C∗ ‖g‖∗ · 2^{−(1−ε)min{Λ/2, βλ}}.*

Here H(λ, Λ) is the set of C^r diffeomorphisms G: V' → V of open subsets of the unit disk D ⊂ E (Euclidean space of dimension 2d+1) that preserve the standard contact form α₀ and the Reeb vector field v₀, expand the unstable cone by factor ≥ 2λ, and satisfy det(DG_z|_Y) ≥ 2Λ for all (d+1)-dimensional subspaces Y in the stable cone.

---

## Corollary 1.2 (exponential decay of correlations)

**Corollary 1.2.** *For any 0 < α < min{Λ₀, (r−1)λ₀}/2, there exist finitely many complex numbers η_i with −α ≤ ℜ(η_i) < 0 and integers k_i ≥ 0 for 1 ≤ i ≤ ℓ such that, for any ψ and ϕ in C^r(M), we have the asymptotic*

> *(1/ω(M)) ∫ ψ · ϕ∘F^t dω − (1/ω(M))² ∫ ψ dω · ∫ ϕ dω = Σᵢ Σⱼ Cᵢⱼ(ϕ, ψ) · t^{kᵢ} · 2^{η_i t} + O(2^{−αt}) as t → ∞,*

*where C_{ij}(ϕ, ψ) are constants depending on ψ and ϕ bilinearly.*

This implies exponential decay of correlations for all C^r observables.

---

## Summary of hypotheses

| Requirement | Content |
|---|---|
| Manifold | Closed (compact, no boundary), orientable, C^r, dim = 2d+1 |
| Regularity | r ≥ 3 (minimum; larger r gives smaller essential spectral radius bound) |
| Contact structure | C^r contact form α preserved by the flow; non-integrability of α is used in the central-part estimate |
| Flow | C^r Anosov flow preserving α; exponential expansion/contraction on E^s, E^u |
| Hilbert space parameter | 0 < β < (r−1)/2 |
| Time | L^t acts on B^β boundedly for sufficiently large t |

---

## Remarks on generalization beyond smooth manifolds

The construction of the Hilbert spaces B^β uses the following manifold-specific tools:
1. A finite system of Darboux charts adapted to the contact structure (Darboux theorem for contact geometry).
2. A C^∞ partition of unity on the cotangent bundle T*D, used to define the Littlewood–Paley decomposition.
3. Sobolev space theory on the local charts D ⊂ E (Euclidean space).
4. The non-integrability condition on the contact form α, used in the oscillatory integral estimates for the central-part operator (Section 12).

The non-integrability of α is used in a very specific way: the transfer operator acting on frequency components aligned along the flow direction (the "central part") produces oscillatory phase cancellations that are controlled by the symplectic structure of dα on E^s ⊕ E^u. This mechanism is inherently tied to the smooth differential-geometric structure.

Items 1–3 depend on the smooth manifold structure and its standard analytic tools (partition of unity, Fourier transform on Euclidean space, Sobolev embedding). It is not obvious how to replace them in a setting without smooth structure (e.g., symbolic dynamical systems, piecewise smooth maps, or infinite-dimensional state spaces). The Darboux theorem in particular is specific to finite-dimensional contact geometry; no analogue exists in general Banach manifolds without additional structure.

The paper's author notes (Remark 1.3) that the argument for the central part is related to his earlier work on expanding semi-flows as a simplified model, and that the specific use of the contact form is what makes the estimate 2^{−Λ₀t/2} tight (matching the bound suggested by Selberg's zeta function for constant curvature surfaces). Any category other than smooth contact manifolds would require new ideas for the central-part estimate.
