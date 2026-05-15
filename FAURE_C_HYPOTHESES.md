# FAURE candidate C — Faure-Tsujii 2013 "Band structure of the Ruelle spectrum of contact Anosov flows" (CRAS)

**Source:** C. R. Acad. Sci. Paris, Ser. I 351 (2013) 385–391. Extracted via pypdf cleanly.

## Verbatim hypotheses

> "If X is a **contact Anosov vector field** on a **smooth compact manifold** M and V ∈ **C^∞(M)**, it is known that the differential operator A = −X + V has some discrete spectrum called Ruelle–Pollicott resonances in specific Sobolev spaces."

> "**Definition 2.1.** On a smooth Riemannian compact manifold (M, g), a smooth vector field X generating a flow φ_t: M → M, t ∈ R, is **Anosov** if there exists a φ_t-invariant decomposition of the tangent bundle TM = E_u ⊕ E_s ⊕ E_0, where E_0 = RX and C > 0, λ > 0 such that for every t ≥ 0:
> ||Dφ_t/E_s|| ≤ Ce^{−λt}, ||Dφ_{-t}/E_u|| ≤ Ce^{−λt}.            (2.1)"

> "**Definition 2.3.** The Anosov one-form α ∈ C(T*M) is defined by Ker α = E_u ⊕ E_s, α(X) = 1. X is a **contact Anosov vector field** if α is a smooth contact one-form, i.e. (dα)|_{E_u ⊕ E_s} is non-degenerate (symplectic)."

## Theorem 5.1 (Band structure) — paraphrased

For |Im z| → ∞, the Ruelle-Pollicott eigenvalues of A = −X + V are contained in the union of vertical bands B_k = {z ∈ C, Re(z) ∈ [γ_k^−, γ_k^+]}, k ≥ 0, with γ_{k+1}^± < γ_k^±. The values γ_k^± are given by max/min of the time average along trajectories of the damping function D = V − (1/2) div X|_{E_u}. In each isolated band, the density of eigenvalues is given by Weyl law: N(b) ≍ b^d where dim M = 2d+1.

## Identification of hypothesis TYPES

(i) **THE OBJECT:** smooth compact Riemannian manifold M, smooth contact Anosov vector field X, smooth potential V ∈ C^∞(M). Continuous-time flow, not a discrete map.
(ii) **THE FUNCTION SPACE:** specific anisotropic Sobolev spaces H_C of distributions on M (Theorem 3.3).
(iii) **THE EXPANSION PROPERTY:** UNIFORM hyperbolicity (Anosov: exponential expansion on E_u, exponential contraction on E_s, neutral on E_0 = RX). CONTACT condition required (α smooth contact 1-form).
(iv) **THE CONCLUSION:** discrete Ruelle resonance spectrum, organized into VERTICAL BANDS (in C, with Re separating bands), each band obeying a Weyl law.

## Match to Syracuse / PADE picture

The band structure with **vertical-band organization in C** is the CATEGORICAL match for PADE's complex-conjugate-pair structure (period 9.2, θ ≈ 0.68 rad). Faure-Tsujii predicts: resonances clustered along vertical lines in C with imaginary parts oscillating — exactly what produces cos(nθ+φ) modulation in coefficient asymptotics.

HOWEVER:
- Syracuse is DISCRETE TIME (Tao recursion μ̂_n → μ̂_{n+1}), not continuous-time flow
- Syracuse acts on a PROFINITE GROUP (Z_3* / (Z/3^n)*), not a smooth manifold
- There is no smooth Anosov vector field; there is no contact one-form; there is no Riemannian metric
- The CATEGORICAL OBJECT is incompatible: Faure-Tsujii is intrinsically smooth-manifold continuous-time
