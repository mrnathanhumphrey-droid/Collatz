# FAURE candidate A — Faure 2009 "Semiclassical origin of the spectral gap for transfer operators of partially expanding maps"

**Source:** arXiv:0903.2747v1 [math.DS] 16 Mar 2009, extracted via pdfminer (pypdf used custom-font encoding that produced glyph names; pdfminer cracked it).
**Cite:** Faure 2009.

## The model (verbatim, with notation cleaned where the PDF dropped characters)

> "Let g: S¹ → S¹ be a C^∞ diffeomorphism (on S¹ := R/Z). g can be written as g: R → R with g(x+1) = g(x) + 1, ∀x ∈ R. Let k ∈ N, k ≥ 2, and let the map E: S¹ → S¹ be defined by E(x) = kg(x) mod 1.
> Let E_min := min_x (dE/dx)(x) = k min_x(dg/dx)(x). We will suppose that the function g is such that **E_min > 1** so that E is a uniform expanding map on S¹.
> The map E is then a k:1 map (i.e. every point y has k previous images x ∈ E⁻¹(y)). Let τ: S¹ → R be a C^∞ function, and define a map f on T² = S¹ × S¹ by:
>
> f: (x, s) ↦ (x' = E(x) = kg(x) mod 1, s' = s + (1/2π)τ(x) mod 1)            (3)
>
> The map f is also a k:1 map. The map f is a very simple example of a compact group extension of the expanding map E. It is also a special example of a **partially hyperbolic map**."

## Transfer operator and Fourier reduction

> "The pull back operator, also called the Koopman operator, or Ruelle transfer operator:
>
> (F̂ψ)(x) = ψ(f(x))            (5)
>
> For a function ψ(x,s) = φ(x) e^{i2πνs} (a Fourier mode in s, ν ∈ Z):
> (F̂ψ)(x,s) = φ(E(x)) e^{iντ(x)} e^{i2πνs}.
>
> Therefore F̂ preserves the decomposition L²(T²) = ⊕_{ν∈Z} H_ν, H_ν := {φ(x) e^{i2πνs} : φ ∈ L²(S¹)}. The operator F̂ restricted to H_ν ≃ L²(S¹) is:
>
> (F̂_ν φ)(x) := φ(E(x)) e^{iντ(x)}            (7)
>
> The parameter ν is a semiclassical parameter; ν → ∞ is the semiclassical limit."

## Theorem 1 (Discrete spectrum of resonances) — verbatim

> "**Theorem 1. Discrete spectrum of resonances.**
> Let m < 0. The operator F̂_ν leaves the Sobolev space H^m(S¹) invariant, and F̂_ν: H^m(S¹) → H^m(S¹) is a bounded operator and can be written
>
> F̂_ν = R̂ + K̂            (9)
>
> where K̂ is a compact operator, and R̂ has a small norm:
>
> ||R̂|| ≤ r_m := 1/E_min^|m| ≤ k/E_min            (10)
>
> Therefore, F̂_ν has an essential spectral radius less than r_m, which means that F̂_ν has discrete (eventually empty) spectrum of generalized eigenvalues λ_i outside the circle of radius r_m. The eigenvalues λ_i are called **Ruelle resonances**."

## Theorem 2 (Spectral gap in the semiclassical limit) — verbatim

> "**Theorem 2. Spectral gap in the semiclassical limit.**
> If the map f is **partially captive** (definition given page 15) (and m small enough), then the spectral radius of the operator F̂_ν: H^m(S¹) → H^m(S¹) does not depend on m and satisfies in the semi-classical limit ν → ∞:
>
> r_s(F̂_ν) ≤ 1/√E_min + o(1)            (11)
>
> which is strictly smaller than 1 from (3)."

## General bound for r_s without the partially captive hypothesis

> "A general bound for r_s(F̂_ν) (with no hypothesis on f) is given by
>
> r_s(F̂_ν) ≤ (1/√E_min) · exp((1/2) lim_{n→∞} log(N(n))/n) + o(1)            (12)
>
> where the function N(n) will be defined in Eq.(30). This bound is similar to the bound given in [Tsu08a, Theorem 1.1] by M. Tsujii."

## Identification of hypothesis TYPES

(i) **THE MAP:**
- f: T² → T² is C^∞
- f is k:1 (k ≥ 2 integer)
- f is a "compact group extension" of the uniformly expanding map E on S¹
- f has explicit form f(x,s) = (kg(x) mod 1, s + τ(x)/(2π) mod 1) — i.e. a SKEW PRODUCT over an expanding base
- The base map E is uniformly expanding: E_min > 1
- The neutral direction is s ∈ S¹ — the fiber direction
- f is a "partially hyperbolic map" with expanding (x) + neutral (s) splitting

(ii) **THE FUNCTION SPACE:**
- Sobolev space H^m(S¹) for m < 0 (Theorem 1)
- The reduction F̂ → ⊕ F̂_ν via Fourier-in-s on the compact neutral direction is LOAD-BEARING — it converts the partially expanding map's transfer operator into a family of weighted transfer operators on the expanding direction only
- The semiclassical parameter is ν ∈ Z (the Fourier mode in s)

(iii) **THE EXPANSION PROPERTY:**
- UNIFORM expansion on the base: E_min > 1, every x has dE/dx ≥ E_min > 1
- NEUTRAL direction has trivial dynamics (s → s + τ(x)/(2π) is an isometric translation in the fiber, no expansion or contraction)
- The expansion is C^∞ and uniform — NOT measure-theoretic
- "Partially captive" hypothesis (Theorem 2): a specific dynamical condition on the trapped set K ⊂ T*S¹ in the cotangent bundle of the base, requiring that "an initial wave packet ϕ₀ represented as a point on K evolves in k wave packets, but in general only one wave packet remains on K and the (k−1) others escape towards infinity" — i.e. the trapped set is "thin" in a quantifiable Hausdorff-dimension sense

(iv) **THE CONCLUSION:**
- Theorem 1: DISCRETE Ruelle resonance spectrum in Sobolev H^m for m < 0
- Theorem 2: Spectral radius ≤ 1/√E_min + o(1) in semiclassical limit ν → ∞ (this IS the spectral gap)
- Equivalently: F̂_ν is QUASI-COMPACT (compact + small remainder)
- The spectral gap location is 1/√E_min — a real number depending only on the base expansion rate
- The o(1) correction goes to 0 as ν → ∞ — this is asymptotic in the SEMICLASSICAL parameter, not in iteration time

## Critical observations for Phase 1

1. The smoothness setting is C^∞. The Tao recursion / Syracuse setup is on the PROFINITE GROUP (Z/3^n)* (or its inverse limit Z_3*). There is no C^∞ structure. This is the central category gap.
2. The expansion structure: E is uniformly expanding (E_min > 1). Syracuse's analog: the "expansion" comes from the Geom(2) random multiplicative factor 2^{-a_i}, NOT a deterministic uniform expansion on a manifold. The (1+3)^u substrate makes 4 = 1+3 a 3-adic unit, so multiplication by 4 acts as a translation on the 3-adic principal-unit group — NOT as expansion.
3. The skew-product structure: Faure's f is a skew-product. Tao's recursion μ̂_{n+1}(ξ) = Σ_v 2^{-v} A_v μ̂_n(ξ · 2^{-v} mod 3^n) is NOT a skew product in the standard sense — it's a renewal-product over an n-fold iid Geom(2) tuple-index (per C1_TAO_RECURSION_FORM.md §1-2).
4. The semiclassical limit: ν → ∞ in Faure corresponds to a Fourier mode in the FIBER. In Syracuse, the analog might be n → ∞ (the level) OR ξ → large (the Fourier variable on Z/3^n). The matching is not 1-to-1: Syracuse's n is the iteration count, which has no direct Faure analog.
5. The PADE multi-spectral picture (complex pair period ≈9.2, asymptotic at z=1.016) and Faure's resonance spectrum: Faure's spectrum is a set of resonances inside the unit disk (since spectral radius < 1). The PADE picture says ε_k's leading singularity at n=13 is at |z| ≈ 1.57 (a power-series radius), with asymptotic at z ≈ 1.016. Converting Faure resonance λ ↔ singularity z: typically z = 1/λ. So |z| = 1.57 ↔ |λ| ≈ 0.64; z = 1.016 ↔ λ ≈ 0.984. These are spectral radii STRICTLY LESS THAN 1 — consistent with Faure's spectral gap.
6. The complex-conjugate-pair structure (period 9.2, θ ≈ 0.68 rad): Faure's resonances are eigenvalues of F̂_ν in C, generically COMPLEX, and "repulse each other like eigenvalues of random complex matrices" (Faure remarks on Figure 2). The PADE complex pair fits this Faure prediction CATEGORICALLY.

## Cross-references to other Faure-school theorems

- The general bound (12) involves N(n) which is the "number of trapped trajectories at time n" on the cotangent dynamics — this connects to the trapped-set Hausdorff-dimension framework that Faure-Tsujii 2013/2021 develops further.
- The "partially captive" hypothesis is verified by Tsujii [Tsu08a, Theorem 1.2] for almost all τ — it's a generic hypothesis, not a structural constraint.
- The semiclassical band-structure prediction of Faure-Tsujii 2013 generalizes Theorem 2: resonances cluster on vertical bands at spectral radii 1/√E_min^k (k = 0, 1, 2, ...), with density obeying Weyl law in each band.

End Phase 0 for candidate A.
