# ADELIC_D — Adelic Poisson summation (Riemann-Roch theorem)

**Source:** Binder Tate-thesis notes §5.3, p. 22-23 (C:/tmp/adelic/Binder_Chicago_REU_Tate_Thesis.txt lines 1264-1438).

## Verbatim statement (Theorem 5.14, p. ~24)

> "Theorem 5.14. Let f(x) be continuous and in L¹(𝔸), with f̂ ∈ ℓ¹(k). Then Σ_{ξ∈k} f̂(ξ/a)/|a| = Σ_{ξ∈k} f(aξ) for all ideles a."

(This is Tate's "Riemann-Roch theorem" — adelic Poisson summation relative to the cocompact discrete subgroup k ⊂ 𝔸.)

## Hypotheses isolated

- **h1 (DOMAIN):** 𝔸 = adeles of number field k.
- **h2 (DISCRETE SUBGROUP):** k embedded in 𝔸 diagonally; k is *cocompact* discrete in 𝔸 (k\𝔸 is compact).
- **h3 (FUNCTION):** f : 𝔸 → ℂ continuous, L¹(𝔸), with f̂ ∈ ℓ¹(k) — i.e., the Fourier coefficients summed over k converge absolutely.
- **CONCLUSION:** Σ_{ξ∈k} f̂(ξ/a)/|a| = Σ_{ξ∈k} f(aξ) — Poisson duality between 𝔸 and k via the idele a.

## Hypothesis × input check

| Hyp | (1) μ_n | (3) R75 Plancherel |
|---|---|---|
| h1 (𝔸 = adeles) | NEEDS_LIFT — μ_n lives on (Z/3^n)*, must lift to 𝔸_ℚ | — |
| h2 (k cocompact in 𝔸) | SATISFIED for ℚ in 𝔸_ℚ in abstract | — |
| h3 (f ∈ L¹(𝔸) and f̂ ∈ ℓ¹(k)) | FAILED — μ_n is a *discrete measure on a finite group*; lifted to 𝔸_ℚ it would be a singular measure on a profinite component (ℤ_3*) tensored trivially elsewhere. It is not an L¹ *function* on 𝔸. The hypothesis on f̂ is also unverifiable. | RELATED — R75 establishes a *finite-group Plancherel* identity on Z/3^n (Theorem 75.1): X_k = Σ_ξ \|μ̂_k(ξ)\|² over high-freq. This is the *finite-group Pontryagin Poisson*, not the *adelic Poisson*. The two are related (adelic Poisson reduces to finite-group Poisson on each level n via approximation), but R75's identity is already as strong as the adelic Poisson would give for this specific measure. |

## Disposition for D

**NO_FIT (with note that R75 already captures what adelic Poisson would give).**

Adelic Poisson summation Σ_{ξ∈k} f̂(ξ/a)/|a| = Σ_{ξ∈k} f(aξ) requires f ∈ L¹(𝔸). Syracuse μ_n is a discrete *measure* on a finite group, not an L¹ function on 𝔸. The finite-group Poisson identity (R75 Theorem 75.1) is the correct, already-proved analog.

The adelic Poisson is the *proof tool inside Tate's FE* — it's not a standalone closure path. Since Tate's FE fails to fit (candidate A NO_FIT), the proof tool inside it doesn't separately fire.

**Structural observation:** R75 is a finite-group Plancherel, not an adelic Poisson. The adelic Poisson reduces to finite-group Poisson after passing through the discrete subgroup quotient k\𝔸 — but in our setting (Z/3^n)* is NOT a k\𝔸 quotient; it's a profinite quotient of ℤ_3*. Different.

## Adelic factorization tag

**ADELIC_FACTORIZATION_INHERENT** — adelic Poisson factors by definition. Doesn't help because the hypothesis class is empty.
