# TAUBERIAN_RESCOPE_A_HYPOTHESES (Flajolet-Sedgewick Ch. VI Singularity Analysis)

**Source PDF:** `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/Flajolet_Sedgewick_Analytic_Combinatorics.pdf`
**Extracted text:** `C:/Collatz/tauberian_extract/A.txt`
**Mode:** E — verbatim from PDF, no inheritance.

---

## Definition VI.1 (Δ-domain) — verbatim (PDF page 405, lines 24332-24340)

> **Definition VI.1.** Given two numbers φ, R with R > 1 and 0 < φ < π/2, the open domain Δ(φ, R) is defined as
>   Δ(φ, R) = { z | |z| < R, z ≠ 1, |arg(z − 1)| > φ }.
> A domain is a Δ–domain at 1 if it is a Δ(φ, R) for some R and φ. For a complex number ζ ≠ 0, a Δ–domain at ζ is the image by the mapping z ↦ ζz of a Δ–domain at 1. A function is Δ–analytic if it is analytic in some Δ–domain.

The Δ-domain is essentially the open disc of radius R > 1 with a slit-and-cone removed near z = 1, with cone half-angle (π − φ) (so the cone fits *outside* the unit disc into the punctured neighborhood of 1). The key property: analytic continuation past the unit circle is required, *except for an acute-angle "spike" at the singular point*.

---

## Theorem VI.1 (Standard function scale) — verbatim (PDF page 397, lines 23623-23649)

> **Theorem VI.1** (Standard function scale). Let α be an arbitrary complex number in ℂ \ ℤ≤0. The coefficient of z^n in
>   f(z) = (1 − z)^{−α}
> admits for large n a complete asymptotic expansion in descending powers of n,
>   [z^n] f(z) ∼ n^{α−1}/Γ(α) · ( 1 + Σ_{k=1}^∞ e_k / n^k ),
> where e_k is a polynomial in α of degree 2k. In particular:
>   [z^n] f(z) ∼ n^{α−1}/Γ(α) · ( 1 + α(α−1)/(2n) + α(α−1)(α−2)(3α−1)/(24 n²) + α²(α−1)²(α−2)(α−3)/(48 n³) + O(1/n^4) )    (13)

(Footnote 2: e_k is a polynomial in α that is divisible by α(α − 1) ⋯ (α − k); 1/Γ(α) vanishes identically when α ∈ ℤ≤0.)

---

## Theorem VI.3 (Transfer, Big-Oh / little-oh) — verbatim (PDF page 406, lines 24362-24381)

> **Theorem VI.3** (Transfer, Big-Oh and little-oh). Let α, β be arbitrary real numbers, α, β ∈ ℝ and let f(z) be a function that is Δ–analytic.
> (i) Assume that f(z) satisfies in the intersection of a neighbourhood of 1 with its Δ–domain the condition
>   f(z) = O( (1 − z)^{−α} (log 1/(1−z))^β ).
> Then one has: [z^n] f(z) = O(n^{α−1} (log n)^β).
> (ii) Assume that f(z) satisfies in the intersection of a neighbourhood of 1 with its Δ–domain the condition
>   f(z) = o( (1 − z)^{−α} (log 1/(1−z))^β ).
> Then one has: [z^n] f(z) = o(n^{α−1} (log n)^β).

---

## Theorem VI.4 (Singularity analysis, single singularity) — verbatim (PDF page 409, lines 24633-24650)

> **Theorem VI.4** (Singularity analysis, single singularity). Let f(z) be function analytic at 0 with a singularity at ζ, such that f(z) can be continued to a domain of the form ζ · Δ_0, for a Δ–domain Δ_0, where ζ · Δ_0 is the image of Δ_0 by the mapping z ↦ ζz. Assume that there exist two functions σ, τ, where σ is a (finite) linear combination of functions in S and τ ∈ S, so that
>   f(z) = σ(z/ζ) + O(τ(z/ζ))   as z → ζ in ζ · Δ_0.
> Then, the coefficients of f(z) satisfy the asymptotic estimate
>   f_n = ζ^{−n} σ_n + O(ζ^{−n} τ*_n),
> where σ_n = [z^n] σ(z) has its coefficients determined by Theorems VI.1, VI.2 and τ*_n = n^{a−1} (log n)^b, if τ(z) = (1 − z)^{−a} λ(z)^b.

Here S is the **standard scale** of functions (from PDF page 408, line 24620-24628):
> S = { (1 − z)^{−α} λ(z)^β | α, β ∈ ℂ }, where λ(z) := (1/z) log(1/(1−z)) ≡ (1/z) L(z).

---

## Theorem VI.5 (multiple singularities) — verbatim (PDF page 412+, lines 24989-25069)

> **Theorem VI.5** (Singularity analysis, multiple singularities). Let f(z) be analytic in |z| < ρ and have a finite number of dominant singularities at ζ_j = ρ e^{iθ_j}, j = 1..r. Assume that there exists a Δ–domain Δ_0 such that f(z) is analytic in the indented disc D = { z : |z| < r, z ≠ ζ_j, |arg(z − ζ_j)| > φ for j = 1..r } for some r > ρ and φ ∈ (0, π/2). [Statement continues with the conclusion that f_n decomposes as sum of contributions from each ζ_j, with controlled error O(ρ^{−n} τ*_n).]

---

## Hypotheses extracted (load-bearing list, for the canonical VI.4 application)

| # | Hypothesis | Source |
|---|---|---|
| h_1 | **f(z) analytic at 0**. | VI.4 lead-in |
| h_2 | **f(z) has a singularity at ζ** (with |ζ| equal to the radius of convergence). | VI.4 |
| h_3 | **f(z) admits analytic continuation to ζ · Δ_0** (a Δ-domain centered at ζ): the function must extend past the disc of convergence into a Pac-Man-shape region around ζ. | VI.4, Def VI.1 |
| h_4 | **Singular expansion** of f exists: ∃ σ in standard scale (finite ℂ-linear combination of (1 − z)^{−α} λ(z)^β functions) and τ ∈ S with τ = o(σ) such that f(z) = σ(z/ζ) + O(τ(z/ζ)) as z → ζ within Δ_0. | VI.4 |

For singular expansion to give n^{-α-1} decay (Δ-analyticity + the standard scale (1-z)^{-α} log^β), the exponent α governs the rate of growth/decay of coefficients.

---

## Notational mapping for our use case

For VI.4 to apply to inputs (1)–(4):

- Construct a generating function f(z) = Σ a_n z^n from input (1) (ε_k) or some derived sequence.
- Find/prove the radius of convergence ρ and dominant singularity ζ on |z| = ρ.
- Prove analytic continuation across |z| = ρ into a Δ-domain at ζ.
- Identify the singular expansion in the standard scale.

The c=7/45 connection would be at the level of the *exponent* α in the singular expansion (1−z/ζ)^{−α}: if α gives rise to a specific n^{α−1} growth/decay matching the 7/45 algebraic structure, that's the SELECTED outcome.

---

## End of A HYPOTHESES extraction.
