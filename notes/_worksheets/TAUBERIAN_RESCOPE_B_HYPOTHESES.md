# TAUBERIAN_RESCOPE_B_HYPOTHESES (Chevalier 2507.15394 Theorem 1.16)

**Source PDF:** `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/arxiv_2507.15394_Tauberian_Square_Root_Singularity.pdf`
**Extracted text:** `C:/Collatz/tauberian_extract/B.txt`
**Mode:** E — verbatim from PDF, no inheritance.

---

## Notational conventions (from earlier pages of the PDF)

- `D` = open unit disc `D(0,1) = {z ∈ ℂ : |z| < 1}` (implicit standard convention; explicit at line 322: "Cauchy's formula applied to g which is continuous on D̄ and holomorphic on D").
- `D̄` = closed unit disc.
- `√` = π_{1/2}(·) = principal-branch square root on ℂ \ ]−∞, 0] (Definition 1.1, p.1 of PDF).
- `D(1, 1)^{1/2}` = the image of D(1, 1) = {|ω − 1| < 1} under the principal square root (compact / closure used when h is required to be holomorphic in a neighborhood of D(1,1)^{1/2} ⊂ ℂ).

---

## Theorem 1.16 — verbatim (PDF lines 665-687)

> **Theorem 1.16.** Let g be a continuous function on D̄ (with value in C ∪ {∞}) whose restriction to D is holomorphic. Let
>   Σ_{n∈N} b_n z^n
> be the power series expansion of g in the neighbourhood of 0. Suppose there exists h_p a meromorphic function defined on a neighbourhood of D(1, 1)^{1/2} ⊂ C that only possesses a pole at 0 with multiplicity M ≥ 1, and suppose that h_p satisfies
>   ∀z ∈ D, g(z) = h_p(√(1 − z))
> Then there exists a unique sequence of constants (d_n)_{n≥1} and a constant D such that for any positive integer K ∈ N \ {0}:
>
>   b_n = D / n^{3/2 − M} · ( 1 + d_1/n^1 + ... + d_K/n^K + O(1/n^{K+1}) )   (1.23)
>
> And if M = 1 then D = − Res(h_p, 0) / √π. Similar formulas can be obtained for the other constants, but we won't calculate them here.

**Proof sketch (PDF lines 688-702):** By assumption on h_p there exists a holomorphic function h in a neighbourhood of D(1, 1)^{1/2} ⊂ ℂ and constants D_M, …, D_1 such that
  h_p(w) = D_M / w^M + ⋯ + D_1 / w + h(w).
In particular we get
  ∀z ∈ D, g(z) = D_M / (√(1−z))^M + ⋯ + D_1 / √(1−z) + h(√(1−z)).
Then applying Lemma 1.10 to the quotients D_j / √(1−z)^j for 1 ≤ j ≤ M and using Theorem 1.14 on the function z ↦ h(√(1−z)), we immediately obtain the desired asymptotic expansion.

---

## Comparison theorem 1.14 — verbatim (PDF lines 483-507)

> **Theorem 1.14.** Tauberian theorem: Asymptotic expansion to any order for square root singularities.
> Let g be a continuous function on D̄ whose restriction to D is holomorphic. Let
>   Σ_{n∈N} a_n z^n
> be the power series expansion of g in a neighbourhood of 0. Suppose there exists a function h that is holomorphic in a neighbourhood of D(1, 1)^{1/2} ⊂ C such that
>   ∀z ∈ D, g(z) = h(√(1 − z))
> Then there exist a unique sequence of constants (c_l)_{l≥1} and a constant C such that for any positive integer K ∈ N \ {0}:
>
>   a_n = (1/n^{3/2}) · ( C + c_1/n + ... + c_K/n^K + O(1/n^{K+1}) ).   (1.15)
>
> And C = − h′(0)/(2√π).

**Difference between 1.14 and 1.16:** 1.14 requires h *holomorphic* (no pole) on a neighborhood of D(1,1)^{1/2}; result is n^{-3/2} asymptotic. 1.16 allows h_p to be *meromorphic with a single pole at 0 of multiplicity M ≥ 1*; result is n^{-(3/2 − M)} asymptotic (i.e., a *growing* sequence for M ≥ 2; growth rate n^{M − 3/2}).

---

## Hypotheses extracted (load-bearing list)

| # | Hypothesis | Location |
|---|---|---|
| h_1 | **g continuous on D̄** (closed unit disc, possibly with value in ℂ ∪ {∞}). | line 665 |
| h_2 | **g holomorphic on D** (open unit disc, i.e. restriction g|_D is holomorphic). | line 666 |
| h_3 | **Σ b_n z^n is the power series expansion of g in a neighborhood of 0**. (i.e. g is determined by its sequence of Taylor coefficients (b_n).) | line 668 |
| h_4 | **∃ h_p meromorphic** on a neighborhood of D(1, 1)^{1/2} ⊂ ℂ. | line 670 |
| h_5 | **h_p has a single pole at 0** with multiplicity M ≥ 1. (No poles elsewhere on D(1,1)^{1/2}.) | line 671 |
| h_6 | **The substitution identity:** ∀z ∈ D, g(z) = h_p(√(1 − z)). | lines 672-674 |

**Conclusion's parameter M:** the multiplicity of h_p's pole at 0. M = 1 corresponds to a square-root singularity (gives n^{-1/2} growth, since 3/2 − M = 1/2 ⇒ n^{−1/2}). M = 2 gives n^{+1/2} (growing). General M gives growth rate n^{M − 3/2}.

---

## Notational mapping for our use case

For Theorem 1.16 to apply to inputs (1)–(4), we need to identify:

- a **target sequence** (b_n)_{n≥1} (a numerical sequence that we want to control asymptotically),
- its **generating function** g(z) = Σ b_n z^n analytic on D and continuous up to D̄,
- a **square-root profile** h_p meromorphic near D(1,1)^{1/2} with the prescribed pole at 0 of multiplicity M, such that g(z) = h_p(√(1−z)).

Candidate target sequences (each tested independently in HYPOTHESIS_CHECK):
  T1. b_n = ε_n (signed exact rationals).
  T2. b_n = |ε_n| · 2^n (normalized magnitudes).
  T3. b_n = |μ̂_n(ξ)|² (the polynomial-in-A target itself).
  T4. b_n = |ε_n|² (squared magnitudes, plain).

The "c=7/45 closure" — Tao's Prop 1.17 / our R77-class spectral threshold — requires control over c=7/45 in the sense that some asymptotic decay coefficient matches the algebraic 7/45 = (4+3)/45. The candidate identification will be evaluated for whether the exponent (M − 3/2) implied by inputs matches the c=7/45 target structurally.

---

## End of B HYPOTHESES extraction.
