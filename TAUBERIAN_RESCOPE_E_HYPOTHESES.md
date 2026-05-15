# TAUBERIAN_RESCOPE_E_HYPOTHESES (Alberts 2508.20814 Theorem 1.1)

**Source PDF:** `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/arxiv_2508.20814_Explicit_Tauberian_Abelian_Number_Fields.pdf`
**Extracted text:** `C:/Collatz/tauberian_extract/E.txt`
**Mode:** E — verbatim from PDF.

---

## Theorem 1.1 — verbatim (PDF lines 53-93)

> **Theorem 1.1.** Let N, N̂ : ℝ_≥0 → ℂ be functions for which N̂ is nondecreasing and |N(X)| ≤ N̂(X). Suppose that each of F = N, N̂ satisfy
>
> (a) L(s, F) converges absolutely on the region Re(s) > σ_a,
>
> (b) L(s, F) has a meromorphic continuation to Re(s) ≥ σ_a − δ with at most finitely many poles in this region, and
>
> (c) For each σ ≥ σ_a − δ, each T ≥ e for which L(s, F)/s does not have a pole on the vertical line [σ + iT, σ + 2iT], and each Z ≥ e/2
>
>   | ∫_T^{2T} L(σ + it, F) Z^{it} dt | ≪_{F, σ} { T^{η̃}   if σ > σ_a − δ ;  T^η (log T)^β   if σ = σ_a − δ }    (1.1)
>
> for some constants η̃ > 0, η > 0, and β ≥ 0.
>
> Then for each X ≥ e
>
>   N(X) = Σ_{Re(z) ≥ σ_a − δ} Res_{s = z} ( L(s, N) · X^s / s ) + O( X^{σ_a − δ/max{η,1}} (log X)^θ ),
>
> where
>   θ = { 0 if η < 1 ;  β + 1 if η = 1 ;  (b − 1)(1 − 1/η) + β/η if η > 1 }
> and b is the order to which s = σ_a is a pole of L(s, N̂).

L(s, F) is defined (PDF lines 34-39):
> L(s, F) = ∫_1^∞ N(x) x^{−s−1} dx · s (Mellin transform of F).

If N(X) = Σ_{n ≤ X} a_n then L(s, N) = Σ a_n n^{-s} (Dirichlet series).

---

## Hypotheses extracted (load-bearing list)

| # | Hypothesis | Source |
|---|---|---|
| h_1 | **N, N̂ : ℝ_≥0 → ℂ** with N̂ nondecreasing. | line 53 |
| h_2 | **|N(X)| ≤ N̂(X)** for all X ≥ 0 (i.e., N̂ is a dominating envelope). | line 53 |
| h_3 | **L(s, F) absolutely convergent** for Re(s) > σ_a, for both F = N and F = N̂. | line 55 |
| h_4 | **L(s, F) meromorphic continuation** to Re(s) ≥ σ_a − δ with finitely many poles. | line 56 |
| h_5 | **Twisted-moment bound (1.1)** on integral of L(σ + it, F) Z^{it} from T to 2T, uniform in σ ∈ [σ_a − δ, σ_a], Z ≥ e/2. | lines 58-66 |
| h_6 | **Existence of constants:** σ_a, δ > 0, η̃ > 0, η > 0, β ≥ 0 satisfying (1.1). | line 68 |
| h_7 | **b** = order of pole of L(s, N̂) at s = σ_a (for the error exponent θ). | line 93 |

Output: an asymptotic for N(X) with explicit error O(X^{σ_a − δ/max(η,1)} (log X)^θ).

---

## Notational mapping for our use case

For Theorem 1.1 to apply to inputs (1)-(4):

- Construct a counting function N: ℝ_≥0 → ℂ from ε_k or |μ̂_n(ξ)|² or related sequence. The natural choice: N(X) = Σ_{k ≤ X} ε_k or Σ_{k ≤ X} |ε_k| · 2^k.
- Build N̂ as a nondecreasing dominating envelope (e.g. N̂(X) = Σ_{k ≤ X} |ε_k| · 2^k since this is nonnegative; nondecreasingness needs the |ε_k| · 2^k pattern to be monotone, which **inputs (1) shows it is NOT** — values are 0.40, 0.038, 0.041, 0.039, 0.037, 0.032, 0.150, 0.191).
- Identify a Dirichlet series L(s, N) = Σ ε_k k^{−s} and prove it has analytic continuation past Re(s) = σ_a.
- Establish twisted-moment bound (h_5) — this is the load-bearing analytic input.

**h_2 mismatch flag (preliminary, before HYPOTHESIS_CHECK):** N(X) = Σ ε_k has |N(X)| dominated by N̂(X) = Σ |ε_k|, but N̂ is nondecreasing iff |ε_k| ≥ 0 which is automatic; **N̂ exists trivially**, so h_1+h_2 are satisfiable.

Load-bearing failures will appear at h_3 (need Dirichlet series convergence, i.e. polynomial bound on |ε_k|) and h_5 (need twisted moment bound). Both rely on first proving an asymptotic property of ε_k itself, which is the closure target — Mode H circularity warning.

---

## End of E HYPOTHESES extraction.
