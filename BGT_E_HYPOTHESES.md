# BGT_E — Bingham-Ostaszewski sequential regular variation (Kendall's Theorem K2)

## Phase 0 — verbatim statement

**Source:** BO 2020 "Sequential regular variation: extensions of Kendall's theorem", Theorem K2 = Kendall's Theorem [BGT 1.9.2]. Lines 182-205 of `Bingham_Ostaszewski_Sequential_Regular_Variation_Kendall_clean.txt`.

> **Theorem K2 (Kendall's Theorem [40, Th. 16], cf. [BGT, 1.9.2]).** For {x_n}_{n∈N} multiplicatively admissible and f : R+ → R+ continuous: if, as n → ∞,
> a_n · f(x_n · λ) → g(λ) (λ ∈ I)
> for some interval I ⊆ (0, ∞), positive sequence {a_n}_{n∈N} and continuous function g : I → R+, then f is regularly varying: for each t > 0,
> K(t) := lim_{x→∞} f(tx) / f(x)
> exists, is finite, multiplicative, and both Baire and measurable. So K(t) = t^ρ for some ρ.

Companion (Theorem 1 of BO 2020 generalizes to Baire functions; Corollary states that if f(x) ~ x^ρ · ℓ(x) with ℓ slowly varying, then a_n ~ c · x_n^ρ / ℓ(x_n)).

## Hypothesis types

- h_1: {x_n} multiplicatively admissible (x_n → ∞, x_{n+1}/x_n → 1).
- h_2 (load-bearing): convergence a_n · f(x_n · λ) → g(λ) for λ in a positive-measure / non-meager set.
- h_3: f continuous (or Baire) and positive-real-valued.
- h_4: g positive and continuous (or Baire).

## Phase 1 — hypothesis × ε_k matrix

Sequential RV is sequence-native, which suits ε_k. Natural choice: x_k = k (multiplicatively admissible? x_{n+1}/x_n = (n+1)/n → 1 ✓), f(k) = |ε_k|·2^k = L(k).

Test h_2: pick λ in some interval, check a_n · L(n·λ) → g(λ).

For n·λ to be integer, restrict to rational λ; or interpret L via a continuous interpolation.

| hyp | check | verdict |
|---|---|---|
| h_1 | x_k = k, x_{k+1}/x_k = (k+1)/k → 1 | SATISFIED. |
| h_3 | f(k) = L(k), positive (taking |·|) | SATISFIED. |
| h_2 | take λ = 2, a_n · L(2n) → g(2). For n=2..3, L(4)=0.039 vs L(2)=0.038, so L(2·2)/L(2) ≈ 1.03. For n=3 → 4, L(6)/L(3) = 0.032/0.041 = 0.78. Across the jump n=3 → 4 with 2n=8: L(8)/L(3) = 0.191/0.041 = 4.66 — emphatically not converging. With a_n trying to normalize: need a_n such that a_n · L(2n) settles, but L(2n) doesn't settle. | **FAILED. h_2 fails because L(2n) doesn't have a stable asymptotic across the k=7 jump.** |
| h_4 | g must be continuous on a non-meager set | UNVERIFIABLE since h_2 fails. |

**Phase 1 verdict: NO_FIT.** Sequential RV's load-bearing convergence hypothesis fails because L(k) doesn't have a stable asymptotic at the index k=7-8 (only 8 coefficients total, and the last two are post-jump).

A variant — restrict the sequence to k ≤ 6 (plateau only) — does provide a clean asymptotic L(k) ~ const ≈ 0.04, and Kendall would fire with K(t) = t^0 = 1, identifying L as a slowly-varying-constant. But this is a finite-range observation; Kendall's theorem requires the asymptotic to hold *as n → ∞*, not just on N ≤ 6. Empirically (with only k ≤ 8 available), we cannot determine whether the post-jump regime is itself a new plateau or a continuing transient.

## Phase 2 — conclusion shape

If Kendall fired (with k restricted or with K ≥ 15-20 coefficients), it would deliver: L(k) = k^ρ · ℓ(k) for slowly-varying ℓ. With ρ = 0 (slow variation) we'd recover the candidate-A representation. ρ ≠ 0 would refine the rate (currently empirically (1/2)^n at the L-level, meaning ρ might be 0 since 2^k is multiplicative not power-law).

Conversion to |μ̂_n(ξ)| polynomial-in-A bound: same as candidate A — slowly-varying multiplicative refinement, no polynomial-in-A delivered directly.

**Phase 2 verdict: SHAPE_PARTIAL** — gives sequence-native RV identification, doesn't directly deliver polynomial-in-A.

## Phase 3 — multi-regime check

Sequential RV is single-regime: a_n · f(x_n·λ) → g(λ) is a SINGLE limit. The k=7 jump is incompatible with this single-limit hypothesis.

**Phase 3 verdict: STRUCTURALLY_BLOCKED at current N=8 data.** A multi-regime extension (e.g., "K_1(t) on a sequence of x_n's ≤ some threshold and K_2(t) afterward") is not in the BO 2020 paper — it would be a research extension.

## Disposition: PARTIAL (within plateau) / NO_FIT (across k=7 jump)

Sequential RV is the categorically cleanest BGT-class candidate for ε_k (sequence-native, requires no Laplace transform or continuous interpolation). It would fire cleanly IF the post-jump regime were either (a) a new plateau (would need to see K ≥ 15-20 with the asymptotic settling), or (b) the jump were a finite-N transient and the true asymptotic was the plateau (would need theoretical chain-side argument that the k=7 jump is transient).

Without (a) or (b), Phase 1 fails at h_2.

**This is the PARTIAL candidate.** Realizable additional input: ε_9..ε_K for K ≥ 15. If the post-jump regime stabilizes, then sequential RV fires with a multi-regime extension that's not in the literature (research-level extension). If the post-jump regime continues to escalate, then sequential RV stays NO_FIT.
