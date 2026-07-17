# Garcia–Young 2023 extraction: disposition

**Date:** 2026-05-11. Author: extraction agent. Files: TRANSLATION + MECHANISM + TRANSLATION_ATTEMPT in this directory.

## Disposition

> **EXTRACTION_FAILS_STRUCTURAL**

Wilson's structural parallel is **correct at the level of size** (√q saturation) and **correct at the level of building block** (Postnikov / truncated p-adic log / Gauss-sum closed form). But the **mechanism producing √q in GY's secondary main term `A`** does not dualize to a bilinear-sum bound on our principal-unit coset.

## One-paragraph rationale

GY's `A` (size `≈ √q`, eq 1.8 and 1.11) arises in the proof of Theorems 1.1/1.2 via Poisson summation in `l` applied to the off-diagonal bilinear `Σ ψ(±n+dl)ψ̄(n) W_{M,N}(±n+dl,n)`, with the inner sum evaluated by Lemma 2.15 (a quadratic Gauss-sum closed form for the Postnikov-coset complete sum `S_{q,d}(ψ,k)`). The √q is the **magnitude of one such Gauss sum** times a Plancherel `1/√q` prefactor (net deterministic factor ≈ 1 in q), with the residual size set by the V-weight Fourier transform `Ŵ` rapid-decay (eq 3.11) restricting the bilinear range to a divisor count `σ_0(|b_ψ|)/√|b_ψ|`. **This is not a bilinear cancellation argument**; it is a Gauss-sum evaluation feeding a divisor-restricted sum whose size is then read off directly. The mechanism cannot transfer to our bilinear `Σ_a 1̂(p·a)·F̂(p·a)` because (i) we have no AFE-style smooth amplitude with rapid Fourier decay — `1̂(p·a)` is a Dirichlet kernel with only `1/|a|` decay; (ii) the hypothesis `d² | q | d³` of Lemma 2.15 forces `r ∈ {1, 2}` at `d = p`, and our empirical range is `r = 8..20`; (iii) attempting to re-derive Lemma 2.15 in the multiplicative coordinate `(1+p)^u` (rather than additive `1+pu`) requires non-trivial new computation that would essentially re-do Theorem 78.4-78.6 at family level. We already have the closed-form evaluation of `F̂(p·a)` (verified F̂ theorem, magnitude `p·√q` on support); what we **need** is a bilinear cancellation among the support points, which GY do not prove (their `A` is the leading non-cancellation term, not a saving).

## Where the mechanism diverges (specific citation)

**Equation (3.11) of Garcia–Young 2023:**
> `|Ŵ^±_{n,d}(x)| ≪_C (M/d) · (1 + |x|M/d)^{-C}` for any `C > 0`.

This rapid decay of the V-weight Fourier transform allows GY to truncate the j-summation at `|j| ≪ q^{1+ε}/M` with negligible error, leaving only a divisor-restricted residue. Our object `1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u)` has Fourier transform decay only `O(min(N, 1/|a|))` — polynomial, not rapid. The truncation step of GY's argument fails on `1̂`.

Consequently, GY's main-term extraction at eq (3.13)–(3.16) cannot be replayed on our bilinear, even after notational translation.

## Adversarial check summary

- **A1 magnitude.** Lemma 2.15 magnitude `√q` matches our F̂ magnitude `p·√q` (with the `p` absorbed into the period-extension factor). **No order-of-magnitude conflict.** But: Lemma 2.15 is not a bilinear bound, so "matching magnitude" doesn't yield C.
- **A2 scope.** Lemma 2.15 requires `(q,3)=1, d²|q|d³`. At `d = p`, this forces `r ∈ {1,2}`. **Our empirical range `r = 8..20` is out of scope.** At `d = p^r` (the natural alternative), only the weaker Lemma 2.16 applies, giving magnitude `q/d = p` (no Gauss-sum √q).
- **A3 hidden infrastructure.** GY use AFE (Lemma 2.17), orthogonality of characters mod `d`, and Heath-Brown's hybrid bilinear (Lemma 2.2). The AFE is **essential** to the (1.9) mechanism: it supplies the `1/√(mn)` amplitudes whose partial-summation against `e_{q/d}(a_ψ l)` truncated at the stationary point `l ≈ q/(d|a_ψ|)` produces `√q`. We have no L-function setup and no AFE; the mechanism doesn't run.
- **A4 extraction difficulty.** Even *attempting* to use Lemma 2.15 for our bilinear requires (i) re-deriving it in the multiplicative coordinate (parametrization mismatch — non-trivial), (ii) extending to `r ≥ 3` (Lemma 2.13(2) truncation fails — requires family-level Theorem 78.4-78.6 analogue, which is open), (iii) actually bounding the residual bilinear (the original ask). **No "just translation" is possible.**

## What the structural parallel *does* tell us (constructively)

- GY's Lemma 2.15 is the **same Gauss-sum building block** that gives our verified F̂ magnitude `p·√q` on support. We already have this closed form at q=3 (Theorem 78.4-78.6); GY's Lemma 2.15 is the family-level template for a clean (closed-form, no Hensel) version — which exists exactly when `q | d³`, i.e., low depth.
- The √q in GY's `A` is **NOT** an analogue of our `√N` saturation. It is an analogue of our F̂ magnitude on a single support point. Our `√N` saturation is a **cardinality factor** from the principal-unit coset size `p^{r-1} = N`, multiplied against the `√q` magnitude — combined effect `√N · √q`. The combined effect is shared between GY's `A` and our `|K|·√q`, but the **decomposition into cardinality × magnitude vs. amplitude × phase** is different.
- The closest published paper in `BURGESS_LITERATURE_FINDINGS.md` is correctly identified by Wilson, but the parallel is at the level of **size of the main term**, not at the level of **machinery producing it**.

## Refs

- [GARCIA_YOUNG_TRANSLATION.md](GARCIA_YOUNG_TRANSLATION.md) — Notation correspondence, parametrization mismatch (1+pu vs (1+p)^u), side-by-side of Lemma 2.15 vs F̂ theorem, side-by-side of eq (1.9) heuristic vs our bilinear.
- [GARCIA_YOUNG_MECHANISM.md](GARCIA_YOUNG_MECHANISM.md) — Step-by-step trace of where `A` arises in §3 of GY: Poisson in `l` (§3.3), Lemma 2.15 application (just before eq 3.13), evaluation at eq (3.13), Ŵ-rapid-decay truncation (eq 3.11), end-state at Lemma 3.7. Both √q contributions identified (Gauss-sum magnitude + Plancherel prefactor); the net √q in `A` is from non-bilinear sources.
- [GARCIA_YOUNG_TRANSLATION_ATTEMPT.md](GARCIA_YOUNG_TRANSLATION_ATTEMPT.md) — Direct translation attempted at `r = 2` (where Lemma 2.15 hypotheses fit). Even allowing the parametrization mismatch to be absorbed by hand, the result is a closed-form evaluation of `F̂(p·a)` — which we already have — not a bilinear bound. Three open paths suggested (A: family-level R78.4-78.6; B: dual F-class Milićević; C: principal-unit decoupling).
