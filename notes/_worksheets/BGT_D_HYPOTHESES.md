# BGT_D — Second-order regular variation (2RV)

## Phase 0 — verbatim statement

**Source:** Hawkes 2RV paper Definition 2.2. Lines 191-203 of `arxiv_2311.02655_Second_Order_Regular_Variation_Hawkes.txt`.

> **Definition 2.2 (Second-order regular variation).** For a function F ∈ RV_∞^α with α ∈ R, if there exist a constant ρ ∈ R and an eventually positive or negative function A on R+ such that
> lim_{t→∞} (F(tx)/F(t) − x^α) / A(t) = x^α · ∫_1^x u^{ρ−1} du, x > 0,
> then F is said to be of second-order regular variation at inﬁnity with first-order index α, second-order index ρ and auxiliary function A. The class of all such functions is denoted 2RV_∞^{α,ρ}(A).

Auxiliary class A_∞^ρ = { A ∈ RV_∞^ρ : lim_{t→∞} A(t) = 0 }.

Companion (Hawkes Theorem 3.2, second-order Karamata representation theorem): for α ∈ R, ρ ≤ 0 and A ∈ A_∞^ρ, F ∈ 2RV iff F has a representation of the form (detail omitted in extract; cited as "Theorem 3.2").

## Hypothesis types

- h_1: F is *first-order* RV with some index α (Definition 2.1).
- h_2: a constant ρ ∈ R and an auxiliary function A ∈ A_∞^ρ exist (eventually one-signed, RV with index ρ ≤ 0, A → 0).
- h_3 (load-bearing): the second-order limit condition (F(tx)/F(t) − x^α) / A(t) → x^α · ∫_1^x u^{ρ−1} du.

## Phase 1 — hypothesis × ε_k matrix

For ε_k, try F(k) := |ε_k|. Test first-order RV with index α = −1 (consistent with 1/2^k ≈ k^{?} not quite, but |ε_k| ~ const·2^{-k} suggests *geometric* not power-law decay — first-order RV is power-law).

Actually, RV is multiplicative (f(tx)/f(t)→x^α): F(tk)/F(k) → x^α. For |ε_k| ~ C·2^{-k}, F(tk)/F(k) → 2^{-(tk-k)} = 2^{k(1-t)}, which depends on k — so |ε_k| is NOT first-order RV (it's *exponentially* decaying, hence regularly varying with index −∞ = "rapidly varying"). 

Try F(k) := |ε_k|·2^k = L(k). Then within plateau, L is bounded between 0.03 and 0.04 — close to constant, slowly varying with index 0. Test F(tk)/F(k) → 1 = t^0. ✓ within plateau. But at k=7 the ratio is 4.72, breaks slow variation. Same as candidate A.

| hyp | check | verdict |
|---|---|---|
| h_1 (F = L slowly varying) | within plateau k=2..6 yes, across k=7 jump no | **FAILED at jump.** |
| h_1 (F = |ε_k| with index ρ_1 = "geometric") | first-order RV is power-law, not geometric — |ε_k| is rapidly varying (de Haan), not RV with finite index | **FAILED categorically.** |
| h_2 (A exists, A ∈ A_∞^ρ) | requires |ε_k| to first BE first-order RV — fails | UNVERIFIABLE because h_1 fails. |
| h_3 (second-order limit) | the deviation (L(k+1)/L(k) − 1) is 0.069, −0.037, −0.061, −0.135 within plateau (drifting toward more-negative), then 3.72 across jump | **FAILED at jump.** |

**Phase 1 verdict: NO_FIT.** Both natural choices of F fail: L=|ε_k|·2^k fails slow-variation at the jump (so 2RV reduces to candidates A-style failure), and |ε_k| itself is rapidly varying (not first-order RV).

## Phase 2 — conclusion shape

If Phase 1 had held, 2RV would deliver: L(k) ~ c + A(k)·c'·something. The auxiliary function A(k) is itself slowly varying with index ρ ≤ 0, decays to 0. This is a second-order *prefactor* refinement on the first-order RV envelope — gives the second-order correction to the slowly-varying L.

For our setup: this would identify a sub-rate at which the plateau itself is *drifting* (L_ratio drifting from 1.07 down to 0.86 across k=2..6 IS suggestive of a slow secular drift in L); 2RV would describe this drift.

But: the rate (1/2)^n in S_n − 7/15 = O((1/2)^n) is *first-order*. 2RV refines the prefactor — it does NOT deliver a stronger rate. So even if 2RV applied with index ρ < 0, the conversion to |μ̂_n(ξ)| via R75/R76/R77 would give a refinement on the slowly-varying multiplicative prefactor, not a polynomial-in-A bound.

**Phase 2 verdict: SHAPE_MISMATCH** (refinement, not polynomial-in-A bound).

## Phase 3 — multi-regime check

2RV is single-regime in the auxiliary A(t). The k=7 jump cannot be encoded in a single A(t) that's regularly varying with index ρ ≤ 0 — A would need to make a 4.72× jump at k=7, which is incompatible with first-order regular variation of A itself.

There's a recent literature on "extended 2RV" (Hawkes Theorem 3.3) that allows broader A, but the load-bearing limit condition (F(tx)/F(t) − x^α)/A(t) → explicit limit is still single-regime.

**Phase 3 verdict: STRUCTURALLY_BLOCKED.**

## Disposition: NO_FIT (at h_1 of first-order RV failure)

2RV theorem requires the underlying function to be first-order RV, which |ε_k| is not (it's rapidly varying) and |ε_k|·2^k is not (slow variation fails at k=7 jump).
