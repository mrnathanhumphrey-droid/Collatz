# BGT_C — de Haan Π-variation

## Phase 0 — verbatim statement

**Source:** Hawkes 2RV paper Definition 2.4 (= BGT 3.0.x). Lines 229-235 of `arxiv_2311.02655_Second_Order_Regular_Variation_Hawkes.txt`.

> **Definition 2.4 (Π-variation).** A measurable function F : R+ ↦→ R is said to belong to the class Π, if there exists an eventually positive or negative function A on R+ such that
> lim_{t→∞} (F(tx) − F(t)) / A(t) = log x, x > 0.
> The class of all such functions is denoted F ∈ Π_∞(A). We refer to A as an auxiliary function for F.

Companion (Hawkes Prop 2.5): F ∈ 2RV_{0,0}(A) iff F ∈ Π_∞(F·A).

## Hypothesis types

- h_1: F measurable real-valued.
- h_2 (load-bearing): there exists an "eventually positive or negative" auxiliary function A such that the *additive* limit (F(tx) − F(t)) / A(t) → log x.
- h_3: by BGT Theorem B.2.7, A must itself be slowly varying.

## Phase 1 — hypothesis × ε_k matrix

For ε_k, choose F(k) := |ε_k|·2^k or F(k) := S_n. Test the additive limit.

| hyp | check | verdict |
|---|---|---|
| h_1 | F real-valued | SATISFIED. |
| h_2 with F = L | L(k+1) − L(k) values: +0.0026, −0.0015, −0.0024, −0.0050, +0.1186, +0.0404. Pattern: small / decreasing in plateau, then *jumps* to large + at k=6→7. Cannot fit (L(tx)−L(t))/A(t) → log(x) for any monotone A: across the k=7 jump, the additive gap explodes 24× while log(x) for x near 7/6 = 1.17 is 0.15. | **FAILED.** |
| h_2 with F = S_n | S_n is convergent, so S_n(k+1) − S_n(k) = ε_{k+1} → 0. The ratio (ε_{k+1})/A(k) → log(x) for x near (k+1)/k. Choose A(k) := log((k+1)/k) · 1/ε_{k+1}? — then A oscillates wildly with sign changes, violating "eventually positive or negative". | FAILED. |
| h_3 | A would need to be slowly varying | **FAILED** by Phase 1A on Karamata representation (same data). |

**Phase 1 verdict: NO_FIT.** The k=7 additive jump kills the Π-class hypothesis directly.

## Phase 2 — conclusion shape

Π-variation gives functional-equation form F(tx) − F(t) = A(t)·log(x) + o(A(t)). This is *additive* slowly-varying correction. If satisfied, it would deliver an asymptotic identification of L(k) of the form L(k) = α + β·log(k) + o(1) for k → ∞ — useful as a structural finding, but the polynomial-in-A Fourier-decay closure target needs more.

**Phase 2 verdict: SHAPE_PARTIAL** — would give log(k) correction structure, not directly a Fourier-decay bound. Conversion via R76/R77 to |μ̂_n| not immediate.

## Phase 3 — multi-regime check

Π-class is fundamentally single-regime (single auxiliary function A). The k=7 jump in L is incompatible with a single A(t).

**Phase 3 verdict: STRUCTURALLY_BLOCKED.**

## Disposition: NO_FIT

Π-variation fails categorically at the k=7 additive jump.
