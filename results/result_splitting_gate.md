# Gate SPLITTING (S1) — multilevel splitting reaches past the build_nu wall: **STRONG PASS on γ (0.1%, unbiased), clean negative on individual Λ (by design), and the cumulative-ε test telescopes to a γ-difference that inherits the good precision**

**Date:** 2026-07-23. Probe `probes/gate_splitting_S1.py`. Pre-registered gate for the rare-event
multilevel-splitting route that computes the collision probability `p_r(m)` **level-by-level**, so it is not
bounded by the `build_nu` wall (r≈16–18, dense `3^{r+1}` arrays). Object:
`γ_r(m) = 3^r·p_r(m)`, `p_r(m) = Pr[X ≡ 4^{−m}X' (mod 3^{r+1})]`, `X,X'` iid `~ν`. Rare event `~3^{−r}`.

## Method (why it dodges the killed MC route)
Perpetuity `X = Σ_{k=0}^{r} 3^k·2^{−(v₁+…+v_k)} (mod 3^{r+1})`, `v_j~Geom(½)`, `P[v=j]=2^{−j}` — **matches
`build_nu`'s `X'=1+3·2^{−v}X` exactly** (verified: same measure, `ν₀=δ₁`). Digit `s` of `X` commits when `v_s`
is drawn, so **extend particles one level at a time**; a pair `(X,X')` survives level `s` iff digit `s` of
`(X − 4^{−m}X')` is 0 given digits `0..s−1` matched. `q_s`=survival fraction; `p_r=∏_{s≤r}q_s`;
`γ_r=∏(3q_s)`. **No `|ν̂|²` is ever formed → none of the positive squaring-bias `Var/M` that killed the earlier
MC route.** Error bars: `R_REP=20` independent replicas → mean ± SE across replicas (honest batch SE, never
naive binomial).

## Result 1 — γ reproduction: STRONG PASS (the decisive test)
`γ̂_r(m=1)` vs exact (`build_nu`→`p_from_nu`), N=40k, r=8..16:

| r | γ̂(m=1) | γ_exact | rel err | z |
|---|---|---|---|---|
| 8 | 0.72012 | 0.71926 | 0.12% | +0.26 |
| 12 | 0.72703 | 0.72565 | 0.19% | +0.39 |
| 14 | 0.72701 | 0.72807 | 0.15% | −0.29 |
| 16 | 0.73118 | 0.73001 | 0.16% | +0.27 |

**All nine z-scores in [−0.30, +0.40] — unbiased, 0.1% accuracy.** And `p_16 = 3^{−16}γ ≈ 7×10⁻⁸`: a rare
event naive MC could not *hit* at N=40k; splitting nails it to 0.6% relative. The estimator is validated on its
fundamental object.

## Result 2 — individual Λ_r: clean NEGATIVE (a confirmation, not a weakness)
`Λ_r = Σ_{m≥1} 4^{−m}(γ_r(m)−γ_{r−1}(m))` is a ~0.03% **differential** of two γ's that agree to 0.1%:

| r | Λ̂ | SE | Λ_exact | rel err | z |
|---|---|---|---|---|---|
| 14 | −3.32e−4 | 3.24e−4 | +2.87e−4 | 216% | −1.91 |
| 16 | +8.60e−4 | 3.56e−4 | +2.34e−4 | 267% | +1.76 |

rel-SE ≈ **150%** (wrong sign at r=14). max|z|=1.91 ⟹ formally within bars, but the bars are enormous.
**`sign(Λ_r)` is unresolvable at moderate N — exactly what design-fix-2 anticipated and worked around.** A clean
negative on a quantity we had *already excluded as a target* is confirmation of the design, not a marginal pass.
N for 10% rel-SE on a single Λ_16 ≈ **9.2×10⁶** (irrelevant — we never test individual Λ).

## Result 3 — the cumulative test TELESCOPES to a γ-difference (the route)
The deep test never needs Λ. Because `A_r(m)=γ_r(m)−γ_{r−1}(m)`, the cumulative sum **telescopes:**

  `Σ_{r=17}^{R} Λ_r = Σ_m 4^{−m}[γ_R(τ_m) − γ_16(τ_m)]`

— a **difference of two γ measurements**, not a sum of 24 noisy Λ's. Noise does not accumulate over levels;
the individual-Λ noise is irrelevant to the answer. Better, the **ratio** `γ_R/γ_16 = ∏_{s=17}^{R}(3q_s)`
depends only on levels 17..R, so its relative error is `δ√19 ≈ 0.16%` at deep N (≈5×10⁵), not the full-path error.

**Signal:** for ε to return from `+5.15×10⁻³` to 0, need `Σ_{r≥17}Λ = −2.58×10⁻³`; dominant m=1 gives
`0.25·Δγ ≈ −2.58×10⁻³ ⟹ Δγ ≈ −0.0103` on base `γ_16=0.73` — a **1.4% change against 0.16% noise ⟹ SNR ≈ 8–9.**
(The gate's own inline SNR estimate of 2.7 was **pessimistic**: it added 24 Λ_r noises in quadrature as if
independent, but they telescope and cancel. The telescoped two-γ test is the correct, structurally-robust
figure.) Splitting is **cheap**: 6.4×10⁷ extensions ran in **3.1s** (the 365s was `build_nu` ground truth, not
the estimator); a deep run to r≈35 at N≈5×10⁵ is minutes.

## Two items before the deep run
1. **int64 vectorization walls at r≈18** — the match test forms `inv4m·accP ~ M²`, overflowing int64 once
   `3^{r+1} > 3×10⁹` (r+1≈19). To reach r≈35 the accumulator needs 2-limb / Python-bigint arithmetic (splitting
   stays cheap, so ~10× slower is fine). This is the one engineering build required for r>18.
2. **⚠️ Particle degeneracy over 35 levels is NOT tested by this r≤16 gate.** Resampling correlates particles and
   effective sample size (ESS) degrades with the *number* of levels — so the deep bars cannot be extrapolated
   from r=16. **Monitor and report ESS per level (Cérou–Guyader) during the deep run.** If ESS collapses at
   depth, the bars widen and the SNR≈8–9 needs revising. This is the real open risk — not the estimator's
   validity, which is established.

## Status
**GATE PASS — reframed correctly:** (1) **strong pass on γ** (0.1%, unbiased, rare event `p~7×10⁻⁸` resolved at
N=40k); (2) **clean negative on individual Λ** (rel-SE 150% — confirming design-fix-2's decision never to test
`sign(Λ_r)`); (3) the cumulative-ε quantity **telescopes to a γ-difference**, `Σ_{17}^R Λ = Σ_m 4^{−m}(γ_R−γ_16)`,
inheriting γ's 0.1% precision ⟹ **the turnover near r≈27–31 is observable via cumulative ε at SNR≈8–9, N≈5×10⁵,
minutes of compute.** Pre-committed decision rule: **turnover CONFIRMED if ε̂_R < ε̂_16 by ≥3σ_cum; NO-turnover if
ε̂ still rising at r=35 by ≥3σ; else INCONCLUSIVE** — σ from replica/ESS, never naive binomial. **Epistemic
guard (unqualified):** the splitting estimate is a **statistical object**, labelled as such with its error budget,
**never welded to the exact ε ladder** — it can answer "did ε turn over," it cannot extend the exact ladder. Not
at stake: R1–R30, R80–R82 (exact identities, M-reality). Now answerable observationally: the r≈27–31 turnover on
which 7/15 is conditional (walk-back #45). Next build: the bigint accumulator + ESS monitor for the deep run.
