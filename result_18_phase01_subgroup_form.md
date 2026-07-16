# Result 18 — PHASE 0 + PHASE 1: ★ `cross(2)` IS a weighted subgroup sum (confirmed exactly). The concentration gate was near-free, as pre-stated — and Phase 1 as scoped **cannot reach the target**.

**Date:** 2026-07-15. **Verdicts: H_P0 ✓ CONFIRMED (cell decomposition exact) / ★ H_P1FORM ✓ CONFIRMED (the reduction is right) / ⚠️ H_CONC — rule fired NOT CONFIRMED (6th mis-specified rule); substance real but near-free, exactly as pre-registered.**

**Headline: the Phase 2/3 literature plan is aimed at the RIGHT object — `cross(2)` reproduces exactly from R13's subgroup condition alone. But two things deflate the plan: Phase 3a does not dodge the character sum, and Phase 1 as scoped is at k=2 where there is no k to be bounded in.**

Probe: `probe_18_phase01_subgroup_form.py`. Log: `result_18_phase01_log.txt`. Runtime: **0.4 s**.

## H_P0 ✓ — Phase 0 done. The general-k cell expression is exact.

Using R14's grading (audited EXACT at k=3 **and k=4** by R17-A2), the `M^k` addresses collapse to cells `(c_1..c_{k−1}, v_k)` with `mass = ∏_{j<k} G^{(j)}_{c_j}·p_{v_k}`, and

> `cross(k) = [Σ_value (Σ_{cells→value} mass)² − Σ_cells mass²] / diag` — **exact, no approximation**

| q,k | cells | cross (cells) | total − within | rel |
|---|---|---|---|---|
| 3,2 | 12 | 0.3976331361 | 0.3976331361 | 1.4e−16 |
| 3,4 | 11,664 | 1.3490563927 | 1.3490563927 | 1.5e−15 |
| 5,3 | 8,000 | 0.0596237751 | 0.0596237751 | 2.9e−15 |
| 7,3 | 9,261 | 0.0680420817 | 0.0680420817 | 1.0e−15 |
| 11,2 | 1,100 | 0.0001137528 | 0.0001137528 | 6.1e−13 |

## ★ H_P1FORM ✓ — the reduction is RIGHT. `cross(2)` is a weighted subgroup sum.

Built from **R13's condition alone** — `H = ⟨2⟩ ⊆ F_q*`, `a := 2^{−A}`, `b := 2^{−v_2}`, collision iff **`a + j·s·b ∈ H`** with `s = (2^d−1)/q mod q` — plus the cell weights `w(h) ∝ 2^{−ind(h)}`:

| q | d | s | cross (subgroup) | cross (cells) | rel |
|---|---|---|---|---|---|
| 5 | 4 | 3 | 0.038875214298 | 0.038875214298 | 5.4e−16 |
| 7 | 3 | 1 | 0.047025340793 | 0.047025340793 | 3.4e−15 |
| 11 | 10 | 5 | 0.000113752827 | 0.000113752827 | 2.4e−13 |
| 13 | 12 | 3 | 0.000186664203 | 0.000186664203 | 8.4e−13 |
| 17 | 8 | 15 | 0.000758521367 | 0.000758521367 | 8.1e−14 |
| 31 | 5 | 1 | 0.000439948615 | 0.000439948615 | 1.2e−13 |

> **`cross` is a MULTIPLICATIVE-SUBGROUP ADDITIVE-SHIFT incidence sum in `F_q`, weighted by `2^{−ind(h)}` — exponential in the discrete log.** The Phase 2 literature target (Konyagin–Shparlinski, *Character Sums with Exponential Functions*) is the right shelf.

## ⚠️ H_CONC — rule fired NOT CONFIRMED. Sixth mis-specified rule.

| q | d | T=5 | T=10 | T=15 | T=20 | T=25 | T=30 | T=40 |
|---|---|---|---|---|---|---|---|---|
| 5 | 4 | 3.57e−1 | 1.72e−2 | 3.48e−4 | 0 | 0 | 0 | 0 |
| 7 | 3 | 3.90e−1 | 9.59e−3 | 7.01e−4 | 1.43e−6 | 0 | 0 | 0 |
| 11 | 10 | **1.00** | **1.00** | 4.77e−3 | 3.80e−3 | 6.50e−6 | 2.80e−6 | 9.74e−10 |
| 13 | 12 | **1.00** | **1.00** | 1.56e−2 | 2.99e−5 | 3.68e−6 | 1.65e−8 | 5.38e−11 |
| 17 | 8 | 1.00 | 7.70e−2 | 1.55e−2 | 9.18e−5 | 7.99e−6 | 2.26e−8 | 2.85e−10 |
| 31 | 5 | 1.00 | 1.14e−1 | 2.27e−3 | 4.35e−4 | 6.41e−7 | 2.17e−7 | 1.66e−9 |

**Why the rule failed — `d`-blind, again.** Truncation captures **nothing** until `T > d`: the smallest cross pair is `(v_2=1, v'_2=1+d)`, hence `err = 1.00` for `T ≤ d` at q=11, 13. And past that the decay is **stepwise in shells of `d`**, not smooth, so "≤0.10 per 5 steps at every step" was the wrong *shape*. *(Prior failures: step<1.5 vs linear growth (R8); ΔR² vs free fits (R9); relative tolerance vs machine-eps noise (R11); d-blind |j|=1 threshold (R13); q-blind flatness threshold (R15).)*

**Substance:** concentration is real and strong (1.0 → 1e−10 as T: 5→40). **But per my own pre-registered caveat** — *"`p_v ~ 2^{−v}` is geometric BY CONSTRUCTION, so a pass proves only the TAIL is geometric, NOT that the HEAD is boundable, which is the actual math. NECESSARY, NOT SUFFICIENT"* — this gate was always going to be nearly free, and it behaved exactly as labelled. **No credit claimed.**

## ★★ Two honest consequences — the plan correcting itself one phase in

**1. Phase 3a does NOT dodge Phase 3b.** The tail is geometric trivially; **the head (`v_2 ≲ 2d`) carries all the mass — and the head is precisely a subgroup-incidence count.** Truncation *localizes* the character sum; it does not avoid it. My "might give the bound with no character theory at all" was too optimistic and is withdrawn.

**2. Phase 1 as scoped CANNOT reach the target.** The target is **boundedness in k**. Phase 1 put **k=2** in subgroup form — and at k=2 there is no k. The k-dependence lives in Phase 0's exact cell expression, which is **not yet in subgroup form**.

> **⇒ Phase 1b required: the general-k subgroup form**, where the level-j contributions `c_j` appear and `c_j ≤ C_q·r_q^j` (`r_q < 1` for q≥5, `r_3 = 1`) is what must be proved. R16 measured `c_j` decaying at ~0.6 at q=5; that decay is the thing to explain.

Catching this at Phase 1 costs a probe; catching it at Phase 3 would have cost the route.

## Where the plan stands

| phase | status |
|---|---|
| 0 — general-k cross expression | **DONE** (exact, gated) |
| 1 — k=2 subgroup form | **DONE** (exact, gated) — but insufficient for the target |
| **1b — general-k subgroup form** | **REQUIRED, next** |
| 2 — literature (Konyagin–Shparlinski et al.) | pending; target confirmed correct |
| 3a truncation / 3b Weil | 3a demoted — it localizes, does not dodge |
| 4 — assembly | pending |

## Not at stake
R10, R11, R13, R14, R15, R16, R17, R5's rate, R6, R7, R12, THEOREM_C_745, Th 78.1–78.3.

_Reporting discipline: both gates were exact identities, not tolerances. H_CONC's limitation ("necessary, not sufficient") was written into the pre-registration BEFORE the run, so its near-free pass earns no credit — and the rule still fired NOT CONFIRMED and is reported as fired, with the d-blindness named rather than the threshold reshaped. The two consequences DEFLATE my own plan: 3a's optimism is withdrawn and Phase 1's scope is admitted insufficient. Author's structural priors this arc: 19-for-29; six mis-specified decision rules, all failing the same way (blind to `d`, to `q`, or to shell structure) — systematic, not incidental._
