# Probe R22 — is f stratum-only — **GATE FAILS: f is NOT stratum-only; the §3 stratum-reduction is void (honest negative)**

**Date:** 2026-07-22  Reuses R7/R9/R10/R21. Probe `probes/probe_strata2_R22.py`. Tests the assumption under the
stratum-reduction framework: does the ratio density f_r(u)=3^r·ρ_r(u) depend only on the stratum j=v₃(u−1)−1, so
the theorem `Σ_m 4^{−m}f(4^{−m})=7/30` collapses to `Σ_j W_j F(j)=7/30`? **It does not.**

## R22-A — STRATUM-ONLY GATE (measurement; the assumption under everything): **f is u-dependent within strata**
Within-stratum spread of f_r(u), by stratum j:

| r | j=0 std/mean | j=0 max/mean | j=1 std/mean | j=2 std/mean |
|---|---|---|---|---|
| 4 | 0.2356 | 1.375 | 0.0626 | 0.0235 |
| 5 | 0.2379 | 1.429 | 0.0669 | 0.0336 |
| 6 | 0.2389 | 1.466 | 0.0688 | 0.0359 |

**The within-stratum spread does NOT shrink with r.** In the dominant j=0 stratum (2/3 of the group) std/mean is
**stable at ≈0.237** and max/mean *grows* (1.37 → 1.47). The higher strata's std/mean is small but *increasing*, not
→0. (The last stratum j=r−1 is always exactly 2 conjugate members — the −S_r/2 pair — so std=0 there is trivial.)
**So f is genuinely u-dependent within a stratum; it is not stratum-only, and the reduction to a stratum profile
F(j) is void as an exact statement.**

## R22-D — F(0) CONVERGENCE (measurement): **the kill — three j=0 orbit points → three different limits**
γ_r(τ_m) for m=1,2,4 (all v₃=0, all stratum j=0), r=1…7:

| r | γ(τ₁) | γ(τ₂) | γ(τ₄) | spread |
|---|---|---|---|---|
| 3 | 0.7030 | 0.4951 | 0.8787 | 0.384 |
| 5 | 0.7108 | 0.4944 | 0.8608 | 0.366 |
| 7 | 0.7171 | 0.4760 | 0.8675 | 0.391 |

**The three j=0 orbit points do NOT agree at any r — they head to three distinct limits (≈0.717, ≈0.476, ≈0.868),
with a stable O(1) spread ≈0.39.** Stratum-only is decisively false *even on the orbit*. (Their mean ≈0.687 is not
the full-stratum mean F(0)=0.667 either — the orbit points are a special subset.) Comparison to Wilson's 0.66841:
the *mean* sits ≈+0.018 above it, but the three values straddle it widely — there is no single "F(0) on the orbit."

## R22-B — THE F(j) TABLE (measurement, NO fit): the stratum **means** are clean and stable
Even though f is not stratum-only, the stratum means F_r(j) are stable across r and clean:

| r | F(0) | F(1) | F(2) | F(3) | F(4) | F(5) |
|---|---|---|---|---|---|---|
| 6 | 0.6667 | 1.4286 | 1.9121 | 2.3723 | 2.8359 | 3.3011 |

with **F(0) = 2/3, F(1) = 10/7** to measured precision. Successive differences dF(j): 0.762, 0.484, 0.460, 0.464,
0.465 — after the first, they **converge to 7/15 = 0.4667** (from just below, creeping up). So the mean profile is
asymptotically **linear with slope 7/15**: F(j) ≈ F(0) + (7/15)j, climbing toward the identity spike f(1)=X_r. (Not
"decreasing toward 7/15 from above" as pre-registered — it undershoots to ~0.46 then approaches from below.)

## R22-C — HAAR CHECK (forced): **PASS** — binning validated
Σ_j (2/3)3^{−j}F_r(j) + f(1)/3^r == 1.000000 exactly, r=4,5,6 (f(1)=X_r = 3.07, 3.53, 4.00). The normalization audit
passes, so R22-A/B's binning is correct — the non-stratum-only verdict is physics, not a binning bug.

## R22-E — WEIGHT/HAAR MISMATCH (measurement): the geometric weight is **95% in j=0**
| j | W_j (geometric) | (2/3)3^{−j} (Haar) | ratio |
|---|---|---|---|
| 0 | 0.31746 | 0.66667 | 0.476 |
| 1 | 0.015869 | 0.22222 | 0.0714 |
| 2 | 3.8e−6 | 0.07407 | 5e−5 |
| ≥3 | ≈0 | … | ≈0 |

**W_0 = 0.3175 is 95% of the total geometric weight Σ W_j = 1/3.** So the theorem `Σ_m 4^{−m}f(4^{−m}) = 7/30` is
95% carried by stratum j=0 (m coprime to 3), 5% by j=1 (m=3,6,…), negligible beyond. And *within* j=0 the weight
4^{−m} concentrates on the smallest m: **m=1 carries ≈75%, m=2 ≈19%** of the whole theorem. So the constant is
dominated by two specific orbit values, **f(τ₁)≈0.717 and f(τ₂)≈0.476** — not stratum means. (The Haar constraint,
by contrast, weights j=0 at only 2/3 and reaches deep into high j; the two equations differ through exactly this
W_j/Haar ratio, which collapses fast — so the second equation carries almost no leverage on the low-j values.)

## Status
**R22: the stratum-only assumption FAILS — honest negative, framework void.** **A** — within-stratum std/mean is
stable at ≈0.24 in the dominant j=0 stratum (not →0); **D** — the three j=0 orbit points γ(τ₁),γ(τ₂),γ(τ₄) converge
to three distinct limits (≈0.717/0.476/0.868, spread ≈0.39). **f is genuinely u-dependent within strata; the
reduction `Σ_j W_j F(j)=7/30` to a stratum profile is void.** The diagnostics survive: **C** binning validated
(Haar sum = 1 exact); **B** the stratum *means* F(j) are stable and clean (F(0)=2/3, diffs→7/15, linear profile);
**E** the geometric weight is 95% in j=0 and ≈75%/19% on m=1/m=2.

**Consequence for the crux (owed to the pen):** R22 **confirms the R20 thin-window reading as the correct one and
retires the stratum-mean reduction.** Because f is not stratum-only, the theorem cannot be reduced to a stratum
profile — it must use the **specific orbit values f(τ_m)**, and R22-E shows those are dominated by **f(τ₁) (≈75%)
and f(τ₂) (≈19%)**. So the owed estimate is the specific values f(τ₁), f(τ₂), and the thin tail — the same object as
R20's window (m≲4), now pinned to two dominant orbit points whose limits (≈0.717, ≈0.476) are what actually set
7/30. Same object as R5's qx+1 step. No fitting; exact Haar/weight gates, labeled numeric density; the stratum-only
assumption reported void as instructed.
