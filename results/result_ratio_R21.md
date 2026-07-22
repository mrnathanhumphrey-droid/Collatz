# Probe R21 — the ratio law — **A gate PASS (independent route), the density is singular at identity but the weight avoids it; C misses (informative)**

**Date:** 2026-07-22  Reuses R7/R9/R10; exact rationals where marked. Probe `probes/probe_ratio_R21.py`. The
plainest reading of γ (no new coordinate — R12-B settled none exists): from the defining congruence
(1+3a)(1+3τ)≡(1+3a′) mod 3^{r+1} with u=1+3τ_m=4^{−m}, **X′/X ≡ 4^{−m}**, so
**γ_r(τ_m) = 3^r·ρ_r(4^{−m})**, ρ_r = law of the ratio R=X′/X, and **f(u) := lim 3^r·ρ_r(u)** is the Haar density
on 1+3ℤ₃ of the ratio of two i.i.d. Syracuse values. The theorem: **Σ_{m≥1} 4^{−m}·f(4^{−m}) = 7/30**.

## R21-A — DIRECT RATIO LAW (forced, a real gate): **PASS** — an independent route to γ
Computed ρ_r(u)=Pr[X′·X^{−1}≡u mod 3^{r+1}] **directly from μ_r⊗μ_r by group division** (integer weights; touches
no τ, no C-table, no engine, no collision formula). Then:
- **γ_r(τ_m) = 3^r·ρ_r(4^{−m}) exactly for all 3^r values of m, r=2…6**, and Σ_u ρ_r(u)=1. This is the first route
  to γ that reuses none of the existing machinery — it **tests** the identity rather than restating it. **#41 not
  incurred.**
- **Weld to 7/30:** Σ_{m=1}^{3^r} 4^{−m}·f_r(4^{−m}) == S_{r+1}/2 exactly, r=2…6 (0.2308 → 0.2331, → 7/30). The
  ratio-density reading reproduces the constant.

So the paper-abstract form is certified: **the weighted value of the Haar ratio-density on the ⟨4⟩-orbit is 7/30.**

## R21-B — THE FULL DENSITY (measurement, NO fit): **near-uniform, with a growing spike at the identity**
f_r(u)=3^r·ρ_r(u) over all 3^r elements of 1+3ℤ/3^{r+1} (⟨4⟩ = the whole group, so this is the full density):

| r | min f | max f | mean | ‖f−1‖²/\|G\| | mean\|f−1\|: v₃(u−1)=1 → r |
|---|---|---|---|---|---|
| 4 | 0.405 | 3.069 | 1.000 | 0.294 | 0.333, 0.429, 0.912, 1.604 |
| 5 | 0.401 | 3.534 | 1.000 | 0.296 | 0.333, 0.429, 0.912, 1.372, 2.069 |
| 6 | 0.395 | 4.000 | 1.000 | 0.297 | 0.333, 0.429, 0.912, 1.372, 1.836, 2.534 |

**f is bounded below (≈0.4) and ≈1 in the bulk** (the v₃(u−1)=1 stratum — most of the group — has mean deviation
0.333, flat in r). **But its max grows with r: max f_r = f_r(1) = X_r** (the accumulation 1+Σ_{j≤r}S_j = 3.07, 3.53,
4.00 ≈ (7/15)r → ∞). **The deviation concentrates near u=1** (mean |f−1| climbs monotonically with v₃(u−1), to 2.53
at v₃=6). So the limiting density is **≈1 almost everywhere with a logarithmically-growing spike at the identity
u=1** — f is singular at the trivial ratio (X′=X), smooth away from it. (u=1 is the DC coincidence X′=X, exactly
γ_r(0)=X_r.)

## R21-D — ORBIT vs BULK (measurement, NO fit): **the geometric weight sits on the smooth part, away from the spike**
| r | geometric-weighted mean\|f−1\| | unweighted mean\|f−1\| | ratio |
|---|---|---|---|
| 4 | 0.111 | 0.444 | 0.250 |
| 5 | 0.109 | 0.444 | 0.246 |
| 6 | 0.109 | 0.444 | 0.244 |

**The `4^{−m}` weighting samples an atypically *smooth* part of f (deviation ≈0.11, a quarter of the bulk 0.44).**
Mechanism: the weight concentrates on small m (m≲4), and for m coprime to 3 (m=1,2,4,5) the point u=4^{−m} has
v₃(u−1)=1 — **far from the identity, in the f≈1 bulk.** The near-identity spike is reached only at m near 3^r (weight
4^{−3^r}). **So the growing singularity at u=1 is weight-suppressed out of the theorem** — the sum Σ4^{−m}f(4^{−m})
lives on the well-behaved region of the density. (Consistent with R20-B: the high-m modes are weight-suppressed.)

## R21-C — ARGMAX AT r=7 (tests the derived prediction): **MISS — the peak is near-degenerate**
Pre-registered from the renewal mechanism: argmax_ξ|μ̂_7(ξ)| = 2⁷=128 (or conjugate 2059), dist-to-trivial (2/3)⁷.
**Measured: argmax = 256 = 2⁸** (|μ̂|=0.07587), with **128 only 4th** (|μ̂|=0.07177, 5% below), dist-to-0 = 0.117 ≈
2·(2/3)⁷, **not (2/3)⁷.** The top five (256, 1931, 2059, 128 — all powers of 2 — within 5% of each other) are
**near-degenerate.** So the clean "argmax = 2^r, dist = (2/3)^r" holds only through r=6; at r=7 a different 2-power
edges ahead. **This is the pre-registered miss, and it is informative exactly as expected:** under Prop 1.17's
superpolynomial decay the additive peak is flattening, its location is not a sharply-determined observable, and the
near-trivial region is a cluster of competing frequencies — not a single migrating mode. The additive-peak thread
is closed by citation (Prop 1.17), and R21-C confirms the peak is too flat to carry a sharp location law.

## Status
**R21: the ratio-density reading is certified and the singularity is shown harmless.** **A GATE PASS** — group
division (independent of all prior machinery) reproduces γ_r(τ_m)=3^r·ρ_r(4^{−m}) for all m, r=2…6, and welds to
S_{r+1}/2 → 7/30 (#41 not incurred). **B** — f is ≈1 in the bulk, bounded below (~0.4), with a
logarithmically-growing spike at the identity u=1 (=X_r); the deviation concentrates near u=1. **D** — the
geometric weight sits on the smooth (far-from-identity) part of f (deviation 0.11 vs bulk 0.44), so **the growing
singularity is weight-suppressed out of the theorem.** **C** — pre-registered miss: argmax_7 = 2⁸ not 2⁷, the peak
is near-degenerate (Prop 1.17 superpoly decay), the (2/3)^r law is r≤6 only.

**Consequence for the crux (owed to the pen):** the theorem is now **`Σ_{m≥1} 4^{−m}f(4^{−m}) = 7/30`**, the
geometric-weighted value of the Haar ratio-density — one sentence a number theorist parses cold. What's owed is
unchanged in substance but cleaner in form: (a) **existence of f at the orbit points that carry weight** — the
*smooth* region (m≲4, u far from identity), where f is bounded and ≈1, **not** the identity spike (which is
weight-suppressed); (b) **evaluation of the weighted sum** to 7/30. This is the same thin window as R20 (m≲4,
O(1)), now read as "f is smooth where the weight lives." Same object as R5's qx+1 step. No fitting; exact
group-division/weld gates, labeled numeric density; C reported as a pre-registered miss.
