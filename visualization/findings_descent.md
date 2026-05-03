# 3x+1 descent visualization — findings

**Setup:** 50,000 odd starting values uniformly sampled from [3, 10⁶]. q=3 (standard 3x+1). All orbits converged (no timeouts). Generation: 4.5s; render: <10s.

**Outputs:**
- [generate_descent_data.py](file:///C:/Collatz/visualization/generate_descent_data.py), [render_descent.py](file:///C:/Collatz/visualization/render_descent.py)
- [viz_outputs/descent_a.csv](file:///C:/Collatz/visualization/viz_outputs/descent_a.csv) (4.63M trajectory points)
- [viz_outputs/descent_b.csv](file:///C:/Collatz/visualization/viz_outputs/descent_b.csv) (50K orbit summaries)
- [descent_3d.blend](file:///C:/Collatz/visualization/descent_3d.blend)
- [snapshots_descent/](file:///C:/Collatz/visualization/snapshots_descent/) (8 PNGs)

---

## Headline finding (Vis 1, the descent channels)

**a★ class predicts σ linearly at roughly +12 steps per class increment, monotonically across all six classes.**

| a★_idx j | a★ = 3ʲ | n_orbits | σ_mean | σ_median | log(peak/n)_mean |
|---|---|---|---|---|---|
| 1 | 3   | 1547  | 107.9 | 100 | 1.12 |
| 2 | 9   | 7853  | 119.1 | 112 | 1.24 |
| 3 | 27  | 15684 | 131.7 | 125 | 1.50 |
| 4 | 81  | 15554 | 143.6 | 138 | 2.04 |
| 5 | 243 | 7808  | 155.7 | 149 | 2.92 |
| 6 | 729 | 1554  | 169.0 | 162 | 3.96 |

The a★ class is determined entirely by `n mod 64` (specifically, by the deterministic prefix from state (a=64, c=r)). Yet it predicts σ — a long-range trajectory statistic — to within ~12 steps mean shift per class.

The naive "prefix length" account doesn't explain this slope: the prefix has length `6 + j`, so prefix length grows by exactly 1 per a★_idx increment. The other ~11 steps of σ-growth happen after the prefix completes. So **a★ class is encoding more than just the deterministic prefix length** — it's biasing the post-prefix dynamics too, presumably through (a★, c★) determining the post-prefix reduced state.

### Sanity-check verification

Three confirmatory checks rule out the obvious confounds:

1. **log(n_start) is constant across classes** (mean = 12.80 ± 0.02 for all six classes, std ~1.0 within each class). a★_idx is uncorrelated with log(n) in the sample, so the σ trend is *not* a log(n) confound. (This is expected analytically — a★ depends only on n mod 64 — but verified empirically.)

2. **Joint OLS:** `σ ≈ -41.84 + 10.68·log(n) + 12.18·a★_idx`, R² = 0.09. The per-class slope of +12.18 steps survives controlling for log(n). Trend is real.

3. **Within-class slopes on log(n) match theory:** all six classes regress σ ≈ ~10.42·log(n) + intercept_j, with intercept_j stepping up ~+12 per class. Slope of 10.42 matches the Crandall heuristic σ ≈ 10.43·ln(n) (note: the ln-form, not log₂; an earlier draft had this constant wrong). So the per-class effect is an **additive intercept shift** on top of the standard σ-vs-log(n) baseline, not a slope change.

**Caveat:** R² ≈ 0.09 at k=6 means (log_n, a★_idx) together explain only 9% of σ variance. The +12 per-class shift is structurally real and statistically clean, but it sits on top of the much larger irregular σ noise that's the famous unsolved feature of 3x+1. R² climbs to ~0.20 at k=16 (more residue resolution explains more variance) — see chase results below.

The log(peak/n_start) trend is steeper still — orbits in class j=6 overshoot their starting value by ~e^4 ≈ 55x on average vs ~e^1 ≈ 2.7x for j=1. Higher a★ class produces both longer-running AND higher-flying orbits.

## Class-size distribution (combinatorial bonus)

The number of odd residue classes mod 64 with a★ = 3ʲ follows the binomial **C(5, j-1) = (1, 5, 10, 10, 5, 1)**:

| j | classes | fraction |
|---|---|---|
| 1 | 1  | 1/32 |
| 2 | 5  | 5/32 |
| 3 | 10 | 10/32 |
| 4 | 10 | 10/32 |
| 5 | 5  | 5/32 |
| 6 | 1  | 1/32 |

This is exactly the binomial-coefficient pattern for paths in the prefix state machine. Combined with the σ trend above, it means the **median 3x+1 orbit (in terms of starting class) is in the j=3 or j=4 class**, which sit at σ_mean ≈ 132 / 144 — the "typical" σ for the typical class.

## Vis 2 (σ surface): structure visible but a★ doesn't crisply stratify

Vis 2 plots one point per starting n with x = log(n_start), y = σ, z = log(peak/n_start), colored by a★ class.

[v2_01_iso.png](file:///C:/Collatz/visualization/snapshots_descent/v2_01_iso.png), [v2_03_side.png](file:///C:/Collatz/visualization/snapshots_descent/v2_03_side.png), [v2_04_top.png](file:///C:/Collatz/visualization/snapshots_descent/v2_04_top.png) show:

- The cloud is a **thick swept surface**, not a thin sheet. So σ ≠ f(log_n, log(peak/n)) — the three carry partly independent info.
- [v2_04_top.png](file:///C:/Collatz/visualization/snapshots_descent/v2_04_top.png) (the σ-vs-log(n) projection) shows the famous **parallel band structure** of 3x+1 σ-distributions clearly. Well-known in the literature; nice visual confirmation but not new.
- Color-by-a★ in this cloud is **mixed at every point** — no spatial clustering visible to the eye. The headline σ-class trend is visible only in the *means* (above table), not as a spatial separation in the cloud.

So Vis 2 confirms:
- σ has structure beyond the (log_n, log(peak/n)) plane
- a★ shifts the cloud's *mean* per class but doesn't stratify it visually (within-class variance dominates between-class shift)

## What this run *adds* over the prior qx+1 run

The prior qx+1 visualization was a null result on the "geometry across q" hypothesis. This descent run produces a real finding because the right axes were chosen:

- **a★ on z (not color)** lets the channels physically separate, removing the color-on-color confound that obscured the prior signal
- **Restricting to q=3** removes the divergence noise that dominated the prior render
- **Per-orbit summary stats by class** surfaces the σ trend that would have been hard to see purely from the cloud

Visually, the descent ribbons in [v1_01_iso.png](file:///C:/Collatz/visualization/snapshots_descent/v1_01_iso.png) and [v1_03_side.png](file:///C:/Collatz/visualization/snapshots_descent/v1_03_side.png) all look qualitatively similar (each band is a sloped fan from log(n_start) down to log(1)=0). The story is in the *mean shift per class*, not the band shape.

## Chase: does the +12 trend hold at higher k? Where does the heuristic break?

### k-sweep (sweep_k.py at 200K orbits, n_max = 10⁶)

The per-class shift is **invariant in k from k=6 to k=16**:

| k | per-class σ shift | SE | heuristic prediction | gap | R² |
|---|---|---|---|---|---|
| 6  | +12.2044 | 0.110 | +12.4975 | -0.293 | 0.089 |
| 8  | +12.2194 | 0.092 | +12.4702 | -0.251 | 0.112 |
| 10 | +12.1937 | 0.080 | +12.4495 | -0.256 | 0.134 |
| 12 | +12.1425 | 0.071 | +12.4359 | -0.294 | 0.155 |
| 14 | +12.1274 | 0.065 | +12.4214 | -0.294 | 0.177 |
| 16 | +12.1067 | 0.060 | +12.4095 | -0.303 | 0.198 |

(Heuristic prediction = `1 + b_logn · ln(3)` using the sample-fitted within-class slope; this is what the Crandall heuristic σ ≈ K·ln(n) predicts for the per-j additive shift.)

Slight downward drift from 12.20 → 12.11 across k=6 → 16 (about 0.1 step over 10 powers of k, basically flat). Within-class slope steady at 10.39–10.47. R² grows monotonically with k (more residue resolution explains more σ variance).

**The heuristic gap is also invariant at -0.28 ± 0.03 across all k.** The Crandall heuristic systematically overpredicts the per-class shift by ~0.27 at every k. This is a sub-leading correction to σ ≈ K·ln(n) when restricted to fixed-residue-class subsamples — same gap regardless of k. Plausibly comes from the c★ correction term or from heuristic-itself slop. Not chased further; resolving 0.28 within ~50K-orbit-class samples requires orders more orbits.

### Per-class deviations from joint-OLS prediction (chase_deviations.py at 200K orbits)

All per-class deviations at k=6, 8, 10 are **within ~2σ** of joint-OLS prediction. No class systematically deviates at the 2.5σ threshold; the largest residual was -1.91σ (k=6 j=5). Endpoint classes (j=1, j=k) show some asymmetry but well within sampling noise. **No structural deviation on offer.**

### Cumulative-v coloring (chase #2, render_cumv.py)

Vis 1 with trajectory points colored by `cumv_resid = Σv_2(3n+1) - 2·n_odd_steps` (the Cramér race deviation). See [snapshots_cumv/](file:///C:/Collatz/visualization/snapshots_cumv/).

- **Visually informative:** the gradient (white → blue → pink) within each trajectory ribbon shows the typical Cramér race — orbit drifts behind heuristic mid-trajectory, recovers at convergence. Useful pedagogical visualization.
- **NOT structurally new.** The per-orbit final cumv_resid is algebraically determined by (σ, log n) via the identity `cum_v_final = log₂(n) + o_total · log₂(3)` (every orbit reaches 1, so total halves are forced). The "Δresid ≈ -1.7 per a★_idx" pattern is just the per-class σ shift expressed in v-coordinates. Same finding, different axes.

## Recursive a★ decomposition (T2.5 follow-up)

The post-prefix value v₁ = a★·m + c★ has its own residue structure, so a recursive decomposition is possible: σ ≈ (k₁+j₁) + (k₂+j₂) + σ(v₂). Heuristic predicts each level should contribute the same +12 per-class shift independently. Tested at depths 1 through 8 with k=6 each.

**Result at first 2-3 levels (clean recursion):**

| depth | per-level coefs | sum | R² |
|---|---|---|---|
| 1 | [+12.17] | 12.17 | 0.087 |
| 2 | [+11.96, +11.71] | 23.67 | 0.132 |
| 3 | [+11.54, +11.17, +11.32] | 34.03 | 0.164 |

Each j_d contributes near-+12 σ shift independently. Recursion works as heuristically predicted.

**Result at deeper levels (recursion degrades):**

| depth | sum-of-coefs / depth | R² |
|---|---|---|
| 4 | +10.68 | 0.179 |
| 5 | +10.02 | **0.186 (peak)** |
| 6 | +9.24  | 0.185 |
| 7 | +8.42  | 0.177 |
| 8 | +7.80  | 0.172 |

**Per-level shift drifts downward; R² peaks at depth=5 then declines.** Cross-level correlations between j_d values explain the deflation:

| depth | corr(j_d, j_{d-1}) |
|---|---|
| 2 | +0.001 |
| 5 | +0.028 |
| 8 | +0.179 |

By depth 8, the `reduce_to_odd` step has shrunk log(v_d) to ~6 (value ≈ 400), so the orbit's prior dynamics have constrained v_d's bits — the deep residue is no longer "random" relative to the orbit, and OLS sees correlated regressors.

**Recursive form is strictly weaker than flat single-level at equivalent k:**

| recursion depth (k=6 each) | recursive R² | single-level R² at equiv k | gap |
|---|---|---|---|
| 2 (eff k=12) | 0.132 | 0.155 | -0.023 |
| 3 (eff k=18) | 0.164 | 0.207 | -0.043 |
| 5 (eff k=30) | 0.186 | (extrapolated ~0.30) | ~-0.11 |

Mixed-depth tests confirmed: (k₁=8, k₂=8) gives R²=0.152 vs single k=16 R²=0.198 — about 23% loss. **The factored decomposition loses information that the flat one captures.** That information is the cross-level correlation structure that develops as v_d shrinks.

**Verdict on T2.5:** the post-prefix residue *does* carry independent σ predictive power for the first 2-3 levels (heuristic-consistent), but the recursion degrades through three compounding mechanisms:

1. **Log-budget exhaustion.** Each recursive level descends mean log(v_d) by ~1.0 (≈ -0.31 from prefix's net log-change + ~ -0.69 from `reduce_to_odd` stripping factors of 2). Heuristic σ ≈ K·log(v_d) shrinks accordingly.
2. **Cumulative c★ contamination of residues.** At depth ≤ 3 the v_d_odd residue is independent of c★ history (R² < 0.001). By depth 8, R² of residue on cumulative c★_1..c★_{d-1} reaches 0.078, and the mean residue has drifted from uniform 32 down to 24. The residue is no longer "free" — it's partially determined by the orbit's prefix path history.
3. **+1 prefix-length floor.** The +1 from prefix-length growth doesn't shrink with log(v_d), so deep per-class shifts decline ~20% slower than naive log-scaling would predict.

Empirical observed/naive-log-scaled ratio is consistent at 1.20 across depths 5 and 8, matching the floor model. The recursion is bounded by all three effects compounding — not a free-lunch refinement, but an honest depth-bounded structural decomposition.

## ε(σ) conjecture: ruled out by k-invariance

Sharpened conjecture (after the recursion arc): **ε(σ) ≈ per_j_gap × ⟨ℓ⟩** where per_j_gap = -0.26 (the invariant gap between observed +12.20 and heuristic +12.46) and ⟨ℓ⟩ = ⟨k+j⟩ = (3k+1)/2 is the mean prefix length. At k=6 this predicts ε(σ) = -0.26·9.5 = -2.47 vs observed -2.45 (in TA.1 data) — a 0.8% match suggesting the gap might be a per-prefix-step heuristic deficit summed over the prefix length.

**Discriminator: k-invariance.** ⟨ℓ⟩ grows linearly with k (9.5 → 24.5 from k=6 to k=16), but ε(σ) is empirically k-invariant. So the simple form predicts ε(σ) should grow with k; if it doesn't, the k=6 match was coincidence.

**Test result (test_epsilon_conjecture.py at 200K orbits, k=6..16):**

| k | ⟨ℓ⟩ | -0.26·⟨ℓ⟩ predicted ε | observed ε | gap |
|---|---|---|---|---|
| 6  | 9.50  | -2.470 | -2.191 | +0.278 |
| 8  | 12.50 | -3.250 | -2.200 | +1.050 |
| 10 | 15.50 | -4.030 | -2.199 | +1.831 |
| 12 | 18.50 | -4.810 | -2.215 | +2.595 |
| 14 | 21.50 | -5.590 | -2.221 | +3.370 |
| 16 | 24.50 | -6.370 | -2.221 | +4.149 |

**Gap grows linearly with k at +0.39 per k unit; observed ε(σ) is invariant at -2.21 ± 0.02.** k=6 match was coincidence. The conjecture is FALSIFIED in its simple form.

Per-class slope β_k of ε(j) on j is also not k-invariant (ranges from -0.07 to -0.40 across k=6..16, no clear pattern), ruling out the "ε(j) is linear in j with slope -0.26" sub-hypothesis.

**Independent confirmation:** ε(σ) is k-invariant at -2.21 ± 0.02 in this dataset, matching the closed_form_findings.md TA.1 result that ε(σ) is k-invariant. (Absolute level differs from TA.1's -2.45 because this dataset is at N=2^20 vs TA.1's N=2^25..2^32; small-N drift accounts for most of the ~0.20 difference.)

**The two gaps are decoupled.** Both are k-invariant, both are open analytically, but they live in different parts of the Tao bridge equation σ ≈ α_det(r) + K_h·log_n + ε(σ):
- The **per-j gap of -0.26** sits inside the OLS slope on j (observed +12.20 vs heuristic +12.46).
- The **global ε(σ) ≈ -2.21** sits in the regression intercept (vs naively-expected 0).

Both gaps might reflect the same "K_h is approximate" finite-N effect at a deeper level, but they are NOT the simple multiplicative relationship -0.26·⟨ℓ⟩.

**Net:** the cleanest "ε(σ) is prefix-algebraic" candidate is ruled out. ε(σ) stays where it was — a post-prefix descent constant requiring renewal-theoretic input, per closed_form_findings.md Result 2.

## Tail-shape test for K_eff (1M orbits at N=2²⁵, 2²⁷)

Hypothesis: the K_h − K_eff = 0.486 deficit reflects shifting tail-weight in the σ distribution. If high-quantile slopes match K_h while median is below, body is the source. If quantile slopes differ from K_h and integrate to K_eff, the tail-shape IS the mechanism.

**Result: third option — clean fan-out of per-quantile slopes around K_h.**

| quantile | slope per ln(n) at N=2²⁵ | slope at N=2²⁷ | gap from K_h |
|---|---|---|---|
| **mean** | 10.30 | 10.47 | -0.04 (matches K_h ✓) |
| q0.5 | 10.82 | 10.20 | -0.21 |
| q0.75 | 10.82 | 11.28 | +0.62 |
| q0.9 | 12.70 | 13.60 | +2.72 |
| q0.95 | 14.43 | 14.68 | +4.13 |
| q0.99 | 17.46 | 17.72 | +7.16 |

Mean slope hits K_h exactly (Crandall asymptotic for the mean), but the body and tail diverge:
- **Median slope BELOW K_h** (10.20-10.82) — central body undershoots heuristic
- **Tail slope WAY ABOVE K_h** (q99 ≈ 17.5)
- They balance to produce mean ≈ K_h

**GPD tail-shape ξ is roughly k-invariant** (-0.06 to -0.15 across octaves, slope on log_2(n) ≈ +0.009), all slightly negative (bounded heavy tail). **The tail SHAPE doesn't change with log(n); the tail's POSITION shifts faster than the median.** q99 − q50 grows from 163 (octave 20-21) to 196 (octave 26-27) at N=2²⁷.

This isn't tail-weight shifting per se — it's the σ distribution scaling its **spread** non-linearly in log(n) while keeping its shape fixed. The K_eff deficit at first-passage threshold f is plausibly inherited from the σ-distribution body (which preferentially samples first-crossing populations). Tested explicitly via σ-quantile-stratified s_mean(f) regression — see below.

## σ-quantile stratified s_mean(f) test (extension of experiment 39)

Setup: 1M orbits at N=2²⁷, walk Collatz from each odd start, record first-passage step s(n;f) for f ∈ {2⁵, 2¹⁰, 2¹⁵, 2²⁰}. Three regressions.

### R1: s_mean(f) on log(n) at fixed f

| threshold | slope per ln(n) | gap K_h | R² |
|---|---|---|---|
| 2²⁰ | 10.4653 | +0.037 | 1.0000 |
| 2¹⁵ | 10.4899 | +0.062 | 0.9998 |
| 2¹⁰ | 10.4465 | +0.018 | 1.0000 |
| 2⁵  | 10.4819 | +0.054 | 0.9999 |

**Slope hits K_h = 10.43 cleanly at every threshold** (R² ≈ 1.000). Body-slope hypothesis as stated (slope ~10.5) FALSIFIED — slope is K_h, not the median's 10.5.

### R2: s_mean(f) on log(f) at fixed log_2(n) bin

| log_2(n) bin | -slope = K_eff | gap K_h |
|---|---|---|
| [20,21] | 10.95 | +0.53 |
| [21,22] | 10.97 | +0.54 |
| [22,23] | 10.97 | +0.54 |
| [23,24] | 10.95 | +0.52 |
| [24,25] | 10.98 | +0.55 |
| [25,26] | 10.97 | +0.54 |
| [26,27] | 10.95 | +0.52 |

**K_eff = 10.97 in this dataset, invariant in log(n) bin** (range ±0.02). +0.55 ABOVE K_h. Discrepancy with closed_form_findings.md's K_eff = 9.94 is from different threshold range (theirs ∈ [√N/log N, N^(2/3)], mine ∈ {2⁵..2²⁰}) and likely different averaging convention. **The bin-invariance is the key:** K_eff doesn't depend on log(n), only on the σ-distribution's interaction with the threshold range.

### R3: σ-quantile-stratified K_eff — the headline

| σ-quantile band | n_orbits | mean_σ | **K_eff_band** | gap K_h |
|---|---|---|---|---|
| [0.00, 0.25] (low-σ) | 246,046 | 108.76 | **6.15** | -4.28 |
| [0.25, 0.50] | 248,592 | 158.79 | **9.29** | -1.13 |
| [0.50, 0.75] | 250,533 | 204.64 | **11.95** | +1.52 |
| [0.75, 0.95] | 204,448 | 261.30 | **15.26** | +4.84 |
| [0.95, 1.00] (tail) | 50,381 | 350.62 | **19.76** | +9.33 |

**K_eff varies by 3.2× across σ-quantile bands.** Low-σ orbits: K_eff ≈ 6.15. High-σ tail: K_eff ≈ 19.76. The orbit-count-weighted average across bands ≈ 10.89, matches my global K_eff = 10.97 within rounding.

**Three corollaries:**

1. **K_eff is NOT a fundamental constant — it's a σ-quantile-weighted average of band-specific K_eff_band values.** The closed_form_findings.md "K_eff = 9.94" and my "K_eff = 10.97" are both correct *for their respective threshold weightings* — they sample different parts of the K_eff_band distribution.

2. **The +0.486 slope on log(f) in ε(f) is generated by σ-distribution body-vs-tail spread.** As threshold f changes, the σ-quantile bands that dominate the global K_eff average re-weight, shifting K_eff by ~0.486 per log(f) unit.

3. **K_h is hit at σ-quantile ~0.6** (boundary between bands 2 and 3 where K_eff_band crosses 10.43). Below median: K_eff < K_h (body undershoots). Above median: K_eff > K_h (tail overshoots). Body-slope hypothesis is partially right for low-σ band; the new finding is that the variation continues monotonically up through the tail.

**Linear-in-quantile fit:** K_eff_band(quantile q) ≈ K_h + 13.6·(q − 0.6) over the 5 bands tested. Any closed-form derivation of K_eff (or of the +0.486 slope on log(f)) needs to predict this 13.6-per-quantile slope, not just the weighted average.

**This is structurally substantive new finding:** K_eff joins ε(σ) as a "constant" that's actually an average over a distribution. The Tao bridge gap on log(f) decomposes into a quantile-dependent K_eff_band function. Closed-form derivation moves from "find K_eff" to "find K_eff_band(σ-quantile)" — harder, but better-posed.

## σ-std linear scaling test (closed-form K_eff_band candidate)

Hypothesis: K_eff_band(q) = K_h + b·z_q, generated by Gaussian σ at fixed log(n) with std φ(log n) = a + b·log(n) growing linearly. Predicted b ≈ 3.5 from K_eff_band slope ÷ z-quartile shift.

**Test: linear fit std(σ) ≈ 27.09 + 2.275·log(n), R² = 0.967.**

Linear scaling holds well, but empirical b = 2.275, not 3.5 (gap -1.22). Predicted K_eff_band(q) underestimates empirical at the tails:

| q | z(q) | pred (Gaussian) | empirical | gap |
|---|---|---|---|---|
| 0.125 | -1.150 | 7.81 | 6.15 | -1.66 |
| 0.375 | -0.319 | 9.70 | 9.29 | -0.41 |
| 0.625 | +0.319 | 11.15 | 11.95 | +0.80 |
| 0.875 | +1.150 | 13.05 | 15.26 | +2.22 |
| 0.975 | +1.960 | 14.89 | 19.76 | +4.87 |

**Why the gap: σ distribution is non-Gaussian with k-invariant shape.**

| octave | skew | excess kurtosis |
|---|---|---|
| 20.50 | +0.562 | +0.064 |
| 23.50 | +0.661 | +0.533 |
| 26.50 | +0.660 | +0.585 |

Skew (~+0.66) and excess kurtosis (~+0.53) are essentially **invariant across octaves**. The σ distribution has a fixed right-skewed heavy-tail shape that just shifts (mean = K_h·log(n)) and scales (std = 27.09 + 2.275·log(n)) with log(n).

**Verdict: closed-form mechanism PARTIAL.** The σ distribution decomposes as:
```
σ(log n) = K_h · log(n) + const + φ(log n) · X
K_eff_band(q) = K_h + b · X_q
```
where X is the standardized log(n)-invariant shape (right-skewed, kurt ~+0.53) and X_q is its q-th quantile. For Gaussian X, X_q = z_q and the prediction matches body well, underestimates tails.

**Three open analytic targets remain after this run:**
1. **K_h** — derived (Crandall, closed form ✓)
2. **b = 2.275** — empirical std growth rate, mechanism: renewal-theoretic variance accumulation per Syracuse step
3. **X-shape (skew +0.66, kurt +0.53)** — empirical, mechanism: asymmetric v_2 distribution propagating through random walk

(2) and (3) are renewal-theoretic. Closed forms exist in principle (renewal theory has variance and central-moment results for generic random walks); applying to the Collatz v_2 ~ Geom(1/2) structure would require derivation. Sanity-check confirmed: per-quantile σ slopes from the octave table match K_h + 2.275·z within ±1.1 — the linear-in-z structure of σ at fixed log(n) is real, just b < 3.5.

## Findings worth keeping

1. **a★ class at any k ≥ 6 predicts σ with a clean +12.20 ± 0.05 additive shift per class** — invariant in k, monotonic across all classes. This is a regime where the Crandall heuristic can be tested between residue classes (not just within), and it passes to within ~2%.
2. **Class-size distribution at level k is exactly C(k-1, j-1)** (binomial). Verified at k = 6, 8, 10, 12, 14, 16. Combinatorial fact about the prefix state machine; deserves a one-page proof.
3. **A systematic ~0.27 gap below the Crandall heuristic** holds across all k tested — sub-leading correction worth understanding analytically if anyone wants to refine the heuristic.
4. **The 3D rendering didn't add structural insight beyond what the per-class statistics tables show.** The visual was useful as a sanity check (channels are physically separated, deviations live in tables not in geometry) but the actual findings came from OLS on per-class means. The cumv visualization is the most useful pedagogically but is algebraically redundant with σ.

## One-line summary

**a★ class at k=6 (and any k ≥ 6) predicts σ with a clean +12.20 additive shift per class — exactly what the Crandall heuristic predicts to within 2% — invariant in k, with no structural deviation hiding under it. The cleanest empirical fact from the run is the binomial C(k-1, j-1) class-size distribution; the cleanest unresolved residue is the ~0.27 systematic gap below heuristic.**
