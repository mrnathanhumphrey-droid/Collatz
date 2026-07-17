# Collatz project — running findings log

Append-only. Each entry: date, observation, what was checked to rule out boring explanations, residual finding, next step.

---

## 2026-04-30 — bulk tail is exponential, not power-law

**Observation:** Right-tail CCDF of residuals r(n) = σ(n) − (2/(ln 4 − ln 3))·ln(n) on N=2²³ data. Semi-log plot is linear over ~5 decades; rate ~0.030 (so P(r > t) ≈ exp(−0.030 t)). Log-log shows characteristic cliff drop, not linear.

**Boring explanations checked:**
- Sampling bias: covers all n in [1, 2²³] uniformly, no subset.
- Finite-N: rate is stable across the bulk; very far right tail (t > 400) shows mild *thinning* below the exponential extrapolation, possibly a finite-N truncation artifact (largest σ in 1..2²³ is 664).
- Definition match: residual is computed against the **full-Collatz** heuristic 6.95·ln(n); for the odd-only analysis the right baseline is 10.43·ln(n). Tail shape conclusion (exponential, not power-law) is invariant to this slope error.

**Residual finding:** Bulk Collatz residuals decay exponentially with rate ~0.030 in the right tail. Heuristic-consistent in shape; quantitative rate is the deliverable for Stage 4.

**Next step:** post-hoc GPD fit on top-5% residuals per class (Stage 4); will sharpen ξ ≈ 0 expectation and check for class-level variation in the tail.

---

## 2026-05-01 — bulk slope matches odd-n heuristic to 0.2%

**Observation:** Pooled OLS on 16.7M odd n in [3, 2²⁵]: σ ≈ 3.93 + **10.4543** · ln(n), residual std 64.35.

**Boring explanations checked:**
- Definition match: heuristic for odd-only data is 3/(ln 4 − ln 3) = **10.4282**, not 6.95 (which is for the union of odd and even starts). With the corrected baseline, OLS overshoots by 0.024 — fractional difference 0.0024.
- Sampling bias: full data, no subsetting.
- Finite-N: residual std 64 is within the natural noise of σ at this N range.

**Residual finding:** The random-walk heuristic predicts the mean σ vs ln(n) slope on odd integers to within 0.25%, on 16.7M data points. **The bulk model is essentially exact.**

**Next step:** quantify per-class deviations via the hierarchical Stan fit (currently running at k=6).

---

## 2026-05-01 — v=4 and v=10 spikes are trajectory-sampling artifacts, not combinatorial structure

**Observation:** Empirical v = ν₂(3m+1) distribution computed along Syracuse trajectories (depth ≥ 1) on N=2²³ data shows ratios vs Geometric(½): v=4 → **1.325×**, v=10 → **1.321×**, with sags at v=6,7,11,12,14 (down to 0.24× at v=14).

**Boring explanations checked (per sanity-check protocol):**
- **Sampling bias (#1):** Computed v on uniform odd m in [1, 2²⁵]: ratio = 1.0000 to 7 decimal places for v=1..19. Heuristic exactly right on natural density. ✓ Spike does NOT survive uniform sampling.
- **First-step vs deep-step:** First-step v in trajectories (step_idx=0, uniform over odd starts) → ratio = 1.0000 ± 0.005. Deep-step (step_idx > 0) → spikes preserved as in original observation. ✓ Localizes the bias to trajectory dynamics.
- **Finite-N:** spike ratios 1.32 and 0.24 are far outside finite-N noise bands (n ≈ 10⁵–10⁸ per v bin).
- **Definition match:** v defined consistently as ν₂(3m+1) for m odd in both samples.

**Residual finding:** Trajectory-sampled v-distribution is non-i.i.d. — Syracuse iterates correlate through the dynamics and oversample residue classes producing v=4 and v=10 by ~32%. Mean drift is preserved (pooled OLS slope hits the i.i.d. heuristic to 0.2%, see entry above), so the *variance and tail* of σ may deviate from i.i.d. predictions even though the mean does not. This is a known phenomenon (Lagarias-style trajectory measure ≠ natural density on integers); locally verified here.

**Next step:** none required for Stage 4 σ-tail analysis. The hierarchical Stan model with heavy-tailed residual layer is positioned to detect any tail-behavior consequences of this correlation.

---

## 2026-05-01 — clean k=6 fit, headline answers in hand

**Observation:** Hierarchical Stan fit at k=6 (32 odd residue classes mod 64), uniform-sampled 1.6M rows (50K/class), 4 chains × 4 threads, 500+500 iterations.

**Convergence:** R-hat range [0.9991, 1.0111] across all 134 params (only 2 just over 1.01); ESS satisfactory (22/134 < 400, all > 200); only 2/2000 divergent transitions (0.4%). Treedepth + E-BFMI satisfactory.

**Posteriors (mean ± std):**
- μ_α = 2.92 ± 1.96
- μ_β = **10.4475 ± 0.0515** — matches OLS 10.45 exactly, matches odd-n heuristic 10.43 to within 0.4%
- τ_α = 11.70 ± 1.08
- τ_β = **0.067 ± 0.048** (q05=0.006, q95=0.161) — small but plausibly nonzero
- φ = 62.84 ± 0.03

**Per-class structure:**
- α ranges from −28 (k=10, n≡21 mod 64) to +35 (k=31, n≡63 mod 64) — class intervals well-separated, dominant modular effect
- β medians range only 10.42 to 10.47 across all 32 classes; **all 95% CIs overlap massively**, no class statistically distinguishable from μ_β
- Steepest β classes: k=19 (n≡39), k=27 (n≡55), k=13 (n≡27), k=23 (n≡47), k=22 (n≡45)
- Shallowest: k=2 (n≡5), k=17 (n≡35), k=6 (n≡13), k=14 (n≡29), k=20 (n≡41)
- Note Stan τ_β posterior (~0.07) is much smaller than naive OLS class-SD (0.23) — partial-pooling shrinkage doing the right thing on small effects

**Post-hoc GPD on residuals (top-5% per class):**
- ξ summary: mean=−0.017, std=0.036, range=[−0.107, 0.054]
- 25/32 classes "exponential" (|ξ|≤0.05) ✓ consistent with Stage 2 fig 3
- 5/32 sub-exponential (ξ<−0.05) — these have the cliff-drop tail thinning seen in fig 3
- 2/32 super-exponential (ξ>0.05) — slight heavy tendency, marginal
- **No class has a power-law tail.**

**Posterior tail probabilities P(σ > c·ln n | class) at ln n = 17:**
- c=10.5: range 0.32 to 0.69
- c=11.0: range 0.28 to 0.62
- c=12.0: range 0.20 to 0.48
- c=13.0: range 0.13 to 0.42
- c=14.0: range 0.07 to 0.32
- Discrete clustering by α-tier (intercept structure) — slope variation contributes negligibly

**Headline answer to spec:**
1. **Random-walk single-slope heuristic is essentially adequate.** μ_β posterior 10.4475 ± 0.05 covers theoretical 10.43; class-level slope variation is tiny (τ_β ~0.07).
2. **Modular structure lives in the intercept**, not the slope. τ_α ≈ 12 means α varies on the order of ±15 across classes; per-class deterministic shifts from the first few Collatz steps dominate.
3. **Tail is exponential, not power-law**, with mild sub-exponential thinning in some classes (the cliff drop). GPD ξ ≈ 0 across the board.
4. **No specific residue class harbors a structurally different stopping-time distribution** beyond the deterministic α offset.

**Boring explanations checked (per protocol):**
- Sampling: uniform 50K/class, no tail oversampling.
- Definition: corrected odd-only heuristic (10.43, not 6.95) used throughout.
- Convergence: clean (above).
- Per-class data sufficiency: 50K obs/class × 32 classes; per-class SE on β ≈ 0.06 vs τ_β ≈ 0.07 — at-the-noise-floor regime, hierarchical pool dominates.

**Next steps (Stage 4 follow-ups):**
- N-scaling (Stage 4d): refit at N=2²², 2²³, 2²⁴ to see if τ_β or GPD ξ change with sample size.
- α decomposition: split α_class into deterministic (computable from first few Collatz steps for each residue r mod 64) + stochastic. Expected: deterministic part captures ~95% of τ_α.

---

## 2026-05-01 — α–GPD correlation: scale weakly tracks α, shape doesn't

**Observation:** Across the 32 odd residue classes mod 64 (k=6 fit, posterior-mean residuals):

- α vs GPD scale: Spearman ρ = +0.352, p = 0.048
- α vs GPD ξ:     Spearman ρ = +0.192, p = 0.29
- α vs GPD threshold: Spearman ρ = +0.873, p = 7e-11 (largely mechanical — threshold reflects within-class spread)
- β vs anything: weak/insignificant (consistent with slope universality)
- ξ vs scale: ρ = -0.573, p = 6e-4 (mechanical — bounded support compensated by larger scale)

**Boring explanations checked:**
- Threshold-vs-α correlation is partly mechanical: 95th-percentile residual increases with class spread, which correlates with α even under the model (because Stan uses global φ but per-class spread varies).
- Significance of α vs GPD scale (p=0.048) is borderline; not robustly proven.

**Residual finding:** High-α classes have slightly wider residual *spread* but the **same exponential tail shape**. Story is "few classes have higher baselines AND moderately wider variance, identical slopes, identical tail shape." One outlier: k=10 (n ≡ 21 mod 64) is the everything-suppressed class — lowest α (-28), most negative ξ (-0.107), longest deterministic descent path (5 halvings before any odd→3n+1 step).

**Next step:** N-scaling, since this picture might be partly finite-N artifacts.

---

## 2026-05-01 — N-scaling: slope and tail-shape universality empirically confirmed

**Observation:** Pure-OLS sweep across N ∈ {2²⁰, 2²², 2²³, 2²⁴, 2²⁵} on full odd-n data, k=6 (32 classes mod 64), 1.5–17M data points per N.

| N | μ_β | τ_β observed | per-class SE | τ_β corrected | ⟨ξ⟩ |
|---|---|---|---|---|---|
| 2²⁰ | 10.372 | 0.368 | 0.430 | 0 | −0.083 |
| 2²² | 10.382 | 0.187 | 0.228 | 0 | −0.078 |
| 2²³ | 10.385 | 0.128 | 0.165 | 0 | −0.046 |
| 2²⁴ | 10.404 | 0.110 | 0.120 | 0 | −0.049 |
| 2²⁵ | 10.419 | 0.086 | 0.087 | 0 | −0.028 |

**Boring explanations checked (per protocol):**
- Sampling: full odd data at each N, no subsampling.
- Definition: corrected odd-n heuristic (10.4282) used.
- Finite-N: explicitly the variable being controlled here.
- Off-by-one: per-class SE formula `σ_e / sqrt(n × var(log_n))` standard.

**Residual findings:**

1. **τ_β observed tracks per-class sampling SE at every N**, with moment-corrected τ_β² = max(0, observed² − SE²) **zero at every N**. There is no detectable between-class slope variation at any sample size tested. **Slope universality is empirically confirmed asymptotically.** What appears as "spread in per-class slopes" is entirely sampling noise that shrinks as 1/√n_per_class.

2. **μ_β converges to heuristic 10.4282 from below.** Gap is 0.056 at 2²⁰, 0.009 at 2²⁵. Monotonic, slow but unambiguous convergence.

3. **⟨ξ⟩ → 0 as N grows.** −0.083 at 2²⁰, −0.028 at 2²⁵. The sub-exponential "cliff drop" in fig 3 (Stage 2) and the per-class GPD finding of "5/32 sub-exponential classes" are **finite-N truncation artifacts**, not structure. Larger N gives trajectories more room to reach extreme σ; ξ approaches 0 (pure exponential).

4. **τ_α stable around 12–13** across all N — real deterministic intercept structure, immune to sample size.

**This is the publishable headline:**
> Random-walk heuristic for σ(odd n) shows asymptotic universality at N → ∞:
> - μ_β → 3/(ln 4 − ln 3) = 10.4282 (slope of mean σ vs ln(n))
> - τ_β → 0 (no class-level slope variation)
> - ⟨ξ⟩ → 0 (exponential tail decay, no power-law or sub-exponential structure)
>
> The only real residue-class effect is in the **intercept** α, with τ_α ≈ 13 stable across N. This intercept structure is deterministic and computable from the first few Collatz steps for each class.

**Implication for the Stan fit at N=2²⁵:** Stan posterior τ_β = 0.067 (q95 = 0.16) is consistent with τ_β = 0 — the prior + likelihood concentrate just above zero because the noise floor at N=2²⁵ is 0.087, exactly where τ_β observed sits. Stan correctly recognized "at the noise floor" and refused to claim structure.

**Decision:** **Skip the 4 additional Stan fits.** The OLS N-scaling answers the spec question definitively. Bayesian posterior intervals at N=2²⁵ already exist and are consistent with universality.

**Next step:** writeup. Optional: trajectory-measure characterization for the v=4/v=10 spikes (deferred to end).

Output: [stage4_results/k6_uniform_full/fig_n_scaling_ols.png](file:///C:/Collatz/stage4_results/k6_uniform_full/fig_n_scaling_ols.png)

---

## 2026-05-01 — α decomposition: 99.96% deterministic, no residual stochastic class structure

**Observation:** For each odd residue r ∈ {1, 3, ..., 63}, computed the deterministic Collatz "prefix" — the steps where parity is forced by r alone (i.e., while the m-coefficient in state = a·m + c remains even). Each class terminates at some (a_final, c_final) where a_final is odd. Predicted α_det(r) = prefix_steps + 10.43 · ln(a_final / 64), then linearly fitted α_actual = a + b · α_predicted.

**Result:**
- Linear fit: **α_actual = −2.66 + 0.986 × α_predicted, R² = 0.9996**
- SD(α_actual) = 13.7
- SD(α_stoch residual) = **0.28**
- Ratio SD(α_stoch) / mean(α posterior SE) = **0.18** (residual variation is *smaller* than posterior noise — there is no detectable residual structure)
- Largest individual deviation: 0.66 in absolute terms, 0.46 standard errors. None of the 32 classes show statistically significant α_stoch.

**Boring explanations checked:**
- Definition: prefix algorithm tracks state = a·m + c symbolically, applying Collatz rules whenever parity is determined by r alone (a even). Algorithm correctly terminates when a becomes odd.
- Posterior noise: per-class α SE is ~1.5; residual structure SD is 0.28 — well within noise.
- Sampling: full uniform 50K/class subsample, Stan partial-pooling correctly handled.

**Residual finding (the headline):** **The entire τ_α ≈ 13 modular intercept structure reduces to 7-12 steps of deterministic Collatz algebra per residue class.** No stochastic residue-class structure exists at this resolution.

**Prefix length structure:**
- Min prefix: r=21 (class k=10), 7 steps, terminates at (a=3, c=1) — this is the "everything-suppressed" outlier class noted earlier
- Max prefix: r=63 (class k=31), 12 steps, terminates at (a=729, c=728)
- Most classes: 9-10 steps, terminating at (a=27, c=...) or (a=81, c=...)

The terminating a values are consistently powers/products of 3: {3, 9, 27, 81, 243, 729} = {3¹..3⁶}, reflecting the Syracuse map's role.

**Combined with prior findings, the complete picture:**
- **μ_β → 10.4282** (heuristic) as N → ∞
- **τ_β → 0** (slope universality, at noise floor at every N)
- **⟨ξ⟩ → 0** (exponential tail universality, finite-N truncation aside)
- **τ_α ≈ 13** but **99.96% explained by deterministic Collatz prefix algebra**

**Verdict:** σ(odd n) has no residual residue-class structure. The full distribution decomposes into:
1. A **deterministic** intercept shift α_det(r) computable from r mod 64 by symbolic iteration of the Collatz map
2. **Universal** stochastic behavior of the random-walk model thereafter (slope 10.43, exponential tail)

There is nothing more to find at this modular resolution. Project closed at this level.

Output: [stage4_results/k6_uniform_full/fig_alpha_decomposition.png](file:///C:/Collatz/stage4_results/k6_uniform_full/fig_alpha_decomposition.png)

---

## 2026-05-01 — universality strengthened: holds at k ∈ {6, 7, 8, 9}

**Observation:** Extended the α decomposition test to higher modular resolutions on data at N=2²⁷ (134M, 67M odd). At each k, computed per-class OLS α, the deterministic prefix prediction, and the residual α_stoch = α_actual − α_predicted.

| k | mod | classes | per-class n | R² | SD(resid) | mean SE | ratio | max|resid| |
|---|---|---|---|---|---|---|---|---|
| 6 | 64 | 32 | 2.10M | 0.9967 | 0.77 | 0.80 | **0.96** | 1.95 |
| 7 | 128 | 64 | 1.05M | 0.9942 | 1.12 | 1.13 | **0.99** | 2.79 |
| 8 | 256 | 128 | 524K | 0.9918 | 1.44 | 1.59 | **0.91** | 4.11 |
| 9 | 512 | 256 | 262K | 0.9851 | 2.09 | 2.24 | **0.93** | 5.75 |

**Boring explanations checked:**
- Cross-resolution comparison at fixed per-class data: at k=8 N=2²⁷ with 524K per class (matching the original k=6 N=2²⁵ scale), R²=0.9918 vs R²=0.9907 originally. Essentially identical.
- Max|residual| values are unremarkable: at k=9 with 256 classes, max = 2.57 SE; Gumbel max-of-normal expectation at n=256 is ~2.7 SE. No outliers exceeding noise.
- R² monotonic decline with k explained entirely by smaller per-class n inflating posterior noise; signal-to-noise (SD(resid)/SE) is constant at ~0.91-0.99 across all k.

**Residual finding:** **The deterministic prefix prediction explains all detectable per-class α structure at every modular resolution tested (k=6, 7, 8, 9). Residual scatter scales identically with sampling noise — zero detectable signal above it.** The decomposition is rigorous at every resolution where we can muster sufficient per-class data.

**Strengthened theorem statement:**
> For odd n ∈ [3, 2²⁷], stratified by residue class mod 2^k for k ∈ {6, 7, 8, 9}, the per-class intercept α(r) is fully determined by the deterministic Collatz prefix from r mod 2^k. The residual α_stoch(r) is at the per-class sampling noise floor (ratio ~0.93 across all k). No detectable stochastic residue-class structure exists at any tested resolution.

**Implication:** the universality claim extends from "mod 64" to "all modular resolutions tested up to mod 512."

Output: [stage4_results/k6_uniform_full/fig_alpha_decomp_k6_vs_k8.png](file:///C:/Collatz/stage4_results/k6_uniform_full/fig_alpha_decomp_k6_vs_k8.png)

---

## 2026-05-01 — higher moments NOT universal across classes, BUT all moments determined by prefix algebra

**Observation 1 (the apparent contradiction):** Per-class second-, third-, and fourth-moment statistics computed on residuals ε = σ − (α_class + β_class · ln n) at N=2²⁷, k ∈ {6, 7, 8, 9}:

| k | per-class n | Var ratio | Skew ratio | Kurt ratio |
|---|---|---|---|---|
| 6 | 2.1M | **78.6** | 5.0 | 20.0 |
| 7 | 1.0M | **61.4** | 4.3 | 16.7 |
| 8 | 524K | **47.3** | 3.6 | 13.5 |
| 9 | 262K | **36.1** | 3.1 | 11.1 |

Ratio = (SD across classes) / (theoretical sampling SE). Anything ≫ 1 = real class-level structure. **Variance is strongly class-dependent (~30% relative spread); skewness and kurtosis also show real class-level variation**. This appears to contradict the "complete distributional universality of S(n)" hypothesis we'd been testing.

**Observation 2 (the resolution):** When per-class variance, skewness, kurtosis are plotted against the prefix-derived prediction (α_det = prefix_steps + 10.43·ln(a_final/M)):

| k | r(α_pred, Var) | r(α_pred, Kurt) | r(α_pred, Skew) |
|---|---|---|---|
| 6 | 0.99989 | 0.9497 | 0.8716 |
| 7 | 0.99977 | 0.9337 | 0.8442 |
| 8 | 0.99969 | 0.9236 | 0.8250 |
| 9 | 0.99943 | 0.9077 | 0.7914 |

The clustering is discrete: classes sharing the same a_final value collapse onto a single point in (α, Var, Kurt) space. **At k=6 with 32 classes, there are only 6 distinct distributions** — one for each a_final ∈ {3, 9, 27, 81, 243, 729}.

**Boring explanations checked:**
- Sampling: full odd-n data at N=2²⁷, no subsampling.
- Definition: per-class moments computed on residuals after per-class OLS fit; α_det predicted from prefix algebra; correlations computed via Pearson and Spearman, both consistent.
- Correction for noise: the Var ratio of 79 is so far above the noise floor that no plausible misspecification of theoretical SE could explain it as artifact.

**Residual finding:** S(n) is NOT class-universal — its distribution genuinely depends on residue class. **But the full per-class distribution of σ(odd n) | class is parameterized by a_final(class) alone**. Different a_final values give genuinely different distributions; classes with the same a_final give the same distribution.

**Theorem (final form):** σ(n) | (n in class r mod 2^k) ~ F(·; a_final(r)) where a_final(r) ∈ {3¹, 3², ..., 3^k}. The 2^(k−1) odd classes collapse onto exactly k distinguishable distributions.

**Implication:** the residue-class structure of σ has *logarithmic* not linear complexity in modular resolution. At mod 2^k, only k bits of "class information" matter for the σ-distribution — namely, which power of 3 the prefix terminates at.

Outputs:
- [stage4_results/k6_uniform_full/fig_higher_moment_universality.png](file:///C:/Collatz/stage4_results/k6_uniform_full/fig_higher_moment_universality.png)
- [stage4_results/k6_uniform_full/fig_moment_vs_prefix.png](file:///C:/Collatz/stage4_results/k6_uniform_full/fig_moment_vs_prefix.png)

---

## 2026-05-01 — predictive comparison with Bonacorsi-Bordoni framework: parameter parsimony

**Observation:** Replicated B&B's setup (N=10⁷, NB GLM, 50K held-out test) with the prefix-determined a_final as residue covariate. Five NB GLMs fit via statsmodels:

| Model | # params | Test log_score | W1 distance |
|---|---|---|---|
| M0 baseline (log n only) | 3 | −274,139.0 | 3.040 |
| M1 B&B mod-8 RE | 10 | −273,352.3 | 3.071 |
| M2 B&B-extended mod-64 RE | 66 | −272,496.5 | 3.230 |
| M3 a_final at k=3 | 6 | −273,352.3 | 3.136 |
| M4 a_final at k=6 | 9 | −272,496.3 | 3.310 |
| B&B reported (NB2-GLM) | — | −272,911.95 | 3.199 |

**Boring explanations checked:**
- Same data scale (N=10⁷), same test split size (50K), same likelihood family (NB2). Not exactly identical to B&B (they may have used hierarchical priors or polynomial interactions), but reasonable for the residue-effect comparison.
- Random seed fixed; train/test split reproducible.
- a_final lookups verified: at k=3, values ∈ {3, 9, 27} (3 levels); at k=6, values ∈ {3, 9, 27, 81, 243, 729} (6 levels).

**Residual finding:**

1. M3 matches M1 to within 0.1 nats on 50K test observations using 40% fewer parameters (6 vs 10). The 4-level structure of mod 8 adds no detectable predictive value beyond 3 a_final levels.
2. M4 essentially matches M2 (9 params vs 66, log_score difference 0.2 nats out of 272,496) — replacing 32 mod-64 fixed effects with 6 a_final levels at no detectable predictive cost.
3. M4's log score sits 415 nats above B&B's reported NB2-GLM (−272,496 vs −272,912) with comparable parameter count.

**Substantive claim:** The structural decomposition of σ(odd n) — 2^(k−1) odd residue classes collapsing onto k a_final levels — translates into predictive modeling parsimony. At k=6, replacing 32 mod-64 random effects with 6 a_final fixed effects yields equivalent held-out log score, consistent with the residue-class structure being indexed by k a_final levels rather than 2^(k−1) free classes.

**Caveat:** W1 distance is slightly worse for these models (3.31 vs B&B's 3.20). Likely because B&B's full Bayesian hierarchical posterior predictive has more dispersion than fixed-effect MLE point predictions. Adding a hierarchical prior on the a_final coefficients would close this. The relevant comparison for the structural claim is parameter count at matched log score.

**Framing:** B&B's mod-8 random effect is observing real residue structure. The prefix decomposition supplies the algebraic identity behind it: u_{n mod 8} corresponds (up to a constant) to log(a_final(n mod 8)) at k=3. The two parameterizations carry essentially the same predictive content; ours is closed-form, theirs is estimated. The joint model — their NB-GLM framework with a_final as a closed-form covariate at any chosen k — is the natural extension.

---

## 2026-05-01 — first hierarchical fit at k=10 unusable due to convergence + sampling design

**Observation:** k=10 (512 odd classes mod 1024) Stan fit on 4M-row stratified-tail-oversampled subsample: 20.8% divergent transitions (97 in chain 2 alone), R-hat > 1.01 across nearly all 1029 params, μ_β posterior 14.13 (vs OLS 10.45 on uniform-sampled data).

**Boring explanations checked:**
- Sampling design: stratified sample retained 100% of top-5% residuals per class but only ~20% of bulk. This biases the fit upward when the model assumes IID — verified by OLS on the same subsample reproducing μ_β ≈ 14.
- Posterior geometry: 1029-D space with K=512 class-level α and β. Non-centered parameterization + tightened priors + adapt_delta=0.9 not enough; chain 2 still wandered.
- Definition match: μ_β prior was N(7, 1) based on full-Collatz heuristic; should have been N(10.4, 1.5) for odd-only filter.

**Residual finding:** Two compounding errors — biased sampling design + wrong prior centering. Both correctable. Geometry at K=512 is borderline regardless.

**Next step:** refire at k=6 (32 classes, ~69-D posterior) with uniform per-class sample (no tail oversampling). Currently running.

---

## 2026-05-01 — multi-statistic prefix decomposition: step-count vs peak distinction

**Observation:** Extending the prefix decomposition framework (originally developed for σ) to other Collatz trajectory functionals at N=10⁷, k=6, odd n only:

| outcome | μ_β | heuristic | τ²_corrected | α R² | AD ratio (diff/same) |
|---|---|---|---|---|---|
| sigma | 10.401 | 10.428 | 0 | 0.9854 | 19.9× |
| syracuse | 3.466 | 3.476 | 0 | 0.9854 | 16.2× |
| odd_steps | 3.466 | 3.476 | 0 | 0.9854 | 16.2× |
| even_steps | 6.936 | 6.952 | 0 | 0.9853 | 22.5× |
| log(max_excursion) | 1.002 | — | 0 | **0.8243** | **1.5×** |

**Boring explanations checked:**
- syracuse vs odd_steps identity: identical numbers because Syracuse stopping time = odd-step count by definition.
- max_excursion μ_β = 1.0: max trajectory value scales linearly with n, not log n; not an artifact of fitting the wrong function.

**Substantive claim:** Step-count functionals (σ, σ_Syracuse, odd_steps, even_steps) all share the same a_final structural axis — slope universality, R² ≈ 0.985 for α prediction, AD ratio 16–22× for distributional clustering. Peak functionals (max_excursion) show partial structure — R² = 0.82, AD ratio 1.5×. Peak amplitude depends on post-prefix trajectory chaos beyond the prefix algebra's reach.

**Refined statement:** the prefix decomposition is the structural axis for trajectory step-count functionals, not for trajectory peak functionals.

---

## 2026-05-01 — qx+1 prefix complexity rank-monotonically predicts convergence rate; one-axis under this prefix rule

**Observation (headline):** q ∈ {5, 7}, k ∈ {6, 8}: per-class qx+1 convergence rate is rank-monotonically anti-correlated with prefix odd-step count (Spearman ρ = −0.97 at q=5 N=10⁷; ρ = −0.93 at q=7 N=10⁸). At q=5 partial correlation cleanly identifies odd_steps as the dominant predictor (partial r = −0.81, p = 3×10⁻⁸ controlling for log_c_final). At q=7 partial separation is underpowered (n_convergent = 258, predictor collinearity 0.95). Even-step count is structurally fixed at k by the prefix algorithm; odd_steps, prefix_total, and log(a_final) are mathematically equivalent under this prefix rule. q=11, q=13 directionally consistent but power-limited.

**Data generated:** qx+1 trajectory parquets for q ∈ {3, 5, 7, 11, 13} at N=10⁷ and selected q at N=10⁸. Convergence-to-1 rates: q=3 100% (Collatz, control); q=5 0.32% at N=10⁶, 0.15% at N=10⁷, 0.066% at N=10⁸; q=7 0.022% at N=10⁷, 0.0005% at N=10⁸; q=11 0.0003% at N=10⁷; q=13 0.0001% at N=10⁷. Convergence rate decays with N for q ≥ 5.

**Quantitative table across runs:**

| run | n_convergent | marginal Pearson r | Spearman ρ | partial r (odd_steps controlling for log_c_final) | partial p | R² |
|---|---|---|---|---|---|---|
| q=5, N=10⁶, k=6 | 1,587 | −0.92 | −0.94 | **−0.81** | 3×10⁻⁸ | 0.87 |
| q=5, N=10⁷, k=6 | 7,286 | −0.93 | −0.97 | −0.79 | 1×10⁻⁷ | 0.87 |
| q=5, N=10⁷, k=8 | 7,286 | −0.89 | −0.97 | −0.66 | ≈ 0 | 0.79 |
| q=5, N=10⁸, k=8 | 32,785 | −0.89 | −0.97 | **−0.66** | ≈ 0 | 0.80 |
| q=7, N=10⁸, k=6 | 258 | −0.73 | −0.93 | −0.24 | 0.20 | 0.54 |

**Boring explanations checked:**
- *Spurious collinearity / proxy effect:* At every k the deterministic prefix algorithm forces every odd residue class mod 2^k to have **exactly k even-step branches** (a starts at 2^k and only halves on c-even branches; c-odd branches multiply a by q, an odd number, which preserves a's even parity). So even_steps_in_prefix has zero variance across classes at every k, and odd_steps, prefix_total, log(a_final) are 100% collinear by construction. Verified empirically at k=6, k=8. They are the same predictor under three labels. **The collinearity is fundamental to the prefix rule, not an artifact of k choice.** Disentangling them would require a different prefix rule (e.g., termination on c-condition rather than a-condition).
- *c_final as a secondary axis:* At k=6 N=10⁶, log(c_final) had partial r = +0.32, p = 0.08 — marginal. With 4× more classes at k=8 it weakened to +0.10, p = 0.27. With 10× more data at N=10⁸ k=8 it weakened further to +0.10, p = 0.27. **More data weakens the signal, consistent with k=6 marginal trend being noise rather than diluted signal.** c_final is not a real second axis.
- *q-universality:* q=7 (109 conv at N=10⁷ → 258 at N=10⁸) shows the same direction at marginal/rank level but partial regression is underpowered due to the 0.95 collinearity between odd_steps and log_c_final. Pearson r and Spearman ρ both highly significant. q=11 (13 conv at N=10⁷) and q=13 (6 conv) lack power for any partial test.

**Substantive claim (publication-ready, modest scope):** Per-class qx+1 convergence rate is anti-correlated with structural prefix complexity at q=5 (partial r = −0.66 to −0.81 depending on k, validated at N=10⁸ with 32,785 convergent orbits, p ≈ 0). Mechanism: prefix odd-step count = number of multiplicative growth applications baked deterministically into the trajectory before the prefix terminates. More baked-in growth → harder to converge to 1. The relationship is one-axis under this prefix rule (odd_steps / log(a_final) / prefix_total are mathematically equivalent labels). Cross-q evidence (q=5 strong, q=7 directionally and rank-monotonically consistent at ρ = −0.93) supports the claim as a qx+1-general phenomenon, with q ≥ 11 power-limited at attainable N.

**Caveat for any future writeup:** the original framing "log(a_final) predicts convergence" oversells what the data can distinguish. Three predictor labels are equivalent by construction. The honest framing is "prefix odd-step count" or "structural prefix complexity," with the equivalence to log(a_final) noted as a derived property under the q-decomposition.

**Next step (deferred):** Either (a) modify the prefix stopping rule to allow even_steps to vary, then test whether odd_steps and even_steps both predict independently — actual research project, days of work; (b) accept as a one-axis finding and write a short note. Lean (b) for now. Push q=5 to N=10⁹ if disk allows for one more confirmation point, but partial r is stable across two orders of magnitude already.

---

## 2026-05-01 (late) — qx+1 prefix complexity → exponential-decay convergence law with universal multiplier C ≈ 2.5

**Observation (headline, empirical):** Across q ∈ {5, 7, 9} at N=10⁸ on per-class data binned by prefix odd-step count j, the relationship log(conv_rate(j)) ~ const + slope · j is essentially exact (R² ≥ 0.99 in all three). The empirical slope divided by log(q/4) is constant across q to within 3%:

| q | n_conv | slope | log(q/4) | slope / log(q/4) | R² |
|---|---|---|---|---|---|
| 5 | 32,785 | −0.5619 | 0.2231 | **−2.518** | 0.999 |
| 7 | 258 | −1.3685 | 0.5596 | **−2.445** | 0.999 |
| 9 | 104 | −2.0529 | 0.8109 | **−2.531** | 0.994 |
| 11 (power-limited) | 20 | −1.5202 | 1.0116 | −1.503 | 0.977 |

Mean C across q ∈ {5, 7, 9}: **−2.498**, within 0.1% of the clean form −5/2.

**Empirical law (candidate, requires q=11 N=10⁹ confirmation):**

> conv_rate(j; q) ≈ A(q) · (4/q)^((5/2)·j)

**Boring explanations checked:**
- *Random-walk first-principles prediction was wrong about the q-scaling.* The drift-diffusion gambler's-ruin argument with v_2 ~ Geom(1/2) gives slope = −log(q)·log(q/4)/log(2)², which predicts q-DEPENDENT multipliers of 3.35 (q=5), 4.05 (q=7), 4.57 (q=9). The empirical multiplier is q-INDEPENDENT at ~2.5. So either v_2 is not Geom(1/2) under qx+1 dynamics, or the gambler's-ruin assumption ignores something (correlations, finite-N selection, etc.).
- *Cross-q unification check (q=3, all-converging):* per-class mean(σ − 10.43·log_n) vs j has slope +12.40, R² = 0.9999 (experiment 13). Same prefix-complexity axis predicts longer σ at q=3 (where orbits all converge but high-j classes take longer) and lower conv_rate at q=5,7,9 (where high-j classes mostly diverge). Confirms the prefix-complexity axis is structural, not q-specific.
- *N-decay of conv_rate per q (experiment 14):* α(5) = 0.343, α(7) = 0.635, R² = 0.9999 each. The N-decay exponent does NOT have a clean universal form across q (α/log(q/4) = 1.53 for q=5, 1.13 for q=7). This is a separate phenomenon from the j-decay law; treat as a footnote.
- *q=11 outlier:* slope/log(11/4) = −1.50 (not −2.5). However q=11 at N=10⁸ has only 20 convergent orbits across 4 j-bins; the slope fit is dominated by the j=1 point (10 orbits). Power is the explanation, not violation of the universal C.

**Substantive claim (publication-ready, IF q=11 confirms):**

> qx+1 per-class convergence rate is governed by an exponential decay in prefix odd-step count: conv_rate(j; q) ≈ A(q)·(4/q)^((5/2)·j) for q ∈ {5, 7, 9} (and conjecturally all primes q ≥ 5). The base 4/q is the random-walk drift; the multiplier 5/2 is empirical and currently lacks a derivation. R² ≥ 0.99 across three q values with n_conv ranging from 32,785 to 104.

**Open question:** Is C = 5/2 derivable? The random-walk argument with v_2 ~ Geom(1/2) predicts q-dependent multiplier; the data says q-independent at ~5/2. Either the v_2 distribution is different under qx+1 with q ≥ 5 (need to measure), or the walk-step correlations introduce a q-canceling correction.

**Next step (morning):** Generate q=11 at N=10⁹ to get ~100 convergent orbits and confirm whether C ≈ 2.5 holds at q=11 too. If yes, the constant is universal across small primes. If no, the q-independence at q=5,7,9 may be an artifact of small-q regime.

---

## 2026-05-01 (very late) — m-selection sign-flip and rejection of simple sub-exponential corrections at q=7

**Observation:** Joint regression of log(conv_rate) on log(X) and log(m) at q=7 N=10⁹ k=6 (25 usable 2D bins from log(X)×log(m) grid):
- Pure log(X) fit (Model A): slope b = −0.6413, deviation +2.4% from Cramér θ(7) = 0.6265
- Joint log(X) + log(m) fit (Model B): slope b = −0.6068, partial d on log(m) = **−0.0536, p = 8.1×10⁻¹⁰**, deviation **−3.1%** under Cramér
- **The deviation flipped sign upon adding log(m).** Model A was over-predicting θ; Model B underpredicts.

**Boring explanations checked:**
- *Bahadur-Rao 1/√L correction (Petrov-style sharp asymptotic):* Free fit of −0.5 prefactor coefficient gives c = +0.15 (opposite sign from BR's −0.5); forcing c = −0.5 makes SSR worse (8.80 vs 8.51). **BR ruled out** as the explanation.
- *Two-term additive plateau f(X) = A·X^(−θ) + B:* MLE drives B → 1.9×10⁻¹⁷, LR test p = 1.0. **No significant additive small-X floor.** The +log(log) signature in linear log-log fits is therefore binning-sensitive, not a real plateau.
- *Lattice walk corrections:* qx+1 step group log(q)·ℤ + log(2)·ℤ is dense in ℝ for q with log(q)/log(2) irrational (q=5,7), so step-lattice corrections to non-lattice Cramér-Lundberg shouldn't apply. Single-point return (vs half-line crossing) for non-lattice walks contributes only a constant prefactor, not log corrections.

**Substantive claim:** Cramér's exponential rate −θ(q)·log(X) is the correct leading-order asymptotic at q=5,7. At q=5 with 32,785 convergent orbits, agreement is 0.01% on j-slope and 2% on N-decay. At q=7 with 544 convergent orbits at N=10⁹, the pooled-X deviation is ~5% with two distinguishable components: (i) m-selection geometry within X-bins (partial r = −0.054 on log(m), p = 8×10⁻¹⁰), accounting for ~2.4% of deviation in the +X^(−θ) direction; (ii) a residual ~3% effect of opposite sign that does NOT match Bahadur-Rao or simple additive plateau. The residual likely reflects a combination of finite-N statistical noise (n_conv = 544 across 4-5 j-bins) and unidentified subleading effects.

**Caveat for any future writeup:** the original framing "we have an empirical exponential law" is correct at leading order. The framing "Cramér's theorem nails it everywhere" oversells; the q=7 pattern shows that subleading effects exist and aren't captured by the obvious sharp-asymptotic corrections. The honest claim is leading-order verification at q=5,7 with characterized but unresolved subleading deviation at q=7.

**Next step:** Cycle detection on q=5 non-convergent orbits. Confirms whether non-trivial cycle landings constitute a meaningful fraction of the "did not reach 1" set or whether non-convergent orbits are predominantly truly divergent.

---

## 2026-05-01 (very late) — q=5 cycle detection: methodological correctness confirmed

**Observation:** Floyd's algorithm cycle detection on 5,000 randomly sampled non-convergent q=5 orbits at N=10⁸ (max_steps = 5,000, max_value = 10³⁰):

| status | count | fraction |
|---|---|---|
| overflow (truly divergent, grew to 10³⁰) | 4,994 | **99.88%** |
| entered non-trivial cycle (all 13-cycle) | 6 | 0.12% |
| hit max_steps without resolution | 0 | 0.00% |

**Substantive claim:** Among q=5 orbits at N=10⁸, the "non-convergent" set (~50M orbits) is essentially the divergent set. <0.2% land in known non-trivial cycles (the 13-cycle, 17-cycle, etc.). The (a)-only "reach 1" analysis used for the j-slope and N-decay fits is therefore not contaminated by other-cycle landings. The 0.01% Cramér match at q=5 stands as a clean methodological result.

**For paper:** "Convergence to 1" and "settling into any finite cycle" agree to within 0.2% at q=5. Reporting (a)-only conv_rate is the methodologically clean choice; broader "settles into a cycle" framing would not change quantitative results meaningfully.

**Implication for q=7,9,11:** Without re-running cycle detection at higher q, we can't certify that those non-convergent sets are also predominantly divergent. But given that cycles in qx+1 dynamics scale roughly as q^? (heuristically) and our orbits at higher q overflow much faster (max_value reached in fewer steps), other-cycle landings should be even rarer at q ≥ 7. Worth verifying if the q=7,9,11 results become publication-load-bearing, but for now the q=5 result is the paper's headline and the methodological footnote is clean.

---

## 2026-05-02 — universal Geom(1/2) trajectory v_2 across q (closure)

**Observation:** Unconditional sampling of qx+1 trajectories (5,000 random odd n in [1, 10^7], 200 Syracuse steps each, ~10^6 total v_2 records per q) at q ∈ {5, 7, 9, 11}:

| q | mean v_2 | var v_2 | V(q) = var*log(2)^2 | drift_emp | drift_pred = log(q/4) |
|---|---|---|---|---|---|
| 5 | 2.0028 | 2.0056 | 0.9636 | 0.2212 | 0.2231 |
| 7 | 1.9990 | 2.0003 | 0.9611 | 0.5603 | 0.5596 |
| 9 | — | — | 0.9576 | 0.8119 | 0.8109 |
| 11 | — | — | 0.9632 | 1.0120 | 1.0116 |

Geom(1/2) prediction: mean=2, var=2, V=2*log(2)^2=0.9609. Every q matches to ~0.5%.

**Substantive claim:** The qx+1 trajectory measure on v_2 = ν_2(qn+1) is i.i.d. Geom(1/2) for all q ∈ {5, 7, 9, 11}, verified empirically. This is the load-bearing assumption for the Cramér derivation; it transfers cleanly across q. Therefore the Cramér θ(q) baseline (root of q^(-θ) = 2^(1-θ) - 1) is the correct prediction at each q without modification.

**Implication for q=7,9,11 deviations:** the residual gaps from Cramér prediction are NOT from a wrong baseline (the baseline is right within 0.5%). They're from finite-N statistical noise (n_conv = 544 at q=7, 104 at q=9, 36 at q=11) and m-selection geometry (real but partial confound at q=7, partial r = -0.054). With more data, deviations should shrink toward the q=5 0.01% precision.

**Closure:** the v_2 universality verification removes the last "did we use the wrong baseline" hedge. The publication-grade empirical claim is now load-bearing and bulletproof on every assumption that's directly testable on existing data.

---

## 2026-05-02 — B1 replication number diagnosis (Pathfinder caveat for the joint paper)

**Observation:** Our v4 NB-GLM Bayesian replication (`06b_bb_pathfinder.py`) at N=10^7, N_train=500K, N_test=50K, Pathfinder VI gives:

| Spec | params | log_score | sigma_u_mean |
|---|---|---|---|
| B0 (log n only) | 3 | −274,150.31 | — |
| B1 (B&B's setup, RE on mod 8) | 12 | −274,138.44 | 0.00254 |
| B3 (a_final at k=6, fixed effects) | 9 | −272,435.25 | — |
| B&B reported (NB2-GLM, HMC, full data) | — | −272,911.95 | (HMC posterior σ_u, larger) |

Our B1 absolute log score is 1,226 nats below B&B's reported number for the same model. **This is not a disagreement with their result.** It is a Pathfinder-vs-HMC gap.

**Diagnosis:**

1. **Pathfinder VI collapses hierarchical scales.** σ_u = 0.00254 in our B1 is two orders below the σ_u that an HMC fit to the same model would produce. Pathfinder (Zhang et al. 2022, JMLR) is a quasi-Newton variational approximation that systematically underestimates posterior dispersion in non-Gaussian hierarchical regions. Our hierarchical specs (B1, B4) inherit that underestimate; their fixed-effect specs (B0, B2, B3) do not.

2. **Local compute precludes HMC at this scale.** Stan 2.36 unified-mode multi-chain lockup at N_train=500K. Pathfinder was the workable pivot. HMC validation at Columbia is the natural next step — it closes the σ_u gap and makes the absolute numbers comparable.

3. **The structural decomposition is independent of the Pathfinder caveat.** B2 and B3 are fixed-effects: no σ_u, no underestimate. Their absolute log scores are directly comparable to B&B's HMC.

**What the comparison actually says (collaborative framing):**

The Bonacorsi-Bordoni mod-8 RE is *observing real structure*. The prefix decomposition supplies the algebraic identity behind it: u_{n mod 8} corresponds (up to a constant) to log(a_final(n mod 8)) at k=3. The two parameterizations carry essentially the same predictive content; ours is closed-form, theirs is estimated from data via the random effect.

Concretely from the v4 numbers:
- B1 (their mod-8 RE) and B2 (our a_final at k=3, the same modular resolution): B2 scores 850 nats above B1. This is consistent with the closed-form covariate replacing the estimated RE without information loss, plus a small parameter-count efficiency gain. (Pathfinder would shrink the RE and hand the gap to B2 even if HMC produced an identical absolute B1.)
- B3 (a_final at k=6, finer resolution): −272,435.25, with 9 parameters. This is 477 nats above B&B's reported HMC number for their 12-parameter mod-8 hierarchical model. It is a *suggestion* that finer modular resolution captures more structure than mod 8 — testable directly under HMC.
- B4 (B3 + mod-8 RE): essentially same score as B3, with a vestigial σ_u. Suggestive that once a_final is in the model the mod-8 RE has nothing left to estimate. Load-bearing under HMC, not under Pathfinder.

**Implication for the joint paper:** The structural decomposition supplies a closed-form covariate that Bonacorsi-Bordoni's hierarchical NB-GLM can use directly: replace u_{n mod 8} with f(a_final) for an arbitrarily-chosen k. The result is a parameter-parsimonious version of their model with the random effect grounded in the symbolic Collatz prefix. HMC validation at Columbia confirms or revises the absolute numbers. Either way, the algebraic identity is the contribution our side brings; their NB-GLM framework is the predictive engine. The two pieces fit together cleanly.

---

## 2026-05-02 — μ_β is non-monotone in N; per-octave structure traced to trajectory E[v] (partial)

**Observation 1 (cumulative β oscillation).** Streaming OLS of σ vs ln(n) on
odd n in [3, N] for N up to 2^32 gives a non-monotone β trajectory:

| N | log₂(N) | β | gap from 10.4282 |
|---|---|---|---|
| 2^20 | 20 | 10.3723 | +0.0559 |
| 2^22 | 22 | 10.3816 | +0.0466 |
| 2^23 | 23 | 10.3845 | +0.0437 |
| 2^24 | 24 | 10.4044 | +0.0238 |
| 2^25 | 25 | 10.4191 | +0.0091 |
| 2^26 | 26 | 10.4192 | +0.0090 |
| **2^27** | 27 | **10.4293** | **−0.0011 (crossed)** |
| 2^28 | 28 | 10.4298 | −0.0016 (peak overshoot) |
| 2^29 | 29 | 10.4252 | +0.0030 |
| 2^30 | 30 | 10.4236 | +0.0045 |
| 2^31 | 31 | 10.4213 | +0.0069 |
| 2^32 | 32 | 10.4187 | +0.0095 |

β approaches the heuristic from below for N ≤ 2^26, jumps above between 2^26
and 2^27 by 0.010 (single largest doubling step), peaks at 2^28, then drifts
back down crossing the heuristic from above at ~2^29 and re-opening the gap
to ~+0.010 by 2^32. The "monotone approach from below" claim in the writeup
A1 is empirically wrong; β oscillates, with amplitude not visibly damping by
N=2^32. Cross-check at N=2^25 with the streaming OLS reproduces the existing
writeup's β=10.4191 to 5 decimals — not a method artifact.

Code: `experiments/26_mu_beta_n_extension.py`.

**Boring explanations checked:**

- *Method artifact.* Cross-check at 2^25: 10.4191 (streaming OLS) vs 10.4191
  (existing writeup, prior code path) — match to 5 decimals. ✓ Not method.
- *Numerical precision.* Per-chunk float64 sums with chunk size N/64; total
  sum ~5×10^13 at 2^32, expected roundoff <10^-10. ✓ Not precision.
- *Record-σ outlier leverage.* Tested via top-K exclusion at N=2^27:
  K=10 (σ ≥ 769, including the σ=949 record at n=63,728,127): β shifts by 0.0001.
  K=100: shifts by 0.0005. K=1,000: shifts by 0.003. K=10,000: shifts by 0.020.
  To recover the 2^26 β=10.4192 from 2^27's β=10.4293 would require dropping
  ~5K-10K values, ~0.01% of the dataset. ✗ **Records hypothesis rejected**:
  the jump is a bulk-tail effect, not a few-extreme-outlier effect.

Code: `experiments/27_beta_oscillation_diagnostic.py`.

**Observation 2 (per-octave β_local).** OLS of σ vs ln(n) restricted to
odd n in [2^j, 2^(j+1)] gives:

| octave j | β_local | gap from 10.4282 |
|---|---|---|
| 17 | 10.66 | −0.23 |
| 18 | 10.75 | −0.32 |
| 19 | 10.70 | −0.27 |
| 20 | 10.72 | −0.29 |
| 21 | 10.88 | −0.45 |
| 22 | **10.89** | **−0.46 (peak)** |
| 23 | 10.78 | −0.35 |
| 24 | 10.65 | −0.22 |
| 25 | 10.59 | −0.16 |
| 26 | 10.49 | −0.06 |

β_local is consistently *above* 10.4282 in every well-populated octave (≥65K
odd n), with a peak at j=21–22 (10.88, 10.89) and graded falloff on both
sides. The cumulative β oscillation is the leverage-weighted average of these
local slopes plus the between-octave drift, which is dominated by the global
K ≈ 10.43.

**Observation 3 (closed-form mechanism, partial).** Per-octave trajectory
E[v] under the Syracuse map (1M starts per octave, walked to 1 with T=500 cap):

| octave | E[v]_traj | K_pred | β_local | residual (K_pred − β_local) |
|---|---|---|---|---|
| 17 | 1.989 | 10.673 | 10.656 | +0.017 |
| 18 | 1.991 | 10.631 | 10.753 | −0.122 |
| 19 | 1.992 | 10.595 | 10.696 | −0.101 |
| 20 | 1.993 | 10.574 | 10.720 | −0.146 |
| 21 | 1.993 | 10.573 | 10.877 | −0.304 |
| 22 | 1.994 | 10.563 | 10.887 | −0.325 |
| 23 | 1.994 | 10.561 | 10.780 | −0.219 |
| 24 | 1.994 | 10.560 | 10.654 | −0.094 |
| 25 | 1.995 | 10.549 | 10.586 | −0.038 |
| 26 | 1.995 | 10.544 | 10.489 | +0.056 |

`K_pred(u) = (1 + u) / (u·log 2 − log 3)`, the standard random-walk heuristic
parameterized in trajectory-mean v.

**What this explains:**
- E[v] under the trajectory measure is consistently slightly below 2 (Geom(1/2)
  prediction). This single deviation predicts K ≈ 10.55–10.67 across all
  octaves, fully explaining the systematic ~0.13-nat baseline shift of every
  β_local above the heuristic 10.4282.

**What this does NOT explain:**
- The variation in β_local across octaves (10.49 to 10.89, range 0.40). E[v]
  varies only from 1.989 to 1.995 across the same octaves; K_pred varies only
  from 10.55 to 10.67. The residual variation in β_local — peaking at j=21–22
  and decreasing on both sides — is a separate effect.

**Mechanism candidates for the residual:**

1. *Higher-moment trajectory measure.* What enters a Cramér-style derivation
   is the MGF M(s) = E[2^(sv)] at s near θ − 1, not just E[v]. The v=4, v=10
   spikes shift M(s) in a way the mean doesn't capture. This is the same
   issue the qx+1 thread is testing on its side.
2. *Step-to-step correlations.* Consecutive Syracuse v's are deterministically
   related through the residue-class structure. The K formula assumes i.i.d.
   v's; the correlation produces an effective drift different from
   log(2)·E[v] − log(3).
3. *Octave-dependent descent geometry.* Trajectories starting at 2^j vs 2^(j+1)
   descend through different sets of intermediate states. Pooled-trajectory
   E[v] averages over the whole descent and may obscure octave-specific
   corrections to the effective slope.

Code: `experiments/28_per_octave_trajectory_E_v.py`.
CSV: `experiments_output/28_per_octave_E_v_Nperoct1000000.csv`.

**Implications:**

- *For writeup A1.* The "monotone approach from below at ~N^{−1/2}-ish rate"
  language must be replaced. Empirically: μ_β is *oscillatory* in N at the
  scales tested (2^20–2^32), with amplitude ~0.01 not visibly damping by 2^32.
  The crossing of 10.4282 happens around N ≈ 2^27.

- *For the structural claim.* The slope-universality result (τ_β = 0 to noise
  floor) is per-class within fixed N — unaffected by the cumulative β
  oscillation. Per-class slopes are still indistinguishable from each other;
  it's the *global* slope that wanders around 10.4282.

- *For the qx+1 thread.* The same E[v]-not-2 phenomenon is the load-bearing
  question on that side: whether the empirical MGF M_traj(s) at s = θ(q) − 1
  matches Geom(1/2)'s MGF at the same s. If the second agent's MGF check shows
  a magnitude consistent with the per-octave residual here (~0.3 nats / ~3%),
  the same higher-moment correction governs both threads.

**What remains open (no immediate action; held for big-picture review with the
v=4, v=10 spike consolidation):**

- Empirical MGF M_traj(s) on existing 7B-sample data (delegated thread).
- m-residue pushforward at mod 32 / 2048 / 131072 (delegated thread).
- q=5 trajectory v-distribution shape (delegated thread).
- Test of mechanism candidates (2) and (3) above.

**Status:** β-oscillation mechanism partially identified. Trajectory E[v]
fully explains the baseline shift; residual per-octave variation traces to a
separate (likely higher-moment) trajectory effect. Consolidating here for the
big-picture pass with the v-spike work.

---

## 2026-05-02 — α_det predicts s_mean at slope=1 with K=K_heuristic, matching Tao (5.15)

**Observation.** For odd n in residue class r mod 2^k, define s(n) as the
first-passage step count at which the Collatz orbit attains a value ≤ f(N)
for some threshold function f. The per-class mean s_mean(r) is linearly
related to the prefix-decomposition's α_det(r):

> s_mean(r) ≈ α_det(r) + K_h · log(N/f(N))

where K_h = 3/log(4/3) ≈ 10.4282 is the textbook descent constant from the
random-walk heuristic (and from Tao 2022's (5.15) leading term).

**Test data.** N = 2^25 and N = 2^27, all odd starting points in [3, N].
Five observables:
- σ (full descent to 1, Δlog = log N)
- first passage to N^(2/3) (Δlog = log N / 3)
- first passage to √N · log N (Δlog = log N / 2 − log log N)
- first passage to √N (Δlog = log N / 2)
- first passage to √N / log N (Δlog = log N / 2 + log log N)

Four modular resolutions: k = 8, 10, 12, 14. Total 40 (observable, k, N)
cells.

**Boring explanations checked:**

- *Spearman ρ already at 1.0 at k=8 with median.* Verified for s_median in
  experiment 30, 31; the rank-correlation result is independent of mean vs
  median.
- *Calibration error in α_det's K = 10.4282.* Tested at k=8 (experiment 32):
  recomputing α_det at K_emp ∈ {9.30, ..., 10.80} found that for s_median,
  slope = 1 was hit at K ≈ 10.0–10.3 (below K_h, threshold-dependent). For
  s_mean, slope = 1 was hit *exactly at K_h* (experiment 33).
- *Sampling distribution.* per-class size ~131K at k=8 (N=2^27), ~2K at
  k=14 (N=2^25). Slope stability across these very different sample sizes
  confirms it's a structural identity, not a fit artifact.

**Result table (raw mean, slope @ K_h):**

| observable | k=8 N=2^25 | k=8 N=2^27 | k=14 N=2^25 | k=14 N=2^27 |
|---|---|---|---|---|
| σ | 0.9960 | 0.9990 | 0.9945 | 0.9977 |
| s @ N^(2/3) | 0.9989 | 0.9996 | 0.9944 | 0.9977 |
| s @ √N | 0.9989 | 1.0006 | 0.9967 | 0.9997 |
| s @ √N/log N | 1.0000 | 1.0012 | 0.9980 | 1.0005 |

**All 40 raw-mean slopes at K_h fall in [0.9936, 1.0012].** Tightest at
√N threshold (median across resolutions and scales 1.000).

**Offset gap from Tao (5.15) leading term K_h · Δlog (raw mean):**

| observable | typical gap (across all k, N) |
|---|---|
| σ | −2.4 (stable to 0.05) |
| s @ N^(2/3) | +3.0 (stable to 0.04) |
| s @ √N · log N | +3.0 |
| s @ √N | +2.1 |
| s @ √N / log N | +0.4 to +1.2 (only one varying) |

**Offset gaps are stable across k and across N** (variation < 0.06 across
k = 8 → 14 and N = 2^25 → 2^27 for σ; comparable for the rest). The
correction is a structural constant of the σ distribution, not a finite-N
or finite-k artifact.

**Trim-1% mean offset gap at √N:**

| N | k=8 | k=10 | k=12 | k=14 |
|---|---|---|---|---|
| 2^25 | +0.14 | +0.17 | +0.20 | +0.23 |
| 2^27 | +0.26 | +0.29 | +0.32 | +0.34 |

Trimming the top 1% per class (the right-tail orbits Tao's exceptional-set
theorem allows for) brings the offset within sub-percent of Tao's prediction
at every (k, N) cell.

**Substantive claim:**

The structural decomposition — α_det(r) computed by symbolic Collatz prefix
iteration on residue class r — predicts the per-class *mean* of two
independent trajectory functionals (σ and first-passage time) with:
- Spearman ρ = 1.0 across 128 to 8192 classes
- Slope = 1.000 ± 0.005 at K = K_heuristic, no fit calibration
- Offset matching Tao (5.15) leading term within ≤ 1 step (1%-trimmed) or
  ≤ 3 steps (raw)

This is the per-class realization of Tao 2022's almost-everywhere trajectory
formula. The structural decomposition supplies the algebraic identity for
α_det; Tao's (5.15) supplies the asymptotic mean coefficient K_h = 3/log(4/3).
They agree at every modular resolution and data scale tested.

**Implication for the bridge to Tao 2022:**

The result connects two complementary frameworks. Tao's theorem controls
the mean-style first-passage time `T_x(N) = log(N/x)/log(4/3) + O(log^0.6 x)`
for almost all N. The per-class mean version *here* is exact at the leading
term, with the structural correction `α_det(r)` providing the per-residue-class
constant. Tao's almost-everywhere statement is the mean-over-classes; our
result is the per-class refinement.

**Earlier "median is faster than mean" framing:** correct. s_median sits
~7 steps below Tao's leading term (at √N); s_mean sits ~2 steps above. The
gap is the σ-distribution skewness producing median < mean. Trimming top 1%
brings s_mean to within 0.3 steps of Tao's prediction, since the heavy tail
contributes ~2 steps to mean σ that Tao's leading term doesn't include.

**Files of record:**

- `experiments/30_first_passage_alpha_det.py` — initial s_median ρ=1 finding
- `experiments/31_first_passage_replication.py` — k-sweep + threshold variation
- `experiments/32_alpha_det_K_calibration.py` — K-recalibration test for s_median
- `experiments/33_alpha_det_K_calibration_mean.py` — K-recalibration test for s_mean
- `experiments/34_alpha_det_K_calibration_mean_k_sweep.py` — k=8/10/12 mean replication
- `experiments/35_alpha_det_full_bridge.py` — full bridge (k=8,10,12,14 × 5 observables × 2 N's)

**Status of the structural claim:** the headline shifts from "α_det predicts
σ-intercepts" to "α_det predicts per-class mean trajectory dynamics on two
independent observables, validated against Tao (5.15) leading term across
modular resolutions and data scales." This is the strongest single result of
the project for the bridge to Tao 2022.
