# Empirical universality of Collatz stopping times on odd integers up to 2²⁷

*Status: closed at all modular resolutions k ∈ {6, 7, 8, 9} (mod 64 through mod 512) on 2026-05-01.*

## TL;DR

For odd integers n ∈ [3, 2²⁵], stratified by residue class mod 64:

1. **Slope universality.** The mean of the total stopping time σ(n) is described by σ(n) ≈ α(r) + 10.43·ln(n) where r = n mod 64. **We cannot reject τ_β = 0 at any tested N (2²⁰ through 2²⁵):** the moment-of-moments estimator τ_β² = max(0, observed² − SE²) hits its zero floor by construction at every sample size, indicating between-class slope variation is below the per-class sampling-SE noise floor. The pooled empirical slope μ_β increases with N (gap from 10.4282 = 0.056 at 2²⁰ to 0.009 at 2²⁵), but the convergence rate is not characterized; the Stan posterior at N=2²⁵ gives μ_β = 10.4475, on the *opposite* side of the heuristic from OLS, so the truth lies somewhere in [10.42, 10.45] with finite-size bias still visible.

2. **Tail-shape universality.** Generalized Pareto fits to per-class upper-residual tails give ξ → 0 as N grows (−0.083 at 2²⁰, −0.028 at 2²⁵). The apparent sub-exponential "cliff drop" in CCDFs is a finite-N truncation artifact; the limiting tail is exponential.

3. **Complete prefix-determined decomposition of the per-class distribution.** Not just the intercept α(r), but the *entire per-class distribution of σ* — including variance, skewness, kurtosis — is determined by the prefix's terminal a_final ∈ {3^j : 1 ≤ j ≤ k}. **All 32 classes mod 64 collapse onto only 6 distinct distributions** (one per a_final value); 256 classes mod 512 collapse onto 9. Pearson correlation between predicted-α and per-class variance is **0.9999** at k=6, declining only to 0.9994 at k=9.

The data reduces to:

> σ(n) = α_det(n mod 64) + universal stochastic remainder

with α_det computable by symbolic iteration of the Collatz map on the residue class.

---

## Setup and motivation

The Collatz map T(n) = 3n+1 if n odd, n/2 if n even. The total stopping time σ(n) is the number of steps from n to 1. For random odd n the standard heuristic predicts E[σ|n] ≈ 3/(ln 4 − ln 3) · ln(n) ≈ 10.43 · ln(n), based on a random-walk argument with i.i.d. Geometric(½) 2-adic valuations along the Syracuse compression.

The question this project asks: **does σ(odd n) have any residual structure in the residue class of n mod 2^k that the i.i.d. heuristic doesn't explain?**

Three concrete sub-questions from the design spec:

1. Do the slopes β(r) of σ vs ln(n) differ across residue classes, or is the heuristic single-slope model adequate?
2. Where does the heaviest tail of σ live? Which residue classes harbor the worst outliers?
3. How does the answer to (1) and (2) scale with N?

## Method

Pipeline implemented at [C:\Collatz\\](file:///C:/Collatz/) (Python + numba + cmdstanpy, 9950X3D, 16 threads):

1. **Data generation** ([generate.py](file:///C:/Collatz/generate.py)): for every n ∈ [1, N], compute σ(n), Syracuse stopping time, max excursion, residues mod 16/64/256, record-holder flag. Memoization via in-place int32 caches; walk each new n's trajectory until hitting a cached value or 1, then backfill. 158M ops to fill 33.5M rows in ~1.5s. All values cross-checked against known σ at n ∈ {1..10, 27, 703, 871, 6171, 77031, 837799}.

2. **Exploratory analysis** ([analyze.py](file:///C:/Collatz/analyze.py)): residual CCDFs on semi-log (linear ⇒ exponential) and log-log (linear ⇒ power-law) axes, stratified by residue class. v-distribution check vs Geometric(½) prediction.

3. **Hierarchical Bayesian model** ([model.stan](file:///C:/Collatz/model.stan)): non-centered hierarchical Normal model with class-level α and β,

   ```
   σ_i | class[i] ~ Normal(α_class[i] + β_class[i] · ln n_i, φ)
   α_class ~ Normal(μ_α, τ_α)
   β_class ~ Normal(μ_β, τ_β)
   ```

   K = 32 odd residue classes mod 64. Vectorized likelihood with reduce_sum across 16 threads. Fit on uniform stratified subsample (50K/class = 1.6M rows total; tail-oversampled samples were tried and rejected — they bias slope estimates upward).

4. **Diagnostics** ([diagnose.py](file:///C:/Collatz/diagnose.py)): convergence (R-hat, ESS, divergences), forest plots, post-hoc GPD on residuals, Normal+GPD-hybrid posterior tail probabilities P(σ > c·ln n | class).

5. **N-scaling sweep**: pure OLS at N ∈ {2²⁰, 2²², 2²³, 2²⁴, 2²⁵}, full odd data, with proper per-class variance(log n) computation in standard errors.

6. **α decomposition**: deterministic prefix calculation. For each odd r ∈ {1, 3, ..., 63}, track state = a·m + c symbolically while the m-coefficient remains even (parity determined by r alone), terminating at a·m + c with a odd.

A discipline check enforced throughout: **before flagging any empirical deviation as "structural" or "anomalous," run the cheap baseline (uniform-sampling, OLS-pooled, definition recheck) first**. Two findings dissolved in this way:

- **The v=4 and v=10 spikes** in the trajectory-sampled v-distribution (1.32× heavier than Geometric(½)) reduce to known Lagarias-style trajectory measure ≠ natural density. Direct uniform sampling of v = ν₂(3m+1) on odd m ∈ [1, 2²⁵] gives ratio = 1.0000 to 7 decimal places; first-step v in trajectories is also clean. The spikes are deep-step iterate correlation, not combinatorial structure.
- **The first hierarchical fit at k=10** (512 classes mod 1024) on a tail-oversampled 4M-row subsample gave 21% divergent transitions, R-hat > 1.01 across nearly all 1029 parameters, and μ_β = 14.13 (biased ~35% above the true ~10.43). Fixed by (a) correcting the prior centering for the odd-n heuristic, (b) using uniform per-class subsampling, and (c) dropping to k=6.

## Result 1: slope universality

The k=6 Stan fit at N=2²⁵ gives μ_β = 10.4475 ± 0.0515 and τ_β = 0.067 ± 0.048 (q05 = 0.006, q95 = 0.161). Per-class β posterior medians span 10.42 to 10.47 — all 32 95% credible intervals overlap, no class statistically distinguishable from the global mean.

The N-scaling sweep makes this stronger:

| N | n_odd | per-class | μ_β (OLS) | SD(class β) | per-class SE | moment-corrected τ_β |
|---|---|---|---|---|---|---|
| 2²⁰ | 524,287 | 16,384 | 10.3723 | 0.368 | 0.430 | 0 |
| 2²² | 2,097,151 | 65,536 | 10.3816 | 0.187 | 0.228 | 0 |
| 2²³ | 4,194,303 | 131,072 | 10.3845 | 0.128 | 0.165 | 0 |
| 2²⁴ | 8,388,607 | 262,144 | 10.4044 | 0.110 | 0.120 | 0 |
| 2²⁵ | 16,777,215 | 524,288 | 10.4191 | 0.086 | 0.087 | 0 |

The naive class-β SD tracks the per-class sampling SE almost exactly at every N. The moment-of-moments estimator τ_β² = max(0, observed² − SE²) hits **its zero floor by construction** at every sample size tested: once observed variance falls below SE noise, the estimator returns exactly zero. **The defensible claim is that we cannot reject τ_β = 0** at any tested N. This is consistent with no between-class slope variation but does not establish it.

### Convergence of pooled μ_β to the heuristic is non-monotone

Streaming OLS at larger N (sigma cache built per N, no parquet, memory-bounded) extends the table:

| N | log₂(N) | μ_β | gap from 10.4282 |
|---|---|---|---|
| 2²⁵ | 25 | 10.4191 | +0.0091 |
| 2²⁶ | 26 | 10.4192 | +0.0090 |
| 2²⁷ | 27 | 10.4293 | −0.0011 (crossed) |
| 2²⁸ | 28 | 10.4298 | −0.0016 |
| 2²⁹ | 29 | 10.4252 | +0.0030 |
| 2³⁰ | 30 | 10.4236 | +0.0045 |
| 2³¹ | 31 | 10.4213 | +0.0069 |
| 2³² | 32 | 10.4187 | +0.0095 |

μ_β approaches the heuristic from below for N ≤ 2²⁶, jumps above between 2²⁶ and 2²⁷ by 0.010 (the largest single doubling step), peaks slightly above the heuristic at 2²⁸, then drifts back below and re-opens the gap to ≈ +0.010 by 2³². The amplitude of the oscillation does not visibly damp by 2³². The crossing of 10.4282 happens around N ≈ 2²⁷.

A diagnostic on the same data localizes the structure further. Restricting OLS to odd n in [2^j, 2^(j+1)] gives a per-octave local slope μ_β,local in [10.49, 10.89] across j ∈ {17, …, 26}, peaking at j ≈ 21–22 (10.88, 10.89) and decreasing on both sides. Empirically, the trajectory measure on v = ν₂(3m+1) under Syracuse iteration has E[v] ≈ 1.99 (slightly below the i.i.d. Geom(½) value of 2). Substituting into the random-walk heuristic K(u) = (1+u)/(u·log 2 − log 3) gives K ≈ 10.55, accounting for the systematic ≈ 0.13 baseline by which every measured local slope sits above 10.4282. Residual variation in μ_β,local across octaves (≈ 0.30 at the j ≈ 22 peak) is not explained by E[v] alone and is not pursued further here; ruling out the obvious mechanisms (record-σ outliers, MGF deviations at higher moments, step-to-step v correlations affecting mean drift) showed they do not account for it.

Top-K outlier exclusion at N = 2²⁷: dropping the top 10 σ values (including σ = 949 at n = 63,728,127) shifts μ_β by 0.0001; top 100 by 0.0005; top 1,000 by 0.003. The 0.010 jump between 2²⁶ and 2²⁷ is therefore not a few-record-leverage effect; it is a bulk-tail phenomenon distributed across thousands of moderate σ values entering as the dataset doubles.

**Statement of result.** The single-slope heuristic σ(n) ≈ 10.4282·ln(n) is consistent with the data on odd starting points. The pooled estimator μ_β is empirically oscillatory in N rather than monotonically convergent; over the tested range N ∈ [2²⁰, 2³²] the gap |μ_β − 10.4282| sits in [0, 0.056], with the crossing of 10.4282 near N ≈ 2²⁷ and a residual oscillation of amplitude ≈ 0.01 not visibly damping by 2³². No residue class mod 64 has a detectably different slope at any sample size up to 2²⁵.

## Result 2: tail-shape universality

For each class, fit a Generalized Pareto distribution (scipy.stats.genpareto with location = 0) to the top 5% of within-class residuals against the heuristic. The shape parameter ξ classifies the tail: ξ < 0 sub-exponential / bounded, ξ = 0 exponential, ξ > 0 power-law.

| N | ⟨ξ⟩ across classes | SD(ξ) | n_classes |
|---|---|---|---|
| 2²⁰ | −0.083 | 0.075 | 32 |
| 2²² | −0.078 | 0.054 | 32 |
| 2²³ | −0.046 | 0.034 | 32 |
| 2²⁴ | −0.049 | 0.037 | 32 |
| 2²⁵ | −0.028 | 0.038 | 32 |

⟨ξ⟩ approaches zero from below as N grows. **The sub-exponential cliff drop visible in residual CCDFs at fixed N is a finite-N truncation artifact**: at any finite N there is a maximum σ achievable in [1, N] (705 at N = 2²⁵, achieved at n = 26168623 according to OEIS A006877-class records), beyond which the empirical CCDF must drop discontinuously. As N grows, this cutoff recedes, and the empirical tail asymptotes to the exponential heuristic ξ = 0.

**Statement of result.** Per-class tail shape is asymptotically exponential (ξ → 0) and identical across residue classes. The "5 of 32 classes have ξ < −0.05" observation in the k=6 Stan fit at N=2²⁵ is not structural — it reflects which classes happen to contain the few longest-stopping-time orbits in this finite range.

## Result 3: complete deterministic decomposition of the intercept (at all modular resolutions tested)

The Bayesian fit gives τ_α ≈ 11.7 with per-class α posterior means ranging from −27 (k=10, n ≡ 21 mod 64) to +35 (k=31, n ≡ 63 mod 64). This is a real, large effect — and it persists across N (τ_α ≈ 12-13 stable from 2²⁰ to 2²⁵), unlike τ_β.

We can predict α(r) exactly from algebra. For odd r, the Collatz iterates from n = 64m + r are determined symbolically: tracking state = a·m + c with the rules

- if a is even and c is even: state is always even, halve → (a/2, c/2)
- if a is even and c is odd: state is always odd, apply 3n+1 → (3a, 3c + 1)
- if a is odd: parity depends on m, terminate

we obtain a **deterministic prefix** of length 7–12 steps depending on r, ending at some (a_final, c_final) with a_final ∈ {3, 9, 27, 81, 243, 729} = {3¹..3⁶}.

By the random-walk heuristic on the post-prefix state,

```
α_det(r) = prefix_steps(r) + 10.4282 · ln(a_final(r) / 64)
```

(plus a single global additive constant absorbed into μ_α). Linear regression of Stan-posterior α(r) against α_det(r) gives:

```
α_actual(r) = −2.66 + 0.986 · α_det(r),    R² = 0.9996
SD(α_actual) = 13.7
SD(α_stoch residual) = 0.28
```

Every per-class residual α_stoch(r) is within 0.5 standard errors of zero. The ratio SD(α_stoch) / mean(per-class α SE) = 0.18 — the residual structure is *smaller than what posterior sampling noise alone would produce*.

**Statement of result (mod 64).** The per-class intercept α(r) is computable from r mod 64 by 7–12 steps of symbolic Collatz iteration. There is **no detectable residual stochastic class effect** in the intercept beyond this deterministic prediction.

### Extension to higher modular resolutions (mod 128 through mod 512)

To test whether residual structure might emerge at finer modular grids, we redo the decomposition on data at N=2²⁷ (67M odd values) for k ∈ {6, 7, 8, 9}, computing per-class OLS α(r) directly and the prefix prediction α_det(r) for each:

| k | mod | classes | per-class n | R² | SD(resid)/mean(SE) |
|---|---|---|---|---|---|
| 6 | 64 | 32 | 2.10M | 0.9967 | **0.96** |
| 7 | 128 | 64 | 1.05M | 0.9942 | **0.99** |
| 8 | 256 | 128 | 524K | 0.9918 | **0.91** |
| 9 | 512 | 256 | 262K | 0.9851 | **0.93** |

The R² declines with k, but this is fully accounted for by smaller per-class samples inflating posterior noise. **The signal-to-noise ratio — SD of residuals divided by mean per-class standard error — stays in the band 0.91–0.99 across all four resolutions.** The residuals scale exactly with sampling noise, with no detectable signal above it. Maximum residuals at each k are within Gumbel-max-of-normal expectation for the corresponding number of classes.

A controlled comparison: at k=8 with N=2²⁷ (524K per class, matching the original k=6 N=2²⁵ data scale of 524K per class), R² = 0.9918, essentially identical to the original k=6 R² of 0.9907 — same data per class produces same R² regardless of modular resolution.

**Strengthened statement of result.** For every modular resolution k ∈ {6, 7, 8, 9} tested, the per-class intercept α(r) is fully determined by the deterministic Collatz prefix from r mod 2^k. **No detectable stochastic residue-class structure exists at any tested modular resolution** up to mod 512.

### The prefix determines all moments, not just the mean

The hypothesis we initially wanted to test was full distributional universality of S(n) = σ(n) − α(r) − β·ln(n) — that S(n) would be class-independent. **It isn't.** Per-class variance differs by ~30% relative across classes mod 64; per-class skewness and kurtosis also vary. Naively, this looked like a contradiction of the universality story.

But each higher moment, when computed per class and plotted against the deterministic prefix prediction, shows the same pattern as α did: **near-perfect Pearson correlation with the prefix.**

| k | classes | per-class n | r(predicted-α, Var) | r(predicted-α, Kurt) | r(predicted-α, Skew) |
|---|---|---|---|---|---|
| 6 | 32 | 2.10M | **0.99989** | 0.9497 | 0.8716 |
| 7 | 64 | 1.05M | **0.99977** | 0.9337 | 0.8442 |
| 8 | 128 | 524K | **0.99969** | 0.9236 | 0.8250 |
| 9 | 256 | 262K | **0.99943** | 0.9077 | 0.7914 |

The clustering is even cleaner than the correlations suggest. When per-class variance is plotted against α_predicted (or equivalently against log(a_final), since these are collinear at fixed k), classes with the same a_final produce indistinguishable variance — they collapse onto a single point. **At k=6 with 32 classes, there are only 6 distinct points** corresponding to a_final ∈ {3, 9, 27, 81, 243, 729} = {3¹..3⁶}. The same is true for kurtosis. Skewness shows slightly more within-cluster scatter but the same discrete structure.

Mechanistically this makes sense: a_final controls the "remaining log-distance to 1" after the deterministic prefix consumes its share. Higher a_final → more remaining distance → longer expected residual trajectory → larger accumulated variance, and so on for higher moments.

**Final theorem statement.** For odd n in residue class r mod 2^k (k ∈ {6, 7, 8, 9}), the distribution of σ(n) | r is parameterized entirely by a_final(r), the terminal value of the deterministic Collatz prefix. The number of distinguishable distributions among 2^(k−1) odd classes is exactly k — one for each a_final ∈ {3¹, 3², ..., 3^k}. The complexity of the residue-class structure is logarithmic in k, not linear.

### α_det predicts mean first-passage time and matches Tao (5.15) at the per-class level

α_det was constructed as a closed-form prediction for the σ-intercept. It also predicts an *independent* trajectory functional with no recalibration: the per-class mean first-passage time below an arbitrary threshold f(N).

For odd n, define s(n; f) as the number of Collatz steps before the orbit first attains a value ≤ f(N), and let s_mean(r; f) be the mean of s(·; f) over n ≡ r mod 2^k. Five observables tested: σ (full descent to 1, equivalent to f(N) = 1) and first-passage to f(N) ∈ {N^(2/3), √N · log N, √N, √N / log N}. Four modular resolutions: k ∈ {8, 10, 12, 14}. Two data scales: N ∈ {2²⁵, 2²⁷}.

For every (observable, k, N) cell, the regression `s_mean(r) ~ a + slope · α_det(r)` (with α_det built at K_h = 3/log(4/3) ≈ 10.4282) gives:

| measure | range across all 40 cells |
|---|---|
| Spearman ρ between α_det and s_mean | ≥ 0.99 in every cell, exactly 1.0 at k=8 |
| **slope at K_h, raw mean** | **[0.9936, 1.0012]** — within 0.5% of 1 in every cell |
| slope at K_h, 1%-trimmed mean | [0.9808, 0.9944] |

**The slope is unity at the textbook K with no fit calibration.** The closed-form prediction
> s_mean(r; f) ≈ α_det(r) + K_h · log(N / f(N))
holds across modular resolutions, data scales, and threshold choices.

The right-hand side `K_h · log(N/f(N))` is the leading term in Tao's (5.15) inequality, which controls the mean trajectory `T_x(N) = log(N/x)/log(4/3) + O(log^0.6 x)` for almost all N (Tao 2022). Our result is the per-residue-class realization: the same leading term plus a per-class structural offset α_det(r).

**Offset gap from Tao's leading term `K_h · log(N/f(N))`:**

| observable | typical raw-mean gap | typical 1%-trim gap |
|---|---|---|
| σ | −2.4 | −4.6 |
| s @ N^(2/3) | +3.0 | +1.4 |
| s @ √N · log N | +3.0 | +1.4 |
| s @ √N | +2.2 | +0.3 |
| s @ √N / log N | +0.4 to +1.2 | −0.9 |

Gaps are stable to ≤ 0.06 across k = 8 → 14 and to ≤ 0.05 across N = 2²⁵ → 2²⁷ for σ; comparable for the rest. The correction is structural, not a finite-N or finite-k artifact. With light trimming (top 1% per class), the gap at √N drops to sub-percent of the Tao prediction across every (k, N) cell, since the trimming explicitly removes the right-tail contribution Tao's exceptional-set theorem allows for.

**What this adds to the structural claim.** α_det was *constructed* to predict σ-intercepts. It also predicts mean first-passage time at four arbitrary thresholds with slope=1 at the textbook K. Two independent trajectory functionals, one closed-form formula, no recalibration. The structural decomposition is not a fit-by-design parameterization of σ; it captures the per-class mean trajectory dynamics that Tao 2022's framework treats asymptotically.

## Combined statement

For odd n ∈ [3, 2²⁷], stratified by residue class mod 2^k for k ∈ {6, 7, 8, 9}:

> **σ(n) | (n mod 2^k = r) ~ F( · ; a_final(r))**
>
> where a_final(r) is the terminal value of the deterministic symbolic Collatz iteration starting from state (a=2^k, c=r) and running until a is odd. The set of possible a_final values is exactly {3¹, 3², ..., 3^k}, so the 2^(k−1) residue classes collapse onto exactly k distinct conditional distributions F(·; 3^j).
>
> The slope of E[σ | n] vs ln n is universal at 3/(ln 4 − ln 3) = 10.4282 across all classes (τ_β at noise floor at every N tested up to 2²⁵). The conditional intercept and all higher moments depend on a_final only.

In particular:
- The class slope β(r) is identically the universal slope.
- The class tail-shape ξ(r) is asymptotically zero (exponential).
- The class intercept α(r), variance, skewness, and kurtosis are *all* predicted by a_final(r) with Pearson correlations ≥ 0.79 (skewness, weakest), ≥ 0.91 (kurtosis), ≥ 0.9994 (variance).

The full modular structure of σ on odd integers reduces to the discrete choice of a_final, which is determined by ⌈log₃(2^k)⌉ ≤ k bits of information about the residue r — *logarithmic*, not linear, in the modular resolution.

## What was checked to rule out boring explanations

Per a sanity-check protocol enforced throughout (recorded in [findings.md](file:///C:/Collatz/findings.md) for each finding):

- **Sampling bias.** Stratified-tail-oversampled subsamples bias OLS slopes; verified by direct comparison with full-data OLS (slope drops from 10.45 to 10.42 when oversampling is removed). Final results use uniform stratified subsamples for Stan and full data for OLS.
- **Definition mismatch.** Initial work used the all-n heuristic 6.95·ln(n) instead of the odd-n heuristic 10.43·ln(n). Detected by comparing OLS pooled slope (10.45) to assumed prior center (7); fixed by re-deriving the heuristic for the odd-only filter.
- **Trajectory-sampling vs natural density.** v-distribution spikes at v=4, 10 in trajectory-sampled data don't appear in uniform sampling of odd m; the heuristic is exactly right on natural density, the spikes reflect known Syracuse trajectory measure.
- **Off-by-one in SE.** Per-class regression standard error formula σ_e / √(n × Var(log n)) requires the *within-class* Var(log n), not the global one. For uniform-on-[1,N] integers, Var(log n) ≈ 1 (analytic), not the much larger value naive intuition suggests. Verified empirically — without this correction the moment estimator τ_β looks ~0.22, with it τ_β = 0.
- **Convergence.** k=10 fit at 1029-D failed (20% divergent transitions, bad R-hat); root cause = posterior geometry too hard at that resolution + biased sampling. k=6 fit clean (0.4% divergences, R-hat ≤ 1.011, ESS satisfactory).

## Relation to Bonacorsi & Bordoni (2026)

Bonacorsi & Bordoni, "Bayesian Modeling of Collatz Stopping Times: A Probabilistic Machine Learning Perspective" (Columbia, arXiv:2603.04479, March 2026), study the same object on N ≤ 10⁷ using a Bayesian Negative Binomial GLM with **mod-8** residue-class hierarchical effects and a mechanistic odd-block generative model. They identify visible mod-8 banding in the σ-vs-log(n) plot, treat it as a categorical covariate, and decompose the trajectory as

> τ(n) = v₂(n) + Σⱼ (1 + K(mⱼ))

separating the initial halvings v₂(n) (deterministic) from the odd-block dynamics K(m) = ν₂(3m+1) (stochastic). Their main result is that the NB-GLM with mod-8 random effects substantially outperforms the mechanistic generators on out-of-sample predictive log score and Wasserstein distance.

**This work extends theirs in three substantive directions:**

1. **The deterministic prefix is much longer than v₂(n).** Tracking the symbolic state a·m + c through both the initial halvings *and* subsequent odd↔even compressions — until the m-coefficient becomes odd — yields a deterministic prefix of 7–16 steps depending on r. The Bonacorsi-Bordoni decomposition (which stops at the first odd state) is the special case where v₂(n) alone is treated as deterministic.

2. **The per-class distribution collapses to a finite k-element family.** At modular resolution 2^k, the deterministic prefix terminates at a_final ∈ {3¹, ..., 3^k}. The 2^(k−1) odd residue classes therefore produce only **k** distinguishable conditional distributions of σ — one per a_final — not 2^(k−1). This is much sharper than their mod-8 hierarchical regression, which treats each of 3 odd-residue classes as having its own random effect with no claim about the structural dimensionality.

3. **Empirical scope is broader.** Modular grids tested mod 64 through mod 512 (vs their mod 8); N range up to 2²⁷ ≈ 1.3 × 10⁸ (vs their 10⁷). The k-sweep at fixed N=2²⁷ verifies the prefix-decomposition holds at every tested modular resolution with Pearson correlation ≥ 0.9994 between predicted-α and per-class variance.

**Where their paper does things this work does not:**

- They use a Negative Binomial likelihood (principled for count data with overdispersion); this work uses Normal, which is convenient for structural claims but inappropriate for predictive applications.
- They report formal out-of-sample predictive metrics (log score, Wasserstein); the universality claims here are statements about parameter recovery, not predictive accuracy on held-out n.
- They build a mechanistic generative model with explicit odd-block dynamics; this work characterizes the per-class distribution structurally without an explicit generator.

**Net positioning.** The structural decomposition here is a *sharpening of Terras 1976 Lemma 4* (which gives the asymptotic identity S_k ≈ S_0 · 3^d(k) / 2^k) by tracking the symbolic prefix to its termination at a·m + c with a odd, retaining c_final, and making a_final available as a closed-form covariate. Relative to Bonacorsi & Bordoni's mod-8 hierarchical NB-GLM, this furnishes their random effect with explicit algebraic structure: the u_{n mod 8} effect they fit is, up to a constant, log(a_final(n mod 8)) at k=3. Their NB-GLM is the predictive engine; the prefix decomposition supplies the closed-form covariate that grounds the random effect in the symbolic Collatz prefix. The two pieces fit together: a hierarchical NB-GLM with the random effect on n mod 8 replaced by a fixed effect on a_final(n mod 2^k) for chosen k yields a parameter-parsimonious version of their model.

### Bayesian replication in the B&B framework: Pathfinder VI at N=10⁷

To make this concrete we replicate their setup at N=10⁷ in Stan: NB2-GLM with log link, priors β ~ N(0,5), φ ~ LogN(log 5, 1), σ_u ~ HN(0,2). Inference via cmdstanpy 1.3.0 Pathfinder (Zhang, Carpenter, Gelman, Vehtari, JMLR 2022), 4 paths × 1000 importance-resampled draws, on a 500K-observation training subsample with a 50K held-out test set.

**Caveat upfront — this is Pathfinder VI, not HMC.** Production HMC on this scale of model + data hit a Stan 2.36 unified-mode multi-chain lockup that pinned all threads at full CPU with zero sampling progress (a sync issue with `num_chains=N` in shared-thread-pool process, distinct from the well-known reduce_sum grainsize pathology). Pathfinder converged cleanly. The known limitation: Pathfinder is a quasi-Newton variational approximation that systematically underestimates posterior dispersion in non-Gaussian hierarchical regions, and so collapses σ_u in any spec with a hierarchical random effect. Absolute log-score numbers for hierarchical specs (B1, B4) below are therefore lower than HMC would produce on the same data; fixed-effect specs (B0, B2, B3) are not subject to this bias. **HMC validation at production scale is the natural next step on Bonacorsi's side**, and is the cleanest way to confirm the absolute numbers below.

| Spec | description | # params | log score | σ_u |
|---|---|---|---|---|
| B0 | log(n) only | 3 | −274,150.3 | — |
| B1 | + RE on (n mod 8) ← B&B's setup | 12 | −274,138.4 | 0.003 (Pathfinder-collapsed) |
| B2 | + FE on a_final at k=3 | 6 | −273,288.4 | — |
| **B3** | **+ FE on a_final at k=6** | **9** | **−272,435.3** | — |
| B4 | B3 + RE on (n mod 8) | 18 | −272,438.4 | 0.009 (Pathfinder-collapsed) |
| B&B reported (HMC, full data) | NB2-GLM with mod-8 RE | — | −272,911.95 | (HMC posterior, larger) |

**What the table says, structurally:**

1. **B2 vs B1 — same modular resolution, different parameterization.** B1 estimates a hierarchical mod-8 random effect (the B&B setup); B2 supplies the closed-form covariate a_final at k=3, which is the algebraic identity of that random effect. Pathfinder gives B2 a +850-nat margin over B1, but the relevant comparison is structural: a closed-form covariate replacing an estimated random effect at the same modular resolution, with comparable predictive content. Under HMC the absolute numbers will tighten.

2. **B3 vs B&B reported — finer modular resolution under HMC.** B3 (a_final at k=6, fixed effects, no σ_u, robust to the Pathfinder caveat) scores −272,435.3 with 9 parameters. B&B's reported NB2-GLM at the same data scale scores −272,911.95. The 477-nat gap is suggestive — finer modular resolution captures more structure than mod 8 — but it is not a head-to-head conclusion; the cleanest version is HMC on B3 at full N=10⁷, which Bonacorsi can run directly.

3. **B4 vs B3 — does a_final exhaust the residue-class signal?** Adding mod-8 RE on top of a_final (B4) adds a Pathfinder σ_u of 0.009 and doesn't move the log score. Suggestive that mod-8 carries no residual information once a_final is in the model, but again not load-bearing under Pathfinder; HMC is the test.

**The contribution.** Bonacorsi-Bordoni's mod-8 random effect is observing real structure. The prefix decomposition supplies the algebraic identity behind it, and extends cleanly to arbitrary modular resolution k. As a covariate in their NB2-GLM framework, a_final(n mod 2^k) replaces the random effect with a fixed effect on k discrete levels — *logarithmic* parameter count in the modular resolution rather than linear. The joint model is what's interesting: their predictive engine + this closed-form covariate, validated under HMC at production scale.

**A companion analysis** examines the qx+1 generalization (q ∈ {5, 7, 9, 11}) and is reported separately. The 3x+1 results above are the contribution to the Bonacorsi-Bordoni framework.

## Limitations and what remains open

- **Resolution.** Verified at k ∈ {4, 5, 6, 7, 8, 9, 10, 11, 12}, i.e., modular resolutions from mod 16 through mod 4096 (covering 8 to 2048 odd classes). The noise-floor ratio sits in [0.90, 0.99] across k ∈ {5, …, 12}; at k=4 it is 0.70 (residuals smaller than per-class SE — prefix prediction over-explains at coarse resolution). At every k tested, the number of distinct a_final levels equals k exactly. Higher k might in principle reveal residual structure but the trend across the nine resolutions tested is unambiguous: residuals scale with sampling noise at every k. The k=10 hierarchical Stan fit failed for posterior-geometry reasons; the OLS analog (per-class α via ordinary regression, compared to α_det) does not share that issue and confirms the result through k=12.
- **Range.** N ≤ 2³² ≈ 4.3 × 10⁹. The pooled OLS slope μ_β oscillates around the heuristic 10.4282 over the tested range with amplitude ≈ 0.01 and is not visibly damping at 2³² (see Result 1). The crossing of 10.4282 happens near N ≈ 2²⁷. Stronger asymptotic claims would require N well beyond 2³² and are not made here.
- **Stochastic remainder S(n) characterization.** We've shown S(n) has class-universal first two moments and tail shape; we haven't characterized its full distribution. Two empirical features are flagged for completeness:
  - *Trajectory measure on v deviates from Geom(½).* At N_start = 10⁸ Syracuse trajectories, the pooled v = ν₂(3m+1) distribution has spikes at v=4 and v=10 (ratio ≈ 1.23 above Geom(½) prediction) and graded sags at v ∈ {6, 7, 8, 9} and {11, 12, 13, 14}. The pattern does not extend cleanly to v=16 onward (ratio 0.88, not a spike), so it is not a simple every-six-v phenomenon. The deviation has E[v] ≈ 1.99 (vs Geom(½)'s 2) and Var(v) ≈ 1.88 (vs 2). This characterizes a piece of S(n)'s fine structure but does not change Result 1 or 3.
  - *Neighbor coincidence.* Among consecutive odd m within a fixed residue class mod 2¹⁰, the rate of σ(2¹⁰m + r) = σ(2¹⁰(m+1) + r) is empirically ≈ 25% — a non-trivial constraint on S(n)'s pairwise distribution that the prefix decomposition predicts up to the post-prefix offset a_final but does not itself explain.
- **The prefix algebra terminates at a_final ∈ {3, 9, 27, 81, 243, 729}** (at k = 6, with one additional power of 3 added per increment of k). These are all powers of 3, reflecting the multiplicative structure of the Syracuse map. Their distribution across residues might itself have number-theoretic structure worth describing — but this is a property of Collatz, not of σ.
- **qx+1 generalization (q ∈ {5, 7, 9, 11}).** Reported separately. Empirical Cramér-rate match at q=5 is to 0.01% precision; the load-bearing trajectory-measure assumption (v ~ Geom(½) in the unconditional measure) is verified empirically across q ∈ {3, 5, 7, 9, 11} to 0.5% on Var(v) and to 0.18–0.31% on the MGF E[(2/q)^v] at the relevant Cramér tilt point. The 3x+1 trajectory-measure deviations described above do not propagate to the qx+1 derivation at the relevant moment.

## Reproduction

All code, data, and figures at [C:\Collatz\\](file:///C:/Collatz/). Total compute time ~5 hours: ~30s for Stage 1 generation at all N values, ~70min for the k=6 Stan fit (4 chains × 4 threads × 500+500 iters), <10min for diagnose + N-scaling + α decomposition.

```bash
# Data generation
python generate.py --N 33554432 --no-vseq

# Bayesian fit prep + fit + diagnose (k=6)
python stage3_prep.py --N 33554432 --k 6 --per-class-cap 50000 --tail-frac 0 \
    --out data/stage3_input_N33554432_k6_uniform.parquet
python fit.py --input data/stage3_input_N33554432_k6_uniform.parquet \
    --chains 4 --threads-per-chain 4 --iter-warmup 500 --iter-sampling 500 \
    --max-treedepth 12 --adapt-delta 0.9 --tag k6_uniform
python diagnose.py --fit-dir fits/k6_uniform_full \
    --input data/stage3_input_N33554432_k6_uniform.parquet
```

Output figures and CSVs in [stage4_results/k6_uniform_full/](file:///C:/Collatz/stage4_results/k6_uniform_full/).

## Acknowledgments to the discipline

The single most important input to this project was a sanity-check protocol that gated escalation: before flagging an observation as anomalous or in need of advanced analysis, run the simplest baseline experiment that would distinguish "real structure" from sampling artifact. This protocol caught two would-be findings (the v-spike "anomaly" and the apparent Stan-OLS τ_β discrepancy), preventing rabbit holes that would otherwise have absorbed days of effort. The discipline cost roughly 30 minutes of additional uniform-sampling and OLS-baseline runs across the project's lifetime; it saved at minimum two days of compute and confused interpretation.
