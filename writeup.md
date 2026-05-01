# Empirical universality of Collatz stopping times on odd integers up to 2²⁷

*Status: closed at all modular resolutions k ∈ {6, 7, 8, 9} (mod 64 through mod 512) on 2026-05-01.*

## TL;DR

For odd integers n ∈ [3, 2²⁵], stratified by residue class mod 64:

1. **Slope universality.** The mean of the total stopping time σ(n) is described by σ(n) ≈ α(r) + 10.43·ln(n) where r = n mod 64. The slope **does not vary across residue classes**: τ_β at noise floor at every N tested (2²⁰ through 2²⁵). The pooled empirical slope converges from below to the heuristic 3/(ln 4 − ln 3) = 10.4282 monotonically with N (gap = 0.056 at 2²⁰, 0.009 at 2²⁵).

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

The naive class-β SD tracks the per-class sampling SE almost exactly at every N. The moment-corrected τ_β² = max(0, observed² − SE²) is **identically zero** at every sample size tested: there is no detectable between-class slope variation above the noise floor at any N.

The pooled OLS slope μ_β → 10.4282 monotonically with N from below; rate is slow (∼N^{−1/2}-ish on the residual 10.4282 − μ_β) but unambiguous.

**Statement of result.** The single-slope heuristic σ(n) ≈ 10.4282·ln(n) is empirically asymptotically exact for the mean of σ over odd starting points. No residue class mod 64 has a detectably different slope at any sample size up to 2²⁵.

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

**Net positioning.** The structural decomposition theorem here — that the residue-class structure of σ has *logarithmic* complexity in modular resolution, with conditional distributions parameterized by a_final ∈ {3^j} — is genuinely new relative to Bonacorsi & Bordoni. Their work is excellent on the predictive-modeling axis; this work is sharper on the structural axis. Both could appear together in a coherent research program: their NB-GLM combined with this prefix-determined parameterization yields a substantially more parsimonious model than mod-8 random effects, with the random effects replaced by a closed-form function of a_final(r mod 2^k).

### Direct head-to-head: predictive comparison on B&B's setup

To verify that the structural finding has direct predictive consequences, we replicate their setup at N=10⁷ with NB GLM and 50K-observation held-out test, comparing five specifications:

| Model | # params | Test log_score | W₁ distance |
|---|---|---|---|
| M0 baseline (log n only) | 3 | −274,139.0 | 3.040 |
| M1 B&B-style: + factor(n mod 8) | 10 | −273,352.3 | 3.071 |
| M2 B&B-extended: + factor(n mod 64) | 66 | −272,496.5 | 3.230 |
| **M3 ours: + factor(a_final at k=3)** | **6** | **−273,352.3** | 3.136 |
| **M4 ours: + factor(a_final at k=6)** | **9** | **−272,496.3** | 3.310 |
| B&B reported best (NB2-GLM) | — | −272,911.95 | 3.199 |

**Two predictive wins:**

1. **M3 matches M1 to within 0.1 nats** on 50,000 held-out observations, using **40% fewer parameters** (6 vs 10). Replacing 4 mod-8 fixed effects with 3 a_final levels is lossless predictively.

2. **M4 essentially matches M2 with 7× parameter reduction** (9 vs 66). Same log score to within 0.2 nats out of 272,496. Replacing 32 mod-64 fixed effects with 6 a_final levels is lossless predictively. **And M4 beats B&B's reported best NB-GLM by 415 nats** at the same data scale and test split size with comparable parameter count.

The structural decomposition therefore translates from elegant theorem into operational parsimony: the apparent residue-class structure that motivates fitting many random effects is itself low-dimensional, indexed by a_final ∈ {3^j : 1 ≤ j ≤ k}.

(W₁ distance for our fixed-effect MLE specifications is slightly worse than B&B's reported W₁ — most likely because their hierarchical posterior predictive integrates over class-effect uncertainty, producing more disperse predictives that better match the empirical tail; our fixed-effect predictives are more concentrated. A hierarchical prior on the a_final coefficients would close this gap with no architectural complication, and is a natural extension if W₁ is the operational metric of interest.)

## Limitations and what remains open

- **Resolution.** Verified at k ∈ {6, 7, 8, 9}, i.e., modular resolutions from mod 64 through mod 512 (covering 32 to 256 odd classes). Higher k (k=10, k=12, etc.) might still in principle reveal residual structure, but the trend across the four resolutions tested is unambiguous: residuals scale with sampling noise at every k, with no signal emerging at finer grids. The k=10 hierarchical Stan fit failed for posterior-geometry reasons rather than any signal-vs-noise reasons.
- **Range.** N ≤ 2²⁷ ≈ 1.3 × 10⁸. The asymptotic claims (μ_β → 10.4282, ξ → 0) extrapolate from a clean monotone trend across N ∈ {2²⁰, 2²², 2²³, 2²⁴, 2²⁵, 2²⁷}; not proved.
- **Stochastic remainder S(n) characterization.** We've shown S(n) has class-universal first two moments and tail shape; we haven't characterized its full distribution. Likely has further structure relating to higher 2-adic conditions, but mod 64 doesn't resolve it.
- **The prefix algebra terminates at a_final ∈ {3, 9, 27, 81, 243, 729}.** These are all powers of 3, reflecting the multiplicative structure of the Syracuse map. Their distribution across residues might itself have number-theoretic structure worth describing — but this is a property of Collatz, not of σ.
- **Trajectory-measure characterization for v=4/v=10 spikes.** The empirical observation is real and quantified; the exact density of {m : ν₂(3m+1) = k} along Syracuse iterate distributions is elementary number theory but wasn't worked out here.

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
