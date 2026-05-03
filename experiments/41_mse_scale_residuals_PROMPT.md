# Experiment 41 — Prompt (saved 2026-05-02)

**Task: Test whether MSE of σ predictions against α_det + K_h · log(n) reveals the K_h − K_eff = 0.486 mechanism**

## Background

The bridge equation `s_mean(r; f) ≈ α_det(r) + K_h · log(N/f) + ε(σ) + (K_h − K_eff) · log(f)` has K_h − K_eff = 0.486 as a structural slope on log(threshold) that's empirically stable (R² = 0.989 across observables) but mechanism-open. Two prior candidates were falsified: trajectory-weighted σ-vs-log-v slope (sub-band instability, slope varies −12 to +38 depending on aggregation) and trajectory-measure correction E[v] = 1.995 (wrong direction, predicts K up not down).

The hypothesis: K_h − K_eff lives in the variance/MSE structure of σ predictions, not in mean-slope behavior. If σ residuals (σ_actual − [α_det + K_h · log(n)]) have scale-dependent variance that grows differently above vs below √N, that scale-dependence might generate the observed slope on log(threshold) through how MSE-driven corrections accumulate during descent.

## Method

1. For each odd n in a sample (use existing parquet at N = 2²⁵ or 2²⁷, sample ≥ 200K orbits), compute the residual:
   ```
   resid(n) = σ(n) − α_det(r(n)) − K_h · log(n)
   ```
   where r(n) = n mod 2^k for k = 6 or 8, and α_det is the closed-form predictor.

2. Walk each orbit, recording (v, resid_partial(v)) pairs where:
   - v is the orbit value at each step
   - `resid_partial(v) = steps_remaining(v) − K_h · log(v)`
   This is the analog of the residual at the orbit value v rather than the starting point.

3. **The MSE-by-scale measurement:**
   - Bin by log(v) into 20 bins from log(1) to log(N).
   - In each bin, compute mean and variance of resid_partial(v).
   - Plot mean(resid_partial) and var(resid_partial) vs log(v).

## What to look for

- Does var(resid_partial) have scale-dependent structure? Specifically, is variance higher at small v (below √N) than at large v?
- Does mean(resid_partial) show a piecewise-linear structure with a kink at log(v) ≈ log(√N)? If yes, the slope above the kink should match K_h asymptote and the slope below should match K_eff.
- Does the integral of var(resid_partial) along the descent recover the empirical slope 0.486 on log(threshold)?

## Specific quantitative test

If MSE-driven corrections are the mechanism:
```
K_h − K_eff = ∫ d(var(resid))/d(log v) · weight(v) dv / ∫ weight(v) dv
```
where weight(v) is the trajectory-visit measure. Compute this integral numerically and compare to 0.486.

## Output

- Plot 1: mean(resid_partial) vs log(v), 20 bins, with K_h slope reference line. Look for kink at log(√N).
- Plot 2: var(resid_partial) vs log(v), 20 bins. Look for scale-dependent structure.
- Plot 3: cumulative residual variance during descent (running sum) vs log(v).
- Numerical: integral of d(var)/d(log v) weighted by trajectory measure, compared to 0.486.

## Decisive outcomes

- If mean(resid_partial) shows a clean kink at √N with slopes K_h above and K_eff below: two-regime descent is real, MSE not needed.
- If var(resid_partial) has scale-dependent structure that integrates to 0.486 against trajectory weights: MSE-by-scale is the mechanism.
- If neither: K_h − K_eff lives in something else (correlation structure across scales, finite-N boundary effects, or pure stopping-time-distribution shape).

Standard procedure: compute, plot, report. If the integral matches 0.486 to within 5%, MSE-by-scale is identified. If it matches within 20%, partial mechanism. If it doesn't match, ruled out.
