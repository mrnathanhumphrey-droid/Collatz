# Tao Bridge Tightening — Consolidated Findings

**Tasks:** TA.1 (N-stability), TA.2 (trim-quantile sweep), TA.3 (parametric form of gap).
**Date:** 2026-05-02
**Data:** existing parquets at N ∈ {2²⁵, 2²⁷}; sigma cache rebuilt at N ∈ {2²⁸, 2³⁰, 2³²}.

The base bridge claim: `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N)) + ε`, with `K_h = 3/log(4/3)`, slope at K_h = 1.000 ± 0.005 across 40 verification cells, ε an observable-dependent structural correction.

The three tasks tighten the characterization of ε.

---

## TA.1 — N-stability of the σ structural offset

**Run.** Sigma cache built and per-class mean σ computed at k = 8, 10, 12 for
each N ∈ {2²⁵, 2²⁷, 2²⁸, 2³⁰, 2³²}. `gap_at_mean = σ̄ − ⟨α_det⟩ − K_h · ⟨log N⟩`.

| N       | log N   | gap (k=8) | gap (k=10) | gap (k=12) | per-class SE (k=8) |
|---------|---------|-----------|------------|------------|---------------------|
| 2²⁵     | 16.329  | −2.4468   | −2.4468    | −2.4468    | 0.174               |
| 2²⁷     | 17.715  | −2.4514   | −2.4514    | −2.4514    | 0.091               |
| 2²⁸     | 18.408  | −2.4492   | −2.4492    | −2.4492    | 0.066               |
| 2³⁰     | 19.794  | −2.4526   | −2.4526    | −2.4526    | 0.034               |
| 2³²     | 21.181  | −2.4574   | −2.4574    | −2.4574    | 0.018               |

**Variation across k at fixed N: identically zero.** All k=8/10/12 gaps
match to 4 decimal places at every N. The gap is a property of the σ
distribution, not the modular grid.

**Variation across N at fixed k: 0.0105 across 7 doublings of N.**
Roughly −0.0022 per unit log N (not strictly monotone — small reversal at
2²⁸). The drift magnitude is ~5× the per-class SE at the largest N, so
detectable but small. Tao's O((log N)^0.6) bound at N=2³² is ~6.4; the
observed drift sits well within it.

**Verdict.** The σ gap is a structural near-constant: ≈ −2.45 across the
entire tested range, with tiny finite-N drift consistent with sub-leading
Tao corrections. Closed-form derivation candidate (constant ε for σ) is
empirically supported.

**Closed-form decomposition (analytical check).** ⟨α_det⟩ at every k equals
`E[ℓ] − K_h · ⟨log(2^k/a_final)⟩`. Computed exactly:

| k | E[prefix_steps] | K_h · ⟨descent⟩ | ⟨α_det⟩ = diff |
|---|----------------|-----------------|-----------------|
| 6  | 9.50  | 3.27 | +6.23 |
| 8  | 12.50 | 6.27 | +6.23 |
| 10 | 15.50 | 9.27 | +6.23 |
| 12 | 18.50 | 12.27 | +6.23 |

⟨α_det⟩ = +6.23 **exactly across all k**, an invariant of the prefix
algebra. The Tao residual at the global mean
`σ̄ − K_h · log N = ⟨α_det⟩ + ε = +6.23 − 2.45 = +3.78`. This decomposes
into the prefix-algebra invariant (+6.23) plus a post-prefix descent
correction (−2.45). The −2.45 says the actual post-prefix descent takes
2.45 fewer Collatz steps than the K_h · log(post-prefix value)
random-walk heuristic predicts.

---

## TA.2 — Trim-quantile sweep, find q* where gap = 0 at √N

**Run.** N = 2²⁷, k = 8, observable = s @ √N. Trim quantile q ∈
{0.05%, 0.1%, 0.2%, 0.3%, 0.5%, 0.75%, 1%, 1.5%, 2%, 3%, 5%}. Within each
class, drop the top-q values and compute mean of the rest.

| q (%)   | offset gap from K_h · log(N)/2 |
|---------|--------------------------------|
| 0       | +2.16 (raw mean)               |
| 0.05    | +2.01                          |
| 0.10    | +1.88                          |
| 0.20    | +1.66                          |
| 0.30    | +1.46                          |
| 0.50    | +1.08                          |
| 0.75    | +0.65                          |
| 1.00    | +0.26                          |
| 1.50    | −0.46                          |
| 2.00    | −1.14                          |
| 3.00    | −2.35                          |
| 5.00    | −4.47                          |

Linear interpolation: **q* ≈ 1.18%** drives gap to 0 (between 1% and 1.5%).

**Comparison to log^(−c) N.** At N = 2²⁷, log N = 17.715:

| c    | log^(−c) N | q* / log^(−c) N |
|------|------------|------------------|
| 0.4  | 0.317      | 0.037           |
| 0.5  | 0.238      | 0.050           |
| 0.6  | 0.178      | 0.066           |
| 0.7  | 0.134      | 0.088           |
| 0.8  | 0.100      | 0.118           |
| 1.0  | 0.0564     | 0.209           |
| 1.5  | 0.0134     | 0.879           |
| 2.0  | 0.00319    | 3.698           |

c such that log^(−c) N = q*: **c = 1.55**. Not a clean Tao exceptional-set
exponent (Tao's main bounds use exponents like 0.6 from O((log N)^0.6)).

**Verdict.** The trim that closes the σ-fiber gap to zero (at √N, k=8,
N=2²⁷) is 1.18%. This does NOT correspond to a recognizable Tao
exceptional-set density at this N. The trim-zero-gap interpretation is
weaker than hoped: the trim fraction looks like a property of the σ
distribution shape (right-tail mass) at this N rather than a structural
Tao quantity. Worth re-running at N = 2³² to test N-dependence; not
pursued in this round.

---

## TA.3 — Parametric form of gap vs threshold / Δlog

**Run.** 13 (observable, N) cells: 5 observables at N ∈ {2²⁵, 2²⁷}, plus
σ at 2²⁸, 2³⁰, 2³². Fit candidate parametric forms via least squares.

| Model                                          | SSE    | R²    |
|------------------------------------------------|--------|-------|
| `gap = a + b · Δlog`                           | 3.28   | 0.955 |
| `gap = a + b · log(Δlog)`                      | 4.23   | 0.942 |
| `gap = a + b / Δlog`                           | 9.37   | 0.872 |
| `gap = a + b · Δlog + c · Δlog²`               | 2.89   | 0.961 |
| **`gap = a + b · log(N) + c · Δlog`**          | **0.74** | **0.990** |
| `gap = a + b · Δlog + c · log(log N)`          | 0.76   | 0.990 |
| `gap = a + b · log(log N)`                     | 53.5   | 0.269 |
| **`gap = a + b · log(threshold)`**             | **0.80** | **0.989** |
| `gap = a + b · Δlog + c · log(Δlog)`           | 3.17   | 0.957 |
| `gap = a + b/Δlog + c · Δlog`                  | 3.24   | 0.956 |

**Best two-parameter form:** `gap ≈ −2.35 + 0.486 · log(threshold)`.
log(threshold) = log(f(N)) where f is the first-passage threshold (or 1 for σ).

**Slope ≈ 0.486 ≈ 1/2.** The "non-monotone in Δlog" appearance was an
artifact of using Δlog instead of log(threshold) as the axis; the
relationship is cleanly linear in log(threshold).

**Predicted gaps under `gap = −2.35 + 0.486 · log(f)`:**

| observable | log(f) | predicted | observed |
|------------|--------|-----------|----------|
| σ                | 0     | −2.35 | −2.45 (off 0.10) |
| s @ √N/log N     | 5.98  | +0.56 | +1.15 (off 0.59) |
| s @ √N           | 8.86  | +1.96 | +2.16 (off 0.20) |
| s @ √N · log N   | 11.73 | +3.36 | +3.03 (off 0.33) |
| s @ N^(2/3)      | 11.81 | +3.39 | +3.01 (off 0.38) |

Fit is good but not perfect. The residual structure may carry the small
finite-N corrections we already see in TA.1's σ-only N-sweep.

**Pure log(log N) fails (R²=0.27).** The gap is not driven by iterated-log
alone; the threshold-depth axis is dominant.

**Verdict.** Gap structure is dominantly linear in `log(threshold)` with
slope ≈ 1/2 and intercept ≈ −2.35. The intercept matches the σ structural
constant (TA.1). No clean Tao-style log-log signature.

---

## Synthesis

1. The bridge is structurally clean at the leading order: `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N))` with slope = 1.000 ± 0.005 at K = K_h, no fit recalibration.

2. The structural correction ε(observable, N) decomposes as:

   `ε ≈ −2.35 + (1/2) · log(f(N))`

   The first term (−2.35) is a near-constant in N (drifts by 0.01 across N=2²⁵ → 2³², well within Tao's O((log N)^0.6) bound). The second term (1/2 · log f) is observable-dependent.

3. At the global mean across classes:
   `σ̄ − K_h · log N = ⟨α_det⟩ + ε(σ) = +6.23 − 2.45 = +3.78`.
   The +6.23 is `E[prefix_steps] − K_h · ⟨descent during prefix⟩`, an exact invariant of the prefix algebra (k-independent). The −2.45 is the post-prefix descent correction.

4. Trim-1.18% mean drives the σ-fiber gap (at √N) to zero. The trim fraction does not match a clean log^(−c) N exponent (implied c = 1.55), so the trim-zero-gap interpretation is closer to "removes the σ right-tail" than to "Tao exceptional-set density." Worth retesting at larger N.

5. **What's now exposed for the bridge:** the −2.35 constant (across all observables) and the +0.486 ≈ 1/2 slope on log(threshold). Both deserve closed-form derivation. The 1/2 slope hints at a √f scaling — the gap grows like log(√f) — but the underlying mechanism is not resolved here.

6. **What's not the answer (negative results):** pure log(log N), 1/Δlog, Δlog², and quadratic Δlog forms all fit worse than the log(threshold) form. The user's hypothesis "gap = −E[ℓ] + heuristic prefix" doesn't match numerically (E[ℓ] − heuristic = +6.23, not −2.45) — but the framing is right: the bridge gap decomposes into a prefix-algebra invariant plus a post-prefix correction, with the −2.45 being the latter only.

## Files

- `experiments/36_TA1_sigma_offset_N_sweep.py` — N=2²⁵..2³² run
- `experiments/37_TA2_trim_quantile_sweep.py` — trim sweep at N=2²⁷
- `experiments/38_TA3_parametric_fit.py` — parametric fits on aggregated data
- `experiments_output/36_TA1_sigma_offset_N_sweep.csv`
- `experiments_output/37_TA2_trim_quantile_sweep_N134217728_k8.csv`
