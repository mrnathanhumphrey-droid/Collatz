# Result 38: Log-trajectory shape characterization — Brownian heuristic FALSIFIED for log_excursion (Gamma not Exp); peak occurs early; descent depth nearly constant

**Date:** 2026-05-03. Per-orbit log-trajectory observables characterized by σ-quantile band at N=2^32 (5 seeds × 1M orbits, 1.8s compute).

Tests four observables:
- **log_excursion** = log(n_max) - log(n_start) [peak above start]
- **log_descent** = log(n_start) - log(n_min) [depth below start]
- **peak_fraction** = t_peak / σ_S [where in orbit lifetime peak occurs]
- **cross-correlations** between these and σ_S

**Verdict per brief decisive outcomes:**
- (a) Brownian heuristic prediction (log_excursion ~ Exp): **FALSIFIED** at ΔAIC = +10K to +21K
- (b) Different parametric form: **CONFIRMED**. log_excursion is Gamma with shape ~1.6-1.9 across all bands
- (c) NA — clean parametric form (gamma) found
- (d) peak_fraction structure: **CONFIRMED**. Peak occurs early (median 3-7% of σ_S) with band-dependent timing

---

## 1. Per-band log_excursion summary

| q-band | mean | sd | p50 | p95 | max |
|---|---|---|---|---|---|
| 0.125 | 0.496 | 0.630 | 0.41 | 1.74 | 7.25 |
| 0.375 | 0.699 | 0.815 | 0.41 | 2.38 | 8.36 |
| 0.625 | 0.836 | 0.944 | 0.52 | 2.75 | 10.27 |
| 0.875 | 1.121 | 1.165 | 0.81 | 3.44 | 12.56 |
| 0.975 | 1.749 | 1.602 | 1.33 | 4.91 | 15.67 |

**log_excursion mean grows monotonically with q-band**: 0.50 → 1.75 nats from bottom to top σ-quartile. Top-σ orbits (slow descent) have ~3.5× the excursion magnitude of bottom-σ orbits.

**Median p50 always smaller than mean** → right-skewed distributions. Max excursion reaches 15.67 nats (4.7M× start magnitude) for top-band orbits.

## 2. Brownian-with-drift heuristic FALSIFIED

Heuristic prediction: for a random walk with negative drift μ and step variance σ², the maximum excursion above start has Exponential distribution with rate 2|μ|/σ². If true, log_excursion | band ~ Exp.

**Empirical:** log_excursion is Gamma with shape ~1.5-1.9, not Exp:

| q-band | best fit | gamma_a | exp_mean | ΔAIC_exp | ΔAIC_lognorm |
|---|---|---|---|---|---|
| 0.125 | gamma | 1.92 | 0.82 | **+20,852** | +2,650 |
| 0.375 | gamma | 1.73 | 1.00 | +15,457 | +4,089 |
| 0.625 | gamma | 1.65 | 1.14 | +12,983 | +3,563 |
| 0.875 | gamma | 1.59 | 1.40 | +11,391 | +6,047 |
| 0.975 | gamma | 1.56 | 1.98 | +10,617 | +10,940 |

Exp is REJECTED at ΔAIC > 10,000 across all bands. The empirical Gamma shape parameter varies smoothly from 1.92 (bottom) to 1.56 (top). For Gamma, shape=1 ⇒ Exponential. Empirical shapes are well above 1, decisively non-exponential.

**Implication:** Brownian motion with negative drift is NOT a good heuristic for the Collatz log-trajectory's peak excursion. The trajectory has structural deviations from Brownian — specifically, longer/heavier tails relative to Exp at the peak observable.

## 3. log_descent is nearly constant; tracks log(N/m_j)

| q-band | mean | sd | p50 | min | max |
|---|---|---|---|---|---|
| 0.125 | 18.98 | 1.61 | 19.57 | 3.06 | 20.57 |
| 0.375 | 19.32 | 1.30 | 19.78 | 7.04 | 20.57 |
| 0.625 | 19.49 | 1.13 | 19.86 | 9.03 | 20.57 |
| 0.875 | 19.62 | 0.99 | 19.94 | 9.37 | 20.57 |
| 0.975 | 19.74 | 0.87 | 20.02 | 10.31 | 20.57 |

**log_descent ≈ log(n_start) - log(m_j_attractor)**, since most orbits have n_min = m_j (the attractor).
- log(N=2^32) - log(m_2=5) = 22.18 - 1.61 = 20.57 — matches max p value ✓
- Mean log_descent ≈ 19-20 nats: ⟨log m_start | absorbing at j=2⟩ ≈ 21.18, log(m_2)=1.61, descent = 19.57 ✓

Best fit is **gamma with very large shape** (123-477) — essentially Gaussian-like (shape → ∞ ⇒ Gaussian). The distribution is narrow because n_min is dominantly m_j (small set of values per absorbing class).

This is structurally trivial — log_descent is determined by attractor selection (Result 30/34) and starting m, not by orbit shape.

## 4. peak_fraction: peaks occur EARLY in orbit lifetime

| q-band | mean | sd | p50 (median) | min | max |
|---|---|---|---|---|---|
| 0.125 | 0.0568 | 0.0807 | **0.0256** | 0 | 0.82 |
| 0.375 | 0.0630 | 0.0872 | 0.0308 | 0 | 0.82 |
| 0.625 | 0.0638 | 0.0883 | 0.0270 | 0 | 0.87 |
| 0.875 | 0.0756 | 0.0965 | 0.0388 | 0 | 0.88 |
| 0.975 | 0.1090 | 0.1205 | **0.0662** | 0 | 0.80 |

**Median peak occurs within first 3-7% of orbit lifetime.** Typical Collatz orbit:
- Rises briefly to peak at ~5% of σ_S (median)
- Then descends throughout the remaining ~95% of orbit
- Peak excursion magnitude is small (median log ~0.4-1.3 nats)

**This is decisively asymmetric.** A symmetric trajectory (rise then equal-magnitude fall) would have peak_fraction ≈ 0.5. Empirical median is 0.026-0.066 — peak is far in the early phase.

Top-σ orbits peak slightly later (median 6.6% vs bottom-σ 2.6%) and have somewhat larger excursions, but all bands show "rise quickly, then descend" structure.

## 5. Cross-correlations

| q-band | ρ(exc, desc) | ρ(exc, σ) | ρ(peak, σ) | **ρ(peak, exc)** |
|---|---|---|---|---|
| 0.125 | -0.04 | +0.16 | +0.06 | **+0.71** |
| 0.375 | -0.06 | +0.04 | -0.00 | +0.68 |
| 0.625 | -0.06 | +0.05 | +0.01 | +0.67 |
| 0.875 | -0.07 | +0.11 | +0.07 | +0.64 |
| 0.975 | -0.08 | +0.22 | +0.17 | +0.58 |

**Excursion and descent are nearly INDEPENDENT** within each band: ρ(exc, desc) ≈ -0.05. Peak excursion height doesn't predict depth of descent. This is consistent with descent depth being "fixed" by attractor (mostly determined by m_j) while excursion is independent variation.

**Peak timing strongly correlates with excursion magnitude**: ρ(peak, exc) = 0.58-0.71. Bigger excursions take longer to reach the peak (relatively). Makes sense: bigger excursion requires more accumulated upward steps before turn-around.

Excursion / σ correlations are weak (0.04-0.22) — orbit length and excursion height are nearly independent within band (after band-conditioning).

## 6. Implications for the framework

**(1) Brownian heuristic falsified for log_excursion:** the trajectory measure deviates from Brownian-motion-with-drift in the peak-excursion observable. This is a fourth structural slice of the trajectory measure (alongside w_q(q), P(q|j), ⟨v|q,j⟩-residual from Result 34).

**(2) Asymmetric trajectory shape:** peaks early, descends most of the time. The "rise" portion is small in both magnitude (excursion <2 nats) and duration (~5% of σ_S), while the "fall" portion dominates. This is consistent with the negative drift but quantifies the asymmetry precisely.

**(3) Independence of excursion and descent:** allows separate parametric characterization. log_descent reduces to attractor-selection (Result 30/34); log_excursion has its own gamma-shaped distribution per band.

**(4) Per-band gamma shape parameter** as a new structural quantity: shape α_q decreases from 1.92 (bottom-q) to 1.56 (top-q). This is itself a function of q-band — analogous to w_q(q) Esscher tilt parameter. Closed form for α_q(q) is the next reducible-to-trajectory-measure piece.

## 7. Verdict

Per brief outcomes:
- **(a) log_excursion ~ Exp:** FALSIFIED at ΔAIC > 10K
- **(b) Different parametric form:** CONFIRMED — Gamma with shape varying smoothly across q-bands
- **(c) No clean form:** NA — gamma is clean
- **(d) peak_fraction structure:** CONFIRMED — peaks early (~5%), with band-dependent timing

**For v3.6 framing:** adds three new per-band structural quantities:
- log_excursion gamma shape α_q(q) and scale σ_q(q)
- peak_fraction distribution per band
- cross-correlation ρ(peak, exc) ≈ 0.6-0.7 (band-stable)

These are **additional empirical slices of the trajectory measure** (not new Lagarias-class observables; they reduce to it). Brownian heuristic is decisively wrong as a leading-order model for log-trajectory shape — informative constraint on the trajectory measure structure.

## 8. Files

- `log_trajectory_shape.py` — walker + per-band parametric fits + correlations
- `log_trajectory_shape_test.md` — this document (Result 38)
- `experiments_output/log_trajectory_shape.csv` — per-band fit parameters
- `experiments_output/log_trajectory_shape_log.txt` — full output

Total compute: 1.8s (5M orbits at N=2^32, numba-parallel walker with trajectory tracking).
