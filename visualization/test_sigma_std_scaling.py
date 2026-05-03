"""
test_sigma_std_scaling.py - Test K_eff_band(q) = K_h + b·z_q closed-form candidate.

Hypothesis: σ at fixed log(n) is approximately Gaussian with std φ(log n) growing
linearly in log(n): φ(log n) = a + b·log(n).

If true:
  X_q(log n) = K_h·log(n) + c + φ(log n)·z_q
  d X_q / d log(n) = K_h + b·z_q
  K_eff_band(q) = K_h + b·z_q

Empirical K_eff_band slope ≈ 2.7 per quartile, z-quartile shift ≈ 0.78,
so predicted b ≈ 3.5.

Method:
1. Per log_2(n) octave: compute mean(σ), std(σ), quantile values
2. Linear fit std(σ) = a + b·log(n)
3. Predict K_eff_band(q) = K_h + b·z_q at q ∈ {0.125, 0.375, 0.625, 0.875, 0.975}
4. Compare to empirical K_eff_band from prior s_mean(f) σ-quantile-stratified run
"""
import math
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

K_H = 3.0 / (math.log(4.0) - math.log(3.0))  # 10.43

OUT_DIR = Path(__file__).parent / "viz_outputs"
SAMPLE_CSV = OUT_DIR / "tail_sample_2pow27_n1000000.csv"

# Empirical K_eff_band values from R3 of test_s_mean_f.py
# (band_low, band_high, n_orbits, mean_sigma, K_eff_band)
EMPIRICAL_BANDS = [
    (0.00, 0.25, 246_046, 108.76,  6.1468),
    (0.25, 0.50, 248_592, 158.79,  9.2945),
    (0.50, 0.75, 250_533, 204.64, 11.9488),
    (0.75, 0.95, 204_448, 261.30, 15.2643),
    (0.95, 1.00,  50_381, 350.62, 19.7607),
]


def main():
    df = pl.read_csv(SAMPLE_CSV)
    n_arr = df["n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)
    log_n = np.log(n_arr.astype(np.float64))
    log2_n = log_n / math.log(2)
    print(f"Loaded {len(df):,} orbits at N=2^27.")
    print(f"K_h = {K_H:.5f}")
    print()

    # ============= Octave table =============
    print(f"{'='*100}")
    print(f"  Octave stats: mean(σ), std(σ), per-quantile σ values")
    print(f"{'='*100}")
    print(f"  {'octave':>9} {'mid':>6} {'n':>9} {'mean_σ':>9} {'std_σ':>9} "
          f"{'q12.5':>7} {'q37.5':>7} {'q62.5':>7} {'q87.5':>7} {'q97.5':>7}")

    octave_edges = np.arange(20, 28)
    octave_mid_log2 = (octave_edges[:-1] + octave_edges[1:]) / 2.0
    octave_mid_ln = octave_mid_log2 * math.log(2)

    means = []
    stds = []
    n_per_oct = []
    quantiles_per_oct = {q: [] for q in [0.125, 0.375, 0.625, 0.875, 0.975]}

    for i in range(len(octave_edges) - 1):
        lo, hi = octave_edges[i], octave_edges[i+1]
        mask = (log2_n >= lo) & (log2_n < hi)
        n_oct = int(mask.sum())
        if n_oct < 100:
            print(f"  [{lo:>2d},{hi:>2d}] {octave_mid_log2[i]:>6.2f} {n_oct:>9d}  (insufficient)")
            means.append(np.nan); stds.append(np.nan); n_per_oct.append(n_oct)
            for q in quantiles_per_oct:
                quantiles_per_oct[q].append(np.nan)
            continue
        sub = sigma[mask]
        mu = sub.mean()
        sd = sub.std()
        qs = {q: np.quantile(sub, q) for q in quantiles_per_oct}
        means.append(mu); stds.append(sd); n_per_oct.append(n_oct)
        for q, v in qs.items():
            quantiles_per_oct[q].append(v)
        print(f"  [{lo:>2d},{hi:>2d}] {octave_mid_log2[i]:>6.2f} {n_oct:>9d} "
              f"{mu:>9.2f} {sd:>9.2f} "
              f"{qs[0.125]:>7.1f} {qs[0.375]:>7.1f} {qs[0.625]:>7.1f} "
              f"{qs[0.875]:>7.1f} {qs[0.975]:>7.1f}")

    means = np.array(means); stds = np.array(stds)
    valid = ~np.isnan(stds)

    # ============= Linear fit std(σ) = a + b·log(n) =============
    print(f"\n{'='*100}")
    print(f"  Linear fit: std(σ) = a + b·log(n)  (using natural log)")
    print(f"{'='*100}")
    slope_b, intercept_a = np.polyfit(octave_mid_ln[valid], stds[valid], 1)
    yhat = slope_b * octave_mid_ln[valid] + intercept_a
    ss_res = np.sum((stds[valid] - yhat)**2)
    ss_tot = np.sum((stds[valid] - stds[valid].mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    print(f"  std(σ) ≈ {intercept_a:.4f} + {slope_b:.4f} · log(n)")
    print(f"  R² = {r2:.4f}")
    print(f"  Predicted b ≈ 3.5 (from K_eff_band slope ≈ 2.7 / quartile-z-shift ≈ 0.78)")
    print(f"  Empirical b = {slope_b:.4f}, gap = {slope_b - 3.5:+.4f}")

    # Per-octave residuals
    print(f"\n  Per-octave fit residuals:")
    print(f"  {'octave':>8} {'std_obs':>9} {'std_fit':>9} {'resid':>9}")
    for i, oct_mid in enumerate(octave_mid_ln):
        if np.isnan(stds[i]): continue
        fit = slope_b * oct_mid + intercept_a
        resid = stds[i] - fit
        print(f"  {octave_mid_log2[i]:>8.2f} {stds[i]:>9.3f} {fit:>9.3f} {resid:>+9.3f}")

    # ============= Predicted vs empirical K_eff_band =============
    print(f"\n{'='*100}")
    print(f"  Predicted vs empirical K_eff_band(q) = K_h + b·z_q")
    print(f"  Using fitted b = {slope_b:.4f}, K_h = {K_H:.4f}")
    print(f"{'='*100}")

    # Need quantiles of σ-band CENTERS within the GLOBAL σ distribution
    # (not per-octave). The empirical bands use global quantiles.
    band_centers = [0.125, 0.375, 0.625, 0.875, 0.975]
    z_at_q = [stats.norm.ppf(q) for q in band_centers]
    print(f"  {'q_center':>10} {'z(q)':>8} {'pred K_eff':>12} {'emp K_eff':>11} {'gap':>10}")
    for q, z, (lo, hi, n, ms, k_emp) in zip(band_centers, z_at_q, EMPIRICAL_BANDS):
        pred = K_H + slope_b * z
        gap = k_emp - pred
        print(f"  {q:>10.3f} {z:>+8.3f} {pred:>12.4f} {k_emp:>11.4f} {gap:>+10.4f}")

    # Also test: linear fit per-quantile?
    print(f"\n{'='*100}")
    print(f"  Per-quantile slopes from octave table (sanity check K_eff_band)")
    print(f"  Predicts: K_eff at quantile q = K_h + b·z_q if Gaussian model holds")
    print(f"{'='*100}")
    print(f"  {'q':>7} {'z(q)':>8} {'fitted slope':>14} {'pred slope (K_h+b·z)':>22} {'gap':>10}")
    for q in band_centers:
        z = stats.norm.ppf(q)
        vals = np.array(quantiles_per_oct[q])
        valid_q = ~np.isnan(vals)
        if valid_q.sum() < 3: continue
        slope_q, _ = np.polyfit(octave_mid_ln[valid_q], vals[valid_q], 1)
        pred_q = K_H + slope_b * z
        print(f"  {q:>7.3f} {z:>+8.3f} {slope_q:>14.4f} {pred_q:>22.4f} {slope_q - pred_q:>+10.4f}")

    # GPD check: is σ Gaussian or heavier-tailed?
    print(f"\n{'='*100}")
    print(f"  Gaussian check: per-octave skew & kurt of σ distribution")
    print(f"  (Gaussian: skew=0, excess_kurt=0; positive excess_kurt = heavier tail)")
    print(f"{'='*100}")
    print(f"  {'octave':>8} {'skew':>8} {'excess_kurt':>13}")
    for i in range(len(octave_edges) - 1):
        lo, hi = octave_edges[i], octave_edges[i+1]
        mask = (log2_n >= lo) & (log2_n < hi)
        if mask.sum() < 100: continue
        sub = sigma[mask]
        skew = stats.skew(sub)
        kurt = stats.kurtosis(sub)  # Fisher's: 0 for Gaussian
        print(f"  {octave_mid_log2[i]:>8.2f} {skew:>+8.3f} {kurt:>+13.3f}")


if __name__ == "__main__":
    main()
