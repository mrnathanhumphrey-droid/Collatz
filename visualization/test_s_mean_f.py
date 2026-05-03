"""
test_s_mean_f.py - First-passage s(n; f) for σ-quantile + body-slope discrimination.

For each odd n in the N=2^27 sample, compute s(n; f) = smallest step t with v_t ≤ f,
for f ∈ {2^5, 2^10, 2^15, 2^20}. Then run two regressions:

(R1) s_mean(f) on log(n) at fixed f
     If body-slope hypothesis holds: predicts ~10.5 (median σ slope from tail-shape test).
     If matches K_h = 10.43 instead: body-slope hypothesis falsified.

(R2) s_mean(f) on log(f) at fixed log(n) bin
     Slope ≈ -K_h. Compare across log(n) bins to see if K_eff = 9.94 is
     body-slope + correction or a separate quantity.
"""
import math
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

K_H = 3.0 / (math.log(4.0) - math.log(3.0))  # 10.43

OUT_DIR = Path(__file__).parent / "viz_outputs"
SAMPLE_CSV = OUT_DIR / "tail_sample_2pow27_n1000000.csv"


@njit(parallel=True, cache=True)
def walk_with_first_passage(starts, thresholds, max_steps):
    """Walk Collatz from each start; record first-passage step for each threshold.
    Thresholds must be sorted descending. Returns (sigma, s_arr).
    s_arr[i, j] = step at which orbit i first satisfies v_t <= thresholds[j].
    sigma[i] = total steps to reach 1.
    """
    n_orbits = len(starts)
    n_thresh = len(thresholds)
    sigmas = np.full(n_orbits, -1, dtype=np.int64)
    s_arr = np.full((n_orbits, n_thresh), -1, dtype=np.int64)

    for i in prange(n_orbits):
        n = starts[i]
        next_idx = 0
        # Initial check
        while next_idx < n_thresh and n <= thresholds[next_idx]:
            s_arr[i, next_idx] = 0
            next_idx += 1
        for step in range(1, max_steps + 1):
            if n == 1:
                while next_idx < n_thresh:
                    s_arr[i, next_idx] = step - 1
                    next_idx += 1
                sigmas[i] = step - 1
                break
            if n % 2 == 0:
                n = n >> 1
            else:
                n = 3 * n + 1
            while next_idx < n_thresh and n <= thresholds[next_idx]:
                s_arr[i, next_idx] = step
                next_idx += 1
    return sigmas, s_arr


def main():
    df = pl.read_csv(SAMPLE_CSV)
    n_arr = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n_arr.astype(np.float64))
    print(f"Loaded {len(n_arr):,} orbits from {SAMPLE_CSV.name}")

    thresholds = np.array([1 << 20, 1 << 15, 1 << 10, 1 << 5], dtype=np.int64)
    log_f = np.log(thresholds.astype(np.float64))
    log2_f = log_f / math.log(2)
    threshold_labels = [f"2^{int(x)}" for x in log2_f]

    # JIT warmup
    _ = walk_with_first_passage(np.array([7], dtype=np.int64), thresholds, 5000)

    print(f"\n[walk] computing s(n; f) for f ∈ {{2^20, 2^15, 2^10, 2^5}} ...")
    t0 = time.perf_counter()
    sigmas, s_arr = walk_with_first_passage(n_arr, thresholds, max_steps=20000)
    print(f"[walk] done in {time.perf_counter()-t0:.1f}s")

    valid = sigmas >= 0
    print(f"  converged: {int(valid.sum()):,} / {len(sigmas):,}")

    # Restrict to converged orbits
    n_v = n_arr[valid]
    log_n_v = log_n[valid]
    sigma_v = sigmas[valid]
    s_v = s_arr[valid]

    # Per-threshold global stats
    print(f"\n{'='*100}")
    print(f"  Per-threshold global stats")
    print(f"{'='*100}")
    print(f"  {'threshold':>12} {'log_2(f)':>10} {'<s(f)>':>9} {'<sigma - s(f)>':>17} {'n_with_s>0':>12}")
    s_means = []
    for j, f in enumerate(thresholds):
        sj = s_v[:, j]
        valid_j = sj > 0  # exclude orbits already below f
        if valid_j.sum() == 0:
            continue
        s_means.append(sj[valid_j].mean())
        rest = sigma_v[valid_j] - sj[valid_j]
        print(f"  {threshold_labels[j]:>12} {log2_f[j]:>10.2f} {sj[valid_j].mean():>9.2f} "
              f"{rest.mean():>17.2f} {int(valid_j.sum()):>12,}")

    # ============= R1: s_mean(f) on log(n) at fixed f =============
    print(f"\n{'='*100}")
    print(f"  R1: s_mean(f) regressed on log(n) at fixed f")
    print(f"      Predicted slopes:")
    print(f"        body-slope hypothesis (from per-quantile q50 ~ 10.2-10.8): ~10.5")
    print(f"        K_h baseline:                                              {K_H:.3f}")
    print(f"{'='*100}")

    # Bin by log_2(n) octave
    log2_n_v = log_n_v / math.log(2)
    octave_lo = 20
    octave_hi = 27
    octave_edges = np.arange(octave_lo, octave_hi + 1)
    octave_mid_log2 = (octave_edges[:-1] + octave_edges[1:]) / 2.0
    octave_mid_ln = octave_mid_log2 * math.log(2)

    print(f"  {'threshold':>12} {'slope/ln(n)':>13} {'gap K_h':>10} {'per-octave R^2':>15}")
    for j, f in enumerate(thresholds):
        sj = s_v[:, j]
        valid_j = sj > 0
        if valid_j.sum() < 1000:
            continue
        s_per_octave = []
        for i in range(len(octave_edges) - 1):
            lo, hi = octave_edges[i], octave_edges[i+1]
            mask = (log2_n_v >= lo) & (log2_n_v < hi) & valid_j
            if mask.sum() < 100:
                s_per_octave.append(np.nan)
            else:
                s_per_octave.append(sj[mask].mean())
        s_per_octave = np.array(s_per_octave)
        valid_oct = ~np.isnan(s_per_octave)
        if valid_oct.sum() < 3:
            continue
        slope, intercept = np.polyfit(octave_mid_ln[valid_oct], s_per_octave[valid_oct], 1)
        # R^2
        yhat = slope * octave_mid_ln[valid_oct] + intercept
        ss_res = np.sum((s_per_octave[valid_oct] - yhat) ** 2)
        ss_tot = np.sum((s_per_octave[valid_oct] - s_per_octave[valid_oct].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        gap = slope - K_H
        print(f"  {threshold_labels[j]:>12} {slope:>13.4f} {gap:>+10.4f} {r2:>15.4f}")

    # ============= R2: s_mean(f) on log(f) at fixed log(n) bin =============
    print(f"\n{'='*100}")
    print(f"  R2: s_mean(f) regressed on log(f) at fixed log_2(n) bin")
    print(f"      Predicted slope ≈ -K_h = -{K_H:.3f}.  K_h - K_eff = +0.486 if K_eff = 9.94.")
    print(f"      If slope is ~constant across bins: K_eff is universal, not bin-dependent.")
    print(f"      If slope varies: K_eff inherits bin-specific σ-distribution structure.")
    print(f"{'='*100}")

    print(f"  {'log_2(n) bin':>14} {'n_bin':>8} {'slope on log(f)':>17} {'-slope':>9} {'gap from K_h':>13}")
    for i in range(len(octave_edges) - 1):
        lo, hi = octave_edges[i], octave_edges[i+1]
        mask = (log2_n_v >= lo) & (log2_n_v < hi)
        if mask.sum() < 100:
            continue
        # For each threshold compute mean s(f) over this bin
        s_means_bin = []
        log_fs_bin = []
        for j, f in enumerate(thresholds):
            sj = s_v[mask, j]
            valid_j = sj > 0
            if valid_j.sum() < 30:
                continue
            s_means_bin.append(sj[valid_j].mean())
            log_fs_bin.append(log_f[j])
        if len(s_means_bin) < 3:
            continue
        s_means_bin = np.array(s_means_bin)
        log_fs_bin = np.array(log_fs_bin)
        slope, intercept = np.polyfit(log_fs_bin, s_means_bin, 1)
        gap = (-slope) - K_H
        print(f"  [{lo},{hi}]{'':>5}  {int(mask.sum()):>8d} {slope:>17.4f} {-slope:>9.4f} {gap:>+13.4f}")

    # ============= R3: σ-quantile-stratified s_mean(f) on log(f) =============
    print(f"\n{'='*100}")
    print(f"  R3: σ-quantile-stratified s_mean(f) regression")
    print(f"      For each σ-quantile band, slope of s_mean(f) on log(f) = K_eff_band")
    print(f"      Tests whether K_eff = 9.94 is body-driven (low-σ band) or tail-driven.")
    print(f"{'='*100}")

    # σ-quantile bands (within each log_2(n) bin to remove log_n confound)
    # Or: globally? Let's do globally for simplicity, then per-bin.
    sigma_quantiles = [0.0, 0.25, 0.5, 0.75, 0.95, 1.001]
    quantile_edges = np.quantile(sigma_v, sigma_quantiles[:-1])
    quantile_edges = np.append(quantile_edges, sigma_v.max() + 1)
    # Bands: (q[i], q[i+1]] for i in 0..len-1

    print(f"  σ-quantile bands (global): {sigma_quantiles}")
    print(f"  σ thresholds: {[f'{q:.0f}' for q in quantile_edges[:-1]]}")
    print(f"  {'band':>14} {'n_orbits':>10} {'mean_sigma':>11} {'slope/log(f)':>14} {'K_eff = -slope':>15} {'gap K_h':>10}")
    for b in range(len(sigma_quantiles) - 1):
        lo_q, hi_q = quantile_edges[b], quantile_edges[b+1]
        mask = (sigma_v >= lo_q) & (sigma_v < hi_q)
        if mask.sum() < 100:
            continue
        s_means_band = []
        log_fs_band = []
        for j, f in enumerate(thresholds):
            sj = s_v[mask, j]
            valid_j = sj > 0
            if valid_j.sum() < 30:
                continue
            s_means_band.append(sj[valid_j].mean())
            log_fs_band.append(log_f[j])
        if len(s_means_band) < 3:
            continue
        s_means_band = np.array(s_means_band)
        log_fs_band = np.array(log_fs_band)
        slope, intercept = np.polyfit(log_fs_band, s_means_band, 1)
        K_eff_band = -slope
        gap = K_eff_band - K_H
        band_label = f"[{sigma_quantiles[b]:.2f},{sigma_quantiles[b+1]:.2f}]"
        print(f"  {band_label:>14} {int(mask.sum()):>10,} {sigma_v[mask].mean():>11.2f} "
              f"{slope:>14.4f} {K_eff_band:>15.4f} {gap:>+10.4f}")


if __name__ == "__main__":
    main()
