"""
Experiment 32 — recalibrate alpha_det's K constant against the slope of
s_median vs alpha_det.

Hypothesis (user): the slope ~0.98 in s_median = slope * alpha_det + offset
is calibration error in alpha_det. alpha_det was built with K = 10.4282
(textbook heuristic). If we substitute the empirical K (~10.55 from the
E[v]-corrected K_pred, or whatever the data implies), the slope tightens
toward 1.

Test: sweep K in {9.5, 10.0, 10.22, 10.43, 10.55, 10.67, 10.80} and refit
s_median vs alpha_det(K) for each threshold at k=8. Report slope, offset,
K_implied (= offset / delta_logN). Find the K that drives slope to exactly 1.

Uses the existing CSV from experiment 31 (per-class prefix_steps, a_final,
log_n_mean, s_median per threshold).
"""
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix(r, a0, max_steps=400):
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return steps, a, c


def fit_slope_offset(x, y):
    x_c = x - x.mean()
    y_c = y - y.mean()
    slope = (x_c * y_c).sum() / (x_c * x_c).sum() if (x_c * x_c).sum() > 0 else 0.0
    offset = y.mean() - slope * x.mean()
    pred = slope * x + offset
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = (y_c ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(offset), float(r2)


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    csv_path = out_dir / "31_first_passage_replication_N134217728.csv"
    df = pl.read_csv(csv_path)
    print(f"[load] {csv_path}   {len(df)} classes", flush=True)

    # Recompute prefix_steps per class (since CSV has alpha_det at K=10.4282
    # but not prefix_steps directly).
    M = 256  # k=8
    prefix_arr = np.zeros(M // 2, dtype=np.int64)
    a_final_arr = np.zeros(M // 2, dtype=np.int64)
    for kk in range(M // 2):
        r = 2 * kk + 1
        s, a_f, _ = deterministic_prefix(r, M)
        prefix_arr[kk] = s
        a_final_arr[kk] = a_f

    # Sanity-check: alpha_det reconstructed at K=10.4282 should match CSV
    K_h = 3.0 / (np.log(4.0) - np.log(3.0))
    log_af_norm = np.log(a_final_arr / float(M))
    alpha_det_h = prefix_arr + K_h * log_af_norm

    csv_alpha = df.sort("class_idx")["alpha_det"].to_numpy()
    diff = np.max(np.abs(alpha_det_h - csv_alpha))
    print(f"[check] reconstructed alpha_det matches CSV to {diff:.2e}", flush=True)

    # log_n_avg same for all classes (since data is uniform over [1, 2^27])
    log_n_avg = float(df["log_n_mean"].mean())
    log_log_n = float(np.log(log_n_avg))
    threshold_names = ["N^(2/3)", "sqrt(N)*log(N)", "sqrt(N)", "sqrt(N)/log(N)"]
    delta_logs = {
        "N^(2/3)": log_n_avg / 3.0,
        "sqrt(N)*log(N)": log_n_avg / 2.0 - log_log_n,
        "sqrt(N)": log_n_avg / 2.0,
        "sqrt(N)/log(N)": log_n_avg / 2.0 + log_log_n,
    }

    # Sweep K
    K_grid = [9.30, 9.59, 9.70, 10.0, 10.22, 10.4282, 10.50, 10.55, 10.67, 10.80]

    df_sorted = df.sort("class_idx")
    s_med_cols = {
        "N^(2/3)": df_sorted["s_median_N^(2/3)"].to_numpy(),
        "sqrt(N)*log(N)": df_sorted["s_median_sqrt(N)*log(N)"].to_numpy(),
        "sqrt(N)": df_sorted["s_median_sqrt(N)"].to_numpy(),
        "sqrt(N)/log(N)": df_sorted["s_median_sqrt(N)/log(N)"].to_numpy(),
    }

    for tname in threshold_names:
        print()
        print(f"=== Threshold: {tname}   (delta_log = {delta_logs[tname]:.4f}) ===")
        print(f"  {'K_used_in_alpha_det':>20} {'slope':>9} {'offset':>9} "
              f"{'R^2':>7} {'K_implied_offset':>17}")
        s_med = s_med_cols[tname]
        for K in K_grid:
            alpha_det_K = prefix_arr + K * log_af_norm
            slope, offset, r2 = fit_slope_offset(alpha_det_K, s_med)
            K_implied = offset / delta_logs[tname]
            mark = "  <- heuristic" if abs(K - K_h) < 1e-3 else ""
            mark = "  <- E[v]≈1.99" if abs(K - 10.55) < 1e-3 else mark
            print(f"  {K:>20.4f} {slope:>9.4f} {offset:>9.3f} "
                  f"{r2:>7.4f} {K_implied:>17.4f}{mark}")

        # Find K* such that slope=1 exactly. Slope is a continuous, monotone
        # function of K (within the regime where Var(prefix) << K^2 Var(log)).
        # Use binary search.
        def slope_at_K(K):
            alpha = prefix_arr + K * log_af_norm
            return fit_slope_offset(alpha, s_med)[0]

        # Bracket from grid: find adjacent K values straddling slope=1
        lo, hi = 5.0, 15.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if slope_at_K(mid) < 1.0:
                lo = mid
            else:
                hi = mid
        K_star = mid
        slope_star, offset_star, r2_star = fit_slope_offset(
            prefix_arr + K_star * log_af_norm, s_med)
        K_implied_star = offset_star / delta_logs[tname]
        print(f"  ---")
        print(f"  K* (slope = 1.0): {K_star:.4f}    "
              f"offset = {offset_star:.3f}   K_implied = {K_implied_star:.4f}   "
              f"slope = {slope_star:.6f}")

    print()
    print("=== Summary ===")
    print("If K_emp lies in a single value across thresholds, the slope-1")
    print("calibration recovers a single empirical descent constant.")
    print("If K_emp varies with threshold, there's threshold-dependent")
    print("behavior (median descent rate depends on how deep we go).")


if __name__ == "__main__":
    main()
