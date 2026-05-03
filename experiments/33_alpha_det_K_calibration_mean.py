"""
Experiment 33 — alpha_det K-calibration test against s_mean (not s_median).

Companion to exp 32. If the slope < 1 in s_median ~ alpha_det reflects the
distribution being right-skewed (median < mean), then redoing with s_mean
should drive the slope toward 1 at K = 10.43 (textbook heuristic) or
K = 10.55 (E[v]-corrected).

Pipeline:
  - Walk all odd N <= 2^27, record first-passage step count for 4 thresholds.
  - Group by N mod 256 (k=8), compute MEAN s per class (not median).
  - For each K in a grid, recompute alpha_det and fit s_mean ~ alpha_det.
  - Report slope, offset, K_implied across K's, find K* where slope=1.

If K* lands near 10.43 or 10.55, the median-vs-mean framing is confirmed
and the slope-not-1 was distribution-shape, not calibration.

Usage:
    python 33_alpha_det_K_calibration_mean.py --N 134217728
"""
import argparse
import gc
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import spearmanr

sys.stdout.reconfigure(encoding="utf-8")


@njit(parallel=True, cache=True)
def first_passage_multi(starts, thresholds_per_n):
    n = len(starts)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    for i in prange(n):
        N = starts[i]
        m = N
        steps = 0
        next_idx = 0
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0
            next_idx += 1
        max_steps = 20000
        while next_idx < 4 and steps < max_steps:
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            steps += 1
            while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                s_arr[i, next_idx] = steps
                next_idx += 1
    return s_arr


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
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1 << 27)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    print(f"[load] N = {args.N:,}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    print(f"        odd rows: {len(n_arr):,}", flush=True)

    # Build thresholds (same as exp 31)
    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)

    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)
    thresholds_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    threshold_names = ["N^(2/3)", "sqrt(N)*log(N)", "sqrt(N)", "sqrt(N)/log(N)"]

    print(f"[run]   walking trajectories ...", flush=True)
    t0 = time.perf_counter()
    s_arr = first_passage_multi(n_arr, thresholds_sorted)
    print(f"[run]   done in {time.perf_counter()-t0:.1f}s", flush=True)

    # Re-route to physical threshold ordering
    s_by_thresh = np.zeros((len(n_arr), 4), dtype=np.int32)
    inv_sort = np.argsort(sort_idx, axis=1)
    for col in range(4):
        rev = inv_sort[:, col]
        s_by_thresh[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()

    del raw, sort_idx, inv_sort, thresholds_sorted, s_arr
    del sqrt_n, n_two_third, sqrt_log, sqrt_div_log
    gc.collect()

    # Groupby at k=8 with both mean and median
    M = 256
    K_classes = M // 2
    res = (n_arr % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    df_g = pl.DataFrame({
        "class_idx": class_idx,
        "log_n": log_n,
        "s0": s_by_thresh[:, 0],
        "s1": s_by_thresh[:, 1],
        "s2": s_by_thresh[:, 2],
        "s3": s_by_thresh[:, 3],
    })

    # Aggregations: median, raw mean, trimmed-99% mean, trimmed-99.9% mean.
    aggs = [
        pl.col("log_n").mean().alias("log_n_mean"),
        pl.len().alias("n_per_class"),
    ]
    for tc in range(4):
        col = f"s{tc}"
        aggs.append(pl.col(col).median().alias(f"{col}_med"))
        aggs.append(pl.col(col).mean().alias(f"{col}_mean"))
        # Trimmed mean at 1% (drop values >= 99th percentile)
        aggs.append(
            pl.col(col).filter(pl.col(col) < pl.col(col).quantile(0.99)).mean()
            .alias(f"{col}_mean_t1")
        )
        # Trimmed mean at 0.1%
        aggs.append(
            pl.col(col).filter(pl.col(col) < pl.col(col).quantile(0.999)).mean()
            .alias(f"{col}_mean_t01")
        )

    agg = df_g.group_by("class_idx").agg(aggs).sort("class_idx")
    del df_g
    gc.collect()

    # Materialize into arrays
    log_n_mean_per = np.zeros(K_classes)
    n_per_class = np.zeros(K_classes, dtype=np.int64)
    s_mean_per = np.zeros((K_classes, 4))
    s_mean_t1_per = np.zeros((K_classes, 4))
    s_mean_t01_per = np.zeros((K_classes, 4))
    s_med_per = np.zeros((K_classes, 4))
    ci = agg["class_idx"].to_numpy()
    log_n_mean_per[ci] = agg["log_n_mean"].to_numpy()
    n_per_class[ci] = agg["n_per_class"].to_numpy()
    for tc in range(4):
        s_med_per[ci, tc] = agg[f"s{tc}_med"].to_numpy()
        s_mean_per[ci, tc] = agg[f"s{tc}_mean"].to_numpy()
        s_mean_t1_per[ci, tc] = agg[f"s{tc}_mean_t1"].to_numpy()
        s_mean_t01_per[ci, tc] = agg[f"s{tc}_mean_t01"].to_numpy()

    # Compute prefix and a_final
    prefix_arr = np.zeros(K_classes, dtype=np.int64)
    a_final_arr = np.zeros(K_classes, dtype=np.int64)
    for kk in range(K_classes):
        r = 2 * kk + 1
        s, a_f, _ = deterministic_prefix(r, M)
        prefix_arr[kk] = s
        a_final_arr[kk] = a_f
    log_af_norm = np.log(a_final_arr / float(M))

    log_n_avg = float(log_n_mean_per.mean())
    log_log_n = float(np.log(log_n_avg))
    delta_logs = {
        "N^(2/3)": log_n_avg / 3.0,
        "sqrt(N)*log(N)": log_n_avg / 2.0 - log_log_n,
        "sqrt(N)": log_n_avg / 2.0,
        "sqrt(N)/log(N)": log_n_avg / 2.0 + log_log_n,
    }

    K_h = 3.0 / (np.log(4.0) - np.log(3.0))
    K_grid = [9.30, 9.59, 9.70, 10.00, 10.22, 10.4282, 10.50, 10.55, 10.67, 10.80, 11.00, 11.50]

    def find_K_star(s_target):
        # Find K such that slope(K) = 1. Slope decreases monotonically with K
        # in this regime, so use binary search with correct direction.
        def slope_at(K):
            return fit_slope_offset(prefix_arr + K * log_af_norm, s_target)[0]
        lo, hi = 5.0, 20.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            sl = slope_at(mid)
            if sl > 1.0:  # slope too high -> K too low -> raise lo
                lo = mid
            else:         # slope too low -> K too high -> lower hi
                hi = mid
        return mid

    # Side-by-side comparison: median vs raw mean vs 1%-trim vs 0.1%-trim
    print()
    print("=" * 100)
    print("Side-by-side: at K_heuristic = 10.4282 and at K* (slope=1)")
    print("=" * 100)
    targets = [
        ("median",        s_med_per),
        ("raw_mean",      s_mean_per),
        ("trim_1%",       s_mean_t1_per),
        ("trim_0.1%",     s_mean_t01_per),
    ]
    for tcol, tname in enumerate(threshold_names):
        delta_log = delta_logs[tname]
        # Tao (5.15) leading term: log(N/threshold) / log(4/3)
        # threshold = N^(2/3): log(N/N^(2/3))/log(4/3) = log(N)/3/log(4/3) = K_h * delta_log / 3? No.
        # Actually: log(N/threshold) = log(N) - log(threshold) = delta_log (we already use this).
        # Tao says T ~ delta_log / log(4/3). And K_h = 3/log(4/3), so
        # Tao_leading = delta_log / log(4/3) = K_h * delta_log / 3.
        # But we measure first-passage in COLLATZ STEPS (= odd + even), which is
        # 3 * # odd-Syracuse steps on Geom(1/2). So:
        # Collatz-step prediction = 3 * delta_log / log(4/3) = K_h * delta_log.
        tao_leading = K_h * delta_log
        print()
        print(f"=== Threshold: {tname}   (delta_log={delta_log:.4f}, "
              f"Tao leading = K_h * delta_log = {tao_leading:.3f}) ===")
        print(f"  {'aggregation':>12} {'slope@K_h':>10} {'offset@K_h':>11} "
              f"{'K*(slope=1)':>12} {'offset@K*':>10} {'K_impl@K*':>11}")
        for label, arr in targets:
            s_target = arr[:, tcol]
            slope_h, offset_h, _ = fit_slope_offset(prefix_arr + K_h * log_af_norm, s_target)
            K_star = find_K_star(s_target)
            slope_star, offset_star, _ = fit_slope_offset(
                prefix_arr + K_star * log_af_norm, s_target)
            K_impl_star = offset_star / delta_log
            print(f"  {label:>12} {slope_h:>10.4f} {offset_h:>11.3f} "
                  f"{K_star:>12.4f} {offset_star:>10.3f} {K_impl_star:>11.4f}")

    # Save
    out_cols = {
        "class_idx": np.arange(K_classes),
        "residue": np.arange(1, M, 2),
        "a_final": a_final_arr,
        "prefix_steps": prefix_arr,
        "n_per_class": n_per_class,
        "log_n_mean": log_n_mean_per,
    }
    for i, tn in enumerate(threshold_names):
        out_cols[f"s_median_{tn}"] = s_med_per[:, i]
        out_cols[f"s_mean_{tn}"] = s_mean_per[:, i]
        out_cols[f"s_mean_t1_{tn}"] = s_mean_t1_per[:, i]
        out_cols[f"s_mean_t01_{tn}"] = s_mean_t01_per[:, i]
    out = pl.DataFrame(out_cols)
    out_dir = here.parent / "experiments_output"
    out_path = out_dir / f"33_first_passage_mean_N{args.N}.csv"
    out.write_csv(out_path)
    print()
    print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
