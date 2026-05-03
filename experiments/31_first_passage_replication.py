"""
Experiment 31 — first-passage / alpha_det replication.

Three sub-tests in one pass on N=2^27 odd data:

  (1) Higher modular resolution: k=8 (mod 256, 128 classes), k=10 (mod 1024,
      512 classes), k=12 (mod 4096, 2048 classes). Does Spearman rho stay
      ~1.0?

  (2) Threshold variation: f(N) in {N^(2/3), sqrt(N)*log(N), sqrt(N),
      sqrt(N)/log(N)}. For each N, walk Collatz once and record first-passage
      step count for each threshold. Tests whether alpha_det predicts
      first-passage timing universally or only at the sqrt(N) point.

  (3) Slope and offset of linear fit s_median(r) ~ slope * alpha_det(r) +
      offset. Is slope exactly 1? Does offset match K * <log N> * delta,
      where delta is the log-fraction implied by the threshold? Implied K
      from offset compared to per-octave beta finding (~10.55).

Usage:
    python 31_first_passage_replication.py --N 134217728
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import spearmanr

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))


@njit(parallel=True, cache=True)
def first_passage_multi(starts, thresholds_per_n):
    """For each odd start N and 4 thresholds (sorted descending), walk Collatz
    until each threshold has been first-crossed. thresholds_per_n shape
    (n, 4): for row i, sorted thresholds for start[i] descending.
    Returns s_arr shape (n, 4) of step counts."""
    n = len(starts)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    for i in prange(n):
        N = starts[i]
        m = N
        steps = 0
        next_idx = 0
        # If N is already <= largest threshold (impossible for our choices on
        # N >= 3), skip ahead.
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


def compute_alpha_det(k):
    M = 1 << k
    K = M // 2
    alpha = np.zeros(K)
    a_final = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r = 2 * kk + 1
        s, a_f, _ = deterministic_prefix(r, M)
        alpha[kk] = s + LOG_FACTOR_ODD * np.log(a_f / float(M))
        a_final[kk] = a_f
    return alpha, a_final


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
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    print(f"        odd rows: {len(n_arr):,}  loaded in "
          f"{time.perf_counter()-t0:.1f}s", flush=True)

    # Build per-n thresholds (sorted descending for the kernel).
    # f(N) ordering for descent (largest -> smallest):
    #   1. N^(2/3)   (largest)
    #   2. sqrt(N) * log(N)
    #   3. sqrt(N)
    #   4. sqrt(N) / log(N)   (smallest)
    print("[setup] building thresholds ...", flush=True)
    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)

    # Sorted descending per N: typically [N^(2/3), sqrt*log, sqrt, sqrt/log]
    # but verify ordering: at log_n=17.7, sqrt_n=11586 etc; N=2^27, N^(2/3)=2^18=262144
    # sqrt*log = 11586*17.7=205100. So sqrt_log < N^(2/3). Order may flip for small N.
    # Stack into (n, 4) and sort each row descending.
    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)  # descending per row
    thresholds_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    # We need to know which physical threshold each column represents per row.
    # Since the typical ordering should be stable for moderate-to-large N,
    # we record the dominant ordering.
    typical_order = np.argsort(-raw[len(n_arr) // 2])
    threshold_names = ["N^(2/3)", "sqrt(N)*log(N)", "sqrt(N)", "sqrt(N)/log(N)"]
    print(f"[setup] dominant threshold ordering (largest -> smallest):")
    for col, idx in enumerate(typical_order):
        print(f"          col {col}: {threshold_names[idx]}")

    # Run first-passage walk
    print(f"[run]   walking trajectories with 4 thresholds each ...", flush=True)
    t0 = time.perf_counter()
    s_arr = first_passage_multi(n_arr, thresholds_sorted)
    t1 = time.perf_counter()
    print(f"[run]   done in {t1 - t0:.1f}s   any uncrossed: "
          f"{int((s_arr == -1).sum()):,}", flush=True)

    # We want results indexed by physical threshold (not sorted-column).
    # Re-route: for each row, s_arr columns correspond to thresholds_sorted's
    # ordering (which we captured in sort_idx). Build s_by_threshold:
    # s_by_threshold[i, raw_index] = s_arr[i, sorted_column].
    # raw_index 0,1,2,3 = [N^(2/3), sqrt*log, sqrt, sqrt/log]
    s_by_thresh = np.zeros((len(n_arr), 4), dtype=np.int32)
    inv_sort = np.argsort(sort_idx, axis=1)
    for col in range(4):
        rev = inv_sort[:, col]
        s_by_thresh[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()

    # Free large intermediates we no longer need
    del raw, sort_idx, inv_sort, thresholds_sorted, s_arr
    del sqrt_n, n_two_third, sqrt_log, sqrt_div_log
    import gc
    gc.collect()

    print()

    # For each k in {8, 10, 12}, for each threshold, group and compute
    # Spearman rho, slope, offset. Use polars groupby for speed.
    for k in [8, 10, 12]:
        M = 1 << k
        K = M // 2
        print(f"=== k = {k}  (mod {M}, {K} odd classes) ===", flush=True)
        t_k = time.perf_counter()
        res = (n_arr % M).astype(np.int32)
        class_idx = ((res - 1) // 2).astype(np.int32)

        alpha_det, a_final = compute_alpha_det(k)

        # Build polars frame and group
        df_g = pl.DataFrame({
            "class_idx": class_idx,
            "log_n": log_n,
            "s0": s_by_thresh[:, 0],
            "s1": s_by_thresh[:, 1],
            "s2": s_by_thresh[:, 2],
            "s3": s_by_thresh[:, 3],
        })
        agg = df_g.group_by("class_idx").agg([
            pl.col("log_n").mean().alias("log_n_mean"),
            pl.len().alias("n_per_class"),
            pl.col("s0").median().alias("s0_med"),
            pl.col("s1").median().alias("s1_med"),
            pl.col("s2").median().alias("s2_med"),
            pl.col("s3").median().alias("s3_med"),
        ]).sort("class_idx")
        del df_g
        gc.collect()

        # Materialize into arrays indexed by class_idx 0..K-1
        log_n_mean_per = np.zeros(K)
        n_per_class = np.zeros(K, dtype=np.int64)
        s_med_per = np.zeros((K, 4), dtype=np.float64)
        ci = agg["class_idx"].to_numpy()
        log_n_mean_per[ci] = agg["log_n_mean"].to_numpy()
        n_per_class[ci] = agg["n_per_class"].to_numpy()
        for tc, col in enumerate(["s0_med", "s1_med", "s2_med", "s3_med"]):
            s_med_per[ci, tc] = agg[col].to_numpy()
        del agg
        gc.collect()
        print(f"  groupby completed in {time.perf_counter()-t_k:.1f}s   "
              f"per-class min {n_per_class[n_per_class > 0].min():,}", flush=True)

        # For each physical threshold compute per-class median s
        print(f"  {'threshold':>22} {'expected delta_logN':>22} "
              f"{'Spearman_rho':>13} {'slope':>9} {'offset':>9} "
              f"{'R^2':>7} {'K_implied':>10}")

        for tcol, tname in enumerate(threshold_names):
            mask_have = n_per_class > 0
            ad = alpha_det[mask_have]
            sm = s_med_per[mask_have, tcol]
            log_n_avg = float(log_n_mean_per[mask_have].mean())

            # Expected slope and intercept under the heuristic:
            # s = alpha_det + K_const * delta_logN
            # where delta_logN depends on threshold.
            # f(N) = N^(2/3): delta_log = log N - log(N^(2/3)) = log N / 3
            # f(N) = sqrt(N)*log N: delta_log = log N / 2 - log log N
            # f(N) = sqrt(N): delta_log = log N / 2
            # f(N) = sqrt(N)/log N: delta_log = log N / 2 + log log N
            log_log_n = np.log(log_n_avg)
            delta_log = {
                "N^(2/3)": log_n_avg / 3.0,
                "sqrt(N)*log(N)": log_n_avg / 2.0 - log_log_n,
                "sqrt(N)": log_n_avg / 2.0,
                "sqrt(N)/log(N)": log_n_avg / 2.0 + log_log_n,
            }[tname]

            rho, _ = spearmanr(ad, sm)
            slope, offset, r2 = fit_slope_offset(ad, sm)
            K_implied = offset / delta_log if delta_log > 0 else float("nan")

            print(f"  {tname:>22} {delta_log:>22.4f} "
                  f"{rho:>13.4f} {slope:>9.4f} {offset:>9.3f} "
                  f"{r2:>7.4f} {K_implied:>10.4f}")
        print()

    # Save consolidated CSV (k=8 detail with all thresholds)
    print("[save] consolidated per-class data at k=8 ...", flush=True)
    M = 256
    K = M // 2
    res = (n_arr % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)
    alpha_det_8, a_final_8 = compute_alpha_det(8)

    rows = []
    for kk in range(K):
        mask = class_idx == kk
        if mask.sum() < 30:
            continue
        rec = {
            "class_idx": kk,
            "residue": 2 * kk + 1,
            "a_final": int(a_final_8[kk]),
            "alpha_det": float(alpha_det_8[kk]),
            "n_per_class": int(mask.sum()),
            "log_n_mean": float(log_n[mask].mean()),
        }
        for tcol, tname in enumerate(threshold_names):
            rec[f"s_median_{tname}"] = float(np.median(s_by_thresh[mask, tcol]))
        rows.append(rec)

    out = pl.DataFrame(rows)
    out_dir = here.parent / "experiments_output"
    out_path = out_dir / f"31_first_passage_replication_N{args.N}.csv"
    out.write_csv(out_path)
    print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
