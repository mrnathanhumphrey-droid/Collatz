"""
Experiment 34 — replicate s_mean / alpha_det match at k=10 and k=12.

The s_mean result at k=8 produced slope = 1.000 +/- 0.001 at K = 10.4282
and offset matching Tao (5.15) leading term within ~1 step (raw mean) or
0.3% (1%-trimmed). This experiment confirms the result at higher modular
resolution.

Same single walk as exp 33; groupby at k=8, 10, 12; report slope/offset/
K_implied for raw mean, trim-1%, and median (for reference).

Usage:
    python 34_alpha_det_K_calibration_mean_k_sweep.py --N 134217728
"""
import argparse
import gc
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

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


def find_K_star(prefix_arr, log_af_norm, s_target):
    def slope_at(K):
        return fit_slope_offset(prefix_arr + K * log_af_norm, s_target)[0]
    lo, hi = 5.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        sl = slope_at(mid)
        if sl > 1.0:
            lo = mid
        else:
            hi = mid
    return mid


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

    s_by_thresh = np.zeros((len(n_arr), 4), dtype=np.int32)
    inv_sort = np.argsort(sort_idx, axis=1)
    for col in range(4):
        rev = inv_sort[:, col]
        s_by_thresh[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()

    del raw, sort_idx, inv_sort, thresholds_sorted, s_arr
    del sqrt_n, n_two_third, sqrt_log, sqrt_div_log
    gc.collect()

    K_h = 3.0 / (np.log(4.0) - np.log(3.0))

    log_n_avg_global = float(log_n.mean())
    log_log_n = float(np.log(log_n_avg_global))
    delta_logs = {
        "N^(2/3)": log_n_avg_global / 3.0,
        "sqrt(N)*log(N)": log_n_avg_global / 2.0 - log_log_n,
        "sqrt(N)": log_n_avg_global / 2.0,
        "sqrt(N)/log(N)": log_n_avg_global / 2.0 + log_log_n,
    }

    for k in [8, 10, 12]:
        M = 1 << k
        K_classes = M // 2
        print()
        print("=" * 100)
        print(f"k = {k}  (mod {M}, {K_classes} odd classes)")
        print("=" * 100)
        t_k = time.perf_counter()
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

        aggs = [
            pl.col("log_n").mean().alias("log_n_mean"),
            pl.len().alias("n_per_class"),
        ]
        for tc in range(4):
            col = f"s{tc}"
            aggs.append(pl.col(col).median().alias(f"{col}_med"))
            aggs.append(pl.col(col).mean().alias(f"{col}_mean"))
            aggs.append(
                pl.col(col).filter(pl.col(col) < pl.col(col).quantile(0.99)).mean()
                .alias(f"{col}_mean_t1")
            )
        agg = df_g.group_by("class_idx").agg(aggs).sort("class_idx")
        del df_g
        gc.collect()
        print(f"  groupby in {time.perf_counter()-t_k:.1f}s", flush=True)

        # Materialize
        ci = agg["class_idx"].to_numpy()
        log_n_mean_per = np.zeros(K_classes)
        log_n_mean_per[ci] = agg["log_n_mean"].to_numpy()
        s_mean_per = np.zeros((K_classes, 4))
        s_mean_t1_per = np.zeros((K_classes, 4))
        s_med_per = np.zeros((K_classes, 4))
        for tc in range(4):
            s_med_per[ci, tc] = agg[f"s{tc}_med"].to_numpy()
            s_mean_per[ci, tc] = agg[f"s{tc}_mean"].to_numpy()
            s_mean_t1_per[ci, tc] = agg[f"s{tc}_mean_t1"].to_numpy()

        # Compute prefix and a_final
        prefix_arr = np.zeros(K_classes, dtype=np.int64)
        a_final_arr = np.zeros(K_classes, dtype=np.int64)
        for kk in range(K_classes):
            r = 2 * kk + 1
            sp, a_f, _ = deterministic_prefix(r, M)
            prefix_arr[kk] = sp
            a_final_arr[kk] = a_f
        log_af_norm = np.log(a_final_arr / float(M))

        for tcol, tname in enumerate(threshold_names):
            delta_log = delta_logs[tname]
            tao_pred = K_h * delta_log
            print()
            print(f"  --- {tname}  (Tao leading = K_h * delta_log = {tao_pred:.3f}) ---")
            print(f"  {'aggregation':>12} {'slope@K_h':>10} {'offset@K_h':>11} "
                  f"{'offset-Tao':>11} {'K*(slope=1)':>12} {'offset@K*':>10}")
            for label, arr in [
                ("median",    s_med_per),
                ("raw_mean",  s_mean_per),
                ("trim_1%",   s_mean_t1_per),
            ]:
                s_target = arr[:, tcol]
                slope_h, offset_h, _ = fit_slope_offset(prefix_arr + K_h * log_af_norm, s_target)
                K_star = find_K_star(prefix_arr, log_af_norm, s_target)
                slope_star, offset_star, _ = fit_slope_offset(
                    prefix_arr + K_star * log_af_norm, s_target)
                print(f"  {label:>12} {slope_h:>10.4f} {offset_h:>11.3f} "
                      f"{offset_h - tao_pred:>+11.3f} {K_star:>12.4f} "
                      f"{offset_star:>10.3f}")


if __name__ == "__main__":
    main()
