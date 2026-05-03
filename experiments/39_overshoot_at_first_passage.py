"""
Experiment 39 — overshoot at first-passage.

Test the decomposition  ε(f) = ε(σ) + K_h * <overshoot(f)>  where
overshoot(n; f) = log(f) - log(v*(n; f))  and  v*(n; f) is the first orbit
value <= f along n's Collatz trajectory.

If this decomposition holds:
  <overshoot(f)> = (ε(f) - ε(σ)) / K_h
                = (ε(f) + 2.45) / K_h    (using ε(σ) ≈ -2.45)

For ε(f) ≈ -2.35 + 0.486 * log(f):  predicted <overshoot> ≈ (0.10 + 0.486*log f) / K_h
                                                          ≈ 0.0096 + 0.0466 * log f

Empirical test: walk all odd n at N=2^27, record v*(n; f) at first crossing
for each of 4 thresholds. Compute <log f - log v*> per threshold and
compare to the predicted overshoot.

If empirical <overshoot> matches predicted to <0.05, the decomposition is
confirmed and the slope on log f is the slope of mean overshoot vs log f.

Usage:
    python 39_overshoot_at_first_passage.py
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))


@njit(parallel=True, cache=True)
def first_passage_with_overshoot(starts, thresholds_per_n):
    """For each odd start, walk Collatz to first crossing of each of 4
    thresholds (sorted descending). Returns:
      s_arr[i, j] = step count at first crossing of j-th threshold
      v_arr[i, j] = orbit value at first crossing of j-th threshold
    """
    n = len(starts)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    v_arr = np.zeros((n, 4), dtype=np.int64)
    for i in prange(n):
        N = starts[i]
        m = N
        steps = 0
        next_idx = 0
        # Initial check
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0
            v_arr[i, next_idx] = m
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
                v_arr[i, next_idx] = m
                next_idx += 1
    return s_arr, v_arr


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    N = 1 << 27
    print(f"[load] N = {N:,}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    print(f"  {len(n_arr):,} odd rows", flush=True)

    # Build thresholds as in exp 35
    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)
    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)
    thresholds_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    threshold_names = ["N^(2/3)", "sqrt(N)*log(N)", "sqrt(N)", "sqrt(N)/log(N)"]

    print("[run] walking trajectories with overshoot recording...", flush=True)
    t0 = time.perf_counter()
    s_arr, v_arr = first_passage_with_overshoot(n_arr, thresholds_sorted)
    print(f"  done in {time.perf_counter()-t0:.1f}s", flush=True)

    # Reorder to physical thresholds
    inv_sort = np.argsort(sort_idx, axis=1)
    v_by_thresh = np.zeros_like(v_arr)
    for col in range(4):
        rev = inv_sort[:, col]
        v_by_thresh[:, col] = np.take_along_axis(v_arr, rev[:, None], axis=1).flatten()
    threshold_by_n = np.zeros_like(thresholds_sorted)
    for col in range(4):
        rev = inv_sort[:, col]
        threshold_by_n[:, col] = np.take_along_axis(thresholds_sorted, rev[:, None], axis=1).flatten()
    del raw, sort_idx, inv_sort, thresholds_sorted, s_arr, v_arr
    del sqrt_n, n_two_third, sqrt_log, sqrt_div_log
    gc.collect()

    log_n_avg = float(log_n.mean())
    log_log_n = float(np.log(log_n_avg))
    delta_logs = {
        "N^(2/3)":         log_n_avg / 3.0,
        "sqrt(N)*log(N)":  log_n_avg / 2.0 - log_log_n,
        "sqrt(N)":         log_n_avg / 2.0,
        "sqrt(N)/log(N)":  log_n_avg / 2.0 + log_log_n,
    }
    # Empirical ε(f) values from exp 35 raw mean at N=2^27, k=8:
    eps_obs = {
        "N^(2/3)":        +3.009,
        "sqrt(N)*log(N)": +3.028,
        "sqrt(N)":        +2.161,
        "sqrt(N)/log(N)": +1.152,
    }
    eps_sigma = -2.45  # from TA.1

    print()
    print(f"  <log N> = {log_n_avg:.4f}")
    print()
    print(f"=== Overshoot at first crossing per threshold ===")
    print(f"  {'threshold':>16} {'log(f)_avg':>11} "
          f"{'<log(f/v*)>':>13} {'SE':>9} "
          f"{'predicted <oversh>':>20} {'matches?':>10}")

    for tcol, tname in enumerate(threshold_names):
        # Per-N threshold and crossing value
        f_per_n = threshold_by_n[:, tcol].astype(np.float64)
        v_star = v_by_thresh[:, tcol].astype(np.float64)
        log_f = np.log(np.maximum(f_per_n, 1.0))
        log_v_star = np.log(np.maximum(v_star, 1.0))
        # Overshoot = log(f) - log(v*)  (positive when v* < f)
        overshoot = log_f - log_v_star
        # Some entries might have v* > f if the threshold wasn't crossed within max_steps
        # (shouldn't happen at N=2^27 / threshold up to N^(2/3) ~ 2^18, but guard)
        valid = (v_star > 0) & (v_star <= f_per_n)
        if valid.sum() < len(overshoot):
            print(f"  WARNING: {(~valid).sum()} entries with invalid crossing for {tname}",
                  flush=True)
        os_valid = overshoot[valid]
        log_f_avg = float(log_f.mean())
        os_mean = float(os_valid.mean())
        os_se = float(os_valid.std() / np.sqrt(len(os_valid)))

        # Predicted from decomposition: (eps_obs - eps_sigma) / K_h
        eps_f = eps_obs[tname]
        os_pred = (eps_f - eps_sigma) / K_H
        match = "YES" if abs(os_mean - os_pred) < 0.05 else "no"

        print(f"  {tname:>16} {log_f_avg:>11.4f} "
              f"{os_mean:>13.5f} {os_se:>9.5f} "
              f"{os_pred:>20.5f} {match:>10}")

    print()
    print(f"=== Mean overshoot vs log(threshold) — fit slope ===")
    log_f_per_thresh = []
    overshoot_per_thresh = []
    for tcol, tname in enumerate(threshold_names):
        f_per_n = threshold_by_n[:, tcol].astype(np.float64)
        v_star = v_by_thresh[:, tcol].astype(np.float64)
        log_f = np.log(np.maximum(f_per_n, 1.0))
        log_v_star = np.log(np.maximum(v_star, 1.0))
        overshoot = log_f - log_v_star
        valid = (v_star > 0) & (v_star <= f_per_n)
        log_f_per_thresh.append(float(log_f.mean()))
        overshoot_per_thresh.append(float(overshoot[valid].mean()))
    log_f_per_thresh = np.array(log_f_per_thresh)
    overshoot_per_thresh = np.array(overshoot_per_thresh)

    # Linear fit overshoot ~ a + b * log(f)
    x_c = log_f_per_thresh - log_f_per_thresh.mean()
    y_c = overshoot_per_thresh - overshoot_per_thresh.mean()
    slope = (x_c * y_c).sum() / (x_c * x_c).sum()
    intercept = overshoot_per_thresh.mean() - slope * log_f_per_thresh.mean()
    print(f"  <overshoot> = {intercept:.4f} + {slope:.4f} * log(f)")
    print(f"  K_h * slope = {K_H * slope:.4f}    (compare to TA.3 slope on log(f) = +0.486)")
    print(f"  K_h * intercept = {K_H * intercept:.4f}    (compare to ε(σ) constant component ~ -2.35 - eps_sigma = +0.10)")


if __name__ == "__main__":
    main()
