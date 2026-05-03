"""
Experiment 37 — k=16 verification of s_mean ≈ α_det + K_h · log(N/f(N)).

Primary agent's exp 35 verified the structural identity at k ∈ {8, 10, 12, 14}
on N = 2^25 and 2^27. This pushes to k=16 (mod 65536, 32768 odd classes) at
N=2^27, where per-class size is ~2048 — marginal but adequate for slope
estimation.

Same observable set:
  σ                full descent to 1
  s @ N^(2/3)      first passage to N^(2/3)
  s @ √N · log N
  s @ √N
  s @ √N / log N

Slope at K_h should remain 1.000 ± 0.005 if universality extends. Offset gap
from Tao (5.15) leading term should remain stable across the resolution.

Usage:
    python 37_alpha_det_k16_verification.py
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
    return float(slope), float(offset)


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    N_val = 1 << 27  # 2^27 = 134,217,728
    print(f"# N = {N_val:,} = 2^27, k=16 (mod 65536, 32768 odd classes)")
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{N_val}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int32)
    print(f"  loaded {len(n_arr):,} odd rows in {time.perf_counter()-t0:.1f}s")

    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)
    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)
    thresholds_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)

    print(f"  walking trajectories ...")
    t0 = time.perf_counter()
    s_arr = first_passage_multi(n_arr, thresholds_sorted)
    print(f"  done in {time.perf_counter()-t0:.1f}s")

    s_by_thresh = np.zeros((len(n_arr), 4), dtype=np.int32)
    inv_sort = np.argsort(sort_idx, axis=1)
    for col in range(4):
        rev = inv_sort[:, col]
        s_by_thresh[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()
    del raw, sort_idx, inv_sort, thresholds_sorted, s_arr
    del sqrt_n, n_two_third, sqrt_log, sqrt_div_log
    gc.collect()

    log_n_avg = float(log_n.mean())
    log_log_n = float(np.log(log_n_avg))
    delta_logs = {
        "sigma":           log_n_avg,
        "s @ N^(2/3)":    log_n_avg / 3.0,
        "s @ sqrt*log":   log_n_avg / 2.0 - log_log_n,
        "s @ sqrt(N)":    log_n_avg / 2.0,
        "s @ sqrt/log":   log_n_avg / 2.0 + log_log_n,
    }

    k = 16
    M = 1 << k
    K_classes = M // 2
    per_class_min = (len(n_arr) // K_classes)
    print(f"\n  k={k}: per-class ~{per_class_min:,}")
    t_k = time.perf_counter()
    res = (n_arr % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    df_g = pl.DataFrame({
        "class_idx": class_idx,
        "log_n": log_n,
        "sigma": sigma_arr,
        "s0": s_by_thresh[:, 0],
        "s1": s_by_thresh[:, 1],
        "s2": s_by_thresh[:, 2],
        "s3": s_by_thresh[:, 3],
    })
    cols_for_aggs = ["sigma", "s0", "s1", "s2", "s3"]
    aggs = [pl.col("log_n").mean().alias("log_n_mean")]
    for col in cols_for_aggs:
        aggs.append(pl.col(col).mean().alias(f"{col}_mean"))
        aggs.append(
            pl.col(col).filter(pl.col(col) < pl.col(col).quantile(0.99)).mean()
            .alias(f"{col}_t1")
        )
    agg = df_g.group_by("class_idx").agg(aggs).sort("class_idx")
    print(f"    groupby in {time.perf_counter()-t_k:.1f}s, K_classes={K_classes}, populated classes: {len(agg)}")
    del df_g; gc.collect()

    ci = agg["class_idx"].to_numpy()
    means = {c: np.zeros(K_classes) for c in cols_for_aggs}
    t1s = {c: np.zeros(K_classes) for c in cols_for_aggs}
    for c in cols_for_aggs:
        means[c][ci] = agg[f"{c}_mean"].to_numpy()
        t1s[c][ci] = agg[f"{c}_t1"].to_numpy()

    # Compute alpha_det for all 32768 odd classes
    print(f"    computing α_det for {K_classes:,} classes ...")
    t_p = time.perf_counter()
    prefix_arr = np.zeros(K_classes, dtype=np.int64)
    a_final_arr = np.zeros(K_classes, dtype=np.int64)
    for kk in range(K_classes):
        r = 2 * kk + 1
        sp, a_f, _ = deterministic_prefix(r, M)
        prefix_arr[kk] = sp
        a_final_arr[kk] = a_f
    log_af_norm = np.log(a_final_arr / float(M))
    alpha_det_h = prefix_arr + K_H * log_af_norm
    print(f"    done in {time.perf_counter()-t_p:.1f}s")

    observables = [
        ("sigma",        "sigma"),
        ("s @ N^(2/3)",  "s0"),
        ("s @ sqrt*log", "s1"),
        ("s @ sqrt(N)",  "s2"),
        ("s @ sqrt/log", "s3"),
    ]
    delta_keys = {
        "sigma": "sigma",
        "s @ N^(2/3)": "s @ N^(2/3)",
        "s @ sqrt*log": "s @ sqrt*log",
        "s @ sqrt(N)": "s @ sqrt(N)",
        "s @ sqrt/log": "s @ sqrt/log",
    }
    print(f"\n  Result table:")
    print(f"    {'observable':>14} {'Tao_pred':>10} "
          f"{'slope_raw':>10} {'gap_raw':>9} "
          f"{'slope_t1%':>10} {'gap_t1%':>9}")
    rows = []
    for obs_label, obs_col in observables:
        delta_log = delta_logs[delta_keys[obs_label]]
        tao = K_H * delta_log
        slope_r, off_r = fit_slope_offset(alpha_det_h, means[obs_col])
        slope_t, off_t = fit_slope_offset(alpha_det_h, t1s[obs_col])
        print(f"    {obs_label:>14} {tao:>10.3f} "
              f"{slope_r:>10.4f} {off_r-tao:>+9.3f} "
              f"{slope_t:>10.4f} {off_t-tao:>+9.3f}")
        rows.append({
            "observable": obs_label,
            "tao_pred": tao,
            "slope_raw": slope_r,
            "gap_raw": off_r - tao,
            "slope_t1": slope_t,
            "gap_t1": off_t - tao,
        })

    out_dir = here.parent / "experiments_output"
    out_csv = out_dir / "37_alpha_det_k16_verification_N134217728.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n  [save] {out_csv}")


if __name__ == "__main__":
    main()
