"""
Experiment 37 (TA.2) — trim-quantile sweep, find trim* where Tao gap = 0.

For each trim quantile q in {0.05%, 0.1%, 0.2%, 0.3%, 0.5%, 0.75%, 1.0%,
1.5%, 2.0%}, compute the per-class mean s for the s @ sqrt(N) observable
with the top-q quantile of values trimmed within each class. Regress
against alpha_det at K_h, report offset, K_implied, and gap from Tao
prediction K_h * log(N)/2.

Identify the trim q* at which gap = 0 (within +/- 0.05). Compare to
log^(-c) N for various c values; if q* matches log^(-c) N for some c,
that's the empirical exceptional-set density at this N.

Done at N = 2^27 (existing parquet) and k = 8.

Usage:
    python 37_TA2_trim_quantile_sweep.py
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
def first_passage_sqrt(starts):
    """First-passage to sqrt(N) per start. Returns step counts."""
    n = len(starts)
    s_arr = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        N = starts[i]
        sqrt_N = np.int64(np.floor(np.sqrt(np.float64(N))))
        m = N
        steps = 0
        while m > sqrt_N:
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            steps += 1
            if steps > 20000:
                break
        s_arr[i] = steps
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
    K_cls = M // 2
    alpha = np.zeros(K_cls)
    for kk in range(K_cls):
        r = 2 * kk + 1
        sp, a_f, _ = deterministic_prefix(r, M)
        alpha[kk] = sp + K_H * np.log(a_f / float(M))
    return alpha


def fit_slope_offset(x, y):
    x_c = x - x.mean()
    y_c = y - y.mean()
    slope = (x_c * y_c).sum() / (x_c * x_c).sum()
    offset = y.mean() - slope * x.mean()
    pred = slope * x + offset
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = (y_c ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(offset), float(r2)


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    N = 1 << 27
    print(f"[load] N=2^27 = {N:,}", flush=True)
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    print(f"  {len(n_arr):,} odd rows in {time.perf_counter()-t0:.1f}s", flush=True)

    print(f"[run] walking trajectories (first-passage to sqrt(N))...", flush=True)
    t0 = time.perf_counter()
    s_arr = first_passage_sqrt(n_arr)
    print(f"  done in {time.perf_counter()-t0:.1f}s", flush=True)

    log_n = np.log(n_arr.astype(np.float64))
    log_n_mean_global = float(log_n.mean())
    log_log_N = float(np.log(log_n_mean_global))
    print(f"  <log N> = {log_n_mean_global:.4f}    log log N = {log_log_N:.4f}", flush=True)

    # k=8 groupby
    M = 256
    K_cls = M // 2
    res = (n_arr % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    # Build polars frame for trimmed-quantile aggregations
    df_g = pl.DataFrame({
        "class_idx": class_idx,
        "log_n": log_n,
        "s": s_arr,
    })

    quantiles = [0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05]

    # Aggregate: mean log_n + raw mean s + trimmed means
    aggs = [
        pl.col("log_n").mean().alias("log_n_mean"),
        pl.col("s").mean().alias("s_mean_raw"),
    ]
    for q in quantiles:
        threshold_q = 1.0 - q
        col_alias = f"s_t_{int(q*1e5):d}"  # encode q in basis points (0.05% = 5)
        aggs.append(
            pl.col("s").filter(pl.col("s") < pl.col("s").quantile(threshold_q)).mean().alias(col_alias)
        )
    print(f"[agg] computing per-class means + {len(quantiles)} trim levels ...", flush=True)
    t0 = time.perf_counter()
    agg = df_g.group_by("class_idx").agg(aggs).sort("class_idx")
    print(f"  done in {time.perf_counter()-t0:.1f}s", flush=True)
    del df_g
    gc.collect()

    ci = agg["class_idx"].to_numpy()
    log_n_per = np.zeros(K_cls)
    log_n_per[ci] = agg["log_n_mean"].to_numpy()
    s_raw = np.zeros(K_cls)
    s_raw[ci] = agg["s_mean_raw"].to_numpy()
    s_trim = {}
    for q in quantiles:
        s_arr_q = np.zeros(K_cls)
        s_arr_q[ci] = agg[f"s_t_{int(q*1e5):d}"].to_numpy()
        s_trim[q] = s_arr_q

    alpha_det = compute_alpha_det(8)
    delta_log = log_n_mean_global / 2.0
    tao_pred = K_H * delta_log
    print(f"  Tao prediction K_h * log(N)/2 = {tao_pred:.4f}", flush=True)

    print()
    print(f"=== Trim quantile sweep at N=2^27, k=8, observable = s @ sqrt(N) ===")
    print(f"  {'trim_q':>10} {'slope':>9} {'offset':>10} {'gap':>10}")
    print(f"  {'(none)':>10} ", end="")
    slope_r, offset_r, _ = fit_slope_offset(alpha_det, s_raw)
    print(f"{slope_r:>9.4f} {offset_r:>+10.4f} {offset_r - tao_pred:>+10.4f}")
    rows = [(0.0, slope_r, offset_r, offset_r - tao_pred)]
    for q in quantiles:
        slope_q, offset_q, _ = fit_slope_offset(alpha_det, s_trim[q])
        gap = offset_q - tao_pred
        print(f"  {100*q:>9.3f}% {slope_q:>9.4f} {offset_q:>+10.4f} {gap:>+10.4f}")
        rows.append((q, slope_q, offset_q, gap))

    # Find trim* where gap = 0 by interpolation
    qs = np.array([r[0] for r in rows])
    gaps = np.array([r[3] for r in rows])
    # gap is a decreasing function of q (more trim -> lower mean)
    # Find the first q where gap crosses 0
    crossings = []
    for i in range(len(qs) - 1):
        if (gaps[i] > 0) and (gaps[i + 1] < 0):
            t = gaps[i] / (gaps[i] - gaps[i + 1])
            q_star = qs[i] + t * (qs[i + 1] - qs[i])
            crossings.append((qs[i], qs[i + 1], q_star))
        elif (gaps[i] < 0) and (gaps[i + 1] > 0):
            t = -gaps[i] / (gaps[i + 1] - gaps[i])
            q_star = qs[i] + t * (qs[i + 1] - qs[i])
            crossings.append((qs[i], qs[i + 1], q_star))

    print()
    if crossings:
        for lo, hi, q_star in crossings:
            print(f"  Zero-gap crossing: q* ≈ {100*q_star:.4f}%   "
                  f"(between {100*lo:.3f}% and {100*hi:.3f}%)")
        q_star = crossings[0][2]
        # Compare to log^(-c) N for various c
        print()
        print(f"  log^(-c) N for various c (at N=2^27, log N = {log_n_mean_global:.4f}):")
        for c in [0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0]:
            pred = log_n_mean_global ** (-c)
            print(f"    c={c}: log^(-c) N = {pred:.6f}   ratio q*/log^(-c)N = "
                  f"{q_star/pred:.4f}")
        # Solve for c that makes log^(-c) N = q*
        c_match = -np.log(q_star) / np.log(log_n_mean_global)
        print(f"  c such that log^(-c) N = q*:  c = {c_match:.4f}")
        print(f"  (Tao's exceptional-set bound is log^(-c) for some c > 0; "
              f"check if c_match has structural meaning.)")
    else:
        print("  No zero-gap crossing in the swept range.")

    # Save
    df_out = pl.DataFrame(
        [(q, sl, off, g) for q, sl, off, g in rows],
        schema=["trim_q", "slope", "offset", "gap_from_tao"],
        orient="row",
    )
    out = out_dir / "37_TA2_trim_quantile_sweep_N134217728_k8.csv"
    df_out.write_csv(out)
    print()
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
