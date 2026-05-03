"""
Experiment 30 — first-passage time below sqrt(N), per-class alpha_det correlation.

For each odd N <= 2^27, walk Collatz steps and record s(N) = step count at
which the orbit first attains a value <= sqrt(N). Group by residue class
r = N mod 256 (k=8, 128 odd classes). Compute s_median(r). Compute alpha_det(r)
from the deterministic prefix at k=8. Compute Spearman correlation between
alpha_det and s_median, and partial Spearman correlation removing per-class
mean log(N).

Hypothesis: alpha_det captures prefix-driven early-descent structure. If so,
s_median (first-passage below sqrt(N), dominated by prefix dynamics) should
correlate with alpha_det beyond what global log(N) explains.

Usage:
    python 30_first_passage_alpha_det.py --N 134217728
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
def first_passage_below_sqrt(starts):
    """For each odd start N, walk Collatz until value <= sqrt(N).
    Returns step count s(N) per N."""
    n = len(starts)
    s_arr = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        N = starts[i]
        sqrt_N = np.int64(np.floor(np.sqrt(np.float64(N))))
        m = N
        steps = 0
        # If N itself <= sqrt(N), we'd never enter; but for odd N>=3, sqrt(N) >= 1.7
        # and N > sqrt(N) always for N > 1. So m starts > sqrt(N).
        while m > sqrt_N:
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            steps += 1
            if steps > 10000:
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
    n = df["n"].to_numpy().astype(np.int64)
    print(f"        odd rows: {len(n):,}  loaded in {time.perf_counter()-t0:.1f}s",
          flush=True)

    # First-passage walks
    print(f"[run]   computing first-passage step counts ...", flush=True)
    t0 = time.perf_counter()
    s = first_passage_below_sqrt(n)
    print(f"[run]   done in {time.perf_counter()-t0:.1f}s   "
          f"s range: [{int(s.min())}, {int(s.max())}]", flush=True)

    # Group by N mod 256
    M = 256
    K = M // 2
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    # Per-class median s and mean log N
    s_median_per = np.zeros(K, dtype=np.float64)
    log_n_mean_per = np.zeros(K, dtype=np.float64)
    n_per_class = np.zeros(K, dtype=np.int64)

    log_n = np.log(n.astype(np.float64))
    for kk in range(K):
        mask = class_idx == kk
        if mask.sum() < 50:
            continue
        s_median_per[kk] = float(np.median(s[mask]))
        log_n_mean_per[kk] = float(log_n[mask].mean())
        n_per_class[kk] = int(mask.sum())

    print(f"[group] {K} odd classes mod {M}; per-class n: "
          f"[{n_per_class.min():,}, {n_per_class.max():,}]", flush=True)

    # Compute alpha_det per class
    alpha_det = np.zeros(K)
    a_final_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r = 2 * kk + 1
        steps, a_f, _ = deterministic_prefix(r, M)
        alpha_det[kk] = steps + LOG_FACTOR_ODD * np.log(a_f / float(M))
        a_final_per[kk] = a_f

    # Spearman correlations
    rho_raw, p_raw = spearmanr(alpha_det, s_median_per)
    rho_logn, p_logn = spearmanr(log_n_mean_per, s_median_per)
    rho_alpha_logn, p_alpha_logn = spearmanr(alpha_det, log_n_mean_per)

    # Partial Spearman: regress alpha_det and s_median on log_n_mean,
    # then correlate the residuals.
    def rank(x):
        return np.argsort(np.argsort(x)).astype(np.float64)

    r_alpha = rank(alpha_det)
    r_s = rank(s_median_per)
    r_logn = rank(log_n_mean_per)

    def lin_resid(y, x):
        x_c = x - x.mean()
        y_c = y - y.mean()
        slope = (x_c * y_c).sum() / (x_c * x_c).sum() if (x_c * x_c).sum() > 0 else 0.0
        return y_c - slope * x_c

    e_alpha = lin_resid(r_alpha, r_logn)
    e_s = lin_resid(r_s, r_logn)
    rho_partial, p_partial = spearmanr(e_alpha, e_s)

    print()
    print(f"=== Results (k=8, 128 odd classes mod 256, N <= 2^{int(round(np.log2(args.N)))}) ===")
    print(f"  Spearman rho(alpha_det, s_median)            = {rho_raw:>+.4f}    p = {p_raw:.3e}")
    print(f"  Spearman rho(<log N>, s_median)              = {rho_logn:>+.4f}    p = {p_logn:.3e}")
    print(f"  Spearman rho(alpha_det, <log N>)             = {rho_alpha_logn:>+.4f}    p = {p_alpha_logn:.3e}")
    print(f"  Partial Spearman(alpha_det, s_median | logN) = {rho_partial:>+.4f}    p = {p_partial:.3e}")

    # Ranges & sanity
    print()
    print(f"=== Ranges ===")
    print(f"  alpha_det:    [{alpha_det.min():.3f}, {alpha_det.max():.3f}]   "
          f"SD = {alpha_det.std():.3f}")
    print(f"  s_median:     [{s_median_per.min():.1f}, {s_median_per.max():.1f}]   "
          f"SD = {s_median_per.std():.3f}")
    print(f"  <log N>:      [{log_n_mean_per.min():.3f}, {log_n_mean_per.max():.3f}]   "
          f"SD = {log_n_mean_per.std():.3f}")

    # Per a_final: pool s_median across classes sharing a_final
    print()
    print(f"=== s_median by a_final group ===")
    a_final_unique = sorted(set(a_final_per.tolist()))
    print(f"  {'a_final':>9} {'classes':>9} {'mean_alpha_det':>15} {'mean_s_median':>14}")
    for af in a_final_unique:
        mask = (a_final_per == af) & (n_per_class > 0)
        if mask.sum() == 0:
            continue
        print(f"  {af:>9} {int(mask.sum()):>9} "
              f"{alpha_det[mask].mean():>15.4f} "
              f"{s_median_per[mask].mean():>14.2f}")

    # Save
    out = pl.DataFrame({
        "class_idx": np.arange(K),
        "residue": np.arange(1, M, 2),
        "a_final": a_final_per,
        "alpha_det": alpha_det,
        "s_median": s_median_per,
        "log_n_mean": log_n_mean_per,
        "n_per_class": n_per_class,
    })
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"30_first_passage_alpha_det_N{args.N}.csv"
    out.write_csv(out_path)
    print()
    print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
