"""
sweep_k.py - Test +12.20 per-class shift invariance from k=6 up to k=16.

Reuses descent_b_enlarged.csv (200K orbits, n_max=10^6).
For each k, builds a* table, runs joint OLS, reports the per-class slope.
"""
import importlib.util
import math
import time
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
EXP29 = ROOT / "experiments" / "29_qx1_cycle_classification.py"
spec = importlib.util.spec_from_file_location("exp29", EXP29)
exp29 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp29)
qx1_prefix = exp29.qx1_prefix


def a_star_idx(r, k, q=3):
    if r == 0:
        return 0
    _, a_star, _ = qx1_prefix(r, k, q)
    j = 0
    a = a_star
    while a > 1 and a % q == 0:
        a //= q
        j += 1
    return j if a == 1 else -1


def precompute_a_star_table(k):
    table = np.zeros(1 << k, dtype=np.int8)
    for r in range(1, 1 << k, 2):
        table[r] = a_star_idx(r, k)
    return table


def chase_k(df, k, ln3=math.log(3)):
    table = precompute_a_star_table(k)
    n_arr = df["n"].to_numpy()
    j_arr = table[n_arr % (1 << k)]
    df_k = df.with_columns(pl.Series("j", j_arr.astype(np.int64)))

    # Joint OLS: sigma = a + b*log(n) + c*j
    X = np.column_stack([
        np.ones(len(df_k)),
        df_k["log_n"].to_numpy(),
        df_k["j"].to_numpy().astype(np.float64),
    ])
    y = df_k["sigma"].to_numpy().astype(np.float64)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n = len(y)
    r2 = 1 - np.sum(resid**2) / np.sum((y - y.mean())**2)
    sigma2 = np.sum(resid**2) / (n - 3)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se_b = math.sqrt(cov[1, 1])
    se_c = math.sqrt(cov[2, 2])

    # Class population spread
    class_counts = df_k.group_by("j").agg(pl.len().alias("n_orbits")).sort("j")
    n_min = class_counts["n_orbits"].min()
    n_max_class = class_counts["n_orbits"].max()
    n_classes_used = len(class_counts)

    # Heuristic prediction with corrected constant
    log_n_slope_predicted = beta[1]  # use observed instead of K=10.43 dogma
    j_shift_predicted = 1 + log_n_slope_predicted * ln3

    return dict(
        k=k,
        beta_intercept=beta[0],
        beta_logn=beta[1],
        beta_j=beta[2],
        se_logn=se_b,
        se_j=se_c,
        r2=r2,
        n_classes_used=n_classes_used,
        n_min=n_min,
        n_max=n_max_class,
        j_shift_predicted=j_shift_predicted,
        gap_obs_minus_pred=beta[2] - j_shift_predicted,
    )


def main():
    df = pl.read_csv(Path(__file__).parent / "viz_outputs" / "descent_b_enlarged.csv")
    print(f"Loaded {len(df):,} orbits. n_max={df['n'].max():,}, log(n_max)={math.log(df['n'].max()):.2f}")
    print()

    results = []
    for k in [6, 8, 10, 12, 14, 16]:
        print(f"  building a* table for k={k} ({1<<k} residues) ...", flush=True)
        t0 = time.perf_counter()
        r = chase_k(df, k)
        elapsed = time.perf_counter() - t0
        print(f"    done in {elapsed:.1f}s", flush=True)
        results.append(r)

    print(f"\n{'='*100}")
    print(f"{'k':>3} {'b_logn':>8} {'b_j(obs)':>10} {'SE_j':>7} {'pred_j':>8} {'gap':>7} "
          f"{'R^2':>7} {'n_cls':>6} {'n_min':>7} {'n_max':>8}")
    print(f"{'-'*100}")
    for r in results:
        print(f"{r['k']:>3} {r['beta_logn']:>8.3f} {r['beta_j']:>+10.4f} {r['se_j']:>7.4f} "
              f"{r['j_shift_predicted']:>8.4f} {r['gap_obs_minus_pred']:>+7.4f} "
              f"{r['r2']:>7.4f} {r['n_classes_used']:>6d} {r['n_min']:>7d} {r['n_max']:>8d}")
    print(f"{'='*100}")
    print(f"  pred_j = 1 + observed_b_logn * ln(3)   (heuristic with sample-fitted slope)")
    print(f"  gap = observed - predicted (negative means heuristic overpredicts)")


if __name__ == "__main__":
    main()
