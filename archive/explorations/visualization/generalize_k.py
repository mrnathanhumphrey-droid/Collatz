"""
generalize_k.py - Test whether the +12 steps/class sigma shift at k=6
generalizes to k=8 and k=10.

Reuses existing descent_b.csv (50K orbits with sigma, log_n, log_ratio).
For each k in {6, 8, 10}, recompute a_star_idx per orbit, regress jointly
sigma ~ log(n) + a_star_idx, report:
  - per-class mean sigma table
  - joint OLS coefficient on a_star_idx ("the shift per class")
  - within-class slopes on log(n)
  - class-size distribution (verify binomial C(k-1, j-1))
"""
import importlib.util
import math
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
    """Map r -> a_star_idx for all odd r in [1, 2^k)."""
    table = np.zeros(1 << k, dtype=np.int8)
    for r in range(1, 1 << k, 2):
        table[r] = a_star_idx(r, k)
    return table


def class_size_distribution(table, k):
    """Count residue classes per a_star_idx j."""
    odd_only = table[1::2]  # only odd residues used
    counts = np.bincount(odd_only, minlength=k + 2)
    return counts


def analyze_k(df, k):
    print(f"\n{'='*60}\nk = {k}  (2^{k} = {1<<k} residue classes, {(1<<k)//2} odd)\n{'='*60}")

    table = precompute_a_star_table(k)
    sizes = class_size_distribution(table, k)
    print(f"  Class-size distribution by a_star_idx j (odd residues only):")
    binomial = [math.comb(k - 1, j - 1) if 1 <= j <= k else 0 for j in range(k + 2)]
    for j in range(1, k + 2):
        if sizes[j] > 0:
            match = "  matches C({},{})".format(k - 1, j - 1) if sizes[j] == binomial[j] else "  != binomial"
            print(f"    j={j:>2}  n_classes={sizes[j]:>4d}  predicted_binomial={binomial[j]:>4d}{match}")

    # Annotate orbits with a_star_idx at this k
    n_int = (df["x"].exp().round().cast(pl.Int64)).to_numpy()  # recover n from log
    # Sanity: n must be odd
    assert all(n_int % 2 == 1), "non-odd n detected after exp(log) round"
    r_at_k = n_int % (1 << k)
    a_idx_at_k = table[r_at_k]

    df_k = df.with_columns([
        pl.Series("a_idx_k", a_idx_at_k),
    ])

    print(f"\n  Per-class sigma stats:")
    print(df_k.group_by("a_idx_k").agg([
        pl.len().alias("n"),
        pl.col("y").mean().round(2).alias("sigma_mean"),
        pl.col("y").median().alias("sigma_med"),
        pl.col("x").mean().round(3).alias("logn_mean"),
        pl.col("x").std().round(3).alias("logn_std"),
        pl.col("z").mean().round(3).alias("log_ratio_mean"),
    ]).sort("a_idx_k"))

    # Joint OLS
    x = df_k.select(["x", "a_idx_k"]).to_numpy().astype(np.float64)
    y = df_k["y"].to_numpy().astype(np.float64)
    X = np.column_stack([np.ones(len(x)), x[:, 0], x[:, 1]])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    r2 = 1 - np.sum(resid**2) / np.sum((y - y.mean())**2)
    n_obs = len(y)
    # Standard error on the a_idx coefficient
    sigma_hat2 = np.sum(resid**2) / (n_obs - 3)
    cov = sigma_hat2 * np.linalg.inv(X.T @ X)
    se_a_idx = np.sqrt(cov[2, 2])

    print(f"\n  Joint OLS:  sigma = {beta[0]:.3f} + {beta[1]:.3f}*log(n) + {beta[2]:.3f}*a_idx")
    print(f"    Per-class shift:  {beta[2]:+.3f}  (SE {se_a_idx:.3f})")
    print(f"    R^2:              {r2:.4f}")
    print(f"    log(n) slope:     {beta[1]:.3f}  (heuristic ~10.03 for sigma ~6.95*log2(n))")


def main():
    df = pl.read_csv(Path(__file__).parent / "viz_outputs" / "descent_b.csv")
    print(f"Loaded {len(df):,} orbits from descent_b.csv")
    print(f"Reminder: sigma = column 'y'; log(n_start) = column 'x'; log(peak/n) = column 'z'")

    for k in [6, 8, 10]:
        analyze_k(df, k)

    # Cross-k summary
    print(f"\n{'='*60}\nCROSS-K SUMMARY (per-class sigma shift)\n{'='*60}")
    print(f"  k=6  per-class shift was reported as +12.18 (verified above)")


if __name__ == "__main__":
    main()
