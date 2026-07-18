"""
chase_deviations.py - Where does the +12.02 heuristic prediction break?

For each (k, j) class:
  predicted sigma_mean(j) = a + b*log(n_mean) + 12.02*j   (b=10.03 from theory)
                           ≈ joint-OLS prediction
  observed sigma_mean(j) = empirical mean

Compute residuals (obs - pred) and SE per class. Flag any class where
|residual| > 2*SE as a candidate deviation worth following.

Also re-walks an enlarged sample (200K orbits) so the smallest classes at
k=10 (94 orbits at 50K) get enough resolution.
"""
import csv
import importlib.util
import math
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
EXP29 = ROOT / "experiments" / "29_qx1_cycle_classification.py"
spec = importlib.util.spec_from_file_location("exp29", EXP29)
exp29 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp29)
qx1_prefix = exp29.qx1_prefix

OUT_DIR = Path(__file__).parent / "viz_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
ENLARGED_CSV = OUT_DIR / "descent_b_enlarged.csv"

N_SAMPLES_ENLARGED = 200_000
N_MAX = 1_000_000
MAX_STEPS = 5_000
SEED = 42


@njit(cache=True)
def sigma_3x1(n_start, max_steps):
    n = n_start
    for step in range(max_steps):
        if n == 1:
            return step
        if n % 2 == 0:
            n = n >> 1
        else:
            n = 3 * n + 1
    return -1


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


def generate_enlarged():
    """200K orbits, sigma + log(n) only."""
    if ENLARGED_CSV.exists():
        print(f"[gen] using cached {ENLARGED_CSV}")
        return
    print(f"[gen] generating {N_SAMPLES_ENLARGED:,} orbits ...")
    rng = np.random.default_rng(SEED)
    odd_pool = np.arange(3, N_MAX + 1, 2, dtype=np.int64)
    starts = rng.choice(odd_pool, size=N_SAMPLES_ENLARGED, replace=False)
    starts.sort()

    # JIT warmup
    _ = sigma_3x1(np.int64(7), MAX_STEPS)
    t0 = time.perf_counter()
    sigmas = np.zeros(N_SAMPLES_ENLARGED, dtype=np.int64)
    for i, n in enumerate(starts):
        sigmas[i] = sigma_3x1(np.int64(n), MAX_STEPS)
    print(f"[gen] walk done in {time.perf_counter()-t0:.1f}s")

    df = pl.DataFrame({"n": starts, "log_n": np.log(starts), "sigma": sigmas})
    df.write_csv(ENLARGED_CSV)
    print(f"[gen] wrote {ENLARGED_CSV}  ({len(df):,} rows)")


def chase(k, df_enlarged):
    print(f"\n{'='*70}\nk = {k}    chasing per-class deviations from heuristic +12.02\n{'='*70}")

    table = precompute_a_star_table(k)
    n_arr = df_enlarged["n"].to_numpy()
    r_at_k = n_arr % (1 << k)
    a_idx = table[r_at_k]

    df_k = df_enlarged.with_columns([pl.Series("j", a_idx)])

    # Joint OLS to get baseline regression
    X = np.column_stack([
        np.ones(len(df_k)),
        df_k["log_n"].to_numpy(),
        df_k["j"].to_numpy().astype(np.float64),
    ])
    y = df_k["sigma"].to_numpy().astype(np.float64)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    sigma_resid2 = np.sum((y - yhat) ** 2) / (len(y) - 3)

    print(f"  Joint OLS:  sigma = {beta[0]:.3f} + {beta[1]:.3f}*log(n) + {beta[2]:.3f}*j")
    print(f"  Per-class shift (data) = {beta[2]:+.3f}")
    print(f"  Heuristic prediction   = {1 + 6.95 * math.log2(3):+.3f}  (= 1 + 6.95*log2(3))")
    print(f"  Within-orbit residual std = {math.sqrt(sigma_resid2):.2f} steps")

    # Per-class deviations from joint OLS prediction
    rows = []
    for j in range(1, k + 2):
        sub = df_k.filter(pl.col("j") == j)
        if len(sub) == 0:
            continue
        n_obs = len(sub)
        sigma_obs = sub["sigma"].mean()
        logn_mean = sub["log_n"].mean()
        sigma_pred = beta[0] + beta[1] * logn_mean + beta[2] * j
        dev = sigma_obs - sigma_pred
        sigma_within = sub["sigma"].std()
        se = sigma_within / math.sqrt(n_obs)
        z = dev / se if se > 0 else 0
        rows.append((j, n_obs, sigma_obs, sigma_pred, dev, se, z))

    print(f"\n  Per-class (k={k}) deviations from joint-OLS prediction:")
    print(f"    {'j':>3} {'n':>7} {'sigma_obs':>10} {'sigma_pred':>11} {'dev':>8} {'SE':>7} {'z':>7}")
    for j, n_obs, so, sp, dv, se, z in rows:
        flag = " ***" if abs(z) > 2.5 else (" *" if abs(z) > 1.96 else "")
        print(f"    {j:>3d} {n_obs:>7d} {so:>10.3f} {sp:>11.3f} {dv:>+8.3f} {se:>7.3f} {z:>+7.2f}{flag}")

    # Extra: the cumulative per-class shift (obs - obs_at_j_minus_1)
    print(f"\n  Adjacent-class shifts (obs sigma_mean(j+1) - obs sigma_mean(j)):")
    print(f"    Heuristic predicts each step = +12.02")
    sigma_means = [r[2] for r in rows]
    js = [r[0] for r in rows]
    for i in range(len(rows) - 1):
        diff = sigma_means[i+1] - sigma_means[i]
        print(f"    j={js[i]}->{js[i+1]}:  +{diff:>6.3f}  (dev from 12.02: {diff-12.02:+.3f})")


def main():
    generate_enlarged()
    df = pl.read_csv(ENLARGED_CSV)
    print(f"\nLoaded {len(df):,} orbits.")

    for k in [6, 8, 10]:
        chase(k, df)


if __name__ == "__main__":
    main()
