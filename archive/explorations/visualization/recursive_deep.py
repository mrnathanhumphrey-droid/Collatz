"""
recursive_deep.py - Push recursive a* decomposition to depth 1..8 at k=6.

For each orbit, iteratively decompose: n -> v1 -> v2 -> ... -> v_D
Each level d gives an a*_idx j_d. Joint OLS at each depth D:
  sigma ~ log(n) + j_1 + j_2 + ... + j_D

Report:
  - per-level coefficients (do they all stay near +12?)
  - R^2 vs depth (where does it saturate?)
  - validity rate at each depth (when does v_d become unusable?)

Also test mixed-depth (k1, k2) configurations to verify level-k invariance.
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


def precompute_table(k, q=3):
    """For each odd r in [1, 2^k), return (j, a_star, c_star)."""
    j_arr = np.zeros(1 << k, dtype=np.int8)
    a_arr = np.zeros(1 << k, dtype=np.int64)
    c_arr = np.zeros(1 << k, dtype=np.int64)
    for r in range(1, 1 << k, 2):
        steps, a_star, c_star = qx1_prefix(r, k, q)
        j = 0
        a = a_star
        while a > 1 and a % q == 0:
            a //= q
            j += 1
        j_arr[r] = j if a == 1 else -1
        a_arr[r] = a_star
        c_arr[r] = c_star
    return j_arr, a_arr, c_arr


def reduce_to_odd(v):
    """Strip factors of 2 from each element of v (in place semantics)."""
    v_odd = v.copy()
    while True:
        even_mask = (v_odd > 0) & (v_odd % 2 == 0)
        if not even_mask.any():
            break
        v_odd[even_mask] //= 2
    return v_odd


def decompose_one_level(v, j_table, a_table, c_table, k):
    """v -> (j_at_v, v_next_odd, valid). v assumed odd >= 2^k."""
    r = v % (1 << k)
    m = v // (1 << k)
    j = j_table[r]
    a = a_table[r]
    c = c_table[r]
    v_next = a * m + c
    v_next_odd = reduce_to_odd(v_next)
    valid = v_next_odd >= (1 << k)
    return j, v_next_odd, valid


def joint_ols(features, target):
    X = np.column_stack([np.ones(len(target))] + list(features))
    y = target
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n = len(y); p = X.shape[1]
    r2 = 1 - np.sum(resid**2) / np.sum((y - y.mean())**2)
    sigma2 = np.sum(resid**2) / (n - p)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return beta, r2, se


def main():
    df = pl.read_csv(Path(__file__).parent / "viz_outputs" / "descent_b_enlarged.csv")
    n_arr = df["n"].to_numpy()
    log_n = df["log_n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)
    print(f"Loaded {len(df):,} orbits.")

    k = 6
    j_table, a_table, c_table = precompute_table(k)

    print(f"\n{'='*80}\nDEPTH SWEEP at k={k}\n{'='*80}")

    # Iteratively build up j arrays per level
    v = n_arr.copy()
    valid = np.ones(len(n_arr), dtype=bool)
    j_levels = []  # list of arrays, one per depth

    print(f"  {'depth':>5} {'n_valid':>9} {'j_d coef':>10} {'R^2':>8} {'log_n_coef':>11}")
    for depth in range(1, 9):
        if depth > 1:
            # Compute next-level v from prior v
            jd, v_next, valid_next = decompose_one_level(v, j_table, a_table, c_table, k)
            j_levels.append(jd)
            v = v_next
            valid = valid & valid_next
        else:
            jd, v_next, valid_next = decompose_one_level(v, j_table, a_table, c_table, k)
            j_levels.append(jd)
            v = v_next
            valid = valid_next

        # Run joint OLS over the current valid set with all j_levels so far
        n_valid = int(valid.sum())
        if n_valid < 1000:
            print(f"  {depth:>5} {n_valid:>9} (insufficient orbits, stopping)")
            break
        log_n_v = log_n[valid]
        sigma_v = sigma[valid]
        feats = [log_n_v]
        for jl in j_levels:
            feats.append(jl[valid].astype(np.float64))
        beta, r2, se = joint_ols(feats, sigma_v)
        # beta = [intercept, log_n_coef, j1_coef, j2_coef, ..., jD_coef]
        log_n_coef = beta[1]
        last_j_coef = beta[-1]
        coef_str = " ".join(f"{c:+.2f}" for c in beta[2:])
        print(f"  {depth:>5} {n_valid:>9} {last_j_coef:>+10.4f} {r2:>8.4f} {log_n_coef:>11.4f}  all_j_coefs=[{coef_str}]")

    # Mixed-depth tests
    print(f"\n{'='*80}\nMIXED-DEPTH TESTS\n{'='*80}")
    for (k1, k2) in [(6, 10), (10, 6), (6, 14), (14, 6), (8, 8)]:
        t1 = precompute_table(k1)
        t2 = precompute_table(k2)
        j1, v1_odd, valid1 = decompose_one_level(n_arr, *t1, k1)
        if not valid1.any():
            print(f"  (k1={k1}, k2={k2}): no valid orbits at level 1, skipping")
            continue
        j2, v2_odd, valid2 = decompose_one_level(v1_odd, *t2, k2)
        valid_both = valid1 & valid2
        n_v = int(valid_both.sum())
        if n_v < 1000:
            print(f"  (k1={k1}, k2={k2}): only {n_v} valid orbits, skipping")
            continue
        log_n_v = log_n[valid_both]
        sigma_v = sigma[valid_both]
        j1f = j1[valid_both].astype(np.float64)
        j2f = j2[valid_both].astype(np.float64)
        beta, r2, se = joint_ols([log_n_v, j1f, j2f], sigma_v)
        print(f"  (k1={k1}, k2={k2}): n_valid={n_v:,}  j1_coef={beta[2]:+.3f} (SE {se[2]:.3f})  "
              f"j2_coef={beta[3]:+.3f} (SE {se[3]:.3f})  R^2={r2:.4f}")


if __name__ == "__main__":
    main()
