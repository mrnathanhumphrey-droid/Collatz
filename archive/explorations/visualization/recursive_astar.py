"""
recursive_astar.py - Test recursive a* decomposition.

For each orbit, after the first prefix at level k1, the post-prefix value is
  v1 = a*_1 * m1 + c*_1   where m1 = n // 2^k1, r1 = n % 2^k1
That v1 has its own residue mod 2^k2, so its own a*_2 = 3^j2.

Hypothesis (from heuristic): joint OLS sigma ~ log(n) + j1 + j2 should give
both j1 and j2 coefficients near +12.20 (additive in the recursion).

Compare to direct single-level OLS at k = k1 + k2 to see whether the factored
(j1, j2) classification explains as much / more variance than the flat-prefix
classification at the combined level.
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


def a_star_full(r, k, q=3):
    """Return (j, a_star, c_star) where a_star = q^j."""
    if r == 0:
        return 0, 1, 0
    _, a_star, c_star = qx1_prefix(r, k, q)
    j = 0
    a = a_star
    while a > 1 and a % q == 0:
        a //= q
        j += 1
    return (j if a == 1 else -1), a_star, c_star


def precompute_table(k):
    """For each odd r in [1, 2^k), store (j, a_star, c_star) at level k."""
    j_arr = np.zeros(1 << k, dtype=np.int8)
    a_arr = np.zeros(1 << k, dtype=np.int64)
    c_arr = np.zeros(1 << k, dtype=np.int64)
    for r in range(1, 1 << k, 2):
        j, a, c = a_star_full(r, k)
        j_arr[r] = j
        a_arr[r] = a
        c_arr[r] = c
    return j_arr, a_arr, c_arr


def recursive_decompose(n_arr, k1, k2):
    """For each n, compute (j1, v1, j2). Returns three arrays.
    v1 may be even — for j2 we use v1 mod 2^k2 directly (decomposes any value)."""
    j1_table, a1_table, c1_table = precompute_table(k1)
    # We need j2 at level k2 for arbitrary v1 — including even v1's residue mod 2^k2.
    # The qx1_prefix function expects ODD residues. If v1 is even, we have to first
    # reduce v1 to its odd part by stripping factors of 2 (matching what Collatz would do).
    # For the recursive decomposition to be meaningful, we follow actual Collatz from v1
    # to the next odd value, then decompose that odd value's residue mod 2^k2.
    j2_table_full, _, _ = precompute_table(k2)

    n = n_arr
    r1 = n % (1 << k1)
    m1 = n // (1 << k1)
    j1 = j1_table[r1]
    a1 = a1_table[r1]
    c1 = c1_table[r1]
    v1 = a1 * m1 + c1

    # Reduce v1 to its odd part (mimicking the halve-only Collatz steps before next odd)
    v1_odd = v1.copy()
    even_mask = (v1_odd > 0) & (v1_odd % 2 == 0)
    while even_mask.any():
        v1_odd[even_mask] //= 2
        even_mask = (v1_odd > 0) & (v1_odd % 2 == 0)

    # Compute j2 from v1_odd's residue mod 2^k2
    valid = v1_odd >= (1 << k2)  # only orbits where second prefix isn't trivially exhausted
    j2 = np.zeros_like(j1)
    r2 = v1_odd % (1 << k2)
    j2[valid] = j2_table_full[r2[valid]]

    return j1, v1, v1_odd, j2, valid


def joint_ols(features, target):
    """Returns beta, R^2, SE on each coefficient."""
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
    print(f"Loaded {len(df):,} orbits.")
    n_arr = df["n"].to_numpy()
    log_n = df["log_n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)

    # Single-level baselines (recompute for sanity)
    for k_single in [6, 12]:
        j1_table, _, _ = precompute_table(k_single)
        j_single = j1_table[n_arr % (1 << k_single)].astype(np.float64)
        beta, r2, se = joint_ols([log_n, j_single], sigma)
        print(f"\n[baseline] single-level k={k_single}:")
        print(f"  sigma = {beta[0]:.3f} + {beta[1]:.3f}*log(n) + {beta[2]:.3f}*j   (SE_j={se[2]:.3f})")
        print(f"  R^2 = {r2:.4f}")

    # Recursive (j1 at k1=6, j2 at k2=6)
    print(f"\n{'='*70}\nRecursive decomposition (k1=6, k2=6)\n{'='*70}")
    j1, v1, v1_odd, j2, valid = recursive_decompose(n_arr, 6, 6)
    n_valid = int(valid.sum())
    print(f"  orbits where second prefix is meaningful (v1_odd >= 2^6): {n_valid:,} / {len(n_arr):,}")

    j1f = j1.astype(np.float64)[valid]
    j2f = j2.astype(np.float64)[valid]
    log_n_v = log_n[valid]
    sigma_v = sigma[valid]

    # Joint OLS sigma ~ log_n + j1 + j2
    beta, r2, se = joint_ols([log_n_v, j1f, j2f], sigma_v)
    print(f"\n  Joint OLS:  sigma = {beta[0]:.3f} + {beta[1]:.3f}*log(n) + {beta[2]:.3f}*j1 + {beta[3]:.3f}*j2")
    print(f"  SE: log(n)={se[1]:.3f}  j1={se[2]:.3f}  j2={se[3]:.3f}")
    print(f"  R^2 = {r2:.4f}")
    print(f"  Heuristic predicts both j1 and j2 coefficients should be ~+12.20")

    # Compare to single-level k=12
    j_table_12, _, _ = precompute_table(12)
    j12_v = j_table_12[n_arr[valid] % (1 << 12)].astype(np.float64)
    beta12, r2_12, se12 = joint_ols([log_n_v, j12_v], sigma_v)
    print(f"\n  COMPARISON: same-orbits single-level k=12:")
    print(f"    sigma = {beta12[0]:.3f} + {beta12[1]:.3f}*log(n) + {beta12[2]:.3f}*j12   (SE_j={se12[2]:.3f})")
    print(f"    R^2 = {r2_12:.4f}")

    # Try summed feature (j1 + j2) to see if it's the same as the level-12 prefix
    j_sum = j1f + j2f
    beta_s, r2_s, se_s = joint_ols([log_n_v, j_sum], sigma_v)
    print(f"\n  COMBINED: sigma ~ log(n) + (j1+j2) :")
    print(f"    sigma = {beta_s[0]:.3f} + {beta_s[1]:.3f}*log(n) + {beta_s[2]:.3f}*(j1+j2)")
    print(f"    R^2 = {r2_s:.4f}")

    # Try level-3 recursion (k1=k2=k3=6 -> "depth 3" decomposition)
    print(f"\n{'='*70}\nRecursive depth-3 decomposition (k1=k2=k3=6)\n{'='*70}")
    j1, v1, v1_odd, j2, valid12 = recursive_decompose(n_arr, 6, 6)
    j2_table, a2_table, c2_table = precompute_table(6)
    r2_v = v1_odd[valid12] % 64
    m2 = v1_odd[valid12] // 64
    a2 = a2_table[r2_v]
    c2 = c2_table[r2_v]
    v2 = a2 * m2 + c2
    v2_odd = v2.copy()
    em = (v2_odd > 0) & (v2_odd % 2 == 0)
    while em.any():
        v2_odd[em] //= 2
        em = (v2_odd > 0) & (v2_odd % 2 == 0)
    valid3 = v2_odd >= 64
    n_valid3 = int(valid3.sum())
    print(f"  orbits where third prefix is meaningful (v2_odd >= 2^6): {n_valid3:,}")

    if n_valid3 > 1000:
        j3_table, _, _ = precompute_table(6)
        j3 = j3_table[v2_odd[valid3] % 64].astype(np.float64)
        j1_3 = j1.astype(np.float64)[valid12][valid3]
        j2_3 = j2.astype(np.float64)[valid12][valid3]
        log_n_3 = log_n[valid12][valid3]
        sigma_3 = sigma[valid12][valid3]

        beta, r2, se = joint_ols([log_n_3, j1_3, j2_3, j3], sigma_3)
        print(f"\n  Joint OLS:  sigma = {beta[0]:.3f} + {beta[1]:.3f}*log(n) + "
              f"{beta[2]:.3f}*j1 + {beta[3]:.3f}*j2 + {beta[4]:.3f}*j3")
        print(f"  SE: log(n)={se[1]:.3f}  j1={se[2]:.3f}  j2={se[3]:.3f}  j3={se[4]:.3f}")
        print(f"  R^2 = {r2:.4f}")

        # Compare to single-level k=18
        j_table_18, _, _ = precompute_table(18)
        n_3 = n_arr[valid12][valid3]
        j18 = j_table_18[n_3 % (1 << 18)].astype(np.float64)
        beta18, r2_18, _ = joint_ols([log_n_3, j18], sigma_3)
        print(f"\n  COMPARISON: same-orbits single-level k=18:")
        print(f"    sigma = {beta18[0]:.3f} + {beta18[1]:.3f}*log(n) + {beta18[2]:.3f}*j18")
        print(f"    R^2 = {r2_18:.4f}")


if __name__ == "__main__":
    main()
