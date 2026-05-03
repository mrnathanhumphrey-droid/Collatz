"""
Experiment 19 — Bahadur-Rao prefactor correction test on q=7 X-binned data.

Cramer asymptotic:           log P(X) ~ const - theta*log(X)
Bahadur-Rao sharp asymptotic: log P(X) ~ const - theta*log(X) - 0.5*log(log(X))

For random walks with non-lattice step distributions, Bahadur-Rao gives the
1/sqrt(L) prefactor on top of the Cramer exponential rate. In our coordinate
L = log(X), this becomes a -0.5 log(log(X)) correction.

Fit both forms to the pooled q=7 X-binned data from experiment 18 and compare:
  - residual SSR
  - theta estimate
  - whether 5.5% deviation collapses

Usage:
    python 19_bahadur_rao.py --q 7 --N 1000000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import brentq, curve_fit

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix_q(r, q, a0, max_steps=400):
    a, c = a0, r
    odd_s = 0; even_s = 0
    while a % 2 == 0 and (odd_s + even_s) < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2; even_s += 1
        else:
            a *= q; c = q*c + 1; odd_s += 1
    return a, c, odd_s, even_s


def cramer_theta(q):
    def f(theta): return q ** (-theta) - (2 ** (1 - theta) - 1)
    return brentq(f, 0.001, 0.99)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=7)
    ap.add_argument("--N", type=int, default=1_000_000_000)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--n_bins", type=int, default=15)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    M = 1 << args.k
    K = M // 2

    a_final_per = np.zeros(K, dtype=np.int64)
    c_final_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r_class = 2 * kk + 1
        a, c, _, _ = deterministic_prefix_q(r_class, args.q, M)
        a_final_per[kk] = a
        c_final_per[kk] = c

    print(f"[load] q={args.q}, N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet")
    n = df["n"].to_numpy().astype(np.int64)
    converged = df["converged"].to_numpy().astype(bool)
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)
    m = (n - res) // M
    a_final_per_n = a_final_per[class_idx]
    c_final_per_n = c_final_per[class_idx]
    X = a_final_per_n.astype(np.float64) * m.astype(np.float64) + c_final_per_n.astype(np.float64)
    log_X = np.log(np.maximum(X, 1.0))

    bin_edges = np.linspace(log_X.min(), log_X.max(), args.n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_idx = np.clip(np.searchsorted(bin_edges[1:-1], log_X, side="right"), 0, args.n_bins - 1)
    total_per_bin = np.bincount(bin_idx, minlength=args.n_bins)
    conv_per_bin = np.bincount(bin_idx[converged], minlength=args.n_bins)
    rate_per_bin = conv_per_bin / np.maximum(total_per_bin, 1)
    keep = (rate_per_bin > 0) & (total_per_bin >= 100) & (bin_centers > 1.0)

    print()
    print(f"[bins] usable bins: {keep.sum()}")

    log_X_v = bin_centers[keep]
    log_rate_v = np.log(rate_per_bin[keep])
    weights = conv_per_bin[keep].astype(float)

    theta_cramer = cramer_theta(args.q)
    print(f"[theta] Cramer prediction theta(q={args.q}) = {theta_cramer:.5f}")

    # Model 1: Cramer  log_rate = a + b * log(X)
    # Weighted OLS with sigma_i ~ 1/sqrt(weights_i) (Poisson-like errors)
    A1 = np.column_stack([np.ones(len(log_X_v)), log_X_v])
    W = np.diag(weights)
    AtWA = A1.T @ W @ A1
    AtWy = A1.T @ W @ log_rate_v
    coef1 = np.linalg.solve(AtWA, AtWy)
    pred1 = A1 @ coef1
    resid1 = log_rate_v - pred1
    ssr1 = np.sum(weights * resid1 ** 2)
    sst = np.sum(weights * (log_rate_v - np.average(log_rate_v, weights=weights)) ** 2)
    r2_1 = 1.0 - ssr1 / sst

    # Model 2: Bahadur-Rao  log_rate = a + b * log(X) - 0.5 * log(log(X))
    # First treat the -0.5 log(log(X)) as a known term: subtract it from log_rate, fit residual
    log_log_X = np.log(np.maximum(log_X_v, 1.0))
    log_rate_corrected = log_rate_v + 0.5 * log_log_X
    AtWy2 = A1.T @ W @ log_rate_corrected
    coef2 = np.linalg.solve(AtWA, AtWy2)
    pred2 = A1 @ coef2 - 0.5 * log_log_X  # back to log_rate prediction
    resid2 = log_rate_v - pred2
    ssr2 = np.sum(weights * resid2 ** 2)
    r2_2 = 1.0 - ssr2 / sst

    # Model 3: free prefactor exponent  log_rate = a + b * log(X) + c * log(log(X))
    A3 = np.column_stack([np.ones(len(log_X_v)), log_X_v, log_log_X])
    AtWA3 = A3.T @ W @ A3
    AtWy3 = A3.T @ W @ log_rate_v
    coef3 = np.linalg.solve(AtWA3, AtWy3)
    pred3 = A3 @ coef3
    resid3 = log_rate_v - pred3
    ssr3 = np.sum(weights * resid3 ** 2)
    r2_3 = 1.0 - ssr3 / sst

    print()
    print(f"=== Model comparison ===")
    print()
    print(f"Model 1 (Cramer pure): log_rate = a + b * log(X)")
    print(f"  fitted: a = {coef1[0]:.4f},  b = {coef1[1]:.4f}")
    print(f"  Cramer prediction:        b = {-theta_cramer:.4f}")
    print(f"  ratio b_emp/b_pred = {coef1[1]/-theta_cramer:.4f}")
    print(f"  weighted SSR = {ssr1:.4f}    R^2 = {r2_1:.4f}")

    print()
    print(f"Model 2 (Bahadur-Rao, -0.5 prefactor fixed): log_rate = a + b * log(X) - 0.5 * log(log(X))")
    print(f"  fitted: a = {coef2[0]:.4f},  b = {coef2[1]:.4f}")
    print(f"  Cramer prediction:        b = {-theta_cramer:.4f}")
    print(f"  ratio b_emp/b_pred = {coef2[1]/-theta_cramer:.4f}")
    print(f"  weighted SSR = {ssr2:.4f}    R^2 = {r2_2:.4f}")

    print()
    print(f"Model 3 (free log(log) coefficient): log_rate = a + b * log(X) + c * log(log(X))")
    print(f"  fitted: a = {coef3[0]:.4f},  b = {coef3[1]:.4f},  c = {coef3[2]:.4f}")
    print(f"  Cramer prediction:        b = {-theta_cramer:.4f}")
    print(f"  Bahadur-Rao prediction:   c = -0.5")
    print(f"  ratio b_emp/b_pred = {coef3[1]/-theta_cramer:.4f}")
    print(f"  weighted SSR = {ssr3:.4f}    R^2 = {r2_3:.4f}")

    print()
    print(f"=== Verdict ===")
    print(f"  Cramer-only b = {coef1[1]:>+.4f}    (theta_emp = {-coef1[1]:.4f}; deviation {100*(coef1[1]/-theta_cramer - 1):+.1f}%)")
    print(f"  BR-fixed   b = {coef2[1]:>+.4f}    (theta_emp = {-coef2[1]:.4f}; deviation {100*(coef2[1]/-theta_cramer - 1):+.1f}%)")
    print(f"  BR-free    b = {coef3[1]:>+.4f}    c = {coef3[2]:>+.4f}    (deviation {100*(coef3[1]/-theta_cramer - 1):+.1f}%)")
    print()
    if abs(coef2[1]/-theta_cramer - 1) < abs(coef1[1]/-theta_cramer - 1) / 2:
        print(f"  *** Bahadur-Rao prefactor cuts the deviation in half or more.")
        print(f"      The 5.5% gap is largely explained by the 1/sqrt(L) correction.")
    elif abs(coef3[2] - (-0.5)) < 0.15:
        print(f"  *** Free fit recovers c ~ -0.5, consistent with Bahadur-Rao theoretical prefactor.")
    else:
        print(f"  *** Bahadur-Rao prefactor does NOT cleanly explain the deviation.")
        print(f"      Free coefficient c = {coef3[2]:.3f} vs theoretical -0.5; mismatch.")


if __name__ == "__main__":
    main()
