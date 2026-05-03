"""
Experiment 18 — Three diagnostics on the q=7 j-slope deviation.

Cramer asymptotic theorem predicts P(reach 1 from log-value V) = exp(-theta * V).
We see:
  q=5: matches to 0.01% on slope-vs-j
  q=7: 14% deviation, persists at N=10^9

Three candidate explanations:
  (1) Finite-j corrections: at same X = q^j*m + c_final, classes with different
      (j, c_final) decompositions give different P(reach 1)
  (2) Boundary effect at integer state 1: P(X) deviates from X^(-theta) at small X
  (3) Step correlation: i.i.d. assumption fails

To discriminate, we bin all q=7 convergent + non-convergent orbits by X (not j),
fit log(conv_rate) vs log(X), and look at:
  A. Pooled slope (all 32 classes pooled): should be -theta(q) = -0.6265 if X is
     the only relevant variable
  B. Per-class slopes: if they vary at fixed X, finite-j corrections (cand #1)
  C. Residuals from pooled fit, plotted against j: if correlated with j,
     class-specific structure (also cand #1). If curvature in X-fit, cand #2.

Usage:
    python 18_q7_x_binning_diagnostic.py --q 7 --N 1000000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import linregress, pearsonr

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=7)
    ap.add_argument("--N", type=int, default=1_000_000_000)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--n_logx_bins", type=int, default=15)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    M = 1 << args.k
    K = M // 2

    # Per-class deterministic prefix lookup
    a_final_per = np.zeros(K, dtype=np.int64)
    c_final_per = np.zeros(K, dtype=np.int64)
    odd_s_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r_class = 2 * kk + 1
        a, c, os_, es_ = deterministic_prefix_q(r_class, args.q, M)
        a_final_per[kk] = a
        c_final_per[kk] = c
        odd_s_per[kk] = os_

    print(f"[load] q={args.q}, N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet")
    n = df["n"].to_numpy().astype(np.int64)
    converged = df["converged"].to_numpy().astype(bool)
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)
    print(f"        rows={len(n):,}, converged={converged.sum():,}", flush=True)

    # Compute m and X = a_final*m + c_final per orbit
    m = (n - res) // M  # integer division (n = res + m*M)
    a_final_per_n = a_final_per[class_idx]
    c_final_per_n = c_final_per[class_idx]
    j_per_n = odd_s_per[class_idx]
    # Use float64 for X since q^j * m can be huge
    X = a_final_per_n.astype(np.float64) * m.astype(np.float64) + c_final_per_n.astype(np.float64)
    log_X = np.log(np.maximum(X, 1.0))

    # log(X) range across data
    lx_min = log_X.min()
    lx_max = log_X.max()
    print(f"[range] log(X) in [{lx_min:.2f}, {lx_max:.2f}]", flush=True)

    # Pooled bin: divide log(X) range into N_bins of equal width, count total + conv per bin
    bin_edges = np.linspace(lx_min, lx_max, args.n_logx_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    bin_idx = np.clip(np.searchsorted(bin_edges[1:-1], log_X, side="right"), 0, args.n_logx_bins - 1)
    total_per_bin = np.bincount(bin_idx, minlength=args.n_logx_bins)
    conv_per_bin = np.bincount(bin_idx[converged], minlength=args.n_logx_bins)
    rate_per_bin = conv_per_bin / np.maximum(total_per_bin, 1)

    print()
    print(f"=== A. Pooled fit (all 32 classes binned by log(X)) ===")
    print(f"{'bin':>4} {'log(X) ctr':>10} {'n_total':>14} {'n_conv':>9} {'conv_rate':>14} {'log(rate)':>11}")
    keep_mask = (rate_per_bin > 0) & (total_per_bin >= 100)
    for i in range(args.n_logx_bins):
        rate = rate_per_bin[i]
        log_rate = np.log(rate) if rate > 0 else float("-inf")
        kept = "*" if keep_mask[i] else " "
        log_rate_str = f"{log_rate:.4f}" if rate > 0 else "      -inf"
        print(f"{i:>4} {bin_centers[i]:>10.3f} {total_per_bin[i]:>14,} {conv_per_bin[i]:>9,} "
              f"{rate:>14.4e} {log_rate_str:>11}{kept}")

    if keep_mask.sum() < 3:
        print("\n[warn] fewer than 3 usable bins; cannot fit")
        return
    log_rate_vals = np.log(rate_per_bin[keep_mask])
    log_X_vals = bin_centers[keep_mask]
    weights = conv_per_bin[keep_mask].astype(float)

    res_pool_u = linregress(log_X_vals, log_rate_vals)
    sum_w = weights.sum()
    x_w = (weights * log_X_vals).sum() / sum_w
    y_w = (weights * log_rate_vals).sum() / sum_w
    sxx = (weights * (log_X_vals - x_w) ** 2).sum()
    sxy = (weights * (log_X_vals - x_w) * (log_rate_vals - y_w)).sum()
    b_w = sxy / sxx
    a_w = y_w - b_w * x_w

    # Cramer prediction for q
    from scipy.optimize import brentq
    def f(theta): return args.q ** (-theta) - (2 ** (1 - theta) - 1)
    theta_cramer = brentq(f, 0.001, 0.99)

    print()
    print(f"  Pooled OLS (unweighted): slope = {res_pool_u.slope:>9.4f}  intercept = {res_pool_u.intercept:>9.4f}  R^2 = {res_pool_u.rvalue**2:.4f}")
    print(f"  Pooled OLS (weighted):   slope = {b_w:>9.4f}  intercept = {a_w:>9.4f}")
    print(f"  Cramer prediction:       slope = {-theta_cramer:>9.4f}  (theta = {theta_cramer:.4f})")
    print(f"  Ratio empirical (weighted) / cramer: {b_w / -theta_cramer:.4f}")

    print()
    print(f"=== B. Per-class fit (each of 32 odd classes separately) ===")
    print(f"  binning each class's orbits by log(X) and fitting log(rate) vs log(X)")
    print()
    print(f"{'class':>5} {'j':>3} {'a_final':>8} {'c_final':>8} {'n_conv':>7} {'slope':>9} {'R^2':>7} {'n_bins':>7}")
    per_class_slopes = []
    per_class_jvals = []
    for kk in range(K):
        m_c = class_idx == kk
        n_c = n[m_c]
        log_X_c = log_X[m_c]
        conv_c = converged[m_c]
        n_conv_c = conv_c.sum()
        if n_conv_c < 10:
            continue
        # Bin within class
        n_bins_c = max(3, min(8, n_conv_c // 5))
        edges_c = np.linspace(log_X_c.min(), log_X_c.max(), n_bins_c + 1)
        bin_idx_c = np.clip(np.searchsorted(edges_c[1:-1], log_X_c, side="right"), 0, n_bins_c - 1)
        tot_c = np.bincount(bin_idx_c, minlength=n_bins_c)
        cnv_c = np.bincount(bin_idx_c[conv_c], minlength=n_bins_c)
        rate_c = cnv_c / np.maximum(tot_c, 1)
        keep_c = (rate_c > 0) & (tot_c >= 100)
        if keep_c.sum() < 3:
            continue
        ctrs_c = 0.5 * (edges_c[:-1] + edges_c[1:])
        try:
            r_lr = linregress(ctrs_c[keep_c], np.log(rate_c[keep_c]))
            per_class_slopes.append(r_lr.slope)
            per_class_jvals.append(int(odd_s_per[kk]))
            print(f"{kk:>5} {int(odd_s_per[kk]):>3} {int(a_final_per[kk]):>8} "
                  f"{int(c_final_per[kk]):>8} {int(n_conv_c):>7} "
                  f"{r_lr.slope:>9.4f} {r_lr.rvalue**2:>7.3f} {keep_c.sum():>7}")
        except Exception as e:
            continue

    if len(per_class_slopes) >= 3:
        slopes = np.array(per_class_slopes)
        print()
        print(f"  Per-class slope statistics ({len(slopes)} classes):")
        print(f"    mean   = {slopes.mean():>9.4f}")
        print(f"    median = {np.median(slopes):>9.4f}")
        print(f"    std    = {slopes.std(ddof=1):>9.4f}")
        print(f"    range  = [{slopes.min():.4f}, {slopes.max():.4f}]")
        print(f"    Cramer prediction: {-theta_cramer:.4f}")

    print()
    print(f"=== C. Residual structure: residuals from pooled fit vs j ===")
    # For each orbit (or each class), compute residual from pooled fit
    # Use per-class mean residual for tractability
    print(f"  Per-class mean residual after subtracting pooled fit:")
    print(f"{'class':>5} {'j':>3} {'log(X)_med':>11} {'mean_resid':>11}")
    resid_per_class = []
    j_per_class = []
    for kk in range(K):
        m_c = class_idx == kk
        log_X_c = log_X[m_c]
        conv_c = converged[m_c]
        n_conv_c = conv_c.sum()
        if n_conv_c < 5:
            continue
        # Take median log(X) of convergent orbits in this class
        log_X_conv_c = log_X_c[conv_c]
        if len(log_X_conv_c) == 0:
            continue
        log_X_med = float(np.median(log_X_conv_c))
        # Predicted log(rate) at this log(X)
        predicted_log_rate = a_w + b_w * log_X_med
        # Empirical log(rate) for this class as a whole
        rate_c = n_conv_c / m_c.sum()
        if rate_c <= 0:
            continue
        emp_log_rate = np.log(rate_c)
        resid = emp_log_rate - predicted_log_rate
        resid_per_class.append(resid)
        j_per_class.append(int(odd_s_per[kk]))
        print(f"{kk:>5} {int(odd_s_per[kk]):>3} {log_X_med:>11.3f} {resid:>11.4f}")

    if len(resid_per_class) >= 3:
        resid = np.array(resid_per_class)
        j_arr = np.array(j_per_class)
        if np.std(j_arr) > 0:
            r_pj, p_pj = pearsonr(j_arr, resid)
            print(f"\n  Pearson r(j, residual) = {r_pj:>+.4f}  p = {p_pj:.3e}")
            if abs(r_pj) < 0.2:
                print(f"  --> j-residuals are flat. X is the dominant variable.")
                print(f"      The j-slope deviation comes from class-X-range correlation, NOT")
                print(f"      from class-specific dynamics at fixed X.")
            elif abs(r_pj) > 0.5:
                print(f"  --> j-residuals correlate strongly with j. Candidate #1 (finite-j /")
                print(f"      class-specific c_final corrections) is the explanation.")
            else:
                print(f"  --> moderate j-correlation ({r_pj:.2f}); some class-specific structure beyond X.")


if __name__ == "__main__":
    main()
