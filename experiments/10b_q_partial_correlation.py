"""
Experiment 10b — Partial correlation analysis for qx+1 convergence prediction.

The Pearson r=-0.92 for r(log(a_final), conv_rate) at q=5 could be driven by
genuine a_final dependence OR by collinearity with prefix length. Both are
deterministically related: log_q(a_final) = odd_steps_in_prefix.

This script regresses per-class convergence rate on three candidate predictors
simultaneously, reports partial regression coefficients and partial correlations:

  predictors:
    1. odd_steps_in_prefix  (== log_q(a_final))
    2. even_steps_in_prefix (= prefix_steps_total - odd_steps_in_prefix)
    3. log(c_final + 1)     (additive offset magnitude, secondary modulation)

If odd_steps survives at high partial r and even_steps drops out, the finding
is "convergence depends on the q-growth count in the prefix, not its length."
If even_steps also predicts independently, it's mixed. If both drop in partial
regression, the finding is weaker than the marginal r=-0.92 suggests.

Usage:
    python 10b_q_partial_correlation.py --q 5 --N 1000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix_q(r, q, a0, max_steps=400):
    """Track odd_steps + even_steps separately."""
    a, c = a0, r
    odd_s = 0; even_s = 0
    while a % 2 == 0 and (odd_s + even_s) < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
            even_s += 1
        else:
            a *= q
            c = q * c + 1
            odd_s += 1
    return a, c, odd_s, even_s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--N", type=int, default=1_000_000)
    ap.add_argument("--k", type=int, default=6)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    M = 1 << args.k
    K = M // 2

    print(f"[load] q={args.q}, N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet")
    n = df["n"].to_numpy().astype(np.int64)
    converged = df["converged"].to_numpy().astype(bool)
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    # Per-class deterministic prefix
    a_final_per = np.zeros(K, dtype=np.int64)
    c_final_per = np.zeros(K, dtype=np.int64)
    odd_s_per = np.zeros(K, dtype=np.int64)
    even_s_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r_class = 2 * kk + 1
        a, c, os_, es_ = deterministic_prefix_q(r_class, args.q, M)
        a_final_per[kk] = a
        c_final_per[kk] = c
        odd_s_per[kk] = os_
        even_s_per[kk] = es_

    # Per-class convergence rate
    n_total_per = np.zeros(K, dtype=np.int64)
    n_conv_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        m = class_idx == kk
        n_total_per[kk] = m.sum()
        n_conv_per[kk] = (m & converged).sum()
    conv_rate = n_conv_per / np.maximum(n_total_per, 1)

    # Build predictor matrix and response
    log_a_final = np.log(a_final_per.astype(np.float64))
    prefix_total = odd_s_per + even_s_per
    log_c_final = np.log(c_final_per.astype(np.float64) + 1.0)
    y = conv_rate.astype(np.float64)

    # Marginal Pearson r
    print()
    print(f"=== Marginal correlations (q={args.q}, k={args.k}, N={args.N:,}) ===")
    for name, x in [("odd_steps_in_prefix", odd_s_per.astype(float)),
                    ("even_steps_in_prefix", even_s_per.astype(float)),
                    ("prefix_total_steps", prefix_total.astype(float)),
                    ("log(a_final)", log_a_final),
                    ("log(c_final+1)", log_c_final)]:
        if np.std(x) < 1e-9:
            print(f"  {name:>22}: var=0, skip")
            continue
        r_, p_ = pearsonr(x, y)
        rs_, ps_ = spearmanr(x, y)
        print(f"  {name:>22}: Pearson r={r_:>7.4f}  p={p_:.3e}    Spearman rs={rs_:>7.4f}  p={ps_:.3e}")

    # Print collinearity matrix among predictors
    preds = {
        "odd_steps": odd_s_per.astype(float),
        "even_steps": even_s_per.astype(float),
        "log_a_final": log_a_final,
        "log_c_final": log_c_final,
    }
    print()
    print(f"=== Predictor collinearity (Pearson r) ===")
    names = list(preds.keys())
    print(f"  {'':>14}", " ".join(f"{n:>12}" for n in names))
    for n1 in names:
        row = [f"{n1:>14}"]
        for n2 in names:
            if np.std(preds[n1]) < 1e-9 or np.std(preds[n2]) < 1e-9:
                row.append(f"{'--':>12}")
            else:
                r_, _ = pearsonr(preds[n1], preds[n2])
                row.append(f"{r_:>12.4f}")
        print(" ".join(row))

    # Multivariate regression with predictors that have variance:
    #   y ~ b0 + b1*odd_steps + b2*even_steps + b3*log_c_final
    # Standardize X for interpretable coefficients.
    keep_preds = []
    keep_names = []
    for name in ["odd_steps", "even_steps", "log_c_final"]:
        x = preds[name]
        if np.std(x) < 1e-9:
            print(f"\n  [skip] {name} has no variance, dropping from regression")
            continue
        keep_preds.append((x - x.mean()) / x.std())
        keep_names.append(name)
    if not keep_preds:
        print("  no usable predictors")
        return

    X = np.column_stack(keep_preds)
    y_std = (y - y.mean()) / y.std() if y.std() > 0 else y - y.mean()

    # OLS via normal equations
    XtX = X.T @ X
    Xty = X.T @ y_std
    try:
        coefs = np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError:
        coefs, *_ = np.linalg.lstsq(X, y_std, rcond=None)

    # Compute residuals + standard errors
    y_pred = X @ coefs
    resid = y_std - y_pred
    n_obs = len(y_std)
    n_par = X.shape[1]
    if n_obs > n_par + 1:
        s2 = (resid ** 2).sum() / (n_obs - n_par - 1)
        try:
            cov = s2 * np.linalg.inv(XtX)
            se = np.sqrt(np.diag(cov))
        except np.linalg.LinAlgError:
            se = np.full(n_par, np.nan)
        t_stats = coefs / se
        from scipy.stats import t as t_dist
        p_vals = 2 * (1 - t_dist.cdf(np.abs(t_stats), df=n_obs - n_par - 1))
    else:
        se = np.full(n_par, np.nan)
        t_stats = np.full(n_par, np.nan)
        p_vals = np.full(n_par, np.nan)

    # R^2
    ss_total = ((y_std - y_std.mean()) ** 2).sum()
    ss_resid = (resid ** 2).sum()
    r_sq = 1.0 - ss_resid / ss_total if ss_total > 0 else float("nan")

    # Partial correlations: r_partial(xi, y | other x's) = sign(coef) * sqrt(t^2 / (t^2 + df))
    df_resid = n_obs - n_par - 1
    partial_r = []
    for ci, ti in zip(coefs, t_stats):
        if not np.isnan(ti):
            sign = 1.0 if ci > 0 else -1.0
            pr = sign * np.sqrt(ti ** 2 / (ti ** 2 + df_resid))
        else:
            pr = np.nan
        partial_r.append(pr)

    print()
    print(f"=== Multivariate regression: conv_rate (std) ~ standardized predictors ===")
    print(f"  R^2 (whole model)        = {r_sq:.4f}")
    print(f"  N classes                = {n_obs}")
    print(f"  df_residual              = {df_resid}")
    print()
    print(f"  {'predictor':>15} {'beta':>9} {'SE':>8} {'t':>8} {'p':>10} {'partial_r':>11}")
    for name, b, s, t, p, pr in zip(keep_names, coefs, se, t_stats, p_vals, partial_r):
        print(f"  {name:>15} {b:>9.4f} {s:>8.4f} {t:>8.3f} {p:>10.3e} {pr:>11.4f}")

    print()
    print("=== Verdict ===")
    for name, ti, pi, pr in zip(keep_names, t_stats, p_vals, partial_r):
        sig = "**INDEPENDENT**" if (pi < 0.01) else ("marginal" if (pi < 0.05) else "not independent")
        print(f"  {name:>15}: partial r = {pr:>+.4f}, p = {pi:.3e}  --> {sig}")


if __name__ == "__main__":
    main()
