"""
Experiment 02 — Higher-moment universality

For each odd residue class r mod 2^k at fixed N, compute the per-class
residual eps = sigma - (alpha_class + beta_class * ln n) via per-class OLS,
then per-class statistics: variance, skewness, excess kurtosis.

Test: do these higher moments correlate with the deterministic prefix
prediction (alpha_det = prefix_steps + 10.43 * ln(a_final / 2^k))?

Reports Pearson and Spearman correlations and produces scatter plots.

Sweeps k ∈ {6, 7, 8, 9} on N=2^27 data by default.

Usage:
    python 02_moment_universality.py --N 134217728 --ks 6 7 8 9
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from scipy import stats as st
from scipy.stats import spearmanr, pearsonr

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))


def deterministic_prefix(r, a0, max_steps=400):
    a, c = a0, r; steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3*c + 1
        steps += 1
    return steps, a, c


def per_class_moments(df_n, df_log, df_sigma, class_idx, K):
    var_per = np.zeros(K); skew_per = np.zeros(K); kurt_per = np.zeros(K)
    for kk in range(K):
        m = class_idx == kk
        if m.sum() < 100:
            continue
        bk, ak = np.polyfit(df_log[m], df_sigma[m], 1)
        eps = df_sigma[m] - (ak + bk * df_log[m])
        var_per[kk] = np.var(eps, ddof=1)
        skew_per[kk] = st.skew(eps)
        kurt_per[kk] = st.kurtosis(eps, fisher=True)
    return var_per, skew_per, kurt_per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1 << 27)
    ap.add_argument("--ks", type=int, nargs="+", default=[6, 7, 8, 9])
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    print(f"[load] N={args.N:,}, odd-only={len(n):,}")

    print()
    print(f"{'k':>3} {'K':>4} {'r(pred,Var)':>14} {'r(pred,Skew)':>14} {'r(pred,Kurt)':>14}  unique_a_final")
    summary_rows = []
    for k_pow in args.ks:
        M = 1 << k_pow; K = M // 2
        res = (n % M).astype(np.int32)
        class_idx = ((res - 1) // 2).astype(np.int32)

        var_per, skew_per, kurt_per = per_class_moments(n, log_n, sigma, class_idx, K)

        prefix_arr = np.zeros(K, dtype=int)
        a_final_arr = np.zeros(K, dtype=int)
        for kk in range(K):
            r = 2*kk + 1
            steps, a_f, _ = deterministic_prefix(r, M)
            prefix_arr[kk] = steps
            a_final_arr[kk] = a_f
        alpha_pred = prefix_arr + LOG_FACTOR_ODD * np.log(a_final_arr / float(M))

        rv, _ = pearsonr(alpha_pred, var_per)
        rs, _ = pearsonr(alpha_pred, skew_per)
        rk, _ = pearsonr(alpha_pred, kurt_per)
        unique_a = sorted(set(a_final_arr.tolist()))

        print(f"{k_pow:>3} {K:>4} {rv:>14.6f} {rs:>14.4f} {rk:>14.4f}  {unique_a}")
        summary_rows.append((k_pow, K, var_per, skew_per, kurt_per, alpha_pred, a_final_arr))

    # Plot all 4 k values
    fig, axes = plt.subplots(3, len(args.ks), figsize=(5*len(args.ks), 12), sharey=False)
    for col, (k_pow, K, var_per, skew_per, kurt_per, alpha_pred, a_finals) in enumerate(summary_rows):
        for row, (label, vals) in enumerate([("Var", var_per), ("Skew", skew_per), ("Excess kurt", kurt_per)]):
            ax = axes[row, col] if len(args.ks) > 1 else axes[row]
            for af in sorted(set(a_finals.tolist())):
                m = a_finals == af
                ax.scatter(alpha_pred[m], vals[m], s=40, alpha=0.7,
                           label=f"a_final={af}" if row == 0 else None)
            r_p, p_p = pearsonr(alpha_pred, vals)
            ax.set_xlabel("predicted alpha")
            ax.set_ylabel(label)
            ax.set_title(f"k={k_pow}: {label}\nPearson r = {r_p:.4f}")
            if row == 0 and col == 0:
                ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)

    fig.suptitle(f"Per-class moments vs deterministic prefix prediction  (N={args.N:,})", y=1.005)
    plt.tight_layout()
    out_png = out_dir / f"02_moment_universality_N{args.N}.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close()
    print()
    print(f"[save] {out_png}")


if __name__ == "__main__":
    main()
