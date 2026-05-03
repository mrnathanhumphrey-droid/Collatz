"""
Experiment 20 — m-selection confound test for q=7 conv_rate signature.

Pooled X-binned fit at q=7 N=10^9 gave a free-coefficient log(log(X)) term
with c = +0.15 (opposite sign from Bahadur-Rao). Three candidate explanations:
  (1) Class-specific c_final corrections — ruled out by tight per-class slopes
  (2) Lattice random-walk discreteness corrections to Cramer
  (3) m-distribution selection: at fixed X across classes, different m values
      are sampled, and if P(reach 1) depends on m beyond just X, the pooled
      X-fit picks up a confound

This script tests #3:
  Bin orbits by log(X) AND by log(m) jointly. If conv_rate factorizes — i.e.,
  P(converge | X, m) = f(X) — then m-selection is not the issue. If conv_rate
  varies systematically with m at fixed X, m-selection is the confound and
  the +log(log) signature is methodological, not a lattice phenomenon.

The cleanest single number: fit log(rate) ~ a + b*log(X) + d*log(m) and report
the partial coefficient d with its standard error. Under the null (X is the
only relevant variable), d ~ 0. Under m-selection, d != 0.

Usage:
    python 20_m_selection_test.py --q 7 --N 1000000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import t as t_dist

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
    ap.add_argument("--n_logx_bins", type=int, default=10)
    ap.add_argument("--n_logm_bins", type=int, default=10)
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
    log_m = np.log(np.maximum(m.astype(np.float64), 1.0))
    print(f"        rows={len(n):,}, conv={converged.sum():,}", flush=True)

    # 2D bin by (log_X, log_m). Equal-width in each dimension.
    lx_min, lx_max = log_X.min(), log_X.max()
    lm_min, lm_max = log_m.min(), log_m.max()
    lx_edges = np.linspace(lx_min, lx_max, args.n_logx_bins + 1)
    lm_edges = np.linspace(lm_min, lm_max, args.n_logm_bins + 1)

    lx_idx = np.clip(np.searchsorted(lx_edges[1:-1], log_X, side="right"), 0, args.n_logx_bins - 1)
    lm_idx = np.clip(np.searchsorted(lm_edges[1:-1], log_m, side="right"), 0, args.n_logm_bins - 1)

    # 2D count and conv_count
    flat_idx = lx_idx * args.n_logm_bins + lm_idx
    total_2d = np.bincount(flat_idx, minlength=args.n_logx_bins * args.n_logm_bins)
    conv_2d = np.bincount(flat_idx[converged], minlength=args.n_logx_bins * args.n_logm_bins)
    rate_2d = (conv_2d / np.maximum(total_2d, 1)).reshape(args.n_logx_bins, args.n_logm_bins)
    total_2d = total_2d.reshape(args.n_logx_bins, args.n_logm_bins)
    conv_2d = conv_2d.reshape(args.n_logx_bins, args.n_logm_bins)

    # Build regression dataset: keep cells with rate > 0 and total >= some minimum
    cells = []
    for ix in range(args.n_logx_bins):
        for im in range(args.n_logm_bins):
            r = rate_2d[ix, im]
            if r <= 0 or total_2d[ix, im] < 100:
                continue
            lx_ctr = 0.5 * (lx_edges[ix] + lx_edges[ix + 1])
            lm_ctr = 0.5 * (lm_edges[im] + lm_edges[im + 1])
            cells.append((lx_ctr, lm_ctr, np.log(r), int(conv_2d[ix, im]), int(total_2d[ix, im])))

    if len(cells) < 4:
        print("[warn] not enough usable 2D cells for regression")
        return

    arr = np.array(cells, dtype=np.float64)
    log_X_v = arr[:, 0]
    log_m_v = arr[:, 1]
    log_rate_v = arr[:, 2]
    n_conv_w = arr[:, 3]

    # Weighted regression: log(rate) ~ a + b*log(X) + d*log(m)
    A = np.column_stack([np.ones(len(log_X_v)), log_X_v, log_m_v])
    W = np.diag(n_conv_w)
    AtWA = A.T @ W @ A
    AtWy = A.T @ W @ log_rate_v
    coef = np.linalg.solve(AtWA, AtWy)
    pred = A @ coef
    resid = log_rate_v - pred
    sw = n_conv_w.sum()
    s2 = (n_conv_w * resid ** 2).sum() / (sw - len(coef))
    cov = s2 * np.linalg.inv(AtWA) * sw / (n_conv_w.sum())  # rough
    se = np.sqrt(np.diag(cov))
    t_stats = coef / se
    df_resid = len(log_X_v) - len(coef)
    p_vals = 2 * (1 - t_dist.cdf(np.abs(t_stats), df=df_resid))

    # Also fit X-only model as baseline
    A1 = np.column_stack([np.ones(len(log_X_v)), log_X_v])
    AtWA1 = A1.T @ W @ A1
    AtWy1 = A1.T @ W @ log_rate_v
    coef1 = np.linalg.solve(AtWA1, AtWy1)
    pred1 = A1 @ coef1
    resid1 = log_rate_v - pred1
    ssr1 = (n_conv_w * resid1 ** 2).sum()
    ssr2 = (n_conv_w * resid ** 2).sum()
    sst = (n_conv_w * (log_rate_v - np.average(log_rate_v, weights=n_conv_w)) ** 2).sum()
    r2_1 = 1 - ssr1 / sst
    r2_2 = 1 - ssr2 / sst

    print()
    print(f"=== Joint (log(X), log(m)) regression at q={args.q} ===")
    print(f"  Usable 2D cells: {len(log_X_v)}  (out of {args.n_logx_bins * args.n_logm_bins})")
    print()
    print(f"Model A: log(rate) = a + b*log(X)")
    print(f"  intercept a = {coef1[0]:.4f}")
    print(f"  slope b     = {coef1[1]:.4f}")
    print(f"  R^2 = {r2_1:.4f}    SSR = {ssr1:.4f}")
    print()
    print(f"Model B: log(rate) = a + b*log(X) + d*log(m)")
    print(f"  intercept a = {coef[0]:.4f}    SE = {se[0]:.4f}    t = {t_stats[0]:>7.3f}    p = {p_vals[0]:.3e}")
    print(f"  slope b     = {coef[1]:.4f}    SE = {se[1]:.4f}    t = {t_stats[1]:>7.3f}    p = {p_vals[1]:.3e}")
    print(f"  partial d   = {coef[2]:.4f}    SE = {se[2]:.4f}    t = {t_stats[2]:>7.3f}    p = {p_vals[2]:.3e}")
    print(f"  R^2 = {r2_2:.4f}    SSR = {ssr2:.4f}")

    # Cramer prediction
    from scipy.optimize import brentq
    def f(theta): return args.q ** (-theta) - (2 ** (1 - theta) - 1)
    theta_cramer = brentq(f, 0.001, 0.99)
    print()
    print(f"Cramer prediction: b = {-theta_cramer:.4f}")
    print(f"  Model A b deviation: {100*(coef1[1]/-theta_cramer - 1):>+.1f}%")
    print(f"  Model B b deviation: {100*(coef[1]/-theta_cramer - 1):>+.1f}%")

    print()
    print(f"=== Verdict on m-selection confound ===")
    if p_vals[2] < 0.05 and abs(coef[2]) > 0.05:
        print(f"  partial d = {coef[2]:+.4f} (p = {p_vals[2]:.3e}) is significant.")
        print(f"  m carries information BEYOND X. Pooled X-fit is confounded by m-selection.")
        if abs(coef[1] / -theta_cramer - 1) < 0.02:
            print(f"  After regressing out m, b matches Cramer to {100*abs(coef[1] / -theta_cramer - 1):.1f}%.")
            print(f"  --> Candidate #3 (m-selection) is the explanation. The 5.5% deviation")
            print(f"      is methodological, not a lattice phenomenon.")
        else:
            print(f"  But Model B b still deviates from Cramer by {100*(coef[1]/-theta_cramer - 1):+.1f}%.")
            print(f"  --> m-selection is part of the story but not the whole story.")
    else:
        print(f"  partial d = {coef[2]:+.4f} (p = {p_vals[2]:.3e}) is not significant.")
        print(f"  At fixed X, m carries no information about conv_rate.")
        print(f"  --> Candidate #3 (m-selection) ruled out. Pursue lattice-correction story (#2).")


if __name__ == "__main__":
    main()
