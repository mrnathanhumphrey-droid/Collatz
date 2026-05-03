"""
Experiment 21 — Two-term fit testing finite-path-floor mechanism for q=7
small-X deviation from pure Cramer.

Hypothesis: P(reach 1 from X) ~ A * X^(-theta) + B * 1[X <= X_*]

For small X, orbits reach 1 via short deterministic paths (B-term floor);
for large X, random walk descent dominates (A * X^(-theta)). The +log(log)
signature in the pure linear log-log fit is then explained as curvature
from the additive small-X floor, not a lattice/Markov correction to Cramer.

We fix theta = theta_Cramer(q) and fit (A, B, X_*) on the pooled X-binned data.

Usage:
    python 21_two_term_fit.py --q 7 --N 1000000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import brentq, minimize

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

    # Build X-binned data
    bin_edges = np.linspace(log_X.min(), log_X.max(), args.n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_idx = np.clip(np.searchsorted(bin_edges[1:-1], log_X, side="right"), 0, args.n_bins - 1)
    total_per_bin = np.bincount(bin_idx, minlength=args.n_bins)
    conv_per_bin = np.bincount(bin_idx[converged], minlength=args.n_bins)
    keep = (total_per_bin >= 100) & (bin_centers > 0.5)

    log_X_v = bin_centers[keep]
    X_v = np.exp(log_X_v)
    n_total_v = total_per_bin[keep]
    n_conv_v = conv_per_bin[keep]
    rate_v = n_conv_v / np.maximum(n_total_v, 1)

    theta = cramer_theta(args.q)
    print(f"[theta] Cramer theta(q={args.q}) = {theta:.5f}")
    print()

    # Model 1: pure Cramer  rate(X) = A * X^(-theta)
    # Model 2: two-term     rate(X) = A * X^(-theta) + B
    # Fit both via Poisson likelihood on (n_conv_v, n_total_v, rate)

    def nll_cramer(params):
        log_A = params[0]
        rate_pred = np.exp(log_A) * X_v ** (-theta)
        rate_pred = np.clip(rate_pred, 1e-30, 1.0)
        # Binomial NLL: -sum(n_conv * log(p) + (n_total - n_conv) * log(1-p))
        ll = (n_conv_v * np.log(rate_pred) +
              (n_total_v - n_conv_v) * np.log(1.0 - rate_pred + 1e-300))
        return -ll.sum()

    def nll_two_term(params):
        log_A = params[0]
        log_B = params[1]
        rate_pred = np.exp(log_A) * X_v ** (-theta) + np.exp(log_B)
        rate_pred = np.clip(rate_pred, 1e-30, 1.0 - 1e-15)
        ll = (n_conv_v * np.log(rate_pred) +
              (n_total_v - n_conv_v) * np.log(1.0 - rate_pred))
        return -ll.sum()

    # Initial guesses
    log_A_init = np.log(rate_v[-3] * X_v[-3] ** theta)  # back out A from large-X bin
    res1 = minimize(nll_cramer, [log_A_init], method="Nelder-Mead", options={"xatol": 1e-6})
    res2 = minimize(nll_two_term, [log_A_init, np.log(1e-10)], method="Nelder-Mead",
                    options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 5000})

    A1 = np.exp(res1.x[0])
    A2 = np.exp(res2.x[0])
    B2 = np.exp(res2.x[1])

    nll1 = res1.fun
    nll2 = res2.fun
    # Likelihood ratio test, df=1
    lr = 2 * (nll1 - nll2)

    print(f"=== Model 1: rate(X) = A * X^(-theta), theta = {theta:.4f} (Cramer) ===")
    print(f"  fitted A = {A1:.4e}    log(A) = {np.log(A1):.4f}")
    print(f"  NLL      = {nll1:.4f}")
    print()
    print(f"=== Model 2: rate(X) = A * X^(-theta) + B  (two-term) ===")
    print(f"  fitted A = {A2:.4e}    log(A) = {np.log(A2):.4f}")
    print(f"  fitted B = {B2:.4e}    log(B) = {np.log(B2):.4f}")
    print(f"  NLL      = {nll2:.4f}")
    print(f"  LR statistic (Model 2 vs Model 1): {lr:.2f}  (df=1)")
    from scipy.stats import chi2
    p_lr = 1 - chi2.cdf(lr, df=1)
    print(f"  p-value for B != 0: {p_lr:.3e}")

    print()
    print(f"=== Per-bin fits (Model 2) ===")
    print(f"{'log(X)':>9} {'n_total':>14} {'n_conv':>9} {'rate_emp':>12} {'rate_pred':>12} {'AX^(-th)':>12} {'B_floor':>11}")
    for i in range(len(log_X_v)):
        rate_pred = A2 * X_v[i] ** (-theta) + B2
        a_term = A2 * X_v[i] ** (-theta)
        print(f"{log_X_v[i]:>9.3f} {n_total_v[i]:>14,} {n_conv_v[i]:>9,} {rate_v[i]:>12.4e} "
              f"{rate_pred:>12.4e} {a_term:>12.4e} {B2:>11.4e}")

    print()
    print(f"=== Verdict ===")
    if p_lr < 0.01:
        print(f"  Two-term model is significantly better (LR p = {p_lr:.2e}).")
        # Where does B dominate vs A*X^(-theta)?
        X_star = (A2 / B2) ** (1.0 / theta)
        print(f"  Crossover X_* (A*X^(-theta) = B):  X = {X_star:.2f}, log(X) = {np.log(X_star):.3f}")
        print(f"  --> Below log(X) ~ {np.log(X_star):.2f}, the floor B dominates.")
        print(f"  --> Above log(X) ~ {np.log(X_star):.2f}, pure Cramer A * X^(-theta) takes over.")
        print(f"  Mechanism: small-X orbits reach 1 via deterministic short paths,")
        print(f"  large-X orbits reach 1 via random-walk descent. The +log(log) signature")
        print(f"  in linear log-log fits is curvature from this additive floor.")
    else:
        print(f"  Two-term model is NOT significantly better than Cramer (LR p = {p_lr:.2e}).")
        print(f"  The +log(log) signature is not explained by a small-X plateau.")


if __name__ == "__main__":
    main()
