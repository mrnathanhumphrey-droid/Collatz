"""
Experiment 12 — Analytical derivation check for qx+1 convergence rate.

Conjectured form (random-walk heuristic):
  Post-prefix orbit walks log(value) with:
    drift     mu = log(q/4)
    variance  sigma^2 ~ 2 * log(2)^2   (from Geometric(1/2) v_2)
  Starting log-value V_j ~ j * log(q) for prefix odd-step count j.
  Gambler's ruin / large-deviation: P(return to 1) ~ exp(-V_j * 2mu/sigma^2)
                                                  = exp(-j * log(q) * mu/log(2)^2)

So log(conv_rate(j)) = const - j * log(q) * log(q/4) / log(2)^2.

Predicted slopes of log(conv_rate) vs j (per unit increment in odd_steps):
  q=5: -log(5)*log(5/4)/log(2)^2 ~= -0.7466
  q=7: -log(7)*log(7/4)/log(2)^2 ~= -2.2671

If empirical slopes match these (modulo a fittable additive constant from the
prefactor), the relationship has a closed form derivable from first principles.

Empirical procedure:
  1. Load qx+1 data at chosen N
  2. For each odd n, compute its odd_steps_in_prefix from its residue mod 2^k
  3. Group n by odd_steps_in_prefix
  4. Per group: compute conv_rate = (# converged) / (# total)
  5. Linear regression of log(conv_rate) vs odd_steps
  6. Compare empirical slope to predicted slope; report ratio.

Usage:
    python 12_q_convrate_analytical.py --q 5 --N 100000000 --k 8
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import linregress

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix_q(r, q, a0, max_steps=400):
    a, c = a0, r
    odd_s = 0; even_s = 0
    while a % 2 == 0 and (odd_s + even_s) < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
            even_s += 1
        else:
            a *= q; c = q*c + 1
            odd_s += 1
    return a, c, odd_s, even_s


def predicted_slope(q):
    """Random-walk prediction for slope of log(conv_rate) vs odd_steps."""
    return -np.log(q) * np.log(q / 4.0) / (np.log(2.0) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--N", type=int, default=100_000_000)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] q={args.q}, N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet")
    n = df["n"].to_numpy().astype(np.int64)
    converged = df["converged"].to_numpy().astype(bool)
    M = 1 << args.k
    res = (n % M).astype(np.int32)

    # Per-class odd_steps lookup (only for ODD residues r = 2k+1)
    K = M // 2
    odd_steps_per_class = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r_class = 2 * kk + 1
        _, _, os_, _ = deterministic_prefix_q(r_class, args.q, M)
        odd_steps_per_class[kk] = os_

    # Map each n's residue to its odd_steps_in_prefix
    class_idx = ((res - 1) // 2).astype(np.int32)
    odd_steps_per_n = odd_steps_per_class[class_idx]

    # Group by odd_steps value, compute convergence rate
    j_values = np.unique(odd_steps_per_class)
    rows = []
    for j in j_values:
        mask = odd_steps_per_n == j
        n_total = mask.sum()
        n_conv = (mask & converged).sum()
        if n_total == 0:
            continue
        rate = n_conv / n_total
        # Number of residue classes with this j
        n_classes = (odd_steps_per_class == j).sum()
        rows.append((int(j), int(n_classes), int(n_total), int(n_conv), float(rate)))

    print()
    print(f"=== conv_rate by prefix odd-step count (q={args.q}, k={args.k}, N={args.N:,}) ===")
    print(f"{'j':>4} {'#classes':>9} {'n_total':>12} {'n_conv':>9} {'conv_rate':>12} {'log_rate':>11}")
    for j, n_cls, n_tot, n_c, rate in rows:
        log_rate_str = f"{np.log(rate):.4f}" if rate > 0 else "      -inf"
        print(f"{j:>4} {n_cls:>9} {n_tot:>12,} {n_c:>9,} {rate:>12.6e} {log_rate_str:>11}")

    # Linear regression of log(conv_rate) vs j (only points with conv_rate > 0)
    rows_pos = [r for r in rows if r[4] > 0]
    if len(rows_pos) < 3:
        print("\n[warn] fewer than 3 j-values with conv_rate > 0; cannot fit slope")
        return

    js = np.array([r[0] for r in rows_pos], dtype=np.float64)
    log_rates = np.array([np.log(r[4]) for r in rows_pos])

    # Weight regression by inverse variance of log_rate ~ 1/(n_conv) for binomial
    weights = np.array([r[3] for r in rows_pos], dtype=np.float64)  # n_conv as weights

    # Unweighted OLS first (simple)
    res_u = linregress(js, log_rates)
    # Weighted OLS: minimize sum(w_i * (y_i - a - b*x_i)^2)
    sum_w = weights.sum()
    x_w = (weights * js).sum() / sum_w
    y_w = (weights * log_rates).sum() / sum_w
    sxx = (weights * (js - x_w) ** 2).sum()
    sxy = (weights * (js - x_w) * (log_rates - y_w)).sum()
    b_w = sxy / sxx
    a_w = y_w - b_w * x_w

    pred_slope = predicted_slope(args.q)

    print()
    print(f"=== Slope analysis ===")
    print(f"  Unweighted OLS:  slope = {res_u.slope:>9.4f}  intercept = {res_u.intercept:>9.4f}  R^2 = {res_u.rvalue**2:.4f}  p = {res_u.pvalue:.3e}")
    print(f"  Weighted OLS:    slope = {b_w:>9.4f}  intercept = {a_w:>9.4f}  (weights = n_conv per j-group)")
    print(f"  Predicted slope: {pred_slope:>9.4f}  (random-walk drift={np.log(args.q/4):.4f}, sigma^2 ~ {2*np.log(2)**2:.4f})")
    print(f"  Ratio empirical/predicted (unweighted): {res_u.slope / pred_slope:.4f}")
    print(f"  Ratio empirical/predicted (weighted):   {b_w / pred_slope:.4f}")

    # Save CSV
    out_csv = out_dir / f"12_q_convrate_analytical_q{args.q}_k{args.k}_N{args.N}.csv"
    df_out = pl.DataFrame([
        {"j_odd_steps": r[0], "n_classes": r[1], "n_total": r[2], "n_conv": r[3], "conv_rate": r[4]}
        for r in rows
    ])
    df_out.write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
