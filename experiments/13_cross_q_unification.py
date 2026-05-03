"""
Experiment 13 — Cross-q unification of the prefix complexity axis.

For q=3 (Collatz, 100% convergence), high prefix complexity (large odd_steps_in_prefix)
predicts LONGER sigma — the orbit takes more total steps to reach 1.
For q=5 (mostly divergent), high prefix complexity predicts LOWER convergence rate.

Both findings should be different manifestations of the same structural axis:
prefix complexity = how many multiplicative growth applications get baked into
the trajectory before the prefix terminates.

Unification target: compute the rate at which mean(sigma) grows per unit j_odd_steps
at q=3, compare to the rate at which log(conv_rate) decays per unit j at q=5.

Procedure:
  q=3 data: existing N=10^7 parquet (main_N10000000.parquet, all-n).
  - For each odd residue class mod 64 (k=6), compute mean(sigma) restricted
    to convergent orbits (which is all of them at q=3).
  - Subtract heuristic baseline 10.4282 * mean(log_n) per class to isolate
    the prefix contribution.
  - Regress the residual mean against odd_steps_in_prefix.
  - That slope = additional sigma per unit j at q=3.

Compare against q=5 slope from experiment 12.

Usage:
    python 13_cross_q_unification.py --N 10000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import linregress

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))  # ~10.4282


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
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--k", type=int, default=6)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    M = 1 << args.k
    K = M // 2  # number of odd residue classes

    print(f"[load] q=3 (Collatz), N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)
    print(f"        odd-only rows = {len(n):,}", flush=True)

    # Per-class odd_steps lookup
    odd_steps_per_class = np.zeros(K, dtype=np.int64)
    a_final_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r_class = 2 * kk + 1
        a, c, os_, es_ = deterministic_prefix_q(r_class, 3, M)
        odd_steps_per_class[kk] = os_
        a_final_per[kk] = a

    # Per-class mean sigma residual (after subtracting heuristic baseline)
    mean_sigma_residual = np.zeros(K)
    n_per_class = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        m = class_idx == kk
        if m.sum() == 0:
            continue
        n_per_class[kk] = m.sum()
        mean_sigma_residual[kk] = (sigma[m] - LOG_FACTOR_ODD * log_n[m]).mean()

    # Group by odd_steps value
    j_values = np.unique(odd_steps_per_class)
    rows = []
    for j in j_values:
        mask = odd_steps_per_class == j
        n_classes = int(mask.sum())
        if n_classes == 0:
            continue
        # Average per-class residual across classes with this j
        # Weight by class size (which is roughly equal at fixed N, k anyway)
        total_n = int(n_per_class[mask].sum())
        # Pool: combine all n in classes with this j and compute (sigma - heuristic*log_n).mean()
        m_all = np.isin(class_idx, np.where(mask)[0])
        pooled_residual = float((sigma[m_all] - LOG_FACTOR_ODD * log_n[m_all]).mean())
        rows.append((int(j), n_classes, total_n, pooled_residual))

    print()
    print(f"=== q=3: mean sigma residual (sigma - 10.43*log_n) by prefix odd-step count ===")
    print(f"{'j':>4} {'#classes':>9} {'n_total':>12} {'mean_residual':>14}")
    for j, n_cls, n_tot, resid in rows:
        print(f"{j:>4} {n_cls:>9} {n_tot:>12,} {resid:>14.4f}")

    js = np.array([r[0] for r in rows], dtype=np.float64)
    residuals = np.array([r[3] for r in rows])
    if len(js) >= 3:
        res_lr = linregress(js, residuals)
        print()
        print(f"=== q=3 unification slope ===")
        print(f"  Slope of mean(sigma_residual) vs j_odd_steps: {res_lr.slope:>8.4f}")
        print(f"  Intercept: {res_lr.intercept:>8.4f}")
        print(f"  R^2: {res_lr.rvalue**2:.4f}")
        print(f"  p: {res_lr.pvalue:.3e}")

        print()
        print(f"=== Comparison to q=5 ===")
        # q=5 emp slope of log(conv_rate) vs j was -0.5619 (k=8 N=10^8)
        # That's per unit j: log(conv_rate) drops by 0.5619 per +1 j
        # Per unit log(q/4): drops by 0.5619/log(5/4) = 2.518
        q5_emp_slope = 0.5619
        q5_emp_per_logq4 = q5_emp_slope / np.log(5/4)
        # For q=3: each odd-step in prefix adds slope_q3 to mean sigma_residual.
        # Translate to "per unit log(q/4)": for q=3, log(3/4) = -0.2877 (negative!)
        # So per +1 j: residual increases by slope_q3
        # The unifying axis: rate per unit prefix complexity scales differently per q
        # because the underlying drift is different sign.
        slope_q3 = res_lr.slope
        log_qfact_q3 = np.log(3/4)
        ratio_q3 = slope_q3 / log_qfact_q3  # additional sigma per unit log(3/4)
        print(f"  q=3 slope (sigma residual per +1 j):       {slope_q3:>9.4f}")
        print(f"  q=3 slope / log(3/4):                       {ratio_q3:>9.4f}")
        print(f"  q=5 -log(conv_rate) slope (decrement per +1 j): {q5_emp_slope:>9.4f}")
        print(f"  q=5 (-slope) / log(5/4):                    {q5_emp_per_logq4:>9.4f}")
        print()
        print(f"  Note: q=3 has log(3/4) < 0 (orbits shrink on average),")
        print(f"        q=5 has log(5/4) > 0 (orbits grow on average).")
        print(f"        The same prefix-complexity axis modulates BOTH the rate of shrinking")
        print(f"        (longer sigma at q=3) AND the rate of growing-out-of-convergence (lower")
        print(f"        conv_rate at q=5).")

    # Save CSV
    out_csv = out_dir / f"13_cross_q_unification_k{args.k}_N{args.N}.csv"
    pl.DataFrame([
        {"j_odd_steps": r[0], "n_classes": r[1], "n_total": r[2], "mean_sigma_residual": r[3]}
        for r in rows
    ]).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
