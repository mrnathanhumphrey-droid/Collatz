"""
Experiment 08 — All-n version of the prefix decomposition

The structural decomposition was developed for ODD n only. Bonacorsi-Bordoni
model all n in [1, 10^7]. This script extends the framework to all n.

Approach:
  - Symbolic prefix algorithm (state = a*m + c) extends trivially to even
    residues r mod 2^k: for r=0, state halves down to (a=1, c=0) without
    encountering an odd state during the prefix.
  - Post-prefix state = a_final * m + c_final with a_final odd.
  - For all-n (m takes integer values, half even half odd), state parity
    depends on m, so post-prefix trajectory follows the all-n heuristic
    slope 6.95 = 2/(ln 4 - ln 3) rather than the odd-only slope 10.43.

Caveat: residues r=0 mod 2^k host n with variable v_2(n) (initial halvings),
so the prefix length itself is not deterministic for that class. We exclude
r=0 from the per-class regressions and treat it separately.

We test:
  - Per-class OLS slope: is it universal across r mod 8 (excluding r=0)?
    If yes, the slope universality theorem extends.
  - Per-class OLS intercept α: is it predicted by α_det = prefix_steps +
    6.95 * ln(a_final / 2^k)?

Usage:
    python 08_all_n_decomposition.py --N 10000000 --k 3
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))   # ~10.4282
# All-n slope is the same 10.4282 because sigma(n) = nu_2(n) + sigma(odd_part(n))
# and slope of sigma vs ln(n) is governed by the odd-part trajectory.
# nu_2(n) contributes only an additive offset (mean ~1 for random integers).
LOG_FACTOR_ALL = LOG_FACTOR_ODD


def deterministic_prefix(r, a0, max_steps=400):
    """Symbolic prefix from state (a0, r). Stop when a is odd."""
    a, c = a0, r; steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3*c + 1
        steps += 1
    return steps, a, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--k", type=int, default=3, help="Modular resolution (default 3 = mod 8 to match B&B)")
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    M = 1 << args.k  # mod
    print(f"[load] N={args.N:,}, mod {M}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(pl.col("n") >= 2)
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    res = (n % M).astype(np.int32)
    print(f"        {len(n):,} rows total")

    # Pooled
    b_p, a_p = np.polyfit(log_n, sigma, 1)
    print(f"[pool] all-n pooled OLS slope = {b_p:.4f}  (heuristic for all n: {LOG_FACTOR_ALL:.4f})")
    print(f"       intercept = {a_p:.4f}")
    print()

    # Per-residue OLS
    print(f"{'r':>3} {'count':>10} {'slope':>9} {'intercept':>10} "
          f"{'prefix':>7} {'a_final':>8} {'c_final':>8} {'pred_a':>9}")
    rows = []
    for r in range(M):
        m_mask = res == r
        n_class = m_mask.sum()
        if n_class < 50:
            continue
        bk, ak = np.polyfit(log_n[m_mask], sigma[m_mask], 1)
        steps, a_f, c_f = deterministic_prefix(r, M)
        # Predicted intercept: with all-n heuristic for post-prefix
        pred_a = steps + LOG_FACTOR_ALL * np.log(a_f / float(M))
        # For r=0 case where a_final=1, the post-prefix is just m+0 with all-n heuristic.
        # But r=0 has variable v_2 so the prefix length itself is not constant — we
        # report whatever the smallest deterministic prefix gives.
        rows.append((r, n_class, bk, ak, steps, a_f, c_f, pred_a))
        print(f"{r:>3} {n_class:>10,} {bk:>9.4f} {ak:>10.3f} "
              f"{steps:>7} {a_f:>8} {c_f:>8} {pred_a:>9.3f}")

    rows = np.array(rows, dtype=[
        ("r", "i4"), ("count", "i8"), ("slope", "f8"), ("intercept", "f8"),
        ("prefix", "i4"), ("a_final", "i8"), ("c_final", "i8"), ("pred_a", "f8")
    ])

    # Slope universality across residues (excluding r=0 which has variable v_2)
    rows_nz = rows[rows["r"] != 0]
    print()
    print(f"=== Slope across residues (excluding r=0, which has variable v_2) ===")
    print(f"  mean slope         = {rows_nz['slope'].mean():.4f}")
    print(f"  SD slope            = {rows_nz['slope'].std():.4f}")
    print(f"  range               = [{rows_nz['slope'].min():.4f}, {rows_nz['slope'].max():.4f}]")
    print(f"  heuristic (all n)   = {LOG_FACTOR_ALL:.4f}")
    print(f"  heuristic (odd n)   = {LOG_FACTOR_ODD:.4f}")
    if rows[rows["r"] == 0].size:
        r0 = rows[rows["r"] == 0][0]
        print(f"  r=0 slope (anomaly) = {r0['slope']:.4f}  (depends on v_2(m), variable)")

    # Intercept decomposition
    print()
    print(f"=== Intercept decomposition: alpha_actual vs alpha_det ===")
    print(f"  using all-n heuristic constant = {LOG_FACTOR_ALL:.4f}")
    if len(rows_nz) >= 2:
        b_fit, a_fit = np.polyfit(rows_nz["pred_a"], rows_nz["intercept"], 1)
        fitted = a_fit + b_fit * rows_nz["pred_a"]
        resid = rows_nz["intercept"] - fitted
        ss_total = np.sum((rows_nz["intercept"] - rows_nz["intercept"].mean())**2)
        ss_resid = np.sum(resid**2)
        r_squared = 1 - ss_resid / ss_total if ss_total > 0 else float("nan")
        print(f"  fit: alpha_actual = {a_fit:.3f} + {b_fit:.4f} * alpha_pred")
        print(f"  R^2 = {r_squared:.4f}")
        print(f"  SD(intercept) across r = {rows_nz['intercept'].std():.3f}")
        print(f"  SD(residual) = {resid.std():.3f}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.bar(rows["r"], rows["slope"], color="C0", alpha=0.7, edgecolor="black")
    ax.axhline(LOG_FACTOR_ALL, color="red", linestyle="--", alpha=0.7,
               label=f"all-n heuristic = {LOG_FACTOR_ALL:.3f}")
    ax.axhline(LOG_FACTOR_ODD, color="green", linestyle=":", alpha=0.6,
               label=f"odd-n heuristic = {LOG_FACTOR_ODD:.3f}")
    ax.set_xlabel("residue r mod 2^k")
    ax.set_ylabel("per-class OLS slope")
    ax.set_title(f"All-n slope per class  (k={args.k}, mod {M}, N={args.N:,})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    if len(rows_nz) >= 2:
        ax.scatter(rows_nz["pred_a"], rows_nz["intercept"], s=80, color="C0", alpha=0.8)
        for r in rows_nz:
            ax.annotate(f"r={r['r']}", (r["pred_a"], r["intercept"]),
                        xytext=(5,5), textcoords="offset points", fontsize=8)
        xs = np.linspace(rows_nz["pred_a"].min(), rows_nz["pred_a"].max(), 100)
        ax.plot(xs, a_fit + b_fit*xs, "r--", lw=1.2,
                label=f"R^2 = {r_squared:.4f}")
        ax.set_xlabel("predicted alpha (prefix + 6.95 ln(a_final/M))")
        ax.set_ylabel("alpha actual (per-class OLS intercept)")
        ax.set_title(f"Intercept vs prefix prediction (all-n, r ≠ 0)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"All-n prefix decomposition  (k={args.k}, N={args.N:,})", y=1.02)
    plt.tight_layout()
    out_png = out_dir / f"08_all_n_decomposition_k{args.k}_N{args.N}.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close()
    print()
    print(f"[save] {out_png}")


if __name__ == "__main__":
    main()
