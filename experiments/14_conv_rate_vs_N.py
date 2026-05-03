"""
Experiment 14 — Convergence rate as a function of N for q in {5, 7}.

At q=5, conv_rate decays gently with N. At q=7, it decays steeply.
Quantify: P(converge | n in [1,N]) ~ N^(-alpha(q)) and find alpha for each q.

For each q, fit log(conv_rate) ~ log(N) on the points already generated and
report the decay exponent alpha. Then compare to potential analytical forms:
  - alpha = log(q/4) / log(2)  (rough scaling from drift)
  - alpha = some other clean q-dependent form

Usage:
    python 14_conv_rate_vs_N.py
"""
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import linregress

sys.stdout.reconfigure(encoding="utf-8")


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    runs = [
        # (q, N) — only those we've actually generated
        (5, 1_000_000),
        (5, 10_000_000),
        (5, 100_000_000),
        (7, 1_000_000),
        (7, 10_000_000),
        (7, 100_000_000),
        (3, 1_000_000),
        (3, 10_000_000),
    ]

    rows = []
    for q, N in runs:
        path = data_dir / f"q_main_q{q}_N{N}.parquet"
        if not path.exists():
            continue
        df = pl.read_parquet(path)
        n_total = len(df)
        n_conv = int(df["converged"].sum())
        rate = n_conv / max(n_total, 1)
        rows.append((q, N, n_total, n_conv, rate))

    # Pretty print
    print(f"\n=== conv_rate vs N (across q) ===")
    print(f"{'q':>3} {'N':>12} {'n_total':>14} {'n_conv':>10} {'conv_rate':>14}")
    for q, N, n_tot, n_c, r in rows:
        print(f"{q:>3} {N:>12,} {n_tot:>14,} {n_c:>10,} {r:>14.6e}")

    # For each q with >= 2 N values, fit log(rate) vs log(N)
    print()
    print(f"=== Decay exponent alpha (conv_rate ~ N^-alpha) per q ===")
    print(f"{'q':>3} {'#points':>9} {'alpha':>10} {'intercept':>12} {'R^2':>8} {'p':>10}")

    by_q = {}
    for q, N, n_tot, n_c, r in rows:
        if r > 0:
            by_q.setdefault(q, []).append((N, r))

    fits = {}
    for q, pts in by_q.items():
        if len(pts) < 2:
            continue
        Ns = np.array([p[0] for p in pts], dtype=float)
        rates = np.array([p[1] for p in pts])
        log_N = np.log(Ns)
        log_r = np.log(rates)
        if len(pts) == 2:
            # exact fit
            slope = (log_r[1] - log_r[0]) / (log_N[1] - log_N[0])
            intercept = log_r[0] - slope * log_N[0]
            r_sq = 1.0
            p_val = float("nan")
        else:
            res = linregress(log_N, log_r)
            slope, intercept, r_sq, p_val = res.slope, res.intercept, res.rvalue ** 2, res.pvalue
        alpha = -slope
        fits[q] = (alpha, intercept, r_sq, p_val)
        p_str = f"{p_val:.3e}" if not np.isnan(p_val) else "       --"
        print(f"{q:>3} {len(pts):>9} {alpha:>10.4f} {intercept:>12.4f} {r_sq:>8.4f} {p_str:>10}")

    print()
    print(f"=== Comparison to candidate analytical forms ===")
    print(f"  Candidate A: alpha = log(q/4) / log(2)")
    print(f"  Candidate B: alpha = log(q/4) / log(2) * some C")
    print(f"  Candidate C: alpha = (log(q) - log(4)) / log(2) (= log_2(q/4))")
    print()
    print(f"{'q':>3} {'alpha_emp':>11} {'log(q/4)':>10} {'log_2(q/4)':>12} {'alpha/log(q/4)':>15} {'alpha/log_2(q/4)':>17}")
    for q, (alpha, _, _, _) in fits.items():
        if q == 3:
            print(f"{q:>3} {alpha:>11.4f} {np.log(q/4):>10.4f} {np.log2(q/4):>12.4f}  (q=3 doesn't decay; orbits all converge)")
            continue
        log_qfact = np.log(q/4)
        log2_qfact = np.log2(q/4)
        if abs(log_qfact) > 1e-9:
            ratio_log = alpha / log_qfact
            ratio_log2 = alpha / log2_qfact
            print(f"{q:>3} {alpha:>11.4f} {log_qfact:>10.4f} {log2_qfact:>12.4f} {ratio_log:>15.4f} {ratio_log2:>17.4f}")

    # Save CSV
    out_csv = out_dir / f"14_conv_rate_vs_N.csv"
    pl.DataFrame([
        {"q": q, "N": N, "n_total": nt, "n_conv": nc, "conv_rate": r}
        for q, N, nt, nc, r in rows
    ]).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
