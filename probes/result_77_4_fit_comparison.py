"""
R77.4 — Operator-shape fit comparison for |eps_n| * 2^n envelope.

Five hypotheses fitted on n = 2..N (n=1 is transient, excluded):

  H0 (constant):    e_n = c
  H1 (Jordan):      e_n = a + b*n
  H2 (Logarithmic): e_n = a + b*log(n)
  H3 (Power-law):   e_n = c * n^{-alpha}
  H4 (Combined):    e_n = c * n^{-alpha} * (log n)^beta

where e_n := |eps_n| * 2^n.

Outputs:
  - C:/Collatz/result_77_4_envelope_data.csv (table with predictions + residuals)
  - prints AIC / Delta-AIC summary used in the writeup.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import numpy as np

DATA_CSV = Path(r"C:/Collatz/experiments_output/S_k_exact_through_6.csv")
K7_CSV = Path(r"C:/Collatz/experiments_output/result_77_4_S_k_exact_through_7.csv")
OUT_CSV = Path(r"C:/Collatz/result_77_4_envelope_data.csv")


def load_envelope() -> tuple[np.ndarray, np.ndarray, bool]:
    """Return (n_array, |eps|*2^n array, has_k7).

    Reads the through-6 CSV as the base, and if the k=7 result file exists,
    appends the n=7 row.
    """
    rows: list[tuple[int, float]] = []
    with DATA_CSV.open() as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            n = int(r["k"])
            env = float(r["abs_eps_times_2k"])
            rows.append((n, env))

    has_k7 = False
    if K7_CSV.exists():
        with K7_CSV.open() as fh:
            reader = csv.DictReader(fh)
            for r in reader:
                n = int(r["k"])
                env = float(r["abs_eps_times_2k"])
                if n == 7:
                    rows.append((n, env))
                    has_k7 = True

    rows.sort()
    ns = np.array([r[0] for r in rows], dtype=float)
    es = np.array([r[1] for r in rows], dtype=float)
    return ns, es, has_k7


def fit_h0(n: np.ndarray, e: np.ndarray):
    """Constant model e = c."""
    c = float(np.mean(e))
    pred = np.full_like(e, c)
    rss = float(np.sum((e - pred) ** 2))
    k_params = 1
    return {"name": "H0_constant", "params": {"c": c}, "pred": pred, "rss": rss, "k": k_params}


def fit_h1(n: np.ndarray, e: np.ndarray):
    """Linear in n: e = a + b*n."""
    A = np.column_stack([np.ones_like(n), n])
    coef, *_ = np.linalg.lstsq(A, e, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    pred = a + b * n
    rss = float(np.sum((e - pred) ** 2))
    ss_tot = float(np.sum((e - e.mean()) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else float("nan")
    return {"name": "H1_jordan", "params": {"a": a, "b": b}, "pred": pred, "rss": rss, "k": 2, "r2": r2}


def fit_h2(n: np.ndarray, e: np.ndarray):
    """Log model: e = a + b*log(n)."""
    ln = np.log(n)
    A = np.column_stack([np.ones_like(n), ln])
    coef, *_ = np.linalg.lstsq(A, e, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    pred = a + b * ln
    rss = float(np.sum((e - pred) ** 2))
    ss_tot = float(np.sum((e - e.mean()) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else float("nan")
    return {"name": "H2_log", "params": {"a": a, "b": b}, "pred": pred, "rss": rss, "k": 2, "r2": r2}


def fit_h3(n: np.ndarray, e: np.ndarray):
    """Power law: e = c * n^{-alpha}, fit via log-log linear regression."""
    ln = np.log(n)
    le = np.log(e)
    A = np.column_stack([np.ones_like(n), ln])
    coef, *_ = np.linalg.lstsq(A, le, rcond=None)
    log_c, neg_alpha = float(coef[0]), float(coef[1])
    c = math.exp(log_c)
    alpha = -neg_alpha  # since regression slope is -alpha
    pred = c * n ** (-alpha)
    rss = float(np.sum((e - pred) ** 2))  # original-scale RSS
    # log-scale R^2 (for the fit quality assessment)
    le_pred = log_c + neg_alpha * ln
    rss_log = float(np.sum((le - le_pred) ** 2))
    ss_tot_log = float(np.sum((le - le.mean()) ** 2))
    r2_log = 1.0 - rss_log / ss_tot_log if ss_tot_log > 0 else float("nan")
    return {"name": "H3_power", "params": {"c": c, "alpha": alpha}, "pred": pred, "rss": rss, "k": 2, "r2_log": r2_log}


def fit_h4(n: np.ndarray, e: np.ndarray):
    """Combined: e = c * n^{-alpha} * (log n)^beta.

    Log-linearizes to log(e) = log(c) - alpha * log(n) + beta * log(log(n)).
    Requires log(log(n)) > 0 i.e. n > e (~2.718). For n=2 log(log(2)) is
    negative but defined; that's fine for the regression.
    """
    ln = np.log(n)
    lln = np.log(ln)
    le = np.log(e)
    A = np.column_stack([np.ones_like(n), ln, lln])
    coef, *_ = np.linalg.lstsq(A, le, rcond=None)
    log_c, neg_alpha, beta = float(coef[0]), float(coef[1]), float(coef[2])
    c = math.exp(log_c)
    alpha = -neg_alpha
    pred = c * n ** (-alpha) * ln ** beta
    rss = float(np.sum((e - pred) ** 2))
    le_pred = log_c + neg_alpha * ln + beta * lln
    rss_log = float(np.sum((le - le_pred) ** 2))
    ss_tot_log = float(np.sum((le - le.mean()) ** 2))
    r2_log = 1.0 - rss_log / ss_tot_log if ss_tot_log > 0 else float("nan")
    return {
        "name": "H4_combined",
        "params": {"c": c, "alpha": alpha, "beta": beta},
        "pred": pred,
        "rss": rss,
        "k": 3,
        "r2_log": r2_log,
    }


def aic(rss: float, n_pts: int, k_params: int) -> float:
    """AIC for least-squares with Gaussian errors (up to additive constant)."""
    return n_pts * math.log(rss / n_pts) + 2 * k_params


def main():
    ns_all, es_all, has_k7 = load_envelope()
    print(f"Loaded {len(ns_all)} points (n = {ns_all.tolist()}); has_k7 = {has_k7}")

    # Fit on n >= 2 (drop transient n=1)
    mask = ns_all >= 2
    n = ns_all[mask]
    e = es_all[mask]
    n_pts = len(n)
    print(f"Fitting on n = {n.astype(int).tolist()}  (n_pts={n_pts})")
    print(f"|eps|*2^n values: {e}")
    print()

    fits = [fit_h0(n, e), fit_h1(n, e), fit_h2(n, e), fit_h3(n, e), fit_h4(n, e)]

    # AIC
    print("Fits:")
    for f in fits:
        f["aic"] = aic(f["rss"], n_pts, f["k"])
        print(f"  {f['name']:14s}  params={f['params']}  RSS={f['rss']:.3e}  k={f['k']}  AIC={f['aic']:.3f}")

    aic_min = min(f["aic"] for f in fits)
    print()
    print(f"Min AIC = {aic_min:.3f}")
    print("Delta AIC (relative to best):")
    for f in fits:
        f["dAIC"] = f["aic"] - aic_min
        print(f"  {f['name']:14s}  Delta AIC = {f['dAIC']:+.3f}")

    # Residuals — check sign patterns
    print()
    print("Residuals (predicted - observed) per model:")
    print("  n :    obs   ", "   ".join(f["name"] for f in fits))
    for i, ni in enumerate(n):
        row = f"  {int(ni)}: {e[i]:.6f}  "
        for f in fits:
            row += f"  {f['pred'][i] - e[i]:+.4e}"
        print(row)

    # Sign pattern (for trend-vs-random check)
    print()
    print("Residual sign sequences (predicted - observed):")
    for f in fits:
        signs = "".join("+" if r > 0 else ("-" if r < 0 else "0") for r in (f["pred"] - e))
        print(f"  {f['name']:14s}  {signs}")

    # H1 direction check
    h1 = next(f for f in fits if f["name"] == "H1_jordan")
    b = h1["params"]["b"]
    print()
    print(f"H1 direction check: b = {b:+.6e}.  b > 0 (Jordan-consistent)? {'YES' if b > 0 else 'NO'}")

    # Write CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.writer(fh)
        header = ["n", "abs_eps_times_2n"] + [f"{f['name']}_pred" for f in fits] + [f"{f['name']}_resid" for f in fits]
        w.writerow(header)
        for i, ni in enumerate(n):
            row = [int(ni), e[i]]
            row.extend(f["pred"][i] for f in fits)
            row.extend(f["pred"][i] - e[i] for f in fits)
            w.writerow(row)
    print(f"\nWrote {OUT_CSV}")

    # Summary block for the writeup (machine-parseable)
    print()
    print("=== SUMMARY ===")
    print(f"n_pts = {n_pts}")
    print(f"has_k7 = {has_k7}")
    for f in fits:
        params = ", ".join(f"{k}={v:+.6e}" for k, v in f["params"].items())
        print(f"{f['name']}: {params}; RSS={f['rss']:.4e}; AIC={f['aic']:.3f}; dAIC={f['dAIC']:+.3f}")


if __name__ == "__main__":
    main()
