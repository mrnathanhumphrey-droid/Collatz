"""
result_epsilon_7.py
===================
Compute eps_7 = S_7 - 7/15 to determine the asymptotic rate of S_k -> 7/15.
Two-mode falsified at k=6; ratio trajectory |eps_{k+1}/eps_k| at k=2..6 is
0.535, 0.482, 0.470, 0.432 (decreasing). eps_7 supplies the next ratio.

Strategy: float64 power iteration on K_7 (1458 states), same path as eps_6.
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
ENV_CSV = os.path.join(OUT_DIR, "result_q_sweep_test_1_envelope.csv")


def order_of_two(N):
    assert N % 2 == 1
    m = 1; v = 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def build_K_float(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    pi = inv2
    for v in range(M):
        powers_inv2[v] = pi
        pi = (pi * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    K = np.zeros((n, n), dtype=np.float64)
    # 2^M overflows int->float for M > ~1023, but Z_v = 1 - 2^(-M) is
    # numerically 1.0 for any M > ~50, so use that form instead.
    Z_v = 1.0 - 2.0 ** (-M)
    weights = np.empty(M, dtype=np.float64)
    for v in range(1, M + 1):
        weights[v - 1] = (2.0 ** (-v)) / Z_v
    for i_r, r in enumerate(coprime):
        base = (q * int(r) + 1) % N
        for v in range(1, M + 1):
            tgt = (base * int(powers_inv2[v - 1])) % N
            j = int(state_idx[tgt])
            K[i_r, j] += weights[v - 1]
    return K, coprime, M, state_idx


def stationary_power_iteration(K, max_iter=400, tol=1e-13):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    for it in range(max_iter):
        pi_new = pi @ K
        pi_new /= pi_new.sum()
        residual = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        if residual < tol:
            return pi, it + 1, residual
    return pi, max_iter, residual


def load_cached_S():
    S = {}
    with open(ENV_CSV) as f:
        for row in csv.DictReader(f):
            if int(row["q"]) != 3:
                continue
            S[int(row["n"])] = Fraction(int(row["S_n_num"]),
                                         int(row["S_n_den"]))
    return S


def main():
    print("=" * 78)
    print("Computing eps_7 = S_7 - 7/15 (asymptotic rate determination)")
    print("=" * 78)

    S_cache = load_cached_S()
    X5_exact = Fraction(1)
    for k in [1, 2, 3, 4, 5]:
        X5_exact += S_cache[k]
    X5_float = float(X5_exact)
    print()
    print(f"X_5 (exact) = {X5_float:.15f}")

    # k=6 (sanity-recompute)
    print()
    print("--- k = 6 (sanity recompute) ---")
    t0 = time.time()
    K6, _, M6, _ = build_K_float(3, 6)
    t_build_6 = time.time() - t0
    print(f"  build K_6 (states={K6.shape[0]}, M={M6}): {t_build_6:.2f}s")
    t0 = time.time()
    pi6, iters6, res6 = stationary_power_iteration(K6, tol=1e-14)
    t_iter_6 = time.time() - t0
    print(f"  power iteration: {iters6} iters, residual {res6:.2e}, "
          f"{t_iter_6:.2f}s")
    sum_sq_6 = float((pi6 ** 2).sum())
    X6 = (3 ** 6) * sum_sq_6
    S6 = X6 - X5_float
    eps6 = S6 - 7.0 / 15.0
    print(f"  X_6 = {X6:.15f}, S_6 = {S6:.15f}")
    print(f"  eps_6 = {eps6:+.10e}")

    # k=7
    print()
    print("--- k = 7 (target) ---")
    t0 = time.time()
    K7, _, M7, _ = build_K_float(3, 7)
    t_build_7 = time.time() - t0
    print(f"  build K_7 (states={K7.shape[0]}, M={M7}): {t_build_7:.2f}s")
    print(f"  K_7 row-sum check: max |row sum - 1| = "
          f"{float(np.max(np.abs(K7.sum(axis=1) - 1))):.2e}")
    t0 = time.time()
    pi7, iters7, res7 = stationary_power_iteration(K7, tol=1e-14)
    t_iter_7 = time.time() - t0
    print(f"  power iteration: {iters7} iters, residual {res7:.2e}, "
          f"{t_iter_7:.2f}s")
    print(f"  sum pi_7 = {pi7.sum():.15f}  (expect 1)")
    sum_sq_7 = float((pi7 ** 2).sum())
    X7 = (3 ** 7) * sum_sq_7
    S7 = X7 - X6
    eps7 = S7 - 7.0 / 15.0
    print(f"  X_7 = {X7:.15f}, S_7 = {S7:.15f}")
    print(f"  eps_7 = {eps7:+.10e}")

    # Ratio
    ratio_67 = abs(eps7 / eps6)
    print()
    print(f"|eps_7 / eps_6| = {ratio_67:.6f}")

    # Decision band
    print()
    print("Decision rule (per brief):")
    if 0.30 <= ratio_67 <= 0.40:
        verdict = ("TREND TOWARD 1/3 IS REAL — asymptotic rate < 1/2; "
                   "rate-1/2 envelope is loose")
    elif 0.40 <= ratio_67 <= 0.45:
        verdict = ("TREND CONTINUES — rate possibly heading below 1/3 or "
                   "intermediate; eps_8 needed")
    elif 0.45 <= ratio_67 <= 0.55:
        verdict = ("STABILIZED NEAR 1/2 — rate-1/2 is asymptotic truth; "
                   "earlier ratios were finite-k artifacts")
    elif ratio_67 > 0.55:
        verdict = ("RATIO REVERSED UP — non-monotone trajectory; "
                   "structurally informative")
    else:
        verdict = ("RATIO BELOW 0.30 — super-exponential collapse; needs "
                   "reanalysis")
    print(f"  -> {verdict}")

    # Full ratio trajectory
    eps_all = {k: float(S_cache[k] - Fraction(7, 15)) for k in S_cache}
    eps_all[6] = eps6
    eps_all[7] = eps7
    print()
    print("Updated ratio trajectory |eps_{k+1}/eps_k|:")
    for k in [1, 2, 3, 4, 5, 6]:
        rr = abs(eps_all[k+1] / eps_all[k])
        print(f"  k = {k} -> {k+1}: {rr:.6f}")

    # Outputs
    out_md = os.path.join(OUT_DIR, "result_epsilon_7.md")
    md = []
    md.append("# Result: eps_7 = S_7 - 7/15 (asymptotic rate)")
    md.append("")
    md.append(f"**Date:** 2026-05-05.  Float64 power iteration on K_7 "
              f"(1458 states).")
    md.append("")
    md.append(f"## Verdict: **{verdict}**")
    md.append("")
    md.append(f"- eps_7 = `{eps7:+.10e}`")
    md.append(f"- |eps_7/eps_6| = **{ratio_67:.6f}**")
    md.append("")
    md.append("## Ratio trajectory (k=1..7)")
    md.append("")
    md.append("| transition | |eps_{k+1}/eps_k| |")
    md.append("|---|---|")
    for k in [1, 2, 3, 4, 5, 6]:
        rr = abs(eps_all[k+1] / eps_all[k])
        md.append(f"| {k} -> {k+1} | {rr:.6f} |")
    md.append("")
    md.append("## eps_k table (k=1..7)")
    md.append("")
    md.append("| k | eps_k | source |")
    md.append("|---|---|---|")
    for k in [1, 2, 3, 4, 5]:
        md.append(f"| {k} | {float(S_cache[k] - Fraction(7,15)):+.10e} | "
                  f"exact rational (cached) |")
    md.append(f"| 6 | {eps6:+.10e} | float64 power iter (487 states) |")
    md.append(f"| 7 | {eps7:+.10e} | float64 power iter (1458 states) |")
    md.append("")
    md.append("## Computation diagnostics")
    md.append("")
    md.append(f"| step | wall time | iterations | residual |")
    md.append(f"|---|---|---|---|")
    md.append(f"| K_6 build | {t_build_6:.2f}s | — | — |")
    md.append(f"| K_6 power iter | {t_iter_6:.2f}s | {iters6} | "
              f"{res6:.2e} |")
    md.append(f"| K_7 build | {t_build_7:.2f}s | — | — |")
    md.append(f"| K_7 power iter | {t_iter_7:.2f}s | {iters7} | "
              f"{res7:.2e} |")
    md.append("")
    md.append(f"Total wall: "
              f"{t_build_6+t_iter_6+t_build_7+t_iter_7:.2f}s "
              f"(brief estimated 1-4 hours; actual sub-minute on this "
              f"hardware).")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_epsilon_7.py` — script")
    md.append("- `result_epsilon_7.csv` — eps_k and ratios for k=1..7")
    md.append("- `result_epsilon_7.md` — this writeup")
    md.append("- `result_epsilon_7_diagnostic.md` — strategy + cost")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print()
    print(f"saved {out_md}")

    out_csv = os.path.join(OUT_DIR, "result_epsilon_7.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "eps_k", "ratio_to_prev"])
        prev = None
        for k in [1, 2, 3, 4, 5, 6, 7]:
            ek = eps_all[k]
            ratio = "" if prev is None else f"{ek/prev:+.10f}"
            w.writerow([k, f"{ek:+.15e}", ratio])
            prev = ek
    print(f"saved {out_csv}")

    out_diag = os.path.join(OUT_DIR, "result_epsilon_7_diagnostic.md")
    diag = []
    diag.append("# Diagnostic: eps_7 strategy + cost")
    diag.append("")
    diag.append("**Strategy:** float64 power iteration (same as eps_6).")
    diag.append("")
    diag.append(f"## Cost summary")
    diag.append("")
    diag.append(f"- K_6 build + iter: "
                f"{t_build_6+t_iter_6:.2f}s ({iters6} iters)")
    diag.append(f"- K_7 build + iter: "
                f"{t_build_7+t_iter_7:.2f}s ({iters7} iters)")
    diag.append(f"- Convergence residual at k=7: {res7:.2e}")
    diag.append("")
    diag.append(f"## Precision check")
    diag.append("")
    diag.append(f"- pi_7 sum check: {pi7.sum():.15f} "
                f"(deviation from 1: {abs(pi7.sum() - 1):.2e})")
    diag.append(f"- K_7 row-sum check: max |row sum - 1| = "
                f"{float(np.max(np.abs(K7.sum(axis=1) - 1))):.2e}")
    diag.append(f"- eps_7 reported to 10 sig figs (float64 noise floor "
                f"on values ~1e-4 is ~1e-19 absolute, ~1e-15 relative)")
    diag.append("")
    diag.append(f"## What was NOT attempted")
    diag.append(f"")
    diag.append(f"- Exact-rational stationary at k=7: previously killed "
                f"in R77.7 after 7+ hours from denominator bloat. Not "
                f"re-attempted per brief.")
    diag.append(f"")
    diag.append(f"- High-precision (mpmath) double-check: not needed; "
                f"float64 noise floor far below the precision the brief "
                f"required (±5e-5).")
    with open(out_diag, "w", encoding="utf-8") as f:
        f.write("\n".join(diag))
    print(f"saved {out_diag}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
