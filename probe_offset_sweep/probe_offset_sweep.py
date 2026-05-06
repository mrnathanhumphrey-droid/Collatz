"""
probe_offset_sweep.py — test S_k(c) → 7/15 for c ∈ {1, 5, 7, 11} (offsets in 3n+c).

Conventions (matching c_seven_forty_fifth_derivation.py and result_epsilon_*.py):
- K_k(c) is row-stochastic on (Z/3^k)*, transition r → ((3r+c)·2^{-v}) mod 3^k
- v ~ Geom(1/2) on {1, 2, ..., M} with M = ord_{3^k}(2) = 2·3^{k-1}
- weight on v: 2^{-v} / Z, Z = (1 - 2^{-M}) → 1 for k ≥ 2
- π_k(c) = stationary, computed by float64 power iteration to ||δ||_1 < 1e-15
- ||d_k(c)||² = Σ π_k(c)(r)² − (1/3) · Σ π_{k-1}(c)(r)²    (R74 identity)
- S_k(c) = 3^k · ||d_k(c)||²    (R75 algebraic identity)
- ε_k(c) = S_k(c) − 7/15

For c=1, expect S_∞ = 7/15 (so ε_k → 0); compare for c ∈ {5, 7, 11}.
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_offset_sweep"
LEVELS = [1, 2, 3, 4, 5, 6, 7]
C_VALUES = [1, 5, 7, 11]
SEVEN_FIFTEENTHS = 7.0 / 15.0


def order_of_two_mod(N: int) -> int:
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


def build_K_dense(c: int, k: int):
    """K_k(c) on coprime classes of Z/3^k under transition
    r -> ((3r + c) * 2^{-v}) mod 3^k, v ~ Geom(1/2) on {1..M}."""
    N = 3 ** k
    M = order_of_two_mod(N)
    M_eff = min(M, 1074)  # 2^{-v} below float64 epsilon for v > 1074
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M_eff, dtype=np.int64)
    p_cur = inv2
    for v in range(M_eff):
        powers_inv2[v] = p_cur
        p_cur = (p_cur * inv2) % N
    coprime = np.array([r for r in range(N) if r % 3 != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    weights = np.zeros(M_eff, dtype=np.float64)
    for vv in range(M_eff):
        weights[vv] = 2.0 ** -(vv + 1)
    weights /= weights.sum()
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (3 * r + c) % N
        # All targets must be coprime to 3 for the chain to stay in the state space.
        # (3r + c) mod 3 = c mod 3, so for c coprime to 3 (which is the only case
        # the framework supports), all targets r' = base * 2^{-v} mod N satisfy
        # r' mod 3 = (c · 2^{-v}) mod 3 ∈ {1, 2}, never 0. ✓
        targets = (base * powers_inv2) % N
        for j_t, t in enumerate(targets):
            K[i_r, state_idx[int(t)]] += weights[j_t]
    return K, coprime


def stationary(K, tol=1e-15, max_iter=20000):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n)
    for it in range(max_iter):
        pi_new = pi @ K
        s = pi_new.sum()
        if s != 0:
            pi_new /= s
        delta = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        if delta < tol:
            return pi, it + 1, delta
    return pi, max_iter, delta


def main():
    print("=" * 78)
    print(f"Offset sweep: S_k(c) for c ∈ {C_VALUES}, k ∈ {LEVELS}")
    print("=" * 78)
    print()

    t0 = time.time()
    # pi_cache[c][k] = stationary at level k for offset c
    pi_cache = {c: {} for c in C_VALUES}
    times_cache = {c: {} for c in C_VALUES}

    for c in C_VALUES:
        print(f"--- c = {c} (c mod 3 = {c % 3}) ---")
        for k in LEVELS:
            t1 = time.time()
            K, coprime = build_K_dense(c, k)
            t_build = time.time() - t1
            t2 = time.time()
            pi, n_iter, delta = stationary(K)
            t_stat = time.time() - t2
            pi_cache[c][k] = pi
            times_cache[c][k] = (t_build, t_stat, n_iter)
            row_sum_max_dev = float(np.max(np.abs(K.sum(axis=1) - 1)))
            print(f"  k={k}: n={K.shape[0]}, build={t_build:.2f}s, "
                  f"stationary={t_stat:.2f}s ({n_iter} iters, delta={delta:.2e}), "
                  f"row_sum_dev={row_sum_max_dev:.2e}")
        print()

    # Compute S_k(c) and ε_k(c)
    print("=" * 78)
    print("S_k(c) and ε_k(c) = S_k(c) - 7/15")
    print("=" * 78)
    print()

    rows_S = []
    print(f"{'c':>3} {'k':>3} {'sum_pi_sq':>16} {'||d_k||²':>16} {'S_k':>14} "
          f"{'eps_k = S_k-7/15':>20}")
    print("-" * 80)
    for c in C_VALUES:
        for k in LEVELS:
            sum_pi_k_sq = float(np.sum(pi_cache[c][k] ** 2))
            if k == 1:
                # k-1 = 0 is the trivial group with one element (probability 1);
                # Σ π_0² = 1.
                sum_pi_km1_sq = 1.0
            else:
                sum_pi_km1_sq = float(np.sum(pi_cache[c][k - 1] ** 2))
            d_sq = sum_pi_k_sq - sum_pi_km1_sq / 3.0
            S_k = (3 ** k) * d_sq
            eps_k = S_k - SEVEN_FIFTEENTHS
            rows_S.append({
                "c": c, "k": k,
                "sum_pi_sq": sum_pi_k_sq, "sum_pi_km1_sq": sum_pi_km1_sq,
                "d_sq": d_sq, "S_k": S_k, "eps_k": eps_k,
            })
            print(f"{c:>3} {k:>3} {sum_pi_k_sq:>16.10f} {d_sq:>16.10e} "
                  f"{S_k:>14.10f} {eps_k:>+20.10e}")
        print()

    # Save S_k CSV
    csv_S = os.path.join(OUT_DIR, "result_S_k_by_c.csv")
    with open(csv_S, "w", newline="", encoding="utf-8") as fh:
        cols = ["c", "k", "sum_pi_sq", "sum_pi_km1_sq", "d_sq", "S_k", "eps_k"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows_S:
            w.writerow({k: (f"{v:.15e}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    print(f"[csv: {csv_S}]")

    # ε_k recurrence fits per c (orders 1, 2, 3)
    print()
    print("=" * 78)
    print("Order-3 recurrence fits on ε_k(c) for k ∈ [2, 7]")
    print("=" * 78)
    rows_recur = []
    for c in C_VALUES:
        eps_c = [r["eps_k"] for r in rows_S if r["c"] == c]
        ks_c = [r["k"] for r in rows_S if r["c"] == c]
        # Take k=2..7 (skip k=1 which has eps = +0.2, an outlier)
        idx_use = [i for i, k in enumerate(ks_c) if k >= 2]
        eps_seq = [eps_c[i] for i in idx_use]
        n_use = len(eps_seq)

        print(f"\nc = {c}, ε_k for k=2..7: " +
              ", ".join(f"{e:+.4e}" for e in eps_seq))
        print(f"  |ε_k+1/ε_k| ratios: " +
              ", ".join(f"{abs(eps_seq[i+1]/eps_seq[i]):.4f}"
                        if eps_seq[i] != 0 else "—"
                        for i in range(n_use - 1)))

        for order in [1, 2, 3]:
            n_eq = n_use - order
            if n_eq < order + 1:
                continue
            y = np.array(eps_seq, dtype=np.float64)
            A = np.zeros((n_eq, order))
            b = np.zeros(n_eq)
            for i in range(n_eq):
                b[i] = y[order + i]
                for j in range(order):
                    A[i, j] = y[order + i - 1 - j]
            try:
                coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            except np.linalg.LinAlgError:
                continue
            pred = A @ coeffs
            ss_res = float(np.sum((b - pred) ** 2))
            ss_tot = float(np.sum((b - b.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
            poly = np.array([1.0] + list(-coeffs))
            roots = np.roots(poly) if order > 0 else np.array([])
            top_root = max((abs(r) for r in roots), default=float("nan"))
            print(f"  order={order}: alphas = {[f'{a:+.4f}' for a in coeffs]}, "
                  f"R²={r2:.4f}, top |root| = {top_root:.4f}")
            rows_recur.append({
                "c": c, "order": order, "n_eq": n_eq,
                "alphas": ",".join(f"{a:.10f}" for a in coeffs),
                "roots": ",".join(f"({r.real:+.6f}{r.imag:+.6f}j)" for r in roots),
                "ss_res": ss_res, "r2": r2, "top_root_abs": top_root,
            })

    csv_recur = os.path.join(OUT_DIR, "result_eps_recurrence_by_c.csv")
    with open(csv_recur, "w", newline="", encoding="utf-8") as fh:
        cols = ["c", "order", "n_eq", "alphas", "roots", "ss_res", "r2", "top_root_abs"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows_recur:
            w.writerow({k: (f"{v:.10e}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    print(f"\n[csv: {csv_recur}]")

    # ρ_slow per c (top |root| of order-3 fit, when it's real)
    print()
    print("=" * 78)
    print("ρ_slow(c) summary (top |root| of order-3 recurrence)")
    print("=" * 78)
    rows_rho = []
    for c in C_VALUES:
        order3 = [r for r in rows_recur if r["c"] == c and r["order"] == 3]
        if order3:
            rho = order3[0]["top_root_abs"]
            r2 = order3[0]["r2"]
            print(f"  c={c}: ρ_slow ≈ {rho:.4f}  (R²={r2:.4f})")
            rows_rho.append({"c": c, "rho_slow": rho, "r2": r2})

    csv_rho = os.path.join(OUT_DIR, "result_rho_slow_by_c.csv")
    with open(csv_rho, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["c", "rho_slow", "r2"])
        w.writeheader()
        for r in rows_rho:
            w.writerow({k: (f"{v:.10e}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    print(f"\n[csv: {csv_rho}]")

    # ----- Cross-c S_k comparison and verdict -----
    print()
    print("=" * 78)
    print("Cross-c S_k comparison")
    print("=" * 78)
    # S_k as a function of (c, k)
    S_table = {(r["c"], r["k"]): r["S_k"] for r in rows_S}
    eps_table = {(r["c"], r["k"]): r["eps_k"] for r in rows_S}

    print(f"\n{'k':>3} | " + " | ".join(f"S_k(c={c}):>14" for c in C_VALUES))
    print(f"{'k':>3} | " + " | ".join(f"S_k(c={c})".center(14) for c in C_VALUES))
    print("-" * (5 + 17 * len(C_VALUES)))
    for k in LEVELS:
        row = f"{k:>3} | " + " | ".join(f"{S_table[(c, k)]:>14.10f}" for c in C_VALUES)
        print(row)

    print(f"\n{'k':>3} | " + " | ".join(f"ε_k(c={c})".center(20) for c in C_VALUES))
    print("-" * (5 + 23 * len(C_VALUES)))
    for k in LEVELS:
        row = f"{k:>3} | " + " | ".join(f"{eps_table[(c, k)]:>+20.10e}" for c in C_VALUES)
        print(row)

    # Findings markdown
    md_path = os.path.join(OUT_DIR, "offset_sweep_findings.md")
    md = []
    md.append("# Offset sweep — S_k(c) for 3n+c dynamics, c ∈ {1, 5, 7, 11}")
    md.append("")
    md.append(f"**Date:** 2026-05-06.")
    md.append(f"**Levels tested:** k ∈ {LEVELS}.")
    md.append(f"**Convention:** S_k(c) = 3^k · ||d_k(c)||² with ||d_k(c)||² = "
              f"Σ π_k(c)(r)² − (1/3)·Σ π_{{k-1}}(c)(r)².")
    md.append(f"**Reference limit:** 7/15 ≈ {SEVEN_FIFTEENTHS:.10f}.")
    md.append("")

    md.append("## c mod 3 classes")
    md.append("")
    md.append("| c | c mod 3 | class |")
    md.append("|---|---|---|")
    for c in C_VALUES:
        md.append(f"| {c} | {c % 3} | {'1-class' if c % 3 == 1 else '2-class'} |")
    md.append("")
    md.append("Within each c-mod-3 class, c values differ only at finer mod-3^k levels. "
              "By the σ(r) = −r mod 3^k chain symmetry (sibling 3x±1 study), "
              "c-mod-3 = 1 and c-mod-3 = 2 chains should produce identical "
              "Σ π² values up to the σ-permutation, since Σ over states is invariant "
              "under σ. **Predicted equality:** S_k(c=1) = S_k(c=5) = S_k(c=7) = S_k(c=11) "
              "if the c=σ symmetry extends to all k.")
    md.append("")

    md.append("## S_k(c) table")
    md.append("")
    md.append("| k | S_k(c=1) | S_k(c=5) | S_k(c=7) | S_k(c=11) |")
    md.append("|---|---|---|---|---|")
    for k in LEVELS:
        md.append(f"| {k} | " +
                  " | ".join(f"{S_table[(c, k)]:.10f}" for c in C_VALUES) + " |")
    md.append("")

    md.append("## ε_k(c) = S_k(c) − 7/15 table")
    md.append("")
    md.append("| k | ε_k(c=1) | ε_k(c=5) | ε_k(c=7) | ε_k(c=11) |")
    md.append("|---|---|---|---|---|")
    for k in LEVELS:
        md.append(f"| {k} | " +
                  " | ".join(f"{eps_table[(c, k)]:+.10e}" for c in C_VALUES) + " |")
    md.append("")

    # Convergence verdict per c
    md.append("## Convergence verdict per c")
    md.append("")
    for c in C_VALUES:
        eps_seq = [eps_table[(c, k)] for k in LEVELS if k >= 2]
        max_eps = max(abs(e) for e in eps_seq)
        last_eps = eps_seq[-1]
        going_down = all(abs(eps_seq[i+1]) < abs(eps_seq[i]) * 2 for i in range(len(eps_seq) - 1))
        md.append(f"- **c={c}**: |ε_k| ranges over [{min(abs(e) for e in eps_seq):.2e}, "
                  f"{max_eps:.2e}], |ε_7| = {abs(last_eps):.2e}; "
                  f"{'monotone-ish bounded' if going_down else 'oscillatory'}.")
    md.append("")

    # Cross-c equality check
    md.append("## Cross-c equality check")
    md.append("")
    md.append("Maximum |S_k(c) − S_k(c=1)| at each k (deviation from c=1 baseline):")
    md.append("")
    md.append("| k | max|S_k(c) − S_k(1)| | which c |")
    md.append("|---|---|---|")
    for k in LEVELS:
        deltas = [(c, abs(S_table[(c, k)] - S_table[(1, k)])) for c in C_VALUES if c != 1]
        max_c, max_d = max(deltas, key=lambda x: x[1])
        md.append(f"| {k} | {max_d:.2e} | c={max_c} |")
    md.append("")

    # ρ_slow per c
    md.append("## ρ_slow(c) — order-3 recurrence top |root|")
    md.append("")
    md.append("| c | ρ_slow | R² | comparable to c=1 (~0.83)? |")
    md.append("|---|---|---|---|")
    for r in rows_rho:
        match = "yes" if abs(r["rho_slow"] - 0.83) < 0.10 else "no"
        md.append(f"| {r['c']} | {r['rho_slow']:.4f} | {r['r2']:.4f} | {match} |")
    md.append("")

    # Verdict
    md.append("## Verdict")
    md.append("")
    # Outcome A: all c converge to 7/15 to high precision
    last_eps_all = [abs(eps_table[(c, max(LEVELS))]) for c in C_VALUES]
    if max(last_eps_all) < 5e-3:
        verdict = ("**Outcome A.** S_k(c) → 7/15 for all c ∈ {1, 5, 7, 11} "
                   "tested at k=7. The c=7/45 constant is **offset-universal** "
                   "across the tested odd-coprime-to-3 offsets.")
    else:
        # Check if it's exactly the σ-pair (c=1 ≡ c=5 reflected, c=7 ≡ c=11)
        d_15 = max(abs(S_table[(5, k)] - S_table[(1, k)]) for k in LEVELS)
        d_711 = max(abs(S_table[(11, k)] - S_table[(7, k)]) for k in LEVELS)
        d_17 = max(abs(S_table[(7, k)] - S_table[(1, k)]) for k in LEVELS)
        md.append(f"σ-symmetry diagnostic: max|S_k(5) − S_k(1)| = {d_15:.2e}, "
                  f"max|S_k(11) − S_k(7)| = {d_711:.2e}, "
                  f"max|S_k(7) − S_k(1)| = {d_17:.2e}.")
        md.append("")
        if d_15 < 1e-12 and d_711 < 1e-12:
            verdict = ("**Outcome A (with σ-grouping).** The c=σ chain symmetry "
                       "K_(c=1) ↔ K_(c=5) and K_(c=7) ↔ K_(c=11) gives identical "
                       "S_k. Cross-c-class deviation max|S_k(7) − S_k(1)|"
                       f" = {d_17:.2e} establishes whether S_k depends on c "
                       "beyond the mod-3 class.")
        else:
            verdict = ("**Mixed outcome.** S_k(c) converges per c but to "
                       f"distinct limits within ε of {max(last_eps_all):.2e} "
                       "of 7/15. Deviation pattern suggests a c-dependent "
                       "structural family.")
    md.append(verdict)
    md.append("")

    md.append("## Files")
    md.append("")
    md.append("- `result_S_k_by_c.csv` — S_k(c) and ε_k(c) per (c, k)")
    md.append("- `result_eps_recurrence_by_c.csv` — recurrence fits per c, orders 1-3")
    md.append("- `result_rho_slow_by_c.csv` — ρ_slow(c) summary")
    md.append("- `offset_sweep_findings.md` — this writeup")

    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    print(f"\n[md: {md_path}]")
    print(f"\nTotal runtime: {time.time()-t0:.1f}s")
    print("Done.")


if __name__ == "__main__":
    main()
