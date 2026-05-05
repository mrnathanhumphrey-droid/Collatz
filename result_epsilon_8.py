"""
result_epsilon_8.py
===================
Compute eps_8 = S_8 - 7/15 to test the bouncing/oscillation hypothesis.

Background: at k=7 we got eps_7 = -1.18e-3, with |eps_7/eps_6| = 2.36 — a
clear reversal of the monotone-decreasing ratio trajectory observed for
k=2..6. This breaks the simple "rate-1/2 envelope decay" reading and opens
three structural possibilities:
  H1: damped oscillation (single rotating mode)
  H2: multi-mode mode-crossing
  H3: continued non-monotone trajectory with no clean asymptotic form yet

Strategy: float64 power iteration on K_8 (4374 states) — same path as eps_7,
plus scipy.sparse.linalg.eigs cross-check on the leading eigenpair. No
exact-rational attempt (denominator bloat killed it at k=7 after 7+ hours).

Decision rule on |eps_8/eps_7|:
  > 0.7  bouncing/oscillation real
  < 0.5  decay resumed
  0.5..0.7  ambiguous, eps_9 needed
"""
from __future__ import annotations

import csv
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
ENV_CSV = os.path.join(OUT_DIR, "result_q_sweep_test_1_envelope.csv")


def order_of_two(N):
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
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


def stationary_power_iteration(K, max_iter=10000, tol=1e-13):
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


def stationary_eigs(K):
    """Cross-check: leading eigenvector of K^T via sparse Arnoldi."""
    n = K.shape[0]

    def mv(v):
        return K.T @ v

    op = LinearOperator((n, n), matvec=mv, dtype=np.float64)
    vals, vecs = eigs(op, k=1, which="LM", maxiter=10000, tol=1e-12)
    val = float(vals[0].real)
    vec = vecs[:, 0].real
    # Normalize so vec sums to 1 (probability)
    s = vec.sum()
    if s < 0:
        vec = -vec
        s = -s
    vec = vec / s
    return val, vec


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
    print("Computing eps_8 = S_8 - 7/15 (bouncing vs decay-resumed)")
    print("=" * 78)

    S_cache = load_cached_S()
    X5_exact = Fraction(1)
    for k in [1, 2, 3, 4, 5]:
        X5_exact += S_cache[k]
    X5_float = float(X5_exact)
    print(f"\nX_5 (exact, cached) = {X5_float:.15f}")

    eps_all = {k: float(S_cache[k] - Fraction(7, 15)) for k in S_cache}

    # k=6 (sanity)
    print("\n--- k = 6 (sanity recompute) ---")
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
    eps_all[6] = eps6
    print(f"  X_6 = {X6:.15f}, S_6 = {S6:.15f}")
    print(f"  eps_6 = {eps6:+.10e}")
    del K6, pi6  # free memory

    # k=7 (sanity)
    print("\n--- k = 7 (sanity recompute) ---")
    t0 = time.time()
    K7, _, M7, _ = build_K_float(3, 7)
    t_build_7 = time.time() - t0
    print(f"  build K_7 (states={K7.shape[0]}, M={M7}): {t_build_7:.2f}s")
    t0 = time.time()
    pi7, iters7, res7 = stationary_power_iteration(K7, tol=1e-14)
    t_iter_7 = time.time() - t0
    print(f"  power iteration: {iters7} iters, residual {res7:.2e}, "
          f"{t_iter_7:.2f}s")
    sum_sq_7 = float((pi7 ** 2).sum())
    X7 = (3 ** 7) * sum_sq_7
    S7 = X7 - X6
    eps7 = S7 - 7.0 / 15.0
    eps_all[7] = eps7
    print(f"  X_7 = {X7:.15f}, S_7 = {S7:.15f}")
    print(f"  eps_7 = {eps7:+.10e}")
    del K7, pi7

    # k=8 (target)
    print("\n--- k = 8 (target) ---")
    t0 = time.time()
    K8, _, M8, _ = build_K_float(3, 8)
    t_build_8 = time.time() - t0
    n8 = K8.shape[0]
    print(f"  build K_8 (states={n8}, M={M8}): {t_build_8:.2f}s")
    print(f"  K_8 row-sum check: max |row sum - 1| = "
          f"{float(np.max(np.abs(K8.sum(axis=1) - 1))):.2e}")

    t0 = time.time()
    pi8, iters8, res8 = stationary_power_iteration(K8, tol=1e-13,
                                                    max_iter=10000)
    t_iter_8 = time.time() - t0
    print(f"  power iteration: {iters8} iters, residual {res8:.2e}, "
          f"{t_iter_8:.2f}s")
    print(f"  sum(pi_8) = {pi8.sum():.15f}  (expect 1)")

    # Cross-check with scipy.eigs
    print("\n  cross-check via scipy.sparse.linalg.eigs (Arnoldi)...")
    t0 = time.time()
    val_eigs, pi8_eigs = stationary_eigs(K8)
    t_eigs = time.time() - t0
    print(f"  eigs: leading eigenvalue = {val_eigs:.12f}  "
          f"(expect 1.0), {t_eigs:.2f}s")
    diff_l1 = float(np.linalg.norm(pi8 - pi8_eigs, ord=1))
    diff_max = float(np.max(np.abs(pi8 - pi8_eigs)))
    print(f"  |pi8_power - pi8_eigs|_1 = {diff_l1:.4e}, "
          f"max-abs = {diff_max:.4e}")

    # Compute S_8, eps_8
    sum_sq_8 = float((pi8 ** 2).sum())
    X8 = (3 ** 8) * sum_sq_8
    S8 = X8 - X7
    eps8 = S8 - 7.0 / 15.0
    eps_all[8] = eps8
    print(f"\n  X_8 = {X8:.15f}, S_8 = {S8:.15f}")
    print(f"  eps_8 = {eps8:+.10e}")

    # Cross-check S_8 via eigs vector
    sum_sq_8_eigs = float((pi8_eigs ** 2).sum())
    X8_eigs = (3 ** 8) * sum_sq_8_eigs
    S8_eigs = X8_eigs - X7
    eps8_eigs = S8_eigs - 7.0 / 15.0
    print(f"  cross-check eps_8 via eigs vector: {eps8_eigs:+.10e}")
    print(f"  agreement: |eps8_power - eps8_eigs| = "
          f"{abs(eps8 - eps8_eigs):.4e}")

    # FFT-based S_8 (full Plancherel sanity)
    # pi_8 lives on (Z/3^8)* (4374 states); embed to Z/3^8 (6561 entries) and FFT
    print("\n  FFT cross-check of S_8 = sum_{xi != 0 mod 3} |pi_hat|^2...")
    N8 = 3 ** 8
    pi_full = np.zeros(N8, dtype=np.float64)
    coprime_8 = np.array([r for r in range(N8) if r % 3 != 0], dtype=np.int64)
    pi_full[coprime_8] = pi8
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N8)
    mask_nontrivial = xi_arr % 3 != 0
    S8_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps8_fft = S8_fft - 7.0 / 15.0
    print(f"  S_8 via FFT = {S8_fft:.15f}, eps_8_fft = {eps8_fft:+.10e}")
    print(f"  agreement: |S_8_X - S_8_FFT| = {abs(S8 - S8_fft):.4e}")

    # Ratio
    ratio_78 = abs(eps8 / eps7)
    print(f"\n|eps_8 / eps_7| = {ratio_78:.6f}")

    # Decision
    if ratio_78 > 0.7:
        verdict = ("BOUNCING/OSCILLATION REAL — |eps_8/eps_7| > 0.7. The "
                   "k=7 reversal is not a one-off; oscillatory structure "
                   "in k-space is genuine. Distinguishes H1 (damped "
                   "oscillation) and H2 (multi-mode crossing) from H3 "
                   "(no clean form). Period and amplitude characterization "
                   "needs eps_9, eps_10.")
    elif ratio_78 < 0.5:
        verdict = ("DECAY RESUMED — |eps_8/eps_7| < 0.5. The k=7 spike was "
                   "a transient; underlying decay continues. The rate-1/2 "
                   "envelope reading survives, modulo the k=7 transient. "
                   "Asymptotic ratio still tracking toward 1/2 or below.")
    else:
        verdict = ("AMBIGUOUS — 0.5 < |eps_8/eps_7| < 0.7. Neither clean "
                   "decay nor sustained bouncing. eps_9 (k=9 chain, ~13k "
                   "states, ~3-5 min compute) would disambiguate.")

    print(f"\n*** VERDICT: {verdict}")

    # Full ratio trajectory
    print("\nUpdated ratio trajectory |eps_{k+1}/eps_k|:")
    for k in [1, 2, 3, 4, 5, 6, 7]:
        rr = abs(eps_all[k + 1] / eps_all[k])
        print(f"  k = {k} -> {k+1}: {rr:.6f}")

    # Outputs
    out_md = os.path.join(OUT_DIR, "result_epsilon_8.md")
    md = []
    md.append("# Result: eps_8 = S_8 - 7/15")
    md.append("")
    md.append(f"**Date:** 2026-05-05.  Float64 power iteration on K_8 "
              f"({n8} states, M={M8}), with scipy.sparse.linalg.eigs "
              f"Arnoldi cross-check.")
    md.append("")
    md.append(f"## Verdict")
    md.append("")
    md.append(verdict)
    md.append("")
    md.append("## Headline numbers")
    md.append("")
    md.append(f"- eps_8 = `{eps8:+.10e}` (power iteration)")
    md.append(f"- eps_8 = `{eps8_eigs:+.10e}` (scipy.eigs cross-check; "
              f"agreement {abs(eps8 - eps8_eigs):.2e})")
    md.append(f"- eps_8 = `{eps8_fft:+.10e}` (FFT cross-check; agreement "
              f"{abs(eps8 - eps8_fft):.2e})")
    md.append(f"- |eps_8/eps_7| = **{ratio_78:.6f}**")
    md.append(f"- |eps_7/eps_6| = {abs(eps7/eps6):.6f} (prior, for context)")
    md.append("")
    md.append("## Ratio trajectory (k=1..8)")
    md.append("")
    md.append("| k → k+1 | |eps_{k+1}/eps_k| |")
    md.append("|---|---|")
    for k in [1, 2, 3, 4, 5, 6, 7]:
        rr = abs(eps_all[k + 1] / eps_all[k])
        md.append(f"| {k} → {k+1} | {rr:.6f} |")
    md.append("")
    md.append("## eps_k table (k=1..8)")
    md.append("")
    md.append("| k | eps_k | source |")
    md.append("|---|---|---|")
    for k in [1, 2, 3, 4, 5]:
        md.append(f"| {k} | {float(S_cache[k] - Fraction(7,15)):+.10e} | "
                  f"exact rational (cached) |")
    md.append(f"| 6 | {eps6:+.10e} | float64 power iter, K_6 ({3**6} states) |")
    md.append(f"| 7 | {eps7:+.10e} | float64 power iter, K_7 (1458 states) |")
    md.append(f"| 8 | {eps8:+.10e} | float64 power iter, K_8 ({n8} states) "
              f"+ eigs cross-check |")
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
    md.append(f"| K_8 build | {t_build_8:.2f}s | — | — |")
    md.append(f"| K_8 power iter | {t_iter_8:.2f}s | {iters8} | "
              f"{res8:.2e} |")
    md.append(f"| K_8 eigs cross-check | {t_eigs:.2f}s | — | — |")
    md.append("")
    md.append(f"**Cross-checks at k=8:**")
    md.append(f"- power iter pi vs eigs vector: |·|_1 diff = {diff_l1:.2e}, "
              f"max-abs = {diff_max:.2e}")
    md.append(f"- S_8 via X_8 - X_7 vs S_8 via FFT: diff = "
              f"{abs(S8 - S8_fft):.2e}")
    md.append(f"- pi_8 sum = {pi8.sum():.15f}")
    md.append(f"- K_8 row-sum max deviation from 1 = "
              f"{float(np.max(np.abs(K8.sum(axis=1) - 1))):.2e}")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_epsilon_8.py` — script")
    md.append("- `result_epsilon_8.csv` — eps_k and ratios for k=1..8")
    md.append("- `result_epsilon_8.md` — this writeup")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"\nsaved {out_md}")

    out_csv = os.path.join(OUT_DIR, "result_epsilon_8.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "eps_k", "abs_eps_k", "ratio_to_prev_signed",
                    "abs_ratio_to_prev"])
        prev = None
        for k in [1, 2, 3, 4, 5, 6, 7, 8]:
            ek = eps_all[k]
            ratio_signed = "" if prev is None else f"{ek/prev:+.10f}"
            ratio_abs = "" if prev is None else f"{abs(ek/prev):.10f}"
            w.writerow([k, f"{ek:+.15e}", f"{abs(ek):.15e}",
                        ratio_signed, ratio_abs])
            prev = ek
    print(f"saved {out_csv}")
    print("\nDone.")


if __name__ == "__main__":
    main()
