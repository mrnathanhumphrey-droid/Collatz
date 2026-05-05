"""
result_epsilon_9.py
===================
Compute eps_9 = S_9 - 7/15 to disambiguate the 0.5..0.7 ratio band at k=8.

At k=8 we got |eps_8/eps_7| = 0.634 — between "decay resumed" (<0.5,
matching pre-spike asymptote of 0.43-0.55) and "sustained bouncing" (>0.7).
eps_9 distinguishes:
  - Returns to ~0.43 → confirms decay resumed; k=7 spike was transient.
  - Stays in 0.55-0.65 band → multi-mode/slow oscillation, asymptotic
    rate is higher than pre-spike trajectory implied.
  - Spikes again → period-2 oscillation in k.

K_9 has 2*3^8 = 13122 states. Float64 dense matrix is 1.38 GB. Need
vectorized bincount build (Python inner loop at this size = minutes).
"""
from __future__ import annotations

import csv
import gc
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


def build_K_float_vectorized(q, k):
    """Build K via vectorized bincount per row. Each row = M weights placed
    on unique target columns (2^(-v) ranges over (Z/3^k)* exactly once)."""
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
    Z_v = 1.0 - 2.0 ** (-M)
    weights = (2.0 ** (-np.arange(1, M + 1))) / Z_v
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (q * r + 1) % N
        targets = (base * powers_inv2) % N
        js = state_idx[targets]
        K[i_r] = np.bincount(js, weights=weights, minlength=n)
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
    n = K.shape[0]

    def mv(v):
        return K.T @ v

    op = LinearOperator((n, n), matvec=mv, dtype=np.float64)
    vals, vecs = eigs(op, k=1, which="LM", maxiter=10000, tol=1e-12)
    val = float(vals[0].real)
    vec = vecs[:, 0].real
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


def compute_X_S_eps(K, X_prev, k):
    pi, iters, res = stationary_power_iteration(K, tol=1e-14)
    sum_sq = float((pi ** 2).sum())
    X = (3 ** k) * sum_sq
    S = X - X_prev
    eps = S - 7.0 / 15.0
    return pi, iters, res, X, S, eps


def main():
    print("=" * 78)
    print("Computing eps_9 = S_9 - 7/15 (disambiguate 0.5..0.7 band at k=8)")
    print("=" * 78)

    S_cache = load_cached_S()
    X = {0: 1.0}
    for k in [1, 2, 3, 4, 5]:
        X[k] = X[k - 1] + float(S_cache[k])
    print(f"\nX_5 (exact, cached) = {X[5]:.15f}")

    eps_all = {k: float(S_cache[k] - Fraction(7, 15)) for k in S_cache}

    # k = 6, 7, 8 — sanity recomputes
    for k in [6, 7, 8]:
        print(f"\n--- k = {k} (sanity recompute) ---")
        t0 = time.time()
        K, _, M, _ = build_K_float_vectorized(3, k)
        t_build = time.time() - t0
        print(f"  build K_{k} (states={K.shape[0]}, M={M}): "
              f"{t_build:.2f}s")
        t0 = time.time()
        pi, iters, res, X[k], S_k, eps_k = compute_X_S_eps(K, X[k - 1], k)
        eps_all[k] = eps_k
        t_iter = time.time() - t0
        print(f"  power iteration: {iters} iters, residual {res:.2e}, "
              f"{t_iter:.2f}s")
        print(f"  X_{k} = {X[k]:.15f}, S_{k} = {S_k:.15f}")
        print(f"  eps_{k} = {eps_k:+.10e}")
        del K, pi
        gc.collect()

    # k = 9 (target)
    print(f"\n--- k = 9 (TARGET) ---")
    print(f"  expected size: 13,122 states, dense float64 ~ 1.38 GB")
    t0 = time.time()
    K9, _, M9, _ = build_K_float_vectorized(3, 9)
    t_build_9 = time.time() - t0
    n9 = K9.shape[0]
    print(f"  build K_9 (states={n9}, M={M9}): {t_build_9:.2f}s")
    print(f"  K_9 row-sum check: max |row sum - 1| = "
          f"{float(np.max(np.abs(K9.sum(axis=1) - 1))):.2e}")

    t0 = time.time()
    pi9, iters9, res9 = stationary_power_iteration(K9, tol=1e-13,
                                                    max_iter=10000)
    t_iter_9 = time.time() - t0
    print(f"  power iteration: {iters9} iters, residual {res9:.2e}, "
          f"{t_iter_9:.2f}s")
    print(f"  sum(pi_9) = {pi9.sum():.15f}  (expect 1)")

    print("\n  cross-check via scipy.sparse.linalg.eigs (Arnoldi)...")
    t0 = time.time()
    val_eigs, pi9_eigs = stationary_eigs(K9)
    t_eigs = time.time() - t0
    print(f"  eigs: leading eigenvalue = {val_eigs:.12f}, "
          f"{t_eigs:.2f}s")
    diff_l1 = float(np.linalg.norm(pi9 - pi9_eigs, ord=1))
    diff_max = float(np.max(np.abs(pi9 - pi9_eigs)))
    print(f"  |pi9_power - pi9_eigs|_1 = {diff_l1:.4e}, "
          f"max-abs = {diff_max:.4e}")

    sum_sq_9 = float((pi9 ** 2).sum())
    X[9] = (3 ** 9) * sum_sq_9
    S9 = X[9] - X[8]
    eps9 = S9 - 7.0 / 15.0
    eps_all[9] = eps9
    print(f"\n  X_9 = {X[9]:.15f}, S_9 = {S9:.15f}")
    print(f"  eps_9 = {eps9:+.10e}")

    sum_sq_9_eigs = float((pi9_eigs ** 2).sum())
    X9_eigs = (3 ** 9) * sum_sq_9_eigs
    S9_eigs = X9_eigs - X[8]
    eps9_eigs = S9_eigs - 7.0 / 15.0
    print(f"  cross-check eps_9 via eigs vector: {eps9_eigs:+.10e}")
    print(f"  agreement: |eps9_power - eps9_eigs| = "
          f"{abs(eps9 - eps9_eigs):.4e}")

    print("\n  FFT cross-check of S_9...")
    N9 = 3 ** 9
    pi_full = np.zeros(N9, dtype=np.float64)
    coprime_9 = np.array([r for r in range(N9) if r % 3 != 0],
                          dtype=np.int64)
    pi_full[coprime_9] = pi9
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N9)
    mask_nontrivial = xi_arr % 3 != 0
    S9_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps9_fft = S9_fft - 7.0 / 15.0
    print(f"  S_9 via FFT = {S9_fft:.15f}, eps_9_fft = "
          f"{eps9_fft:+.10e}")
    print(f"  agreement: |S_9_X - S_9_FFT| = {abs(S9 - S9_fft):.4e}")

    eps7 = eps_all[7]
    eps8 = eps_all[8]
    ratio_89 = abs(eps9 / eps8)
    ratio_78 = abs(eps8 / eps7)
    print(f"\n|eps_9 / eps_8| = {ratio_89:.6f}  (compare |eps_8/eps_7| "
          f"= {ratio_78:.6f})")

    # Decision
    if 0.40 <= ratio_89 <= 0.55:
        verdict = ("DECAY RESUMED — |eps_9/eps_8| in 0.40..0.55, matching "
                   "pre-spike asymptote 0.43-0.55. The k=7 spike was a "
                   "transient; underlying decay continues. Rate-1/2 "
                   "envelope reading survives.")
    elif ratio_89 > 0.7:
        verdict = ("SUSTAINED BOUNCING — |eps_9/eps_8| > 0.7. Multi-step "
                   "structure in k-space is genuine; the 0.634 at k=8 was "
                   "not a return to decay. Period and amplitude need "
                   "k=10..12 to characterize.")
    elif ratio_89 < 0.40:
        verdict = ("ACCELERATED DECAY — |eps_9/eps_8| < 0.40, faster than "
                   "pre-spike. The k=7 spike was a one-off and post-spike "
                   "trajectory is faster than the k=2..6 trend. "
                   "Asymptotic rate may be > 1/2.")
    elif 0.55 < ratio_89 <= 0.7:
        verdict = ("ELEVATED RATIO — |eps_9/eps_8| in 0.55..0.70. Stays "
                   "above pre-spike band; multi-mode/slow oscillation is "
                   "consistent. The post-k=6 trajectory is structurally "
                   "different from the k=2..6 trajectory.")
    else:
        verdict = (f"UNEXPECTED — |eps_9/eps_8| = {ratio_89:.4f} doesn't "
                   "fit known decision bands. Inspect.")

    print(f"\n*** VERDICT: {verdict}")

    print("\nUpdated ratio trajectory |eps_{k+1}/eps_k|:")
    for k in [1, 2, 3, 4, 5, 6, 7, 8]:
        rr = abs(eps_all[k + 1] / eps_all[k])
        print(f"  k = {k} -> {k+1}: {rr:.6f}")

    # Outputs
    out_md = os.path.join(OUT_DIR, "result_epsilon_9.md")
    md = []
    md.append("# Result: eps_9 = S_9 - 7/15")
    md.append("")
    md.append(f"**Date:** 2026-05-05.  Float64 power iteration on K_9 "
              f"({n9} states, M={M9}), with scipy.sparse.linalg.eigs "
              f"Arnoldi cross-check and FFT verification.")
    md.append("")
    md.append(f"## Verdict")
    md.append("")
    md.append(verdict)
    md.append("")
    md.append("## Headline numbers")
    md.append("")
    md.append(f"- eps_9 = `{eps9:+.10e}` (power iteration)")
    md.append(f"- eps_9 = `{eps9_eigs:+.10e}` (eigs cross-check; "
              f"agreement {abs(eps9 - eps9_eigs):.2e})")
    md.append(f"- eps_9 = `{eps9_fft:+.10e}` (FFT cross-check; agreement "
              f"{abs(eps9 - eps9_fft):.2e})")
    md.append(f"- |eps_9/eps_8| = **{ratio_89:.6f}**")
    md.append(f"- |eps_8/eps_7| = {ratio_78:.6f} (prior, for context)")
    md.append("")
    md.append("## Ratio trajectory (k=1..9)")
    md.append("")
    md.append("| k → k+1 | |eps_{k+1}/eps_k| |")
    md.append("|---|---|")
    for k in [1, 2, 3, 4, 5, 6, 7, 8]:
        rr = abs(eps_all[k + 1] / eps_all[k])
        md.append(f"| {k} → {k+1} | {rr:.6f} |")
    md.append("")
    md.append("## eps_k table (k=1..9)")
    md.append("")
    md.append("| k | eps_k | source |")
    md.append("|---|---|---|")
    for k in [1, 2, 3, 4, 5]:
        md.append(f"| {k} | {float(S_cache[k] - Fraction(7,15)):+.10e} | "
                  f"exact rational (cached) |")
    md.append(f"| 6 | {eps_all[6]:+.10e} | float64 power iter (486 states) |")
    md.append(f"| 7 | {eps_all[7]:+.10e} | float64 power iter (1458 states) |")
    md.append(f"| 8 | {eps_all[8]:+.10e} | float64 power iter (4374 states) |")
    md.append(f"| 9 | {eps9:+.10e} | float64 power iter ({n9} states) "
              f"+ eigs + FFT cross-checks |")
    md.append("")
    md.append("## Computation diagnostics (k=9)")
    md.append("")
    md.append(f"- K_9 build time (vectorized bincount): {t_build_9:.2f}s")
    md.append(f"- Power iter: {iters9} iterations, residual = {res9:.2e}, "
              f"{t_iter_9:.2f}s")
    md.append(f"- eigs Arnoldi cross-check: {t_eigs:.2f}s, leading eval "
              f"= {val_eigs:.12f}")
    md.append(f"- |pi9_power - pi9_eigs|_1 = {diff_l1:.2e}, max-abs = "
              f"{diff_max:.2e}")
    md.append(f"- S_9 via X_9 - X_8 vs S_9 via FFT: diff = "
              f"{abs(S9 - S9_fft):.2e}")
    md.append(f"- pi_9 sum = {pi9.sum():.15f}")
    md.append(f"- K_9 row-sum max deviation from 1 = "
              f"{float(np.max(np.abs(K9.sum(axis=1) - 1))):.2e}")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_epsilon_9.py` — script")
    md.append("- `result_epsilon_9.csv` — eps_k and ratios for k=1..9")
    md.append("- `result_epsilon_9.md` — this writeup")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"\nsaved {out_md}")

    out_csv = os.path.join(OUT_DIR, "result_epsilon_9.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "eps_k", "abs_eps_k", "ratio_to_prev_signed",
                    "abs_ratio_to_prev"])
        prev = None
        for k in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
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
