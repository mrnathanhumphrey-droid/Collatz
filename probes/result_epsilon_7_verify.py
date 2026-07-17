"""
Verification: cross-check the float64 power-iteration path used in
result_epsilon_7.py against (a) exact rationals at k=5, (b) scipy
eigsolver at k=6 and k=7. If these agree, the surprising eps_7 result
is real and not a numerical artifact.
"""
from __future__ import annotations
import csv
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix

sys.stdout.reconfigure(encoding="utf-8")

ENV_CSV = r"C:\Collatz\result_q_sweep_test_1_envelope.csv"


def order_of_two(N):
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
    return K


def power_iter(K, max_iter=2000, tol=1e-14):
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


def eigs_left_dominant(K):
    """Find left dominant eigenvector via scipy.eigs on K^T."""
    KT = csr_matrix(K.T)
    vals, vecs = eigs(KT, k=1, which="LM", maxiter=10000, tol=1e-13)
    v = np.array(vecs[:, 0]).real
    if v.sum() < 0:
        v = -v
    v /= v.sum()
    return v, vals[0].real


def load_S():
    S = {}
    with open(ENV_CSV) as f:
        for row in csv.DictReader(f):
            if int(row["q"]) != 3:
                continue
            S[int(row["n"])] = Fraction(int(row["S_n_num"]),
                                          int(row["S_n_den"]))
    return S


def main():
    S_cache = load_S()
    print("=" * 78)
    print("Verification: float64 power iter vs exact rationals + scipy eigs")
    print("=" * 78)

    # Sanity at k=5
    print()
    print("--- k=5 sanity ---")
    print("Exact cached values:")
    X4_exact = Fraction(1) + sum(S_cache[k] for k in [1, 2, 3, 4])
    X5_exact = X4_exact + S_cache[5]
    print(f"  X_4 (exact) = {float(X4_exact):.15f}")
    print(f"  X_5 (exact) = {float(X5_exact):.15f}")
    print(f"  S_5 (exact) = {float(S_cache[5]):.15f}")

    print("\nMy float64 K_5 power iter:")
    t0 = time.time()
    K5 = build_K_float(3, 5)
    t_b5 = time.time() - t0
    print(f"  build K_5: {t_b5:.2f}s")
    pi5, iters5, res5 = power_iter(K5)
    sum_sq_5 = float((pi5 ** 2).sum())
    X5_float = (3 ** 5) * sum_sq_5
    S5_float = X5_float - float(X4_exact)
    print(f"  iters: {iters5}, residual: {res5:.2e}")
    print(f"  X_5 (float64) = {X5_float:.15f}")
    print(f"  S_5 (float64) = {S5_float:.15f}")
    print(f"  diff X_5: {X5_float - float(X5_exact):.3e}")
    print(f"  diff S_5: {S5_float - float(S_cache[5]):.3e}")
    if abs(X5_float - float(X5_exact)) < 1e-12:
        print("  ** k=5 float matches exact rational to 1e-12: float code is CORRECT")
    else:
        print("  ** k=5 float DOES NOT match exact: code has a bug, halt")
        return

    # k=6 sanity via two methods
    print()
    print("--- k=6 cross-check: power iter vs scipy.eigs ---")
    t0 = time.time()
    K6 = build_K_float(3, 6)
    print(f"  build K_6: {time.time()-t0:.2f}s")
    pi6_pwr, iters6, res6 = power_iter(K6)
    print(f"  power iter: {iters6} iters, res {res6:.2e}, "
          f"X_6 = {(3**6)*float((pi6_pwr**2).sum()):.15f}")
    pi6_eig, lam6 = eigs_left_dominant(K6)
    print(f"  scipy.eigs: lam={lam6:.10f}, "
          f"X_6 = {(3**6)*float((pi6_eig**2).sum()):.15f}")
    print(f"  L1 diff between methods: "
          f"{float(np.linalg.norm(pi6_pwr - pi6_eig, ord=1)):.3e}")

    # k=7 cross-check
    print()
    print("--- k=7 cross-check: power iter vs scipy.eigs ---")
    t0 = time.time()
    K7 = build_K_float(3, 7)
    print(f"  build K_7: {time.time()-t0:.2f}s")
    pi7_pwr, iters7, res7 = power_iter(K7)
    print(f"  power iter: {iters7} iters, res {res7:.2e}, "
          f"X_7 = {(3**7)*float((pi7_pwr**2).sum()):.15f}")
    pi7_eig, lam7 = eigs_left_dominant(K7)
    print(f"  scipy.eigs: lam={lam7:.10f}, "
          f"X_7 = {(3**7)*float((pi7_eig**2).sum()):.15f}")
    print(f"  L1 diff between methods: "
          f"{float(np.linalg.norm(pi7_pwr - pi7_eig, ord=1)):.3e}")

    # Now compute eps_5, eps_6, eps_7 from each method
    X5_p = (3**5)*float((pi5**2).sum())
    X6_p = (3**6)*float((pi6_pwr**2).sum())
    X6_e = (3**6)*float((pi6_eig**2).sum())
    X7_p = (3**7)*float((pi7_pwr**2).sum())
    X7_e = (3**7)*float((pi7_eig**2).sum())

    eps5 = float(S_cache[5] - Fraction(7, 15))   # exact
    eps6_p = (X6_p - X5_p) - 7/15
    eps6_e = (X6_e - X5_p) - 7/15
    eps7_p = (X7_p - X6_p) - 7/15
    eps7_e = (X7_e - X6_e) - 7/15

    print()
    print("--- final eps comparison ---")
    print(f"  eps_5 (exact):           {eps5:+.10e}")
    print(f"  eps_6 (power iter):      {eps6_p:+.10e}")
    print(f"  eps_6 (scipy.eigs):      {eps6_e:+.10e}")
    print(f"  eps_7 (power iter):      {eps7_p:+.10e}")
    print(f"  eps_7 (scipy.eigs):      {eps7_e:+.10e}")
    print()
    print(f"  |eps_7/eps_6| (power):  {abs(eps7_p/eps6_p):.6f}")
    print(f"  |eps_7/eps_6| (eigs):   {abs(eps7_e/eps6_e):.6f}")


if __name__ == "__main__":
    main()
