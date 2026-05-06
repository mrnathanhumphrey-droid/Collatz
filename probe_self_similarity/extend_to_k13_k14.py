"""
extend_to_k13_k14.py — compute π_13 and π_14, then re-run the entropy-deficit
analysis to see whether ρ_Δ saturates near 0.83 or keeps drifting upward.

Strategy: full dense K_k is intractable beyond k=10 (k=13 would need ~9 TB).
Use sparse K_k truncated at v_max = 60 (since 2^{-60} ≈ 8.7e-19, well below
double-precision threshold; tail contribution is negligible). For k=13,
sparse K has ~50 nonzeros per row × 1M rows = 50M nonzeros, manageable.

Time budget per k:
  k=13: build ~30s, power-iter ~1 min
  k=14: build ~2 min, power-iter ~3 min  (3.2M rows, 160M nnz)
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_self_similarity"
PROFINITE_DIR = r"C:\Collatz\probe_profinite"

V_MAX = 60  # truncation: weight(v) = 2^{-v} drops below machine epsilon by v=60


def order_of_two(N):
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


def build_K_sparse_truncated(k, v_max=V_MAX):
    """Build K_k as scipy.sparse.csr_matrix using only v ∈ {1..v_max}.
    Truncation error per matvec ~ 2^{-v_max} (below 1e-18 for v_max=60).
    """
    N = 3 ** k
    M = order_of_two(N)
    v_eff = min(M, v_max)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(v_eff, dtype=np.int64)
    p = inv2
    for v in range(v_eff):
        powers_inv2[v] = p
        p = (p * inv2) % N
    # coprime states (skip multiples of 3)
    coprime_mask = np.ones(N, dtype=bool)
    coprime_mask[::3] = False  # r=0,3,6,...
    coprime_idx_in_N = np.where(coprime_mask)[0]
    n = len(coprime_idx_in_N)
    state_idx = -np.ones(N, dtype=np.int64)
    state_idx[coprime_idx_in_N] = np.arange(n)
    weights = np.array([2.0 ** -(v + 1) for v in range(v_eff)], dtype=np.float64)
    weights /= weights.sum()  # normalize on truncated support

    # COO assembly
    rows = np.empty(n * v_eff, dtype=np.int64)
    cols = np.empty(n * v_eff, dtype=np.int64)
    vals = np.empty(n * v_eff, dtype=np.float64)
    base = (3 * coprime_idx_in_N + 1) % N  # length n
    # For each v: target = base * powers_inv2[v] mod N → state_idx
    for v in range(v_eff):
        targets = (base * powers_inv2[v]) % N
        # All targets coprime to 3 since (3r+1)·2^{-v} mod 3 = 1·2^{-v} mod 3 ≠ 0
        rows[v * n:(v + 1) * n] = np.arange(n)
        cols[v * n:(v + 1) * n] = state_idx[targets]
        vals[v * n:(v + 1) * n] = weights[v]
    K = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    K.sum_duplicates()
    return K, coprime_idx_in_N


def stationary_sparse(K, tol=1e-13, max_iter=200):
    """Power iteration on sparse K acting on row-vectors via π_new = π @ K."""
    n = K.shape[0]
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    K_T = K.T.tocsr()  # for column-vector form: π_new^T = K^T @ π^T
    for it in range(max_iter):
        pi_new = K_T @ pi  # equivalent to π @ K
        s = pi_new.sum()
        if s != 0:
            pi_new /= s
        delta = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        if delta < tol:
            return pi, it + 1, delta
    return pi, max_iter, delta


def shannon_entropy(pi):
    p = pi[pi > 0]
    return float(-np.sum(p * np.log(p)))


def main():
    print("=" * 78)
    print(f"Extending entropy series to k=13, 14 (sparse K_k, v truncated at "
          f"v_max={V_MAX})")
    print("=" * 78)
    print()

    # Existing entropy values from cache
    H_cache = {
        5: 4.607200, 6: 5.656623, 7: 6.715097, 8: 7.780191,
        9: 8.850412, 10: 9.924671, 11: 11.002234, 12: 12.082527,
    }

    new_entropies = {}
    for k in [13, 14]:
        print(f"--- k = {k} ---")
        # Resource projection
        n_proj = 2 * 3 ** (k - 1)
        nnz_proj = n_proj * V_MAX
        gb_K = nnz_proj * 12 / 1e9  # ~12 bytes per nnz (8 val + 4 col)
        print(f"  n = {n_proj:,}, projected nnz ≈ {nnz_proj:,}, "
              f"K storage ~{gb_K:.2f} GB")
        t0 = time.time()
        K, coprime = build_K_sparse_truncated(k)
        print(f"  K built: shape {K.shape}, nnz {K.nnz:,}, "
              f"t = {time.time()-t0:.1f}s")
        # Sanity: row sums should be 1.0
        rsum = np.array(K.sum(axis=1)).flatten()
        print(f"  row sums: min {rsum.min():.10f}, max {rsum.max():.10f}, "
              f"mean {rsum.mean():.10f}")
        t0 = time.time()
        pi, iters, residual = stationary_sparse(K, tol=1e-13, max_iter=300)
        print(f"  stationary: {iters} iters, ||δ||_1 = {residual:.2e}, "
              f"t = {time.time()-t0:.1f}s")
        # Entropy
        H = shannon_entropy(pi)
        H_uni = math.log(len(pi))
        print(f"  H(π_{k}) = {H:.6f}, H_uniform = {H_uni:.6f}, "
              f"D_KL = {H_uni - H:.6f}")
        new_entropies[k] = H
        # Cache π_k
        out_pi = os.path.join(OUT_DIR, f"pi_{k}_truncated.npz")
        np.savez_compressed(out_pi, pi=pi, coprime=coprime, k=k, v_max=V_MAX)
        print(f"  saved {out_pi}")
        print()

    # Combined entropy series
    H_all = dict(H_cache)
    H_all.update(new_entropies)
    print("=" * 78)
    print("Extended entropy series (k=5..14)")
    print("=" * 78)
    print()
    LOG3 = math.log(3.0)
    ks = sorted(H_all.keys())
    print(f"  {'k':>3} {'H(π_k)':>14} {'ΔH_k':>14} "
          f"{'Δ_k = log3 - ΔH_k':>20} {'ratio':>10}")
    Delta = {}
    for i, k in enumerate(ks[:-1]):
        kp1 = ks[i + 1]
        dH = H_all[kp1] - H_all[k]
        delta_k = LOG3 - dH
        Delta[k] = delta_k
        ratio_str = ""
        if i > 0:
            prev_delta = Delta[ks[i - 1]]
            if prev_delta > 0:
                ratio_str = f"{delta_k / prev_delta:.6f}"
        print(f"  {k:>3} {H_all[k]:>14.6f} {dH:>14.10f} "
              f"{delta_k:>20.10e} {ratio_str:>10}")
    print(f"  {ks[-1]:>3} {H_all[ks[-1]]:>14.6f} {'—':>14} {'—':>20} {'—':>10}")

    # Re-fit ρ_Δ on extended series
    print()
    print("=" * 78)
    print(f"Geometric fit log Δ_k = a + b·k on k = {ks[0]}..{ks[-2]}")
    print("=" * 78)
    print()
    Δ_arr = np.array([Delta[k] for k in sorted(Delta.keys())])
    k_arr = np.array(sorted(Delta.keys()), dtype=np.float64)
    log_Δ = np.log(Δ_arr)
    A = np.column_stack([np.ones_like(k_arr), k_arr])
    coeffs, _, _, _ = np.linalg.lstsq(A, log_Δ, rcond=None)
    a_fit, b_fit = coeffs
    rho_full = math.exp(b_fit)
    pred = A @ coeffs
    ss_res = float(np.sum((log_Δ - pred) ** 2))
    ss_tot = float(np.sum((log_Δ - log_Δ.mean()) ** 2))
    r2_full = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    print(f"  Full fit (k=5..{ks[-2]}, {len(k_arr)} points):  "
          f"ρ_Δ = {rho_full:.6f}, R² = {r2_full:.6f}")

    # Late fit (last 5 points only) to see if rate has stabilized
    if len(k_arr) >= 6:
        late_k = k_arr[-5:]
        late_log_Δ = log_Δ[-5:]
        A_late = np.column_stack([np.ones_like(late_k), late_k])
        c_late, _, _, _ = np.linalg.lstsq(A_late, late_log_Δ, rcond=None)
        rho_late = math.exp(c_late[1])
        pred_late = A_late @ c_late
        r2_late = 1.0 - np.sum((late_log_Δ - pred_late) ** 2) / np.sum(
            (late_log_Δ - late_log_Δ.mean()) ** 2)
        print(f"  Late fit (last 5 transitions, k={int(late_k[0])}..{int(late_k[-1])}):  "
              f"ρ_Δ = {rho_late:.6f}, R² = {r2_late:.6f}")

    # Compare to ρ_slow
    rho_slow = 0.826934
    print()
    print(f"  ρ_slow (order-3 recurrence) = {rho_slow:.6f}")
    print(f"  ρ_Δ (full fit) = {rho_full:.6f}, gap = {rho_full - rho_slow:+.6f}")
    if len(k_arr) >= 6:
        print(f"  ρ_Δ (late fit) = {rho_late:.6f}, gap = {rho_late - rho_slow:+.6f}")
    # Trend in ratios
    ratios = [Delta[ks[i + 1]] / Delta[ks[i]]
              for i in range(len(ks) - 2)
              if Delta[ks[i]] > 0]
    print()
    print(f"  Ratio sequence Δ_{{k+1}}/Δ_k:  "
          + ", ".join(f"{r:.4f}" for r in ratios))
    print(f"  Trend: {'rising monotonically' if all(ratios[i+1] > ratios[i] for i in range(len(ratios)-1)) else 'mixed/stabilizing'}")


if __name__ == "__main__":
    main()
