"""
compute_pi15_epsilon15.py — build pi_15 via sparse power iteration, then compute
epsilon_15 via FFT (same path as epsilon_13, epsilon_14).

Scale: n = 2*3^14 = 9,565,938 states; nnz <= n*V_MAX with V_MAX=60.
K storage ~ 6.9 GB (csr float64). Power iteration ~ matvecs of 574M nnz each;
expect ~5 min build + 5-15 min iter.

Pi saved to probe_self_similarity/pi_15_truncated.npz for future re-use.
"""
from __future__ import annotations
import math, os, sys, time, json
import numpy as np
import scipy.sparse as sp

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_epsilon_15"
PI_OUT = r"C:\Collatz\probe_self_similarity\pi_15_truncated.npz"
os.makedirs(OUT_DIR, exist_ok=True)
V_MAX = 60

EPS_KNOWN = {
    1: +2.0000000000e-01, 2: +9.5238095238e-03, 3: -5.0919863259e-03,
    4: -2.4522582483e-03, 5: -1.1517469151e-03, 6: -4.9790566522e-04,
    7: -1.1752368304e-03, 8: -7.4554636729e-04, 9: -7.5202571564e-06,
    10: +7.2075091711e-04, 11: +1.5019670121e-03, 12: +2.2747137206e-03,
    13: +2.9482473172e-03, 14: +3.5876674275e-03,
}


def order_of_two(N):
    assert N % 2 == 1
    m, v = 1, 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def build_K_sparse_truncated(k, v_max=V_MAX):
    N = 3 ** k
    M = order_of_two(N)
    v_eff = min(M, v_max)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(v_eff, dtype=np.int64)
    p = inv2
    for v in range(v_eff):
        powers_inv2[v] = p
        p = (p * inv2) % N
    coprime_mask = np.ones(N, dtype=bool)
    coprime_mask[::3] = False
    coprime_idx_in_N = np.where(coprime_mask)[0]
    n = len(coprime_idx_in_N)
    state_idx = -np.ones(N, dtype=np.int64)
    state_idx[coprime_idx_in_N] = np.arange(n)
    weights = np.array([2.0 ** -(v + 1) for v in range(v_eff)], dtype=np.float64)
    weights /= weights.sum()

    rows = np.empty(n * v_eff, dtype=np.int64)
    cols = np.empty(n * v_eff, dtype=np.int64)
    vals = np.empty(n * v_eff, dtype=np.float64)
    base = (3 * coprime_idx_in_N + 1) % N
    for v in range(v_eff):
        targets = (base * powers_inv2[v]) % N
        rows[v * n:(v + 1) * n] = np.arange(n)
        cols[v * n:(v + 1) * n] = state_idx[targets]
        vals[v * n:(v + 1) * n] = weights[v]
    K = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    K.sum_duplicates()
    return K, coprime_idx_in_N


def stationary_sparse(K, tol=1e-13, max_iter=300):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    K_T = K.T.tocsr()
    for it in range(max_iter):
        pi_new = K_T @ pi
        s = pi_new.sum()
        if s != 0:
            pi_new /= s
        delta = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        if delta < tol:
            return pi, it + 1, delta
        if (it + 1) % 10 == 0:
            print(f"    iter {it+1}: ||delta||_1 = {delta:.3e}")
    return pi, max_iter, delta


def main():
    print("=" * 78)
    print("Computing pi_15 + epsilon_15 (sparse, v_max=60)")
    print("=" * 78)

    k = 15
    n_proj = 2 * 3 ** (k - 1)
    nnz_proj = n_proj * V_MAX
    print(f"  n = {n_proj:,}, projected nnz ≈ {nnz_proj:,}, ~{nnz_proj*12/1e9:.1f} GB")

    t0 = time.time()
    print("  Building K_15 ...")
    K, coprime = build_K_sparse_truncated(k)
    print(f"  K built: shape {K.shape}, nnz {K.nnz:,}, t = {time.time()-t0:.1f}s")

    rsum = np.array(K.sum(axis=1)).flatten()
    print(f"  row sums: min {rsum.min():.10f}, max {rsum.max():.10f}")

    t0 = time.time()
    print("  Power iteration on K_15 ...")
    pi, iters, residual = stationary_sparse(K, tol=1e-13, max_iter=400)
    print(f"  stationary: {iters} iters, ||delta||_1 = {residual:.2e}, t = {time.time()-t0:.1f}s")

    print(f"  Saving pi_15 to {PI_OUT}")
    np.savez_compressed(PI_OUT, pi=pi, coprime=coprime, k=k, v_max=V_MAX)

    # === S_15 via FFT ===
    print("\n  Computing S_15 via FFT ...")
    t0 = time.time()
    N15 = 3 ** k
    pi_full = np.zeros(N15, dtype=np.float64)
    pi_full[coprime] = pi
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N15)
    mask_nontrivial = xi_arr % 3 != 0
    S15_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps15_fft = S15_fft - 7.0 / 15.0
    t_fft = time.time() - t0
    print(f"  S_15 (FFT) = {S15_fft:.15f}  ({t_fft:.1f}s)")
    print(f"  epsilon_15 (FFT) = {eps15_fft:+.12e}")

    # Compare with k=13, 14
    print(f"\n  eps_13 = {EPS_KNOWN[13]:+.10e}")
    print(f"  eps_14 = {EPS_KNOWN[14]:+.10e}")
    print(f"  eps_15 = {eps15_fft:+.10e}")
    print(f"  ratio eps_15/eps_14 = {eps15_fft/EPS_KNOWN[14]:+.6f}")
    print(f"  ratio eps_14/eps_13 = {EPS_KNOWN[14]/EPS_KNOWN[13]:+.6f}")
    print(f"  sign change k=14 -> 15: {'YES' if (eps15_fft > 0) != (EPS_KNOWN[14] > 0) else 'NO'}")

    # Hadamard radii
    had_13 = 1.0 / abs(EPS_KNOWN[13]) ** (1/13)
    had_14 = 1.0 / abs(EPS_KNOWN[14]) ** (1/14)
    had_15 = 1.0 / abs(eps15_fft) ** (1/15)
    print(f"\n  Hadamard radius |eps_k|^(-1/k):")
    print(f"    k=13: {had_13:.4f}")
    print(f"    k=14: {had_14:.4f}")
    print(f"    k=15: {had_15:.4f}  (predicted asymptote 1.016)")

    out = {
        "k": 15, "S_15": S15_fft, "epsilon_15": eps15_fft,
        "ratio_eps15_eps14": eps15_fft / EPS_KNOWN[14],
        "hadamard_at_15": had_15,
        "envelope_2to15": abs(eps15_fft) * 2**15,
        "eps_known_through_14": EPS_KNOWN,
        "power_iter_iters": iters,
        "power_iter_residual": residual,
        "fft_time_sec": t_fft,
        "K_nnz": int(K.nnz),
        "n_states": int(K.shape[0]),
        "v_max_truncation": V_MAX,
    }
    with open(os.path.join(OUT_DIR, "epsilon_15_result.json"), "w") as f:
        json.dump(out, f, indent=2)
    with open(os.path.join(OUT_DIR, "S_15_epsilon_15.txt"), "w") as f:
        f.write(f"k=15\nS_15 = {S15_fft:.20e}\nepsilon_15 = {eps15_fft:+.20e}\n")
        f.write(f"|eps_15| * 2^15 = {abs(eps15_fft)*2**15:.6e}\n")
        f.write(f"Hadamard radius at k=15: {had_15:.6f}\n")
    print(f"\nSaved: {OUT_DIR}/epsilon_15_result.json")


if __name__ == "__main__":
    main()
