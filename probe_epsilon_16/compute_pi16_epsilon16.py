"""
compute_pi16_epsilon16.py — build pi_16 via sparse power iteration, then compute
epsilon_16 via FFT.

Scale: n = 2*3^15 = 28,697,814 states; nnz <= n*V_MAX with V_MAX=60.
K storage ~ 20.7 GB (csr float64). Expect build ~2-3 min, iter ~5-10 min.

Pi saved to probe_self_similarity/pi_16_truncated.npz for future re-use.
"""
from __future__ import annotations
import math, os, sys, time, json
import numpy as np
import scipy.sparse as sp

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_epsilon_16"
PI_OUT = r"C:\Collatz\probe_self_similarity\pi_16_truncated.npz"
os.makedirs(OUT_DIR, exist_ok=True)
V_MAX = 60

EPS_KNOWN = {
    1: +2.0000000000e-01, 2: +9.5238095238e-03, 3: -5.0919863259e-03,
    4: -2.4522582483e-03, 5: -1.1517469151e-03, 6: -4.9790566522e-04,
    7: -1.1752368304e-03, 8: -7.4554636729e-04, 9: -7.5202571564e-06,
    10: +7.2075091711e-04, 11: +1.5019670121e-03, 12: +2.2747137206e-03,
    13: +2.9482473172e-03, 14: +3.5876674275e-03, 15: +4.1611139460e-03,
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


def stationary_sparse(K, tol=1e-13, max_iter=400):
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
        if (it + 1) % 5 == 0:
            print(f"    iter {it+1}: ||delta||_1 = {delta:.3e}", flush=True)
    return pi, max_iter, delta


def main():
    print("=" * 78)
    print("Computing pi_16 + epsilon_16 (sparse, v_max=60)")
    print("=" * 78)

    k = 16
    n_proj = 2 * 3 ** (k - 1)
    nnz_proj = n_proj * V_MAX
    print(f"  n = {n_proj:,}, projected nnz ≈ {nnz_proj:,}, ~{nnz_proj*12/1e9:.1f} GB", flush=True)

    t0 = time.time()
    print("  Building K_16 ...", flush=True)
    K, coprime = build_K_sparse_truncated(k)
    print(f"  K built: shape {K.shape}, nnz {K.nnz:,}, t = {time.time()-t0:.1f}s", flush=True)

    rsum = np.array(K.sum(axis=1)).flatten()
    print(f"  row sums: min {rsum.min():.10f}, max {rsum.max():.10f}", flush=True)

    t0 = time.time()
    print("  Power iteration on K_16 ...", flush=True)
    pi, iters, residual = stationary_sparse(K, tol=1e-13, max_iter=400)
    print(f"  stationary: {iters} iters, ||delta||_1 = {residual:.2e}, t = {time.time()-t0:.1f}s", flush=True)

    print(f"  Saving pi_16 to {PI_OUT}", flush=True)
    np.savez_compressed(PI_OUT, pi=pi, coprime=coprime, k=k, v_max=V_MAX)

    # === S_16 via FFT ===
    print("\n  Computing S_16 via FFT ...", flush=True)
    t0 = time.time()
    N16 = 3 ** k
    pi_full = np.zeros(N16, dtype=np.float64)
    pi_full[coprime] = pi
    pi_hat = np.fft.fft(pi_full)
    xi_arr = np.arange(N16)
    mask_nontrivial = xi_arr % 3 != 0
    S16_fft = float(np.sum(np.abs(pi_hat[mask_nontrivial]) ** 2))
    eps16_fft = S16_fft - 7.0 / 15.0
    t_fft = time.time() - t0
    print(f"  S_16 (FFT) = {S16_fft:.15f}  ({t_fft:.1f}s)", flush=True)
    print(f"  epsilon_16 (FFT) = {eps16_fft:+.12e}", flush=True)

    print(f"\n  eps_13 = {EPS_KNOWN[13]:+.10e}", flush=True)
    print(f"  eps_14 = {EPS_KNOWN[14]:+.10e}", flush=True)
    print(f"  eps_15 = {EPS_KNOWN[15]:+.10e}", flush=True)
    print(f"  eps_16 = {eps16_fft:+.10e}", flush=True)
    print(f"  ratio eps_16/eps_15 = {eps16_fft/EPS_KNOWN[15]:+.6f}", flush=True)
    print(f"  ratio eps_15/eps_14 = {EPS_KNOWN[15]/EPS_KNOWN[14]:+.6f}", flush=True)
    print(f"  ratio eps_14/eps_13 = {EPS_KNOWN[14]/EPS_KNOWN[13]:+.6f}", flush=True)
    print(f"  sign change k=15 -> 16: {'YES' if (eps16_fft > 0) != (EPS_KNOWN[15] > 0) else 'NO'}", flush=True)

    had_13 = 1.0 / abs(EPS_KNOWN[13]) ** (1/13)
    had_14 = 1.0 / abs(EPS_KNOWN[14]) ** (1/14)
    had_15 = 1.0 / abs(EPS_KNOWN[15]) ** (1/15)
    had_16 = 1.0 / abs(eps16_fft) ** (1/16)
    print(f"\n  Hadamard radius |eps_k|^(-1/k):", flush=True)
    print(f"    k=13: {had_13:.4f}", flush=True)
    print(f"    k=14: {had_14:.4f}  (inward {had_13-had_14:.4f})", flush=True)
    print(f"    k=15: {had_15:.4f}  (inward {had_14-had_15:.4f})", flush=True)
    print(f"    k=16: {had_16:.4f}  (inward {had_15-had_16:.4f})  -- predicted asymptote 1.016", flush=True)

    out = {
        "k": 16, "S_16": S16_fft, "epsilon_16": eps16_fft,
        "ratio_eps16_eps15": eps16_fft / EPS_KNOWN[15],
        "hadamard_at_16": had_16,
        "envelope_2to16": abs(eps16_fft) * 2**16,
        "eps_known_through_15": EPS_KNOWN,
        "power_iter_iters": iters,
        "power_iter_residual": residual,
        "fft_time_sec": t_fft,
        "K_nnz": int(K.nnz),
        "n_states": int(K.shape[0]),
        "v_max_truncation": V_MAX,
    }
    with open(os.path.join(OUT_DIR, "epsilon_16_result.json"), "w") as f:
        json.dump(out, f, indent=2)
    with open(os.path.join(OUT_DIR, "S_16_epsilon_16.txt"), "w") as f:
        f.write(f"k=16\nS_16 = {S16_fft:.20e}\nepsilon_16 = {eps16_fft:+.20e}\n")
        f.write(f"|eps_16| * 2^16 = {abs(eps16_fft)*2**16:.6e}\n")
        f.write(f"Hadamard radius at k=16: {had_16:.6f}\n")
    print(f"\nSaved: {OUT_DIR}/epsilon_16_result.json", flush=True)


if __name__ == "__main__":
    main()
