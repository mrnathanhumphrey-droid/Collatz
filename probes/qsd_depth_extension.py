"""
qsd_depth_extension.py — Result 51: extend QSD computation to deeper
state spaces K > 6 (Chang's depth) and test convergence v_K → D_empirical.

For each K in {6, 8, 10, 12, 14, 16, 18}:
  1. Build cylinder-averaged kernel P_K on odd residues mod 2^K
     (state size 2^(K-1)), using 128 lifts for fine cylinder averaging
     within each state.
  2. Identify absorbing states: r mod 32 = 21 (m_3 cylinder).
  3. Eigendecompose P_K_sub^T -> leading eigenvalue lambda_PF and
     eigenvector v_K (the QSD).
  4. Project v_K to mod 32: D_K(r) = sum_{m ≡ r mod 32} v_K(m) / pi_32(r).
  5. Compare to empirical D_avg from Result 50.

Validation: at K=6 the construction should reproduce Chang's pi exactly
and yield Result 50's (b) {21, 53} QSD: lambda_PF = 0.9375, near-uniform v.
"""
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit, prange
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigs

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


# ----------------------------------------------------------------------
# Numba-parallel kernel construction
# ----------------------------------------------------------------------
@njit(parallel=True, cache=True)
def build_P_K_pairs(K, N_LIFTS):
    """
    For state space mod 2^K (n_states = 2^(K-1) odd residues), generate
    transition pairs (row_idx, col_idx) from 128 lifts per state.
    Returns flat arrays of length n_states * N_LIFTS (with duplicates).

    State indexing: state i corresponds to residue r = 2*i + 1 (odd).
    Mod 2^K mask: M_mask = (1 << K) - 1.
    """
    M = np.int64(1) << K            # 2^K
    M_mask = M - 1                  # for fast mod
    n_states = M >> 1               # 2^(K-1)

    rows = np.empty(n_states * N_LIFTS, dtype=np.int64)
    cols = np.empty(n_states * N_LIFTS, dtype=np.int64)

    for i in prange(n_states):
        r = np.int64(2 * i + 1)
        offset = i * N_LIFTS
        for k in range(N_LIFTS):
            m = r + M * k          # lift modulo 2^(K+log2(N_LIFTS))
            three_m = 3 * m + 1
            # Strip trailing zeros (Syracuse: divide by 2^v_2(3m+1))
            while three_m & 1 == 0:
                three_m >>= 1
            r_next = three_m & M_mask  # mod 2^K
            j = (r_next - 1) >> 1
            rows[offset + k] = i
            cols[offset + k] = j

    return rows, cols


# ----------------------------------------------------------------------
# Build sparse P_K and absorbing mask
# ----------------------------------------------------------------------
def build_P_K(K, N_LIFTS=128):
    n_states = 1 << (K - 1)
    log(f"  K={K}: building kernel ({n_states:,} states, {N_LIFTS} lifts each)")
    t0 = time.perf_counter()
    rows, cols = build_P_K_pairs(K, N_LIFTS)
    t_build = time.perf_counter() - t0

    data = np.full(len(rows), 1.0 / N_LIFTS, dtype=np.float64)
    P_coo = coo_matrix(
        (data, (rows, cols)), shape=(n_states, n_states)
    )
    P_csr = P_coo.tocsr()
    P_csr.sum_duplicates()
    nnz = P_csr.nnz

    # Absorbing residues: r mod 32 == 21 -> i mod 16 == 10 (since r = 2i+1)
    # r mod 32 = 21 means r = 21, 53, 85, ... = 21 + 32k for k = 0, 1, 2, ...
    # In state index: i = (r - 1) / 2 = (21 + 32k - 1) / 2 = 10 + 16k
    absorbing_idx = np.where(np.arange(n_states) % 16 == 10)[0]
    n_absorbing = len(absorbing_idx)

    surviving_idx = np.where(np.arange(n_states) % 16 != 10)[0]
    n_surv = len(surviving_idx)

    log(f"    build_pairs: {t_build:.2f}s")
    log(f"    nnz after dedupe: {nnz:,}")
    log(f"    absorbing states (r mod 32 = 21): {n_absorbing:,}")
    log(f"    surviving states: {n_surv:,}")

    return P_csr, surviving_idx, absorbing_idx, n_states


# ----------------------------------------------------------------------
# Eigendecomposition (sparse leading eigenpair)
# ----------------------------------------------------------------------
def compute_qsd(P_csr, surviving_idx, n_top=2):
    """
    P_sub = P_csr restricted to surviving states.
    Compute leading n_top eigenpairs of P_sub^T.
    Return lambda_PF, v_PF (probability normalized), lambda_2, |lambda_2|/lambda_PF.
    """
    P_sub = P_csr[surviving_idx, :][:, surviving_idx]
    n_surv = P_sub.shape[0]

    log(f"    eigs (P_sub^T, k={n_top}, n={n_surv:,})...")
    t0 = time.perf_counter()
    # Largest magnitude eigenvalues. Use sigma=1 shift-invert? No — for k=2
    # leading eigenvalues, which='LM' is fine.
    if n_surv <= 2:
        # eigs requires k < n - 1
        from scipy.sparse.linalg import eigs as _eigs
        # Fall back to dense
        from numpy.linalg import eig
        Pdense = P_sub.toarray().T
        vals_d, vecs_d = eig(Pdense)
        order = np.argsort(-np.abs(vals_d))
        vals = vals_d[order][:n_top]
        vecs = vecs_d[:, order[:n_top]]
    else:
        vals, vecs = eigs(P_sub.T, k=n_top, which='LM', maxiter=20000, tol=1e-10)
        order = np.argsort(-np.abs(vals))
        vals = vals[order]
        vecs = vecs[:, order]
    t_eigs = time.perf_counter() - t0
    log(f"    eigs: {t_eigs:.2f}s")

    lambda_PF = float(vals[0].real)
    v_PF = vecs[:, 0].real
    if v_PF.sum() < 0:
        v_PF = -v_PF
    v_PF = v_PF / v_PF.sum()

    if n_top >= 2:
        lambda_2 = vals[1]
        spec_gap = abs(lambda_2) / lambda_PF
        v_2 = vecs[:, 1].real
        # v_2 has zero sum; normalize by max absolute value
        if abs(v_2).sum() > 0:
            v_2 = v_2 / np.max(np.abs(v_2))
    else:
        lambda_2 = 0
        spec_gap = 0
        v_2 = None

    return lambda_PF, v_PF, lambda_2, spec_gap, v_2


# ----------------------------------------------------------------------
# Project v_K to mod 32 (and compute D_K)
# ----------------------------------------------------------------------
def project_to_mod32(v_PF, surviving_idx, n_states, pi_32):
    """
    v_PF is the QSD on surviving states (probability).
    Embed back to full state space (zeros at absorbing), then sum
    v_full[i] over all i with state-residue r=2i+1 having same r mod 32.
    Divide by pi_32(r mod 32) to get D_K(r mod 32).
    """
    v_full = np.zeros(n_states, dtype=np.float64)
    v_full[surviving_idx] = v_PF

    # Group by r mod 32. r = 2i+1, so i mod 16 indexes the 16 odd residues
    # r = 1, 3, 5, ..., 31.
    # i mod 16 = 0  -> r = 1
    # i mod 16 = 1  -> r = 3
    # i mod 16 = j  -> r = 2j + 1
    bin_idx = np.arange(n_states) % 16
    v32 = np.bincount(bin_idx, weights=v_full, minlength=16)

    odd_r32 = np.array([2 * j + 1 for j in range(16)])
    D_K = {}
    for j, r in enumerate(odd_r32):
        if pi_32[r] > 0 and v32[j] > 0:
            D_K[r] = v32[j] / pi_32[r]
        elif pi_32[r] > 0:
            D_K[r] = 0.0  # absorbing residues
        else:
            D_K[r] = None

    return D_K, v32


# ----------------------------------------------------------------------
# Chang's pi at depth K=6 (mod 32)
# ----------------------------------------------------------------------
def chang_pi_32():
    """Result 50's pi mod 32, which is uniform 1/16 to within ~1.6%."""
    # Recompute exactly using fractions for precision
    from fractions import Fraction
    odd_residues = list(range(1, 64, 2))
    idx = {r: i for i, r in enumerate(odd_residues)}
    P = [[Fraction(0)] * 32 for _ in range(32)]
    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(128):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm & 1 == 0:
                mm >>= 1
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], 128)
    n = 32
    A = [[P[j][i] - (Fraction(1) if i == j else Fraction(0))
          for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = next(row for row in range(col, n) if A[row][col] != 0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    pi = b
    pi_32 = {}
    for r32 in range(1, 32, 2):
        pi_32[r32] = float(pi[idx[r32]] + pi[idx[r32 + 32]])
    return pi_32


# ----------------------------------------------------------------------
# Empirical D_avg from Result 50's late-t data
# ----------------------------------------------------------------------
def load_empirical_D_avg():
    """Load late-t D_avg from qsd_late_t_avg.csv."""
    D_avg = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(",")
            r = int(parts[0])
            D_avg[r] = float(parts[2])
    return D_avg


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    log("=" * 78)
    log("QSD depth-extension: v_K convergence to empirical trajectory measure")
    log("=" * 78)

    # Reference data
    pi_32 = chang_pi_32()
    D_avg = load_empirical_D_avg()
    log("\nChang pi_32 (at K=6, projected):")
    for r in sorted(pi_32.keys()):
        log(f"  r={r:>2}: pi = {pi_32[r]:.6f}, pi*16 = {pi_32[r] * 16:.4f}")

    log("\nEmpirical D_avg (late t, from Result 50):")
    for r in sorted(D_avg.keys()):
        log(f"  r={r:>2}: D_avg = {D_avg[r]:.4f}")

    # Test depths
    K_list = [6, 8, 10, 12, 14, 16, 18]
    results = {}

    for K in K_list:
        log(f"\n{'=' * 78}")
        log(f"K = {K} (state space mod 2^{K} = {1 << K}, "
            f"odd states = {1 << (K - 1):,})")
        log(f"{'=' * 78}")

        try:
            P_csr, surv_idx, abs_idx, n_states = build_P_K(K, N_LIFTS=128)
            lambda_PF, v_PF, lambda_2, spec_gap, v_2 = compute_qsd(
                P_csr, surv_idx, n_top=2
            )
            D_K, v32 = project_to_mod32(v_PF, surv_idx, n_states, pi_32)

            # Compare to empirical
            common_r = [r for r in D_K if D_K[r] is not None
                        and D_K[r] > 0 and r in D_avg]
            total_dev = sum(abs(D_avg[r] - D_K[r]) for r in common_r)
            mean_dev = total_dev / len(common_r) if common_r else float("nan")

            # Also project second eigenmode
            D_K_2, v32_2 = project_to_mod32(
                v_2 if v_2 is not None else np.zeros_like(v_PF),
                surv_idx, n_states, pi_32
            ) if v_2 is not None else ({}, None)

            log(f"\n  lambda_PF = {lambda_PF:.6f}")
            log(f"  |lambda_2| = {abs(lambda_2):.6f}")
            log(f"  spectral gap |lambda_2|/lambda_PF = {spec_gap:.4f}")
            log(f"  total |D_K - D_avg| = {total_dev:.4f} ({len(common_r)} residues)")
            log(f"  mean |D_K - D_avg| per residue = {mean_dev:.4f}")

            log(f"\n  D_K(r) vs empirical D_avg(r):")
            log(f"  {'r':>3}  {'D_K':>7}  {'D_avg':>7}  {'diff':>8}  "
                f"{'D_avg/D_K':>10}")
            for r in sorted(D_K.keys()):
                d_k = D_K[r]
                d_emp = D_avg.get(r)
                if d_k > 0 and d_emp is not None:
                    diff = d_emp - d_k
                    ratio = d_emp / d_k if d_k > 0 else float("nan")
                    log(f"  {r:>3}  {d_k:>7.4f}  {d_emp:>7.4f}  "
                        f"{diff:>+8.4f}  {ratio:>10.4f}")
                elif d_k == 0:
                    log(f"  {r:>3}  {'absorb':>7}  {d_emp:>7.4f}  "
                        f"{'—':>8}  {'—':>10}")

            results[K] = {
                "lambda_PF": lambda_PF,
                "lambda_2": lambda_2,
                "spec_gap": spec_gap,
                "D_K": D_K,
                "D_K_2": D_K_2,
                "total_dev": total_dev,
                "n_states": n_states,
                "n_surv": len(surv_idx),
            }

            # Free memory
            del P_csr, v_PF
            if v_2 is not None:
                del v_2
            import gc
            gc.collect()
        except MemoryError as e:
            log(f"  MEMORY ERROR at K={K}: {e}")
            log(f"  Stopping at K={K-2}")
            break
        except Exception as e:
            log(f"  ERROR at K={K}: {e}")
            log(f"  Stopping at K={K-2}")
            break

    # ------------------------------------------------------------------
    # Convergence diagnostics
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Convergence diagnostics")
    log("=" * 78)

    log("\nlambda_PF evolution:")
    log(f"  {'K':>3}  {'n_states':>10}  {'lambda_PF':>10}  "
        f"{'spec_gap':>9}  {'total_dev':>10}")
    for K, r in sorted(results.items()):
        log(f"  {K:>3}  {r['n_states']:>10,}  {r['lambda_PF']:>10.6f}  "
            f"{r['spec_gap']:>9.4f}  {r['total_dev']:>10.4f}")

    log("\nPer-residue D_K(r) evolution across K:")
    K_done = sorted(results.keys())
    log(f"  {'r':>3}  {'D_avg':>7}  " +
        "  ".join(f"K={K:>2}" for K in K_done))
    for r in sorted(D_avg.keys()):
        line = f"  {r:>3}  {D_avg[r]:>7.4f}"
        for K in K_done:
            d_k = results[K]["D_K"].get(r)
            if d_k is None or d_k == 0:
                line += f"  {'—':>5}"
            else:
                line += f"  {d_k:>5.3f}"
        log(line)

    # Save CSVs
    log("\n[Save] CSVs")
    with open(OUT / "qsd_depth_extension.csv", "w") as f:
        f.write("K,r,D_K,D_avg,abs_diff\n")
        for K in sorted(results.keys()):
            for r in sorted(D_avg.keys()):
                d_k = results[K]["D_K"].get(r)
                d_emp = D_avg.get(r)
                d_str = f"{d_k:.6f}" if d_k is not None else ""
                diff_str = (f"{abs(d_emp - d_k):.6f}"
                            if d_k is not None and d_k > 0 else "")
                f.write(f"{K},{r},{d_str},{d_emp:.6f},{diff_str}\n")

    with open(OUT / "qsd_lambda_evolution.csv", "w") as f:
        f.write("K,n_states,n_surviving,lambda_PF,lambda_2,spec_gap,total_dev\n")
        for K, r in sorted(results.items()):
            f.write(f"{K},{r['n_states']},{r['n_surv']},"
                    f"{r['lambda_PF']:.10f},{r['lambda_2']:.10f},"
                    f"{r['spec_gap']:.6f},{r['total_dev']:.6f}\n")

    with open(OUT / "qsd_second_eigenmode.csv", "w") as f:
        f.write("K,r,D_K_2nd_mode\n")
        for K in sorted(results.keys()):
            if results[K]["D_K_2"]:
                for r in sorted(results[K]["D_K_2"].keys()):
                    d_2 = results[K]["D_K_2"][r]
                    if d_2 is not None:
                        f.write(f"{K},{r},{d_2:.6f}\n")

    (OUT / "qsd_depth_extension_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("[wrote] qsd_depth_extension.csv, qsd_lambda_evolution.csv, "
        "qsd_second_eigenmode.csv, qsd_depth_extension_log.txt")


if __name__ == "__main__":
    main()
