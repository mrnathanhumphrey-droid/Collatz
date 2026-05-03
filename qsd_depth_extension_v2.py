"""
qsd_depth_extension_v2.py — Round 2: test alternative absorption conventions.

Round 1 finding: cylinder absorption at "r mod 32 = 21" gives D_K invariant
in K (D_K projected to mod 32 = Chang's depth-6 QSD at every K). This is
because the absorption is at a depth-6 cylinder; deeper kernels project
back to the same depth-6 chain.

To test the depth-extension framework with non-trivial K-dependence, this
round uses SPECIFIC-VALUE absorption (m_j attractor specific residues at
depth K, NOT the depth-6 cylinder).

Three absorption conventions:
  (A) "values": absorb at specific m_j values {1, 5, 21, 85, ..., m_j_max}
      where m_j < 2^K. As K grows, more m_j fit; absorption set grows
      slowly with K.
  (B) "value21": absorb only at residue 21 mod 2^K (ONE specific residue
      at each depth). Single absorbing state at deepest cylinder.
  (C) "value1": absorb only at residue 1 mod 2^K.

For each, project v_K to mod 32 and check whether D_K depends on K and
whether it converges to D_avg.

Decisive outcomes:
  (alpha) D_K -> D_avg as K grows: depth-extension framework correct
          with value-based absorption
  (beta)  D_K depends on K but doesn't reach D_avg
  (gamma) D_K invariant in K for these conventions too: framework wrong
"""
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit, prange
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


@njit(parallel=True, cache=True)
def build_P_K_pairs(K, N_LIFTS):
    M = np.int64(1) << K
    M_mask = M - 1
    n_states = M >> 1
    rows = np.empty(n_states * N_LIFTS, dtype=np.int64)
    cols = np.empty(n_states * N_LIFTS, dtype=np.int64)
    for i in prange(n_states):
        r = np.int64(2 * i + 1)
        offset = i * N_LIFTS
        for k in range(N_LIFTS):
            m = r + M * k
            three_m = 3 * m + 1
            while three_m & 1 == 0:
                three_m >>= 1
            r_next = three_m & M_mask
            j = (r_next - 1) >> 1
            rows[offset + k] = i
            cols[offset + k] = j
    return rows, cols


def build_P_K(K, N_LIFTS=128):
    n_states = 1 << (K - 1)
    log(f"  K={K}: building kernel ({n_states:,} states, {N_LIFTS} lifts)")
    t0 = time.perf_counter()
    rows, cols = build_P_K_pairs(K, N_LIFTS)
    data = np.full(len(rows), 1.0 / N_LIFTS, dtype=np.float64)
    P_csr = coo_matrix(
        (data, (rows, cols)), shape=(n_states, n_states)
    ).tocsr()
    P_csr.sum_duplicates()
    log(f"    nnz: {P_csr.nnz:,}, build_time: {time.perf_counter() - t0:.2f}s")
    return P_csr, n_states


def get_absorbing_set(K, convention):
    """Return list of absorbing state indices at depth K."""
    M = 1 << K
    n_states = 1 << (K - 1)

    if convention == "values":
        # Absorb at specific m_j values that fit in mod 2^K
        m_j = []
        m = 1
        j = 0
        while m < M:
            m_j.append(m)
            j += 1
            m = (4 ** (j + 1) - 1) // 3
        # State index: r odd, i = (r-1)//2
        absorbing_idx = [(m - 1) // 2 for m in m_j if m % 2 == 1]
        label = f"values{{{','.join(str(m) for m in m_j[:3])},...}} ({len(absorbing_idx)} states)"

    elif convention == "value21":
        # Absorb only at residue 21 mod 2^K (single state)
        absorbing_idx = [(21 - 1) // 2] if 21 < M else []
        label = "value21"

    elif convention == "value1":
        # Absorb only at residue 1 mod 2^K (single state)
        absorbing_idx = [(1 - 1) // 2]
        label = "value1"

    elif convention == "cylinder21":
        # Absorb at all r mod 32 = 21 (Round 1 default)
        absorbing_idx = [i for i in range(n_states) if (2 * i + 1) % 32 == 21]
        label = "cylinder21"

    elif convention == "cylinder5":
        # Absorb at all r mod 32 = 5
        absorbing_idx = [i for i in range(n_states) if (2 * i + 1) % 32 == 5]
        label = "cylinder5"

    elif convention == "fine_cylinder":
        # Absorb at residue 21 mod 2^min(K, 14): a finer cylinder than mod 32
        # As K grows beyond 14, this becomes specific value 21 mod 16384.
        K_cyl = min(K, 14)
        M_cyl = 1 << K_cyl
        absorbing_idx = [
            i for i in range(n_states)
            if (2 * i + 1) % M_cyl == 21
        ]
        label = f"residue21mod2^{K_cyl}"

    else:
        raise ValueError(f"Unknown convention: {convention}")

    return absorbing_idx, label


def compute_qsd(P_csr, surv_idx, n_top=2):
    P_sub = P_csr[surv_idx, :][:, surv_idx]
    n_surv = P_sub.shape[0]
    if n_surv <= 3:
        from numpy.linalg import eig
        Pdense = P_sub.toarray().T
        vals_d, vecs_d = eig(Pdense)
        order = np.argsort(-np.abs(vals_d))
        vals = vals_d[order][:n_top]
        vecs = vecs_d[:, order[:n_top]]
    else:
        vals, vecs = eigs(P_sub.T, k=min(n_top, n_surv - 2),
                          which='LM', maxiter=20000, tol=1e-10)
        order = np.argsort(-np.abs(vals))
        vals = vals[order]
        vecs = vecs[:, order]
    lambda_PF = float(vals[0].real)
    v_PF = vecs[:, 0].real
    if v_PF.sum() < 0:
        v_PF = -v_PF
    if v_PF.sum() > 0:
        v_PF = v_PF / v_PF.sum()
    if len(vals) >= 2:
        lambda_2 = vals[1]
    else:
        lambda_2 = 0
    spec_gap = abs(lambda_2) / lambda_PF if lambda_PF > 0 else 0
    return lambda_PF, v_PF, lambda_2, spec_gap


def project_to_mod32(v_PF, surv_idx, n_states, pi_32):
    v_full = np.zeros(n_states, dtype=np.float64)
    v_full[surv_idx] = v_PF
    bin_idx = np.arange(n_states) % 16
    v32 = np.bincount(bin_idx, weights=v_full, minlength=16)
    odd_r32 = np.array([2 * j + 1 for j in range(16)])
    D_K = {}
    for j, r in enumerate(odd_r32):
        if pi_32[r] > 0 and v32[j] > 0:
            D_K[r] = v32[j] / pi_32[r]
        else:
            D_K[r] = 0.0
    return D_K, v32


def chang_pi_32():
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
    pi_32 = {}
    for r32 in range(1, 32, 2):
        pi_32[r32] = float(b[idx[r32]] + b[idx[r32 + 32]])
    return pi_32


def load_empirical_D_avg():
    D_avg = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            r = int(parts[0])
            D_avg[r] = float(parts[2])
    return D_avg


def main():
    log("=" * 78)
    log("QSD depth-extension v2: test alternative absorption conventions")
    log("=" * 78)

    pi_32 = chang_pi_32()
    D_avg = load_empirical_D_avg()

    # Test multiple conventions across K
    K_list = [6, 8, 10, 12, 14]
    conventions = ["cylinder21", "values", "value21", "value1",
                   "fine_cylinder"]

    log("\nTest matrix: each cell shows lambda_PF, total |D_K - D_avg|")
    log("  (absorbing set varies by convention)")

    all_results = {}
    for K in K_list:
        try:
            P_csr, n_states = build_P_K(K, N_LIFTS=128)
            for convention in conventions:
                absorbing_idx, label = get_absorbing_set(K, convention)
                if not absorbing_idx:
                    log(f"  K={K} {convention}: no absorbing states; skip")
                    continue
                surv_idx = np.array([
                    i for i in range(n_states) if i not in set(absorbing_idx)
                ])
                lambda_PF, v_PF, lambda_2, spec_gap = compute_qsd(
                    P_csr, surv_idx
                )
                D_K, v32 = project_to_mod32(v_PF, surv_idx, n_states, pi_32)
                common_r = [r for r in D_K if D_K[r] > 0 and r in D_avg]
                total_dev = sum(abs(D_avg[r] - D_K[r]) for r in common_r)
                all_results[(K, convention)] = {
                    "lambda_PF": lambda_PF,
                    "lambda_2": lambda_2,
                    "spec_gap": spec_gap,
                    "D_K": D_K,
                    "total_dev": total_dev,
                    "n_absorbing": len(absorbing_idx),
                    "label": label,
                }
                log(f"  K={K:>2}  {convention:<14}  "
                    f"abs_states={len(absorbing_idx):>5}  "
                    f"lambda_PF={lambda_PF:.6f}  "
                    f"|lam2|/PF={spec_gap:.4f}  "
                    f"sum|D-D_avg|={total_dev:.4f}")
            del P_csr
            import gc
            gc.collect()
        except Exception as e:
            log(f"  ERROR at K={K}: {e}")

    # Per-residue evolution for the most informative convention
    log("\n" + "=" * 78)
    log("D_K(r) evolution across K for each convention")
    log("=" * 78)

    for convention in conventions:
        log(f"\n--- Convention: {convention} ---")
        log(f"  {'r':>3}  {'D_avg':>7}  " +
            "  ".join(f"K={K:>2}" for K in K_list))
        for r in sorted(D_avg.keys()):
            line = f"  {r:>3}  {D_avg[r]:>7.4f}"
            for K in K_list:
                key = (K, convention)
                if key in all_results:
                    d_k = all_results[key]["D_K"].get(r)
                    if d_k is None or d_k == 0:
                        line += f"  {'absorb' if d_k == 0 else '—':>5}"
                    else:
                        line += f"  {d_k:>5.3f}"
                else:
                    line += f"  {'—':>5}"
            log(line)

    # Save CSV
    with open(OUT / "qsd_depth_v2.csv", "w") as f:
        f.write("K,convention,r,D_K,D_avg,abs_diff,lambda_PF,n_absorbing\n")
        for (K, conv), r in all_results.items():
            for resid in sorted(D_avg.keys()):
                d_k = r["D_K"].get(resid)
                d_emp = D_avg[resid]
                d_str = f"{d_k:.6f}" if d_k is not None else ""
                diff_str = (f"{abs(d_emp - d_k):.6f}"
                            if d_k is not None and d_k > 0 else "")
                f.write(f"{K},{conv},{resid},{d_str},{d_emp:.6f},"
                        f"{diff_str},{r['lambda_PF']:.6f},"
                        f"{r['n_absorbing']}\n")

    (OUT / "qsd_depth_v2_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("[wrote] qsd_depth_v2.csv, qsd_depth_v2_log.txt")


if __name__ == "__main__":
    main()
