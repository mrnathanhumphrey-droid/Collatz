"""
size_stratified_markov.py — Result 54 candidate.

Test Family D: Markov kernel on joint state (residue mod 32, log-size bin).
Compute Perron-Frobenius eigenvector v_PF of empirical sub-stochastic K,
marginalize to residue, compare to D_avg.

Motivation: residue-only kernels (cylinder R51, inverse tree R52) systematically
miss D_avg. Size-stratified explicitly tracks value information through the
log-size coordinate b = floor(log_2(m)).

State: (r, b)  with r in {1, 3, ..., 31} (16 odd residues), b in {0, ..., B-1}.
B = 64 covers m up to 2^64 (well past int64 risk for our N=2^32 orbits).

Kernel: K[(r, b) -> (r', b')] = #{Syracuse transitions (r,b)->(r',b')} / #visits(r,b).
Sub-stochastic because m=1 is absorbing; row sums leak when last step before m=1.

Outputs:
  size_stratified_kernel.npz       -- sparse K_sub
  size_stratified_eigvec.csv       -- v_PF on (r, b) state space
  size_stratified_residue_marginal.csv -- D_predicted vs D_empirical
  size_stratified_markov.md        -- writeup
"""
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit, prange
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)

LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)

R = 16   # odd residues mod 32
B = 64   # log-size bins (0..63)
N_STATES = R * B   # 1024


@njit(parallel=True, cache=True)
def walk_count_transitions(starts, max_T, n_chunks, R, B):
    """
    Walk Collatz orbits via Syracuse map. For each (odd m -> next odd m')
    transition, increment K_counts[from_state, to_state] where
        state = ((r-1)//2) * B + b,  r = m mod 32,  b = min(floor(log_2 m), B-1)

    Each thread chunk maintains its own count matrix; reduce on host.

    Returns:
      K_counts (n_chunks, R*B, R*B) int64
      sigma_arr (n,)     -- total Collatz step count if orbit completed
      T_arr    (n,)      -- total Syracuse step count
      sigma_log_n (n, 2) -- (sigma, log n_start) for sigma-band assignment later
    """
    n = len(starts)
    n_states = R * B
    chunk_size = (n + n_chunks - 1) // n_chunks
    K_counts = np.zeros((n_chunks, n_states, n_states), dtype=np.int64)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    log_n_arr = np.zeros(n, dtype=np.float64)

    for chunk in prange(n_chunks):
        i_start = chunk * chunk_size
        i_end = min((chunk + 1) * chunk_size, n)
        for i in range(i_start, i_end):
            m = np.int64(starts[i])
            log_n_arr[i] = np.log(np.float64(m))
            sigma = 0
            T = 0
            failed = False

            # Strip leading even factors so we land on odd m
            while (m & 1) == 0 and m > 1:
                m >>= 1
                sigma += 1
            if m == 1:
                sigma_arr[i] = sigma
                continue

            # Compute initial state (r, b)
            r_idx = (np.int64(m & 31) - 1) >> 1
            mm = m
            b = 0
            while mm > 1:
                mm >>= 1
                b += 1
            if b >= B:
                b = B - 1
            s_curr = r_idx * B + b

            while m != 1 and T < max_T:
                if m > MAX_VAL // 3:
                    failed = True
                    break
                x = 3 * m + 1
                v = 0
                while (x & 1) == 0:
                    x >>= 1
                    v += 1

                # Compute next state (r', b')
                r_idx_n = (np.int64(x & 31) - 1) >> 1
                xx = x
                bn = 0
                while xx > 1:
                    xx >>= 1
                    bn += 1
                if bn >= B:
                    bn = B - 1
                s_next = r_idx_n * B + bn

                K_counts[chunk, s_curr, s_next] += 1

                sigma += 1 + v
                T += 1
                m = x
                s_curr = s_next

            if not failed and m == 1:
                sigma_arr[i] = sigma
                T_arr[i] = T
    return K_counts, sigma_arr, T_arr, log_n_arr


def chang_pi_32():
    """Reproduce Chang's stationary pi on residues mod 32 (depth-13 cylinder)."""
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
    pi_32 = {r32: float(b[idx[r32]] + b[idx[r32 + 32]])
             for r32 in range(1, 32, 2)}
    return pi_32


def load_empirical_D_avg():
    D_avg = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            D_avg[int(parts[0])] = float(parts[2])
    return D_avg


def main():
    log("=" * 78)
    log("Size-stratified Markov framework (Result 54 candidate)")
    log("State = (r mod 32, b) with b = floor(log_2 m), B=64")
    log("=" * 78)

    pi_32 = chang_pi_32()
    D_avg = load_empirical_D_avg()
    odd_r32 = list(range(1, 32, 2))

    # ---------------------------------------------------------------
    # Step 1: Walk orbits, count (r,b) -> (r',b') transitions
    # ---------------------------------------------------------------
    log2N = 32
    N = 1 << log2N
    n_orbits_per_seed = 500_000
    seeds = [42, 137, 271]
    total_orbits = len(seeds) * n_orbits_per_seed
    log(f"\nStep 1: walking {total_orbits:,} orbits at N=2^{log2N}, B={B}")
    log(f"        state space size = {R}*{B} = {N_STATES}")

    n_chunks = 12
    log(f"        n_chunks = {n_chunks} (parallel accumulators)")

    K_total_counts = np.zeros((N_STATES, N_STATES), dtype=np.int64)
    all_sigma = []
    all_log_n = []

    t0 = time.time()
    for seed_i, seed in enumerate(seeds):
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2,
                                  size=n_orbits_per_seed,
                                  dtype=np.int64) + 1
        K_counts, sigma_arr, T_arr, log_n_arr = walk_count_transitions(
            starts, 600, n_chunks, R, B
        )
        K_seed = K_counts.sum(axis=0)
        K_total_counts += K_seed
        ok = sigma_arr > 0
        all_sigma.append(sigma_arr[ok].astype(np.float64))
        all_log_n.append(log_n_arr[ok])
        log(f"  seed {seed_i+1}/{len(seeds)}: ok={int(ok.sum()):,}, "
            f"transitions={int(K_seed.sum()):,}, "
            f"elapsed={time.time()-t0:.1f}s")

    sigma = np.concatenate(all_sigma)
    log_n = np.concatenate(all_log_n)
    n_orbits_ok = len(sigma)
    total_transitions = int(K_total_counts.sum())
    log(f"\n  walk: {time.time()-t0:.1f}s, n_orbits ok = {n_orbits_ok:,}, "
        f"total transitions = {total_transitions:,}")

    # ---------------------------------------------------------------
    # Step 2: Build sub-stochastic kernel K with row normalization
    # ---------------------------------------------------------------
    visits = K_total_counts.sum(axis=1)  # outgoing-counts per state
    inflows = K_total_counts.sum(axis=0)  # incoming-counts per state

    n_states_with_visits = int((visits > 0).sum())
    n_states_total_used = int(((visits > 0) | (inflows > 0)).sum())
    log(f"\nStep 2: row-normalize K (sub-stochastic, m=1 absorbs)")
    log(f"  states with outgoing visits: {n_states_with_visits} / {N_STATES}")
    log(f"  states with any flow:        {n_states_total_used} / {N_STATES}")

    # Row-normalize: K[s, :] /= visits[s] (set zero rows to zero)
    K = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    nz = visits > 0
    K[nz, :] = K_total_counts[nz, :] / visits[nz][:, None]

    # Row sums (should be 1 for non-absorbing rows that completed all paths,
    # < 1 if some orbits ended at m=1 before next Syracuse step)
    row_sums = K.sum(axis=1)
    log(f"  row-sum stats over states with visits>0:")
    log(f"    min={row_sums[nz].min():.6f}, max={row_sums[nz].max():.6f}, "
         f"mean={row_sums[nz].mean():.6f}")

    # Mass leakage (= survival rate per Syracuse step)
    weighted_survival = ((row_sums[nz] * visits[nz]).sum() /
                         visits[nz].sum())
    log(f"  visit-weighted mean row-sum (survival rate per step) = {weighted_survival:.6f}")
    log(f"  (Result 50 reported empirical asymptotic survival ~ 0.94)")

    # ---------------------------------------------------------------
    # Step 3: Restrict to states with sufficient mass; compute leading
    # right eigenpair of K_sub
    # ---------------------------------------------------------------
    log("\nStep 3: leading right eigenpair of K (Perron-Frobenius)")

    # Filter: keep states with at least min_visits outgoing OR incoming
    min_visits = 50
    keep_mask = (visits >= min_visits) | (inflows >= min_visits)
    keep_idx = np.where(keep_mask)[0]
    n_kept = len(keep_idx)
    log(f"  keeping {n_kept} states with visits or inflows >= {min_visits}")

    K_sub = K[np.ix_(keep_idx, keep_idx)]
    K_sub_sparse = sp.csr_matrix(K_sub)

    # Leading right eigenvector: K v = lambda v
    # K is row-stochastic (~), so leading right eigvec is uniform (lambda~1)
    # — but we want the QSD interpretation: leading LEFT eigenvector pi K = lambda pi
    # which is the stationary distribution of mass on states.
    #
    # For sub-stochastic K, leading left eigenvector with eigenvalue lambda < 1
    # is the QSD (quasi-stationary distribution).
    #
    # Computing leading LEFT eigenpair = leading right of K^T
    log("  computing leading right eigenpair of K^T (= left of K = QSD)")
    try:
        vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=3,
                               which='LM', maxiter=5000, tol=1e-10)
    except spla.ArpackNoConvergence as e:
        log(f"  ARPACK no convergence: {e}; retrying with k=2")
        vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=2,
                               which='LM', maxiter=10000, tol=1e-8)

    # Sort by |lambda| descending
    order = np.argsort(-np.abs(vals))
    vals = vals[order]
    vecs = vecs[:, order]
    lam_PF = vals[0].real
    v_PF_sub = vecs[:, 0].real
    if v_PF_sub.sum() < 0:
        v_PF_sub = -v_PF_sub
    v_PF_sub = v_PF_sub / v_PF_sub.sum()

    log(f"  Top eigenvalues: " + ", ".join(f"{v:.6f}" for v in vals[:3]))
    log(f"  lambda_PF = {lam_PF:.6f}")
    log(f"  spectral gap = {1 - abs(vals[1])/abs(vals[0]):.4f}")

    # Embed back to full state space
    v_PF = np.zeros(N_STATES, dtype=np.float64)
    v_PF[keep_idx] = v_PF_sub

    # ---------------------------------------------------------------
    # Step 4: Marginalize v_PF to residue space
    # ---------------------------------------------------------------
    log("\nStep 4: marginalize v_PF over log-size bins to residue")

    rho_pred = np.zeros(16, dtype=np.float64)
    for r_idx in range(16):
        rho_pred[r_idx] = v_PF[r_idx * B:(r_idx + 1) * B].sum()
    rho_pred = rho_pred / rho_pred.sum()

    # D_predicted(r) = rho_pred(r) / pi_32(r)
    D_pred = {}
    for i, r in enumerate(odd_r32):
        D_pred[r] = rho_pred[i] / pi_32[r]

    # ---------------------------------------------------------------
    # Step 5: Compare D_predicted vs D_empirical
    # ---------------------------------------------------------------
    log("\nStep 5: D_predicted vs D_avg")
    log(f"\n  {'r':>3}  {'pi_32':>8}  {'D_avg':>8}  {'D_pred':>8}  {'diff':>8}")
    total_dev = 0.0
    diffs = []
    for r in odd_r32:
        d_pred = D_pred[r]
        d_emp = D_avg[r]
        diff = d_pred - d_emp
        total_dev += abs(diff)
        diffs.append(diff)
        log(f"  {r:>3}  {pi_32[r]:>8.5f}  {d_emp:>8.4f}  {d_pred:>8.4f}  "
            f"{diff:>+8.4f}")
    mad = total_dev / 16
    log(f"\n  Total |D_pred - D_avg| = {total_dev:.4f}")
    log(f"  MAD per residue = {mad:.4f}")

    D_avg_arr = np.array([D_avg[r] for r in odd_r32])
    D_pred_arr = np.array([D_pred[r] for r in odd_r32])
    pearson = float(np.corrcoef(D_avg_arr, D_pred_arr)[0, 1])
    spearman_r = D_avg_arr.argsort().argsort()
    spearman_p = D_pred_arr.argsort().argsort()
    spearman = float(np.corrcoef(spearman_r, spearman_p)[0, 1])
    log(f"  Pearson rho = {pearson:.4f}")
    log(f"  Spearman rho = {spearman:.4f}")

    # ---------------------------------------------------------------
    # Step 6: Examine joint structure
    # ---------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 6: joint structure of v_PF(r, b)")
    log("=" * 78)

    v_PF_2d = v_PF.reshape(R, B)  # (16, 64)

    # Marginal over r and over b
    pr = v_PF_2d.sum(axis=1)
    pr = pr / pr.sum()
    pb = v_PF_2d.sum(axis=0)
    pb_norm = pb / pb.sum() if pb.sum() > 0 else pb

    log("\n  Marginal P(b) (log-size distribution under QSD):")
    log(f"  {'b':>3}  {'P(b)':>8}  hist")
    for b_idx in range(B):
        if pb_norm[b_idx] > 1e-4:
            bar_len = int(pb_norm[b_idx] * 200)
            log(f"  {b_idx:>3}  {pb_norm[b_idx]:>8.5f}  "
                f"{'#' * min(bar_len, 60)}")

    # Test factorization: is v_PF(r, b) ~ f(r) * g(b)?
    f_r = pr
    g_b = pb_norm
    v_pred_factor = f_r[:, None] * g_b[None, :]
    factor_resid = v_PF_2d - v_pred_factor
    rms_factor = np.sqrt((factor_resid**2).sum() / (v_PF_2d**2).sum())
    log(f"\n  Test factorization v_PF(r,b) ?= f(r) * g(b):")
    log(f"  RMS residual relative to ||v_PF|| = {rms_factor:.4f}")
    log(f"  ({'separable' if rms_factor < 0.05 else 'NOT separable'} at 5% threshold)")

    # P(b | r): conditional log-size distribution per residue
    log("\n  P(b | r=21) and P(b | r=5) — peak log-size bins:")
    for r in [5, 21]:
        r_idx = (r - 1) // 2
        row = v_PF_2d[r_idx, :]
        if row.sum() > 0:
            row_norm = row / row.sum()
            top5 = np.argsort(-row_norm)[:5]
            log(f"  r={r}: top bins = " +
                ", ".join(f"b={b_}({row_norm[b_]:.3f})" for b_ in top5))

    # ---------------------------------------------------------------
    # Step 7: Connection to band-conditional structure
    # ---------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 7: sigma-band <-> log-size bin correspondence")
    log("=" * 78)

    # Compute sigma-band membership for orbits
    sigma_resid = sigma - K_H * log_n
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']

    # log_n-bin per orbit (start size only — orbit visits many bins during run)
    log_n_bin = np.minimum((log_n / LOG_2).astype(int), B - 1)
    log("\n  Orbit-start log_2(n) distribution by band:")
    log(f"  {'band':>6}  " + "  ".join(f"b={b}" for b in [25, 27, 29, 31]))
    for b_id in range(5):
        mask = band == b_id
        if mask.sum() == 0:
            continue
        starts_b = log_n_bin[mask]
        line = f"  {band_names[b_id]:>6}"
        for tb in [25, 27, 29, 31]:
            line += f"   {(starts_b == tb).mean():.3f}"
        log(line)
    log("\n  (Note: orbit visits MANY size bins during evolution; start-size")
    log("   binning here is only diagnostic. v_PF lives on visit-weighted bins.)")

    # ---------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Saving outputs")
    log("=" * 78)

    # Save sparse kernel
    K_sparse = sp.csr_matrix(K_sub)
    sp.save_npz(OUT / "size_stratified_kernel.npz", K_sparse)
    np.save(OUT / "size_stratified_keep_idx.npy", keep_idx)
    log("  [wrote] size_stratified_kernel.npz, size_stratified_keep_idx.npy")

    # Save eigenvector on (r, b)
    with open(OUT / "size_stratified_eigvec.csv", "w") as f:
        f.write("r,b,v_PF,P_r,P_b_given_r\n")
        for i, r in enumerate(odd_r32):
            r_idx = i
            row = v_PF_2d[r_idx, :]
            row_sum = row.sum()
            for b_idx in range(B):
                p_b_g_r = (row[b_idx] / row_sum) if row_sum > 0 else 0.0
                f.write(f"{r},{b_idx},{v_PF_2d[r_idx, b_idx]:.8e},"
                        f"{f_r[r_idx]:.6f},{p_b_g_r:.6f}\n")
    log("  [wrote] size_stratified_eigvec.csv")

    # Save residue marginal
    with open(OUT / "size_stratified_residue_marginal.csv", "w") as f:
        f.write("r,pi,D_avg,D_pred,diff,abs_diff\n")
        for r in odd_r32:
            d_pred = D_pred[r]
            diff = d_pred - D_avg[r]
            f.write(f"{r},{pi_32[r]:.6f},{D_avg[r]:.6f},{d_pred:.6f},"
                    f"{diff:+.6f},{abs(diff):.6f}\n")
    log("  [wrote] size_stratified_residue_marginal.csv")

    # ---------------------------------------------------------------
    # Final verdict
    # ---------------------------------------------------------------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    log(f"  Total |D_pred - D_avg| = {total_dev:.4f}")
    log(f"  Pearson rho = {pearson:.4f}")
    log(f"  Spearman rho = {spearman:.4f}")
    log(f"  lambda_PF = {lam_PF:.4f} (vs empirical 0.94 per-step survival)")
    log("")
    log(f"  Reference points:")
    log(f"    Cylinder QSD (R51, all 25 cells):   total dev 5.39 - 7.71")
    log(f"    Inverse-tree (R52, all 4 schemes):  total dev 13.2 - 27")
    log(f"    Renewal kernel exit (R53):          total dev 5.04, rho 0.58")
    log(f"    Trivial null (D_full this dataset): total dev 4.72")
    log(f"")
    if total_dev < 4.5 and pearson > 0.7:
        log(f"  Outcome: (alpha) — size-stratified Markov reproduces D_avg.")
    elif total_dev < 5.0 or pearson > 0.5:
        log(f"  Outcome: (beta) — partial; captures shape, missing pieces.")
    else:
        log(f"  Outcome: (gamma) — size-stratified misses despite tracking value info.")

    (OUT / "size_stratified_markov_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] size_stratified_markov_log.txt")


if __name__ == "__main__":
    main()
