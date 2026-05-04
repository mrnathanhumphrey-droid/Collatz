"""
Validation protocol for Result 60 (size-stratified Markov framework).

Stress-tests Result 60's central claim — D_avg is residue marginal of leading
left eigenvector v_PF of the empirical 1024-state kernel K — across:

  Concern 1: Sample size (N), bin count (B), min_visits threshold
  Concern 2: Train/test split holdout
  Concern 3: Bin boundary shift; alternate log base
  Concern 4: Markov-property diagnostic (conditional MI)
  Concern 5: Apples-to-apples framework comparison
  Concern 6: lambda_PF auto-match check (scrambled-kernel null)
  Concern 7: 79% factorization residual vs random-vector null

Outputs
  experiments_output/82_validation_results.csv
  experiments_output/82_markov_diagnostics.csv
  experiments_output/82_validation_log.txt
"""
import sys
import io
import time
from pathlib import Path

import numpy as np
import polars as pl
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import njit, prange

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path("C:/Collatz")
EXP_OUT = Path("C:/Collatz/experiments_output")

LOG_2 = np.log(2.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)

R_RES = 16  # 16 odd residues mod 32

results_log = []
def log(s):
    print(s, flush=True)
    results_log.append(s)


# ============================================================
# Walker — captures ALL transitions with full state history (m, r, b)
# Builds K_counts[from_state, to_state] AND records transitions in a flat array
# for higher-order MI analysis.
# ============================================================

@njit(parallel=True, cache=True)
def walk_count(starts, max_T, B, n_chunks):
    """
    State: idx = ((r-1)//2) * B + b,  r = m mod 32 (odd),  b = floor(log2 m).
    Returns K_counts (n_chunks, R*B, R*B) int64 + per-orbit completion indicator.
    """
    n = len(starts)
    n_states = R_RES * B
    chunk_size = (n + n_chunks - 1) // n_chunks
    K_counts = np.zeros((n_chunks, n_states, n_states), dtype=np.int64)
    completed = np.zeros(n, dtype=np.int8)

    for chunk in prange(n_chunks):
        i_lo = chunk * chunk_size
        i_hi = min((chunk + 1) * chunk_size, n)
        for i in range(i_lo, i_hi):
            m = np.int64(starts[i])
            T = 0
            while (m & 1) == 0 and m > 1:
                m >>= 1
            if m == 1:
                completed[i] = 1
                continue

            r_idx = (np.int64(m & 31) - 1) >> 1
            mm = m; b = 0
            while mm > 1:
                mm >>= 1; b += 1
            if b >= B: b = B - 1
            s_curr = r_idx * B + b

            failed = False
            while m != 1 and T < max_T:
                if m > MAX_VAL // 3:
                    failed = True; break
                x = 3 * m + 1
                v = 0
                while (x & 1) == 0:
                    x >>= 1; v += 1
                r_idx_n = (np.int64(x & 31) - 1) >> 1
                xx = x; bn = 0
                while xx > 1:
                    xx >>= 1; bn += 1
                if bn >= B: bn = B - 1
                s_next = r_idx_n * B + bn
                K_counts[chunk, s_curr, s_next] += 1
                T += 1; m = x
                s_curr = s_next
            if not failed and m == 1:
                completed[i] = 1
    return K_counts, completed


@njit(parallel=True, cache=True)
def walk_record_flat(starts, max_T, max_steps_per_orbit, B):
    """Record (s_curr, s_next) pairs for every transition. For higher-order MI."""
    n = len(starts)
    flat_curr = np.full(n * max_steps_per_orbit, -1, dtype=np.int32)
    flat_next = np.full(n * max_steps_per_orbit, -1, dtype=np.int32)
    flat_orbit = np.full(n * max_steps_per_orbit, -1, dtype=np.int32)
    n_records = np.zeros(n, dtype=np.int32)

    for i in prange(n):
        m = np.int64(starts[i])
        while (m & 1) == 0 and m > 1:
            m >>= 1
        if m == 1: continue

        r_idx = (np.int64(m & 31) - 1) >> 1
        mm = m; b = 0
        while mm > 1: mm >>= 1; b += 1
        if b >= B: b = B - 1
        s_curr = r_idx * B + b

        T = 0
        offset = i * max_steps_per_orbit
        while m != 1 and T < max_T and T < max_steps_per_orbit:
            if m > MAX_VAL // 3: break
            x = 3 * m + 1
            v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1
            r_idx_n = (np.int64(x & 31) - 1) >> 1
            xx = x; bn = 0
            while xx > 1: xx >>= 1; bn += 1
            if bn >= B: bn = B - 1
            s_next = r_idx_n * B + bn
            flat_curr[offset + T] = s_curr
            flat_next[offset + T] = s_next
            flat_orbit[offset + T] = i
            T += 1; m = x
            s_curr = s_next
        n_records[i] = T

    # Truncate to actual records
    total = int(n_records.sum())
    out_curr = np.empty(total, dtype=np.int32)
    out_next = np.empty(total, dtype=np.int32)
    out_orb = np.empty(total, dtype=np.int32)
    pos = 0
    for i in range(n):
        cnt = n_records[i]
        offset = i * max_steps_per_orbit
        out_curr[pos:pos+cnt] = flat_curr[offset:offset+cnt]
        out_next[pos:pos+cnt] = flat_next[offset:offset+cnt]
        out_orb[pos:pos+cnt] = flat_orbit[offset:offset+cnt]
        pos += cnt
    return out_curr, out_next, out_orb


# ============================================================
# Helpers: Chang's pi_32, D_avg, framework eval
# ============================================================

def chang_pi_32():
    """Reproduce Chang's stationary pi on residues mod 32."""
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
    return {r32: float(b[idx[r32]] + b[idx[r32 + 32]]) for r32 in range(1, 32, 2)}


def load_D_avg():
    df = pl.read_csv(OUT / "qsd_late_t_avg.csv")
    return {row['r']: row['D_avg'] for row in df.iter_rows(named=True)}


def evaluate_kernel(K_counts, B, min_visits, pi_32, D_avg, label=""):
    """Build K, restrict, compute v_PF, marginalize, compare to D_avg."""
    odd_r32 = list(range(1, 32, 2))
    N_STATES = R_RES * B
    visits = K_counts.sum(axis=1)
    inflows = K_counts.sum(axis=0)
    nz = visits > 0
    K = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    K[nz, :] = K_counts[nz, :] / visits[nz][:, None]

    keep_mask = (visits >= min_visits) | (inflows >= min_visits)
    keep_idx = np.where(keep_mask)[0]
    n_kept = len(keep_idx)
    if n_kept < 50:
        return None
    K_sub = K[np.ix_(keep_idx, keep_idx)]
    K_sub_sparse = sp.csr_matrix(K_sub)

    try:
        vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=3,
                               which='LM', maxiter=5000, tol=1e-10)
    except spla.ArpackNoConvergence:
        try:
            vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=2,
                                   which='LM', maxiter=10000, tol=1e-8)
        except Exception:
            return None
    order = np.argsort(-np.abs(vals))
    vals = vals[order]; vecs = vecs[:, order]
    lam_PF = float(vals[0].real)
    v_sub = vecs[:, 0].real
    if v_sub.sum() < 0:
        v_sub = -v_sub
    if v_sub.sum() == 0:
        return None
    v_sub = v_sub / v_sub.sum()
    v_PF = np.zeros(N_STATES); v_PF[keep_idx] = v_sub

    rho_pred = np.zeros(R_RES)
    for r_idx in range(R_RES):
        rho_pred[r_idx] = v_PF[r_idx * B:(r_idx + 1) * B].sum()
    if rho_pred.sum() == 0:
        return None
    rho_pred /= rho_pred.sum()

    D_pred = np.array([rho_pred[i] / pi_32[r] for i, r in enumerate(odd_r32)])
    D_emp = np.array([D_avg[r] for r in odd_r32])

    diff = D_pred - D_emp
    total_dev = float(np.abs(diff).sum())
    pearson = float(np.corrcoef(D_emp, D_pred)[0, 1])
    spearman = float(np.corrcoef(D_emp.argsort().argsort(),
                                  D_pred.argsort().argsort())[0, 1])

    return dict(
        label=label, n_kept=n_kept, lam_PF=lam_PF,
        total_dev=total_dev, pearson=pearson, spearman=spearman,
        v_PF=v_PF, D_pred=D_pred, D_emp=D_emp,
    )


def walk_orbits(N, n_orbits, seed, max_T, B, n_chunks=12):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    K_counts, completed = walk_count(starts, max_T, B, n_chunks)
    return K_counts.sum(axis=0), completed.sum(), starts


def main():
    pi_32 = chang_pi_32()
    D_avg = load_D_avg()
    odd_r32 = list(range(1, 32, 2))
    rows = []

    log("=" * 80)
    log("RESULT 60 VALIDATION PROTOCOL")
    log("=" * 80)
    log(f"Reference (R60 baseline): N=2^32, 1.5M orbits, B=64, min_visits=50")
    log(f"  -> total_dev=3.40, Pearson=0.80, lambda_PF=0.9566")

    # ===========================================================
    # CONCERN 1: Sample-size dependence
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 1: Sample-size, B, threshold dependence")
    log("=" * 80)

    log("\n  1a. Vary N at fixed B=64, min_visits=50, 500K orbits/seed * 3 seeds:")
    log(f"  {'N':>10}  {'n_kept':>7}  {'total_dev':>10}  {'Pearson':>8}  {'lam_PF':>8}")
    for log2N in [28, 30, 32, 34]:
        N = 1 << log2N
        K_total = np.zeros((R_RES * 64, R_RES * 64), dtype=np.int64)
        t0 = time.time()
        for seed in [42, 137, 271]:
            K_seed, _, _ = walk_orbits(N, 500_000, seed, 600, 64)
            K_total += K_seed
        elapsed = time.time() - t0
        res = evaluate_kernel(K_total, 64, 50, pi_32, D_avg, label=f"N=2^{log2N}")
        log(f"  2^{log2N:>3}  {res['n_kept']:>7}  {res['total_dev']:>10.4f}  {res['pearson']:>8.4f}  {res['lam_PF']:>8.4f}  ({elapsed:.0f}s)")
        rows.append({'concern': '1a_sample_size', 'param': f'N=2^{log2N}',
                    'n_kept': res['n_kept'], 'total_dev': res['total_dev'],
                    'pearson': res['pearson'], 'lam_PF': res['lam_PF']})

    # 1b: Vary B at fixed N=2^32, min_visits=50
    log("\n  1b. Vary B at fixed N=2^32:")
    log(f"  {'B':>4}  {'n_kept':>7}  {'total_dev':>10}  {'Pearson':>8}  {'lam_PF':>8}")
    for B_val in [32, 64, 128]:
        K_total = np.zeros((R_RES * B_val, R_RES * B_val), dtype=np.int64)
        for seed in [42, 137, 271]:
            K_seed, _, _ = walk_orbits(1 << 32, 500_000, seed, 600, B_val)
            K_total += K_seed
        res = evaluate_kernel(K_total, B_val, 50, pi_32, D_avg, label=f"B={B_val}")
        log(f"  {B_val:>4}  {res['n_kept']:>7}  {res['total_dev']:>10.4f}  {res['pearson']:>8.4f}  {res['lam_PF']:>8.4f}")
        rows.append({'concern': '1b_B', 'param': f'B={B_val}',
                    'n_kept': res['n_kept'], 'total_dev': res['total_dev'],
                    'pearson': res['pearson'], 'lam_PF': res['lam_PF']})

    # 1c: Vary min_visits at N=2^32, B=64
    log("\n  1c. Vary min_visits threshold at N=2^32, B=64 (reuse from 1a):")
    K_total = np.zeros((R_RES * 64, R_RES * 64), dtype=np.int64)
    for seed in [42, 137, 271]:
        K_seed, _, _ = walk_orbits(1 << 32, 500_000, seed, 600, 64)
        K_total += K_seed
    log(f"  {'thresh':>7}  {'n_kept':>7}  {'total_dev':>10}  {'Pearson':>8}  {'lam_PF':>8}")
    for thr in [25, 50, 100, 200]:
        res = evaluate_kernel(K_total, 64, thr, pi_32, D_avg, label=f"thresh={thr}")
        log(f"  {thr:>7}  {res['n_kept']:>7}  {res['total_dev']:>10.4f}  {res['pearson']:>8.4f}  {res['lam_PF']:>8.4f}")
        rows.append({'concern': '1c_threshold', 'param': f'thresh={thr}',
                    'n_kept': res['n_kept'], 'total_dev': res['total_dev'],
                    'pearson': res['pearson'], 'lam_PF': res['lam_PF']})

    # ===========================================================
    # CONCERN 2: Train-test split
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 2: Train-test split (build K on half, eval on other half's D)")
    log("=" * 80)
    log("\n  Plan: walk orbits with same seed but split by start range.")
    log("        Train half: starts in lower N/2; Test half: upper N/2.")

    log2N = 32; N = 1 << log2N; B_val = 64
    n_per_half = 500_000

    rng_train = np.random.default_rng(11111)
    rng_test = np.random.default_rng(99999)
    starts_train = 2 * rng_train.integers(1, N // 4 - 1, size=n_per_half, dtype=np.int64) + 1
    starts_test = 2 * rng_test.integers(N // 4, (N - 1) // 2, size=n_per_half, dtype=np.int64) + 1

    log(f"  Train: {n_per_half:,} orbits in [3, N/2-1], Test: {n_per_half:,} in [N/2+1, N-1]")

    # For test set, we need D_test (residue marginal of empirical visits)
    # Easiest: build a kernel from test, use its visit counts to compute
    # D_test_residue = (sum of all visits in residue r) / pi_32(r)
    K_train, _ = walk_count(starts_train, 600, B_val, 12)
    K_train = K_train.sum(axis=0)
    K_test, _ = walk_count(starts_test, 600, B_val, 12)
    K_test = K_test.sum(axis=0)

    # Train kernel -> v_PF, marginalize, compare to D_avg (which is independent global truth)
    res_train = evaluate_kernel(K_train, B_val, 50, pi_32, D_avg, label="train_kernel_vs_D_avg")
    log(f"\n  Train K -> v_PF marginal vs global D_avg:")
    log(f"    total_dev={res_train['total_dev']:.4f}  Pearson={res_train['pearson']:.4f}  lam_PF={res_train['lam_PF']:.4f}")
    rows.append({'concern': '2_train_test', 'param': 'train_v_PF_vs_global_D_avg',
                'n_kept': res_train['n_kept'], 'total_dev': res_train['total_dev'],
                'pearson': res_train['pearson'], 'lam_PF': res_train['lam_PF']})

    # Build D_test from TEST set's visits — independent measurement of residue marginal
    visits_test = K_test.sum(axis=1)  # per-state visits
    rho_test = np.zeros(R_RES)
    for r_idx in range(R_RES):
        rho_test[r_idx] = visits_test[r_idx * B_val:(r_idx + 1) * B_val].sum()
    rho_test = rho_test / rho_test.sum() if rho_test.sum() > 0 else rho_test
    D_test = np.array([rho_test[i] / pi_32[r] for i, r in enumerate(odd_r32)])

    # Train v_PF predicted residue marginal, compare to TEST-derived D
    diff_train_vs_test = res_train['D_pred'] - D_test
    total_dev_tt = float(np.abs(diff_train_vs_test).sum())
    pearson_tt = float(np.corrcoef(D_test, res_train['D_pred'])[0, 1])
    log(f"\n  Train K v_PF residue marginal vs Test-derived D_test:")
    log(f"    total_dev={total_dev_tt:.4f}  Pearson={pearson_tt:.4f}")
    rows.append({'concern': '2_train_test', 'param': 'train_v_PF_vs_test_D',
                'n_kept': res_train['n_kept'], 'total_dev': total_dev_tt,
                'pearson': pearson_tt, 'lam_PF': res_train['lam_PF']})

    # Sanity: D_test vs D_avg (test set's marginal vs reference)
    diff_test_vs_avg = D_test - res_train['D_emp']
    pearson_t_a = float(np.corrcoef(res_train['D_emp'], D_test)[0, 1])
    log(f"\n  Sanity: D_test (independent measurement) vs global D_avg:")
    log(f"    total_dev={float(np.abs(diff_test_vs_avg).sum()):.4f}  Pearson={pearson_t_a:.4f}")
    log(f"    (high correlation = D_avg is well-defined and reproduces across samples)")

    # ===========================================================
    # CONCERN 3: Bin boundary sensitivity
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 3: Bin boundary sensitivity")
    log("=" * 80)

    # 3a: Shifted bins via different log base / offset.
    # We test by remapping b post-hoc inside walker can't easily — just rerun
    # with a shifted size scaling by walking with a modified B definition.
    # Cleanest: write a small variant walker that computes b = floor(log_2(m) + 0.5)
    log("\n  3a. Shifted bins: b = floor(log_2 m + 0.5):")

    @njit(parallel=True, cache=True)
    def walk_count_shifted(starts, max_T, B, n_chunks):
        n = len(starts)
        n_states = R_RES * B
        chunk_size = (n + n_chunks - 1) // n_chunks
        K_counts = np.zeros((n_chunks, n_states, n_states), dtype=np.int64)
        for chunk in prange(n_chunks):
            i_lo = chunk * chunk_size
            i_hi = min((chunk + 1) * chunk_size, n)
            for i in range(i_lo, i_hi):
                m = np.int64(starts[i])
                T = 0
                while (m & 1) == 0 and m > 1: m >>= 1
                if m == 1: continue
                # b = floor(log2(m) + 0.5) = floor(log2(2^0.5 * m)) = floor(log2(m * sqrt2))
                # equivalently: m * sqrt(2) approximated via m << 1 sqrt approx
                # Easier: compute using float log
                r_idx = (np.int64(m & 31) - 1) >> 1
                b = int(np.floor(np.log2(np.float64(m)) + 0.5))
                if b >= B: b = B - 1
                if b < 0: b = 0
                s_curr = r_idx * B + b
                while m != 1 and T < max_T:
                    if m > MAX_VAL // 3: break
                    x = 3 * m + 1
                    while (x & 1) == 0: x >>= 1
                    r_idx_n = (np.int64(x & 31) - 1) >> 1
                    bn = int(np.floor(np.log2(np.float64(x)) + 0.5))
                    if bn >= B: bn = B - 1
                    if bn < 0: bn = 0
                    s_next = r_idx_n * B + bn
                    K_counts[chunk, s_curr, s_next] += 1
                    T += 1; m = x
                    s_curr = s_next
        return K_counts

    K_total = np.zeros((R_RES * 64, R_RES * 64), dtype=np.int64)
    for seed in [42, 137, 271]:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, ((1 << 32) - 1) // 2, size=500_000, dtype=np.int64) + 1
        K_seed = walk_count_shifted(starts, 600, 64, 12).sum(axis=0)
        K_total += K_seed
    res = evaluate_kernel(K_total, 64, 50, pi_32, D_avg, label="b_shifted_0.5")
    log(f"    total_dev={res['total_dev']:.4f}  Pearson={res['pearson']:.4f}  lam_PF={res['lam_PF']:.4f}")
    rows.append({'concern': '3a_bin_shift', 'param': 'b=floor(log2 m + 0.5)',
                'n_kept': res['n_kept'], 'total_dev': res['total_dev'],
                'pearson': res['pearson'], 'lam_PF': res['lam_PF']})

    # 3b: Different log base (log_3, log_e, log_1.5)
    log("\n  3b. Different log base (B chosen to keep state count similar):")

    @njit(parallel=True, cache=True)
    def walk_count_logbase(starts, max_T, B, log_base, n_chunks):
        n = len(starts)
        n_states = R_RES * B
        chunk_size = (n + n_chunks - 1) // n_chunks
        K_counts = np.zeros((n_chunks, n_states, n_states), dtype=np.int64)
        log_bs = np.log(log_base)
        for chunk in prange(n_chunks):
            i_lo = chunk * chunk_size
            i_hi = min((chunk + 1) * chunk_size, n)
            for i in range(i_lo, i_hi):
                m = np.int64(starts[i])
                T = 0
                while (m & 1) == 0 and m > 1: m >>= 1
                if m == 1: continue
                r_idx = (np.int64(m & 31) - 1) >> 1
                b = int(np.floor(np.log(np.float64(m)) / log_bs))
                if b >= B: b = B - 1
                if b < 0: b = 0
                s_curr = r_idx * B + b
                while m != 1 and T < max_T:
                    if m > MAX_VAL // 3: break
                    x = 3 * m + 1
                    while (x & 1) == 0: x >>= 1
                    r_idx_n = (np.int64(x & 31) - 1) >> 1
                    bn = int(np.floor(np.log(np.float64(x)) / log_bs))
                    if bn >= B: bn = B - 1
                    if bn < 0: bn = 0
                    s_next = r_idx_n * B + bn
                    K_counts[chunk, s_curr, s_next] += 1
                    T += 1; m = x
                    s_curr = s_next
        return K_counts

    log(f"  {'base':>5}  {'B_chosen':>9}  {'n_kept':>7}  {'total_dev':>10}  {'Pearson':>8}")
    # log base x with B such that B*log(x) ~= 64*log(2): B = 64*log2/log(x)
    for base in [3.0, np.e, 1.5]:
        B_val = max(8, int(np.round(64 * LOG_2 / np.log(base))))
        K_total = np.zeros((R_RES * B_val, R_RES * B_val), dtype=np.int64)
        for seed in [42, 137, 271]:
            rng = np.random.default_rng(seed)
            starts = 2 * rng.integers(1, ((1 << 32) - 1) // 2, size=500_000, dtype=np.int64) + 1
            K_seed = walk_count_logbase(starts, 600, B_val, base, 12).sum(axis=0)
            K_total += K_seed
        res = evaluate_kernel(K_total, B_val, 50, pi_32, D_avg, label=f"log_base={base:.2f}")
        log(f"  {base:>5.2f}  {B_val:>9}  {res['n_kept']:>7}  {res['total_dev']:>10.4f}  {res['pearson']:>8.4f}")
        rows.append({'concern': '3b_log_base', 'param': f'base={base:.2f},B={B_val}',
                    'n_kept': res['n_kept'], 'total_dev': res['total_dev'],
                    'pearson': res['pearson'], 'lam_PF': res['lam_PF']})

    # ===========================================================
    # CONCERN 4: Markov-property diagnostic
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 4: Markov assumption diagnostic")
    log("=" * 80)
    log("\n  Compute mutual info I(s_{t+1}; s_{t-1} | s_t) — should be ~0 for Markov.")
    log("  Compare to I(s_{t+1}; s_t) — main coupling.")

    # Walk smaller sample with full record at B=64
    rng = np.random.default_rng(42)
    starts_mi = 2 * rng.integers(1, ((1 << 32) - 1) // 2, size=50_000, dtype=np.int64) + 1
    flat_curr, flat_next, flat_orb = walk_record_flat(starts_mi, 600, 600, 64)

    # Build adjacent-pair state sequence for orbits with >= 3 records
    # For each orbit, if records are s0, s1, s2, ..., we have pairs (s_t -> s_{t+1})
    # in flat_curr/flat_next. To get triples, group by orbit.
    log(f"  Records: {len(flat_curr):,} from {len(starts_mi):,} orbits")

    # Triple (s_{t-1}, s_t, s_{t+1}) construction
    # In the flat array: for orbit i with k records, we have s_curr indexed by 0..k-1
    # and s_next indexed by 1..k. So s_{t-1} = flat_curr[j-1], s_t = flat_curr[j] = flat_next[j-1],
    # s_{t+1} = flat_next[j] for j = 1..k-1.
    # Easier: triples come from positions j and j+1 in same orbit.
    # Use flat_orb to detect orbit boundaries.
    s_prev = flat_curr[:-1]
    s_t = flat_curr[1:]
    s_next = flat_next[1:]
    same_orb = (flat_orb[:-1] == flat_orb[1:])
    s_prev = s_prev[same_orb]
    s_t = s_t[same_orb]
    s_next = s_next[same_orb]
    n_triples = len(s_prev)
    log(f"  Triples (s_{{t-1}}, s_t, s_{{t+1}}): {n_triples:,}")

    # Mutual information via entropy decomposition
    def entropy_int(x):
        _, counts = np.unique(x, return_counts=True)
        p = counts.astype(np.float64) / counts.sum()
        return -float((p * np.log2(p)).sum())

    def joint_entropy(x, y):
        # bijection: combine into single int (assumes x, y < 2^16)
        xy = x.astype(np.int64) << 16 | y.astype(np.int64)
        return entropy_int(xy)

    H_t = entropy_int(s_t)
    H_next = entropy_int(s_next)
    H_t_next = joint_entropy(s_t, s_next)
    I_t_next = H_t + H_next - H_t_next
    log(f"  H(s_t) = {H_t:.3f} bits, H(s_next) = {H_next:.3f}, I(s_t; s_next) = {I_t_next:.3f}")

    # I(s_prev; s_next | s_t) = H(s_prev | s_t) + H(s_next | s_t) - H(s_prev, s_next | s_t)
    # = H(s_prev, s_t) + H(s_next, s_t) - H(s_t) - H(s_prev, s_t, s_next)
    H_prev_t = joint_entropy(s_prev, s_t)
    H_next_t = joint_entropy(s_next, s_t)
    # Triple entropy: combine 3 ints — need 64-bit safety; states < 1024 so 10 bits each fits in 30 bits
    pst = (s_prev.astype(np.int64) << 20) | (s_t.astype(np.int64) << 10) | s_next.astype(np.int64)
    H_pst = entropy_int(pst)
    I_prev_next_given_t = H_prev_t + H_next_t - H_t - H_pst
    log(f"  I(s_prev; s_next | s_t) = {I_prev_next_given_t:.4f} bits")
    ratio = I_prev_next_given_t / I_t_next if I_t_next > 0 else float('nan')
    log(f"  Ratio I(prev;next|t) / I(t;next) = {ratio:.4f}")
    log(f"  Markov assumption: ratio should be << 1 (ideally < 0.05)")

    rows.append({'concern': '4_markov_MI', 'param': 'I(prev;next|t)',
                'n_kept': n_triples, 'total_dev': I_prev_next_given_t,
                'pearson': ratio, 'lam_PF': float('nan')})

    # ===========================================================
    # CONCERN 5: Apples-to-apples framework comparison
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 5: Apples-to-apples framework comparison")
    log("=" * 80)
    log("\n  All competitors at same N=2^32, same D_avg target, same evaluation metric.")
    log("  Pulling baselines from existing files for inverse-tree, renewal, R51, null.")

    # Trivial null = pi_32 distribution (uniform-on-residue-frequency = D = 1 everywhere)
    D_emp = np.array([D_avg[r] for r in odd_r32])
    null_pred = np.ones_like(D_emp)
    null_dev = float(np.abs(null_pred - D_emp).sum())
    log(f"  Trivial null (D_pred = 1 everywhere): total_dev = {null_dev:.4f}")
    rows.append({'concern': '5_compare', 'param': 'null_D=1', 'n_kept': 16,
                'total_dev': null_dev, 'pearson': 0.0, 'lam_PF': float('nan')})

    # Reload R60 baseline at N=2^32 for direct comparison
    K_total = np.zeros((R_RES * 64, R_RES * 64), dtype=np.int64)
    for seed in [42, 137, 271]:
        K_seed, _, _ = walk_orbits(1 << 32, 500_000, seed, 600, 64)
        K_total += K_seed
    res_r60 = evaluate_kernel(K_total, 64, 50, pi_32, D_avg, label="R60_baseline")
    log(f"  R60 baseline N=2^32: total_dev={res_r60['total_dev']:.4f}  Pearson={res_r60['pearson']:.4f}")

    # Inverse tree at SAME conditions: load from existing inverse_tree eigvec at mod 32
    # Result 23's v_max(odd) ratio is the prediction.
    df_inv = pl.read_csv(OUT / "inverse_tree" / "inverse_tree_eigvec_mod32.csv")
    eigvec_full = {row['residue_mod_32']: row['predicted_density'] for row in df_inv.iter_rows(named=True)}
    eigvec_odd = np.array([eigvec_full[r] for r in odd_r32])
    rho_inv = eigvec_odd / eigvec_odd.sum()
    D_inv = np.array([rho_inv[i] / pi_32[r] for i, r in enumerate(odd_r32)])
    inv_dev = float(np.abs(D_inv - D_emp).sum())
    inv_pearson = float(np.corrcoef(D_emp, D_inv)[0, 1])
    log(f"  Inverse tree (R23, no truncation): total_dev = {inv_dev:.4f}  Pearson={inv_pearson:.4f}")
    rows.append({'concern': '5_compare', 'param': 'inverse_tree_R23',
                'n_kept': 16, 'total_dev': inv_dev, 'pearson': inv_pearson,
                'lam_PF': float('nan')})

    # ===========================================================
    # CONCERN 6: lambda_PF auto-match (scrambled kernel)
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 6: lambda_PF auto-match check (scrambled kernel)")
    log("=" * 80)
    log("\n  Hypothesis: lambda_PF ~ 0.94 is automatic for any kernel with")
    log("  realistic per-state survival rates, NOT informative of correctness.")
    log("\n  Test: scramble the column structure of K (preserve row sums)")
    log("        and check if lambda_PF still ~ 0.94.")

    visits = K_total.sum(axis=1)
    inflows = K_total.sum(axis=0)
    keep = (visits >= 50) | (inflows >= 50)
    keep_idx = np.where(keep)[0]
    K = np.zeros_like(K_total, dtype=np.float64)
    nz = visits > 0
    K[nz, :] = K_total[nz, :] / visits[nz][:, None]
    K_sub = K[np.ix_(keep_idx, keep_idx)]
    row_sums = K_sub.sum(axis=1)

    # Real K eigenvalue
    vals, _ = spla.eigs(sp.csr_matrix(K_sub.T), k=3, which='LM', maxiter=5000, tol=1e-10)
    lam_real = float(np.sort(np.abs(vals))[-1])
    log(f"  Real K_sub lambda_PF = {lam_real:.4f}")

    # Scramble: in each row, permute the column entries (preserves row sums)
    rng = np.random.default_rng(11)
    K_scram = K_sub.copy()
    for i in range(K_scram.shape[0]):
        perm = rng.permutation(K_scram.shape[1])
        K_scram[i, :] = K_scram[i, perm]
    vals_s, _ = spla.eigs(sp.csr_matrix(K_scram.T), k=3, which='LM', maxiter=5000, tol=1e-10)
    lam_scram = float(np.sort(np.abs(vals_s))[-1])
    log(f"  Scrambled K (row entries permuted) lambda_PF = {lam_scram:.4f}")

    # Strict: replace columns with uniform redistribution preserving row sums
    K_uniform = np.zeros_like(K_sub)
    n_kept = K_sub.shape[0]
    for i in range(n_kept):
        K_uniform[i, :] = row_sums[i] / n_kept
    vals_u, _ = spla.eigs(sp.csr_matrix(K_uniform.T), k=3, which='LM', maxiter=5000, tol=1e-10)
    lam_unif = float(np.sort(np.abs(vals_u))[-1])
    log(f"  Uniform-column K (row sums preserved): lambda_PF = {lam_unif:.4f}")
    log(f"\n  If scrambled and uniform also give lambda ~ 0.94, the eigenvalue match")
    log(f"  is automatic and NOT informative evidence for the framework.")

    rows.append({'concern': '6_lambda_auto', 'param': 'real_K',
                'n_kept': len(keep_idx), 'total_dev': float('nan'), 'pearson': float('nan'),
                'lam_PF': lam_real})
    rows.append({'concern': '6_lambda_auto', 'param': 'scrambled_K',
                'n_kept': len(keep_idx), 'total_dev': float('nan'), 'pearson': float('nan'),
                'lam_PF': lam_scram})
    rows.append({'concern': '6_lambda_auto', 'param': 'uniform_col_K',
                'n_kept': len(keep_idx), 'total_dev': float('nan'), 'pearson': float('nan'),
                'lam_PF': lam_unif})

    # ===========================================================
    # CONCERN 7: 79% factorization residual null distribution
    # ===========================================================
    log("\n" + "=" * 80)
    log("CONCERN 7: Factorization residual relative to null distribution")
    log("=" * 80)

    v_PF_2d = res_r60['v_PF'].reshape(R_RES, 64)
    pr = v_PF_2d.sum(axis=1)
    pb = v_PF_2d.sum(axis=0)
    f_r = pr / pr.sum() if pr.sum() > 0 else pr
    g_b = pb / pb.sum() if pb.sum() > 0 else pb
    v_factored = (pr.sum() * f_r[:, None]) * g_b[None, :] / pr.sum()  # = pr.sum * f_r * g_b
    factor_resid = v_PF_2d - v_factored
    rms_real = float(np.sqrt((factor_resid**2).sum() / (v_PF_2d**2).sum()))
    log(f"  Real v_PF factorization residual: {rms_real:.4f} (= {rms_real*100:.1f}%)")

    # Null distribution: random vectors with same support
    support = (v_PF_2d > 0).astype(np.float64)
    n_support = int(support.sum())
    log(f"  v_PF support size: {n_support}/{R_RES*64} = {n_support/(R_RES*64)*100:.1f}%")

    rng = np.random.default_rng(7)
    null_resids = []
    for trial in range(200):
        v_random = rng.exponential(scale=1.0, size=(R_RES, 64)) * support
        if v_random.sum() == 0: continue
        v_random = v_random / v_random.sum()
        pr_r = v_random.sum(axis=1)
        pb_r = v_random.sum(axis=0)
        f_r_r = pr_r / pr_r.sum() if pr_r.sum() > 0 else pr_r
        g_b_r = pb_r / pb_r.sum() if pb_r.sum() > 0 else pb_r
        v_fact = pr_r.sum() * f_r_r[:, None] * g_b_r[None, :]
        resid = v_random - v_fact
        null_resids.append(np.sqrt((resid**2).sum() / (v_random**2).sum()))

    null_resids = np.array(null_resids)
    log(f"  Null (200 random vectors with same support, exponential weights):")
    log(f"    mean RMS residual = {null_resids.mean():.4f}")
    log(f"    std RMS residual  = {null_resids.std():.4f}")
    log(f"    range            = [{null_resids.min():.4f}, {null_resids.max():.4f}]")
    log(f"  Real value {rms_real:.4f} z-score vs null = {(rms_real - null_resids.mean()) / null_resids.std():.2f}")

    rows.append({'concern': '7_factor_resid', 'param': 'real_v_PF',
                'n_kept': n_support, 'total_dev': rms_real,
                'pearson': (rms_real - null_resids.mean()) / null_resids.std(),
                'lam_PF': float('nan')})
    rows.append({'concern': '7_factor_resid', 'param': 'null_mean',
                'n_kept': 200, 'total_dev': float(null_resids.mean()),
                'pearson': float(null_resids.std()), 'lam_PF': float('nan')})

    # ===========================================================
    # SAVE + Final summary
    # ===========================================================
    log("\n" + "=" * 80)
    log("SAVING")
    log("=" * 80)

    df_out = pl.DataFrame(rows)
    df_out.write_csv(EXP_OUT / "82_validation_results.csv")
    log(f"  [save] 82_validation_results.csv")

    (EXP_OUT / "82_validation_log.txt").write_text("\n".join(results_log), encoding="utf-8")
    log(f"  [save] 82_validation_log.txt")


if __name__ == "__main__":
    main()
