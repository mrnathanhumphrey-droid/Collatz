"""
Result 60 v2 outcome-delta verification: validate Pearson 0.99 at B=150,
log_base=1.159 (or similar fine binning that spans 2^32 exactly).

Concern: is the 0.99 Pearson genuine structural identification, or
overfitting to kernel-estimation noise? Apply rigorous validation:
  1. Train-test split at B=150
  2. Sample-size stability at B=150
  3. Spectral gap analysis (very tight gap at fine binning may be artifact)
  4. Markov assumption (does it improve further?)
  5. Sweep B more carefully: where does Pearson saturate? What about 200, 250?
     (Try smaller min_visits threshold for finer bins.)
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
R_RES = 16

results_log = []
def log(s):
    print(s, flush=True)
    results_log.append(s)


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
            while mm & 1 == 0: mm >>= 1
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], 128)
    n = 32
    A = [[P[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n-1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n-1] = Fraction(1)
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


def evaluate_kernel(K_counts, B, min_visits, pi_32, D_avg, return_full=False):
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
    if n_kept < 50: return None
    K_sub = K[np.ix_(keep_idx, keep_idx)]
    K_sub_sparse = sp.csr_matrix(K_sub)
    try:
        vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=3,
                               which='LM', maxiter=20000, tol=1e-10)
    except Exception:
        try:
            vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=2,
                                   which='LM', maxiter=30000, tol=1e-7)
        except Exception:
            return None
    order = np.argsort(-np.abs(vals))
    vals = vals[order]; vecs = vecs[:, order]
    lam_PF = float(vals[0].real)
    lam_2 = float(vals[1].real) if len(vals) > 1 else 0.0
    v_sub = vecs[:, 0].real
    if v_sub.sum() < 0: v_sub = -v_sub
    if v_sub.sum() == 0: return None
    v_sub = v_sub / v_sub.sum()
    v_PF = np.zeros(N_STATES); v_PF[keep_idx] = v_sub

    rho_pred = np.zeros(R_RES)
    for r_idx in range(R_RES):
        rho_pred[r_idx] = v_PF[r_idx * B:(r_idx + 1) * B].sum()
    if rho_pred.sum() == 0: return None
    rho_pred /= rho_pred.sum()

    D_pred = np.array([rho_pred[i] / pi_32[r] for i, r in enumerate(odd_r32)])
    D_emp = np.array([D_avg[r] for r in odd_r32])
    diff = D_pred - D_emp
    total_dev = float(np.abs(diff).sum())
    pearson = float(np.corrcoef(D_emp, D_pred)[0, 1])

    res = dict(n_kept=n_kept, lam_PF=lam_PF, lam_2=lam_2,
               spectral_gap=1 - abs(lam_2)/abs(lam_PF) if abs(lam_PF) > 0 else 0.0,
               total_dev=total_dev, pearson=pearson)
    if return_full:
        res['v_PF'] = v_PF; res['D_pred'] = D_pred; res['D_emp'] = D_emp
    return res


def walk(N, n_orbits, seed, B, log_base):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    return walk_count_logbase(starts, 600, B, log_base, 12).sum(axis=0)


def main():
    pi_32 = chang_pi_32()
    df = pl.read_csv(OUT / "qsd_late_t_avg.csv")
    D_avg = {row['r']: row['D_avg'] for row in df.iter_rows(named=True)}

    log("=" * 80)
    log("OUTCOME (delta) VERIFICATION: Pearson 0.99 at B=150, log_base=2^(32/B)")
    log("=" * 80)
    log("\n  At B=109/log_base=1.226: scaling sweep showed Pearson 0.9829")
    log("  At B=150/log_base=1.159: scaling sweep showed Pearson 0.9888")
    log("  Need to verify: train-test, sample-size, Markov, spectral gap.")

    rows = []

    # ==========================
    # Reference: B=150, log_base = 2^(32/150)
    # ==========================
    log("\n=== Reference build at B=150, log_base=2^(32/150)=1.1594 ===\n")
    B_ref = 150
    log_base_ref = 2.0 ** (32.0 / B_ref)
    log(f"  log_base = {log_base_ref:.6f}")
    log(f"  bin range: 1.5^0=1 to {log_base_ref**B_ref:.2e}  (covers N=2^32 = {2**32:.2e})")

    K_total = np.zeros((R_RES * B_ref, R_RES * B_ref), dtype=np.int64)
    t0 = time.time()
    for seed in [42, 137, 271]:
        K_total += walk(1 << 32, 500_000, seed, B_ref, log_base_ref)
    log(f"  walk: {time.time()-t0:.1f}s, total transitions = {K_total.sum():,}")

    res = evaluate_kernel(K_total, B_ref, 50, pi_32, D_avg, return_full=True)
    log(f"\n  n_kept: {res['n_kept']} / {R_RES * B_ref}")
    log(f"  lambda_PF = {res['lam_PF']:.6f}, lambda_2 = {res['lam_2']:.6f}")
    log(f"  spectral gap = {res['spectral_gap']:.6f}")
    log(f"  total_dev = {res['total_dev']:.4f}, Pearson = {res['pearson']:.4f}")
    rows.append({'config': 'reference_B150', 'B': B_ref, 'log_base': log_base_ref,
                'n_kept': res['n_kept'], 'total_dev': res['total_dev'],
                'pearson': res['pearson'], 'lam_PF': res['lam_PF'],
                'spectral_gap': res['spectral_gap']})

    # Per-residue table
    log(f"\n  Per-residue:")
    log(f"  {'r':>3}  {'D_avg':>8}  {'D_pred':>8}  {'diff':>9}")
    odd_r32 = list(range(1, 32, 2))
    for i, r in enumerate(odd_r32):
        log(f"  {r:>3}  {res['D_emp'][i]:>8.4f}  {res['D_pred'][i]:>8.4f}  {res['D_pred'][i]-res['D_emp'][i]:>+9.4f}")

    ref_pearson = res['pearson']
    ref_total_dev = res['total_dev']

    # ==========================
    # Concern 1: Sample-size stability at B=150
    # ==========================
    log(f"\n=== Concern 1: Sample-size stability at B={B_ref}, log_base={log_base_ref:.4f} ===\n")
    log(f"  {'N':>10}  {'n_kept':>7}  {'total_dev':>10}  {'Pearson':>8}  {'lam_PF':>8}")
    for log2N in [28, 30, 32, 34]:
        N = 1 << log2N
        K_t = np.zeros((R_RES * B_ref, R_RES * B_ref), dtype=np.int64)
        for seed in [42, 137, 271]:
            K_t += walk(N, 500_000, seed, B_ref, log_base_ref)
        r = evaluate_kernel(K_t, B_ref, 50, pi_32, D_avg)
        if r is None:
            log(f"  N=2^{log2N}: eval failed")
            continue
        log(f"  2^{log2N:>3}  {r['n_kept']:>7}  {r['total_dev']:>10.4f}  {r['pearson']:>8.4f}  {r['lam_PF']:>8.4f}")
        rows.append({'config': f'sample_size_N=2^{log2N}', 'B': B_ref, 'log_base': log_base_ref,
                    'n_kept': r['n_kept'], 'total_dev': r['total_dev'],
                    'pearson': r['pearson'], 'lam_PF': r['lam_PF'],
                    'spectral_gap': r['spectral_gap']})

    # ==========================
    # Concern 2: Train-test split at B=150
    # ==========================
    log(f"\n=== Concern 2: Train-test split at B={B_ref} ===\n")
    log(f"  Train: starts in [3, N/4]; Test: starts in [N/4, N/2-1]")
    N = 1 << 32
    rng_tr = np.random.default_rng(11111)
    rng_te = np.random.default_rng(99999)
    starts_tr = 2 * rng_tr.integers(1, N // 4 - 1, size=500_000, dtype=np.int64) + 1
    starts_te = 2 * rng_te.integers(N // 4, (N - 1) // 2, size=500_000, dtype=np.int64) + 1
    K_tr = walk_count_logbase(starts_tr, 600, B_ref, log_base_ref, 12).sum(axis=0)
    K_te = walk_count_logbase(starts_te, 600, B_ref, log_base_ref, 12).sum(axis=0)

    res_tr = evaluate_kernel(K_tr, B_ref, 50, pi_32, D_avg, return_full=True)
    log(f"  Train K v_PF marginal vs global D_avg:")
    log(f"    total_dev = {res_tr['total_dev']:.4f}, Pearson = {res_tr['pearson']:.4f}")

    visits_te = K_te.sum(axis=1)
    rho_te = np.zeros(R_RES)
    for r_idx in range(R_RES):
        rho_te[r_idx] = visits_te[r_idx * B_ref:(r_idx + 1) * B_ref].sum()
    rho_te = rho_te / rho_te.sum() if rho_te.sum() > 0 else rho_te
    D_test = np.array([rho_te[i] / pi_32[r] for i, r in enumerate(odd_r32)])

    diff_tt = res_tr['D_pred'] - D_test
    pearson_tt = float(np.corrcoef(D_test, res_tr['D_pred'])[0, 1])
    log(f"  Train K v_PF vs Test-derived D_test:")
    log(f"    total_dev = {float(np.abs(diff_tt).sum()):.4f}, Pearson = {pearson_tt:.4f}")
    pearson_t_a = float(np.corrcoef(res_tr['D_emp'], D_test)[0, 1])
    log(f"    (Sanity: D_test vs D_avg Pearson={pearson_t_a:.4f})")

    rows.append({'config': 'train_vs_global_D_avg', 'B': B_ref, 'log_base': log_base_ref,
                'n_kept': res_tr['n_kept'], 'total_dev': res_tr['total_dev'],
                'pearson': res_tr['pearson'], 'lam_PF': res_tr['lam_PF'],
                'spectral_gap': float('nan')})
    rows.append({'config': 'train_vs_test_D', 'B': B_ref, 'log_base': log_base_ref,
                'n_kept': res_tr['n_kept'], 'total_dev': float(np.abs(diff_tt).sum()),
                'pearson': pearson_tt, 'lam_PF': res_tr['lam_PF'],
                'spectral_gap': float('nan')})

    # ==========================
    # Concern 3: Spectral gap analysis (very tight gap = overfitting?)
    # ==========================
    log(f"\n=== Concern 3: Spectral structure at B={B_ref} ===\n")
    log(f"  lambda_PF = {res['lam_PF']:.6f}")
    log(f"  lambda_2  = {res['lam_2']:.6f}")
    log(f"  spectral gap = 1 - |lam_2/lam_PF| = {res['spectral_gap']:.6f}")
    if res['spectral_gap'] < 0.01:
        log(f"  WARNING: spectral gap < 0.01 — leading eigenvector identification may be unstable")
    elif res['spectral_gap'] < 0.05:
        log(f"  NOTE: spectral gap < 0.05 — moderate, ARPACK should still resolve cleanly")
    else:
        log(f"  spectral gap >= 0.05 — well-separated leading eigenvector")

    # Run a smaller-tolerance ARPACK on the same kernel to cross-check
    visits = K_total.sum(axis=1); inflows = K_total.sum(axis=0)
    keep = (visits >= 50) | (inflows >= 50)
    keep_idx_full = np.where(keep)[0]
    K_full = np.zeros_like(K_total, dtype=np.float64)
    nz = visits > 0
    K_full[nz, :] = K_total[nz, :] / visits[nz][:, None]
    K_sub_arr = K_full[np.ix_(keep_idx_full, keep_idx_full)]

    # Cross-check with stricter tolerance
    try:
        vals_check, _ = spla.eigs(sp.csr_matrix(K_sub_arr.T), k=5,
                                   which='LM', maxiter=50000, tol=1e-12)
        sorted_vals = sorted(vals_check, key=lambda v: -abs(v))
        log(f"  Top 5 eigenvalues (stricter tol=1e-12):")
        for v in sorted_vals[:5]:
            log(f"    |lambda| = {abs(v):.6f}  (real={v.real:.6f}, imag={v.imag:.4e})")
    except Exception as e:
        log(f"  Stricter check failed: {e}")

    # ==========================
    # Concern 4: B sweep more fine, with min_visits adjusted
    # ==========================
    log(f"\n=== Concern 4: Broader B sweep (lowered min_visits for finer bins) ===\n")
    log(f"  {'B':>4}  {'log_base':>9}  {'thr':>4}  {'n_kept':>7}  {'total_dev':>10}  {'Pearson':>8}  {'lam_PF':>8}  {'gap':>7}")

    for B_t in [50, 75, 100, 125, 150, 175, 200, 250, 300]:
        log_base_t = 2.0 ** (32.0 / B_t)
        # Adjust threshold: at finer bins, fewer visits per bin, lower threshold
        # ~equally allocate visit budget across states
        avg_visits_per_state = 112_000_000 / (R_RES * B_t)
        thr = max(20, int(avg_visits_per_state * 0.005))  # 0.5% of avg
        K_t = np.zeros((R_RES * B_t, R_RES * B_t), dtype=np.int64)
        t1 = time.time()
        for seed in [42, 137, 271]:
            K_t += walk(1 << 32, 500_000, seed, B_t, log_base_t)
        r = evaluate_kernel(K_t, B_t, thr, pi_32, D_avg)
        elapsed = time.time() - t1
        if r is None:
            log(f"  {B_t:>4}  {log_base_t:>9.4f}  {thr:>4}  -- eval failed (--)")
            continue
        log(f"  {B_t:>4}  {log_base_t:>9.4f}  {thr:>4}  {r['n_kept']:>7}  {r['total_dev']:>10.4f}  {r['pearson']:>8.4f}  {r['lam_PF']:>8.4f}  {r['spectral_gap']:>7.4f}  ({elapsed:.0f}s)")
        rows.append({'config': f'B_sweep_B={B_t}', 'B': B_t, 'log_base': log_base_t,
                    'n_kept': r['n_kept'], 'total_dev': r['total_dev'],
                    'pearson': r['pearson'], 'lam_PF': r['lam_PF'],
                    'spectral_gap': r['spectral_gap']})

    pl.DataFrame(rows).write_csv(EXP_OUT / "result60_v2_outcome_delta.csv")
    log(f"\n  [save] result60_v2_outcome_delta.csv")

    # ==========================
    # Verdict
    # ==========================
    log(f"\n=== VERDICT ===\n")
    if ref_pearson >= 0.95 and pearson_tt >= 0.85:
        log(f"  Outcome (delta) CONFIRMED: Pearson {ref_pearson:.4f} at B={B_ref} survives train-test {pearson_tt:.4f}")
    elif ref_pearson >= 0.95 and pearson_tt < 0.85:
        log(f"  PARTIAL: Pearson {ref_pearson:.4f} but train-test drops to {pearson_tt:.4f} — overfitting")
    elif ref_pearson >= 0.91:
        log(f"  Outcome (alpha) at finer-than-1.5: Pearson {ref_pearson:.4f}")
    else:
        log(f"  Below 0.91 at B={B_ref}: investigation needed")

    (EXP_OUT / "result60_v2_outcome_delta_log.txt").write_text("\n".join(results_log), encoding="utf-8")


if __name__ == "__main__":
    main()
