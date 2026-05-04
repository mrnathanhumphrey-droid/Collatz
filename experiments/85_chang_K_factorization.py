"""
Operator factorization: Chang's transfer operator P vs trajectory measure
kernel K (R60 v2, B=109, log_base=1.5).

Tests:
  Hypothesis A: K projects through to P via residue marginalization
  Hypothesis B: K factors as K_residue x K_size (separable)
  Hypothesis C: no clean operator relation

Steps:
  1. Build Chang's P at mod 64 (32x32 row-stochastic), project to mod 32 (16x16)
  2. Load K at B=109/log_base=1.5 from result60_v2_kernel.npz
  3. K_residue := v_PF-weighted projection over (b, b') -> 16x16 row-stochastic
  4. Element-wise comparison of K_residue vs P_chang_mod32 (Frobenius, max gap)
  5. Stationary eigvec comparison: rho_K vs pi_chang_mod32
  6. SVD of K's joint (r,b) -> (r',b') structure: separable rank?
  7. Survivor-conditioning test: build K_unconditioned (uniform start), project,
     compare to P_chang
  8. Character decomposition (DFT on residues): identify character classes

Outputs:
  experiments_output/chang_K_factorization.csv
  experiments_output/chang_K_kernel_diff.csv
  experiments_output/K_svd_spectrum.csv
  experiments_output/chang_K_factorization_log.txt
"""
import sys
import io
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import polars as pl
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import njit, prange

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path("C:/Collatz")
EXP_OUT = Path("C:/Collatz/experiments_output")

R_RES = 16
MAX_VAL = np.int64(2**62)

results_log = []
def log(s):
    print(s, flush=True)
    results_log.append(s)


# ============================================================
# Chang's P at mod 64 (Definition C.5, depth 13, 128 lifts/residue)
# ============================================================

def build_chang_P_mod64():
    """32x32 row-stochastic. odd_residues[i] = 2i+1, indexed 0..31."""
    odd_residues = list(range(1, 64, 2))
    idx = {r: i for i, r in enumerate(odd_residues)}
    N_LIFTS = 128
    P = np.zeros((32, 32), dtype=np.float64)
    for i, r in enumerate(odd_residues):
        for k in range(N_LIFTS):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm % 2 == 0: mm //= 2
            P[i, idx[mm % 64]] += 1
        P[i, :] /= N_LIFTS
    return P, odd_residues, idx


def project_mod64_to_mod32(P_mod64, weights=None):
    """Project 32x32 mod-64 kernel to 16x16 mod-32 kernel.

    Aggregating odd r mod 64 -> r mod 32 by pairs (r, r+32).
    Weighted projection if `weights` provided (e.g., stationary).
    """
    odd_r32 = list(range(1, 32, 2))
    if weights is None:
        weights = np.ones(32) / 32.0
    P_mod32 = np.zeros((16, 16), dtype=np.float64)
    for i32, r1 in enumerate(odd_r32):
        # source class = {r1, r1+32}
        sources = []
        for j64, r in enumerate(range(1, 64, 2)):
            if r % 32 == r1: sources.append(j64)
        w_src = weights[sources]
        w_total = w_src.sum()
        if w_total == 0: continue
        for j32, r2 in enumerate(odd_r32):
            targets = []
            for k64, r in enumerate(range(1, 64, 2)):
                if r % 32 == r2: targets.append(k64)
            for s in sources:
                for t in targets:
                    P_mod32[i32, j32] += weights[s] * P_mod64[s, t]
            # nothing — we'll renormalize below
        P_mod32[i32, :] /= w_total
    return P_mod32


def stationary_left(M):
    """Leading left eigvec of M, normalized to sum to 1."""
    eigvals, eigvecs = np.linalg.eig(M.T)
    order = np.argsort(-eigvals.real)
    lam = float(eigvals.real[order[0]])
    v = eigvecs[:, order[0]].real
    if v.sum() < 0: v = -v
    if v.sum() == 0: return None, None
    return v / v.sum(), lam


# ============================================================
# K from R60 v2 (B=109, log_base=1.5)
# ============================================================

@njit(parallel=True, cache=True)
def walk_K_logbase(starts, max_T, B, log_base, n_chunks):
    """K_counts on (r mod 32, log-bin) state space.
    Returns also row-marginals for projection weighting."""
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


@njit(parallel=True, cache=True)
def walk_K_mod64(starts, max_T, B, log_base, n_chunks):
    """Same but with r mod 64 (32 odd residues) for direct Chang comparison."""
    n = len(starts)
    R_64 = 32  # 32 odd residues mod 64
    n_states = R_64 * B
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
            r_idx = (np.int64(m & 63) - 1) >> 1   # 0..31 for odd r mod 64
            b = int(np.floor(np.log(np.float64(m)) / log_bs))
            if b >= B: b = B - 1
            if b < 0: b = 0
            s_curr = r_idx * B + b
            while m != 1 and T < max_T:
                if m > MAX_VAL // 3: break
                x = 3 * m + 1
                while (x & 1) == 0: x >>= 1
                r_idx_n = (np.int64(x & 63) - 1) >> 1
                bn = int(np.floor(np.log(np.float64(x)) / log_bs))
                if bn >= B: bn = B - 1
                if bn < 0: bn = 0
                s_next = r_idx_n * B + bn
                K_counts[chunk, s_curr, s_next] += 1
                T += 1; m = x
                s_curr = s_next
    return K_counts


def evaluate_residue_kernel(K_counts, B, R_size, min_visits=50):
    """Build sub-stochastic K, leading left eigvec, marginalize to residue."""
    N_STATES = R_size * B
    visits = K_counts.sum(axis=1)
    inflows = K_counts.sum(axis=0)
    nz = visits > 0
    K = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    K[nz, :] = K_counts[nz, :] / visits[nz][:, None]
    keep = (visits >= min_visits) | (inflows >= min_visits)
    keep_idx = np.where(keep)[0]
    n_kept = len(keep_idx)
    if n_kept < 50: return None
    K_sub = K[np.ix_(keep_idx, keep_idx)]
    K_sub_sparse = sp.csr_matrix(K_sub)
    try:
        vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=3,
                               which='LM', maxiter=10000, tol=1e-10)
    except Exception:
        try:
            vals, vecs = spla.eigs(K_sub_sparse.T.astype(np.float64), k=2,
                                   which='LM', maxiter=20000, tol=1e-7)
        except Exception:
            return None
    order = np.argsort(-np.abs(vals))
    vals = vals[order]; vecs = vecs[:, order]
    lam_PF = float(vals[0].real)
    v_sub = vecs[:, 0].real
    if v_sub.sum() < 0: v_sub = -v_sub
    if v_sub.sum() == 0: return None
    v_sub = v_sub / v_sub.sum()
    v_PF = np.zeros(N_STATES); v_PF[keep_idx] = v_sub

    rho_marg = np.zeros(R_size)
    for r_idx in range(R_size):
        rho_marg[r_idx] = v_PF[r_idx * B:(r_idx + 1) * B].sum()
    if rho_marg.sum() == 0: return None
    rho_marg /= rho_marg.sum()

    return dict(K_sub=K_sub, K_full=K, keep_idx=keep_idx, v_PF=v_PF,
                rho=rho_marg, lam_PF=lam_PF, visits=visits)


def project_K_to_residue(K_full, v_PF, B, R_size):
    """v_PF-weighted projection: K_residue[r, r'] s.t. rho K_residue = lam rho.

    K_residue[r, r'] = (1/rho(r)) * sum_{b, b'} v_PF(r, b) * K[(r,b), (r',b')]
    """
    K_res = np.zeros((R_size, R_size), dtype=np.float64)
    rho = np.zeros(R_size)
    for r_idx in range(R_size):
        rho[r_idx] = v_PF[r_idx * B:(r_idx + 1) * B].sum()
    for r_idx in range(R_size):
        if rho[r_idx] == 0: continue
        for r_idx_n in range(R_size):
            mass = 0.0
            for b in range(B):
                vb = v_PF[r_idx * B + b]
                if vb == 0: continue
                for bn in range(B):
                    mass += vb * K_full[r_idx * B + b, r_idx_n * B + bn]
            K_res[r_idx, r_idx_n] = mass / rho[r_idx]
    return K_res


def main():
    # ===========================================================
    # Step 1: Chang's P at mod 64 and projection to mod 32
    # ===========================================================
    log("=" * 80)
    log("STEP 1: Chang's transfer operator P (Def C.5, depth 13, mod 64)")
    log("=" * 80)
    P_chang_64, odd_r64, idx_64 = build_chang_P_mod64()
    log(f"\n  P_chang shape: {P_chang_64.shape}")
    log(f"  Row-stochastic: max|row_sum-1| = {abs(P_chang_64.sum(axis=1) - 1).max():.2e}")
    pi_chang_64, lam = stationary_left(P_chang_64)
    log(f"  Leading eigval: {lam:.10f} (should be 1)")
    log(f"  pi_chang_64 sum: {pi_chang_64.sum():.6f}")

    # Projection to mod 32 (uniform-weighted, then stationary-weighted)
    log("\n  Project P_chang from mod 64 to mod 32:")
    P_chang_32_uniform = project_mod64_to_mod32(P_chang_64, weights=None)
    log(f"    Uniform-weighted projection row sums: max|sum-1| = {abs(P_chang_32_uniform.sum(axis=1) - 1).max():.2e}")
    pi_chang_32_unif, lam_u = stationary_left(P_chang_32_uniform)
    log(f"    Leading eigval (uniform proj): {lam_u:.10f}")

    P_chang_32_stat = project_mod64_to_mod32(P_chang_64, weights=pi_chang_64)
    log(f"    Stationary-weighted projection row sums: max|sum-1| = {abs(P_chang_32_stat.sum(axis=1) - 1).max():.2e}")
    pi_chang_32_stat, lam_s = stationary_left(P_chang_32_stat)
    log(f"    Leading eigval (stat proj): {lam_s:.10f}")

    # Direct: aggregate pi_chang_64 to mod 32
    pi_chang_32_direct = np.zeros(16)
    odd_r32 = list(range(1, 32, 2))
    for i32, r in enumerate(odd_r32):
        for j64 in range(32):
            if odd_r64[j64] % 32 == r:
                pi_chang_32_direct[i32] += pi_chang_64[j64]
    log(f"\n  pi_chang_32 (direct aggregation of mod-64 stationary):")
    log(f"  {'r':>3}  {'pi_64_aggr':>10}  {'pi_unif_proj':>13}  {'pi_stat_proj':>13}")
    for i, r in enumerate(odd_r32):
        log(f"  {r:>3}  {pi_chang_32_direct[i]:>10.6f}  {pi_chang_32_unif[i]:>13.6f}  {pi_chang_32_stat[i]:>13.6f}")

    # ===========================================================
    # Step 2-3: Load D_avg and our R60 v2 kernel; build K_residue
    # ===========================================================
    log("\n" + "=" * 80)
    log("STEP 2-3: K from R60 v2 (B=109, log_base=1.5), build K_residue")
    log("=" * 80)

    # Walk fresh K (using same seeds as R60 v2 reference)
    log_base = 1.5
    B_K = 109
    log(f"  Walking 1.5M orbits at N=2^32, B={B_K}, log_base={log_base}")
    t0 = time.time()
    K_total = np.zeros((R_RES * B_K, R_RES * B_K), dtype=np.int64)
    for seed in [42, 137, 271]:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, ((1 << 32) - 1) // 2, size=500_000, dtype=np.int64) + 1
        K_seed = walk_K_logbase(starts, 600, B_K, log_base, 12).sum(axis=0)
        K_total += K_seed
    log(f"  Walk: {time.time()-t0:.1f}s, transitions = {K_total.sum():,}")

    res_K = evaluate_residue_kernel(K_total, B_K, R_RES, min_visits=50)
    log(f"\n  K (R60 v2): n_kept={len(res_K['keep_idx'])}, lam_PF={res_K['lam_PF']:.6f}")
    log(f"  rho (residue marginal of v_PF, normalized to 1):")
    log(f"  {'r':>3}  {'rho':>9}  {'pi_chang_32':>12}  {'rho/pi':>8}")
    for i, r in enumerate(odd_r32):
        ratio = res_K['rho'][i] / pi_chang_32_direct[i]
        log(f"  {r:>3}  {res_K['rho'][i]:>9.6f}  {pi_chang_32_direct[i]:>12.6f}  {ratio:>8.4f}")

    # Build K_residue via v_PF-weighted projection
    log("\n  Building K_residue (16x16) via v_PF-weighted projection over bins...")
    t0 = time.time()
    K_residue = project_K_to_residue(res_K['K_full'], res_K['v_PF'], B_K, R_RES)
    log(f"  K_residue shape: {K_residue.shape}, build time: {time.time()-t0:.1f}s")
    log(f"  K_residue row sums: min={K_residue.sum(axis=1).min():.4f}, max={K_residue.sum(axis=1).max():.4f}")

    # Verify: rho K_residue = lam_PF rho?
    rho_K_residue = res_K['rho'] @ K_residue
    rho_K_residue_normalized = rho_K_residue / rho_K_residue.sum() if rho_K_residue.sum() > 0 else rho_K_residue
    log(f"\n  Self-consistency: rho @ K_residue ?= lam_PF * rho")
    log(f"  ||rho @ K_residue - lam_PF * rho|| = {np.linalg.norm(rho_K_residue - res_K['lam_PF'] * res_K['rho']):.6f}")
    log(f"  Pearson(rho, rho @ K_residue) = {np.corrcoef(res_K['rho'], rho_K_residue_normalized)[0,1]:.6f}")

    # ===========================================================
    # Step 4: Element-wise comparison K_residue vs P_chang_32
    # ===========================================================
    log("\n" + "=" * 80)
    log("STEP 4: K_residue vs P_chang_32 element-wise")
    log("=" * 80)

    diff_unif = K_residue - P_chang_32_uniform
    diff_stat = K_residue - P_chang_32_stat
    log(f"\n  K_residue vs uniform-projected P_chang:")
    log(f"    Frobenius ||K - P||_F = {np.linalg.norm(diff_unif):.4f}")
    log(f"    Max element gap = {np.abs(diff_unif).max():.4f}")
    log(f"    Mean abs element gap = {np.abs(diff_unif).mean():.4f}")
    log(f"\n  K_residue vs stationary-projected P_chang:")
    log(f"    Frobenius ||K - P||_F = {np.linalg.norm(diff_stat):.4f}")
    log(f"    Max element gap = {np.abs(diff_stat).max():.4f}")
    log(f"    Mean abs element gap = {np.abs(diff_stat).mean():.4f}")

    # Per-row Pearson
    pearson_per_row = []
    for i in range(16):
        if K_residue[i].sum() == 0 or P_chang_32_stat[i].sum() == 0: continue
        try:
            p = float(np.corrcoef(K_residue[i], P_chang_32_stat[i])[0, 1])
            pearson_per_row.append(p)
        except Exception:
            pass
    log(f"\n  Per-row Pearson(K_residue[i,:], P_chang[i,:]):")
    log(f"    mean={np.mean(pearson_per_row):.4f}, min={min(pearson_per_row):.4f}, max={max(pearson_per_row):.4f}")

    # Stationary comparison
    log(f"\n  Compare stationary distributions:")
    log(f"    pi_chang_32_direct (Chang stationary, aggregated)")
    log(f"    rho_K (residue marginal of K's v_PF)")
    log(f"  Pearson(pi_chang_32, rho_K) = {np.corrcoef(pi_chang_32_direct, res_K['rho'])[0,1]:.4f}")
    log(f"  ||pi_chang - rho_K||_1 = {np.abs(pi_chang_32_direct - res_K['rho']).sum():.4f}")

    rows_kernel = []
    for i in range(16):
        for j in range(16):
            rows_kernel.append({'r_from': odd_r32[i], 'r_to': odd_r32[j],
                              'K_residue': K_residue[i, j],
                              'P_chang_uniform': P_chang_32_uniform[i, j],
                              'P_chang_stat': P_chang_32_stat[i, j],
                              'diff_K_minus_Pstat': K_residue[i, j] - P_chang_32_stat[i, j]})
    pl.DataFrame(rows_kernel).write_csv(EXP_OUT / "chang_K_kernel_diff.csv")
    log(f"  [save] chang_K_kernel_diff.csv")

    # ===========================================================
    # Step 5: SVD / separability of K
    # ===========================================================
    log("\n" + "=" * 80)
    log("STEP 5: SVD of K's joint structure (Hypothesis B test)")
    log("=" * 80)

    # Reshape K_full as (R, B, R, B) -> matricize as (R*B, R*B) and compare to
    # tensor-product approximation K_r ⊗ K_b
    K_4d = res_K['K_full'].reshape(R_RES, B_K, R_RES, B_K)
    log(f"  K reshaped as (R, B, R, B) = {K_4d.shape}")

    # Build K_r (residue-only) as marginalization
    K_r_visit = K_4d.sum(axis=(1, 3))  # 16x16, sum over (b, b')
    visits_r = res_K['visits'].reshape(R_RES, B_K).sum(axis=1)  # row visits per r
    nz_r = visits_r > 0
    K_r = np.zeros((R_RES, R_RES))
    K_r[nz_r] = K_r_visit[nz_r] / visits_r[nz_r, None] * B_K  # rough normalization

    # Build K_b (size-only) as marginalization (sum over r, r')
    K_b_visit = K_4d.sum(axis=(0, 2))  # B x B
    visits_b = res_K['visits'].reshape(R_RES, B_K).sum(axis=0)
    nz_b = visits_b > 0
    K_b = np.zeros((B_K, B_K))
    K_b[nz_b] = K_b_visit[nz_b] / visits_b[nz_b, None] * R_RES

    # Tensor product approximation
    K_tensor = np.einsum('ij,kl->ikjl', K_r, K_b).reshape(R_RES * B_K, R_RES * B_K)
    K_tensor_sub = K_tensor[np.ix_(res_K['keep_idx'], res_K['keep_idx'])]
    K_actual_sub = res_K['K_full'][np.ix_(res_K['keep_idx'], res_K['keep_idx'])]

    # Normalize tensor approximation row-wise to match scale
    rs = K_tensor_sub.sum(axis=1)
    nz = rs > 0
    K_tensor_norm = np.zeros_like(K_tensor_sub)
    K_tensor_norm[nz] = K_tensor_sub[nz] / rs[nz, None] * K_actual_sub[nz].sum(axis=1)[:, None]

    diff_sep = K_actual_sub - K_tensor_norm
    rms_sep = np.sqrt((diff_sep**2).mean())
    log(f"\n  K vs K_r ⊗ K_b separable approx:")
    log(f"    RMS difference = {rms_sep:.4f}")
    log(f"    Frobenius ratio ||K - K_tensor||_F / ||K||_F = {np.linalg.norm(diff_sep)/np.linalg.norm(K_actual_sub):.4f}")

    # SVD spectrum: rank of K_actual_sub vs separable
    # Reshape K_actual_sub back conceptually: it's already a matrix
    # SVD of K_actual_sub itself gives rank info
    log("\n  SVD spectrum of K_sub (top 10 singular values):")
    K_sparse = sp.csr_matrix(K_actual_sub)
    try:
        u, s, vt = spla.svds(K_sparse.astype(np.float64), k=10)
        s_sorted = sorted(s, reverse=True)
        for i, sv in enumerate(s_sorted[:10]):
            log(f"    sigma_{i+1} = {sv:.6f}")
        log(f"  Total Frobenius: ||K||_F = {np.linalg.norm(K_actual_sub):.4f}")
        log(f"  Rank-1 approx captures: {s_sorted[0]/np.linalg.norm(K_actual_sub):.4f}")
    except Exception as e:
        log(f"  SVD failed: {e}")

    rows_svd = []
    try:
        for i, sv in enumerate(s_sorted[:10]):
            rows_svd.append({'rank': i+1, 'singular_value': float(sv)})
        pl.DataFrame(rows_svd).write_csv(EXP_OUT / "K_svd_spectrum.csv")
    except Exception:
        pass

    # ===========================================================
    # Step 6: Survivor-conditioning test
    # ===========================================================
    log("\n" + "=" * 80)
    log("STEP 6: Survivor-conditioning test")
    log("=" * 80)

    # Build K' using uniform residue counts (i.e., P_chang's lift convention).
    # Here we approximate "unconditioned" by walking just the FIRST step from
    # uniform starts: lift uniform per residue mod 64 with 128 lifts each,
    # see where they land. This is exactly Chang's P!
    log("\n  Note: 'Unconditioned' = first-step image under uniform lift.")
    log("        That IS Chang's P by construction. So:")
    log("          K_unconditioned_residue ≡ P_chang")
    log("        And: K_residue (survivor-conditioned, multi-step composite) ≠ P_chang")
    log("        is expected by definition.\n")

    log("  Quantifying the deviation:")
    log(f"    Mean abs element gap K_residue vs P_chang_stat = {np.abs(diff_stat).mean():.4f}")
    log(f"    Stationary Pearson(rho_K, pi_chang_32) = {np.corrcoef(pi_chang_32_direct, res_K['rho'])[0,1]:.4f}")

    # Quantify per-residue: which residues' marginals diverge most?
    log("\n  Per-residue stationary divergence:")
    log(f"  {'r':>3}  {'pi_chang_32':>12}  {'rho_K':>10}  {'rho/pi':>8}  {'sign':>6}")
    for i, r in enumerate(odd_r32):
        ratio = res_K['rho'][i] / pi_chang_32_direct[i]
        sign = '+' if ratio > 1.05 else '-' if ratio < 0.95 else '='
        log(f"  {r:>3}  {pi_chang_32_direct[i]:>12.6f}  {res_K['rho'][i]:>10.6f}  {ratio:>8.4f}  {sign:>6}")

    # ===========================================================
    # Step 7: Character decomposition via DFT on residues
    # ===========================================================
    log("\n" + "=" * 80)
    log("STEP 7: Character decomposition (DFT on Z/16Z odd residues)")
    log("=" * 80)

    # Stationary distributions as functions on odd residues mod 32 (16 values).
    # Compute their DFT on Z/16Z (since there are 16 odd residues, indexed 0..15
    # by odd_r32[i] = 2i+1).
    pi_chang_arr = pi_chang_32_direct.copy()
    rho_K_arr = res_K['rho'].copy()

    # Center
    pi_centered = pi_chang_arr - pi_chang_arr.mean()
    rho_centered = rho_K_arr - rho_K_arr.mean()

    F_chang = np.fft.fft(pi_centered)
    F_K = np.fft.fft(rho_centered)

    log("\n  DFT of stationary distributions (centered):")
    log(f"  {'k':>3}  {'|F_chang(k)|':>12}  {'|F_K(k)|':>10}  {'arg_chang':>9}  {'arg_K':>7}")
    for k in range(16):
        log(f"  {k:>3}  {abs(F_chang[k]):>12.6f}  {abs(F_K[k]):>10.6f}  {np.angle(F_chang[k]):>9.4f}  {np.angle(F_K[k]):>7.4f}")

    # Identify dominant Fourier modes
    chang_modes = sorted(range(16), key=lambda k: -abs(F_chang[k]))[:5]
    K_modes = sorted(range(16), key=lambda k: -abs(F_K[k]))[:5]
    log(f"\n  Chang's top 5 modes by |F|: {chang_modes}")
    log(f"  K's top 5 modes by |F|:     {K_modes}")

    # ===========================================================
    # Step 8: Synthesis verdict
    # ===========================================================
    log("\n" + "=" * 80)
    log("STEP 8: Synthesis verdict")
    log("=" * 80)

    pearson_stat = float(np.corrcoef(pi_chang_32_direct, res_K['rho'])[0,1])
    rms_kernel = float(np.linalg.norm(diff_stat)) / float(np.linalg.norm(P_chang_32_stat))
    log(f"\n  Stationary Pearson: {pearson_stat:.4f}")
    log(f"  Kernel RMS rel diff: {rms_kernel:.4f}")
    log(f"  Separability rel diff: {np.linalg.norm(diff_sep)/np.linalg.norm(K_actual_sub):.4f}")

    if pearson_stat > 0.95 and rms_kernel < 0.10:
        log("\n  VERDICT (alpha): explicit factorization holds")
    elif pearson_stat > 0.7 or rms_kernel < 0.30:
        log("\n  VERDICT (beta): partial relation; Chang's P approximates K_residue")
    else:
        log("\n  VERDICT (gamma): no clean operator factorization; complementary frameworks")

    # Save summary
    rows_summary = [
        {'metric': 'stationary Pearson(pi_chang, rho_K)', 'value': pearson_stat},
        {'metric': 'kernel Frobenius rel diff (K_res vs P_chang_stat)', 'value': rms_kernel},
        {'metric': 'kernel max element gap', 'value': float(np.abs(diff_stat).max())},
        {'metric': 'kernel mean element gap', 'value': float(np.abs(diff_stat).mean())},
        {'metric': 'separability RMS', 'value': float(rms_sep)},
        {'metric': 'separability rel Frobenius', 'value': float(np.linalg.norm(diff_sep)/np.linalg.norm(K_actual_sub))},
        {'metric': 'lambda_K', 'value': res_K['lam_PF']},
        {'metric': 'lambda_chang', 'value': lam},
    ]
    pl.DataFrame(rows_summary).write_csv(EXP_OUT / "chang_K_factorization.csv")
    log(f"\n  [save] chang_K_factorization.csv")

    (EXP_OUT / "chang_K_factorization_log.txt").write_text("\n".join(results_log), encoding="utf-8")


if __name__ == "__main__":
    main()
