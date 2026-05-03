"""
Experiment 45 — Absorbing Markov chain analysis of the Syracuse map.

DELIVERABLE 1 (PRIMARY): derive P(j) and W_j for entry classes m_j = (4^j-1)/3
on the chain over odd m in [3, M] for M ∈ {10^5, 10^6}, via sparse linear
algebra on (I - Q)^{-1}. Match against empirical 50M-orbit targets at N=2^36:

  j  m_j     P(j) emp   W_j emp
  2  5       0.93787    +7.156
  4  85      0.024      -4.755
  5  341     0.038      +4.59
  ε_S asymptote = 1.374

Method:
  - Each odd m ∈ [3, M] is a state. Transition m → T(m) = (3m+1)/2^{v_2(3m+1)}.
  - Absorbing: m_j ∈ [5, M] for j = 2..max_j. One extra "escape" column for
    transient → successor > M.
  - The chain is DETERMINISTIC (each row of Q has at most 1 nonzero), so:
      B = (I-Q)^{-1} R   gives 0/1 absorption indicators.
      τ = (I-Q)^{-1} 1   gives steps to absorption.
    For deterministic chains, the conditional fundamental matrix N_j reduces
    to τ restricted to ancestors of j. Same answer; computed once globally.
  - Aggregate over starts in [3, M] under three priors.

Reference: Kemeny & Snell, *Finite Markov Chains* (1960), ch. 3.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG_4_3 = np.log(4) - np.log(3)


def entry_classes(M):
    """List of (j, m_j) with m_j = (4^j-1)/3 in [5, M]."""
    out = []
    j = 2
    while True:
        m = (4**j - 1) // 3
        if m > M:
            break
        out.append((j, m))
        j += 1
    return out


@njit(cache=True, parallel=True)
def syracuse_succ(starts, M):
    n = len(starts)
    succ = np.empty(n, dtype=np.int64)
    for i in prange(n):
        x = 3 * starts[i] + 1
        while x & 1 == 0:
            x >>= 1
        succ[i] = -1 if x > M else x
    return succ


def build_chain(M):
    odd_arr = np.arange(3, M + 1, 2, dtype=np.int64)

    ec = entry_classes(M)
    abs_m = np.array([m for _, m in ec], dtype=np.int64)
    abs_j_arr = np.array([j for j, _ in ec], dtype=np.int32)
    K_abs = len(ec)

    abs_map = -np.ones(M + 2, dtype=np.int64)
    abs_map[abs_m] = np.arange(K_abs)

    is_abs_in_odd = abs_map[odd_arr] >= 0
    transient_m = odd_arr[~is_abs_in_odd]
    n_T = len(transient_m)

    transient_map = -np.ones(M + 2, dtype=np.int64)
    transient_map[transient_m] = np.arange(n_T)

    print(f"  building Q on n_T={n_T:,}, K_abs={K_abs}")
    t0 = time.perf_counter()
    succ = syracuse_succ(transient_m, M)
    print(f"    successors in {time.perf_counter()-t0:.2f}s")

    s_clip = np.where(succ < 0, 0, succ)
    is_escape = succ == -1
    a_idx = abs_map[s_clip]
    is_abs_succ = (~is_escape) & (a_idx >= 0)
    is_trans_succ = (~is_escape) & (~is_abs_succ)

    rows_Q = np.where(is_trans_succ)[0]
    cols_Q = transient_map[succ[rows_Q]]
    rows_R_abs = np.where(is_abs_succ)[0]
    cols_R_abs = a_idx[rows_R_abs]
    rows_R_esc = np.where(is_escape)[0]
    cols_R_esc = np.full(len(rows_R_esc), K_abs, dtype=np.int64)

    rows_R = np.concatenate([rows_R_abs, rows_R_esc])
    cols_R = np.concatenate([cols_R_abs, cols_R_esc])

    Q = sp.csr_matrix((np.ones(len(rows_Q)), (rows_Q, cols_Q)),
                      shape=(n_T, n_T), dtype=np.float64)
    R = sp.csr_matrix((np.ones(len(rows_R)), (rows_R, cols_R)),
                      shape=(n_T, K_abs + 1), dtype=np.float64)

    return {
        'transient_m': transient_m, 'abs_m': abs_m, 'abs_j': abs_j_arr,
        'K_abs': K_abs, 'n_T': n_T, 'Q': Q, 'R': R,
    }


def solve_chain(chain):
    n_T, Q, R = chain['n_T'], chain['Q'], chain['R']
    I_minus_Q = (sp.eye(n_T, format='csc') - Q.tocsc()).tocsc()
    print(f"  splu factor (I-Q), nnz={I_minus_Q.nnz:,}")
    t0 = time.perf_counter()
    lu = spla.splu(I_minus_Q)
    print(f"    factor in {time.perf_counter()-t0:.2f}s")

    print(f"  solve B = LU^{{-1}} R_dense ({R.shape[1]} columns)")
    t0 = time.perf_counter()
    B = lu.solve(R.toarray())
    print(f"    B in {time.perf_counter()-t0:.2f}s")

    print(f"  solve τ = LU^{{-1}} 1")
    t0 = time.perf_counter()
    tau = lu.solve(np.ones(n_T))
    print(f"    τ in {time.perf_counter()-t0:.2f}s")

    return B, tau


def aggregate(chain, B, tau):
    transient_m = chain['transient_m']
    abs_m = chain['abs_m']
    K_abs = chain['K_abs']
    M = transient_m.max()

    log_m_T = np.log(transient_m.astype(np.float64))

    priors_raw = {
        'uniform_top_half': (transient_m >= M // 2).astype(np.float64),
        'uniform_all':      np.ones_like(transient_m, dtype=np.float64),
        'log_uniform':      1.0 / transient_m.astype(np.float64),
    }

    log_mj = np.log(abs_m.astype(np.float64))

    results = {}
    for name, p_raw in priors_raw.items():
        s = p_raw.sum()
        if s == 0:
            continue
        p = p_raw / s

        P_marg = p @ B
        P_j = P_marg[:K_abs]
        P_escape = P_marg[K_abs]

        sigma_w = p * (tau + 1)  # +1 for the m_j → 1 step
        sigma_w_sq = p * (tau + 1)**2
        log_w = p * log_m_T

        sigma_per_j = sigma_w @ B[:, :K_abs]
        sigma_sq_per_j = sigma_w_sq @ B[:, :K_abs]
        log_per_j = log_w @ B[:, :K_abs]

        denom = np.maximum(P_j, 1e-300)
        msS_j = sigma_per_j / denom
        msS2_j = sigma_sq_per_j / denom
        var_sS_j = msS2_j - msS_j**2
        sd_sS_j = np.sqrt(np.maximum(var_sS_j, 0))
        mlog_j = log_per_j / denom

        # Conditional Wald W_j
        W_j = msS_j - mlog_j / LOG_4_3 - 1 + log_mj / LOG_4_3

        # ε_S over terminating orbits
        non_esc = B[:, :K_abs].sum(axis=1)
        P_term = P_j.sum()
        if P_term > 0:
            sig_tot = (sigma_w * non_esc).sum() / P_term
            log_tot = (log_w * non_esc).sum() / P_term
            eps_S = sig_tot - log_tot / LOG_4_3
        else:
            eps_S = float('nan')

        results[name] = {
            'P_j': P_j, 'P_escape': P_escape,
            'msS_j': msS_j, 'msS2_j': msS2_j,
            'var_sS_j': var_sS_j, 'sd_sS_j': sd_sS_j,
            'mlog_j': mlog_j, 'W_j': W_j,
            'P_term': P_term, 'eps_S': eps_S,
        }

    return results


EMPIR = {
    2:  {'P': 0.93787, 'W': +7.156, 'm': 5},
    4:  {'P': 0.02400, 'W': -4.755, 'm': 85},
    5:  {'P': 0.03800, 'W': +4.590, 'm': 341},
    7:  {'P': 0.00008, 'W': -12.61, 'm': 5461},
    8:  {'P': 0.00049, 'W': +2.17,  'm': 21845},
    10: {'P': 0.00003, 'W': +1.79,  'm': 349525},
}


@njit(parallel=True, cache=True)
def empirical_per_j_walker(starts, max_value, max_steps):
    """Walk N orbits to completion, accumulate σ_S and σ_S^2 per entry class j.

    Returns (sum_sS_j, sum_sS2_j, cnt_j, sum_log_j, overflow) where j-index 0..12.
    """
    n = len(starts)
    n_chunks = 32
    chunk = (n + n_chunks - 1) // n_chunks
    max_j = 12
    sum_sS_j = np.zeros((n_chunks, max_j+1), dtype=np.float64)
    sum_sS2_j = np.zeros((n_chunks, max_j+1), dtype=np.float64)
    sum_log_j = np.zeros((n_chunks, max_j+1), dtype=np.float64)
    cnt_j = np.zeros((n_chunks, max_j+1), dtype=np.int64)
    overflow = np.zeros(n_chunks, dtype=np.int64)
    for c in prange(n_chunks):
        lo = c*chunk; hi = min(lo+chunk, n)
        for i in range(lo, hi):
            x = np.int64(starts[i])
            log_x = np.log(np.float64(x))
            sigma_S = 0
            steps = 0
            last_j = np.int64(-1)
            failed = False
            while x != 1 and steps < max_steps:
                if x % 2 == 1:
                    if x > max_value // 3:
                        failed = True; break
                    threex_p1 = 3*x + 1
                    if threex_p1 > 0 and (threex_p1 & (threex_p1-1)) == 0:
                        k = 0; tmp = threex_p1
                        while tmp > 1: tmp >>= 1; k += 1
                        if k % 2 == 0:
                            last_j = k // 2
                    x = threex_p1
                    sigma_S += 1
                else:
                    x = x // 2
                steps += 1
            if failed or x != 1:
                overflow[c] += 1
                continue
            if 0 <= last_j <= max_j:
                sum_sS_j[c, last_j] += sigma_S
                sum_sS2_j[c, last_j] += sigma_S*sigma_S
                sum_log_j[c, last_j] += log_x
                cnt_j[c, last_j] += 1
    return (sum_sS_j.sum(axis=0), sum_sS2_j.sum(axis=0),
            cnt_j.sum(axis=0), sum_log_j.sum(axis=0), overflow.sum())


def empirical_var_at_N(N_log2, n_starts):
    """Run walker at N=2^N_log2 with n_starts uniform odd starts. Return per-j Var, mean, SD."""
    N = 1 << N_log2
    rng = np.random.default_rng(42)
    starts = 2 * rng.integers(1, (N-1)//2, size=n_starts, dtype=np.int64) + 1
    print(f"  walking {n_starts:,} orbits at N=2^{N_log2}...")
    t0 = time.perf_counter()
    sum_sS, sum_sS2, cnt, sum_log, of = empirical_per_j_walker(starts, np.int64(2**62), 200_000)
    print(f"    walk done in {time.perf_counter()-t0:.1f}s, overflow={of}")
    out = {}
    for j in range(2, 13):
        if cnt[j] > 100:
            mean_sS = sum_sS[j] / cnt[j]
            mean_sS2 = sum_sS2[j] / cnt[j]
            var = mean_sS2 - mean_sS**2
            mean_log = sum_log[j] / cnt[j]
            out[j] = {
                'mean_sS': mean_sS,
                'var_sS': var,
                'sd_sS': np.sqrt(max(var, 0)),
                'mean_log': mean_log,
                'cnt': int(cnt[j]),
                'P': cnt[j] / cnt.sum(),
            }
    return out


def report(M, chain, results):
    abs_j = chain['abs_j']
    abs_m = chain['abs_m']

    print(f"\n=== Results at M = {M:,} ===")
    print(f"  Empirical asymptote at N=2^36: ε_S = 1.374")

    for name, r in results.items():
        # Normalize P(j) by P_term for empirical comparison (which has 0 escape)
        P_cond = r['P_j'] / max(r['P_term'], 1e-300)
        print(f"\n  Prior: {name}")
        print(f"    P_escape = {r['P_escape']:.5f}")
        print(f"    P_terminate = {r['P_term']:.5f}")
        print(f"    ε_S (chain) = {r['eps_S']:.4f}  (asymp emp: 1.374, gap {r['eps_S']-1.374:+.4f})")
        print(f"    {'j':>3} {'m_j':>9} {'P(j)|term':>10} {'P_emp':>9} "
              f"{'gap_P':>9} {'W_j':>9} {'W_emp':>9} {'gap_W':>9} "
              f"{'⟨σ_S|j⟩':>10} {'SD[σ_S|j]':>10}")
        for k, j in enumerate(abs_j):
            E = EMPIR.get(int(j), {})
            P_e = E.get('P', float('nan'))
            W_e = E.get('W', float('nan'))
            gap_W = r['W_j'][k] - W_e if not np.isnan(W_e) else float('nan')
            gap_P = P_cond[k] - P_e if not np.isnan(P_e) else float('nan')
            print(f"    {j:>3} {abs_m[k]:>9} {P_cond[k]:>10.5f} {P_e:>9.5f} "
                  f"{gap_P:>+9.5f} {r['W_j'][k]:>+9.4f} {W_e:>+9.3f} {gap_W:>+9.4f} "
                  f"{r['msS_j'][k]:>10.3f} {r['sd_sS_j'][k]:>10.3f}")


def to_csv_rows(M, chain, results):
    rows = []
    abs_j = chain['abs_j']
    abs_m = chain['abs_m']
    for name, r in results.items():
        for k, j in enumerate(abs_j):
            rows.append({
                'M': M, 'prior': name,
                'j': int(j), 'm_j': int(abs_m[k]),
                'P_j': float(r['P_j'][k]),
                'W_j': float(r['W_j'][k]),
                'msS_j': float(r['msS_j'][k]),
                'sd_sS_j': float(r['sd_sS_j'][k]),
                'var_sS_j': float(r['var_sS_j'][k]),
                'mlog_j': float(r['mlog_j'][k]),
                'P_term': float(r['P_term']),
                'P_escape': float(r['P_escape']),
                'eps_S': float(r['eps_S']),
            })
    return rows


def report_var_compare(M, chain, results, emp_table, N_log2_emp):
    """Compare per-j chain Var to empirical Var at N=2^N_log2_emp."""
    abs_j = chain['abs_j']
    abs_m = chain['abs_m']
    print(f"\n=== Per-j Var[σ_S|j]: chain (M={M:,}) vs empirical (N=2^{N_log2_emp}) ===")
    for name in ['log_uniform', 'uniform_top_half']:
        if name not in results:
            continue
        r = results[name]
        print(f"\n  Chain prior: {name}")
        print(f"    {'j':>3} {'m_j':>9} {'⟨σ_S⟩_chn':>10} {'⟨σ_S⟩_emp':>10} "
              f"{'SD_chn':>9} {'SD_emp':>9} {'SD ratio':>10}")
        for k, j in enumerate(abs_j):
            jj = int(j)
            if jj not in emp_table:
                continue
            e = emp_table[jj]
            ratio = r['sd_sS_j'][k] / e['sd_sS'] if e['sd_sS'] > 0 else float('nan')
            print(f"    {jj:>3} {abs_m[k]:>9} {r['msS_j'][k]:>10.3f} {e['mean_sS']:>10.3f} "
                  f"{r['sd_sS_j'][k]:>9.3f} {e['sd_sS']:>9.3f} {ratio:>10.3f}")


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    all_rows = []
    chain_results_by_M = {}
    for M in [10**5, 10**6]:
        print(f"\n========== M = {M:,} ==========")
        t0 = time.perf_counter()
        chain = build_chain(M)
        print(f"  chain built in {time.perf_counter()-t0:.2f}s")

        B, tau = solve_chain(chain)
        results = aggregate(chain, B, tau)
        report(M, chain, results)
        all_rows.extend(to_csv_rows(M, chain, results))
        chain_results_by_M[M] = (chain, results)

    out_csv = out_dir / "45_absorbing_chain_results.csv"
    pl.DataFrame(all_rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")

    # === Per-j variance: empirical at N=2^24 (4M starts) for ground-truth comparison ===
    print(f"\n========== Empirical per-j Var[σ_S] at N=2^24 ==========")
    N_log2 = 24
    n_starts = 4_000_000
    emp_table = empirical_var_at_N(N_log2, n_starts)
    print(f"\n  Per-j empirical:")
    print(f"    {'j':>3} {'count':>10} {'⟨σ_S⟩':>10} {'SD':>9} {'⟨log m⟩':>9}")
    for j in sorted(emp_table.keys()):
        e = emp_table[j]
        print(f"    {j:>3} {e['cnt']:>10,} {e['mean_sS']:>10.3f} {e['sd_sS']:>9.3f} "
              f"{e['mean_log']:>9.3f}")

    for M in [10**5, 10**6]:
        chain, results = chain_results_by_M[M]
        report_var_compare(M, chain, results, emp_table, N_log2)

    # Save empirical table to CSV
    emp_rows = []
    for j, e in emp_table.items():
        emp_rows.append({
            'N_log2': N_log2,
            'j': int(j),
            'cnt': int(e['cnt']),
            'P_j': float(e['P']),
            'mean_sS': float(e['mean_sS']),
            'sd_sS': float(e['sd_sS']),
            'var_sS': float(e['var_sS']),
            'mean_log': float(e['mean_log']),
        })
    emp_csv = out_dir / "45_empirical_per_j_variance.csv"
    pl.DataFrame(emp_rows).write_csv(emp_csv)
    print(f"\n[save] {emp_csv}")


if __name__ == "__main__":
    main()
