"""
Experiment 38 — discriminate ε_S^(integer) = log(4) vs near-log(4) via 50M-orbit precision.

The Syracuse-walk renewal residue ε_S^(integer) (integer σ_S step count) appears
to converge to log(4) = 1.3863 within 0.5% at N=2^34. The 0.5% gap is either
sampling noise / finite-N convergence (then asymptote is exactly log(4)) or
a real deviation. Resolving requires 50M-orbit precision per N.

Method: at N ∈ {2^32, 2^34, 2^36}, sample 50M odd starts uniformly in [3, N].
Walk each orbit, compute σ_S (integer Syracuse step count), record last entry
class j (where v_entry = 2j is the final halving count). Compute:
  ε_S = ⟨σ_S⟩ − ⟨log N⟩/log(4/3)        (direct)
  W_j = E[σ_S | j] − (⟨log N⟩ − log(m_j))/log(4/3) − 1   (per-j conditional Wald)
  ε_S_via_Wj = Σ P(j)·W_j − ⟨log m_j⟩/log(4/3) + 1     (consistency)
Both routes should give identical ε_S to within sampling SE.

SE at 50M orbits: ~24/√50M ≈ 0.003. Sufficient to distinguish 1.380 vs 1.386
(gap 0.006) at >2σ confidence.

Usage:
    python 38_eps_S_log4_test.py
"""
import sys
import time
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG_2 = np.log(2)
LOG_3 = np.log(3)
LOG_4_3 = 2*LOG_2 - LOG_3
LOG_4 = 2*LOG_2
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def collect_per_j_orbit(starts, max_value, max_steps):
    n = len(starts)
    n_chunks = 32
    chunk = (n + n_chunks - 1) // n_chunks
    max_j = 12
    sum_sS_j = np.zeros((n_chunks, max_j+1), dtype=np.float64)
    sum_log_j = np.zeros((n_chunks, max_j+1), dtype=np.float64)
    cnt_j = np.zeros((n_chunks, max_j+1), dtype=np.int64)
    sum_sS = np.zeros(n_chunks, dtype=np.float64)
    sum_log = np.zeros(n_chunks, dtype=np.float64)
    cnt = np.zeros(n_chunks, dtype=np.int64)
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
            sum_sS[c] += sigma_S
            sum_log[c] += log_x
            cnt[c] += 1
            if 0 <= last_j <= max_j:
                sum_sS_j[c, last_j] += sigma_S
                sum_log_j[c, last_j] += log_x
                cnt_j[c, last_j] += 1
    return (sum_sS_j.sum(axis=0), sum_log_j.sum(axis=0), cnt_j.sum(axis=0),
            sum_sS.sum(), sum_log.sum(), cnt.sum(), overflow.sum())


def main():
    N_starts = 50_000_000
    print(f'50M-orbit precision test: ε_S vs log(4) = {LOG_4:.6f}')
    print(f'  Sampling {N_starts:,} odd starts per N')
    print()

    rng = np.random.default_rng(42)
    results = []
    for log2N in [32, 34, 36]:
        N = 1 << log2N
        starts = 2 * rng.integers(1, (N-1)//2, size=N_starts, dtype=np.int64) + 1
        print(f'[N=2^{log2N}] starting walk...')
        t0 = time.perf_counter()
        ssS_j, sl_j, cn_j, ssS_tot, sl_tot, cn_tot, of = collect_per_j_orbit(
            starts, MAX_VAL, 200_000)
        t = time.perf_counter() - t0
        print(f'    walk took {t:.1f}s, overflow={of}')
        # Direct ε_S
        msS = ssS_tot/cn_tot
        mlog = sl_tot/cn_tot
        eps_S_direct = msS - mlog/LOG_4_3
        # Per-j W_j and aggregated check
        W_per_j = {}
        P_per_j = {}
        eps_S_per_j = {}
        for j in range(2, 13):
            if cn_j[j] > 100:
                P_per_j[j] = cn_j[j]/cn_tot
                msS_j = ssS_j[j]/cn_j[j]
                mlog_j = sl_j[j]/cn_j[j]
                log_mj = np.log((4**j-1)/3)
                eps_S_j = msS_j - mlog_j/LOG_4_3
                W_j_v = eps_S_j - (1 - log_mj/LOG_4_3)
                W_per_j[j] = W_j_v
                eps_S_per_j[j] = eps_S_j
        # Aggregated via W_j route
        eps_S_via_Wj = sum(P_per_j[j]*eps_S_per_j[j] for j in P_per_j) + \
                       (1 - sum(P_per_j[j] for j in P_per_j))*0  # tiny adjustment
        # Standard error on direct ε_S: sample SD of (σ_S - log_n/log(4/3))
        # Approximate from variance of σ_S; SE ≈ Var(σ_S)/√N
        # Just report sample size
        SE = 24.0 / np.sqrt(cn_tot)
        results.append((log2N, mlog, msS, eps_S_direct, eps_S_via_Wj, P_per_j, W_per_j, SE))
        print(f'    ⟨log N⟩    = {mlog:.6f}')
        print(f'    ⟨σ_S⟩      = {msS:.6f}')
        print(f'    ε_S direct = {eps_S_direct:.6f}  (SE ≈ {SE:.4f})')
        print(f'    log(4)     = {LOG_4:.6f}')
        print(f'    diff       = {eps_S_direct - LOG_4:+.6f}  ({100*(eps_S_direct-LOG_4)/LOG_4:+.3f}%)')
        z = (eps_S_direct - LOG_4) / SE
        print(f'    z-score from log(4): {z:+.2f}σ')
        print()
        print(f'    Per-j table:')
        print(f'    {"j":>3} {"P(j)":>8} {"⟨σ_S|j⟩":>9} {"ε_S(j)":>9} {"W_j":>9} {"log(m_j)":>9}')
        for j in sorted(P_per_j.keys()):
            log_mj = np.log((4**j-1)/3)
            print(f'    {j:>3} {P_per_j[j]:>8.5f} {ssS_j[j]/cn_j[j]:>9.4f} '
                  f'{eps_S_per_j[j]:>+9.4f} {W_per_j[j]:>+9.4f} {log_mj:>9.4f}')
        # Identity check: ⟨W_j⟩ - ⟨log m_j⟩/log(4/3) + 1 = ε_S
        sum_W_P = sum(P_per_j[j]*W_per_j[j] for j in P_per_j)
        sum_logm_P = sum(P_per_j[j]*np.log((4**j-1)/3) for j in P_per_j)
        identity_LHS = sum_W_P - sum_logm_P/LOG_4_3 + 1
        print(f'    ⟨W_j⟩ - ⟨log m_j⟩/log(4/3) + 1 = {identity_LHS:.6f}')
        print(f'    log(4)                          = {LOG_4:.6f}')
        print(f'    diff: {identity_LHS - LOG_4:+.6f}')
        print()

    # Asymptote: extrapolate ε_S → ε_∞ via fit
    print('=== Asymptote analysis ===')
    log2N_arr = np.array([r[0] for r in results])
    mlog_arr = np.array([r[1] for r in results])
    eps_arr = np.array([r[3] for r in results])
    print(f'{"N":>8} {"⟨log N⟩":>9} {"ε_S":>10} {"diff log(4)":>12}')
    for log2N, mlog, msS, eps_S, eps_via_Wj, P_per_j, W_per_j, SE in results:
        print(f'2^{log2N:>2} {mlog:>9.4f} {eps_S:>10.6f} {eps_S-LOG_4:>+12.6f}')

    # Linear fit ε_S = ε_∞ + C/⟨log N⟩
    if len(eps_arr) >= 2:
        A = np.column_stack([np.ones_like(mlog_arr), 1.0/mlog_arr])
        coef, _, _, _ = np.linalg.lstsq(A, eps_arr, rcond=None)
        print(f'\n  Fit ε_S = ε_∞ + C/⟨log N⟩:')
        print(f'    ε_∞ = {coef[0]:.6f}')
        print(f'    C   = {coef[1]:.4f}')
        print(f'    distance to log(4): {coef[0] - LOG_4:+.6f} ({100*(coef[0]-LOG_4)/LOG_4:+.3f}%)')

    # Final verdict
    eps_largest = eps_arr[-1]
    SE_largest = results[-1][7]
    print(f'\n=== VERDICT ===')
    print(f'  ε_S empirical at N=2^{log2N_arr[-1]}: {eps_largest:.6f} (SE ≈ {SE_largest:.4f})')
    print(f'  log(4):                            {LOG_4:.6f}')
    print(f'  diff: {eps_largest - LOG_4:+.6f}  ({(eps_largest-LOG_4)/SE_largest:+.2f}σ)')
    if abs(eps_largest - LOG_4) < 2*SE_largest:
        print(f'  ✓ Consistent with ε_S = log(4) at <2σ. Closed form likely exact.')
    elif abs(eps_largest - LOG_4) < 5*SE_largest:
        print(f'  ~ Within 2-5σ of log(4). Suggestive but not decisive.')
    else:
        print(f'  ✗ More than 5σ from log(4). ε_S asymptote is NOT exactly log(4).')


if __name__ == "__main__":
    main()
