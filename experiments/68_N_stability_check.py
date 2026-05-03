"""
N-stability re-measurement of P(q|j) Gibbs form (Result 37).

Walk N ∈ {2^32, 2^34, 2^36, 2^38}, 5 seeds × 1M orbits each.
For each N: compute P(q|j) at quartile bands × j ∈ {2, 4, 5}.
Fit Gibbs form: log P(q|j_a)/P(q|j_b) = a + α_diff·q.
Extract α(j) per N (anchor α(2)=0); track drift.

Verdict:
- If α(j) values stable within ±1% (or bootstrap): structural
- If drift > 5%: approximate at N tested

Also re-measure w_q at quartile bands for cross-check.
"""
import sys
import io
import time
import numpy as np
from numba import njit, prange
from scipy.stats import norm
import polars as pl
from pathlib import Path

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG2 = np.log(2.0); LOG3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk(starts, max_steps=400000):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    sumv_arr = np.zeros(n, dtype=np.int64)
    nodd_arr = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        m = starts[i]
        steps = 0
        T = 0
        sumv = 0
        j_attr = -1
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            v = 0
            while (three_m & 1) == 0:
                three_m >>= 1
                v += 1
            steps += 1 + v
            sumv += v
            T += 1
            if three_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = three_m
        if m == 1:
            sigma_arr[i] = steps
            j_arr[i] = j_attr
            sumv_arr[i] = sumv
            nodd_arr[i] = T
    return sigma_arr, j_arr, sumv_arr, nodd_arr


def E_v_from_w_pathc(w):
    return 1.0 / (1.0 - 2.0**(-(1.0 + w)))


def w_from_E_v_pathc(E_v):
    if E_v <= 1: return float('inf')
    r = 1.0 - 1.0/E_v
    if r <= 0: return float('inf')
    return -np.log2(r) - 1.0


def measure_one_N(log2N, n_per_seed=200_000, seeds=(42,137,271,314,1729)):
    N = 1 << log2N
    j_targets = [2, 4, 5]
    q_labels = [0.125, 0.375, 0.625, 0.875]

    P_qj_per_seed = {jt: {q: [] for q in q_labels} for jt in j_targets}
    E_v_band_per_seed = {q: [] for q in q_labels}
    w_q_per_seed = {q: [] for q in q_labels}

    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        sigma, j_arr, sumv, nodd = walk(starts)
        valid = sigma > 0
        sigma = sigma[valid]; j_arr = j_arr[valid]
        sumv = sumv[valid]; nodd = nodd[valid]

        # σ-quartile bands using global sigma distribution edges
        edges = np.quantile(sigma, [0.25, 0.50, 0.75])
        band_idx = np.digitize(sigma, edges)  # 0..3

        # P(q | j) per band
        for jt in j_targets:
            mask_j = j_arr == jt
            n_j = int(mask_j.sum())
            if n_j == 0: continue
            for b_idx, q in enumerate(q_labels):
                P_qj_per_seed[jt][q].append(((band_idx == b_idx) & mask_j).sum() / n_j)

        # Per-band E[v_t] = sumv / nodd within band, weighted (per-step rather than per-orbit)
        # Use sum(sumv[band]) / sum(nodd[band]) — per-step v_t mean within band
        for b_idx, q in enumerate(q_labels):
            mask_b = band_idx == b_idx
            sv = sumv[mask_b].sum()
            nt = nodd[mask_b].sum()
            E_v = float(sv) / float(nt) if nt > 0 else 0.0
            E_v_band_per_seed[q].append(E_v)
            # w_q via Path C inversion
            w = w_from_E_v_pathc(E_v) if E_v > 1.001 else float('nan')
            w_q_per_seed[q].append(w)

    walk_time = time.time() - t0

    # Aggregate
    P_qj = {jt: {q: (np.mean(P_qj_per_seed[jt][q]),
                     np.std(P_qj_per_seed[jt][q], ddof=1)/np.sqrt(len(P_qj_per_seed[jt][q])))
                 for q in q_labels} for jt in j_targets}
    E_v_band = {q: (np.mean(E_v_band_per_seed[q]),
                    np.std(E_v_band_per_seed[q], ddof=1)/np.sqrt(len(E_v_band_per_seed[q])))
                for q in q_labels}
    w_q = {q: (np.mean(w_q_per_seed[q]),
               np.std(w_q_per_seed[q], ddof=1)/np.sqrt(len(w_q_per_seed[q])))
           for q in q_labels}

    # Fit Gibbs form: log P(q|j_a)/P(q|j_b) = a + alpha_diff·q
    q_arr = np.array(q_labels)
    pair_results = {}
    for ja, jb in [(2,4), (4,5), (2,5)]:
        log_ratios = []
        for q in q_labels:
            log_ratios.append(np.log(P_qj[ja][q][0] / P_qj[jb][q][0]))
        log_ratios = np.array(log_ratios)
        # OLS in q
        Xc = q_arr - q_arr.mean(); Yc = log_ratios - log_ratios.mean()
        slope = (Xc*Yc).sum() / (Xc*Xc).sum()
        intercept = log_ratios.mean() - slope*q_arr.mean()
        pred = intercept + slope*q_arr
        ss_res = ((log_ratios - pred)**2).sum()
        ss_tot = ((log_ratios - log_ratios.mean())**2).sum()
        R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 1.0
        pair_results[(ja, jb)] = {'slope': float(slope), 'intercept': float(intercept), 'R2': float(R2)}

    # alpha(j) anchored: α(2) = 0
    # slope_24 = α(2) - α(4) = 3.02 → α(4) = -slope_24
    # slope_25 = α(2) - α(5) = 2.30 → α(5) = -slope_25
    alpha = {
        2: 0.0,
        4: -pair_results[(2,4)]['slope'],
        5: -pair_results[(2,5)]['slope'],
    }

    return {
        'log2N': log2N, 'walk_time_s': walk_time,
        'P_qj': P_qj, 'E_v_band': E_v_band, 'w_q': w_q,
        'pairs': pair_results, 'alpha': alpha,
    }


def main():
    out_dir = Path("C:/Collatz/experiments_output")
    print(f"# N-stability check of P(q|j) Gibbs form (Result 37 follow-up)", flush=True)
    print(f"# 5 seeds × 200k orbits × 4 N values = 4M orbits total\n", flush=True)

    results = {}
    for log2N in [32, 34, 36, 38]:
        print(f"# N = 2^{log2N}", flush=True)
        r = measure_one_N(log2N)
        results[log2N] = r
        print(f"  walked in {r['walk_time_s']:.1f}s", flush=True)

    # ===== α(j) drift across N =====
    print(f"\n\n=== α(j) Gibbs parameters across N ===", flush=True)
    print(f"  {'log2N':>6}  {'α(2)':>8}  {'α(4)':>10}  {'α(5)':>10}  "
          f"{'slope_24 R²':>13}  {'slope_45 R²':>13}  {'slope_25 R²':>13}", flush=True)
    for log2N in sorted(results.keys()):
        r = results[log2N]
        print(f"  {log2N:>6}  {r['alpha'][2]:>+7.4f}  {r['alpha'][4]:>+10.4f}  {r['alpha'][5]:>+10.4f}  "
              f"{r['pairs'][(2,4)]['R2']:>13.4f}  {r['pairs'][(4,5)]['R2']:>13.4f}  "
              f"{r['pairs'][(2,5)]['R2']:>13.4f}", flush=True)

    # Drift analysis
    alpha_4_arr = np.array([results[N]['alpha'][4] for N in sorted(results.keys())])
    alpha_5_arr = np.array([results[N]['alpha'][5] for N in sorted(results.keys())])
    drift_4 = (alpha_4_arr.max() - alpha_4_arr.min()) / abs(alpha_4_arr.mean()) * 100
    drift_5 = (alpha_5_arr.max() - alpha_5_arr.min()) / abs(alpha_5_arr.mean()) * 100
    print(f"\n  α(4) drift: {alpha_4_arr.min():+.4f} to {alpha_4_arr.max():+.4f} = {drift_4:.2f}% of |mean|", flush=True)
    print(f"  α(5) drift: {alpha_5_arr.min():+.4f} to {alpha_5_arr.max():+.4f} = {drift_5:.2f}% of |mean|", flush=True)

    # ===== P(q|j) raw values per N =====
    print(f"\n\n=== Raw P(q|j) per N — same q, j across rows ===", flush=True)
    for log2N in sorted(results.keys()):
        r = results[log2N]
        print(f"\n  N=2^{log2N}:", flush=True)
        print(f"    {'q':>6}  " + "  ".join(f"P(q|j={j})" for j in [2,4,5]), flush=True)
        for q in [0.125, 0.375, 0.625, 0.875]:
            line = f"    {q:>6.3f}  "
            for jt in [2,4,5]:
                m, se = r['P_qj'][jt][q]
                line += f"{m:.4f}±{se:.4f}  "
            print(line, flush=True)

    # ===== w_q drift =====
    print(f"\n\n=== w_q per N (per-step E[v] inversion) ===", flush=True)
    print(f"  {'log2N':>6}  " + "  ".join(f"{'w_q@'+str(q):>14}" for q in [0.125, 0.375, 0.625, 0.875]), flush=True)
    for log2N in sorted(results.keys()):
        r = results[log2N]
        line = f"  {log2N:>6}  "
        for q in [0.125, 0.375, 0.625, 0.875]:
            m, se = r['w_q'][q]
            line += f"{m:>+8.4f}±{se:.4f}  "
        print(line, flush=True)

    # Drift on w_q
    print(f"\n  w_q drift across N:", flush=True)
    for q in [0.125, 0.375, 0.625, 0.875]:
        wq_arr = np.array([results[N]['w_q'][q][0] for N in sorted(results.keys())])
        rng = wq_arr.max() - wq_arr.min()
        mean = abs(wq_arr.mean())
        pct = rng/mean*100 if mean > 1e-9 else float('nan')
        print(f"    q={q}: min={wq_arr.min():+.4f}  max={wq_arr.max():+.4f}  range={rng:.4f}  "
              f"({pct:.2f}% of |mean|)", flush=True)

    # ===== VERDICT =====
    print(f"\n\n=== VERDICT ===", flush=True)
    max_alpha_drift = max(drift_4, drift_5)
    if max_alpha_drift < 1.0:
        print(f"  α(j) drift < 1% → Gibbs form parameters STRUCTURAL across N", flush=True)
    elif max_alpha_drift < 5.0:
        print(f"  α(j) drift {max_alpha_drift:.2f}% — APPROXIMATELY structural; small N-drift", flush=True)
    else:
        print(f"  α(j) drift {max_alpha_drift:.2f}% > 5% — approximate at N tested, drifts with N", flush=True)

    # ===== Save =====
    rows = []
    for log2N in sorted(results.keys()):
        r = results[log2N]
        for q in [0.125, 0.375, 0.625, 0.875]:
            for jt in [2,4,5]:
                rows.append({
                    'log2N': log2N, 'q': q, 'j': jt,
                    'P_qj': r['P_qj'][jt][q][0], 'P_qj_se': r['P_qj'][jt][q][1],
                })
        for q in [0.125, 0.375, 0.625, 0.875]:
            rows.append({
                'log2N': log2N, 'q': q, 'j': -1,  # marker for w_q row
                'P_qj': r['w_q'][q][0], 'P_qj_se': r['w_q'][q][1],
            })
    pl.DataFrame(rows).write_csv(out_dir / "68_N_stability_check.csv")

    rows_alpha = []
    for log2N in sorted(results.keys()):
        r = results[log2N]
        rows_alpha.append({
            'log2N': log2N,
            'alpha_2': r['alpha'][2], 'alpha_4': r['alpha'][4], 'alpha_5': r['alpha'][5],
            'R2_24': r['pairs'][(2,4)]['R2'], 'R2_45': r['pairs'][(4,5)]['R2'], 'R2_25': r['pairs'][(2,5)]['R2'],
        })
    pl.DataFrame(rows_alpha).write_csv(out_dir / "68_alpha_per_N.csv")
    print(f"\n[save] CSVs written", flush=True)


if __name__ == "__main__":
    main()
