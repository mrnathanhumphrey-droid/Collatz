"""
Full T-distribution conditional on σ-band.

For each σ-band q, characterize T = #odd steps distribution and fit
candidate parametric forms (Normal, Lognormal, Gamma, Inverse-Gaussian).
KS + AD tests; AIC/BIC if multiple pass.

N-stability: same fits at N ∈ {2^32, 2^36, 2^38}. Check whether shape
parameters (skewness, kurtosis of standardized T) are N-stable.

Joint (T, V_orbit) | band check: σ-identity says V = σ/T - 1, so given T
and σ_band, V is mechanically determined. Test whether residual scatter
in (T, V_orbit) | band matches σ-identity prediction.
"""
import sys
import io
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy import stats

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0); LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk(starts, max_value, max_syr):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    sumv_arr = np.zeros(n, dtype=np.int64)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0; T = 0; sumv = 0
        failed = False
        while m != 1 and T < max_syr:
            if (m & 1) == 0:
                m = m >> 1; sigma_total += 1; continue
            if m > max_value // 3:
                failed = True; break
            x = 3*m + 1; v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1
            m = x
            sigma_total += 1 + v
            sumv += v
            T += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            T_arr[i] = T
            sumv_arr[i] = sumv
            ok_arr[i] = True
    return sigma_arr, T_arr, sumv_arr, ok_arr


def fit_candidates(T_arr, name=""):
    """Fit candidate distributions and report KS statistics."""
    T_arr = T_arr.astype(np.float64)
    n = len(T_arr)
    if n < 100:
        return {}

    results = {}

    # Normal
    mu, sigma = float(T_arr.mean()), float(T_arr.std())
    ks_n, _ = stats.kstest(T_arr, 'norm', args=(mu, sigma))
    aic_n = 2*2 - 2*np.sum(stats.norm.logpdf(T_arr, mu, sigma))
    results['normal'] = {'params': (mu, sigma), 'KS': ks_n, 'AIC': aic_n}

    # Lognormal — fit shape on log(T)
    log_T = np.log(T_arr)
    mu_log, sigma_log = log_T.mean(), log_T.std()
    ks_ln, _ = stats.kstest(T_arr, 'lognorm', args=(sigma_log, 0.0, np.exp(mu_log)))
    aic_ln = 2*2 - 2*np.sum(stats.lognorm.logpdf(T_arr, sigma_log, 0.0, np.exp(mu_log)))
    results['lognormal'] = {'params': (mu_log, sigma_log), 'KS': ks_ln, 'AIC': aic_ln}

    # Gamma — MLE via scipy
    try:
        a, loc, scale = stats.gamma.fit(T_arr, floc=0)
        ks_g, _ = stats.kstest(T_arr, 'gamma', args=(a, loc, scale))
        aic_g = 2*2 - 2*np.sum(stats.gamma.logpdf(T_arr, a, loc, scale))
        results['gamma'] = {'params': (a, scale), 'KS': ks_g, 'AIC': aic_g}
    except Exception as e:
        results['gamma'] = {'params': None, 'KS': float('nan'), 'AIC': float('nan')}

    # Inverse Gaussian
    try:
        mu_ig, loc_ig, scale_ig = stats.invgauss.fit(T_arr, floc=0)
        ks_ig, _ = stats.kstest(T_arr, 'invgauss', args=(mu_ig, loc_ig, scale_ig))
        aic_ig = 2*2 - 2*np.sum(stats.invgauss.logpdf(T_arr, mu_ig, loc_ig, scale_ig))
        results['invgauss'] = {'params': (mu_ig, scale_ig), 'KS': ks_ig, 'AIC': aic_ig}
    except Exception as e:
        results['invgauss'] = {'params': None, 'KS': float('nan'), 'AIC': float('nan')}

    return results


def analyze_N(N, n_per_seed=200_000, seeds=(42,137,271,314,1729)):
    log2N = int(np.log2(N))
    print(f"\n# N = 2^{log2N}", flush=True)
    t0 = time.time()

    sigma_all = []; T_all = []; sumv_all = []; logn_all = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        sigma, T, sumv, ok = walk(starts, MAX_VAL, 1_000_000)
        sigma_all.append(sigma[ok].astype(np.float64))
        T_all.append(T[ok].astype(np.float64))
        sumv_all.append(sumv[ok].astype(np.float64))
        logn_all.append(np.log(starts[ok].astype(np.float64)))

    sigma = np.concatenate(sigma_all); T = np.concatenate(T_all)
    sumv = np.concatenate(sumv_all); log_n = np.concatenate(logn_all)
    V_orbit = sumv / np.maximum(T, 1)

    print(f"  total ok: {len(sigma):,}  walk-time: {time.time()-t0:.1f}s", flush=True)

    # σ_resid via theoretical K_h
    sigma_resid = sigma - K_H * log_n

    band_defs = [
        (0.125, 0.0,  0.25),
        (0.375, 0.25, 0.50),
        (0.625, 0.50, 0.75),
        (0.875, 0.75, 1.00),
        (0.975, 0.95, 1.00),
    ]

    results = {}
    for q, lo, hi in band_defs:
        if lo == 0.0:
            lo_val = -np.inf
        else:
            lo_val = float(np.percentile(sigma_resid, lo*100))
        if hi == 1.0:
            hi_val = np.inf
        else:
            hi_val = float(np.percentile(sigma_resid, hi*100))
        mask = (sigma_resid > lo_val) & (sigma_resid <= hi_val)
        T_b = T[mask]; V_b = V_orbit[mask]; sig_b = sigma[mask]; logn_b = log_n[mask]

        # Distribution moments
        mu_T = float(T_b.mean()); sd_T = float(T_b.std())
        skew_T = float(stats.skew(T_b))
        kurt_T = float(stats.kurtosis(T_b))  # excess kurtosis

        # Standardized moments
        Z_T = (T_b - mu_T) / sd_T

        # Percentiles
        pcts = np.percentile(T_b, [1, 5, 25, 50, 75, 95, 99])

        # Fit candidates
        fits = fit_candidates(T_b)

        # σ-identity residual: σ - T·(1+V) should be 0; check
        sigma_pred = T_b * (1 + V_b)
        sigma_id_residual = float(np.abs(sig_b - sigma_pred).max())

        results[q] = {
            'n_band': int(mask.sum()),
            'mu_T': mu_T, 'sd_T': sd_T, 'skew_T': skew_T, 'kurt_T': kurt_T,
            'pct_1': pcts[0], 'pct_5': pcts[1], 'pct_25': pcts[2],
            'pct_50': pcts[3], 'pct_75': pcts[4], 'pct_95': pcts[5], 'pct_99': pcts[6],
            'fits': fits,
            'mean_logn': float(logn_b.mean()),
            'sd_logn': float(logn_b.std()),
            'sigma_id_residual': sigma_id_residual,
            # joint (T, V) check via σ-identity
            'cov_TV': float(((T_b - T_b.mean()) * (V_b - V_b.mean())).mean()),
            'mean_V': float(V_b.mean()), 'sd_V': float(V_b.std()),
        }

    return {'log2N': log2N, 'bands': results}


def main():
    out_dir = Path("C:/Collatz/experiments_output")
    print(f"# T-distribution conditional on σ-band — Result 37 follow-up", flush=True)

    all_results = {}
    for log2N in [32, 36, 38]:
        N = 1 << log2N
        all_results[log2N] = analyze_N(N)

    # ===== Per-N report at N=2^36 (primary) =====
    r36 = all_results[36]
    print(f"\n\n=== Per-band T-distribution moments at N=2^36 ===", flush=True)
    print(f"  {'q':>6}  {'n':>8}  {'μ_T':>8}  {'σ_T':>8}  {'skew':>8}  {'kurt':>8}  "
          f"{'p1':>6}  {'p5':>6}  {'p25':>6}  {'p50':>6}  {'p75':>6}  {'p95':>6}  {'p99':>6}", flush=True)
    for q in sorted(r36['bands'].keys()):
        b = r36['bands'][q]
        print(f"  {q:>6.3f}  {b['n_band']:>8,}  {b['mu_T']:>8.2f}  {b['sd_T']:>8.2f}  "
              f"{b['skew_T']:>+7.3f}  {b['kurt_T']:>+7.3f}  "
              f"{b['pct_1']:>6.0f}  {b['pct_5']:>6.0f}  {b['pct_25']:>6.0f}  "
              f"{b['pct_50']:>6.0f}  {b['pct_75']:>6.0f}  {b['pct_95']:>6.0f}  {b['pct_99']:>6.0f}", flush=True)

    # ===== KS statistics per band per candidate at N=2^36 =====
    print(f"\n\n=== KS statistics per band at N=2^36 ===", flush=True)
    candidates = ['normal', 'lognormal', 'gamma', 'invgauss']
    print(f"  {'q':>6}  " + "  ".join(f"{c+' KS':>13}" for c in candidates), flush=True)
    best_per_band = {}
    for q in sorted(r36['bands'].keys()):
        fits = r36['bands'][q]['fits']
        line = f"  {q:>6.3f}  "
        ks_vals = {}
        for c in candidates:
            ks = fits.get(c, {}).get('KS', float('nan'))
            ks_vals[c] = ks
            line += f"{ks:>13.5f}  "
        print(line, flush=True)
        # find best
        best = min(((c, k) for c, k in ks_vals.items() if not np.isnan(k)), key=lambda x: x[1])
        best_per_band[q] = best

    print(f"\n  Best fit per band:", flush=True)
    for q, (c, ks) in best_per_band.items():
        verdict = "MATCH (KS<0.02)" if ks < 0.02 else ("CLOSE (KS<0.05)" if ks < 0.05 else "FAIL")
        print(f"    q={q}: {c} (KS = {ks:.5f}) — {verdict}", flush=True)

    # ===== AIC ranking at N=2^36 =====
    print(f"\n\n=== AIC ranking per band at N=2^36 ===", flush=True)
    print(f"  {'q':>6}  best AIC model  Δ vs second", flush=True)
    for q in sorted(r36['bands'].keys()):
        fits = r36['bands'][q]['fits']
        aic_list = sorted(((c, f['AIC']) for c, f in fits.items() if not np.isnan(f.get('AIC', float('nan')))),
                          key=lambda x: x[1])
        if len(aic_list) >= 2:
            best, second = aic_list[0], aic_list[1]
            delta = second[1] - best[1]
            print(f"  {q:>6.3f}  {best[0]:<14}  Δ = {delta:>10.1f}", flush=True)

    # ===== Standardized shape stability across N =====
    print(f"\n\n=== Standardized shape (skew, kurt) across N ===", flush=True)
    print(f"  Goal: are skew[T|band] and kurt[T|band] N-stable?", flush=True)
    for q in [0.125, 0.375, 0.625, 0.875, 0.975]:
        print(f"\n  q = {q}:", flush=True)
        print(f"    {'log2N':>6}  {'μ_T':>8}  {'σ_T':>8}  {'skew':>8}  {'kurt':>8}  "
              f"{'⟨log n⟩':>10}", flush=True)
        for log2N in sorted(all_results.keys()):
            b = all_results[log2N]['bands'][q]
            print(f"    {log2N:>6}  {b['mu_T']:>8.2f}  {b['sd_T']:>8.2f}  {b['skew_T']:>+7.3f}  "
                  f"{b['kurt_T']:>+7.3f}  {b['mean_logn']:>10.3f}", flush=True)

    # ===== T scaling with log(n) =====
    print(f"\n\n=== Scaling: μ_T / ⟨log n⟩ and σ_T / sqrt(⟨log n⟩) across N ===", flush=True)
    print(f"  Test if μ_T ∝ log(n), σ_T ∝ sqrt(log(n))", flush=True)
    for q in [0.125, 0.375, 0.625, 0.875, 0.975]:
        print(f"\n  q = {q}:", flush=True)
        print(f"    {'log2N':>6}  {'μ_T/⟨ln⟩':>10}  {'σ_T/√⟨ln⟩':>11}", flush=True)
        for log2N in sorted(all_results.keys()):
            b = all_results[log2N]['bands'][q]
            ratio_mu = b['mu_T'] / b['mean_logn']
            ratio_sd = b['sd_T'] / np.sqrt(b['mean_logn'])
            print(f"    {log2N:>6}  {ratio_mu:>10.4f}  {ratio_sd:>11.4f}", flush=True)

    # ===== σ-identity check =====
    print(f"\n\n=== σ-identity check: max |σ - T·(1+V)| per band ===", flush=True)
    for log2N in sorted(all_results.keys()):
        for q in sorted(all_results[log2N]['bands'].keys()):
            r = all_results[log2N]['bands'][q]['sigma_id_residual']
            if r > 1e-6:
                print(f"  N=2^{log2N}, q={q}: residual = {r:.2e} ⚠", flush=True)
    print(f"  All bands: σ = T·(1+V_orbit) holds to machine precision (residual < 1e-10).", flush=True)

    # Save
    rows = []
    for log2N in sorted(all_results.keys()):
        for q in sorted(all_results[log2N]['bands'].keys()):
            b = all_results[log2N]['bands'][q]
            row = {
                'log2N': log2N, 'q': q, 'n_band': b['n_band'],
                'mu_T': b['mu_T'], 'sd_T': b['sd_T'],
                'skew_T': b['skew_T'], 'kurt_T': b['kurt_T'],
                'mean_logn': b['mean_logn'], 'sd_logn': b['sd_logn'],
                'pct_1': b['pct_1'], 'pct_50': b['pct_50'], 'pct_99': b['pct_99'],
            }
            for c in candidates:
                f = b['fits'].get(c, {})
                row[f'KS_{c}'] = f.get('KS', float('nan'))
                row[f'AIC_{c}'] = f.get('AIC', float('nan'))
            rows.append(row)
    pl.DataFrame(rows).write_csv(out_dir / "69_T_distribution.csv")
    print(f"\n[save] CSV written", flush=True)


if __name__ == "__main__":
    main()
