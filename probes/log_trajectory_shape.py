"""
log_trajectory_shape.py — characterize per-orbit log-trajectory observables
per σ-quantile band: log_excursion, log_descent, peak_fraction.

For each orbit at N=2^32, record:
  - log_n_start, log_n_max, log_n_min, t_peak, σ_S
  - log_excursion = log(n_max) - log(n_start)
  - log_descent = log(n_start) - log(n_min)
  - peak_fraction = t_peak / σ_S

Per σ-band: fit parametric candidates (exp, gamma, lognormal) to log_excursion
and log_descent. Test Brownian-with-drift prediction (log_excursion ~ Exp).
Compute cross-observable correlations.

5 seeds × 1M orbits at N=2^32. Compute: <30s.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_with_trajectory(starts, max_steps=400000):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    n_max = np.zeros(n, dtype=np.int64)
    n_min = np.zeros(n, dtype=np.int64)
    t_peak = np.zeros(n, dtype=np.int32)
    failed = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = starts[i]
        n_max_i = m
        n_min_i = m
        t_peak_i = 0
        steps = 0
        j_attr = -1
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                failed[i] = True
                break
            three_m = 3 * m + 1
            v = 0
            while (three_m & 1) == 0:
                three_m >>= 1
                v += 1
            steps += 1
            if three_m > n_max_i:
                n_max_i = three_m
                t_peak_i = steps
            if three_m != 1 and three_m < n_min_i:
                n_min_i = three_m
            if three_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = three_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
        n_max[i] = n_max_i
        n_min[i] = n_min_i
        t_peak[i] = t_peak_i
    return sigma_arr, j_arr, n_max, n_min, t_peak, failed


def fit_distributions(data, max_n=100000):
    """Fit Exp, Gamma, Lognormal to data. Return log-likelihoods + params."""
    if len(data) > max_n:
        data = np.random.choice(data, size=max_n, replace=False)
    data = data[data > 0]  # positive values only
    if len(data) < 20:
        return None
    out = {}
    # Exp(loc=0)
    try:
        loc, scale = stats.expon.fit(data, floc=0)
        ll_exp = stats.expon.logpdf(data, loc, scale).sum()
        out['exp'] = {'rate': 1/scale, 'mean': scale, 'll': ll_exp}
    except Exception:
        out['exp'] = None
    # Gamma(loc=0)
    try:
        a, loc, scale = stats.gamma.fit(data, floc=0)
        ll_gamma = stats.gamma.logpdf(data, a, loc, scale).sum()
        out['gamma'] = {'a': a, 'scale': scale, 'mean': a*scale, 'll': ll_gamma}
    except Exception:
        out['gamma'] = None
    # Lognormal
    try:
        s, loc, scale = stats.lognorm.fit(data, floc=0)
        ll_logn = stats.lognorm.logpdf(data, s, loc, scale).sum()
        out['lognorm'] = {'s': s, 'scale': scale, 'mu': np.log(scale), 'sigma': s, 'll': ll_logn}
    except Exception:
        out['lognorm'] = None
    out['empirical_mean'] = float(data.mean())
    out['empirical_sd'] = float(data.std(ddof=1))
    out['n'] = len(data)
    return out


def main():
    log2N = 32
    N = 1 << log2N
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000
    q_edges = [0.0, 0.25, 0.50, 0.75, 0.95, 1.0]
    q_labels = [0.125, 0.375, 0.625, 0.875, 0.975]

    print(f"# Log-trajectory shape characterization at N=2^{log2N}")
    print(f"# {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# Bands at σ-quantile midpoints {q_labels}")
    print()

    # Aggregate across seeds
    band_data = {q: {'log_exc': [], 'log_desc': [], 'peak_frac': [], 'sigma': [], 'log_m': []}
                  for q in q_labels}

    t0_total = time.perf_counter()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, j_arr, n_max, n_min, t_peak, failed = walk_with_trajectory(starts)
        elapsed = time.perf_counter() - t0
        valid = ~failed
        starts_v = starts[valid]
        sigma = sigma[valid].astype(np.float64)
        n_max = n_max[valid].astype(np.float64)
        n_min = n_min[valid].astype(np.float64)
        t_peak = t_peak[valid].astype(np.float64)
        log_n_start = np.log(starts_v.astype(np.float64))
        log_n_max = np.log(n_max)
        log_n_min = np.log(n_min)
        log_excursion = log_n_max - log_n_start
        log_descent = log_n_start - log_n_min
        peak_fraction = t_peak / np.maximum(sigma, 1.0)

        # σ-quantile bands at q_edges
        edges = np.quantile(sigma, q_edges[1:-1])
        band_idx = np.digitize(sigma, edges)  # 0..4

        for b_idx, q in enumerate(q_labels):
            mask_b = band_idx == b_idx
            band_data[q]['log_exc'].append(log_excursion[mask_b])
            band_data[q]['log_desc'].append(log_descent[mask_b])
            band_data[q]['peak_frac'].append(peak_fraction[mask_b])
            band_data[q]['sigma'].append(sigma[mask_b])
            band_data[q]['log_m'].append(log_n_start[mask_b])

        print(f"  seed={seed:>5} ({elapsed:.1f}s)  valid={valid.sum():,}, failed={failed.sum()}")

    print(f"\n# Total compute: {time.perf_counter()-t0_total:.1f}s\n")

    # Concatenate per band
    for q in q_labels:
        for k in band_data[q]:
            band_data[q][k] = np.concatenate(band_data[q][k])

    # Per-band summary stats
    print(f"# Per-band log_excursion (= log(n_max) - log(n_start)) summary:")
    print(f"#   {'q':>6}  {'mean':>9}  {'sd':>9}  {'min':>9}  {'max':>9}  {'p50':>9}  {'p95':>9}")
    for q in q_labels:
        x = band_data[q]['log_exc']
        print(f"#   {q:>6.3f}  {x.mean():>9.4f}  {x.std(ddof=1):>9.4f}  {x.min():>9.4f}  "
              f"{x.max():>9.4f}  {np.percentile(x, 50):>9.4f}  {np.percentile(x, 95):>9.4f}")

    print(f"\n# Per-band log_descent (= log(n_start) - log(n_min)) summary:")
    print(f"#   {'q':>6}  {'mean':>9}  {'sd':>9}  {'min':>9}  {'max':>9}  {'p50':>9}  {'p95':>9}")
    for q in q_labels:
        x = band_data[q]['log_desc']
        print(f"#   {q:>6.3f}  {x.mean():>9.4f}  {x.std(ddof=1):>9.4f}  {x.min():>9.4f}  "
              f"{x.max():>9.4f}  {np.percentile(x, 50):>9.4f}  {np.percentile(x, 95):>9.4f}")

    print(f"\n# Per-band peak_fraction (= t_peak / σ_S) summary:")
    print(f"#   {'q':>6}  {'mean':>9}  {'sd':>9}  {'min':>9}  {'max':>9}  {'p50':>9}")
    for q in q_labels:
        x = band_data[q]['peak_frac']
        print(f"#   {q:>6.3f}  {x.mean():>9.4f}  {x.std(ddof=1):>9.4f}  {x.min():>9.4f}  "
              f"{x.max():>9.4f}  {np.percentile(x, 50):>9.4f}")

    # Parametric fits per band for log_excursion
    print(f"\n# Parametric fits for log_excursion (sample size capped at 100K per band):")
    print(f"#   Best fit by AIC (lower = better, ΔAIC = 2(ll_best - ll_x))")
    print(f"#   {'q':>6}  {'best':>9}  {'exp_mean':>9}  {'gamma_a':>9}  {'lognorm_μ':>10}  {'lognorm_σ':>10}  ΔAIC_exp  ΔAIC_lognorm")
    fit_results = {}
    for q in q_labels:
        x = band_data[q]['log_exc']
        # log_excursion can be 0 (orbit doesn't excurse) — filter for fit
        x_pos = x[x > 0.001]
        f = fit_distributions(x_pos, max_n=100000)
        fit_results[q] = f
        if f is None:
            continue
        lls = {'exp': f['exp']['ll'] if f['exp'] else -1e9,
                'gamma': f['gamma']['ll'] if f['gamma'] else -1e9,
                'lognorm': f['lognorm']['ll'] if f['lognorm'] else -1e9}
        # Penalize by # params: exp=1, gamma=2, lognorm=2
        n_params = {'exp': 1, 'gamma': 2, 'lognorm': 2}
        aics = {k: 2*n_params[k] - 2*v for k, v in lls.items()}
        best = min(aics, key=aics.get)
        d_aic_exp = aics['exp'] - aics[best]
        d_aic_logn = aics['lognorm'] - aics[best]
        e_str = f"{f['exp']['mean']:.4f}" if f['exp'] else 'NA'
        g_str = f"{f['gamma']['a']:.4f}" if f['gamma'] else 'NA'
        ln_mu = f"{f['lognorm']['mu']:.4f}" if f['lognorm'] else 'NA'
        ln_s = f"{f['lognorm']['sigma']:.4f}" if f['lognorm'] else 'NA'
        print(f"#   {q:>6.3f}  {best:>9}  {e_str:>9}  {g_str:>9}  "
              f"{ln_mu:>10}  {ln_s:>10}  {d_aic_exp:>+8.0f}  {d_aic_logn:>+10.0f}")

    # Parametric fits per band for log_descent
    print(f"\n# Parametric fits for log_descent:")
    print(f"#   {'q':>6}  {'best':>9}  {'exp_mean':>9}  {'gamma_a':>9}  {'lognorm_μ':>10}  {'lognorm_σ':>10}")
    desc_fit_results = {}
    for q in q_labels:
        x = band_data[q]['log_desc']
        x_pos = x[x > 0.001]
        f = fit_distributions(x_pos, max_n=100000)
        desc_fit_results[q] = f
        if f is None:
            continue
        lls = {'exp': f['exp']['ll'] if f['exp'] else -1e9,
                'gamma': f['gamma']['ll'] if f['gamma'] else -1e9,
                'lognorm': f['lognorm']['ll'] if f['lognorm'] else -1e9}
        n_params = {'exp': 1, 'gamma': 2, 'lognorm': 2}
        aics = {k: 2*n_params[k] - 2*v for k, v in lls.items()}
        best = min(aics, key=aics.get)
        e_str = f"{f['exp']['mean']:.4f}" if f['exp'] else 'NA'
        g_str = f"{f['gamma']['a']:.4f}" if f['gamma'] else 'NA'
        ln_mu = f"{f['lognorm']['mu']:.4f}" if f['lognorm'] else 'NA'
        ln_s = f"{f['lognorm']['sigma']:.4f}" if f['lognorm'] else 'NA'
        print(f"#   {q:>6.3f}  {best:>9}  {e_str:>9}  {g_str:>9}  {ln_mu:>10}  {ln_s:>10}")

    # Cross-observable correlations per band
    print(f"\n# Cross-observable correlations per band (Pearson):")
    print(f"#   {'q':>6}  {'ρ(exc,desc)':>12}  {'ρ(exc,σ)':>10}  {'ρ(peak,σ)':>10}  {'ρ(peak,exc)':>12}")
    for q in q_labels:
        e = band_data[q]['log_exc']
        d = band_data[q]['log_desc']
        s = band_data[q]['sigma']
        p = band_data[q]['peak_frac']
        if len(e) < 100:
            continue
        r_ed = np.corrcoef(e, d)[0, 1]
        r_es = np.corrcoef(e, s)[0, 1]
        r_ps = np.corrcoef(p, s)[0, 1]
        r_pe = np.corrcoef(p, e)[0, 1]
        print(f"#   {q:>6.3f}  {r_ed:>+12.4f}  {r_es:>+10.4f}  {r_ps:>+10.4f}  {r_pe:>+12.4f}")

    # Save summary
    rows = []
    for q in q_labels:
        f_exc = fit_results.get(q)
        f_desc = desc_fit_results.get(q)
        rows.append({
            'q': q,
            'n': int(len(band_data[q]['log_exc'])),
            'log_exc_mean': float(band_data[q]['log_exc'].mean()),
            'log_exc_sd': float(band_data[q]['log_exc'].std(ddof=1)),
            'log_desc_mean': float(band_data[q]['log_desc'].mean()),
            'log_desc_sd': float(band_data[q]['log_desc'].std(ddof=1)),
            'peak_frac_mean': float(band_data[q]['peak_frac'].mean()),
            'peak_frac_sd': float(band_data[q]['peak_frac'].std(ddof=1)),
            'exc_exp_rate': 1.0 / f_exc['exp']['mean'] if f_exc and f_exc['exp'] else None,
            'exc_lognorm_mu': f_exc['lognorm']['mu'] if f_exc and f_exc['lognorm'] else None,
            'exc_lognorm_sigma': f_exc['lognorm']['sigma'] if f_exc and f_exc['lognorm'] else None,
            'desc_exp_rate': 1.0 / f_desc['exp']['mean'] if f_desc and f_desc['exp'] else None,
            'desc_lognorm_mu': f_desc['lognorm']['mu'] if f_desc and f_desc['lognorm'] else None,
            'desc_lognorm_sigma': f_desc['lognorm']['sigma'] if f_desc and f_desc['lognorm'] else None,
        })
    pl.DataFrame(rows).write_csv("experiments_output/log_trajectory_shape.csv")
    print(f"\n[save] experiments_output/log_trajectory_shape.csv")


if __name__ == "__main__":
    main()
