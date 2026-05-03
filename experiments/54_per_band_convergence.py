"""
Experiment 54 — per-band K_eff_band(q, N) convergence-law characterization.

Computes K_eff_band(q) bootstrap (5 seeds) at N ∈ {2²⁵, 2²⁷, 2²⁸, 2³⁰, 2³², 2³⁴, 2³⁶}
for q ∈ {0.125, 0.375, 0.625, 0.875, 0.975}.

For each q: fit power-law (gap = A·N^(-α)), log-power (A·(log N)^(-α)),
and plateau (K = K_∞ - a·N^(-β)) to K_eff_band(q, N) vs N. Compare R² and
extracted exponents across q. If exponents cluster, K_full convergence is
band-weighted average of per-band convergence with shared dynamics. If they
differ across q, per-band dynamics are heterogeneous and K_full is shaped
by both per-band convergence and band-weight evolution.

Bands at user's 5 q midpoints:
  0.125 = 0-25 percentile
  0.375 = 25-50 percentile
  0.625 = 50-75 percentile
  0.875 = 75-100 percentile (whole top quartile)
  0.975 = 95-100 percentile (top 5% sub-band)

The 5th band is a sub-band of the 4th — they're not disjoint, but each gives
an independent K_eff_band measurement.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import norm
from scipy.optimize import minimize

sys.stdout.reconfigure(encoding="utf-8")

K_H_USER = 3.0 / (np.log(4.0) - np.log(3.0))
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_K(starts, thresholds_per_n, max_value, max_steps):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0; next_idx = 0
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0; next_idx += 1
        failed = False
        while m != 1 and steps < max_steps:
            if m & 1:
                if m > max_value // 3:
                    failed = True; break
                m = 3 * m + 1
            else:
                m = m >> 1
            steps += 1
            while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                s_arr[i, next_idx] = steps; next_idx += 1
        if not failed and m == 1:
            sigma_arr[i] = steps
            ok_arr[i] = True
    return sigma_arr, s_arr, ok_arr


def K_eff_slope(sigma, s_per_thresh, thresh_phys):
    log_f = np.log(np.maximum(thresh_phys, 1.0))
    R = sigma[:, None] - s_per_thresh.astype(np.float64)
    x = log_f.mean(axis=0); y = R.mean(axis=0)
    xc = x - x.mean(); yc = y - y.mean()
    return float((xc * yc).sum() / (xc * xc).sum())


def per_band_K_one_seed(N, seed, n_starts):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
    log_n = np.log(starts.astype(np.float64))
    sqrt_n = np.sqrt(starts.astype(np.float64))
    nt = starts.astype(np.float64) ** (2.0 / 3.0)
    sl = sqrt_n * log_n
    sdl = sqrt_n / np.maximum(log_n, 1.0)
    raw = np.column_stack([nt, sl, sqrt_n, sdl])
    sort_idx = np.argsort(-raw, axis=1)
    th = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    sigma, s, ok = walk_K(starts, th, MAX_VAL, 200_000)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros((len(starts), 4), dtype=np.int32)
    for col in range(4):
        s_phys[:, col] = np.take_along_axis(s, inv_sort[:, col:col+1], axis=1).flatten()
    sigma_ok = sigma[ok].astype(np.float64)
    s_ok = s_phys[ok]
    raw_ok = raw[ok]

    # K_full
    K_full = K_eff_slope(sigma_ok, s_ok, raw_ok)

    # 5 bands at user's midpoints
    edges = np.percentile(sigma_ok, [25, 50, 75, 95])  # 4 edges
    masks = [
        sigma_ok <= edges[0],                          # 0-25 → q=0.125
        (sigma_ok > edges[0]) & (sigma_ok <= edges[1]),# 25-50 → q=0.375
        (sigma_ok > edges[1]) & (sigma_ok <= edges[2]),# 50-75 → q=0.625
        sigma_ok > edges[2],                           # 75-100 → q=0.875
        sigma_ok > edges[3],                           # 95-100 → q=0.975
    ]
    K_band = np.array([K_eff_slope(sigma_ok[m], s_ok[m], raw_ok[m]) for m in masks])
    return K_full, K_band


def fit_power(log_N, K_arr):
    K_h = K_H_USER
    gap = K_h - K_arr
    pos = gap > 0
    if pos.sum() < 3: return None
    log_gap = np.log(gap[pos])
    lN = log_N[pos]
    xc = lN - lN.mean(); yc = log_gap - log_gap.mean()
    slope = float((xc*yc).sum() / (xc*xc).sum())
    alpha = -slope
    A = float(np.exp(log_gap.mean() - slope*lN.mean()))
    pred = K_h - A*np.exp(-alpha*log_N)
    R2 = 1 - float(np.sum((K_arr - pred)**2) / np.sum((K_arr - K_arr.mean())**2))
    return {'alpha': alpha, 'A': A, 'R2': R2, 'pred': pred,
            'resid_sd': float(np.std(K_arr - pred, ddof=1))}


def fit_logpower(log_N, K_arr):
    K_h = K_H_USER
    gap = K_h - K_arr
    pos = gap > 0
    if pos.sum() < 3: return None
    log_gap = np.log(gap[pos])
    lN = log_N[pos]
    log_lN = np.log(lN)
    xc = log_lN - log_lN.mean(); yc = log_gap - log_gap.mean()
    slope = float((xc*yc).sum() / (xc*xc).sum())
    alpha = -slope
    A = float(np.exp(log_gap.mean() - slope*log_lN.mean()))
    pred = K_h - A*np.exp(-alpha*np.log(log_N))
    R2 = 1 - float(np.sum((K_arr - pred)**2) / np.sum((K_arr - K_arr.mean())**2))
    return {'alpha': alpha, 'A': A, 'R2': R2, 'pred': pred,
            'resid_sd': float(np.std(K_arr - pred, ddof=1))}


def fit_plateau(log_N, K_arr):
    """K = K_inf - a · N^(-β)  with K_inf possibly above OR below K_h."""
    def loss(params):
        K_inf, a, beta = params
        if beta <= 0 or beta > 2.0: return 1e9
        pred = K_inf - a * np.exp(-beta * log_N)
        return float(np.sum((K_arr - pred)**2))
    # Try multiple starts to avoid local minima
    best = None
    K_h = K_H_USER
    for K_inf_init in [K_arr.max() + 0.1, K_h, K_h + 0.1]:
        for a_init in [10, 100, 1000]:
            for b_init in [0.1, 0.3, 0.5, 0.8]:
                res = minimize(loss, x0=[K_inf_init, a_init, b_init],
                               method='Nelder-Mead',
                               options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 5000})
                if best is None or res.fun < best.fun:
                    best = res
    K_inf, a, beta = best.x
    pred = K_inf - a*np.exp(-beta*log_N)
    R2 = 1 - float(np.sum((K_arr - pred)**2) / np.sum((K_arr - K_arr.mean())**2))
    return {'K_inf': float(K_inf), 'a': float(a), 'beta': float(beta),
            'R2': R2, 'pred': pred,
            'resid_sd': float(np.std(K_arr - pred, ddof=1))}


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N_list = [25, 27, 28, 30, 32, 34, 36]
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 500_000
    q_list = [0.125, 0.375, 0.625, 0.875, 0.975]
    z_q = np.array([norm.ppf(q) for q in q_list])

    print(f"# Per-band K_eff_band bootstrap")
    print(f"# {len(log2N_list)} N × {len(seeds)} seeds × {len(q_list)} bands")
    print(f"# n_starts per cell = {n_starts:,}")

    # Storage: K_band[N_idx][seed_idx][q_idx]
    K_band_all = np.zeros((len(log2N_list), len(seeds), len(q_list)))
    K_full_all = np.zeros((len(log2N_list), len(seeds)))

    t_total0 = time.perf_counter()
    for i, log2N in enumerate(log2N_list):
        N = 1 << log2N
        for j, seed in enumerate(seeds):
            t0 = time.perf_counter()
            K_full, K_band = per_band_K_one_seed(N, seed, n_starts)
            K_full_all[i, j] = K_full
            K_band_all[i, j, :] = K_band
            elapsed = time.perf_counter() - t0
            print(f"  2^{log2N:>2} seed={seed:>5}: K_full={K_full:.4f}  "
                  f"K_band={[f'{k:.3f}' for k in K_band]}  ({elapsed:.1f}s)")
    print(f"\n# Total compute: {time.perf_counter()-t_total0:.1f}s")

    # Bootstrap means and SDs
    K_band_mean = K_band_all.mean(axis=1)  # (N, q)
    K_band_sd = K_band_all.std(axis=1, ddof=1)
    K_full_mean = K_full_all.mean(axis=1)
    K_full_sd = K_full_all.std(axis=1, ddof=1)

    print(f"\n# Bootstrap-mean K_eff_band(q, N):")
    print(f"#   {'log2N':>6}  {'K_full':>9}", end="")
    for q in q_list:
        print(f"  {'q='+f'{q:.3f}':>9}", end="")
    print()
    for i, log2N in enumerate(log2N_list):
        print(f"#   {log2N:>6}  {K_full_mean[i]:>9.4f}", end="")
        for k in range(len(q_list)):
            print(f"  {K_band_mean[i, k]:>9.4f}", end="")
        print()

    print(f"\n# Single-seed SDs (across 5 seeds):")
    print(f"#   {'log2N':>6}  {'K_full':>9}", end="")
    for q in q_list:
        print(f"  {'q='+f'{q:.3f}':>9}", end="")
    print()
    for i, log2N in enumerate(log2N_list):
        print(f"#   {log2N:>6}  {K_full_sd[i]:>9.4f}", end="")
        for k in range(len(q_list)):
            print(f"  {K_band_sd[i, k]:>9.4f}", end="")
        print()

    # Per-band fits
    log_N = np.array(log2N_list) * np.log(2)

    # Also fit K_full for reference
    fits_full_pl = fit_power(log_N, K_full_mean)
    fits_full_lp = fit_logpower(log_N, K_full_mean)
    fits_full_pt = fit_plateau(log_N, K_full_mean)

    print(f"\n# K_full fits (7 pts) for reference:")
    print(f"#   power-law:  α={fits_full_pl['alpha']:.4f}  R²={fits_full_pl['R2']:.4f}  resSD={fits_full_pl['resid_sd']:.4f}")
    print(f"#   log-power:  α={fits_full_lp['alpha']:.4f}  R²={fits_full_lp['R2']:.4f}  resSD={fits_full_lp['resid_sd']:.4f}")
    print(f"#   plateau:    K_∞={fits_full_pt['K_inf']:.4f}  β={fits_full_pt['beta']:.4f}  R²={fits_full_pt['R2']:.4f}  resSD={fits_full_pt['resid_sd']:.4f}")

    # Per-band fits
    print(f"\n# Per-band fits across 7 N values:")
    band_fits = {}
    for k, q in enumerate(q_list):
        K_band_q = K_band_mean[:, k]
        fit_pl = fit_power(log_N, K_band_q)
        fit_lp = fit_logpower(log_N, K_band_q)
        fit_pt = fit_plateau(log_N, K_band_q)
        band_fits[q] = {'pl': fit_pl, 'lp': fit_lp, 'pt': fit_pt, 'K': K_band_q}
        print(f"\n#   q = {q}  (z_q = {z_q[k]:+.3f}):")
        print(f"#     K_band trajectory: " +
              "  ".join([f'2^{log2N_list[i]}={K_band_q[i]:.3f}' for i in range(len(log2N_list))]))
        if fit_pl:
            print(f"#     power-law:  α={fit_pl['alpha']:+.4f}  R²={fit_pl['R2']:+.4f}  resSD={fit_pl['resid_sd']:.4f}")
        else:
            print(f"#     power-law:  cannot fit (K_band > K_h or insufficient points)")
        if fit_lp:
            print(f"#     log-power:  α={fit_lp['alpha']:+.4f}  R²={fit_lp['R2']:+.4f}  resSD={fit_lp['resid_sd']:.4f}")
        else:
            print(f"#     log-power:  cannot fit")
        print(f"#     plateau:    K_∞={fit_pt['K_inf']:+.4f}  β={fit_pt['beta']:+.4f}  R²={fit_pt['R2']:+.4f}  resSD={fit_pt['resid_sd']:.4f}")

    # Cross-band exponent comparison
    print(f"\n# Per-band fitted exponents (do they cluster or differ?):")
    print(f"#   {'q':>6}  {'PL α':>9}  {'LG α':>9}  {'plat β':>9}  {'plat K_∞':>10}")
    for q in q_list:
        f = band_fits[q]
        pl_alpha = f['pl']['alpha'] if f['pl'] else float('nan')
        lp_alpha = f['lp']['alpha'] if f['lp'] else float('nan')
        pt_beta = f['pt']['beta']
        pt_K_inf = f['pt']['K_inf']
        print(f"#   {q:>6.3f}  {pl_alpha:>+9.4f}  {lp_alpha:>+9.4f}  {pt_beta:>+9.4f}  {pt_K_inf:>+10.4f}")

    # Verdict
    plat_betas = np.array([band_fits[q]['pt']['beta'] for q in q_list])
    plat_K_infs = np.array([band_fits[q]['pt']['K_inf'] for q in q_list])
    print(f"\n# Cross-band β spread: {plat_betas.max() - plat_betas.min():.4f}, mean = {plat_betas.mean():.4f}")
    print(f"# Cross-band K_∞ spread: {plat_K_infs.max() - plat_K_infs.min():.4f}, mean = {plat_K_infs.mean():.4f}")

    # Save
    rows = []
    for i, log2N in enumerate(log2N_list):
        for j, seed in enumerate(seeds):
            for k, q in enumerate(q_list):
                rows.append({
                    'log2N': log2N, 'N': 1 << log2N,
                    'seed': seed, 'q': float(q),
                    'K_band': float(K_band_all[i, j, k]),
                    'K_full': float(K_full_all[i, j]),
                })
    out_csv = out_dir / "54_per_band_convergence.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")

    # Save fit summary
    fit_rows = []
    for k, q in enumerate(q_list):
        f = band_fits[q]
        fit_rows.append({
            'q': float(q),
            'pl_alpha': f['pl']['alpha'] if f['pl'] else None,
            'pl_R2': f['pl']['R2'] if f['pl'] else None,
            'pl_resSD': f['pl']['resid_sd'] if f['pl'] else None,
            'lp_alpha': f['lp']['alpha'] if f['lp'] else None,
            'lp_R2': f['lp']['R2'] if f['lp'] else None,
            'lp_resSD': f['lp']['resid_sd'] if f['lp'] else None,
            'plat_K_inf': f['pt']['K_inf'],
            'plat_beta': f['pt']['beta'],
            'plat_a': f['pt']['a'],
            'plat_R2': f['pt']['R2'],
            'plat_resSD': f['pt']['resid_sd'],
        })
    pl.DataFrame(fit_rows).write_csv(out_dir / "54_per_band_fits.csv")
    print(f"[save] {out_dir / '54_per_band_fits.csv'}")


if __name__ == "__main__":
    main()
