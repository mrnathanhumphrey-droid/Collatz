"""
Experiment 53 — K_full + ξ_X at N=2^36, 5-seed bootstrap. ONE additional
data point that constrains the K_full(N) and ξ_X(N) trajectories.

This is NOT validation of any model. It does:
  - Tighten plateau-fit CI on K_∞
  - Discriminate plateau vs power-law more clearly (different predictions here)
  - Either confirm or shift β estimate away from 1/2
  - Extend ξ_X rotation trajectory past 2^34's -0.10

Predictions at 2^36 from 6-point fits:
  Plateau:    10.268 ± 0.01
  Power-law:  10.279
  Log-power:  10.277
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

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


def gpd_factor(p, xi):
    if abs(xi) < 1e-9:
        return -np.log(1.0 - p)
    return ((1.0 - p) ** (-xi) - 1.0) / xi


def fit_xi_NLS(q_thresh, q_tail, excess):
    from scipy.optimize import minimize
    p = (q_tail - q_thresh) / (1.0 - q_thresh)
    def loss(params):
        xi, sigma = params
        if sigma <= 0 or abs(xi) > 0.5: return 1e20
        pred = sigma * np.array([gpd_factor(pi, xi) for pi in p])
        return float(np.sum((pred - excess) ** 2))
    res = minimize(loss, x0=[0.0, np.median(excess)], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-10})
    return float(res.x[0]), float(res.x[1])


def analyze_one_seed(N, seed, n_starts):
    from scipy.stats import norm
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
    s_phys_ok = s_phys[ok]
    raw_ok = raw[ok]
    K_full = K_eff_slope(sigma_ok, s_phys_ok, raw_ok)

    # Bands at midpoints {0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995}
    edges = np.percentile(sigma_ok, [25, 50, 75, 90, 95, 99])
    masks = [
        sigma_ok <= edges[0],
        (sigma_ok > edges[0]) & (sigma_ok <= edges[1]),
        (sigma_ok > edges[1]) & (sigma_ok <= edges[2]),
        (sigma_ok > edges[2]) & (sigma_ok <= edges[3]),
        (sigma_ok > edges[3]) & (sigma_ok <= edges[4]),
        (sigma_ok > edges[4]) & (sigma_ok <= edges[5]),
        sigma_ok > edges[5],
    ]
    q_arr = np.array([0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995])
    z_q = np.array([float(norm.ppf(q)) for q in q_arr])
    K_emp = np.array([K_eff_slope(sigma_ok[m], s_phys_ok[m], raw_ok[m]) for m in masks])
    body = (q_arr >= 0.30) & (q_arr <= 0.70)
    z_b = z_q[body]; K_b = K_emp[body]
    z_c = z_b - z_b.mean(); K_c = K_b - K_b.mean()
    b_fit = float((z_c * K_c).sum() / (z_c * z_c).sum())
    K_h_fit = float(K_b.mean() - b_fit * z_b.mean())
    X_thresh = norm.ppf(0.75)
    tail = q_arr >= 0.75
    excess_user = (K_emp[tail] - K_H_USER) / 2.275 - X_thresh
    excess_body = (K_emp[tail] - K_h_fit) / b_fit - X_thresh
    xi_user, _ = fit_xi_NLS(0.75, q_arr[tail], excess_user)
    xi_body, _ = fit_xi_NLS(0.75, q_arr[tail], excess_body)
    return K_full, xi_user, xi_body, int(ok.sum())


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N = 36
    N = 1 << log2N
    n_starts = 500_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# K_full + ξ_X at 2^{log2N} bootstrap, {len(seeds)} seeds × {n_starts:,} orbits")
    K_seeds = []
    xi_user_seeds = []
    xi_body_seeds = []
    for seed in seeds:
        t0 = time.perf_counter()
        K, xi_u, xi_b, n_ok = analyze_one_seed(N, seed, n_starts)
        K_seeds.append(K); xi_user_seeds.append(xi_u); xi_body_seeds.append(xi_b)
        print(f"  seed={seed:>5}: K_full={K:.5f}  ξ_user={xi_u:+.4f}  ξ_body={xi_b:+.4f}  "
              f"({time.perf_counter()-t0:.1f}s, {n_ok:,} ok)")

    K_arr = np.array(K_seeds)
    mean = K_arr.mean(); sd = K_arr.std(ddof=1); se = sd / np.sqrt(len(seeds))
    xi_user_arr = np.array(xi_user_seeds)
    xi_body_arr = np.array(xi_body_seeds)
    xi_user_mean = xi_user_arr.mean(); xi_body_mean = xi_body_arr.mean()
    xi_user_sd = xi_user_arr.std(ddof=1); xi_body_sd = xi_body_arr.std(ddof=1)

    print(f"\n# Bootstrap-mean K_full(2^{log2N}) = {mean:.5f}  ± {sd:.5f} (single-seed SD)")
    print(f"# SE of mean (5 seeds): {se:.5f}")
    print(f"# Gap K_h - K_full: {K_H_USER - mean:.5f}")
    print(f"# ξ_X(2^{log2N}) USER: {xi_user_mean:+.4f} ± {xi_user_sd:.4f}")
    print(f"# ξ_X(2^{log2N}) BODY: {xi_body_mean:+.4f} ± {xi_body_sd:.4f}")

    # Compare to candidate predictions
    pred_plateau = 10.268; pred_power = 10.279; pred_logpow = 10.277
    print(f"\n# 2^{log2N} K_full vs predictions from 6-point fits at 2^25..2^34:")
    print(f"#   {'model':<12}  {'pred':>8}  {'empir':>8}  {'gap':>9}  {'σ-dist':>8}")
    for label, pred in [('plateau', pred_plateau), ('power-law', pred_power),
                         ('log-power', pred_logpow)]:
        print(f"#   {label:<12}  {pred:>8.4f}  {mean:>8.4f}  {mean-pred:>+8.4f}  "
              f"{(mean-pred)/se:>+7.2f}")
    print(f"# (one out-of-sample point — narrows weight on candidates, not validation)")

    # Refit ALL THREE models on 7 points to update parameters
    log2N_arr = np.array([25, 27, 28, 30, 32, 34, 36])
    K_arr_full = np.array([9.98124, 10.11622, 10.17708, 10.22815, 10.23989, 10.26312, mean])
    log_N = log2N_arr * np.log(2)
    gap_arr = K_H_USER - K_arr_full
    log_gap = np.log(gap_arr)
    ss_tot = float(np.sum((K_arr_full - K_arr_full.mean())**2))
    ss_tot_log = float(np.sum((log_gap - log_gap.mean())**2))

    # Power-law: gap = A · N^(-α)
    xc = log_N - log_N.mean(); yc = log_gap - log_gap.mean()
    slope_pl = float((xc*yc).sum() / (xc*xc).sum())
    alpha_pl = -slope_pl
    intercept_pl = log_gap.mean() - slope_pl * log_N.mean()
    A_pl = float(np.exp(intercept_pl))
    pred_pl = K_H_USER - A_pl * np.exp(-alpha_pl * log_N)
    R2_pl = 1 - float(np.sum((K_arr_full - pred_pl)**2) / ss_tot)

    # Log-power: gap = A · (log N)^(-α)
    log_logN = np.log(log_N)
    xc = log_logN - log_logN.mean(); yc = log_gap - log_gap.mean()
    slope_lg = float((xc*yc).sum() / (xc*xc).sum())
    alpha_lg = -slope_lg
    intercept_lg = log_gap.mean() - slope_lg * log_logN.mean()
    A_lg = float(np.exp(intercept_lg))
    pred_lg = K_H_USER - A_lg * np.exp(-alpha_lg * log_logN)
    R2_lg = 1 - float(np.sum((K_arr_full - pred_lg)**2) / ss_tot)

    # Plateau: K = K_inf - a · N^(-β)
    from scipy.optimize import minimize
    def plateau_loss(params):
        K_inf, a, beta = params
        if beta <= 0 or a <= 0 or K_inf > K_H_USER + 0.1: return 1e9
        pred = K_inf - a * np.exp(-beta * log_N)
        return float(np.sum((K_arr_full - pred)**2))
    res = minimize(plateau_loss, x0=[10.275, 1500, 0.495], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 20000})
    K_inf, a_pl, beta_pl = res.x
    pred_plateau_arr = K_inf - a_pl * np.exp(-beta_pl * log_N)
    R2_plateau = 1 - float(np.sum((K_arr_full - pred_plateau_arr)**2) / ss_tot)

    print(f"\n# Refits on 7 points (was 6 in exp 52):")
    print(f"#   {'model':<14}  {'param updates':<55}  {'R²':>7}  {'resid SD':>10}")
    print(f"#   {'power-law':<14}  α: 0.1534 → {alpha_pl:.4f}, A: 5.544 → {A_pl:.4f}".ljust(75)
          + f"  {R2_pl:>7.4f}  {np.std(K_arr_full - pred_pl, ddof=1):>10.5f}")
    print(f"#   {'log-power':<14}  α: 3.158 → {alpha_lg:.4f}, A: 3261 → {A_lg:.4f}".ljust(75)
          + f"  {R2_lg:>7.4f}  {np.std(K_arr_full - pred_lg, ddof=1):>10.5f}")
    print(f"#   {'plateau':<14}  K_∞: 10.275 → {K_inf:.4f}, β: 0.495 → {beta_pl:.4f}".ljust(75)
          + f"  {R2_plateau:>7.4f}  {np.std(K_arr_full - pred_plateau_arr, ddof=1):>10.5f}")
    print(f"# K_h reference: {K_H_USER:.4f}; gap K_h - K_∞ (7-pt): {K_H_USER - K_inf:+.4f}")
    print(f"# β change vs exact 1/2: |β - 0.5| = {abs(beta_pl - 0.5):.4f}")

    # Per-point residuals (all three)
    print(f"\n# Per-point fit residuals:")
    print(f"#   {'log2N':>6}  {'empir':>9}  {'res_PL':>+8}  {'res_LG':>+8}  {'res_plat':>+9}")
    for i, log2N_i in enumerate(log2N_arr):
        print(f"#   {log2N_i:>6}  {K_arr_full[i]:>9.5f}  "
              f"{(K_arr_full - pred_pl)[i]:>+8.5f}  "
              f"{(K_arr_full - pred_lg)[i]:>+8.5f}  "
              f"{(K_arr_full - pred_plateau_arr)[i]:>+9.5f}")

    # ξ_X trajectory recap
    print(f"\n# ξ_X(N) trajectory across N (USER frame):")
    xi_traj = {25: -0.2575, 27: -0.2672, 28: -0.3040, 30: -0.2145, 32: -0.1226,
               34: -0.0993, 36: xi_user_mean}
    for k in sorted(xi_traj.keys()):
        marker = "  ← new" if k == log2N else ""
        print(f"#   2^{k}: ξ_user = {xi_traj[k]:+.4f}{marker}")

    # Save
    rows = []
    for i, seed in enumerate(seeds):
        rows.append({
            'log2N': log2N, 'N': N, 'seed': seed,
            'K_full': float(K_seeds[i]),
            'xi_user': float(xi_user_seeds[i]),
            'xi_body': float(xi_body_seeds[i]),
            'K_full_mean': mean, 'K_full_sd': sd,
            'xi_user_mean': xi_user_mean, 'xi_body_mean': xi_body_mean,
            'K_inf_7pt': K_inf, 'beta_7pt': beta_pl, 'R2_plateau_7pt': R2_plateau,
            'alpha_pl_7pt': alpha_pl, 'R2_pl_7pt': R2_pl,
            'alpha_lg_7pt': alpha_lg, 'R2_lg_7pt': R2_lg,
        })
    out_csv = out_dir / "53_K_eff_at_2_36_bootstrap.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
