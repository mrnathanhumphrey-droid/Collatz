"""
Experiment 50 — push K_eff_band(q) and K_full to N=2^32 via direct orbit
sampling (avoids the 64 GB memory wall of full enumeration at 2^32).

Method: sample 500K random odd starts in [3, 2^32]; walk each orbit to 1;
record σ, first-passage steps to 4 thresholds (N^(2/3), √n·log n, √n,
√n/log n). Identical methodology to exp 49 except orbits are walked
directly rather than read from a precomputed parquet.

Tests: confirm post-2^25 power-law convergence gap = 31.9·N^(-0.247),
extend to 6 post-minimum points (2^25..2^32). Also extend ξ_X drift
trajectory.
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
B_USER = 2.275
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_orbit_with_passage(starts, thresholds_per_n, max_value, max_steps):
    """For each start, walk Collatz to 1, record σ and first-passage steps to
    4 thresholds (in DESCENDING order of threshold value).

    Returns:
      sigma_arr (int32, n)         total stopping time
      s_arr (int32, n×4)           steps to first hit threshold[k]
      ok_arr (bool, n)             True if orbit completed within max_steps
    """
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0
        next_idx = 0
        # Initial check: starts always above all thresholds (m > N/2 typically),
        # so don't pre-mark; however protect against m already ≤ thresh[0]:
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0
            next_idx += 1
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
                s_arr[i, next_idx] = steps
                next_idx += 1
        if not failed and m == 1:
            sigma_arr[i] = steps
            ok_arr[i] = True
    return sigma_arr, s_arr, ok_arr


def build_thresholds(n_arr):
    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)
    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)
    thresh_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    return thresh_sorted, sort_idx, raw


def K_eff_slope(sigma, s_per_thresh, thresh_phys):
    log_f = np.log(np.maximum(thresh_phys, 1.0))
    R = sigma[:, None] - s_per_thresh.astype(np.float64)
    log_f_mean = log_f.mean(axis=0)
    R_mean = R.mean(axis=0)
    x = log_f_mean - log_f_mean.mean()
    y = R_mean - R_mean.mean()
    return float((x * y).sum() / (x * x).sum())


def gpd_factor(p, xi):
    if abs(xi) < 1e-9:
        return -np.log(1.0 - p)
    return ((1.0 - p) ** (-xi) - 1.0) / xi


def fit_xi_NLS(q_thresh, q_tail, excess):
    p = (q_tail - q_thresh) / (1.0 - q_thresh)
    def loss(params):
        xi, sigma = params
        if sigma <= 0 or abs(xi) > 0.5: return 1e20
        pred = sigma * np.array([gpd_factor(pi, xi) for pi in p])
        return float(np.sum((pred - excess) ** 2))
    res = minimize(loss, x0=[0.0, np.median(excess)], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-10})
    return float(res.x[0]), float(res.x[1])


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    log2N = 32
    N = 1 << log2N
    n_starts = 500_000
    print(f"# N = 2^{log2N} = {N:,}")
    print(f"# Sampling {n_starts:,} odd starts in [3, N], walking to 1.")

    rng = np.random.default_rng(42)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1

    # Build thresholds
    print(f"# building thresholds ...", flush=True)
    thresh_sorted, sort_idx, raw_thresh = build_thresholds(starts)

    # Walk orbits
    print(f"# walking {n_starts:,} orbits ...", flush=True)
    t0 = time.perf_counter()
    sigma_arr, s_arr, ok_arr = walk_orbit_with_passage(
        starts, thresh_sorted, MAX_VAL, 100_000)
    print(f"#   walked in {time.perf_counter()-t0:.1f}s, "
          f"completed={int(ok_arr.sum()):,}/{n_starts:,}, "
          f"failed={int((~ok_arr).sum())}")

    if ok_arr.sum() < n_starts * 0.99:
        print(f"# WARNING: only {ok_arr.sum()}/{n_starts} completed; "
              f"check max_value/max_steps.")

    # Filter to OK orbits, restore physical-threshold order
    n_arr = starts[ok_arr]
    sigma_arr = sigma_arr[ok_arr]
    s_arr_ok = s_arr[ok_arr]
    sort_idx_ok = sort_idx[ok_arr]
    raw_thresh_ok = raw_thresh[ok_arr]
    inv_sort = np.argsort(sort_idx_ok, axis=1)
    s_phys = np.zeros((len(n_arr), 4), dtype=np.int32)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(s_arr_ok, rev[:, None], axis=1).flatten()

    # K_full
    K_full = K_eff_slope(sigma_arr.astype(np.float64), s_phys, raw_thresh_ok)
    print(f"\n# K_full at N=2^{log2N}: {K_full:.4f}")
    print(f"# K_h reference: {K_H_USER:.4f}")
    print(f"# Gap K_h - K_full: {K_H_USER - K_full:+.4f}")

    # Bands at midpoints {0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995}
    edges = np.percentile(sigma_arr, [25, 50, 75, 90, 95, 99])
    masks = [
        sigma_arr <= edges[0],
        (sigma_arr > edges[0]) & (sigma_arr <= edges[1]),
        (sigma_arr > edges[1]) & (sigma_arr <= edges[2]),
        (sigma_arr > edges[2]) & (sigma_arr <= edges[3]),
        (sigma_arr > edges[3]) & (sigma_arr <= edges[4]),
        (sigma_arr > edges[4]) & (sigma_arr <= edges[5]),
        sigma_arr > edges[5],
    ]
    q_arr = np.array([0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995])
    z_q = np.array([float(norm.ppf(q)) for q in q_arr])
    K_emp = np.array([K_eff_slope(sigma_arr[m].astype(np.float64),
                                   s_phys[m], raw_thresh_ok[m]) for m in masks])

    # Body fit q ∈ {0.375, 0.625}
    body = (q_arr >= 0.30) & (q_arr <= 0.70)
    z_b = z_q[body]; K_b = K_emp[body]
    z_c = z_b - z_b.mean(); K_c = K_b - K_b.mean()
    b_fit = float((z_c * K_c).sum() / (z_c * z_c).sum())
    K_h_fit = float(K_b.mean() - b_fit * z_b.mean())

    # Tail GPD fits (4 points)
    X_thresh = norm.ppf(0.75)
    tail = q_arr >= 0.75
    q_tail = q_arr[tail]
    excess_user = (K_emp[tail] - K_H_USER) / B_USER - X_thresh
    excess_body = (K_emp[tail] - K_h_fit) / b_fit - X_thresh
    xi_user, sig_user = fit_xi_NLS(0.75, q_tail, excess_user)
    xi_body, sig_body = fit_xi_NLS(0.75, q_tail, excess_body)

    print(f"\n# Per-band K_eff at N=2^{log2N}:")
    print(f"#   {'q':>6} {'z_q':>7} {'K_emp':>9}")
    for i in range(len(q_arr)):
        print(f"#   {q_arr[i]:>6.3f} {z_q[i]:>+7.3f} {K_emp[i]:>9.4f}")
    print(f"\n# Body fit: K_h_fit = {K_h_fit:.4f}, b_fit = {b_fit:.4f}")
    print(f"# Tail GPD: ξ_USER = {xi_user:+.4f}, ξ_BODY = {xi_body:+.4f}")

    # Combine with previous results from exp 49 for slow-decay extrapolation test
    prev = {
        25: 9.9906, 26: 10.0434, 27: 10.1144, 28: 10.1849, 30: 10.2349,
    }
    prev_xi_user = {25: -0.2575, 26: -0.2884, 27: -0.2672, 28: -0.3040, 30: -0.2145}
    prev_xi_body = {25: -0.4270, 26: -0.4399, 27: -0.3712, 28: -0.3883, 30: -0.2722}
    prev_b = {25: 1.8604, 26: 2.0315, 27: 2.3690, 28: 2.4419, 30: 2.6138}
    prev_Kh = {25: 9.7837, 26: 9.7001, 27: 9.7213, 28: 9.7774, 30: 9.8072}

    # Add 2^32
    Kf_all = {**prev, log2N: K_full}
    xi_user_all = {**prev_xi_user, log2N: xi_user}
    xi_body_all = {**prev_xi_body, log2N: xi_body}
    b_all = {**prev_b, log2N: b_fit}
    Kh_all = {**prev_Kh, log2N: K_h_fit}

    log2_arr = np.array(sorted(Kf_all.keys()))
    K_arr = np.array([Kf_all[k] for k in log2_arr])
    gap_arr = K_H_USER - K_arr
    print(f"\n# Slow-decay branch (post-minimum) with 2^{log2N} added:")
    print(f"#   {'log2N':>6} {'K_full':>9} {'gap':>8}")
    for i, ln in enumerate(log2_arr):
        print(f"#   {ln:>6} {K_arr[i]:>9.4f} {gap_arr[i]:>8.4f}")

    # Re-fit power law on extended data
    logN = log2_arr * np.log(2)
    log_gap = np.log(gap_arr)
    xc = logN - logN.mean(); yc = log_gap - log_gap.mean()
    slope = float((xc*yc).sum() / (xc*xc).sum())
    intercept = float(log_gap.mean() - slope*logN.mean())
    A = float(np.exp(intercept))
    alpha = -slope
    pred_log = intercept + slope*logN
    ss_res = float(np.sum((log_gap - pred_log)**2))
    ss_tot = float(np.sum((log_gap - log_gap.mean())**2))
    R2 = 1 - ss_res/ss_tot

    print(f"\n# Power law fit gap = A · N^(-α) on {len(log2_arr)} post-minimum points:")
    print(f"#   α = {alpha:.4f}  (5-point fit gave 0.2473)")
    print(f"#   A = {A:.4f}     (5-point gave 31.90)")
    print(f"#   R² = {R2:.4f}   (5-point gave 0.980)")
    print(f"#   Predicted at 2^{log2N}: gap = {A*np.exp(slope*log2N*np.log(2)):.4f}, "
          f"K_full = {K_H_USER - A*np.exp(slope*log2N*np.log(2)):.4f}")
    print(f"#   Actual at 2^{log2N}:    gap = {gap_arr[-1]:.4f}, "
          f"K_full = {K_arr[-1]:.4f}")

    # Discrimination test: is α closer to 0.247 (drift) or 0.250 (rational 1/4)?
    # Restrict to N ≥ 2^27 (cleaner asymptotic regime)
    mask_late = log2_arr >= 27
    if mask_late.sum() >= 2:
        logN_late = logN[mask_late]; log_gap_late = log_gap[mask_late]
        xc_l = logN_late - logN_late.mean(); yc_l = log_gap_late - log_gap_late.mean()
        slope_late = float((xc_l*yc_l).sum() / (xc_l*xc_l).sum())
        alpha_late = -slope_late
        print(f"\n# Late-N (log2N ≥ 27) sub-fit α = {alpha_late:.4f}")
        print(f"# Distance to α=1/4: {abs(alpha_late - 0.25):.4f}")

    # ξ_X trend with 2^32 included
    print(f"\n# ξ_X trajectory (post-minimum branch):")
    print(f"#   {'log2N':>6} {'ξ_USER':>9} {'ξ_BODY':>9} {'b_fit':>8} {'K_h_fit':>9}")
    for ln in log2_arr:
        print(f"#   {ln:>6} {xi_user_all[ln]:>+9.4f} {xi_body_all[ln]:>+9.4f} "
              f"{b_all[ln]:>8.4f} {Kh_all[ln]:>9.4f}")

    # Save row
    rows = []
    for i, q in enumerate(q_arr):
        rows.append({
            'log2N': log2N, 'N': N,
            'q': float(q), 'z_q': float(z_q[i]),
            'K_emp': float(K_emp[i]),
            'K_full': K_full, 'K_h_fit': K_h_fit, 'b_fit': b_fit,
            'xi_user': xi_user, 'xi_body': xi_body,
            'alpha_powerlaw_6pt': alpha, 'A_powerlaw_6pt': A,
            'R2_powerlaw_6pt': R2,
        })
    out_csv = out_dir / "50_K_eff_at_2_32.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
