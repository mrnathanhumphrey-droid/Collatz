"""
Experiment 51 — K_full and ξ_X at N=2^34 with 5-seed bootstrap, to
discriminate three candidate convergence laws for K_h - K_full(N):

  (1) Logarithmic: gap ~ C/log(N)            → 32→34 Δ ~ +0.005
  (2) Stretched exp: gap ~ C·exp(-γ·log(N)^β) → 32→34 Δ ~ +0.0003 (10× smaller)
  (3) Plateau: K_full → K_∞ < K_h            → 32→34 Δ ~ +0.001 (matches 30→32)

Also reports ξ_X to test whether the Weibull→Gumbel drift continues:
  - ξ_X stays at -0.12 or grows back negative → finite-N rotation
  - ξ_X keeps moving toward 0 → asymptotic Gumbel limit

Reference values from prior runs:
  K_full(2^30) = 10.2349  (gap 0.1933)
  K_full(2^32) = 10.2399 ± 0.008  (gap 0.1883)
  ξ_X USER:    2^30 = -0.21,  2^32 = -0.12
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
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0
        next_idx = 0
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


def analyze_one_seed(N, seed, n_starts):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
    thresh_sorted, sort_idx, raw_thresh = build_thresholds(starts)
    sigma_arr, s_arr, ok_arr = walk_orbit_with_passage(
        starts, thresh_sorted, MAX_VAL, 200_000)
    n_arr = starts[ok_arr]
    sigma_arr_ok = sigma_arr[ok_arr]
    s_arr_ok = s_arr[ok_arr]
    sort_idx_ok = sort_idx[ok_arr]
    raw_thresh_ok = raw_thresh[ok_arr]
    inv_sort = np.argsort(sort_idx_ok, axis=1)
    s_phys = np.zeros((len(n_arr), 4), dtype=np.int32)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(s_arr_ok, rev[:, None], axis=1).flatten()

    K_full = K_eff_slope(sigma_arr_ok.astype(np.float64), s_phys, raw_thresh_ok)

    edges = np.percentile(sigma_arr_ok, [25, 50, 75, 90, 95, 99])
    masks = [
        sigma_arr_ok <= edges[0],
        (sigma_arr_ok > edges[0]) & (sigma_arr_ok <= edges[1]),
        (sigma_arr_ok > edges[1]) & (sigma_arr_ok <= edges[2]),
        (sigma_arr_ok > edges[2]) & (sigma_arr_ok <= edges[3]),
        (sigma_arr_ok > edges[3]) & (sigma_arr_ok <= edges[4]),
        (sigma_arr_ok > edges[4]) & (sigma_arr_ok <= edges[5]),
        sigma_arr_ok > edges[5],
    ]
    q_arr = np.array([0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995])
    z_q = np.array([float(norm.ppf(q)) for q in q_arr])
    K_emp = np.array([K_eff_slope(sigma_arr_ok[m].astype(np.float64),
                                   s_phys[m], raw_thresh_ok[m]) for m in masks])

    body = (q_arr >= 0.30) & (q_arr <= 0.70)
    z_b = z_q[body]; K_b = K_emp[body]
    z_c = z_b - z_b.mean(); K_c = K_b - K_b.mean()
    b_fit = float((z_c * K_c).sum() / (z_c * z_c).sum())
    K_h_fit = float(K_b.mean() - b_fit * z_b.mean())

    X_thresh = norm.ppf(0.75)
    tail = q_arr >= 0.75
    q_tail = q_arr[tail]
    excess_user = (K_emp[tail] - K_H_USER) / B_USER - X_thresh
    excess_body = (K_emp[tail] - K_h_fit) / b_fit - X_thresh
    xi_user, sig_user = fit_xi_NLS(0.75, q_tail, excess_user)
    xi_body, sig_body = fit_xi_NLS(0.75, q_tail, excess_body)

    return {
        'seed': seed, 'K_full': K_full,
        'K_h_fit': K_h_fit, 'b_fit': b_fit,
        'xi_user': xi_user, 'xi_body': xi_body,
        'q_arr': q_arr, 'K_emp': K_emp,
        'n_completed': int(ok_arr.sum()),
    }


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N = 34
    N = 1 << log2N
    n_starts = 500_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# N = 2^{log2N} = {N:,}")
    print(f"# Sampling {n_starts:,} odd starts per seed, {len(seeds)} seeds.")

    results = []
    for seed in seeds:
        t0 = time.perf_counter()
        r = analyze_one_seed(N, seed, n_starts)
        elapsed = time.perf_counter() - t0
        print(f"  seed={seed:>5}: K_full={r['K_full']:.5f}  ξ_user={r['xi_user']:+.4f}  "
              f"ξ_body={r['xi_body']:+.4f}  K_h_fit={r['K_h_fit']:.4f}  "
              f"b_fit={r['b_fit']:.4f}  ({elapsed:.1f}s, {r['n_completed']:,} ok)")
        results.append(r)

    K_arr = np.array([r['K_full'] for r in results])
    xi_user_arr = np.array([r['xi_user'] for r in results])
    xi_body_arr = np.array([r['xi_body'] for r in results])
    b_fit_arr = np.array([r['b_fit'] for r in results])
    K_h_fit_arr = np.array([r['K_h_fit'] for r in results])

    K_full_mean = K_arr.mean()
    K_full_sd = K_arr.std(ddof=1)
    print(f"\n# Mean K_full(2^{log2N}) = {K_full_mean:.5f}  ± {K_full_sd:.5f} (single-seed SD)")
    print(f"# K_h - mean K_full      = {K_H_USER - K_full_mean:.5f}")
    print(f"# Reference K_full(2^32) = 10.24011 (5-seed mean), gap 0.18829")
    print(f"# Reference K_full(2^30) = 10.2349, gap 0.1933")

    # Per-octave deltas
    K_30 = 10.2349
    K_32 = 10.24011
    K_34 = K_full_mean
    delta_30_32 = (K_32 - K_30) / 2  # per-octave (over 2 octaves)
    delta_32_34 = (K_34 - K_32) / 2
    print(f"\n# Per-octave Δ K_full:")
    print(f"#   28→30: +0.0250")
    print(f"#   30→32: {delta_30_32:+.5f}")
    print(f"#   32→34: {delta_32_34:+.5f}")

    # Discriminate candidates
    print(f"\n# Candidate predictions (per-octave Δ at 32→34):")
    print(f"#   (1) Logarithmic       ~ +0.005")
    print(f"#   (2) Stretched-exp     ~ +0.0003 (10× smaller than 30→32)")
    print(f"#   (3) Plateau           ~ +0.001 (similar to 30→32)")
    print(f"#   Empirical 32→34       = {delta_32_34:+.5f}")
    if delta_32_34 < 0.0006:
        verdict = "(2) STRETCHED-EXPONENTIAL favored"
    elif delta_32_34 < 0.002:
        verdict = "(3) PLATEAU favored"
    elif delta_32_34 < 0.008:
        verdict = "(1) LOGARITHMIC favored"
    else:
        verdict = "Unexpected — slow climb actually accelerated"
    print(f"#   → {verdict}")

    # ξ_X trajectory
    xi_user_mean = xi_user_arr.mean()
    xi_body_mean = xi_body_arr.mean()
    print(f"\n# ξ_X(2^{log2N}):")
    print(f"#   ξ_user mean = {xi_user_mean:+.4f}  ± {xi_user_arr.std(ddof=1):.4f}")
    print(f"#   ξ_body mean = {xi_body_mean:+.4f}  ± {xi_body_arr.std(ddof=1):.4f}")
    print(f"#   Reference:")
    print(f"#     ξ_user(2^30) = -0.2145, ξ_user(2^32) = -0.1226")
    print(f"#     ξ_body(2^30) = -0.2722, ξ_body(2^32) = -0.1900")
    if xi_user_mean > -0.07:
        xi_verdict = "ξ_X continuing toward 0 — asymptotic Gumbel rotation real"
    elif xi_user_mean > -0.15:
        xi_verdict = "ξ_X still drifting toward 0 (consistent with prior trajectory)"
    elif xi_user_mean > -0.20:
        xi_verdict = "ξ_X moving back toward negative — possible finite-N reversal"
    else:
        xi_verdict = "ξ_X pulled back hard toward negative — Weibull regime restored"
    print(f"#   → {xi_verdict}")

    # Save CSV
    rows = []
    for r in results:
        for i, q in enumerate(r['q_arr']):
            rows.append({
                'log2N': log2N, 'N': N, 'seed': r['seed'],
                'q': float(q), 'K_emp': float(r['K_emp'][i]),
                'K_full': r['K_full'], 'K_h_fit': r['K_h_fit'],
                'b_fit': r['b_fit'],
                'xi_user': r['xi_user'], 'xi_body': r['xi_body'],
            })
    out_csv = out_dir / "51_K_eff_at_2_34_bootstrap.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
