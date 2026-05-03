"""
Experiment 49 — extended N-sweep of K_eff_band(q) and K_eff_full to
discriminate slow-decay vs bounded-oscillation, plus limit-cycle hypothesis
at higher resolution.

Adds N ∈ {2^21, 2^26, 2^28, 2^30} to the existing {2^20, 2^22, 2^23, 2^24,
2^25, 2^27}, giving 10 N values for phase analysis.

Tests:
  1. Slow-decay vs bounded-oscillation: does post-2^25 rise continue
     monotonically (slow decay; converging to K_h) or reverse (bounded
     oscillation around some asymptote)? 2^28 and 2^30 datapoints decisive.
  2. Limit-cycle hypothesis: at 8-10 point resolution, does K_eff(N) cluster
     by log2(N) mod log2(3) classes (the Diophantine residue indexing
     conjectured but ruled out at 6 points in Result 9)?
  3. ξ_X(N) trend across 10 N's: does the GPD shape parameter stabilize
     or continue drifting?
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


@njit(parallel=True, cache=True)
def first_passage_4thresh(starts, thresholds_per_n):
    n = len(starts)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    max_steps = 20000
    for i in prange(n):
        m = starts[i]
        steps = 0
        next_idx = 0
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0
            next_idx += 1
        while next_idx < 4 and steps < max_steps:
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            steps += 1
            while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                s_arr[i, next_idx] = steps
                next_idx += 1
    return s_arr


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


def analyze_N(N_val, data_dir, n_sample=500_000, rng=None):
    log2N = int(np.log2(N_val))
    print(f"\n[N=2^{log2N}] loading parquet ...", end=" ", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{N_val}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int32)
    print(f"{len(n_arr):,} odd rows", flush=True)

    if len(n_arr) > n_sample:
        idx = rng.choice(len(n_arr), size=n_sample, replace=False)
        n_arr = n_arr[idx]; sigma_arr = sigma_arr[idx]

    thresh_sorted, sort_idx, raw_thresh = build_thresholds(n_arr)
    s_arr = first_passage_4thresh(n_arr, thresh_sorted)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros((len(n_arr), 4), dtype=np.int32)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()

    # Full K_eff
    K_full = K_eff_slope(sigma_arr.astype(np.float64), s_phys, raw_thresh)

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
                                   s_phys[m], raw_thresh[m]) for m in masks])

    # Body fit on q ∈ {0.375, 0.625}
    body = (q_arr >= 0.30) & (q_arr <= 0.70)
    z_b = z_q[body]; K_b = K_emp[body]
    z_c = z_b - z_b.mean(); K_c = K_b - K_b.mean()
    b_fit = float((z_c * K_c).sum() / (z_c * z_c).sum())
    K_h_fit = float(K_b.mean() - b_fit * z_b.mean())

    # Tail GPD fits (4 points, q ≥ 0.825)
    X_thresh = norm.ppf(0.75)
    tail = q_arr >= 0.75
    q_tail = q_arr[tail]
    excess_user = (K_emp[tail] - K_H_USER) / B_USER - X_thresh
    excess_body = (K_emp[tail] - K_h_fit) / b_fit - X_thresh
    xi_user, sig_user = fit_xi_NLS(0.75, q_tail, excess_user)
    xi_body, sig_body = fit_xi_NLS(0.75, q_tail, excess_body)

    print(f"  K_full={K_full:.4f}  K_h_fit={K_h_fit:.4f}  b_fit={b_fit:.4f}  "
          f"ξ_user={xi_user:+.3f}  ξ_body={xi_body:+.3f}")

    return {
        'N': N_val, 'log2N': log2N,
        'K_full': K_full, 'K_h_fit': K_h_fit, 'b_fit': b_fit,
        'xi_user': xi_user, 'xi_body': xi_body,
        'sig_user': sig_user, 'sig_body': sig_body,
        'q_arr': q_arr, 'K_emp': K_emp, 'z_q': z_q,
    }


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    rng = np.random.default_rng(42)
    Ns_log2 = [20, 21, 22, 23, 24, 25, 26, 27, 28, 30]
    results = []
    for log2N in Ns_log2:
        N = 1 << log2N
        r = analyze_N(N, data_dir, n_sample=500_000, rng=rng)
        results.append(r)

    # Phase analysis: K_full vs N
    print(f"\n{'='*70}")
    print(f"K_full(N) phase analysis")
    print(f"{'='*70}")
    print(f"  {'log2N':>6}  {'K_full':>9}  {'Δ from prev':>12}  "
          f"{'mod log2(3)':>13}")
    log2_3 = np.log2(3.0)
    for i, r in enumerate(results):
        log2N = r['log2N']
        K = r['K_full']
        dprev = K - results[i-1]['K_full'] if i > 0 else 0.0
        mod_val = log2N % log2_3
        print(f"  {log2N:>6}  {K:>9.4f}  {dprev:>+12.4f}  {mod_val:>13.4f}")

    K_arr = np.array([r['K_full'] for r in results])
    log2N_arr = np.array([r['log2N'] for r in results])

    # Slow-decay vs oscillation: monotone after 2^25?
    print(f"\n  After-2^25 trajectory:")
    for i, r in enumerate(results):
        if r['log2N'] >= 25:
            print(f"    2^{r['log2N']}: K_full = {r['K_full']:.4f}")
    diffs = np.diff([r['K_full'] for r in results if r['log2N'] >= 25])
    sign_changes = np.sum(np.diff(np.sign(diffs)) != 0)
    monotone = sign_changes == 0
    print(f"  Sign changes in ΔK after 2^25: {sign_changes}")
    print(f"  → {'MONOTONE (slow-decay)' if monotone else 'NON-MONOTONE (oscillation)'}")

    # Limit-cycle test: K_full vs log2N mod log2(3)
    mod_arr = log2N_arr % log2_3
    sin_arg = 2 * np.pi * mod_arr / log2_3
    cos_arg = sin_arg
    A = np.column_stack([np.ones_like(K_arr), np.sin(sin_arg), np.cos(cos_arg)])
    coef, _, _, _ = np.linalg.lstsq(A, K_arr, rcond=None)
    pred = A @ coef
    ss_res = np.sum((K_arr - pred) ** 2)
    ss_tot = np.sum((K_arr - K_arr.mean()) ** 2)
    R2_lim = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    print(f"\n  Limit-cycle (sinusoid in log2N mod log2(3)) fit on {len(K_arr)} N values:")
    print(f"    coefs (mean, sin, cos): {coef[0]:.4f}, {coef[1]:.4f}, {coef[2]:.4f}")
    print(f"    R² = {R2_lim:.4f}")
    print(f"    Threshold: R²>0.5 to support limit-cycle (Result 9 used 0.25 at 6 points)")
    print(f"    → {'SUPPORTS limit-cycle' if R2_lim > 0.5 else 'limit-cycle NOT supported'}")

    # ξ_X stability across all N
    print(f"\n  ξ_X across {len(results)} N values:")
    print(f"  {'log2N':>6}  {'ξ_user':>9}  {'ξ_body':>9}  {'b_fit':>8}  {'K_h_fit':>9}")
    for r in results:
        print(f"  {r['log2N']:>6}  {r['xi_user']:>+9.4f}  {r['xi_body']:>+9.4f}  "
              f"{r['b_fit']:>8.4f}  {r['K_h_fit']:>9.4f}")
    xi_user_arr = np.array([r['xi_user'] for r in results])
    xi_body_arr = np.array([r['xi_body'] for r in results])
    print(f"  ξ_user spread: {np.ptp(xi_user_arr):.4f}, mean = {xi_user_arr.mean():+.4f}")
    print(f"  ξ_body spread: {np.ptp(xi_body_arr):.4f}, mean = {xi_body_arr.mean():+.4f}")

    # Save CSV
    rows = []
    for r in results:
        for i, q in enumerate(r['q_arr']):
            rows.append({
                'log2N': r['log2N'], 'N': r['N'],
                'q': float(q), 'z_q': float(r['z_q'][i]),
                'K_emp': float(r['K_emp'][i]),
                'K_full': r['K_full'],
                'K_h_fit': r['K_h_fit'], 'b_fit': r['b_fit'],
                'xi_user': r['xi_user'], 'xi_body': r['xi_body'],
            })
    out_csv = out_dir / "49_K_eff_phase_extended.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")

    # Plot
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        axes[0,0].plot(log2N_arr, K_arr, 'ko-', markersize=8)
        axes[0,0].axhline(K_H_USER, color='gray', linestyle=':', label=f'K_h={K_H_USER:.3f}')
        axes[0,0].set_xlabel('log2(N)'); axes[0,0].set_ylabel('K_full')
        axes[0,0].set_title('K_full vs N (phase analysis)')
        axes[0,0].legend(); axes[0,0].grid(alpha=0.3)

        axes[0,1].plot(log2N_arr, xi_user_arr, 'o-', label='ξ (USER frame)')
        axes[0,1].plot(log2N_arr, xi_body_arr, 's-', label='ξ (BODY frame)')
        axes[0,1].axhline(-0.075, color='red', linestyle=':', label='σ-records ξ=-0.075')
        axes[0,1].set_xlabel('log2(N)'); axes[0,1].set_ylabel('GPD ξ')
        axes[0,1].set_title('Tail GPD ξ vs N')
        axes[0,1].legend(fontsize=8); axes[0,1].grid(alpha=0.3)

        # Limit-cycle: K_full vs (log2N mod log2(3))
        order = np.argsort(mod_arr)
        axes[1,0].scatter(mod_arr, K_arr, color='black', s=80)
        for i, ln in enumerate(log2N_arr):
            axes[1,0].annotate(f'2^{ln}', (mod_arr[i], K_arr[i]),
                               fontsize=8, xytext=(3, 3), textcoords='offset points')
        x_fine = np.linspace(0, log2_3, 200)
        sin_f = np.sin(2*np.pi*x_fine/log2_3); cos_f = np.cos(2*np.pi*x_fine/log2_3)
        y_fine = coef[0] + coef[1]*sin_f + coef[2]*cos_f
        axes[1,0].plot(x_fine, y_fine, 'r-', alpha=0.6, label=f'sinusoid R²={R2_lim:.3f}')
        axes[1,0].set_xlabel('log2(N) mod log2(3)')
        axes[1,0].set_ylabel('K_full')
        axes[1,0].set_title('Limit-cycle test')
        axes[1,0].legend(); axes[1,0].grid(alpha=0.3)

        axes[1,1].plot(log2N_arr, [r['b_fit'] for r in results], 'o-', label='b_fit')
        axes[1,1].plot(log2N_arr, [r['K_h_fit'] for r in results], 's-', label='K_h_fit')
        axes[1,1].axhline(B_USER, color='gray', linestyle=':', alpha=0.5, label=f'b_user={B_USER}')
        axes[1,1].axhline(K_H_USER, color='gray', linestyle='--', alpha=0.5, label=f'K_h_user={K_H_USER:.3f}')
        axes[1,1].set_xlabel('log2(N)'); axes[1,1].set_ylabel('value')
        axes[1,1].set_title('Body-fit (K_h, b) drift')
        axes[1,1].legend(fontsize=8); axes[1,1].grid(alpha=0.3)

        plt.tight_layout()
        plot_path = out_dir / "49_K_eff_phase_extended.png"
        plt.savefig(plot_path, dpi=120)
        print(f"[save] {plot_path}")
    except Exception as e:
        print(f"  (plot skipped: {e})")


if __name__ == "__main__":
    main()
