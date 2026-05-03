"""
Experiment 59 — Push K_eff_band(q=0.875) and ξ_X to N ∈ {2^38, 2^40, 2^42}.

Tests whether the q=0.875 band reversal (peak at 2^30, decreasing through 2^36)
stabilizes at ~14.5 (finite-N transition) or continues decreasing structurally.

Per N ∈ {2^25, 2^27, 2^28, 2^30, 2^32, 2^34, 2^36, 2^38, 2^40, 2^42} × 5 seeds × 500K orbits:
  - K_eff_band(q=0.875): top-quartile K_eff regression slope on log(threshold)
  - ξ_X (USER frame): GPD shape fit to standardized K_emp residuals at q ≥ 0.75
  - E[v]_q875: mean v_per_orbit (n_even/n_odd) over top-quartile orbits
  - Plus K_full, K_eff_band(q=0.125, 0.375, 0.625, 0.975) for context

Decisive outcomes:
  (a) K_eff_band(q=0.875) stabilizes at 14.5 ± 0.1 across 2^38-2^42 → finite-N
      Weibull-to-Gumbel transition; v3.4 K_full → K_h estimate holds.
  (b) Continues decreasing → structural sub-K_h asymptote; revise constant 4.
  (c) Mixed signals → identify what's still moving.

Compute scaling: 500K orbit walks at 2^42 in ~0.2s/seed (numba parallel),
so total ~10s for full sweep. Brief's "1-3 hours/N" estimate was off by ~10⁴×.
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
LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
MAX_VAL = np.int64(2**62)


def K_of_v(v):
    return (1 + v) / (v * LOG_2 - LOG_3)


@njit(parallel=True, cache=True)
def walk_orbit(starts, thresholds_per_n, max_value, max_steps):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)
    n_odd = np.zeros(n, dtype=np.int32)
    n_even = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0; next_idx = 0
        odd_count = 0; even_count = 0
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0; next_idx += 1
        failed = False
        while m != 1 and steps < max_steps:
            if m & 1:
                if m > max_value // 3:
                    failed = True; break
                m = 3 * m + 1
                odd_count += 1
            else:
                m = m >> 1
                even_count += 1
            steps += 1
            while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                s_arr[i, next_idx] = steps; next_idx += 1
        if not failed and m == 1:
            sigma_arr[i] = steps
            ok_arr[i] = True
            n_odd[i] = odd_count
            n_even[i] = even_count
    return sigma_arr, s_arr, ok_arr, n_odd, n_even


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
    p = (q_tail - q_thresh) / (1.0 - q_thresh)
    def loss(params):
        xi, sigma = params
        if sigma <= 0 or abs(xi) > 0.5: return 1e20
        pred = sigma * np.array([gpd_factor(pi, xi) for pi in p])
        return float(np.sum((pred - excess) ** 2))
    res = minimize(loss, x0=[0.0, np.median(np.abs(excess))], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-10})
    return float(res.x[0]), float(res.x[1])


def analyze_one_seed(N, seed, n_starts, max_steps=1_000_000):
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
    sigma, s, ok, n_odd, n_even = walk_orbit(starts, th, MAX_VAL, max_steps)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros((len(starts), 4), dtype=np.int32)
    for col in range(4):
        s_phys[:, col] = np.take_along_axis(s, inv_sort[:, col:col+1], axis=1).flatten()
    sigma_ok = sigma[ok].astype(np.float64)
    s_ok = s_phys[ok]
    raw_ok = raw[ok]
    n_odd_ok = n_odd[ok].astype(np.float64)
    n_even_ok = n_even[ok].astype(np.float64)
    v_per_orbit = n_even_ok / np.maximum(n_odd_ok, 1)

    # K_full
    K_full = K_eff_slope(sigma_ok, s_ok, raw_ok)

    # 5 user-midpoint bands for K_eff_band trajectory
    edges5 = np.percentile(sigma_ok, [25, 50, 75, 95])
    masks5 = [
        sigma_ok <= edges5[0],
        (sigma_ok > edges5[0]) & (sigma_ok <= edges5[1]),
        (sigma_ok > edges5[1]) & (sigma_ok <= edges5[2]),
        sigma_ok > edges5[2],          # q=0.875 (top quartile)
        sigma_ok > edges5[3],          # q=0.975 (top 5%)
    ]
    q5 = [0.125, 0.375, 0.625, 0.875, 0.975]
    K_band5 = np.array([K_eff_slope(sigma_ok[m], s_ok[m], raw_ok[m]) for m in masks5])

    # E[v]_q875: mean v over top quartile
    mask_top = masks5[3]
    E_v_q875 = float(v_per_orbit[mask_top].mean())
    K_pred_q875_from_E_v = K_of_v(E_v_q875)

    # E[v]_q125 for parity with experiment 56
    mask_low = masks5[0]
    E_v_q125 = float(v_per_orbit[mask_low].mean())

    # ξ_X (USER frame): GPD fit on q ≥ 0.75 standardized residuals
    edges7 = np.percentile(sigma_ok, [25, 50, 75, 90, 95, 99])
    masks7 = [
        sigma_ok <= edges7[0],
        (sigma_ok > edges7[0]) & (sigma_ok <= edges7[1]),
        (sigma_ok > edges7[1]) & (sigma_ok <= edges7[2]),
        (sigma_ok > edges7[2]) & (sigma_ok <= edges7[3]),
        (sigma_ok > edges7[3]) & (sigma_ok <= edges7[4]),
        (sigma_ok > edges7[4]) & (sigma_ok <= edges7[5]),
        sigma_ok > edges7[5],
    ]
    q7 = np.array([0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995])
    K7 = np.array([K_eff_slope(sigma_ok[m], s_ok[m], raw_ok[m]) for m in masks7])
    z_q = np.array([float(norm.ppf(q)) for q in q7])
    X_thresh = norm.ppf(0.75)
    tail = q7 >= 0.75
    excess_user = (K7[tail] - K_H_USER) / 2.275 - X_thresh
    xi_user, _ = fit_xi_NLS(0.75, q7[tail], excess_user)

    return {
        'K_full': K_full,
        'K_q125': float(K_band5[0]),
        'K_q375': float(K_band5[1]),
        'K_q625': float(K_band5[2]),
        'K_q875': float(K_band5[3]),
        'K_q975': float(K_band5[4]),
        'E_v_q875': E_v_q875,
        'K_pred_from_E_v_q875': float(K_pred_q875_from_E_v),
        'E_v_q125': E_v_q125,
        'xi_X_user': xi_user,
        'n_total': int(ok.sum()),
        'n_failed': int((~ok).sum()),
    }


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N_list = [25, 27, 28, 30, 32, 34, 36, 38, 40, 42]
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 500_000

    print(f"# Experiment 59: K_eff_band(q=0.875) + ξ_X + E[v]_q875 to N=2^42")
    print(f"# {len(log2N_list)} N × {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# K_h = {K_H_USER:.4f}")
    print()

    rows = []
    summary = {}
    t_total0 = time.perf_counter()
    for log2N in log2N_list:
        N = 1 << log2N
        seed_results = []
        for seed in seeds:
            t0 = time.perf_counter()
            r = analyze_one_seed(N, seed, n_starts)
            elapsed = time.perf_counter() - t0
            seed_results.append(r)
            rows.append({
                'log2N': log2N, 'N': N, 'seed': seed, **r,
            })
            print(f"  2^{log2N:>2} seed={seed:>5}: K_q875={r['K_q875']:.4f}  ξ_X={r['xi_X_user']:+.4f}  "
                  f"E[v]_q875={r['E_v_q875']:.4f}  K_full={r['K_full']:.4f}  "
                  f"failed={r['n_failed']}  ({elapsed:.1f}s)")
        # Aggregate
        K_q875_arr = np.array([r['K_q875'] for r in seed_results])
        xi_arr = np.array([r['xi_X_user'] for r in seed_results])
        E_v_arr = np.array([r['E_v_q875'] for r in seed_results])
        K_full_arr = np.array([r['K_full'] for r in seed_results])
        K_q125_arr = np.array([r['K_q125'] for r in seed_results])
        K_q975_arr = np.array([r['K_q975'] for r in seed_results])
        summary[log2N] = {
            'K_q875_mean': K_q875_arr.mean(), 'K_q875_sd': K_q875_arr.std(ddof=1),
            'K_q875_se': K_q875_arr.std(ddof=1) / np.sqrt(len(seeds)),
            'xi_mean': xi_arr.mean(), 'xi_sd': xi_arr.std(ddof=1),
            'xi_se': xi_arr.std(ddof=1) / np.sqrt(len(seeds)),
            'E_v_mean': E_v_arr.mean(), 'E_v_sd': E_v_arr.std(ddof=1),
            'K_full_mean': K_full_arr.mean(),
            'K_q125_mean': K_q125_arr.mean(),
            'K_q975_mean': K_q975_arr.mean(),
            'n_failed_total': sum(r['n_failed'] for r in seed_results),
        }
    print(f"\n# Total compute: {time.perf_counter()-t_total0:.1f}s")

    # Trajectory tables
    print(f"\n# K_eff_band(q=0.875) trajectory across N (5-seed bootstrap mean ± SE):")
    print(f"#   {'log2N':>6}  {'K_q875':>8}  {'SE':>7}  {'xi_X':>8}  {'xi_SE':>7}  "
          f"{'E[v]_q875':>10}  {'K_full':>8}  {'failed':>7}")
    for log2N in log2N_list:
        s = summary[log2N]
        print(f"#   {log2N:>6}  {s['K_q875_mean']:>8.4f}  {s['K_q875_se']:>7.4f}  "
              f"{s['xi_mean']:>+8.4f}  {s['xi_se']:>7.4f}  "
              f"{s['E_v_mean']:>10.4f}  {s['K_full_mean']:>8.4f}  {s['n_failed_total']:>7}")

    # Per-octave Δ K_q875
    print(f"\n# Per-octave Δ K_q875 with bootstrap SEs:")
    print(f"#   {'step':>10}  {'Δ K_q875':>10}  {'SE on Δ':>9}  {'z':>6}")
    for i in range(1, len(log2N_list)):
        a = log2N_list[i-1]; b = log2N_list[i]
        d_octave = b - a
        delta = summary[b]['K_q875_mean'] - summary[a]['K_q875_mean']
        se = np.sqrt(summary[a]['K_q875_se']**2 + summary[b]['K_q875_se']**2)
        z = delta / se if se > 0 else 0
        print(f"   {a}→{b}:  {delta:>+10.4f}  {se:>9.4f}  {z:>+6.2f}")

    # Verdict logic
    K_q875_38 = summary[38]['K_q875_mean']
    K_q875_40 = summary[40]['K_q875_mean']
    K_q875_42 = summary[42]['K_q875_mean']
    K_q875_36 = summary[36]['K_q875_mean']
    se_38_42 = max(summary[k]['K_q875_se'] for k in [38, 40, 42])

    delta_38_42 = K_q875_42 - K_q875_38
    se_delta = np.sqrt(summary[38]['K_q875_se']**2 + summary[42]['K_q875_se']**2)

    print(f"\n# Verdict (per brief decisive outcomes):")
    print(f"#   K_q875(2^36) = {K_q875_36:.4f}  (existing data)")
    print(f"#   K_q875(2^38) = {K_q875_38:.4f}  (NEW)")
    print(f"#   K_q875(2^40) = {K_q875_40:.4f}  (NEW)")
    print(f"#   K_q875(2^42) = {K_q875_42:.4f}  (NEW)")
    print(f"#   38→42 Δ = {delta_38_42:+.4f} ± {se_delta:.4f}  (z = {delta_38_42/se_delta:+.2f})")

    # Check (a): stable in 14.4-14.6 band across 38-42, |Δ| < 2σ
    stable_38_42 = abs(delta_38_42) < 2 * se_delta
    plateau_band = all(14.4 <= summary[k]['K_q875_mean'] <= 14.6 for k in [38, 40, 42])

    # Check (b): monotone-decreasing 38→42 with |Δ| > 3σ
    decreasing = (delta_38_42 < 0) and (abs(delta_38_42) > 3 * se_delta)

    if stable_38_42 and plateau_band:
        print(f"#   VERDICT (a): K_q875 stabilizes at ~14.5 across 2^38-2^42")
        print(f"#                Reversal is finite-N Weibull-to-Gumbel transition")
        print(f"#                v3.4 K_full → K_h estimate holds")
    elif decreasing:
        print(f"#   VERDICT (b): K_q875 continues decreasing through 2^42")
        print(f"#                Reversal is structural; band asymptote below 14.5")
        print(f"#                Constant 4 framing reopens — K_full → K_∞ < K_h")
    else:
        print(f"#   VERDICT (c): mixed signals; document trajectory shape")

    # ξ_X trajectory
    xi_36 = summary[36]['xi_mean']
    xi_42 = summary[42]['xi_mean']
    print(f"\n# ξ_X trajectory:")
    print(f"#   ξ_X(2^36) = {xi_36:+.4f}, ξ_X(2^42) = {xi_42:+.4f}, change = {xi_42-xi_36:+.4f}")
    print(f"#   Linear extrapolation: ξ_X = 0 at log2N ≈ ", end="")
    # Fit linear ξ_X vs log2N over 2^36..2^42
    log2_pts = np.array([36, 38, 40, 42])
    xi_pts = np.array([summary[k]['xi_mean'] for k in log2_pts])
    p = np.polyfit(log2_pts, xi_pts, 1)
    if p[0] > 1e-9:
        crossing_log2N = -p[1] / p[0]
        print(f"{crossing_log2N:.1f}  (slope {p[0]:+.4f}/octave)")
    else:
        print(f"never (slope {p[0]:+.4f} non-positive)")

    # E[v]_q875 trajectory
    print(f"\n# E[v]_q875 trajectory (does upper band have its own structural deviation?):")
    print(f"#   {'log2N':>6}  {'E[v]_q875':>10}")
    for log2N in log2N_list:
        print(f"#   {log2N:>6}  {summary[log2N]['E_v_mean']:>10.4f}")

    # Save CSV
    out_csv = out_dir / "59_q875_extension.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")

    # Save summary CSV
    summary_rows = []
    for log2N in log2N_list:
        s = summary[log2N]
        summary_rows.append({
            'log2N': log2N, 'N': 1 << log2N,
            'K_q875_mean': s['K_q875_mean'],
            'K_q875_sd': s['K_q875_sd'],
            'K_q875_se': s['K_q875_se'],
            'xi_X_mean': s['xi_mean'],
            'xi_X_sd': s['xi_sd'],
            'xi_X_se': s['xi_se'],
            'E_v_q875_mean': s['E_v_mean'],
            'E_v_q875_sd': s['E_v_sd'],
            'K_full_mean': s['K_full_mean'],
            'K_q125_mean': s['K_q125_mean'],
            'K_q975_mean': s['K_q975_mean'],
            'n_failed': s['n_failed_total'],
        })
    out_csv2 = out_dir / "59_q875_extension_summary.csv"
    pl.DataFrame(summary_rows).write_csv(out_csv2)
    print(f"[save] {out_csv2}")


if __name__ == "__main__":
    main()
