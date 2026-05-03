"""
Experiment 56 — push K_q125 + E[v]_q125 to N=2^42, verify K(v) identity.

For each N ∈ {2^25, 2^27, 2^28, 2^30, 2^32, 2^34, 2^36, 2^38, 2^40, 2^42}:
  - Walk 500K random odd starts × 5 seeds
  - For each orbit, record: σ, n_odd (Syracuse steps), n_even (halvings),
    s_to_4_thresholds
  - Per-orbit v = n_even / n_odd (mean v across orbit's Syracuse steps)
  - Bottom σ-quartile band: q=0.125
  - K_q125 = first-passage slope of (σ̄ - s̄(f)) on log(f), 4-threshold
  - E[v]_q125 = mean v over bottom-quartile orbits
  - K_predicted(v) = (1 + E[v]) / (E[v]·log(2) - log(3))   ← K(mean v)
  - K_avg_per_orbit = mean over band of K(v_per_orbit)     ← <K(v)> direct

Identity check: does K_q125 (regression-based) match K(E[v]_q125) (algebraic)?
For homogeneous v across band, they should match. Jensen gap measures
within-band v-variance.

Also report E[v]_q125 trajectory and fit to test asymptote question.

Memory note: at 2^42, max orbit excursion may approach int64 limit. Walker
flags failed orbits via `m > max_value // 3` check before 3m+1.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

K_H_USER = 3.0 / (np.log(4.0) - np.log(3.0))
LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
MAX_VAL = np.int64(2**62)


def K_of_v(v):
    """K(v) = (1+v) / (v·log2 - log3)"""
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

    # Bottom σ-quartile = q=0.125 band
    q25 = np.percentile(sigma_ok, 25)
    mask_low = sigma_ok <= q25

    K_q125 = K_eff_slope(sigma_ok[mask_low], s_ok[mask_low], raw_ok[mask_low])
    E_v_q125 = float(v_per_orbit[mask_low].mean())
    K_predicted = K_of_v(E_v_q125)
    K_avg_per_orbit = float(K_of_v(v_per_orbit[mask_low]).mean())

    K_full = K_eff_slope(sigma_ok, s_ok, raw_ok)
    E_v_all = float(v_per_orbit.mean())

    return {
        'K_q125': K_q125,
        'E_v_q125': E_v_q125,
        'K_predicted_from_E_v': K_predicted,
        'K_avg_per_orbit_q125': K_avg_per_orbit,
        'K_full': K_full,
        'E_v_all': E_v_all,
        'n_total': int(ok.sum()),
        'n_failed': int((~ok).sum()),
        'n_low': int(mask_low.sum()),
    }


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N_list = [25, 27, 28, 30, 32, 34, 36, 38, 40, 42]
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 500_000

    print(f"# K_q125 + E[v]_q125 trajectory + identity check")
    print(f"# {len(log2N_list)} N × {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# K_h = {K_H_USER:.4f}; K(v→∞) = {1/LOG_2:.4f}")
    print()

    summary = {}
    for log2N in log2N_list:
        N = 1 << log2N
        seed_results = []
        for seed in seeds:
            t0 = time.perf_counter()
            r = analyze_one_seed(N, seed, n_starts)
            elapsed = time.perf_counter() - t0
            seed_results.append(r)
            print(f"  2^{log2N:>2} seed={seed:>5}: K_q125={r['K_q125']:.4f}  E[v]_q125={r['E_v_q125']:.4f}  "
                  f"K(E[v])={r['K_predicted_from_E_v']:.4f}  ⟨K(v)⟩={r['K_avg_per_orbit_q125']:.4f}  "
                  f"K_full={r['K_full']:.4f}  failed={r['n_failed']}  ({elapsed:.1f}s)")
        K_q125_arr = np.array([r['K_q125'] for r in seed_results])
        E_v_q125_arr = np.array([r['E_v_q125'] for r in seed_results])
        K_pred_arr = np.array([r['K_predicted_from_E_v'] for r in seed_results])
        K_avg_arr = np.array([r['K_avg_per_orbit_q125'] for r in seed_results])
        K_full_arr = np.array([r['K_full'] for r in seed_results])
        n_failed = sum(r['n_failed'] for r in seed_results)
        summary[log2N] = {
            'K_q125_mean': K_q125_arr.mean(), 'K_q125_sd': K_q125_arr.std(ddof=1),
            'E_v_q125_mean': E_v_q125_arr.mean(), 'E_v_q125_sd': E_v_q125_arr.std(ddof=1),
            'K_pred_mean': K_pred_arr.mean(),
            'K_avg_per_orbit_mean': K_avg_arr.mean(),
            'K_full_mean': K_full_arr.mean(), 'K_full_sd': K_full_arr.std(ddof=1),
            'n_failed_total': n_failed,
        }

    # Trajectory tables
    print(f"\n# K_q125 trajectory:")
    print(f"#   {'log2N':>6}  {'K_q125':>9}  {'K_q125_SD':>9}  {'E[v]_q125':>10}  {'K(E[v])':>9}  {'⟨K(v)⟩':>9}  "
          f"{'identity gap':>13}  {'K_full':>9}  {'failed':>7}")
    for ln in log2N_list:
        s = summary[ln]
        gap_identity = s['K_q125_mean'] - s['K_avg_per_orbit_mean']
        print(f"#   {ln:>6}  {s['K_q125_mean']:>9.4f}  {s['K_q125_sd']:>9.4f}  "
              f"{s['E_v_q125_mean']:>10.4f}  {s['K_pred_mean']:>9.4f}  {s['K_avg_per_orbit_mean']:>9.4f}  "
              f"{gap_identity:>+13.4f}  {s['K_full_mean']:>9.4f}  {s['n_failed_total']:>7}")

    # Identity verification
    print(f"\n# Identity check: K_q125 (regression) vs ⟨K(v)⟩ (per-orbit avg)")
    print(f"#   If equal: K_q125 = mean of K(v_orbit) over band orbits")
    print(f"#   If gap: residual structural deviation (regression slope ≠ orbit-avg slope)")

    # E[v] vs log(N) fits
    log2N_arr = np.array(log2N_list)
    log_N = log2N_arr * LOG_2
    E_v_arr = np.array([summary[ln]['E_v_q125_mean'] for ln in log2N_list])
    K_q125_arr = np.array([summary[ln]['K_q125_mean'] for ln in log2N_list])

    print(f"\n# E[v]_q125 vs log(N) fits:")
    # Linear: E[v] = a + b·log N
    xc = log_N - log_N.mean(); yc = E_v_arr - E_v_arr.mean()
    slope_lin = float((xc*yc).sum() / (xc*xc).sum())
    intercept_lin = float(E_v_arr.mean() - slope_lin*log_N.mean())
    pred_lin = intercept_lin + slope_lin*log_N
    R2_lin = 1 - float(np.sum((E_v_arr - pred_lin)**2) / np.sum((E_v_arr - E_v_arr.mean())**2))
    print(f"#   linear:  E[v] = {intercept_lin:.4f} + {slope_lin:.5f}·log(N)   R²={R2_lin:.4f}")
    if slope_lin < 0:
        v_at_inf = "asymptote unbounded below"
        # Linear extrapolation to E[v]=2:
        if intercept_lin > 2:
            log_N_at_2 = (2.0 - intercept_lin) / slope_lin
            log2N_at_2 = log_N_at_2 / LOG_2
            print(f"#   linear extrap E[v]→2: at log2N ≈ {log2N_at_2:.1f}")

    # Plateau form: E[v] = E_inf + a·N^(-β)
    from scipy.optimize import minimize
    def plateau_loss(params):
        E_inf, a, beta = params
        if beta <= 0 or beta > 2: return 1e9
        pred = E_inf + a * np.exp(-beta * log_N)
        return float(np.sum((E_v_arr - pred)**2))
    best = None
    for E_init in [2.0, 2.1, 2.2, E_v_arr.min()]:
        for a_init in [0.5, 1.0, 5.0, 50]:
            for b_init in [0.05, 0.1, 0.3, 0.5]:
                res = minimize(plateau_loss, x0=[E_init, a_init, b_init],
                               method='Nelder-Mead',
                               options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 5000})
                if best is None or res.fun < best.fun:
                    best = res
    E_inf, a_pl, b_pl = best.x
    pred_pl = E_inf + a_pl * np.exp(-b_pl * log_N)
    R2_pl = 1 - float(np.sum((E_v_arr - pred_pl)**2) / np.sum((E_v_arr - E_v_arr.mean())**2))
    print(f"#   plateau: E[v] = {E_inf:.4f} + {a_pl:.4f}·N^(-{b_pl:.4f})   R²={R2_pl:.4f}")
    print(f"#     E[v]_∞ = {E_inf:.4f}  (typical Geom(1/2) E[v] = 2.0)")
    if E_inf > 0.001:
        print(f"#     Predicted K_q125 asymptote = K(E[v]_∞) = {K_of_v(E_inf):.4f}")
    print(f"#     gap E[v]_∞ - 2.0 = {E_inf - 2.0:+.4f}")

    # Per-octave Δ E[v]_q125
    print(f"\n# Per-octave Δ E[v]_q125 (does it slow?):")
    for i in range(1, len(log2N_list)):
        d = (E_v_arr[i] - E_v_arr[i-1]) / (log2N_list[i] - log2N_list[i-1])
        print(f"#   {log2N_list[i-1]}→{log2N_list[i]}: ΔE[v]/oct = {d:+.5f}")

    # Save
    rows = []
    for ln in log2N_list:
        s = summary[ln]
        rows.append({
            'log2N': ln, 'N': 1 << ln,
            'K_q125_mean': float(s['K_q125_mean']),
            'K_q125_sd': float(s['K_q125_sd']),
            'E_v_q125_mean': float(s['E_v_q125_mean']),
            'E_v_q125_sd': float(s['E_v_q125_sd']),
            'K_pred_from_E_v': float(s['K_pred_mean']),
            'K_avg_per_orbit': float(s['K_avg_per_orbit_mean']),
            'identity_gap': float(s['K_q125_mean'] - s['K_avg_per_orbit_mean']),
            'K_full_mean': float(s['K_full_mean']),
            'K_full_sd': float(s['K_full_sd']),
            'n_failed': int(s['n_failed_total']),
        })
    out_csv = out_dir / "56_q125_to_2_42_with_v.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
