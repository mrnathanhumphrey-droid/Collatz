"""
Experiment 55 — push K_eff_band(q=0.125) to higher N to track the band's
monotone climb until it plateaus or crosses K_h.

Per-band tracking only (we still walk all orbits to identify the bottom
quartile; the savings vs full K_full bootstrap is in output dimensionality
and analytical clarity).

Bootstrap (5 seeds) at N ∈ {2^38, 2^40, 2^42, 2^44}.

Empirical structural lower bound argument:
  K(v) = (1+v) / (v·log2 - log3),  E[v]→∞ gives K→1/log2 = 1.443
  Empirical q=0.125 = 6.97 ⇔ E[v] ≈ 2.26 (vs typical 2.0)
  Plateau target unknown; expected in (4, 9) range based on realistic
  v-distributions for bottom-σ-quartile orbits.
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


def K_from_v(v_arr):
    """K(v) from random-walk heuristic. Bound: v→∞ gives 1/log(2)."""
    return (1 + v_arr) / (v_arr * LOG_2 - LOG_3)


@njit(parallel=True, cache=True)
def walk_K_with_v(starts, thresholds_per_n, max_value, max_steps):
    """Walk orbits, also accumulate v-step distribution for each orbit.

    Returns sigma, s_per_thresh, ok, sum_v_per_orbit (= total halvings = even
    steps), n_syracuse_per_orbit (= odd steps).
    """
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


def analyze_one_seed(N, seed, n_starts, max_steps=500_000):
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
    sigma, s, ok, n_odd, n_even = walk_K_with_v(starts, th, MAX_VAL, max_steps)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros((len(starts), 4), dtype=np.int32)
    for col in range(4):
        s_phys[:, col] = np.take_along_axis(s, inv_sort[:, col:col+1], axis=1).flatten()
    sigma_ok = sigma[ok].astype(np.float64)
    s_ok = s_phys[ok]
    raw_ok = raw[ok]
    n_odd_ok = n_odd[ok].astype(np.float64)
    n_even_ok = n_even[ok].astype(np.float64)

    # Bottom σ-quartile = q=0.125 band
    q25 = np.percentile(sigma_ok, 25)
    mask_low = sigma_ok <= q25
    K_q125 = K_eff_slope(sigma_ok[mask_low], s_ok[mask_low], raw_ok[mask_low])

    # Effective E[v] in low-σ band: E[v] = (n_odd + n_even) / n_odd... wait
    # Each orbit: sum of v_i over Syracuse steps = n_even (total halvings).
    # n_Syracuse = n_odd. So E[v | orbit] = n_even / n_odd.
    # E[v] over band = mean of n_even / n_odd over band orbits.
    v_per_orbit = n_even_ok / np.maximum(n_odd_ok, 1)
    E_v_low = float(v_per_orbit[mask_low].mean())
    E_v_all = float(v_per_orbit.mean())

    # K_full for reference
    K_full = K_eff_slope(sigma_ok, s_ok, raw_ok)

    return {
        'K_q125': K_q125, 'K_full': K_full,
        'E_v_low': E_v_low, 'E_v_all': E_v_all,
        'n_total': int(ok.sum()),
        'n_low': int(mask_low.sum()),
    }


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N_list = [38, 40]
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 500_000

    print(f"# K_eff_band(q=0.125) push to higher N")
    print(f"# {len(log2N_list)} N × {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# K_h = {K_H_USER:.4f}; K_min(v→∞) = {1/LOG_2:.4f}")
    print()

    results = {}
    for log2N in log2N_list:
        N = 1 << log2N
        seed_results = []
        for seed in seeds:
            t0 = time.perf_counter()
            r = analyze_one_seed(N, seed, n_starts)
            elapsed = time.perf_counter() - t0
            seed_results.append(r)
            print(f"  2^{log2N} seed={seed:>5}: K_q125={r['K_q125']:.4f}  K_full={r['K_full']:.4f}  "
                  f"E[v]_low={r['E_v_low']:.4f}  E[v]_all={r['E_v_all']:.4f}  ({elapsed:.1f}s)")
        K_q125_arr = np.array([r['K_q125'] for r in seed_results])
        K_full_arr = np.array([r['K_full'] for r in seed_results])
        E_v_low_arr = np.array([r['E_v_low'] for r in seed_results])
        E_v_all_arr = np.array([r['E_v_all'] for r in seed_results])
        results[log2N] = {
            'K_q125_mean': K_q125_arr.mean(), 'K_q125_sd': K_q125_arr.std(ddof=1),
            'K_full_mean': K_full_arr.mean(), 'K_full_sd': K_full_arr.std(ddof=1),
            'E_v_low_mean': E_v_low_arr.mean(), 'E_v_low_sd': E_v_low_arr.std(ddof=1),
            'E_v_all_mean': E_v_all_arr.mean(),
            'seeds': seed_results,
        }

    # Combined trajectory
    prev = {25: 6.266, 27: 6.404, 28: 6.482, 30: 6.589, 32: 6.725, 34: 6.867, 36: 6.971}
    print(f"\n# K_eff_band(q=0.125) trajectory (combining prior + new):")
    print(f"#   {'log2N':>6}  {'K_q125':>9}  {'Δ from prev':>13}  {'Δ/oct':>8}  {'gap to K_h':>11}")
    all_log2N = sorted(list(prev.keys()) + list(results.keys()))
    K_traj = []
    for ln in all_log2N:
        K = prev[ln] if ln in prev else results[ln]['K_q125_mean']
        K_traj.append(K)
    for i, ln in enumerate(all_log2N):
        K = K_traj[i]
        if i == 0:
            print(f"#   {ln:>6}  {K:>9.4f}  {'—':>13}  {'—':>8}  {K_H_USER - K:>11.4f}")
        else:
            d = K - K_traj[i-1]
            d_oct = d / (ln - all_log2N[i-1])
            print(f"#   {ln:>6}  {K:>9.4f}  {d:>+13.5f}  {d_oct:>+8.5f}  {K_H_USER - K:>11.4f}")

    # E[v] trajectory
    print(f"\n# E[v] in bottom-σ-quartile band:")
    print(f"#   {'log2N':>6}  {'E[v]_low':>9}  {'E[v]_all':>9}  {'K_from_E[v]_low':>16}")
    for ln in log2N_list:
        r = results[ln]
        K_from_v = (1 + r['E_v_low_mean']) / (r['E_v_low_mean']*LOG_2 - LOG_3)
        print(f"#   {ln:>6}  {r['E_v_low_mean']:>9.4f}  {r['E_v_all_mean']:>9.4f}  {K_from_v:>16.4f}")

    # Solve for empirical-implied E[v] from K_q125
    print(f"\n# Empirical-K → implied E[v] (inverting K(v) formula):")
    print(f"#   {'log2N':>6}  {'K_q125':>9}  {'implied E[v]':>13}")
    for ln in all_log2N:
        K = K_traj[all_log2N.index(ln)]
        # K = (1+v) / (v·log2 - log3), solve for v:
        # K·v·log2 - K·log3 = 1 + v
        # v(K·log2 - 1) = 1 + K·log3
        # v = (1 + K·log3) / (K·log2 - 1)
        v_implied = (1 + K * LOG_3) / (K * LOG_2 - 1)
        print(f"#   {ln:>6}  {K:>9.4f}  {v_implied:>13.4f}")

    # Per-octave Δ trend (does it slow?)
    print(f"\n# Per-octave Δ K_q125 trend (does the climb decelerate?):")
    for i in range(1, len(all_log2N)):
        ln_prev = all_log2N[i-1]
        ln = all_log2N[i]
        d_oct = (K_traj[i] - K_traj[i-1]) / (ln - ln_prev)
        marker = " ← new" if ln in results else ""
        print(f"#   {ln_prev}→{ln}: {d_oct:+.5f}/oct{marker}")

    # Save
    rows = []
    for ln, r in results.items():
        for j, seed_r in enumerate(r['seeds']):
            rows.append({
                'log2N': ln, 'N': 1 << ln, 'seed': seeds[j],
                'K_q125': float(seed_r['K_q125']),
                'K_full': float(seed_r['K_full']),
                'E_v_low': float(seed_r['E_v_low']),
                'E_v_all': float(seed_r['E_v_all']),
                'K_q125_mean': float(r['K_q125_mean']),
                'K_q125_sd': float(r['K_q125_sd']),
            })
    out_csv = out_dir / "55_q125_high_N.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
