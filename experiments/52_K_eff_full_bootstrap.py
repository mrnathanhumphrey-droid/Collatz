"""
Experiment 52 — re-bootstrap K_full at N ∈ {2^25, 2^27, 2^28, 2^30, 2^32, 2^34}
with fresh random-orbit sampling, 5 seeds each, to make all per-octave Δ
estimates comparable.

Prior K_full estimates (2^25..2^30) were single-seed subsamples from
precomputed parquet — the seed variance was implicit. The 32→34 vs 30→32
deceleration question requires apples-to-apples bootstrap at every N.
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


def K_full_one_seed(N, seed, n_starts):
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
    log_f = np.log(np.maximum(raw_ok, 1.0))
    R = sigma_ok[:, None] - s_phys_ok.astype(np.float64)
    x = log_f.mean(axis=0); y = R.mean(axis=0)
    xc = x - x.mean(); yc = y - y.mean()
    return float((xc * yc).sum() / (xc * xc).sum())


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N_list = [25, 27, 28, 30, 32, 34]
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 500_000

    print(f"# Bootstrap K_full at N ∈ {{2^{log2N_list[0]}, ..., 2^{log2N_list[-1]}}}")
    print(f"# {n_starts:,} orbits per seed × {len(seeds)} seeds")
    print()

    summary = {}
    for log2N in log2N_list:
        N = 1 << log2N
        K_seeds = []
        for seed in seeds:
            t0 = time.perf_counter()
            K = K_full_one_seed(N, seed, n_starts)
            elapsed = time.perf_counter() - t0
            K_seeds.append(K)
        K_arr = np.array(K_seeds)
        mean = K_arr.mean()
        sd = K_arr.std(ddof=1)
        se = sd / np.sqrt(len(seeds))
        summary[log2N] = {'mean': mean, 'sd': sd, 'se': se, 'seeds': K_seeds}
        print(f"  2^{log2N}: mean = {mean:.5f}  sd = {sd:.5f}  se = {se:.5f}  "
              f"seeds = {[f'{k:.4f}' for k in K_seeds]}")

    print(f"\n# Per-octave Δ K_full (with bootstrap CIs):")
    print(f"#   {'step':>10}  {'Δ mean':>9}  {'Δ/octave':>10}  {'SE on Δ':>9}")
    for i in range(1, len(log2N_list)):
        a = log2N_list[i-1]; b = log2N_list[i]
        d = b - a
        delta = summary[b]['mean'] - summary[a]['mean']
        se_a = summary[a]['se']; se_b = summary[b]['se']
        se_delta = np.sqrt(se_a**2 + se_b**2)
        print(f"  {a}→{b}:  {delta:+9.5f}  {delta/d:+10.5f}  {se_delta:>9.5f}")

    # Re-fit candidate functional forms
    log2N_arr = np.array(log2N_list)
    K_arr = np.array([summary[k]['mean'] for k in log2N_list])
    gap_arr = K_H_USER - K_arr
    log_N = log2N_arr * np.log(2)

    print(f"\n# Bootstrap-mean K_full and gap to K_h={K_H_USER:.4f}:")
    print(f"#   {'log2N':>6}  {'K_full':>9}  {'gap':>8}")
    for i, log2N in enumerate(log2N_list):
        print(f"#   {log2N:>6}  {K_arr[i]:>9.5f}  {gap_arr[i]:>8.5f}")

    # Power-law fit
    log_gap = np.log(gap_arr)
    xc = log_N - log_N.mean(); yc = log_gap - log_gap.mean()
    slope_pl = float((xc*yc).sum() / (xc*xc).sum())
    alpha_pl = -slope_pl
    intercept_pl = log_gap.mean() - slope_pl * log_N.mean()
    A_pl = float(np.exp(intercept_pl))
    pred_pl = intercept_pl + slope_pl * log_N
    R2_pl = 1 - float(np.sum((log_gap - pred_pl)**2) / np.sum((log_gap - log_gap.mean())**2))
    print(f"\n# Power-law fit gap = A · N^(-α) (6-pt bootstrap means):")
    print(f"#   α = {alpha_pl:.4f}  A = {A_pl:.4f}  R² = {R2_pl:.4f}")

    # 1/log(N) fit
    log_logN = np.log(log_N)
    xc = log_logN - log_logN.mean(); yc = log_gap - log_gap.mean()
    slope_lg = float((xc*yc).sum() / (xc*xc).sum())
    alpha_lg = -slope_lg
    intercept_lg = log_gap.mean() - slope_lg * log_logN.mean()
    A_lg = float(np.exp(intercept_lg))
    pred_lg = intercept_lg + slope_lg * log_logN
    R2_lg = 1 - float(np.sum((log_gap - pred_lg)**2) / np.sum((log_gap - log_gap.mean())**2))
    print(f"\n# (log N)^(-α) fit:")
    print(f"#   α = {alpha_lg:.4f}  A = {A_lg:.4f}  R² = {R2_lg:.4f}")

    # Plateau fit: gap = (K_h - K_inf) + (decaying)·exp(-...) — try simpler:
    # gap = constant + power-law: K_full = K_inf + a·N^(-beta); fit (K_inf, a, beta)
    from scipy.optimize import minimize
    def plateau_loss(params):
        K_inf, a, beta = params
        if beta <= 0 or a <= 0 or K_inf > K_H_USER + 0.1: return 1e9
        pred = K_inf - a * np.exp(-beta * log_N)  # K_full predicted
        return float(np.sum((K_arr - pred)**2))
    res = minimize(plateau_loss, x0=[10.30, 5.0, 0.1], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 10000})
    K_inf, a_pl_fit, beta_pl_fit = res.x
    pred_pl_full = K_inf - a_pl_fit * np.exp(-beta_pl_fit * log_N)
    R2_plateau = 1 - float(np.sum((K_arr - pred_pl_full)**2) /
                            np.sum((K_arr - K_arr.mean())**2))
    print(f"\n# Plateau fit K_full = K_inf - a·N^(-β):")
    print(f"#   K_inf = {K_inf:.5f}  (K_h = {K_H_USER:.5f}; gap K_h - K_inf = {K_H_USER - K_inf:.5f})")
    print(f"#   a = {a_pl_fit:.4f}  β = {beta_pl_fit:.5f}  R² = {R2_plateau:.5f}")
    print(f"\n# Predictions per fit:")
    print(f"#   {'log2N':>6}  {'empir':>8}  {'pred PL':>9}  {'pred 1/log':>10}  {'pred plateau':>13}")
    for i, log2N in enumerate(log2N_list):
        K_emp = K_arr[i]
        K_pl_pred = K_H_USER - A_pl * np.exp(-alpha_pl * log_N[i])
        K_lg_pred = K_H_USER - A_lg * np.exp(-alpha_lg * np.log(log_N[i]))
        K_pl_p_pred = pred_pl_full[i]
        print(f"#   {log2N:>6}  {K_emp:>8.4f}  {K_pl_pred:>9.4f}  {K_lg_pred:>10.4f}  "
              f"{K_pl_p_pred:>13.4f}")

    # Save
    rows = []
    for log2N in log2N_list:
        s = summary[log2N]
        for i, K in enumerate(s['seeds']):
            rows.append({
                'log2N': log2N, 'N': 1 << log2N,
                'seed': seeds[i], 'K_full': float(K),
                'K_full_mean': s['mean'], 'K_full_sd': s['sd'],
            })
    out_csv = out_dir / "52_K_eff_full_bootstrap.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
