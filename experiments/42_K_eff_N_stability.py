"""
Experiment 42 — K_eff N-stability discriminator.

Tests whether K_eff = K_h − slope_on_log(f) is stable across N (rules out
finite-N as the mechanism) or drifts with N (mechanism identified).

Background: Result 3 in closed_form_findings.md identified K_eff ≈ 9.94
(mixed) / 10.11 (first-passage only) at N = 2²⁷ via TA.3 fit. Mechanism is
hypothesized as finite-N μ_β oscillation. If K_eff drifts toward K_h with
N, finite-N is confirmed. If K_eff is N-stable, finite-N is ruled out.

Method (per user, "quick test on existing data"):
  For each parquet at N ∈ {2²⁰, 2²², 2²³, 2²⁴, 2²⁵, 2²⁷}:
    1. Sample 100K odd starts.
    2. Walk each orbit to first-passage of 4 thresholds:
       f ∈ {N^(2/3), √N · log N, √N, √N / log N}.
    3. s_mean(f) = mean steps from start to first-crossing of f.
    4. R(f) = σ_mean − s_mean(f) (segment "below f"), or equivalently
       slope on log(f) of (s_mean(f) + K_h · log(f) − K_h · log(N))
       gives K_h − K_eff.
    5. Linear fit slope on log(f) of either R(f) or ε(f); K_eff = K_h − slope.
  Print K_eff vs log(N).

If K_eff drifts upward with N (toward K_h = 10.4282): finite-N μ_β confirmed.
If K_eff is flat: finite-N ruled out, look elsewhere.

Usage:
    python 42_K_eff_N_stability.py
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))


@njit(parallel=True, cache=True)
def first_passage_4thresh(starts, thresholds_per_n):
    """For each start n, walk Collatz to first-crossing of each of 4
    thresholds. thresholds_per_n[i] is sorted descending so we can record
    in order. Returns s_arr[i, j] = step count at first crossing of
    j-th (descending-sorted) threshold."""
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


def fit_slope(x, y):
    x_, y_ = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    xm, ym = x_.mean(), y_.mean()
    den = ((x_ - xm) ** 2).sum()
    if den < 1e-20:
        return float("nan"), float("nan"), float("nan")
    slope = ((x_ - xm) * (y_ - ym)).sum() / den
    intercept = ym - slope * xm
    pred = slope * x_ + intercept
    ss_res = ((y_ - pred) ** 2).sum()
    ss_tot = ((y_ - ym) ** 2).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-20 else 1.0
    return float(slope), float(intercept), float(r2)


def compute_K_eff_at_N(data_dir, N, n_samples=100_000, seed=42):
    """Compute K_eff at one N value using 4-threshold first-passage."""
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int64)
    n_total = len(n_arr)

    # Use all rows if total is small; else random sample
    rng = np.random.default_rng(seed)
    if n_total > n_samples:
        idx = rng.choice(n_total, size=n_samples, replace=False)
        n_arr = n_arr[idx]
        sigma_arr = sigma_arr[idx]

    n_used = len(n_arr)
    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)

    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)
    thresholds_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    threshold_names = ["N^(2/3)", "sqrt(N)*log(N)", "sqrt(N)", "sqrt(N)/log(N)"]

    s_sorted = first_passage_4thresh(n_arr, thresholds_sorted)

    # De-permute back to physical thresholds
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros_like(s_sorted)
    thresh_phys = np.zeros_like(thresholds_sorted)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(
            s_sorted, rev[:, None], axis=1).flatten()
        thresh_phys[:, col] = np.take_along_axis(
            thresholds_sorted, rev[:, None], axis=1).flatten()

    # Drop unreached rows for safety
    valid_all = np.all(s_phys >= 0, axis=1)
    s_phys = s_phys[valid_all]
    thresh_phys = thresh_phys[valid_all]
    sigma_arr_v = sigma_arr[valid_all]
    log_n_v = log_n[valid_all]

    # For each threshold, compute s_mean and log(f)_mean
    s_mean_per_thresh = np.zeros(4)
    log_f_per_thresh = np.zeros(4)
    R_per_thresh = np.zeros(4)
    for col in range(4):
        s_mean_per_thresh[col] = s_phys[:, col].mean()
        log_f_per_thresh[col] = np.log(thresh_phys[:, col].astype(np.float64)).mean()
        # R(f) = sigma_mean - s_mean(f)  (segment from f to 1)
        R_per_thresh[col] = sigma_arr_v.mean() - s_mean_per_thresh[col]

    # K_eff via slope of R(f) on log(f). Per Result 3 line 104:
    #   R(f) ≈ K_eff · log(f) + const
    K_eff_R, R_int, R_r2 = fit_slope(log_f_per_thresh, R_per_thresh)

    # Cross-check via ε(f) = s_mean(f) - K_h · log(N/f) - <α_det>
    #   ε(f) ≈ ε(σ) + (K_h - K_eff) · log(f)
    # using <α_det> = log(6)/log(4/3) = 6.228263 (Result 1)
    alpha_det_mean = np.log(6.0) / np.log(4.0 / 3.0)
    log_N = np.log(float(N))
    eps_per_thresh = (s_mean_per_thresh - K_H * (log_N - log_f_per_thresh)
                      - alpha_det_mean)
    slope_eps, eps_int, eps_r2 = fit_slope(log_f_per_thresh, eps_per_thresh)
    K_eff_eps = K_H - slope_eps

    return {
        "N": N,
        "log_N": log_N,
        "n_used": n_used,
        "sigma_mean": float(sigma_arr_v.mean()),
        "s_mean": s_mean_per_thresh,
        "log_f": log_f_per_thresh,
        "R": R_per_thresh,
        "eps": eps_per_thresh,
        "K_eff_R": K_eff_R,
        "K_eff_R_R2": R_r2,
        "K_eff_eps": K_eff_eps,
        "K_eff_eps_R2": eps_r2,
        "slope_eps": slope_eps,
    }


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    N_values = [
        1 << 20,   # 2^20 = 1,048,576
        1 << 22,   # 2^22 = 4,194,304
        1 << 23,   # 2^23 = 8,388,608
        1 << 24,   # 2^24 = 16,777,216
        1 << 25,   # 2^25 = 33,554,432
        1 << 27,   # 2^27 = 134,217,728
    ]

    print(f"K_eff N-stability sweep")
    print(f"K_h = {K_H:.6f}")
    print(f"<α_det> = log(6)/log(4/3) = {np.log(6.0)/np.log(4.0/3.0):.6f}")
    print(f"Method: K_eff_R via slope of R(f) = σ̄ − s̄(f) on log(f)")
    print(f"        K_eff_eps via K_h − slope of ε(f) on log(f)")
    print()

    results = []
    for N in N_values:
        log_N = np.log(float(N))
        print(f"[N = 2^{int(round(log_N/np.log(2))):>2}] = {N:>14,d}  ", end="", flush=True)
        t0 = time.perf_counter()
        try:
            r = compute_K_eff_at_N(data_dir, N)
            elapsed = time.perf_counter() - t0
            print(f"({elapsed:.1f}s)")
            print(f"  used {r['n_used']:>8,d} odd starts, σ̄ = {r['sigma_mean']:.4f}")
            print(f"  s̄(f) per threshold:  {[f'{x:.2f}' for x in r['s_mean']]}")
            print(f"  log(f) per threshold: {[f'{x:.3f}' for x in r['log_f']]}")
            print(f"  R(f) per threshold:   {[f'{x:.3f}' for x in r['R']]}")
            print(f"  ε(f) per threshold:   {[f'{x:.3f}' for x in r['eps']]}")
            print(f"  -> K_eff (via R fit):  {r['K_eff_R']:.4f}    R^2 = {r['K_eff_R_R2']:.4f}")
            print(f"  -> K_eff (via ε fit):  {r['K_eff_eps']:.4f}    R^2 = {r['K_eff_eps_R2']:.4f}")
            print(f"  -> slope ε on log(f) = {r['slope_eps']:.4f} (≡ K_h − K_eff)")
            results.append(r)
        except Exception as e:
            print(f"FAILED: {e}")
        print()
        gc.collect()

    # Summary table
    print(f"=== K_eff vs log(N) ===")
    print(f"  {'N':>14}  {'log N':>7}  {'K_eff(R)':>9}  {'K_eff(ε)':>9}  "
          f"{'slope ε':>8}  {'K_h−K_eff':>10}")
    for r in results:
        N = r["N"]
        kpow = int(round(r["log_N"] / np.log(2)))
        K_eff_avg = (r["K_eff_R"] + r["K_eff_eps"]) / 2
        gap = K_H - K_eff_avg
        print(f"  2^{kpow:<3}={N:>10,d}  {r['log_N']:>7.4f}  "
              f"{r['K_eff_R']:>9.4f}  {r['K_eff_eps']:>9.4f}  "
              f"{r['slope_eps']:>+8.4f}  {gap:>10.4f}")

    # Slope of K_eff vs log(N) (drift test)
    if len(results) >= 3:
        log_Ns = np.array([r["log_N"] for r in results])
        K_effs_R = np.array([r["K_eff_R"] for r in results])
        K_effs_eps = np.array([r["K_eff_eps"] for r in results])
        K_effs_avg = (K_effs_R + K_effs_eps) / 2

        slope_R, int_R, r2_R = fit_slope(log_Ns, K_effs_R)
        slope_eps, int_eps, r2_eps = fit_slope(log_Ns, K_effs_eps)
        slope_avg, int_avg, r2_avg = fit_slope(log_Ns, K_effs_avg)

        print(f"\n=== Drift of K_eff with log(N) ===")
        print(f"  K_eff(R)   = {int_R:.4f} + {slope_R:>+.5f} · log(N)    R^2 = {r2_R:.3f}")
        print(f"  K_eff(ε)   = {int_eps:.4f} + {slope_eps:>+.5f} · log(N)    R^2 = {r2_eps:.3f}")
        print(f"  K_eff(avg) = {int_avg:.4f} + {slope_avg:>+.5f} · log(N)    R^2 = {r2_avg:.3f}")
        print(f"  K_h        = {K_H:.4f}")
        print(f"  Range of K_eff: [{K_effs_avg.min():.3f}, {K_effs_avg.max():.3f}]   "
              f"spread = {K_effs_avg.max() - K_effs_avg.min():.3f}")
        print()
        # Verdict
        spread = K_effs_avg.max() - K_effs_avg.min()
        slope_mag = abs(slope_avg)
        if slope_mag > 0.05 and r2_avg > 0.5:
            verdict = "DRIFTS with N (finite-N mechanism CONFIRMED)"
        elif spread < 0.10:
            verdict = "STABLE across N (finite-N mechanism RULED OUT)"
        else:
            verdict = "AMBIGUOUS (drift small, R² low)"
        print(f"  Verdict: {verdict}")
        print(f"  Direction: {'toward K_h' if (slope_avg > 0 and K_effs_avg.mean() < K_H) or (slope_avg < 0 and K_effs_avg.mean() > K_H) else 'away from K_h'}")

    print(f"\n[done]")


if __name__ == "__main__":
    main()
