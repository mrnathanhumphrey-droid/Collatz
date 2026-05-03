"""
Experiment 40 — decompose K_eff into odd-σ slope and halving-chain compensation.

User's hypothesis: K_eff (the slope on log(threshold) − K_h ≈ 9.94 in older
fit, 10.11 in cleaner first-passage-only fit) splits into:

  K_eff = slope_odd_only + halving_chain_compensation

where slope_odd_only is the trajectory-weighted regression of σ(v) on log(v)
restricted to odd v (≈ 9.31 empirical), and halving_chain_compensation is
the contribution from even visits in the orbit (≈ 0.63 empirical).

Test:
  1. Walk many Collatz trajectories from random odd starts. Record each
     visited value and its parity. For v ≤ threshold, count visits.
  2. Look up σ(v) from the main parquet.
  3. Compute weighted regression slope of σ(v) on log(v) under the
     trajectory-visit measure π(v):
        - all visits: slope_all
        - odd visits: slope_odd
        - even visits: slope_even
        - "halving compensation": slope_all − slope_odd
  4. Compare to the K_eff value from the first-passage analysis.

Then test the closed-form prediction for the halving-chain compensation:
under Geom(1/2) on v_2 along trajectories, mean halving-chain length per
Syracuse step is 2. Each halving step adds log(2) to log(v) and 1 to the
remaining-σ count, with no descent.

Threshold: v ≤ 5000 (per user's framing). Also report at v ≤ √N for the
direct comparison.

Usage:
    python 40_K_eff_decomposition.py
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
def collect_visit_counts(starts, threshold):
    """For each start, walk Collatz trajectory. Count visits to each value
    v in [1, threshold]. Returns visit_count[v] for v=0..threshold (with
    visit_count[0] unused). Chunk-parallel to avoid races."""
    n = len(starts)
    chunk = max(1, n // 256)
    n_chunks = (n + chunk - 1) // chunk
    counts_local = np.zeros((n_chunks, threshold + 1), dtype=np.int64)
    max_steps = 20000
    for c in prange(n_chunks):
        start_i = c * chunk
        end_i = min(start_i + chunk, n)
        for i in range(start_i, end_i):
            m = starts[i]
            steps = 0
            while m > 1 and steps < max_steps:
                if m <= threshold:
                    counts_local[c, m] += 1
                if m & 1:
                    m = 3 * m + 1
                else:
                    m = m >> 1
                steps += 1
            # Final visit at m=1
            if m == 1:
                counts_local[c, 1] += 1
    return counts_local.sum(axis=0)


def fit_weighted_slope(v_vals, sigma_vals, weights):
    """Weighted linear regression of sigma_vals on log(v_vals) with weights.
    Returns slope, intercept, R^2."""
    log_v = np.log(v_vals.astype(np.float64))
    w = weights.astype(np.float64)
    w_sum = w.sum()
    if w_sum < 10:
        return float("nan"), float("nan"), float("nan")
    log_v_mean = (w * log_v).sum() / w_sum
    sigma_mean = (w * sigma_vals).sum() / w_sum
    log_v_c = log_v - log_v_mean
    sigma_c = sigma_vals - sigma_mean
    cov = (w * log_v_c * sigma_c).sum() / w_sum
    var = (w * log_v_c * log_v_c).sum() / w_sum
    if var < 1e-20:
        return float("nan"), float("nan"), float("nan")
    slope = cov / var
    intercept = sigma_mean - slope * log_v_mean
    pred = slope * log_v + intercept
    ss_res = (w * (sigma_vals - pred) ** 2).sum()
    ss_tot = (w * sigma_c * sigma_c).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(intercept), float(r2)


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    # Load sigma values from N=2^25 parquet
    N_data = 1 << 25
    print(f"[load] sigma values from main_N{N_data}.parquet", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{N_data}.parquet")
    sigma_lookup = df["sigma"].to_numpy().astype(np.int32)
    # sigma_lookup[n-1] = sigma(n), for n=1..N_data. (Parquet starts at n=1.)
    print(f"  sigma lookup table size: {len(sigma_lookup):,}", flush=True)

    # Sample random odd starts. We want trajectories that pass through small
    # values, so use starts uniformly in [3, N_max] for some N_max.
    rng = np.random.default_rng(42)
    N_max = 10**7
    n_starts = 1_000_000
    odd_starts = (2 * rng.integers(1, N_max // 2, size=n_starts, dtype=np.int64) + 1)
    print(f"[setup] {n_starts:,} odd starts in [3, {N_max:,}]", flush=True)

    for threshold in [5000, 11586]:  # 11586 ~ sqrt(2^27); 5000 user's framing
        print()
        print(f"=== Threshold v ≤ {threshold} ===", flush=True)
        t0 = time.perf_counter()
        visit_counts = collect_visit_counts(odd_starts, threshold)
        print(f"  walked trajectories in {time.perf_counter()-t0:.1f}s", flush=True)
        total_visits = int(visit_counts.sum())
        print(f"  total visits to v ≤ {threshold}: {total_visits:,}", flush=True)

        # Build arrays of (v, sigma(v), weight) for v=1..threshold with weight > 0
        v_arr = np.arange(1, threshold + 1, dtype=np.int64)
        sigma_arr = sigma_lookup[v_arr - 1]  # parquet rows start at n=1
        weights_arr = visit_counts[1:threshold + 1].astype(np.int64)
        valid = weights_arr > 0

        v_v = v_arr[valid]
        sigma_v = sigma_arr[valid].astype(np.float64)
        w_v = weights_arr[valid]

        odd_mask = (v_v % 2 == 1)
        even_mask = ~odd_mask

        # Slope: all visits
        slope_all, int_all, r2_all = fit_weighted_slope(v_v, sigma_v, w_v)
        slope_odd, int_odd, r2_odd = fit_weighted_slope(
            v_v[odd_mask], sigma_v[odd_mask], w_v[odd_mask])
        slope_even, int_even, r2_even = fit_weighted_slope(
            v_v[even_mask], sigma_v[even_mask], w_v[even_mask])

        print(f"  unique values visited: {valid.sum():,} of {threshold}")
        print(f"  odd values visited:    {odd_mask.sum():,}   "
              f"weight sum: {int(w_v[odd_mask].sum()):,}")
        print(f"  even values visited:   {even_mask.sum():,}   "
              f"weight sum: {int(w_v[even_mask].sum()):,}")
        # Mean weight per visited odd value vs per visited even value
        if odd_mask.sum() > 0 and even_mask.sum() > 0:
            print(f"  mean weight per odd value:  {w_v[odd_mask].sum() / odd_mask.sum():.1f}")
            print(f"  mean weight per even value: {w_v[even_mask].sum() / even_mask.sum():.1f}")

        print()
        print(f"  Trajectory-weighted regression of σ(v) on log(v):")
        print(f"    slope_all   = {slope_all:.4f}    intercept = {int_all:>+.4f}    R^2 = {r2_all:.4f}")
        print(f"    slope_odd   = {slope_odd:.4f}    intercept = {int_odd:>+.4f}    R^2 = {r2_odd:.4f}")
        print(f"    slope_even  = {slope_even:.4f}    intercept = {int_even:>+.4f}    R^2 = {r2_even:.4f}")
        print(f"    halving compensation = slope_all − slope_odd = {slope_all - slope_odd:>+.4f}")

        # Closed-form prediction for halving-chain compensation:
        # Per Syracuse step under Geom(1/2): 1 odd visit (start) + E[v]=2 even visits + 1 odd visit (next start, but counted next syracuse step).
        # So per Syracuse step: 1 odd + 2 even = 3 visits, with log values shifted.
        # Halving chain visits: at log(start_value/2^j) for j=0..k-1 (k = halving chain length).
        # Wait, the chain starts at 3m+1 (which has log = log(3m+1) ≈ log(3m)). Then halves.
        # Actually after odd→3m+1: visit value 3m+1. Then halve: visit (3m+1)/2. Etc.
        # For Syracuse step from m to m' = (3m+1)/2^k: chain values are 3m+1, (3m+1)/2, ..., (3m+1)/2^(k-1) = 2·m'.
        # log values: log(3m+1) ≈ log(3) + log(m), log(3) + log(m) - log(2), ..., log(3)+log(m)-(k-1)log(2).
        # The k-th value, (3m+1)/2^k = m', is the next odd start.
        print()
        print(f"  Closed-form check for halving-chain compensation:")
        print(f"  Under Geom(1/2): per Syracuse step, mean halving chain length = 2.")
        print(f"  Halving visits at log values shifted by {{0, log(2), 2·log(2)}} relative to log(3m+1).")
        print(f"  σ values at those visits decrease by 1 each (just halving steps).")
        print(f"  Whether this gives +{slope_all - slope_odd:.3f} requires the full integral over π.")


if __name__ == "__main__":
    main()
