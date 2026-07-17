"""
inverse_tree_growth.py — measure per-attractor inverse-tree growth rate λ_j
via full forward enumeration of odd m ∈ [1, N].

For each odd m in [1, N], walk Syracuse to absorption. Record (σ_S, j_attr).
Aggregate count(σ_S, j) = # of m with that σ_S absorbing at j. Fit linear
log(count(σ_S)) vs σ_S to extract λ_j = exp(slope).

Test:
  Spectral closure: ⟨σ_S | j⟩ ≈ log(N/m_j) / log(λ_j)?
  Equivalently: λ_j has clean form?

This bypasses BFS-with-N-cap issue (which culls excursions). Forward
enumeration gives exact ancestor count by depth without that bias.

Cap < 200 lines.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_to_absorption(starts, max_steps=400000):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0
        j_attr = -1
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            v = 0
            while (three_m & 1) == 0:
                three_m >>= 1
                v += 1
            steps += 1
            if three_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = three_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
    return sigma_arr, j_arr


def fit_growth(counts_per_d, j_label, log_N, log_mj):
    """log count(d) vs d. Fit linear in growth regime."""
    items = sorted(counts_per_d.items())
    d_arr = np.array([d for d, c in items if c >= 5])
    c_arr = np.array([c for d, c in items if c >= 5])
    if len(d_arr) < 5:
        return None
    log_c = np.log(c_arr)
    # Find peak depth (where counts are max)
    peak_d = d_arr[np.argmax(c_arr)]
    # Growth regime: d in [3, peak_d - 2] (rising part)
    growth_mask = (d_arr >= 3) & (d_arr < peak_d - 2)
    if growth_mask.sum() < 5:
        # Fit on lower half before peak
        growth_mask = (d_arr >= 3) & (d_arr <= peak_d)
    d_fit = d_arr[growth_mask]
    c_fit = log_c[growth_mask]
    if len(d_fit) < 3:
        return None
    p = np.polyfit(d_fit, c_fit, 1)
    slope, intercept = p
    pred = slope * d_fit + intercept
    R2 = 1 - np.sum((c_fit - pred)**2) / np.sum((c_fit - c_fit.mean())**2)
    lam_j = float(np.exp(slope))
    # Predicted ⟨σ_S | j⟩ via spectral: log(N/m_j) / log(λ_j)
    sigma_pred = (log_N - log_mj) / slope if slope > 0 else None
    return {
        'slope': float(slope), 'lambda_j': lam_j, 'R2': float(R2),
        'd_fit_min': int(d_fit[0]), 'd_fit_max': int(d_fit[-1]),
        'peak_d': int(peak_d), 'sigma_pred': sigma_pred,
    }


def main():
    N_values = [1<<22, 1<<24, 1<<26]  # 4.2M, 16.8M, 67M odd integers max
    j_targets = [2, 4, 5, 7, 8]

    print(f"# Per-attractor inverse-tree growth rate λ_j via forward enumeration")
    print(f"# Hypothesis: ⟨σ_S | j⟩ ≈ log(N/m_j) / log(λ_j)")
    print()

    for N in N_values:
        log_N = np.log(float(N))
        # Enumerate all odd m in [1, N]
        # For N=2^26 = 67M: 33M odd integers. ~100 steps avg. ~3.3 billion ops.
        # numba parallel handles it in seconds.
        n_odd = N // 2
        print(f"## N = 2^{int(np.log2(N))} = {N:,} (enumerating {n_odd:,} odd integers)")
        starts = np.arange(1, N, 2, dtype=np.int64)
        if len(starts) > 100_000_000:
            print(f"  too large, sampling 50M instead")
            rng = np.random.default_rng(42)
            starts = 2 * rng.integers(1, (N - 1) // 2, size=50_000_000, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma_arr, j_arr = walk_to_absorption(starts)
        elapsed = time.perf_counter() - t0
        print(f"  walked in {elapsed:.1f}s")

        # Bucket per j and depth
        for jt in j_targets:
            mask = (j_arr == jt)
            n_j = int(mask.sum())
            if n_j < 100:
                continue
            sigma_j = sigma_arr[mask]
            log_mj = np.log(float((4**jt - 1) // 3))
            sigma_mean = float(sigma_j.mean())
            sigma_sd = float(sigma_j.std(ddof=1))
            # Per-depth counts
            counts = defaultdict(int)
            for s in sigma_j:
                counts[int(s)] += 1
            fit = fit_growth(counts, jt, log_N, log_mj)
            mode = max(counts, key=counts.get)
            if fit is not None:
                lam = fit['lambda_j']
                sigma_spectral = fit['sigma_pred']
                R2 = fit['R2']
                print(f"  j={jt}: N={n_j:>9,}, ⟨σ_S|j⟩_emp = {sigma_mean:.3f} ± {sigma_sd/np.sqrt(n_j):.3f}, "
                      f"σ-mode = {mode}, λ_j = {lam:.5f} "
                      f"(slope={fit['slope']:+.5f}, fit d∈[{fit['d_fit_min']}..{fit['d_fit_max']}], R²={R2:.4f}), "
                      f"σ_pred = {sigma_spectral:.2f}")
            else:
                print(f"  j={jt}: N={n_j:>9,}, ⟨σ_S|j⟩_emp = {sigma_mean:.3f}, "
                      f"σ-mode = {mode}, fit FAILED")

    print()
    print("# Spectral closure check: does λ_j extracted at one N predict ⟨σ_S | j⟩ at larger N?")
    print("# If λ_j is structural and N-stable: yes. Check by comparing λ_j across N.")


if __name__ == "__main__":
    main()
