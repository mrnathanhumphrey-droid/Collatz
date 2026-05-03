"""
joint_q_j_analysis.py — test whether σ-quantile q and absorbing-attractor j
are functionally equivalent observables or genuinely 2D.

For each orbit, record (σ_S, j_attr, log m_start, ⟨v⟩_orbit). Stratify by
σ-quantile band (5 bands at midpoints {0.125, 0.375, 0.625, 0.875, 0.975}).
Within each (q, j) cell, compute count, ⟨v⟩, ⟨σ_S⟩.

Test: does P(j | q-band) concentrate (q ↔ j highly redundant) or distribute
broadly (genuinely 2D)?

Compute: 5 seeds × 1M orbits at N=2^32. ~10s.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk(starts, max_steps=400000):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    v_sum_arr = np.zeros(n, dtype=np.int64)
    for i in prange(n):
        m = starts[i]
        steps = 0
        j_attr = -1
        v_total = 0
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            v = 0
            while (three_m & 1) == 0:
                three_m >>= 1
                v += 1
            steps += 1
            v_total += v
            if three_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = three_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
        v_sum_arr[i] = v_total
    return sigma_arr, j_arr, v_sum_arr


def main():
    log2N = 32
    N = 1 << log2N
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000
    j_targets = [2, 4, 5]
    q_edges = [0.0, 0.25, 0.50, 0.75, 1.0]
    q_labels = [0.125, 0.375, 0.625, 0.875]

    print(f"# Joint (q-band, j) analysis at N=2^{log2N}")
    print(f"# {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# σ-quantile bands at midpoints {q_labels}\n")

    # Aggregate across seeds
    cell = {(q, j): {'n': 0, 'v_sum_orbit_sum': 0.0, 'sigma_sum': 0.0, 'logm_sum': 0.0,
                      'sigma_sum_orbits': 0, 'v_count_orbits': 0}
            for q in q_labels for j in j_targets}

    P_j_per_q = {q: {j: [] for j in j_targets} for q in q_labels}
    Ev_per_q = {q: [] for q in q_labels}  # bulk mean v in band
    Ev_per_qj = {q: {j: [] for j in j_targets} for q in q_labels}

    t0_total = time.perf_counter()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, j_arr, v_sum = walk(starts)
        elapsed = time.perf_counter() - t0
        valid = sigma > 0
        sigma = sigma[valid]
        j_arr = j_arr[valid]
        v_sum = v_sum[valid]
        starts_v = starts[valid]
        log_m = np.log(starts_v.astype(np.float64))
        v_avg = v_sum.astype(np.float64) / np.maximum(sigma, 1).astype(np.float64)

        # σ-quantile bands
        edges = np.quantile(sigma, q_edges[1:-1])  # 0.25, 0.50, 0.75
        band_idx = np.digitize(sigma, edges)  # 0..3
        # 0 -> q=0.125, 1 -> 0.375, 2 -> 0.625, 3 -> 0.875

        # Per (band, j)
        for b_idx, q in enumerate(q_labels):
            mask_b = band_idx == b_idx
            n_b = int(mask_b.sum())
            if n_b == 0:
                continue
            # Bulk band stats
            Ev_per_q[q].append(float(v_avg[mask_b].mean()))
            for j_target in j_targets:
                mask_bj = mask_b & (j_arr == j_target)
                n_bj = int(mask_bj.sum())
                P_j_per_q[q][j_target].append(n_bj / n_b)
                if n_bj >= 5:
                    Ev_per_qj[q][j_target].append(float(v_avg[mask_bj].mean()))

        print(f"  seed={seed:>5} ({elapsed:.1f}s)  valid={valid.sum():,}")

    print(f"\n# Total: {time.perf_counter()-t0_total:.1f}s\n")

    # P(j | q-band) table
    print(f"# P(j | q-band) — does σ-quantile concentrate the j-distribution?")
    print(f"#   {'q':>6}  " + "  ".join(f"P(j={j}|q)" for j in j_targets) + "  ⟨v|q⟩_bulk")
    for q in q_labels:
        line = f"#   {q:>6.3f}  "
        for j in j_targets:
            P = np.mean(P_j_per_q[q][j])
            line += f"{P:>9.4f}  "
        Ev = np.mean(Ev_per_q[q])
        line += f"{Ev:>9.4f}"
        print(line)

    # ⟨v | q, j⟩ table
    print(f"\n# ⟨v | q, j⟩ joint")
    print(f"#   {'q':>6}  " + "  ".join(f"v|j={j}" for j in j_targets))
    for q in q_labels:
        line = f"#   {q:>6.3f}  "
        for j in j_targets:
            arr = Ev_per_qj[q][j]
            if arr:
                Ev = np.mean(arr)
                line += f"{Ev:>7.4f}  "
            else:
                line += f"{'--':>7}  "
        print(line)

    # Reduction test: is ⟨v | q, j⟩ ≈ ⟨v | j⟩ alone (j-only conditioning sufficient)?
    print(f"\n# Reduction test: does conditioning on j make q irrelevant?")
    print(f"# (Compare ⟨v | q, j⟩ across q for fixed j; if constant in q, then j alone suffices)")
    for j in j_targets:
        vals = []
        for q in q_labels:
            arr = Ev_per_qj[q][j]
            if arr:
                vals.append((q, np.mean(arr)))
        if len(vals) >= 2:
            v_arr = np.array([v for q, v in vals])
            spread = v_arr.max() - v_arr.min()
            print(f"  j={j}: ⟨v | q, j⟩ across q: {[f'{v:.3f}' for q, v in vals]}, spread = {spread:.4f}")

    # Reverse: does conditioning on q make j irrelevant?
    print(f"\n# Reverse reduction: does conditioning on q make j irrelevant?")
    print(f"# (Compare ⟨v | q, j⟩ across j for fixed q; if constant in j, then q alone suffices)")
    for q in q_labels:
        vals = []
        for j in j_targets:
            arr = Ev_per_qj[q][j]
            if arr:
                vals.append((j, np.mean(arr)))
        if len(vals) >= 2:
            v_arr = np.array([v for j, v in vals])
            spread = v_arr.max() - v_arr.min()
            print(f"  q={q}: ⟨v | q, j⟩ across j: {[f'(j={j},v={v:.3f})' for j, v in vals]}, spread = {spread:.4f}")

    # Verdict
    print(f"\n# Verdict:")
    print(f"#   Two observables functionally equivalent if: ⟨v|q,j⟩ depends on only one of (q, j)")
    print(f"#   Genuinely 2D if: ⟨v|q,j⟩ varies with both")


if __name__ == "__main__":
    main()
