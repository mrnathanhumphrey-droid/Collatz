"""
P_q_given_j.py — does P(q | j) match across j?

If P(q | j=2) ≈ P(q | j=4) ≈ P(q | j=5): ⟨σ_S | j⟩ reduces cleanly to
w_q(q) — same q-distribution per j class, different (log m_anc - log m_j)
gives different mean σ_S but no structural q-dependence on j.

If P(q | j) differs across j: orbits absorbing at different attractors
have different q-distributions, meaning j carries information beyond q.
Third structural object beyond w_q(q).

Method: 5 seeds × 1M orbits at N=2^32. Stratify by j_attr, then compute
P(q-band | j) using same q-edges as Result 33.
"""
import sys
import time
import numpy as np
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
    for i in prange(n):
        m = starts[i]
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


def main():
    log2N = 32
    N = 1 << log2N
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000
    j_targets = [2, 4, 5]
    q_edges_quantile = [0.0, 0.25, 0.50, 0.75, 1.0]
    q_labels = [0.125, 0.375, 0.625, 0.875]

    print(f"# P(q-band | j) at N=2^{log2N}, {len(seeds)} seeds × {n_starts:,} orbits\n")

    # P(q | j) per seed
    P_q_given_j_seeds = {jt: {q: [] for q in q_labels} for jt in j_targets}

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
        sigma, j_arr = walk(starts)
        valid = sigma > 0
        sigma = sigma[valid]
        j_arr = j_arr[valid]

        # σ-quantile bands using GLOBAL distribution edges (same for all j to make P(q|j) comparable)
        edges = np.quantile(sigma, q_edges_quantile[1:-1])
        band_idx = np.digitize(sigma, edges)  # 0..3

        for jt in j_targets:
            mask_j = j_arr == jt
            n_j = int(mask_j.sum())
            if n_j == 0:
                continue
            for b_idx, q in enumerate(q_labels):
                P_qj = float(((band_idx == b_idx) & mask_j).sum()) / n_j
                P_q_given_j_seeds[jt][q].append(P_qj)

    # Aggregate
    print(f"# P(q-band | j) — same q-distribution across j → reduction holds")
    print(f"#   {'q-band':>8}  " + "  ".join(f"P(q|j={jt})" for jt in j_targets))
    for q in q_labels:
        line = f"#   {q:>8.3f}  "
        for jt in j_targets:
            arr = P_q_given_j_seeds[jt][q]
            mean = np.mean(arr)
            se = np.std(arr, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0
            line += f"{mean:>7.4f}±{se:.4f}  "
        print(line)

    # Pairwise spread
    print(f"\n# Pairwise spread of P(q | j) across j (per q-band):")
    for q in q_labels:
        means = {jt: np.mean(P_q_given_j_seeds[jt][q]) for jt in j_targets}
        spread = max(means.values()) - min(means.values())
        ses = {jt: np.std(P_q_given_j_seeds[jt][q], ddof=1)/np.sqrt(len(P_q_given_j_seeds[jt][q])) for jt in j_targets}
        # max z-score on pairwise differences
        max_z = 0
        for j1 in j_targets:
            for j2 in j_targets:
                if j1 < j2:
                    z = abs(means[j1] - means[j2]) / np.sqrt(ses[j1]**2 + ses[j2]**2)
                    if z > max_z:
                        max_z = z
        print(f"   q={q}: spread = {spread:.5f}, max pairwise |z| = {max_z:.2f}")

    # Verdict
    print(f"\n# Verdict:")
    print(f"#   If max |z| < 3 across all q: P(q|j) ≈ same across j → ⟨σ_S|j⟩ reduces cleanly to w_q(q)")
    print(f"#   If any max |z| > 5 with structural pattern: third object beyond w_q(q)")


if __name__ == "__main__":
    main()
