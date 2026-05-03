"""
inverse_tree_v_distribution.py — derive closed-form ⟨v | j⟩ via the per-j
v-distribution.

Hypothesis: P(v | j) is Esscher-tilted Geom(1/2) with j-specific tilt w_j.
Geom-tilt structure: P_w(v=k) = (1 - 2^{-(1+w)}) · 2^{-(1+w)(k-1)}
Equivalently: log P(v=k) = const - k · (1+w) · log(2). Linear in k with slope.

Test:
  1. Walk N=2^32 orbits, record per-orbit v-count histogram.
  2. Aggregate per-j: total count of v=k across orbits absorbing at j.
  3. Compute log P(v=k | j) vs k, fit linear, extract slope.
  4. Slope = -(1+w_j) · log(2), so w_j = -slope/log(2) - 1.
  5. Check w_j(j) for closed form.

Closed-form candidates to test:
  - w_j = log_2(P_j) for some structural quantity P_j
  - w_j = -1/log(m_j) · constant
  - w_j = -log(m_j)/log(2) · constant
  - w_j related to depth in inverse Collatz tree
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
def walk_with_v_hist(m_starts, max_steps=400000, v_max=20):
    n = len(m_starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    v_count = np.zeros((n, v_max + 1), dtype=np.int32)
    failed_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(m_starts[i])
        steps = 0
        j_attr = -1
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                failed_arr[i] = True
                break
            three_m_plus_1 = 3 * m + 1
            v = 0
            tmp = three_m_plus_1
            while (tmp & 1) == 0:
                tmp >>= 1
                v += 1
            new_m = tmp
            steps += 1
            if v <= v_max:
                v_count[i, v] += 1
            else:
                v_count[i, v_max] += 1  # bucket overflow
            if new_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = new_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
    return sigma_arr, j_arr, v_count, failed_arr


def fit_geom_tilt(v_hist, kmax=12):
    """Fit log P(v=k) vs k linear, return slope, intercept, w, and R²."""
    total = v_hist[1:kmax+1].sum()
    if total < 100:
        return None
    P = v_hist[1:kmax+1].astype(np.float64) / total
    valid = P > 0
    k_arr = np.arange(1, kmax+1)[valid]
    logP = np.log(P[valid])
    if len(logP) < 3:
        return None
    # Linear fit log P = a + b · k
    p = np.polyfit(k_arr, logP, 1)
    slope, intercept = p[0], p[1]
    pred = slope * k_arr + intercept
    R2 = 1 - np.sum((logP - pred)**2) / np.sum((logP - logP.mean())**2)
    # Geom-tilt: slope = -(1+w)·log 2, so w = -slope/log 2 - 1
    w = -slope / LOG2 - 1
    p_w = 1 - 2 ** (-(1 + w))
    E_v_w = 1 / p_w if p_w > 0 else np.inf
    return {
        'slope': slope, 'intercept': intercept, 'R2': R2,
        'w': w, 'p_w': p_w, 'E_v_predicted': E_v_w,
        'k_fit_range': (k_arr[0], k_arr[-1]),
        'P': P, 'k_arr': np.arange(1, kmax+1),
    }


def main():
    N = 1 << 32
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000
    v_max = 20
    j_targets = [2, 4, 5, 7, 8]

    print(f"# Per-j v-distribution structural analysis at N=2^32")
    print(f"# {len(seeds)} seeds × {n_starts:,} orbits, v tracked up to {v_max}")
    print()

    # Aggregate v-counts per j across all seeds
    j_counts = {jt: np.zeros(v_max + 1, dtype=np.int64) for jt in j_targets}
    j_orbit_counts = {jt: 0 for jt in j_targets}
    j_logm = {jt: [] for jt in j_targets}
    j_sigma = {jt: [] for jt in j_targets}
    overall_counts = np.zeros(v_max + 1, dtype=np.int64)

    t0 = time.perf_counter()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
        sigma, j_arr, v_count, failed = walk_with_v_hist(starts, v_max=v_max)
        valid = ~failed
        starts_v = starts[valid]
        sigma = sigma[valid]
        j = j_arr[valid]
        v_count = v_count[valid]
        log_m = np.log(starts_v.astype(np.float64))
        # Overall
        overall_counts += v_count.sum(axis=0)
        # Per j
        for jt in j_targets:
            mask = j == jt
            if mask.sum() < 10:
                continue
            j_counts[jt] += v_count[mask].sum(axis=0)
            j_orbit_counts[jt] += int(mask.sum())
            j_logm[jt].append(float(log_m[mask].mean()))
            j_sigma[jt].append(float(sigma[mask].mean()))
    elapsed = time.perf_counter() - t0
    print(f"# Compute: {elapsed:.1f}s\n")

    # Empirical P(v=k) overall (Geom(1/2) baseline)
    print(f"# Geom(1/2) baseline check (overall P(v=k) across all orbits):")
    total_steps = overall_counts.sum()
    print(f"#   {'k':>3}  {'P_emp(v=k)':>12}  {'2^(-k) target':>14}  {'gap':>10}")
    for k in range(1, 8):
        P_emp = overall_counts[k] / total_steps
        P_geom = 2.0 ** (-k)
        print(f"#   {k:>3}  {P_emp:>12.6f}  {P_geom:>14.6f}  {P_emp - P_geom:>+10.6f}")
    fit_overall = fit_geom_tilt(overall_counts, kmax=12)
    print(f"\n#   Geom-tilt fit: slope = {fit_overall['slope']:+.5f}, "
          f"w_overall = {fit_overall['w']:+.5f} (Geom(1/2) → w=0), "
          f"E[v]_predicted = {fit_overall['E_v_predicted']:.4f}, R² = {fit_overall['R2']:.6f}")

    # Per-j tilt extraction
    print(f"\n# Per-j Geom-tilt extraction (P(v|j) ≈ Geom_w_j):")
    print(f"#   {'j':>3}  {'orbits':>8}  {'⟨log m|j⟩':>10}  {'⟨σ_S|j⟩':>9}  {'slope':>8}  "
          f"{'w_j':>8}  {'E[v]_pred':>10}  {'E[v]_emp':>10}  {'R²':>7}")
    fits = {}
    for jt in j_targets:
        if j_orbit_counts[jt] < 50:
            continue
        fit = fit_geom_tilt(j_counts[jt], kmax=12)
        if fit is None:
            continue
        n_orbits = j_orbit_counts[jt]
        # Empirical ⟨v|j⟩ from v-counts
        Ev_emp = sum(k * j_counts[jt][k] for k in range(1, v_max + 1)) / j_counts[jt].sum()
        mlm = np.mean(j_logm[jt])
        ms = np.mean(j_sigma[jt])
        fits[jt] = {**fit, 'E_v_empirical': Ev_emp, 'orbits': n_orbits,
                    'mean_log_m': mlm, 'mean_sigma': ms}
        print(f"#   {jt:>3}  {n_orbits:>8,}  {mlm:>10.4f}  {ms:>9.3f}  "
              f"{fit['slope']:>+8.5f}  {fit['w']:>+8.5f}  {fit['E_v_predicted']:>10.4f}  "
              f"{Ev_emp:>10.4f}  {fit['R2']:>7.5f}")

    # Closed-form candidate tests
    print(f"\n# Test closed-form candidates for w_j(j):")
    j_list = sorted([j for j in fits.keys()])
    w_arr = np.array([fits[j]['w'] for j in j_list])
    j_arr_np = np.array(j_list)
    log_mj = np.array([np.log((4**j - 1) / 3) for j in j_list])
    log2_mj = log_mj / LOG2

    print(f"\n#  Empirical w_j: " + ", ".join(f"w_{j}={w:.5f}" for j, w in zip(j_list, w_arr)))

    # Candidate 1: w_j = a · log(m_j)
    if len(j_list) >= 2:
        print(f"\n#  Candidate: w_j = a · log(m_j) (linear in log m_j)")
        # Fit a
        a = np.dot(w_arr, log_mj) / np.dot(log_mj, log_mj)
        print(f"#    Best fit a = {a:+.5f}")
        for j, log_mj_val, w in zip(j_list, log_mj, w_arr):
            pred = a * log_mj_val
            print(f"#    j={j}: w_emp={w:+.5f}, pred={pred:+.5f}, gap={w-pred:+.5f}")

    # Candidate 2: w_j = -a / log(m_j)
    print(f"\n#  Candidate: w_j = -a / log(m_j)")
    inv_log_mj = 1.0 / log_mj
    a = -np.dot(w_arr, inv_log_mj) / np.dot(inv_log_mj, inv_log_mj)
    print(f"#    Best fit a = {a:+.5f}")
    for j, ilm, w in zip(j_list, inv_log_mj, w_arr):
        pred = -a * ilm
        print(f"#    j={j}: w_emp={w:+.5f}, pred={pred:+.5f}, gap={w-pred:+.5f}")

    # Candidate 3: w_j = -log_2(P_j)·c (related to absorption probability)
    print(f"\n#  Candidate: w_j = a + b · log_2(m_j)  (affine in log_2(m_j))")
    if len(j_list) >= 2:
        # OLS: w_j = a + b · log2(m_j)
        x = log2_mj
        b, a = np.polyfit(x, w_arr, 1)
        print(f"#    Fit: w_j = {a:+.5f} + {b:+.5f} · log_2(m_j)")
        for j, x_val, w in zip(j_list, x, w_arr):
            pred = a + b * x_val
            print(f"#    j={j} (log_2 m_j={x_val:.4f}): w_emp={w:+.5f}, pred={pred:+.5f}, gap={w-pred:+.5f}")

    # Candidate 4: simple 1/j scaling
    print(f"\n#  Candidate: w_j = -a/j or affine in 1/j")
    inv_j = 1.0 / j_arr_np
    if len(j_list) >= 2:
        b, a = np.polyfit(inv_j, w_arr, 1)
        print(f"#    Fit: w_j = {a:+.5f} + {b:+.5f}/j")
        for j, inv, w in zip(j_list, inv_j, w_arr):
            pred = a + b * inv
            print(f"#    j={j}: w_emp={w:+.5f}, pred={pred:+.5f}, gap={w-pred:+.5f}")

    # Candidate 5: w_j ∝ ⟨log m | j⟩ - log(m_j) (excess descent)
    print(f"\n#  Candidate: w_j = -a / (⟨log m|j⟩ - log m_j) (inverse of descent magnitude)")
    excess = np.array([fits[j]['mean_log_m'] - np.log((4**j - 1) / 3) for j in j_list])
    if len(j_list) >= 2:
        # w = a + b · 1/excess
        x = 1.0 / excess
        b, a = np.polyfit(x, w_arr, 1)
        print(f"#    Fit: w_j = {a:+.5f} + {b:+.5f}/(⟨log m|j⟩ - log m_j)")
        for j, e, w in zip(j_list, excess, w_arr):
            pred = a + b / e
            print(f"#    j={j}: excess={e:.4f}, w_emp={w:+.5f}, pred={pred:+.5f}, gap={w-pred:+.5f}")

    # SAVE per-v distributions
    rows = []
    for j in j_list:
        f = fits[j]
        for k in range(1, v_max + 1):
            P = j_counts[j][k] / j_counts[j].sum()
            rows.append({
                'j': j, 'v': k, 'count': int(j_counts[j][k]),
                'P_v_given_j': P, 'log_P': np.log(P) if P > 0 else None,
                'orbits_at_j': j_orbit_counts[j],
            })
    pl.DataFrame(rows).write_csv("experiments_output/inverse_tree_v_dist.csv")
    print(f"\n[save] experiments_output/inverse_tree_v_dist.csv")


if __name__ == "__main__":
    main()
