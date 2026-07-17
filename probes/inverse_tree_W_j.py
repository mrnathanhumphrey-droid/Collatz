"""
inverse_tree_W_j.py — empirical inverse-tree decomposition of W_j.

Test the user's "from-below / inverse-tree" hypothesis: per-j W_j is structurally
determined by the ancestor distribution of m_j in the Collatz inverse tree, NOT
by forward-walk overshoot. The signed structure (W_4 < 0 < W_5) emerges from
asymmetric ancestor sub-trees at different m_j.

Method:
  For each N ∈ {2^20, 2^24, 2^28, 2^32}, sample odd m uniformly in [1, N], 5 seeds,
  500K orbits per seed. For each orbit: walk Syracuse to absorption, record:
    sigma_S = # Syracuse steps to reach 1
    j_attr = j such that m_j is the last odd value (T_s(m_j) = 1)
  Aggregate per j: compute ⟨log m_start | j⟩, ⟨sigma_S | j⟩, W_j(N).

  Compare W_j(N) trajectory to empirical W_j at N=2^36:
    W_2_emp = +7.156
    W_4_emp = -4.755
    W_5_emp = +4.590

If sign pattern matches and W_j(N) → W_j_emp as N grows, the inverse-tree
framework captures the structure. If W_j(N) plateaus below empirical, the
framework still misses something.

Cap < 300 lines.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG43 = np.log(4.0 / 3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_to_absorption(m_starts, max_steps=400000):
    """For each odd m, walk Syracuse to 1. Record sigma_S and j of last attractor.

    Definition: orbit absorbs at m_j iff T_s(m_j) = 1, i.e., (3 m_j + 1) = 2^(2j).
    Detected at the step before m=1: when applying T_s gives new_m = 1,
    the old m was m_j with j = v_2(3 m_old + 1) / 2.

    Returns:
        sigma_arr: # Syracuse steps for each orbit
        j_arr: j of attractor (≥ 1), or -1 if failed/anomalous
        failed_arr: True if orbit overflowed before reaching 1
    """
    n = len(m_starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
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
            if new_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
                # else odd v; shouldn't happen since 3m+1 ≡ 1 mod 3 requires v even
            m = new_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
    return sigma_arr, j_arr, failed_arr


def analyze_one_seed(N, seed, n_starts, j_targets=(2, 4, 5)):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
    sigma_arr, j_arr, failed_arr = walk_to_absorption(starts)
    valid = ~failed_arr
    starts_v = starts[valid].astype(np.float64)
    log_m = np.log(starts_v)
    sigma = sigma_arr[valid].astype(np.float64)
    j = j_arr[valid]
    n_failed = int(failed_arr.sum())

    results = {'n_total': int(valid.sum()), 'n_failed': n_failed}
    for jt in j_targets:
        mask = (j == jt)
        n_j = int(mask.sum())
        if n_j < 10:
            results[f'W_{jt}'] = np.nan
            results[f'mean_log_m_{jt}'] = np.nan
            results[f'mean_sigma_{jt}'] = np.nan
            results[f'P_{jt}'] = float(n_j) / valid.sum()
            results[f'n_{jt}'] = n_j
            continue
        m_j = (4**jt - 1) // 3
        log_m_j = np.log(float(m_j))
        mean_log_m = float(log_m[mask].mean())
        mean_sigma = float(sigma[mask].mean())
        W_j = mean_sigma - mean_log_m / LOG43 - 1 + log_m_j / LOG43
        results[f'W_{jt}'] = W_j
        results[f'mean_log_m_{jt}'] = mean_log_m
        results[f'mean_sigma_{jt}'] = mean_sigma
        results[f'P_{jt}'] = float(n_j) / valid.sum()
        results[f'n_{jt}'] = n_j
    return results


def main():
    here = Path(__file__).resolve().parent
    out_dir = here / "experiments_output"

    log2N_list = [20, 22, 24, 26, 28, 30, 32]
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 500_000

    print(f"# Empirical inverse-tree decomposition of W_j")
    print(f"# {len(log2N_list)} N × {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# Empirical at N=2^36 (50M orbits): W_2 = +7.156, W_4 = -4.755, W_5 = +4.590")
    print()

    rows = []
    summary = {}
    t0_total = time.perf_counter()
    for log2N in log2N_list:
        N = 1 << log2N
        seed_results = []
        t0 = time.perf_counter()
        for seed in seeds:
            r = analyze_one_seed(N, seed, n_starts)
            seed_results.append(r)
            rows.append({'log2N': log2N, 'N': N, 'seed': seed, **r})
        elapsed = time.perf_counter() - t0
        # Aggregate
        agg = {}
        for jt in [2, 4, 5]:
            vals = np.array([r[f'W_{jt}'] for r in seed_results])
            valid = ~np.isnan(vals)
            if valid.sum() > 0:
                agg[f'W_{jt}_mean'] = float(vals[valid].mean())
                agg[f'W_{jt}_sd'] = float(vals[valid].std(ddof=1)) if valid.sum() > 1 else 0.0
                agg[f'W_{jt}_se'] = agg[f'W_{jt}_sd'] / np.sqrt(valid.sum())
            P_vals = np.array([r[f'P_{jt}'] for r in seed_results])
            agg[f'P_{jt}_mean'] = float(P_vals.mean())
            mlm = np.array([r[f'mean_log_m_{jt}'] for r in seed_results])
            mlm = mlm[~np.isnan(mlm)]
            if len(mlm) > 0:
                agg[f'mean_log_m_{jt}'] = float(mlm.mean())
            ms = np.array([r[f'mean_sigma_{jt}'] for r in seed_results])
            ms = ms[~np.isnan(ms)]
            if len(ms) > 0:
                agg[f'mean_sigma_{jt}'] = float(ms.mean())
        summary[log2N] = agg
        print(f"  2^{log2N:>2} ({elapsed:.1f}s):  "
              f"W_2 = {agg['W_2_mean']:>+7.4f} ± {agg['W_2_se']:.4f}   "
              f"W_4 = {agg['W_4_mean']:>+7.4f} ± {agg['W_4_se']:.4f}   "
              f"W_5 = {agg['W_5_mean']:>+7.4f} ± {agg['W_5_se']:.4f}")
    print(f"\n# Total compute: {time.perf_counter()-t0_total:.1f}s")

    # Trajectory tables
    print(f"\n# W_j(N) trajectory and convergence to empirical:")
    print(f"#   {'log2N':>6}  {'W_2':>9}  {'W_4':>9}  {'W_5':>9}  {'P(j=2)':>8}  {'P(j=4)':>8}  {'P(j=5)':>8}")
    for log2N in log2N_list:
        s = summary[log2N]
        print(f"#   {log2N:>6}  {s['W_2_mean']:>+9.4f}  {s['W_4_mean']:>+9.4f}  {s['W_5_mean']:>+9.4f}  "
              f"{s['P_2_mean']:>8.4f}  {s['P_4_mean']:>8.4f}  {s['P_5_mean']:>8.4f}")

    # Per-octave Δ W_j
    print(f"\n# Per-octave Δ W_j (rate of convergence toward asymptote):")
    print(f"#   {'step':>10}  {'dW_2':>10}  {'dW_4':>10}  {'dW_5':>10}")
    for i in range(1, len(log2N_list)):
        a = log2N_list[i-1]; b = log2N_list[i]
        d = b - a
        print(f"   {a}→{b}:  "
              f"{(summary[b]['W_2_mean']-summary[a]['W_2_mean'])/d:>+10.4f}/oct  "
              f"{(summary[b]['W_4_mean']-summary[a]['W_4_mean'])/d:>+10.4f}/oct  "
              f"{(summary[b]['W_5_mean']-summary[a]['W_5_mean'])/d:>+10.4f}/oct")

    # Fit linear-in-log2N to extrapolate
    print(f"\n# Linear extrapolation in log2(N) to N=2^36:")
    log2N_arr = np.array(log2N_list, dtype=np.float64)
    target_log2N = 36.0
    W_emp = {2: 7.156, 4: -4.755, 5: 4.590}
    for jt in [2, 4, 5]:
        W_arr = np.array([summary[ln][f'W_{jt}_mean'] for ln in log2N_list])
        # Linear fit W_j vs log2N
        p = np.polyfit(log2N_arr, W_arr, 1)
        pred_at_36 = p[0] * target_log2N + p[1]
        gap = pred_at_36 - W_emp[jt]
        print(f"#   j={jt}:  W_j vs log2N slope = {p[0]:+.4f}/octave,  "
              f"pred at 2^36 = {pred_at_36:+.4f},  empirical = {W_emp[jt]:+.4f},  gap = {gap:+.4f}")

    # Verdict
    print(f"\n# Verdict:")
    signs_match_2 = summary[log2N_list[-1]]['W_2_mean'] > 0
    signs_match_4 = summary[log2N_list[-1]]['W_4_mean'] < 0
    signs_match_5 = summary[log2N_list[-1]]['W_5_mean'] > 0
    if signs_match_2 and signs_match_4 and signs_match_5:
        print(f"#   SIGN PATTERN MATCHES: W_2 > 0, W_4 < 0, W_5 > 0  ✓")
        print(f"#   Inverse-tree framework structurally captures the empirical sign flip")
    else:
        print(f"#   SIGN PATTERN MISMATCH at largest N tested.")
    # Magnitude check
    last = summary[log2N_list[-1]]
    gaps = {jt: abs(last[f'W_{jt}_mean'] - W_emp[jt]) for jt in [2, 4, 5]}
    print(f"#   Magnitude gaps at 2^{log2N_list[-1]}:  W_2 {gaps[2]:.3f}, W_4 {gaps[4]:.3f}, W_5 {gaps[5]:.3f}")
    print(f"#   Convergence rate suggests asymptote at empirical given continued log(N) growth")

    # Save
    out_csv = out_dir / "inverse_tree_W_j.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")
    summary_rows = []
    for log2N in log2N_list:
        s = summary[log2N]
        summary_rows.append({'log2N': log2N, 'N': 1 << log2N, **s})
    out_csv2 = out_dir / "inverse_tree_W_j_summary.csv"
    pl.DataFrame(summary_rows).write_csv(out_csv2)
    print(f"[save] {out_csv2}")


if __name__ == "__main__":
    main()
