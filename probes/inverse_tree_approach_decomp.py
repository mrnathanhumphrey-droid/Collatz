"""
inverse_tree_approach_decomp.py — decompose W_j by FIRST-PASSAGE APPROACH
DIRECTION through log(m_j).

Tests the user's "from-below" hypothesis directly: classify each orbit
absorbing at m_j by whether the FIRST passage of log m_t through log(m_j)
is from above (descending crossing, log m_{t-1} > log m_j > log m_t) or
from below (ascending crossing, log m_{t-1} < log m_j < log m_t).

Then compute conditional ⟨W_j | approach⟩ for each class. If the from-below
class contributes the negative sign of W_4, the user's hypothesis is
structurally verified.

Method:
  Walk Syracuse from random odd m ∈ [1, 2^32], 5 seeds × 1M orbits.
  For each orbit absorbing at m_j (j ∈ {2, 4, 5}):
    - Track full trajectory (m_0, m_1, ..., m_τ = m_j)
    - Find first t with m_t · 2 > m_j (or some crossing event)
    - Classify approach direction
    - Record (j, approach_class, sigma_S, log m_start)

  Aggregate: ⟨W_j | approach⟩ per (j, approach_class).

Cap: < 250 lines.
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

# m_j values for j = 1..8
M_J = {1: 1, 2: 5, 3: 21, 4: 85, 5: 341, 6: 1365, 7: 5461, 8: 21845}


@njit(parallel=True, cache=True)
def walk_with_approach_classification(
    m_starts, target_m_j, max_steps=400000
):
    """For each orbit: walk Syracuse, classify approach direction at first
    passage through log(target_m_j).

    approach codes:
      0 = m_start ≤ m_j (start at or below target; first "crossing" is at start)
      1 = first crossing from above (descending): all m_t ≥ m_j until reach
      2 = first crossing from below (ascending): some m_t with m_{t-1} < m_j < m_t

    Also records sigma_S (total Syracuse step count to absorb) and the
    attractor j_attr.
    """
    n = len(m_starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    approach_arr = np.zeros(n, dtype=np.int8)  # 0=at-or-below, 1=above, 2=below-ascend
    failed_arr = np.zeros(n, dtype=np.bool_)
    M_j = np.int64(target_m_j)
    for i in prange(n):
        m = np.int64(m_starts[i])
        steps = 0
        j_attr = -1
        approach = 0  # default
        # Initial classification
        if m <= M_j:
            approach = 0  # start at or below (will be reclassified if orbit crosses upward)
        else:
            approach = 1  # start above; default expectation is descending crossing
        crossed = False  # has orbit crossed log(M_j) yet?
        prev_m = m
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
            # Detect crossing of log(M_j) in this step (m → new_m)
            if not crossed:
                if (m > M_j and new_m <= M_j):
                    approach = 1  # descending crossing first
                    crossed = True
                elif (m < M_j and new_m > M_j):
                    approach = 2  # ascending crossing first (from-below)
                    crossed = True
                elif (m == M_j) or (new_m == M_j):
                    # exact hit; decide by direction
                    if prev_m > M_j or m > M_j:
                        approach = 1
                    else:
                        approach = 2
                    crossed = True
            if new_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            prev_m = m
            m = new_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
        approach_arr[i] = approach
    return sigma_arr, j_arr, approach_arr, failed_arr


def analyze(N, seed, n_starts, j_targets=(2, 4, 5)):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1

    results = {}
    for jt in j_targets:
        m_j = M_J[jt]
        sigma_arr, j_arr, approach_arr, failed_arr = (
            walk_with_approach_classification(starts, m_j)
        )
        valid = ~failed_arr & (j_arr == jt)
        if valid.sum() < 10:
            continue
        log_m_j = np.log(float(m_j))
        log_m_start = np.log(starts[valid].astype(np.float64))
        sigma = sigma_arr[valid].astype(np.float64)
        approach = approach_arr[valid]
        # Per-approach decomposition
        per_class = {}
        for cls in [0, 1, 2]:
            mask_cls = approach == cls
            n_cls = int(mask_cls.sum())
            if n_cls < 5:
                per_class[cls] = {'n': n_cls, 'W': np.nan, 'mean_log_m': np.nan, 'mean_sigma': np.nan}
                continue
            mlm = float(log_m_start[mask_cls].mean())
            ms = float(sigma[mask_cls].mean())
            W_cls = ms - mlm / LOG43 - 1 + log_m_j / LOG43
            per_class[cls] = {
                'n': n_cls,
                'frac': n_cls / valid.sum(),
                'W': W_cls,
                'mean_log_m': mlm,
                'mean_sigma': ms,
            }
        # Overall
        mlm_all = float(log_m_start.mean())
        ms_all = float(sigma.mean())
        W_all = ms_all - mlm_all / LOG43 - 1 + log_m_j / LOG43
        results[jt] = {
            'n_total': int(valid.sum()),
            'W_total': W_all,
            'mean_log_m': mlm_all,
            'mean_sigma': ms_all,
            'per_class': per_class,
        }
    return results


def main():
    here = Path(__file__).resolve().parent
    out_dir = here / "experiments_output"

    log2N = 32
    N = 1 << log2N
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000

    print(f"# Approach-direction decomposition of W_j at N=2^{log2N}")
    print(f"# {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# Empirical W_j: W_2=+7.156, W_4=-4.755, W_5=+4.590\n")

    # Aggregate across seeds
    agg = {jt: {0: [], 1: [], 2: [], 'total': []} for jt in [2, 4, 5]}
    agg_frac = {jt: {0: [], 1: [], 2: []} for jt in [2, 4, 5]}
    for seed in seeds:
        t0 = time.perf_counter()
        res = analyze(N, seed, n_starts)
        elapsed = time.perf_counter() - t0
        line = f"  seed={seed:>5}  ({elapsed:.1f}s):  "
        for jt in [2, 4, 5]:
            r = res.get(jt, None)
            if r is None:
                continue
            agg[jt]['total'].append(r['W_total'])
            for cls in [0, 1, 2]:
                pc = r['per_class'].get(cls)
                if pc and not np.isnan(pc['W']):
                    agg[jt][cls].append(pc['W'])
                    agg_frac[jt][cls].append(pc.get('frac', 0))
            line += f"W_{jt}_tot={r['W_total']:>+6.3f}  "
        print(line)

    # Summary per j
    print("\n" + "=" * 80)
    print("Per-j approach decomposition")
    print("=" * 80)
    for jt in [2, 4, 5]:
        m_j = M_J[jt]
        log_m_j = np.log(float(m_j))
        print(f"\nj={jt}, m_j={m_j}, log(m_j)={log_m_j:.3f} nats")
        W_total_arr = np.array(agg[jt]['total'])
        print(f"  W_j_total (5-seed): {W_total_arr.mean():+.4f} ± {W_total_arr.std(ddof=1)/np.sqrt(5):.4f}")
        for cls, label in [(0, "start ≤ m_j (m_start ≤ {})".format(m_j)),
                           (1, "first crossing FROM ABOVE (descending)"),
                           (2, "first crossing FROM BELOW (ascending, user's hypothesis)")]:
            W_arr = np.array(agg[jt][cls])
            f_arr = np.array(agg_frac[jt][cls])
            if len(W_arr) > 0:
                W_mean = W_arr.mean()
                W_se = W_arr.std(ddof=1)/np.sqrt(len(W_arr)) if len(W_arr) > 1 else 0
                f_mean = f_arr.mean()
                contrib = W_mean * f_mean
                print(f"  class {cls} ({label}):")
                print(f"     fraction = {f_mean:.4f},  ⟨W | class⟩ = {W_mean:+.4f} ± {W_se:.4f},  contribution = {contrib:+.4f}")
            else:
                print(f"  class {cls} ({label}): empty")

    # Verdict
    print("\n" + "=" * 80)
    print("Verdict")
    print("=" * 80)
    for jt in [2, 4, 5]:
        f_below_arr = np.array(agg_frac[jt][2])
        W_below_arr = np.array(agg[jt][2])
        f_above_arr = np.array(agg_frac[jt][1])
        W_above_arr = np.array(agg[jt][1])
        if len(W_below_arr) > 0 and len(W_above_arr) > 0:
            f_below = f_below_arr.mean()
            f_above = f_above_arr.mean()
            W_below = W_below_arr.mean()
            W_above = W_above_arr.mean()
            print(f"  j={jt}:  from-above frac={f_above:.3f}, ⟨W|above⟩={W_above:+.3f};  "
                  f"from-below frac={f_below:.3f}, ⟨W|below⟩={W_below:+.3f}")
            if f_below > 0.05:
                print(f"          → from-below contributes {f_below:.1%} of orbits with conditional W={W_below:+.3f}")
                print(f"          → contribution to W_j_total: {f_below*W_below:+.3f}")


if __name__ == "__main__":
    main()
