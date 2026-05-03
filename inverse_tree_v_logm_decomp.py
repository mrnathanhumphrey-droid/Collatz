"""
inverse_tree_v_logm_decomp.py — decompose W_j by per-orbit ⟨v⟩ and log m_start
to identify the structural mechanism for the sign pattern (W_4 < 0 < W_5 < W_2).

Hypothesis: for orbits absorbing at attractor m_j, the conditional distributions
⟨v_avg | j⟩ and ⟨log m_start | j⟩ vary with j. Higher ⟨v_avg | j⟩ → faster
descent → smaller σ_S → smaller W_j (potentially negative).

Identity:
  σ_S(orbit) = (log m_start - log m_j) / (E[v_orbit]·log 2 - log 3) + 1 + finite-orbit residual

W_j = ⟨σ_S | j⟩ - ⟨log m_start | j⟩/log(4/3) - 1 + log(m_j)/log(4/3)
    ≈ ⟨log m_start | j⟩ · [1/(⟨v|j⟩·log2 - log3) - 1/log(4/3)] - log(m_j) · [same difference] + residual

If ⟨v|j⟩ = 2 exactly: 1/(2·log2 - log3) = 1/log(4/3), bracket = 0, W_j = residual only.
If ⟨v|j⟩ > 2: bracket < 0 (fastdescentcoef less than Wald), so W_j < (Wald baseline).
If ⟨v|j⟩ < 2: bracket > 0, W_j > Wald baseline.

Method:
  Walk orbits at N=2^32, 5 seeds × 1M orbits. For each orbit:
    σ_S, log m_start, j_attr, sum of v's along orbit (for E[v])
  Aggregate per j: ⟨v | j⟩, ⟨log m_start | j⟩, ⟨σ_S | j⟩.
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
def walk_with_v_track(m_starts, max_steps=400000):
    """For each orbit: walk Syracuse, record σ_S, j_attr, sum of v's."""
    n = len(m_starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    v_sum_arr = np.zeros(n, dtype=np.int64)
    failed_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(m_starts[i])
        steps = 0
        v_sum = 0
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
            v_sum += v
            if new_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = new_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
        v_sum_arr[i] = v_sum
    return sigma_arr, j_arr, v_sum_arr, failed_arr


def analyze(N, seed, n_starts, j_targets=(2, 4, 5)):
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
    sigma_arr, j_arr, v_sum_arr, failed_arr = walk_with_v_track(starts)
    valid = ~failed_arr
    log_m = np.log(starts[valid].astype(np.float64))
    sigma = sigma_arr[valid].astype(np.float64)
    v_sum = v_sum_arr[valid].astype(np.float64)
    v_avg = v_sum / np.maximum(sigma, 1.0)  # mean v per orbit
    j = j_arr[valid]
    n_total = int(valid.sum())

    # Overall (unconditional)
    overall = {
        'mean_v': float(v_avg.mean()),
        'mean_log_m': float(log_m.mean()),
        'mean_sigma': float(sigma.mean()),
    }

    # Per j
    per_j = {}
    for jt in j_targets:
        mask = (j == jt)
        n_j = int(mask.sum())
        if n_j < 10:
            continue
        m_j = (4**jt - 1) // 3
        log_m_j = np.log(float(m_j))
        ml_m = float(log_m[mask].mean())
        m_s = float(sigma[mask].mean())
        m_v = float(v_avg[mask].mean())
        # SDs
        sd_v = float(v_avg[mask].std(ddof=1))
        sd_log_m = float(log_m[mask].std(ddof=1))
        # W_j
        W_j = m_s - ml_m / LOG43 - 1 + log_m_j / LOG43
        # Predicted W_j if E[v] is the only deviation:
        # σ_S_pred(v) = (log m_start - log m_j) / (v·log2 - log3)
        denom = m_v * LOG2 - LOG3
        if denom > 0:
            sigma_pred_from_v = (ml_m - log_m_j) / denom
            W_pred_from_v = sigma_pred_from_v - ml_m / LOG43 - 1 + log_m_j / LOG43
        else:
            W_pred_from_v = np.nan
        per_j[jt] = {
            'n': n_j,
            'P_j': n_j / n_total,
            'mean_v': m_v,
            'sd_v': sd_v,
            'mean_log_m': ml_m,
            'sd_log_m': sd_log_m,
            'mean_sigma': m_s,
            'W_j': W_j,
            'W_pred_from_v': W_pred_from_v,
        }
    return overall, per_j


def main():
    here = Path(__file__).resolve().parent
    out_dir = here / "experiments_output"

    log2N = 32
    N = 1 << log2N
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000

    print(f"# v-and-log_m decomposition of W_j at N=2^{log2N}")
    print(f"# {len(seeds)} seeds × {n_starts:,} orbits")
    print(f"# Empirical W_j: W_2=+7.156, W_4=-4.755, W_5=+4.590")
    print(f"# K_h conversion: 1/log(4/3) = {1/LOG43:.4f}\n")

    overall_arr = []
    per_j_arr = {jt: [] for jt in [2, 4, 5]}
    for seed in seeds:
        t0 = time.perf_counter()
        overall, per_j = analyze(N, seed, n_starts)
        elapsed = time.perf_counter() - t0
        overall_arr.append(overall)
        for jt in [2, 4, 5]:
            if jt in per_j:
                per_j_arr[jt].append(per_j[jt])
        print(f"  seed={seed:>5}  ({elapsed:.1f}s):  "
              f"E[v]_all={overall['mean_v']:.4f}  "
              f"E[v|2]={per_j[2]['mean_v']:.4f}  "
              f"E[v|4]={per_j[4]['mean_v']:.4f}  "
              f"E[v|5]={per_j[5]['mean_v']:.4f}")

    # Aggregate
    print("\n" + "=" * 80)
    print("Per-j conditional means (5-seed bootstrap)")
    print("=" * 80)

    print(f"\n{'j':>3}  {'P(j)':>7}  {'⟨v|j⟩':>10}  {'⟨log m|j⟩':>12}  {'⟨σ_S|j⟩':>12}  "
          f"{'W_j_emp':>10}  {'W_j_predicted_from_⟨v|j⟩':>26}")
    for jt in [2, 4, 5]:
        arr = per_j_arr[jt]
        if not arr:
            continue
        P_j = np.array([r['P_j'] for r in arr]).mean()
        v_arr = np.array([r['mean_v'] for r in arr])
        lm_arr = np.array([r['mean_log_m'] for r in arr])
        sg_arr = np.array([r['mean_sigma'] for r in arr])
        W_arr = np.array([r['W_j'] for r in arr])
        Wpred_arr = np.array([r['W_pred_from_v'] for r in arr])
        v_mean = v_arr.mean(); v_se = v_arr.std(ddof=1)/np.sqrt(5)
        lm_mean = lm_arr.mean(); lm_se = lm_arr.std(ddof=1)/np.sqrt(5)
        sg_mean = sg_arr.mean(); sg_se = sg_arr.std(ddof=1)/np.sqrt(5)
        W_mean = W_arr.mean(); W_se = W_arr.std(ddof=1)/np.sqrt(5)
        Wpred_mean = Wpred_arr.mean()
        print(f"{jt:>3}  {P_j:>7.4f}  "
              f"{v_mean:>6.4f}±{v_se:.4f}  "
              f"{lm_mean:>7.4f}±{lm_se:.4f}  "
              f"{sg_mean:>7.4f}±{sg_se:.4f}  "
              f"{W_mean:>+6.3f}±{W_se:.3f}  "
              f"{Wpred_mean:>+22.3f}")

    # Cross-class comparison
    print(f"\n# Overall (across all orbits, no j-conditioning):")
    Eov_v = np.array([o['mean_v'] for o in overall_arr]).mean()
    Eov_lm = np.array([o['mean_log_m'] for o in overall_arr]).mean()
    Eov_sg = np.array([o['mean_sigma'] for o in overall_arr]).mean()
    print(f"   E[v]_all = {Eov_v:.4f}  (Geom(1/2) baseline = 2.0)")
    print(f"   E[log m]_all = {Eov_lm:.4f}  (uniform on [1, 2^{log2N}] expected ≈ {log2N*LOG2 - 1:.3f})")
    print(f"   E[σ_S]_all = {Eov_sg:.4f}")

    # Sign-flip analysis
    print(f"\n# Sign-flip analysis: which conditional moment drives W_4 < 0?")
    arr2 = per_j_arr[2]; arr4 = per_j_arr[4]; arr5 = per_j_arr[5]
    v2 = np.array([r['mean_v'] for r in arr2]).mean()
    v4 = np.array([r['mean_v'] for r in arr4]).mean()
    v5 = np.array([r['mean_v'] for r in arr5]).mean()
    lm2 = np.array([r['mean_log_m'] for r in arr2]).mean()
    lm4 = np.array([r['mean_log_m'] for r in arr4]).mean()
    lm5 = np.array([r['mean_log_m'] for r in arr5]).mean()
    print(f"   ⟨v|j=2⟩ = {v2:.4f}    ⟨v|j=4⟩ = {v4:.4f}    ⟨v|j=5⟩ = {v5:.4f}")
    print(f"   ⟨log m|j=2⟩ = {lm2:.4f}    ⟨log m|j=4⟩ = {lm4:.4f}    ⟨log m|j=5⟩ = {lm5:.4f}")
    print(f"   sign of W_j ↔ sign of (Wald deviation from log(4/3) descent rate)")
    print(f"   For W_4 < 0: ⟨v|j=4⟩ should be > 2 (fast descent) OR ⟨log m|j=4⟩ ≪ overall")


if __name__ == "__main__":
    main()
