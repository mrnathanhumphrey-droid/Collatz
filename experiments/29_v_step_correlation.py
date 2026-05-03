"""
Experiment 29 — step-to-step correlations of v along Syracuse trajectories.

The K(E[v]) formula assumes i.i.d. v values along trajectories. If
consecutive v's are correlated, the effective drift differs from
log(2)*E[v] - log(3) and per-octave beta_local picks up a residual
deviation from K_pred(E[v]) that the mean alone cannot explain.

Test: walk many trajectories, record (v_t, v_{t+1}) pairs, compute:
  - Lag-1 autocorrelation rho_1 = Cov(v_t, v_{t+1}) / Var(v_t).
  - Conditional E[v_{t+1} | v_t = k] for k = 1, 2, ..., 6.
  - Joint P(v_t, v_{t+1}) heatmap on small grid.

Expected under i.i.d.: rho_1 = 0; conditional E[v_{t+1} | v_t = k] = E[v]
independent of k. Departures from these are the i.i.d.-violation signature.

Usage:
    python 29_v_step_correlation.py --N_start 1000000 --T 200
"""
import argparse
import sys
import time

import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")


@njit(parallel=True, cache=True)
def collect_pair_stats(starts, T, max_v):
    """Chunk-parallel: record (v_t, v_{t+1}) joint histogram and pooled v sums.

    Returns:
      joint:  (max_v+1) x (max_v+1) joint histogram of (v_t, v_{t+1})
      sum_v:  sum of v_t across all transitions counted in joint
      sum_v2: sum of v_t^2 across all transitions counted in joint
      cnt:    number of transitions counted
    """
    n = len(starts)
    chunk = max(1, n // 256)
    n_chunks = (n + chunk - 1) // chunk
    joint = np.zeros((n_chunks, max_v + 1, max_v + 1), dtype=np.int64)
    sum_v = np.zeros(n_chunks, dtype=np.float64)
    sum_v2 = np.zeros(n_chunks, dtype=np.float64)
    sum_vt1 = np.zeros(n_chunks, dtype=np.float64)
    sum_vt1_2 = np.zeros(n_chunks, dtype=np.float64)
    sum_v_vt1 = np.zeros(n_chunks, dtype=np.float64)
    cnt = np.zeros(n_chunks, dtype=np.int64)

    for c in prange(n_chunks):
        start_i = c * chunk
        end_i = min(start_i + chunk, n)
        for i in range(start_i, end_i):
            m = starts[i]
            prev_v = -1
            for _ in range(T):
                if m == 1:
                    break
                tmp = 3 * m + 1
                v = 0
                while tmp % 2 == 0:
                    tmp //= 2
                    v += 1
                if prev_v >= 0 and prev_v <= max_v and v <= max_v:
                    joint[c, prev_v, v] += 1
                    sum_v[c] += prev_v
                    sum_v2[c] += prev_v * prev_v
                    sum_vt1[c] += v
                    sum_vt1_2[c] += v * v
                    sum_v_vt1[c] += prev_v * v
                    cnt[c] += 1
                prev_v = v
                m = tmp

    return (joint.sum(axis=0),
            sum_v.sum(), sum_v2.sum(),
            sum_vt1.sum(), sum_vt1_2.sum(),
            sum_v_vt1.sum(),
            cnt.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N_start", type=int, default=1_000_000)
    ap.add_argument("--T", type=int, default=200)
    ap.add_argument("--N_max", type=int, default=10**9)
    ap.add_argument("--max_v", type=int, default=15)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    starts = (2 * rng.integers(1, args.N_max // 2, size=args.N_start, dtype=np.int64) + 1)

    print(f"[setup] N_start = {args.N_start:,}  T = {args.T}  N_max = {args.N_max:,}",
          flush=True)
    print(f"[run]   collecting (v_t, v_t+1) pairs ...", flush=True)
    t0 = time.perf_counter()
    joint, sv, sv2, svt1, svt1_2, svvt1, cnt = collect_pair_stats(
        starts, args.T, args.max_v)
    t1 = time.perf_counter()
    print(f"[run]   done in {t1 - t0:.1f}s   pairs counted: {cnt:,}",
          flush=True)

    n = float(cnt)
    mean_v = sv / n
    mean_vt1 = svt1 / n
    var_v = sv2 / n - mean_v * mean_v
    var_vt1 = svt1_2 / n - mean_vt1 * mean_vt1
    cov_v_vt1 = svvt1 / n - mean_v * mean_vt1
    rho_1 = cov_v_vt1 / np.sqrt(var_v * var_vt1)

    print()
    print(f"=== Pooled v-step statistics (lag 1) ===")
    print(f"  E[v_t]      = {mean_v:.6f}")
    print(f"  E[v_{{t+1}}]    = {mean_vt1:.6f}")
    print(f"  Var(v_t)    = {var_v:.6f}")
    print(f"  Var(v_t+1)  = {var_vt1:.6f}")
    print(f"  Cov(v_t, v_t+1) = {cov_v_vt1:.6f}")
    print(f"  rho_1       = {rho_1:.6f}    "
          f"(under i.i.d. expected 0; SE ~ {1/np.sqrt(n):.2e})")
    print(f"  (Geom(1/2): E[v]=2, Var=2)")

    print()
    print(f"=== Conditional E[v_{{t+1}} | v_t = k] ===")
    print(f"  {'v_t = k':>9} {'count':>12} {'E[v_t+1|k]':>12} {'minus marginal':>15}")
    for k in range(1, args.max_v + 1):
        row_count = joint[k, :].sum()
        if row_count < 100:
            continue
        cond_mean = sum(j * joint[k, j] for j in range(args.max_v + 1)) / row_count
        diff = cond_mean - mean_vt1
        print(f"  {k:>9} {int(row_count):>12,} {cond_mean:>12.5f} "
              f"{diff:>+15.5f}")

    print()
    print(f"=== Joint P(v_t = i, v_{{t+1}} = j) heatmap (top-left 8x8) ===")
    print(f"  v_t \\ v_t+1   ", end="")
    for j in range(1, 9):
        print(f"{j:>9}", end="")
    print()
    print(f"  marg_indep    ", end="")
    marg_v = joint.sum(axis=1) / n
    marg_vt1 = joint.sum(axis=0) / n
    for j in range(1, 9):
        # under independence, P(v_t=i, v_t+1=j) = P(v_t=i) * P(v_t+1=j); we
        # show the full heatmap of empirical / independent
        pass
    # Print empirical / independent ratio heatmap
    print()
    print(f"=== Empirical / independent ratio (top-left 8x8) ===")
    print(f"  i.i.d. would give all entries = 1.0")
    print(f"  v_t \\ v_t+1   ", end="")
    for j in range(1, 9):
        print(f"{j:>9}", end="")
    print()
    for i in range(1, 9):
        print(f"  v_t={i:>2}        ", end="")
        for j in range(1, 9):
            empirical = joint[i, j] / n
            independent = marg_v[i] * marg_vt1[j]
            ratio = empirical / independent if independent > 0 else float("nan")
            print(f"{ratio:>9.4f}", end="")
        print()

    print()
    print(f"=== Implication for K(E[v]) deviation ===")
    print(f"  Under i.i.d., the random-walk drift mean is exactly")
    print(f"    E[delta_log m per step] = log(3) - log(2) * E[v]")
    print(f"  Step correlations preserve this MEAN (each step's marginal v has")
    print(f"  the same expectation). They affect only the VARIANCE of the partial")
    print(f"  sum, which doesn't enter K. So even if rho_1 != 0, K(E[v]) is")
    print(f"  unchanged at the mean level.")
    print()
    print(f"  *Conditional* E[v_t+1 | v_t = k] varying with k indicates that the")
    print(f"  trajectory's effective sampling distribution differs from the marginal,")
    print(f"  but the long-run mean v is still the marginal E[v]. So this CANNOT")
    print(f"  be the mechanism for the residual beta_local variation.")
    print(f"  The result here is therefore a *negative* test: if rho_1 != 0 and")
    print(f"  conditional means vary, that's expected (Syracuse is deterministic);")
    print(f"  it's not the mechanism we're after.")


if __name__ == "__main__":
    main()
