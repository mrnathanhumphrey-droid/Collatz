"""
Experiment 43 — Long-lag autocovariance of v_2 along Collatz orbits.

Tests whether step-correlation structure across scales is the mechanism for
K_eff finite-N oscillation. Prior exp 29 found Σ_{k=1..20} Cov(v_0, v_k) =
−4.84 (anti-correlated short-range). This extends to long lags (up to 500)
and across N to test:

  - Does the autocovariance decay to 0 at large lag (short-range only) or
    plateau (long-range)?
  - Does the cumulative sum Σ_{k=1..L} Cov stabilize as L grows (finite
    total correction) or diverge?
  - Is the cumulative-sum value N-dependent? If yes, correlation structure
    is a candidate mechanism for K_eff(N) drift.

Method:
  - For each N ∈ {2²⁰, 2²⁵, 2²⁷}, sample 200K odd starts from the parquet.
  - For each orbit, build v_t = ν_2(3*m_t + 1) where m_t is the t-th odd
    value visited (Syracuse step view).
  - Pool pairs (v_t, v_{t+k}) across all orbits for k = 1, 5, 10, 25, 50,
    100, 250, 500. Compute Cov(v_t, v_{t+k}) = E[v_t · v_{t+k}] − E[v]².
  - Cumulative sum at L = 20, 50, 100, 250, 500.

Decisive outcomes:
  - cumsum stabilizes near same value across N → not the K_eff mechanism.
  - cumsum diverges or scales with N → correlation structure is the
    mechanism, predict K_eff(N) from cumsum.

Usage:
    python 43_v_autocovariance_long_lag.py
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")


@njit(cache=True)
def _v2(x):
    """Number of trailing zeros in x. Assumes x > 0."""
    c = 0
    while (x & 1) == 0:
        x >>= 1
        c += 1
    return c


@njit(parallel=True, cache=True)
def walk_and_accumulate_lag_stats(starts, lags, max_v_per_orbit):
    """For each odd start, walk Collatz and build v_2 sequence (Syracuse-step
    view). For each requested lag k, accumulate sum of v_t · v_{t+k} and
    count of pairs across all orbits.

    Also accumulate global E[v] sum/count.

    Args:
      starts: int64 array of odd starts.
      lags: int32 array of lag values (sorted ascending).
      max_v_per_orbit: maximum v_2 sequence length to allocate per orbit.

    Returns:
      sum_pairs[k_idx], cnt_pairs[k_idx]: pair statistics per lag
      v_sum, v_sum2, v_cnt: global v statistics
      n_orbits_completed: count of orbits that finished within max_steps
    """
    n = len(starts)
    n_lags = len(lags)
    n_chunks = 32
    chunk = (n + n_chunks - 1) // n_chunks

    sum_pairs_local = np.zeros((n_chunks, n_lags), dtype=np.float64)
    cnt_pairs_local = np.zeros((n_chunks, n_lags), dtype=np.int64)
    v_sum_local = np.zeros(n_chunks, dtype=np.float64)
    v_sum2_local = np.zeros(n_chunks, dtype=np.float64)
    v_cnt_local = np.zeros(n_chunks, dtype=np.int64)
    completed_local = np.zeros(n_chunks, dtype=np.int64)

    max_steps = 200_000

    for c in prange(n_chunks):
        lo = c * chunk
        hi = min(lo + chunk, n)
        v_buf = np.zeros(max_v_per_orbit, dtype=np.int32)
        for i in range(lo, hi):
            m = starts[i]
            steps = 0
            v_len = 0
            failed = False
            while m > 1 and steps < max_steps:
                if m & 1:
                    # Syracuse step: m -> (3m+1)/2^v
                    threex_p1 = 3 * m + 1
                    v = _v2(threex_p1)
                    if v_len < max_v_per_orbit:
                        v_buf[v_len] = v
                        v_len += 1
                    m = threex_p1 >> v
                    steps += 1 + v  # 1 odd step + v halvings
                else:
                    # shouldn't enter here from odd start, but guard
                    m = m >> 1
                    steps += 1
            if m != 1:
                failed = True
            if failed:
                continue
            completed_local[c] += 1

            # Accumulate v stats and pair stats
            for t in range(v_len):
                v_sum_local[c] += v_buf[t]
                v_sum2_local[c] += v_buf[t] * v_buf[t]
                v_cnt_local[c] += 1
            for k_idx in range(n_lags):
                k = lags[k_idx]
                if v_len > k:
                    for t in range(v_len - k):
                        sum_pairs_local[c, k_idx] += float(v_buf[t]) * float(v_buf[t + k])
                        cnt_pairs_local[c, k_idx] += 1

    sum_pairs = sum_pairs_local.sum(axis=0)
    cnt_pairs = cnt_pairs_local.sum(axis=0)
    v_sum = v_sum_local.sum()
    v_sum2 = v_sum2_local.sum()
    v_cnt = v_cnt_local.sum()
    completed = completed_local.sum()
    return sum_pairs, cnt_pairs, v_sum, v_sum2, v_cnt, completed


def autocov_at_N(data_dir, N, n_samples=200_000, seed=42):
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    n_total = len(n_arr)
    rng = np.random.default_rng(seed)
    if n_total > n_samples:
        idx = rng.choice(n_total, size=n_samples, replace=False)
        n_arr = n_arr[idx]
    n_used = len(n_arr)

    lags = np.array([1, 5, 10, 25, 50, 100, 250, 500], dtype=np.int32)
    # Allocate enough buffer; max σ_S at N=2^27 is ~150, but tail can reach 500+.
    # Be safe with 2000.
    max_v_per_orbit = 2000

    t0 = time.perf_counter()
    sum_pairs, cnt_pairs, v_sum, v_sum2, v_cnt, completed = \
        walk_and_accumulate_lag_stats(n_arr, lags, max_v_per_orbit)
    elapsed = time.perf_counter() - t0

    Ev = v_sum / v_cnt
    Var_v = v_sum2 / v_cnt - Ev * Ev
    cov_per_lag = np.zeros(len(lags))
    for k_idx in range(len(lags)):
        if cnt_pairs[k_idx] > 0:
            E_pair = sum_pairs[k_idx] / cnt_pairs[k_idx]
            cov_per_lag[k_idx] = E_pair - Ev * Ev
        else:
            cov_per_lag[k_idx] = float("nan")

    return {
        "N": N,
        "n_used": n_used,
        "completed": completed,
        "elapsed_s": elapsed,
        "Ev": Ev,
        "Var_v": Var_v,
        "v_cnt": v_cnt,
        "lags": lags,
        "cov_per_lag": cov_per_lag,
        "cnt_pairs_per_lag": cnt_pairs,
    }


@njit(parallel=True, cache=True)
def cumsum_via_full_lag_walk(starts, L_max):
    """Compute Σ_{k=1..L} Cov(v_0, v_k) for L = 20, 50, 100, 250, 500, by
    accumulating sum_{lag=1..L_max} (v_t · v_{t+lag}) per orbit.

    Returns sum_per_lag (length L_max+1, sum_per_lag[k] = sum over orbits of
    sum_{t} v_t · v_{t+k}), cnt_per_lag (number of pairs at each lag),
    and global v stats.
    """
    n = len(starts)
    n_chunks = 32
    chunk = (n + n_chunks - 1) // n_chunks
    L1 = L_max + 1
    sum_per_lag = np.zeros((n_chunks, L1), dtype=np.float64)
    cnt_per_lag = np.zeros((n_chunks, L1), dtype=np.int64)
    v_sum = np.zeros(n_chunks, dtype=np.float64)
    v_cnt = np.zeros(n_chunks, dtype=np.int64)

    max_steps = 200_000

    for c in prange(n_chunks):
        lo = c * chunk
        hi = min(lo + chunk, n)
        v_buf = np.zeros(2000, dtype=np.int32)
        for i in range(lo, hi):
            m = starts[i]
            steps = 0
            v_len = 0
            while m > 1 and steps < max_steps:
                if m & 1:
                    threex_p1 = 3 * m + 1
                    v_val = _v2(threex_p1)
                    if v_len < 2000:
                        v_buf[v_len] = v_val
                        v_len += 1
                    m = threex_p1 >> v_val
                    steps += 1 + v_val
                else:
                    m = m >> 1
                    steps += 1
            if m != 1:
                continue
            for t in range(v_len):
                v_sum[c] += v_buf[t]
                v_cnt[c] += 1
            # lag 0 is variance, but we want lags 1..L_max
            for k in range(1, min(L1, v_len)):
                for t in range(v_len - k):
                    sum_per_lag[c, k] += float(v_buf[t]) * float(v_buf[t + k])
                    cnt_per_lag[c, k] += 1

    return (sum_per_lag.sum(axis=0), cnt_per_lag.sum(axis=0),
            v_sum.sum(), v_cnt.sum())


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "outputs" / "exp43"
    out_dir.mkdir(parents=True, exist_ok=True)

    N_values = [1 << 20, 1 << 25, 1 << 27]
    lag_targets = np.array([1, 5, 10, 25, 50, 100, 250, 500])
    cumsum_Ls = [20, 50, 100, 250, 500]

    print(f"v_2 long-lag autocovariance across N")
    print(f"Lag targets: {lag_targets.tolist()}")
    print(f"Cumsum L values: {cumsum_Ls}")
    print()

    results_by_N = {}

    for N in N_values:
        log2N = int(round(np.log2(N)))
        print(f"=== N = 2^{log2N} = {N:,} ===", flush=True)
        # Sample-based per-lag autocovariance
        r = autocov_at_N(data_dir, N)
        results_by_N[N] = r
        print(f"  used {r['n_used']:,} odd starts, "
              f"completed {r['completed']:,} orbits in {r['elapsed_s']:.1f}s")
        print(f"  E[v] = {r['Ev']:.5f}, Var[v] = {r['Var_v']:.5f}, "
              f"v-count = {r['v_cnt']:,}")
        print(f"  Per-lag autocovariance:")
        print(f"    {'lag':>5}  {'pairs':>14}  {'Cov(v_t,v_t+k)':>15}")
        for k_idx, k in enumerate(r["lags"]):
            print(f"    {int(k):>5}  {r['cnt_pairs_per_lag'][k_idx]:>14,d}  "
                  f"{r['cov_per_lag'][k_idx]:>+15.6f}")

        # Full lag walk for cumsum analysis
        print(f"\n  computing cumsum via full lag walk (L_max=500)...", flush=True)
        df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
            (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
        n_arr_full = df["n"].to_numpy().astype(np.int64)
        rng = np.random.default_rng(42)
        if len(n_arr_full) > 200_000:
            idx = rng.choice(len(n_arr_full), size=200_000, replace=False)
            n_arr_full = n_arr_full[idx]
        t0 = time.perf_counter()
        sum_per_lag, cnt_per_lag, v_sum_g, v_cnt_g = \
            cumsum_via_full_lag_walk(n_arr_full, 500)
        elapsed = time.perf_counter() - t0
        Ev = v_sum_g / v_cnt_g
        cov_full = np.full(501, np.nan)
        for k in range(1, 501):
            if cnt_per_lag[k] > 0:
                cov_full[k] = sum_per_lag[k] / cnt_per_lag[k] - Ev * Ev
        print(f"    full walk in {elapsed:.1f}s")

        # Cumulative sums
        print(f"  Cumulative Σ_{{k=1..L}} Cov(v_0, v_k):")
        for L in cumsum_Ls:
            cum = 0.0
            n_finite = 0
            for k in range(1, L + 1):
                if not np.isnan(cov_full[k]):
                    cum += cov_full[k]
                    n_finite += 1
            n_pairs_at_L = cnt_per_lag[L] if L <= 500 else 0
            print(f"    L = {L:>3}: cumsum = {cum:>+9.4f}    "
                  f"(over {n_finite} finite-cov lags; "
                  f"pair count at L={L}: {n_pairs_at_L:,})")
        results_by_N[N]["cov_full"] = cov_full
        results_by_N[N]["cnt_per_lag"] = cnt_per_lag
        gc.collect()
        print()

    # Cross-N comparison table
    print("=== Cross-N comparison ===")
    header = f"  {'N':>8}  {'E[v]':>7}  {'Var[v]':>7}  " + "  ".join(
        [f"{'Cov[' + str(k) + ']':>12}" for k in lag_targets])
    print(header)
    for N in N_values:
        r = results_by_N[N]
        log2N = int(round(np.log2(N)))
        row = f"  2^{log2N:<6}  {r['Ev']:>7.4f}  {r['Var_v']:>7.4f}  "
        row += "  ".join([f"{r['cov_per_lag'][i]:>+12.5f}" for i in range(len(lag_targets))])
        print(row)

    print()
    print("=== Cumulative sum cross-N table ===")
    print(f"  {'N':>8}  " + "  ".join([f"{'cumsum L=' + str(L):>13}" for L in cumsum_Ls]))
    cumsums_by_N = {}
    for N in N_values:
        r = results_by_N[N]
        log2N = int(round(np.log2(N)))
        cov_full = r["cov_full"]
        cums = []
        for L in cumsum_Ls:
            cum = 0.0
            for k in range(1, L + 1):
                if not np.isnan(cov_full[k]):
                    cum += cov_full[k]
            cums.append(cum)
        cumsums_by_N[N] = cums
        print(f"  2^{log2N:<6}  " + "  ".join([f"{c:>+13.4f}" for c in cums]))

    # Verdict logic
    print()
    print("=== Verdict logic ===")
    cumsum_at_500 = np.array([cumsums_by_N[N][-1] for N in N_values])
    spread = cumsum_at_500.max() - cumsum_at_500.min()
    print(f"  cumsum at L=500: {cumsum_at_500}")
    print(f"  spread across N = {spread:.4f}")
    print(f"  K_eff spread across same N range (Result 6) = ~1.35 (full 2^20-2^27)")

    # Stability: does cumsum(L) plateau within a single N?
    for N in N_values:
        log2N = int(round(np.log2(N)))
        cums = cumsums_by_N[N]
        delta_late = abs(cums[-1] - cums[-2])  # 500 vs 250
        print(f"  N=2^{log2N}: |cumsum(500) - cumsum(250)| = {delta_late:.4f}   "
              f"(plateau if small)")

    # Plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Plot 1: autocov vs lag (3 N values overlaid, log-x)
        fig, ax = plt.subplots(figsize=(8, 5))
        for N in N_values:
            r = results_by_N[N]
            log2N = int(round(np.log2(N)))
            ax.plot(r["lags"], r["cov_per_lag"], "o-", label=f"N=2^{log2N}")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xscale("log")
        ax.set_xlabel("lag k")
        ax.set_ylabel("Cov(v_t, v_{t+k})")
        ax.set_title("Plot 1: autocovariance of v_2 vs lag")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot1_autocov_vs_lag.png", dpi=120)
        plt.close()

        # Plot 2: cumsum vs L
        fig, ax = plt.subplots(figsize=(8, 5))
        Ls_full = list(range(1, 501))
        for N in N_values:
            log2N = int(round(np.log2(N)))
            cov_full = results_by_N[N]["cov_full"]
            cum_arr = np.full(500, np.nan)
            running = 0.0
            for k in range(1, 501):
                if not np.isnan(cov_full[k]):
                    running += cov_full[k]
                cum_arr[k - 1] = running
            ax.plot(Ls_full, cum_arr, "-", label=f"N=2^{log2N}")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xlabel("L (upper limit)")
        ax.set_ylabel("Σ_{k=1..L} Cov(v_t, v_{t+k})")
        ax.set_title("Plot 2: cumulative autocovariance vs L")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot2_cumsum_vs_L.png", dpi=120)
        plt.close()
        print(f"  plots saved to {out_dir}")
    except Exception as e:
        print(f"  [warn] plotting failed: {e}")

    print(f"\n[done]")


if __name__ == "__main__":
    main()
