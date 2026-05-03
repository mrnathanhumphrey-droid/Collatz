"""
Experiment 27 — beta-oscillation diagnostic.

Experiment 26 found beta = sigma vs ln(n) slope is non-monotone in N:
approaches 10.4282 from below for N <= 2^26, jumps above between 2^26 and
2^27, peaks at 2^28, drifts back below for 2^29 through 2^32 with the gap
re-opening to ~+0.01 by 2^32. The jump 2^26 -> 2^27 is the largest single
step.

Three diagnostics on the same sigma cache:

  (A) Cumulative beta on odd n in [3, N] (matches exp 26 output).
  (B) Per-octave local beta: for each j with 2^j <= N, fit beta on odd n
      in [2^j, 2^(j+1)]. Tells us whether local beta systematically shifts
      at the octave 2^26 -> 2^27 or whether oscillation is a weight-shift
      artifact of cumulative averaging across octaves.
  (C) Drop-top-K cumulative beta: refit beta after excluding the top K
      sigma values. If a few extreme outliers in [2^26, 2^27] are pulling
      cumulative beta up, K=100 or K=1000 should pull it back toward
      the smaller-N trend.

For (B), local beta per octave should be ~10.4282 if the heuristic holds at
each scale. Departures localize where the deviation lives.

For (C), if dropping top-K stabilizes beta, outliers are the mechanism.
If beta is unchanged, the deviation is bulk, not tail-driven.

Usage:
    python 27_beta_oscillation_diagnostic.py --N 134217728
    python 27_beta_oscillation_diagnostic.py --N 268435456 --drop_top 0,100,1000,10000
"""
import argparse
import sys
import time

import numpy as np
import numba
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

PATH_MAX = 200_000


@njit(cache=True, boundscheck=False)
def compute_sigma_cache(N):
    sigma = np.full(N + 1, -1, dtype=np.int32)
    sigma[1] = 0
    path = np.empty(PATH_MAX, dtype=np.int64)

    for n in range(2, N + 1):
        if sigma[n] != -1:
            continue
        m = np.int64(n)
        path_len = 0
        overflow = False
        while True:
            if m <= N and sigma[m] != -1:
                break
            if path_len >= PATH_MAX:
                overflow = True
                break
            path[path_len] = m
            path_len += 1
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            if m == 1:
                break

        if overflow:
            m = np.int64(n)
            steps = 0
            while m > 1:
                if m & 1:
                    m = 3 * m + 1
                else:
                    m = m >> 1
                steps += 1
            sigma[n] = np.int32(steps)
            continue

        cur = np.int32(0) if m == 1 else sigma[m]
        for i in range(path_len - 1, -1, -1):
            v = path[i]
            cur += 1
            if v <= N:
                sigma[v] = cur
    return sigma


@njit(cache=True, parallel=True, boundscheck=False)
def streaming_ols_range(sigma, lo, hi):
    """OLS sums over odd n in [max(lo,3), hi)."""
    n_threads = 64
    span = hi - lo
    chunk = max(2, span // n_threads)
    n_chunks = (span + chunk - 1) // chunk

    sx = np.zeros(n_chunks, dtype=np.float64)
    sy = np.zeros(n_chunks, dtype=np.float64)
    sxy = np.zeros(n_chunks, dtype=np.float64)
    sx2 = np.zeros(n_chunks, dtype=np.float64)
    cnt = np.zeros(n_chunks, dtype=np.int64)

    for c in prange(n_chunks):
        a = lo + c * chunk
        b = min(a + chunk, hi)
        if a < 3:
            a = 3
        if a % 2 == 0:
            a += 1
        for n in range(a, b, 2):
            x = np.log(np.float64(n))
            y = np.float64(sigma[n])
            sx[c] += x
            sy[c] += y
            sxy[c] += x * y
            sx2[c] += x * x
            cnt[c] += 1

    return sx.sum(), sy.sum(), sxy.sum(), sx2.sum(), cnt.sum()


def fit_beta(sx, sy, sxy, sx2, cnt):
    n_obs = float(cnt)
    denom = n_obs * sx2 - sx * sx
    beta = (n_obs * sxy - sx * sy) / denom
    alpha = (sy - beta * sx) / n_obs
    return alpha, beta, int(cnt)


@njit(cache=True, parallel=True, boundscheck=False)
def streaming_ols_drop_top(sigma, N, threshold):
    """OLS over odd n in [3, N] excluding sigma >= threshold."""
    n_threads = 64
    chunk = max(2, (N // n_threads) | 1)
    n_chunks = (N + chunk - 1) // chunk

    sx = np.zeros(n_chunks, dtype=np.float64)
    sy = np.zeros(n_chunks, dtype=np.float64)
    sxy = np.zeros(n_chunks, dtype=np.float64)
    sx2 = np.zeros(n_chunks, dtype=np.float64)
    cnt = np.zeros(n_chunks, dtype=np.int64)
    excluded = np.zeros(n_chunks, dtype=np.int64)

    for c in prange(n_chunks):
        lo = c * chunk
        hi = min(lo + chunk, N + 1)
        if lo < 3:
            lo = 3
        if lo % 2 == 0:
            lo += 1
        for n in range(lo, hi, 2):
            s = sigma[n]
            if s >= threshold:
                excluded[c] += 1
                continue
            x = np.log(np.float64(n))
            y = np.float64(s)
            sx[c] += x
            sy[c] += y
            sxy[c] += x * y
            sx2[c] += x * x
            cnt[c] += 1

    return sx.sum(), sy.sum(), sxy.sum(), sx2.sum(), cnt.sum(), excluded.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--drop_K", type=str, default="0,10,100,1000,10000",
                    help="Drop top-K sigma values (comma-separated).")
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    numba.set_num_threads(args.threads)

    print(f"[setup] N = {args.N:,}  threads = {numba.get_num_threads()}",
          flush=True)
    print(f"[mem]   sigma cache: {(args.N + 1) * 4 / 1e9:.2f} GB", flush=True)

    t0 = time.perf_counter()
    sigma = compute_sigma_cache(args.N)
    t1 = time.perf_counter()
    max_sigma = int(sigma[1:].max())
    print(f"[run]   cache built in {t1 - t0:.1f}s   max(sigma) = {max_sigma}",
          flush=True)

    heuristic = 3.0 / (np.log(4.0) - np.log(3.0))

    # (A) Cumulative beta in [3, N]
    print()
    print("=" * 72)
    print("(A) Cumulative beta on odd n in [3, N]")
    print("=" * 72)
    sx, sy, sxy, sx2, cnt = streaming_ols_range(sigma, 0, args.N + 1)
    a, b, n_obs = fit_beta(sx, sy, sxy, sx2, cnt)
    gap_cum = heuristic - b
    print(f"  n_obs = {n_obs:,}   beta = {b:.6f}   gap = {gap_cum:>+.6f}")

    # (B) Per-octave local beta
    print()
    print("=" * 72)
    print("(B) Per-octave local beta on odd n in [2^j, 2^(j+1))")
    print("=" * 72)
    print(f"  {'j':>3} {'lo':>14} {'hi':>14} {'n_obs':>13} "
          f"{'beta_local':>11} {'gap':>9} {'mean(sigma)':>12} {'max(sigma)':>11}")
    j = 1
    octave_data = []
    while (1 << (j + 1)) <= args.N + 1:
        lo = 1 << j
        hi = 1 << (j + 1)
        sx, sy, sxy, sx2, cnt = streaming_ols_range(sigma, lo, hi)
        if cnt > 100:
            a_loc, b_loc, n_obs_loc = fit_beta(sx, sy, sxy, sx2, cnt)
            gap_loc = heuristic - b_loc
            mean_sig = sy / cnt
            # Compute local max sigma cheaply
            max_sig = int(sigma[lo:hi].max())
            print(f"  {j:>3} {lo:>14,} {hi - 1:>14,} {n_obs_loc:>13,} "
                  f"{b_loc:>11.4f} {gap_loc:>+9.4f} "
                  f"{mean_sig:>12.4f} {max_sig:>11}")
            octave_data.append((j, lo, hi, n_obs_loc, b_loc, gap_loc, mean_sig, max_sig))
        j += 1

    # (C) Drop top-K cumulative beta
    print()
    print("=" * 72)
    print("(C) Cumulative beta excluding top-K sigma values")
    print("=" * 72)

    Ks = sorted({int(k) for k in args.drop_K.split(",")})
    if 0 not in Ks:
        Ks = [0] + Ks

    # Need top-K threshold. Sort sigma values once.
    print("  sorting sigma to find top-K thresholds ...", flush=True)
    t0 = time.perf_counter()
    odd_mask_check = (np.arange(args.N + 1) % 2 == 1)
    odd_mask_check[0] = False
    odd_sigma_sample = sigma[odd_mask_check]
    odd_sigma_sample = odd_sigma_sample[odd_sigma_sample >= 0]
    # Use partition for top-K thresholds
    sigma_sorted_desc = np.sort(odd_sigma_sample)[::-1]
    print(f"  sort done in {time.perf_counter()-t0:.1f}s   "
          f"top-10 sigma values: {sigma_sorted_desc[:10].tolist()}", flush=True)

    print()
    print(f"  {'K':>8} {'threshold':>11} {'n_obs':>13} "
          f"{'beta':>11} {'gap':>9} {'excluded':>10}")

    for K in Ks:
        if K == 0:
            sx, sy, sxy, sx2, cnt = streaming_ols_range(sigma, 0, args.N + 1)
            ex = 0
            threshold = max_sigma + 1
        else:
            if K > len(sigma_sorted_desc):
                continue
            threshold = int(sigma_sorted_desc[K - 1])
            sx, sy, sxy, sx2, cnt, ex = streaming_ols_drop_top(sigma, args.N, threshold)
        a, b, n_obs = fit_beta(sx, sy, sxy, sx2, cnt)
        gap = heuristic - b
        print(f"  {K:>8} {threshold:>11} {n_obs:>13,} "
              f"{b:>11.6f} {gap:>+9.6f} {ex:>10}")


if __name__ == "__main__":
    main()
