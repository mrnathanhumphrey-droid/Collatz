"""
Experiment 26 — mu_beta convergence rate at N = 2^30 and 2^32.

Existing N-scaling table covers N in {2^20, 2^22, 2^23, 2^24, 2^25} with
gap from heuristic 10.4282 closing from 0.056 to 0.009. Two more data points
at N = 2^30, 2^32 distinguish O(N^{-1/2}) from logarithmic alternatives.

Approach: build sigma cache for n in [1, N] via memoized Collatz iteration
(numba JIT, sequential — depends on smaller-n cache hits), then a streaming
OLS pass on odd n only:

    sum_x  = sum  ln(n)            for odd n in [3, N]
    sum_y  = sum  sigma(n)
    sum_xy = sum  sigma(n) * ln(n)
    sum_x2 = sum  ln(n)^2
    n_obs  = count of odd n

    beta_OLS = (n_obs * sum_xy - sum_x * sum_y) /
               (n_obs * sum_x2 - sum_x * sum_x)

No parquet writing — just the slope. Memory: only the int32 sigma cache.
For 2^30: ~4 GB. For 2^32: ~17 GB. Both fit on 64 GB RAM.

Usage:
    python 26_mu_beta_n_extension.py --N 1073741824
    python 26_mu_beta_n_extension.py --N 4294967296
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
    """Sequentially fill sigma[n] for n in 1..N via memoized trajectory walk."""
    sigma = np.full(N + 1, -1, dtype=np.int32)
    sigma[1] = 0
    path = np.empty(PATH_MAX, dtype=np.int64)

    for n in range(2, N + 1):
        if sigma[n] != -1:
            continue
        m = np.int64(n)
        path_len = 0
        while True:
            if m <= N and sigma[m] != -1:
                break
            if path_len >= PATH_MAX:
                # Path overflow; fall back to running until m=1
                path_len = -1
                break
            path[path_len] = m
            path_len += 1
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            if m == 1:
                break

        if path_len < 0:
            # Overflow path: walk again without caching, just to end
            m = np.int64(n)
            steps_to_end = 0
            while m > 1:
                if m & 1:
                    m = 3 * m + 1
                else:
                    m = m >> 1
                steps_to_end += 1
            sigma[n] = np.int32(steps_to_end)
            continue

        if m == 1:
            cur_sigma = np.int32(0)
        else:
            cur_sigma = sigma[m]

        for i in range(path_len - 1, -1, -1):
            v = path[i]
            cur_sigma += 1
            if v <= N:
                sigma[v] = cur_sigma

    return sigma


@njit(cache=True, parallel=True, boundscheck=False)
def streaming_ols_odd(sigma, N):
    """Parallel streaming OLS sums over odd n in [3, N]."""
    n_threads = 64
    chunk = max(2, (N // n_threads) | 1)  # ensure chunk parity doesn't matter
    n_chunks = (N + chunk - 1) // chunk

    sx = np.zeros(n_chunks, dtype=np.float64)
    sy = np.zeros(n_chunks, dtype=np.float64)
    sxy = np.zeros(n_chunks, dtype=np.float64)
    sx2 = np.zeros(n_chunks, dtype=np.float64)
    cnt = np.zeros(n_chunks, dtype=np.int64)

    for c in prange(n_chunks):
        lo = c * chunk
        hi = min(lo + chunk, N + 1)
        if lo < 3:
            lo = 3
        if lo % 2 == 0:
            lo += 1
        s = sigma  # local
        for n in range(lo, hi, 2):
            x = np.log(np.float64(n))
            y = np.float64(s[n])
            sx[c] += x
            sy[c] += y
            sxy[c] += x * y
            sx2[c] += x * x
            cnt[c] += 1

    return sx.sum(), sy.sum(), sxy.sum(), sx2.sum(), cnt.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True, help="Upper bound (inclusive).")
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    numba.set_num_threads(args.threads)

    print(f"[setup] N = {args.N:,}  threads = {numba.get_num_threads()}",
          flush=True)
    print(f"[mem]   sigma cache: {(args.N + 1) * 4 / 1e9:.2f} GB", flush=True)
    print(f"[run]   building sigma cache ...", flush=True)

    t0 = time.perf_counter()
    sigma = compute_sigma_cache(args.N)
    t1 = time.perf_counter()
    print(f"[run]   cache built in {t1 - t0:.1f}s  "
          f"max(sigma) = {int(sigma[1:].max())}", flush=True)

    # Sanity check
    known = {1: 0, 27: 111, 703: 170, 871: 178, 6171: 261, 77031: 350, 837799: 524}
    for nk, sk in known.items():
        if nk <= args.N:
            got = int(sigma[nk])
            assert got == sk, f"sigma({nk}) = {got}, expected {sk}"
    print(f"[ok]    sigma matches known checkpoints", flush=True)

    print(f"[run]   streaming OLS on odd n in [3, {args.N:,}] ...", flush=True)
    t0 = time.perf_counter()
    sx, sy, sxy, sx2, cnt = streaming_ols_odd(sigma, args.N)
    t1 = time.perf_counter()
    print(f"[run]   OLS sums in {t1 - t0:.1f}s  (n_obs = {cnt:,})", flush=True)

    n_obs = float(cnt)
    denom = n_obs * sx2 - sx * sx
    beta = (n_obs * sxy - sx * sy) / denom
    alpha = (sy - beta * sx) / n_obs
    heuristic = 3.0 / (np.log(4.0) - np.log(3.0))

    print()
    print(f"=== Result at N = {args.N:,} = 2^{int(round(np.log2(args.N)))} ===")
    print(f"  n_obs (odd integers): {int(cnt):,}")
    print(f"  alpha (intercept)   : {alpha:.6f}")
    print(f"  beta  (slope)       : {beta:.6f}")
    print(f"  heuristic (3/ln(4/3)): {heuristic:.6f}")
    print(f"  gap = heuristic - beta : {heuristic - beta:>+.6f}")

    # Append to reference table
    print()
    print("Reference (from existing writeup):")
    print(f"  {'N':>14} {'log2(N)':>9} {'beta':>10} {'gap':>10}")
    ref = [
        (2 ** 20, 10.3723, heuristic - 10.3723),
        (2 ** 22, 10.3816, heuristic - 10.3816),
        (2 ** 23, 10.3845, heuristic - 10.3845),
        (2 ** 24, 10.4044, heuristic - 10.4044),
        (2 ** 25, 10.4191, heuristic - 10.4191),
    ]
    for Nv, b, g in ref:
        print(f"  {Nv:>14,} {int(round(np.log2(Nv))):>9} {b:>10.4f} {g:>10.4f}")
    print(f"  {args.N:>14,} {int(round(np.log2(args.N))):>9} {beta:>10.4f} "
          f"{heuristic - beta:>10.4f}    [this run]")


if __name__ == "__main__":
    main()
