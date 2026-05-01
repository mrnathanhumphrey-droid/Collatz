"""
Stage 1: Data generation for hierarchical Bayesian Collatz analysis.

For every integer n in [1, N], computes:
  - sigma(n):         total stopping time under standard Collatz
  - syracuse(n):      odd-step count = stopping time under the compressed map
                      n -> (3n+1) / 2^v  for odd n, with v = nu_2(3n+1)
  - max_excursion(n): largest value reached along the orbit from n
  - odd_steps(n) / even_steps(n)
  - residue n mod 2^k for each k in --k (default 4, 6, 8 i.e. mod 16/64/256)
  - is_record(n):     True iff sigma(n) > sigma(n') for all n' < n

Outputs two parquet files:
  - data/main_N{N}.parquet   one row per n, scalar columns
  - data/v_seq_N{N}.parquet  long format (n, step_idx, v) for every Syracuse
                             step in every trajectory (sum of syracuse[n] rows)

Algorithm. Cache pass: walk each new n's trajectory until we hit either 1 or
a value <= N already in the cache, eagerly caching every value <= N visited
along the way. This keeps the per-n marginal work to about path_len ~ O(1)
once the cache fills up. JIT-compiled with numba; cache pass is sequential
(dependency on smaller n), v-seq pass is parallelized with prange.
"""

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
import numba
from numba import njit, prange


KNOWN_SIGMA = {
    1: 0,    2: 1,    3: 7,    4: 2,    5: 5,    6: 8,    7: 16,   8: 3,
    9: 19,   10: 6,   27: 111, 703: 170, 871: 178, 6171: 261, 77031: 350,
    837799: 524,
}

PATH_MAX = 100_000


@njit(cache=True, boundscheck=False)
def compute_caches(N):
    """
    Sequentially fill sigma / syracuse / max_excursion for n in 1..N.

    For each uncached n, walk the standard Collatz trajectory until we hit
    a cached value m <= N (or m == 1). Then backfill the three accumulators
    in reverse trajectory order, eagerly caching every value <= N seen.
    """
    sigma    = np.full(N + 1, -1, dtype=np.int32)
    syracuse = np.full(N + 1, -1, dtype=np.int32)
    max_exc  = np.full(N + 1, -1, dtype=np.int64)

    sigma[1]    = 0
    syracuse[1] = 0
    max_exc[1]  = 1

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
                raise RuntimeError("path buffer overflow; raise PATH_MAX")
            path[path_len] = m
            path_len += 1
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            if m == 1:
                break

        if m == 1:
            cur_sigma = np.int32(0)
            cur_syr   = np.int32(0)
            cur_exc   = np.int64(1)
        else:
            cur_sigma = sigma[m]
            cur_syr   = syracuse[m]
            cur_exc   = max_exc[m]

        for i in range(path_len - 1, -1, -1):
            v = path[i]
            cur_sigma += 1
            if v & 1:
                cur_syr += 1
            if v > cur_exc:
                cur_exc = v
            if v <= N:
                sigma[v]    = cur_sigma
                syracuse[v] = cur_syr
                max_exc[v]  = cur_exc

    return sigma, syracuse, max_exc


@njit(cache=True, parallel=True, boundscheck=False)
def collect_v_sequences(N, syracuse, offsets, n_arr, step_arr, v_arr):
    """
    For each n in 1..N, walk its Syracuse trajectory and write v-values
    (= nu_2(3m+1) for each odd m visited) into pre-sized flat arrays at
    pre-computed offsets[n]. Embarrassingly parallel across n.
    """
    for n in prange(1, N + 1):
        m = np.int64(n)
        while m > 1 and (m & 1) == 0:
            m = m >> 1

        idx = offsets[n]
        step = np.int32(0)
        while m > 1:
            tmp = 3 * m + 1
            v = np.int8(0)
            while (tmp & 1) == 0:
                v += 1
                tmp = tmp >> 1
            n_arr[idx]    = n
            step_arr[idx] = step
            v_arr[idx]    = v
            idx += 1
            m = tmp
            step += 1


def compute_residues(N, k_list):
    n = np.arange(1, N + 1, dtype=np.int64)
    return {k: (n & ((1 << k) - 1)).astype(np.int16) for k in k_list}


@njit(cache=True, boundscheck=False)
def compute_records(sigma):
    """is_record[n] = True iff sigma[n] > sigma[n'] for every 1 <= n' < n."""
    L = sigma.shape[0]
    records = np.zeros(L, dtype=np.bool_)
    running_max = np.int32(-1)
    for n in range(1, L):
        if sigma[n] > running_max:
            records[n] = True
            running_max = sigma[n]
    return records


def sanity_check(sigma):
    failures = []
    for n, expected in KNOWN_SIGMA.items():
        if n < len(sigma):
            got = int(sigma[n])
            if got != expected:
                failures.append((n, expected, got))
    if failures:
        msg = "; ".join(f"sigma({n})={g} expected {e}" for n, e, g in failures)
        raise AssertionError(f"Sanity check failed: {msg}")
    largest = max(k for k in KNOWN_SIGMA if k < len(sigma))
    print(f"[sanity]  sigma matches all known checkpoints up to n={largest:,}",
          flush=True)


def main():
    ap = argparse.ArgumentParser(description="Stage 1 Collatz data generation.")
    ap.add_argument("--N", type=int, default=1 << 20,
                    help="Compute for n in [1, N]. Default 2**20 = 1,048,576.")
    ap.add_argument("--k", type=int, nargs="+", default=[4, 6, 8],
                    help="Residue powers: emit n mod 2^k for each k. Default 4 6 8.")
    ap.add_argument("--out", type=str, default=None,
                    help="Output dir. Default: <script_dir>/data")
    ap.add_argument("--no-vseq", action="store_true",
                    help="Skip v-sequence parquet (useful for quick smoke runs).")
    ap.add_argument("--threads", type=int, default=16,
                    help="Numba thread cap. Default 16 (leaves headroom on a "
                         "9950X3D for other workloads).")
    args = ap.parse_args()

    numba.set_num_threads(args.threads)
    print(f"[threads] numba using {numba.get_num_threads()} threads "
          f"(of {numba.config.NUMBA_DEFAULT_NUM_THREADS} available)", flush=True)

    N = args.N
    out_dir = Path(args.out) if args.out else Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    main_path = out_dir / f"main_N{N}.parquet"
    vseq_path = out_dir / f"v_seq_N{N}.parquet"

    print(f"[generate] N = {N:,} ({hex(N)})", flush=True)
    print(f"[generate] out = {out_dir}", flush=True)

    t0 = time.perf_counter()
    sigma, syracuse, max_exc = compute_caches(N)
    t1 = time.perf_counter()
    print(f"[caches]  built in {t1 - t0:.2f}s "
          f"(max sigma = {int(sigma[1:].max())}, "
          f"max excursion = {int(max_exc[1:].max()):,})", flush=True)

    sanity_check(sigma)

    residues = compute_residues(N, args.k)
    records  = compute_records(sigma)

    n_index    = np.arange(1, N + 1, dtype=np.int64)
    odd_steps  = syracuse[1:].astype(np.int32)
    even_steps = (sigma[1:] - syracuse[1:]).astype(np.int32)

    main_cols = {
        "n":              n_index,
        "sigma":          sigma[1:].astype(np.int32),
        "syracuse":       syracuse[1:].astype(np.int32),
        "odd_steps":      odd_steps,
        "even_steps":     even_steps,
        "max_excursion":  max_exc[1:].astype(np.int64),
        "is_record":      records[1:],
    }
    for k in args.k:
        main_cols[f"res_mod_{1 << k}"] = residues[k]

    df_main = pl.DataFrame(main_cols)
    df_main.write_parquet(main_path, compression="zstd", compression_level=3)
    print(f"[main]    {df_main.shape[0]:,} rows -> {main_path.name} "
          f"({main_path.stat().st_size / 1e6:.1f} MB)", flush=True)

    if args.no_vseq:
        print("[v-seq]   skipped (--no-vseq)", flush=True)
    else:
        cs = np.cumsum(syracuse[1:N + 1].astype(np.int64))
        offsets = np.empty(N + 2, dtype=np.int64)
        offsets[0] = 0
        offsets[1] = 0
        offsets[2:N + 2] = cs
        total = int(cs[-1])

        n_arr    = np.empty(total, dtype=np.int64)
        step_arr = np.empty(total, dtype=np.int32)
        v_arr    = np.empty(total, dtype=np.int8)

        t2 = time.perf_counter()
        collect_v_sequences(N, syracuse, offsets, n_arr, step_arr, v_arr)
        t3 = time.perf_counter()
        print(f"[v-seqs]  {total:,} rows generated in {t3 - t2:.2f}s "
              f"(parallel)", flush=True)

        df_vseq = pl.DataFrame({"n": n_arr, "step_idx": step_arr, "v": v_arr})
        df_vseq.write_parquet(vseq_path, compression="zstd", compression_level=3)
        print(f"[v-seq]   -> {vseq_path.name} "
              f"({vseq_path.stat().st_size / 1e6:.1f} MB)", flush=True)

    print(f"[done]    total {time.perf_counter() - t0:.2f}s", flush=True)


if __name__ == "__main__":
    main()
