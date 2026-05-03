"""
qx+1 generalization data generator.

For each odd integer n in [1, N], computes the trajectory under the map
  T_q(n) = n/2          if n even
  T_q(n) = q*n + 1      if n odd
and records, with a max-step and max-value cutoff:
  - sigma_q(n):       steps to reach 1 (NaN if cutoff hit)
  - converged(n):     True if reached 1, False if cutoff
  - max_excursion(n)
  - odd_steps(n), even_steps(n)
  - residues n mod 2^k for k in {3, 4, 5, 6}

For q = 3 this is the standard Collatz; for q >= 5 the heuristic factor is
q/4 > 1 and many orbits are believed to diverge. Cutoff parameters are
configurable.

Outputs to data/q_main_q{q}_N{N}.parquet.

Usage:
    python generate_q.py --q 3 --N 1000000
    python generate_q.py --q 5 --N 1000000 --max_steps 100000 --max_value 1e18
"""
import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
import numba
from numba import njit, prange


@njit(cache=True, boundscheck=False, parallel=True)
def compute_q_trajectories(N, q, max_steps, max_value):
    """For each odd n in [1, N], compute qx+1 trajectory stats."""
    sigma = np.full(N + 1, -1, dtype=np.int64)        # -1 = not computed
    converged = np.zeros(N + 1, dtype=np.bool_)
    max_exc = np.zeros(N + 1, dtype=np.int64)
    odd_steps = np.zeros(N + 1, dtype=np.int64)
    even_steps = np.zeros(N + 1, dtype=np.int64)

    for n in prange(1, N + 1):
        if n % 2 == 0:
            continue
        x = n
        steps = 0
        os = 0
        es = 0
        m = x
        while x != 1 and steps < max_steps and x < max_value:
            if x > m:
                m = x
            if x % 2 == 0:
                x = x // 2
                es += 1
            else:
                x = q * x + 1
                os += 1
            steps += 1
        if x == 1:
            converged[n] = True
            sigma[n] = steps
        else:
            converged[n] = False
            sigma[n] = -1  # sentinel
        max_exc[n] = m
        odd_steps[n] = os
        even_steps[n] = es
    return sigma, converged, max_exc, odd_steps, even_steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--N", type=int, default=1_000_000)
    ap.add_argument("--max_steps", type=int, default=100_000)
    ap.add_argument("--max_value", type=float, default=1e18)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = Path(args.out) if args.out else here / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[gen] q={args.q}, N={args.N:,}, max_steps={args.max_steps:,}, max_value={args.max_value:.2e}", flush=True)
    t0 = time.perf_counter()
    sigma, converged, max_exc, odd_s, even_s = compute_q_trajectories(
        args.N, args.q, args.max_steps, int(args.max_value)
    )
    dt = time.perf_counter() - t0
    print(f"[gen] done in {dt:.2f}s", flush=True)

    odd_mask = np.arange(args.N + 1) % 2 == 1
    odd_mask[0] = False
    n_arr = np.arange(args.N + 1)[odd_mask].astype(np.int64)
    sigma_arr = sigma[odd_mask].astype(np.int64)
    conv_arr = converged[odd_mask]
    max_exc_arr = max_exc[odd_mask].astype(np.int64)
    odd_s_arr = odd_s[odd_mask].astype(np.int64)
    even_s_arr = even_s[odd_mask].astype(np.int64)

    print(f"[gen] {len(n_arr):,} odd n's processed", flush=True)
    n_conv = conv_arr.sum()
    print(f"[gen] converged: {n_conv:,} ({100.0*n_conv/len(n_arr):.2f}%)", flush=True)
    if args.q == 3:
        # sanity check
        check_pairs = [(1, 0), (3, 7), (5, 5), (7, 16), (27, 111), (703, 170), (871, 178), (6171, 261)]
        for n_check, expected in check_pairs:
            if n_check <= args.N:
                idx = (n_arr == n_check).argmax()
                if n_arr[idx] == n_check:
                    actual = sigma_arr[idx]
                    status = "OK" if actual == expected else "MISMATCH"
                    print(f"  sanity sigma({n_check}) = {actual}  expected {expected}  [{status}]", flush=True)

    df = pl.DataFrame({
        "n": n_arr,
        "sigma_q": sigma_arr,            # -1 if did not converge
        "converged": conv_arr,
        "max_excursion": max_exc_arr,
        "odd_steps": odd_s_arr,
        "even_steps": even_s_arr,
        "res_mod_8": (n_arr % 8).astype(np.int16),
        "res_mod_16": (n_arr % 16).astype(np.int16),
        "res_mod_32": (n_arr % 32).astype(np.int16),
        "res_mod_64": (n_arr % 64).astype(np.int16),
    })
    out_path = out_dir / f"q_main_q{args.q}_N{args.N}.parquet"
    df.write_parquet(out_path)
    sz_mb = out_path.stat().st_size / (1 << 20)
    print(f"[save] {out_path}  ({sz_mb:.2f} MB)", flush=True)


if __name__ == "__main__":
    main()
