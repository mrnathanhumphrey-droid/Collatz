"""
Experiment 28 — per-octave trajectory E[v], closed-form prediction of beta_local.

The standard random-walk heuristic for sigma:
    sigma(m) ~ (1 + E[v]) * # Syracuse steps
    # Syracuse steps ~ log(m) / (log(2) * E[v] - log(3))    (for E[v] > log(3)/log(2))
    => sigma(m) ~ K(E[v]) * log(m)
    where  K(u) = (1 + u) / (u * log(2) - log(3))

For Geom(1/2)  u = E[v] = 2  =>  K = 3 / log(4/3) = 10.4282.

But trajectory measure on v has structural deviations from Geom(1/2)
(see exp 25). If E[v] under the trajectory measure shifts with m's scale,
the local slope K_local shifts accordingly. dK/du ~ -21.7 at u=2, so a
0.01 shift in E[v] moves K by ~0.22.

This experiment measures trajectory E[v] separately within each octave
of starting values [2^j, 2^(j+1)] and predicts K_j. Compares to beta_local
from exp 27.

If K_predicted matches beta_local octave-by-octave, the per-octave deviation
from heuristic is fully explained by one trajectory-measure parameter.

Usage:
    python 28_per_octave_trajectory_E_v.py --N_per_octave 1000000 --T 500
"""
import argparse
import sys
import time

import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")


@njit(parallel=True, cache=True)
def collect_v_per_octave(starts, T, max_v):
    """Chunk-parallel: each chunk owns its own histogram row.
    starts: array of odd integer starting values.
    Returns hist: histogram of v values, total step count, count of trajectories
    that reached 1 within T steps.
    """
    n = len(starts)
    chunk = max(1, n // 256)
    n_chunks = (n + chunk - 1) // chunk
    hist_local = np.zeros((n_chunks, max_v + 2), dtype=np.int64)
    cnt_steps_local = np.zeros(n_chunks, dtype=np.int64)
    cnt_reached_local = np.zeros(n_chunks, dtype=np.int64)
    for c in prange(n_chunks):
        start_i = c * chunk
        end_i = min(start_i + chunk, n)
        for i in range(start_i, end_i):
            m = starts[i]
            for _ in range(T):
                if m == 1:
                    cnt_reached_local[c] += 1
                    break
                tmp = 3 * m + 1
                v = 0
                while tmp % 2 == 0:
                    tmp //= 2
                    v += 1
                if v <= max_v:
                    hist_local[c, v] += 1
                else:
                    hist_local[c, max_v + 1] += 1
                cnt_steps_local[c] += 1
                m = tmp
    return hist_local.sum(axis=0), cnt_steps_local.sum(), cnt_reached_local.sum()


def K_pred(E_v):
    """Heuristic K(E[v]) = (1 + E[v]) / (log(2)*E[v] - log(3))."""
    return (1.0 + E_v) / (E_v * np.log(2.0) - np.log(3.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N_per_octave", type=int, default=1_000_000)
    ap.add_argument("--T", type=int, default=500)
    ap.add_argument("--max_v", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--octaves", type=str, default="17,18,19,20,21,22,23,24,25,26")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    octaves = [int(j) for j in args.octaves.split(",")]

    # beta_local from exp 27 (run on N=2^27 cache).
    # Reproducing here so we have a single side-by-side report.
    beta_local_exp27 = {
        17: 10.6561, 18: 10.7531, 19: 10.6955, 20: 10.7198,
        21: 10.8765, 22: 10.8873, 23: 10.7797, 24: 10.6539,
        25: 10.5863, 26: 10.4888,
    }

    geom_E_v = 2.0
    geom_K = K_pred(geom_E_v)

    print(f"[setup] N_per_octave = {args.N_per_octave:,}  T = {args.T}",
          flush=True)
    print(f"[ref]   Geom(1/2): E[v] = {geom_E_v}, K = {geom_K:.6f}", flush=True)
    print()
    print(f"  {'octave':>6} {'lo':>14} {'hi':>14} {'reached_1':>10} "
          f"{'n_v':>14} {'E[v]':>9} {'K_pred':>10} "
          f"{'beta_local':>11} {'pred-meas':>11}")

    summary = []
    for j in octaves:
        lo = 1 << j
        hi = 1 << (j + 1)
        # sample N_per_octave odd integers uniformly from [lo, hi)
        # odd integers in [lo, hi): (hi - lo) / 2 of them, starting at lo+1 if lo even
        # since lo = 2^j is even, odd starts are lo+1, lo+3, ..., hi-1
        n_odd_in_octave = (hi - lo) // 2
        sample_size = min(args.N_per_octave, n_odd_in_octave)
        # Pick random odd m's
        offsets = rng.integers(0, n_odd_in_octave, size=sample_size, dtype=np.int64)
        starts = (lo + 1 + 2 * offsets).astype(np.int64)

        t0 = time.perf_counter()
        hist, n_steps, n_reached = collect_v_per_octave(starts, args.T, args.max_v)
        t1 = time.perf_counter()

        # E[v] from histogram
        n_v_total = int(hist[1:args.max_v + 1].sum())
        if n_v_total == 0:
            print(f"  {j:>6} {lo:>14,} {hi - 1:>14,} {n_reached:>10,} "
                  f"{n_v_total:>14,} (no v samples)")
            continue
        E_v = sum(v * hist[v] for v in range(1, args.max_v + 1)) / n_v_total
        K = K_pred(E_v)

        beta_local = beta_local_exp27.get(j, float("nan"))
        diff = K - beta_local

        print(f"  {j:>6} {lo:>14,} {hi - 1:>14,} {n_reached:>10,} "
              f"{n_v_total:>14,} {E_v:>9.5f} {K:>10.4f} "
              f"{beta_local:>11.4f} {diff:>+11.4f}    [{t1-t0:.1f}s]")

        summary.append((j, E_v, K, beta_local, diff))

    print()
    print("Interpretation:")
    print(f"  K_pred(E[v]) = (1 + E[v]) / (E[v]*log(2) - log(3))")
    print(f"  If trajectory measure E[v] in octave j matches the empirical")
    print(f"  beta_local in that octave (via the K formula), per-octave beta")
    print(f"  is fully explained by per-octave trajectory mean v.")
    print()
    print("  Magnitude of expected effect: dK/dE[v] = -log(6)/(E[v]*log(2)-log(3))^2")
    print(f"                              = {-np.log(6) / (np.log(2)*2 - np.log(3))**2:.2f}")
    print(f"  i.e., 0.01 shift in E[v] moves K by ~0.22.")

    # Save
    import polars as pl
    if summary:
        df = pl.DataFrame(
            summary, schema=["octave", "E_v", "K_pred", "beta_local_exp27", "K_pred_minus_beta"],
            orient="row",
        )
        from pathlib import Path
        out = Path(__file__).resolve().parent.parent / "experiments_output"
        out.mkdir(parents=True, exist_ok=True)
        df.write_csv(out / f"28_per_octave_E_v_Nperoct{args.N_per_octave}.csv")
        print(f"[save] {out / f'28_per_octave_E_v_Nperoct{args.N_per_octave}.csv'}")


if __name__ == "__main__":
    main()
