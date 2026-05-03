"""
Experiment 25 — high-resolution trajectory measure of v = nu_2(3m+1).

The natural-density measure on odd m gives v ~ Geom(1/2): P(v=k) = 2^{-k}.
Along trajectories of the Syracuse map T(m) = (3m+1)/2^{v(m)}, however, v has
a deviation: previous low-resolution measurement (N=2^23) showed v=4, v=10
spiked at ~1.32x and v=6,7,11,12,14 sagged. The question is whether this
deviation has 2-adic / 3-adic structure that points to a closed-form
trajectory measure.

This script samples N_start uniform odd starts, walks T Syracuse steps from
each, and pools all observed v's. With N_start = 10^6, T = 200 we get
~2 x 10^8 v-samples — about 10000x the original measurement. Per-bin SE
on the ratio: ~1/sqrt(predicted_count). For v=10 we expect ~2x10^5
counts, SE ~0.2%; for v=20, ~200 counts, SE ~7%.

Then look for structure in the ratio: residue patterns of 2^v mod 3^k for
small k, etc.

Usage:
    python 25_trajectory_measure.py --N_start 1000000 --T 200 --max_v 35
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")


@njit(parallel=True, cache=True)
def collect_v_pooled(starts, T, max_v):
    """Chunk-parallel: each prange iteration owns its own histogram row.
    No races possible because chunk c writes only to hist_local[c, :].
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


@njit(parallel=True, cache=True)
def collect_v_first_step(starts, max_v):
    """First-step v only: v(m) for each start. Sanity check — should be Geom(1/2)."""
    n = len(starts)
    chunk = max(1, n // 256)
    n_chunks = (n + chunk - 1) // chunk
    hist_local = np.zeros((n_chunks, max_v + 2), dtype=np.int64)
    for c in prange(n_chunks):
        start_i = c * chunk
        end_i = min(start_i + chunk, n)
        for i in range(start_i, end_i):
            tmp = 3 * starts[i] + 1
            v = 0
            while tmp % 2 == 0:
                tmp //= 2
                v += 1
            if v <= max_v:
                hist_local[c, v] += 1
            else:
                hist_local[c, max_v + 1] += 1
    return hist_local.sum(axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N_start", type=int, default=1_000_000,
                    help="Number of random odd starting integers.")
    ap.add_argument("--T", type=int, default=200,
                    help="Syracuse steps per trajectory.")
    ap.add_argument("--N_max", type=int, default=10**9,
                    help="Upper bound for odd starts.")
    ap.add_argument("--max_v", type=int, default=35)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    # Generate odd starts in [3, N_max]
    starts = 2 * rng.integers(1, args.N_max // 2, size=args.N_start, dtype=np.int64) + 1

    print(f"[setup] N_start={args.N_start:,}  T={args.T}  N_max={args.N_max:,}",
          flush=True)

    # First-step (sanity check)
    print("[run] first-step v (should match Geom(1/2)) ...", flush=True)
    t0 = time.perf_counter()
    hist_first = collect_v_first_step(starts, args.max_v)
    t_first = time.perf_counter() - t0
    print(f"       done in {t_first:.1f}s", flush=True)

    # Trajectory v
    print("[run] trajectory v (T Syracuse steps per start) ...", flush=True)
    t0 = time.perf_counter()
    hist_traj, n_steps, n_reached_1 = collect_v_pooled(starts, args.T, args.max_v)
    t_traj = time.perf_counter() - t0
    print(f"       done in {t_traj:.1f}s", flush=True)
    print(f"       total v-samples: {hist_traj.sum():,}",
          f"trajectories that reached 1: {n_reached_1:,} / {args.N_start:,}",
          flush=True)
    print()

    # Compute ratios
    total_first = hist_first[:args.max_v + 1].sum()
    total_traj = hist_traj[:args.max_v + 1].sum()
    # Geom(1/2) on positive integers: P(v=k) = 2^(-k) for k = 1, 2, 3, ...
    # No normalization — the geometric sum from v=1 to infinity already equals 1.
    pred_geom = np.zeros(args.max_v + 1)
    for v in range(1, args.max_v + 1):
        pred_geom[v] = 2.0 ** (-v)
    pred_geom_norm = pred_geom

    expected_first = pred_geom_norm * total_first
    expected_traj = pred_geom_norm * total_traj

    print(f"{'v':>3} {'2^-v':>10} {'first_emp':>10} {'first_ratio':>11} "
          f"{'traj_emp':>14} {'traj_ratio':>11} {'SE_traj':>9} "
          f"{'2^v_mod9':>9} {'2^v_mod27':>10}")

    rows_out = []
    for v in range(1, args.max_v + 1):
        first_count = hist_first[v]
        traj_count = hist_traj[v]
        first_emp = first_count / total_first
        traj_emp = traj_count / total_traj
        first_ratio = first_emp / pred_geom_norm[v] if pred_geom_norm[v] > 0 else float("nan")
        traj_ratio = traj_emp / pred_geom_norm[v] if pred_geom_norm[v] > 0 else float("nan")
        # SE for ratio: ~ 1/sqrt(expected_count) under Poisson
        se_traj = 1.0 / np.sqrt(expected_traj[v]) if expected_traj[v] > 0 else float("nan")
        # 2^v mod small powers of 3 — for pattern hunting
        mod9 = pow(2, v, 9)
        mod27 = pow(2, v, 27)
        print(f"{v:>3} {pred_geom_norm[v]:>10.6f} {first_count:>10,} "
              f"{first_ratio:>11.4f} {traj_count:>14,} {traj_ratio:>11.4f} "
              f"{se_traj:>9.4f} {mod9:>9} {mod27:>10}")
        rows_out.append({
            "v": v,
            "geom_pred": pred_geom_norm[v],
            "first_count": int(first_count),
            "first_ratio": first_ratio,
            "traj_count": int(traj_count),
            "traj_ratio": traj_ratio,
            "se_traj": se_traj,
            "two_pow_v_mod9": mod9,
            "two_pow_v_mod27": mod27,
        })

    # Pattern-hunting: group ratios by 2^v mod 9, mod 27
    print()
    print("=== Pattern: trajectory ratio grouped by 2^v mod 9 ===")
    rows_arr = rows_out
    print(f"{'2^v mod 9':>10} {'count_v':>9} {'mean_ratio':>11} {'SD_ratio':>10}")
    by_mod9 = {}
    for r in rows_arr[:25]:  # restrict to v=1..25 where SE is small
        by_mod9.setdefault(r["two_pow_v_mod9"], []).append(r["traj_ratio"])
    for k, vals in sorted(by_mod9.items()):
        print(f"{k:>10} {len(vals):>9} {np.mean(vals):>11.4f} {np.std(vals):>10.4f}")

    print()
    print("=== Pattern: trajectory ratio grouped by 2^v mod 27 ===")
    print(f"{'2^v mod 27':>11} {'count_v':>9} {'mean_ratio':>11} {'SD_ratio':>10}")
    by_mod27 = {}
    for r in rows_arr[:25]:
        by_mod27.setdefault(r["two_pow_v_mod27"], []).append(r["traj_ratio"])
    for k, vals in sorted(by_mod27.items()):
        print(f"{k:>11} {len(vals):>9} {np.mean(vals):>11.4f} {np.std(vals):>10.4f}")

    # Save CSV
    import polars as pl
    out = pl.DataFrame(rows_out)
    out_csv = out_dir / f"25_trajectory_measure_Nstart{args.N_start}_T{args.T}.csv"
    out.write_csv(out_csv)
    print()
    print(f"[save] {out_csv}")


if __name__ == "__main__":
    main()
