"""
Experiment 29 — qx+1 cycle classification for q in {5, 7, 11, 13}.

For each odd n in the existing parquet for that q, walk T_q under the
full map T_q(x) = (qx+1) if odd else x/2, with Floyd's tortoise-and-hare
cycle detection. Classify each orbit:
  status = 0  reached trivial cycle (smallest member = 1)
  status = 1  reached a non-trivial cycle (record smallest member, cycle length)
  status = 2  divergent (exceeded max_value)
  status = 3  timeout (exceeded max_steps without cycling or diverging)

Then tabulate by residue class r mod 2^k and by a★(r mod 2^k) where
a★ is the qx+1 prefix decomposition's terminal odd a-coefficient.

Outputs:
  experiments_output/29_qx1_cycles_q{q}_N{N}.parquet  per-orbit classification
  experiments_output/29_qx1_cycle_catalog_q{q}.csv     unique cycles found
  experiments_output/29_qx1_per_class_outcomes_q{q}_k{k}.csv  fraction tables

Usage:
  python 29_qx1_cycle_classification.py --q 5  --parquet data/q_main_q5_N100000000.parquet
  python 29_qx1_cycle_classification.py --q 7  --parquet data/q_main_q7_N100000000.parquet
  python 29_qx1_cycle_classification.py --q 11 --parquet data/q_main_q11_N1000000000.parquet
  python 29_qx1_cycle_classification.py --q 13 --parquet data/q_main_q13_N10000000.parquet
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

# Status codes
S_TRIVIAL  = 0   # reached 1 (full-map trivial cycle {1,6,3,16,8,4,2} for q=5, etc.)
S_NONTRIV  = 1   # reached a cycle with smallest != 1
S_DIVERGE  = 2   # exceeded max_value
S_TIMEOUT  = 3   # exceeded max_steps


@njit(cache=True)
def step_q(x, q):
    """Full-map T_q step."""
    if x % 2 == 0:
        return x // 2
    else:
        return q * x + 1


@njit(cache=True)
def classify_one(start, q, max_value, max_steps):
    """Floyd cycle detection. Returns (status, smallest, cycle_len, steps_taken).
    Uses int64 with pre-step bounds check to avoid overflow.
    """
    # Floyd's: tortoise advances 1, hare advances 2
    t = start
    h = start
    steps = 0
    # Advance both until t == h (cycle), or h diverges, or timeout.
    safe_cap = max_value // q  # to avoid q*x+1 overflow at extreme x
    while steps < max_steps:
        # Tortoise: 1 step
        if (t % 2 == 1) and (t > safe_cap):
            return S_DIVERGE, np.int64(0), np.int64(0), steps
        t = step_q(t, q)
        # Hare: 2 steps
        if (h % 2 == 1) and (h > safe_cap):
            return S_DIVERGE, np.int64(0), np.int64(0), steps
        h = step_q(h, q)
        if h > max_value:
            return S_DIVERGE, np.int64(0), np.int64(0), steps
        if (h % 2 == 1) and (h > safe_cap):
            return S_DIVERGE, np.int64(0), np.int64(0), steps
        h = step_q(h, q)
        if h > max_value:
            return S_DIVERGE, np.int64(0), np.int64(0), steps
        steps += 1
        if t == h:
            break
    if t != h:
        return S_TIMEOUT, np.int64(0), np.int64(0), steps

    # Found cycle: re-locate cycle entry by resetting tortoise
    t2 = start
    while t2 != h:
        t2 = step_q(t2, q)
        h  = step_q(h, q)
        if h > max_value:
            # shouldn't happen since we're inside a cycle, but safe-out
            return S_DIVERGE, np.int64(0), np.int64(0), steps
    cycle_entry = t2

    # Walk one full cycle, find smallest member and cycle length
    smallest = cycle_entry
    cur = step_q(cycle_entry, q)
    cycle_len = np.int64(1)
    while cur != cycle_entry:
        if cur < smallest:
            smallest = cur
        cur = step_q(cur, q)
        cycle_len += 1
        if cycle_len > max_steps:
            return S_TIMEOUT, np.int64(0), np.int64(0), steps

    if smallest == 1:
        return S_TRIVIAL, np.int64(1), cycle_len, steps
    else:
        return S_NONTRIV, smallest, cycle_len, steps


@njit(parallel=True, cache=True)
def classify_all(starts, q, max_value, max_steps):
    n = len(starts)
    status   = np.zeros(n, dtype=np.int8)
    smallest = np.zeros(n, dtype=np.int64)
    cyclen   = np.zeros(n, dtype=np.int64)
    nsteps   = np.zeros(n, dtype=np.int64)
    for i in prange(n):
        s, sm, cl, st = classify_one(starts[i], q, max_value, max_steps)
        status[i]   = s
        smallest[i] = sm
        cyclen[i]   = cl
        nsteps[i]   = st
    return status, smallest, cyclen, nsteps


def qx1_prefix(r, k, q):
    """qx+1 prefix decomposition: starting state (a=2^k, c=r), apply rules
    until a is odd. Returns (prefix_steps, a_star, c_star)."""
    a, c = 1 << k, r
    steps = 0
    while a % 2 == 0 and steps < 4 * k + 50:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a = q * a; c = q * c + 1
        steps += 1
    return steps, a, c


def build_per_class_table(df, q, k, total_count_per_class, ks_residues):
    """For each odd residue r mod 2^k, count outcomes."""
    out_rows = []
    for r in range(1, 1 << k, 2):
        sub = df.filter(pl.col(f"res_mod_{1<<k}") == r) if f"res_mod_{1<<k}" in df.columns else None
        if sub is None:
            # fall back to mask
            sub = df.filter((pl.col("n") % (1 << k)) == r)
        n_total = len(sub)
        if n_total == 0:
            continue
        n_trivial  = sub.filter(pl.col("status") == S_TRIVIAL).height
        n_nontriv  = sub.filter(pl.col("status") == S_NONTRIV).height
        n_diverge  = sub.filter(pl.col("status") == S_DIVERGE).height
        n_timeout  = sub.filter(pl.col("status") == S_TIMEOUT).height
        # Per-cycle breakdown for non-trivial
        cycle_breakdown = {}
        if n_nontriv > 0:
            sub_nt = sub.filter(pl.col("status") == S_NONTRIV)
            for sm in sub_nt["smallest"].unique().to_list():
                cycle_breakdown[int(sm)] = sub_nt.filter(pl.col("smallest") == sm).height
        # a★ for this class
        prefix_steps, a_star, c_star = qx1_prefix(r, k, q)
        out_rows.append(dict(
            r=r,
            n=n_total,
            p_trivial=n_trivial / n_total,
            p_nontriv=n_nontriv / n_total,
            p_diverge=n_diverge / n_total,
            p_timeout=n_timeout / n_total,
            a_star=a_star,
            c_star=c_star,
            prefix_steps=prefix_steps,
            cycle_breakdown=str(cycle_breakdown) if cycle_breakdown else "",
        ))
    return pl.DataFrame(out_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--parquet", type=str, required=True)
    ap.add_argument("--max_value", type=float, default=1e18)
    ap.add_argument("--max_steps", type=int, default=1_000_000)
    ap.add_argument("--max_n", type=int, default=None,
                    help="Optional cap on n (use only n <= max_n)")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.parquet}", flush=True)
    df = pl.read_parquet(args.parquet)
    if args.max_n is not None:
        df = df.filter(pl.col("n") <= args.max_n)
    starts = df["n"].to_numpy().astype(np.int64)
    N_orig = int(starts.max())
    print(f"        {len(starts):,} orbits  q={args.q}  max_value={args.max_value:.1e}  N_max={N_orig:,}", flush=True)

    # Warm up jit
    _ = classify_all(np.array([1, 3, 5], dtype=np.int64), args.q, int(args.max_value), 100)

    print(f"[run] Floyd cycle detection ...", flush=True)
    t0 = time.perf_counter()
    status, smallest, cyclen, nsteps = classify_all(
        starts, args.q, int(args.max_value), args.max_steps)
    t = time.perf_counter() - t0
    print(f"        done in {t:.1f}s", flush=True)

    # Status histogram
    n_total = len(starts)
    n_trivial = int((status == S_TRIVIAL).sum())
    n_nontriv = int((status == S_NONTRIV).sum())
    n_diverge = int((status == S_DIVERGE).sum())
    n_timeout = int((status == S_TIMEOUT).sum())
    print(f"\n  status counts:")
    print(f"    trivial cycle: {n_trivial:>12,}  ({100*n_trivial/n_total:.4f}%)")
    print(f"    other cycle:   {n_nontriv:>12,}  ({100*n_nontriv/n_total:.4f}%)")
    print(f"    divergent:     {n_diverge:>12,}  ({100*n_diverge/n_total:.4f}%)")
    print(f"    timeout:       {n_timeout:>12,}  ({100*n_timeout/n_total:.4f}%)")

    # Cycle catalog
    if n_nontriv > 0:
        cat_rows = []
        sm_arr = smallest[status == S_NONTRIV]
        cl_arr = cyclen[status == S_NONTRIV]
        unique_smallest = np.unique(sm_arr)
        print(f"\n  cycle catalog (non-trivial):")
        print(f"    {'smallest':>12} {'cycle_len':>10} {'count':>10}  example walk (first 10 members from smallest)")
        for sm in sorted(unique_smallest.tolist()):
            mask = sm_arr == sm
            count = int(mask.sum())
            cl = int(cl_arr[mask][0])
            # Trace cycle from smallest
            walk = [int(sm)]
            x = int(sm)
            for _ in range(min(cl - 1, 9)):
                if x % 2 == 0:
                    x = x // 2
                else:
                    x = args.q * x + 1
                walk.append(int(x))
            print(f"    {sm:>12} {cl:>10} {count:>10}  {walk}")
            cat_rows.append(dict(q=args.q, smallest=int(sm), cycle_length=cl, count=count, walk=str(walk)))
        cat_df = pl.DataFrame(cat_rows)
        cat_df.write_csv(out_dir / f"29_qx1_cycle_catalog_q{args.q}.csv")

    # Save per-orbit results
    classified = df.with_columns([
        pl.Series("status", status),
        pl.Series("smallest", smallest),
        pl.Series("cycle_length", cyclen),
        pl.Series("n_steps_floyd", nsteps),
    ])
    out_pq = out_dir / f"29_qx1_cycles_q{args.q}_N{N_orig}.parquet"
    classified.write_parquet(out_pq)
    print(f"\n[save] {out_pq}")

    # Per-class outcome tables
    for k in [6, 8, 10] if args.q == 5 else ([6, 8] if args.q == 7 else [6]):
        if (1 << k) > 65536:
            continue
        # Build residue column if not present
        residue_col = f"res_mod_{1<<k}"
        if residue_col not in classified.columns:
            classified_k = classified.with_columns(
                (pl.col("n") % (1 << k)).cast(pl.Int32).alias(residue_col)
            )
        else:
            classified_k = classified
        per_class = build_per_class_table(classified_k, args.q, k, None, None)
        if len(per_class) == 0:
            continue
        out_pc = out_dir / f"29_qx1_per_class_outcomes_q{args.q}_k{k}.csv"
        per_class.write_csv(out_pc)
        print(f"[save] {out_pc}  ({len(per_class)} classes)")
        # a★ pooled summary
        a_star_summary = per_class.group_by("a_star").agg([
            pl.col("n").sum().alias("n_total"),
            (pl.col("p_trivial") * pl.col("n")).sum().alias("trivial_count"),
            (pl.col("p_nontriv") * pl.col("n")).sum().alias("nontriv_count"),
            (pl.col("p_diverge") * pl.col("n")).sum().alias("diverge_count"),
            pl.col("r").count().alias("n_classes"),
        ]).with_columns([
            (pl.col("trivial_count") / pl.col("n_total")).alias("p_trivial"),
            (pl.col("nontriv_count") / pl.col("n_total")).alias("p_nontriv"),
            (pl.col("diverge_count") / pl.col("n_total")).alias("p_diverge"),
        ]).sort("a_star")
        print(f"\n  k={k} a★ pooled outcome (a★ collapse test):")
        print(a_star_summary.select(["a_star", "n_classes", "n_total", "p_trivial", "p_nontriv", "p_diverge"]))


if __name__ == "__main__":
    main()
