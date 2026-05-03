"""
Experiment 22 — Cycle detection on q=5 non-convergent orbits.

For 5x+1 there's the trivial 1-cycle (containing 1) plus at least two known
non-trivial cycles:
  - 13-cycle: {13, 26, 33, 52, 66, 83, 104, 166, 208, 416}, length 11
  - 17-cycle: {17, 27, 34, 43, 54, 68, 86, 108, 136, 216}, length 11

Our generation classifies orbits as "converged" iff they reach exactly 1 within
max_steps. Orbits that enter non-trivial cycles cycle forever and are marked
non-convergent (they hit max_steps).

Question: among non-convergent orbits, what fraction lands in each known
non-trivial cycle vs grows unboundedly?

If a meaningful fraction lands in non-trivial cycles, the (a)-only "reach 1"
classification might miss structurally similar orbits that "settle in any
cycle". This would refine the publication-grade claim.

We sample non-convergent orbits and run Floyd's algorithm on each.

Usage:
    python 22_q5_cycle_detection.py --N 100000000 --sample 5000 --max_steps 5000
"""
import argparse
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")


# Known 5x+1 cycle "minimum members" — used to classify cycles
KNOWN_CYCLES_Q5 = {
    1: {1, 2, 3, 4, 6, 8, 16},                                    # trivial
    13: {13, 26, 33, 52, 66, 83, 104, 166, 208, 416},             # 13-cycle
    17: {17, 27, 34, 43, 54, 68, 86, 108, 136, 216},              # 17-cycle
}


def step_q(x, q):
    return x // 2 if x % 2 == 0 else q * x + 1


def floyd_detect(n0, q, max_steps, max_value):
    """Return (status, info). status in {'reach_1', 'cycle', 'overflow', 'maxsteps'}."""
    if n0 == 1:
        return "reach_1", None
    tortoise = step_q(n0, q)
    if tortoise == 1:
        return "reach_1", None
    hare = step_q(step_q(n0, q), q) if step_q(n0, q) != 1 else 1
    if hare == 1:
        return "reach_1", None
    steps = 1
    while tortoise != hare:
        tortoise = step_q(tortoise, q)
        if tortoise == 1:
            return "reach_1", None
        if tortoise > max_value:
            return "overflow", tortoise
        hare = step_q(hare, q)
        if hare == 1:
            return "reach_1", None
        if hare > max_value:
            return "overflow", hare
        hare = step_q(hare, q)
        if hare == 1:
            return "reach_1", None
        if hare > max_value:
            return "overflow", hare
        steps += 1
        if steps >= max_steps:
            return "maxsteps", None
    # Cycle detected; find a cycle member and return the minimum
    cycle_start = tortoise
    members = {cycle_start}
    x = step_q(cycle_start, q)
    while x != cycle_start:
        members.add(x)
        x = step_q(x, q)
        if len(members) > 10000:
            break
    cycle_min = min(members)
    return "cycle", cycle_min


def classify_cycle(cycle_min):
    """Map cycle_min to canonical cycle name."""
    for name, mems in KNOWN_CYCLES_Q5.items():
        if cycle_min in mems:
            return name
    return f"unknown_{cycle_min}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--N", type=int, default=100_000_000)
    ap.add_argument("--sample", type=int, default=5000)
    ap.add_argument("--max_steps", type=int, default=5000)
    ap.add_argument("--max_value", type=float, default=1e30)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    print(f"[load] q={args.q}, N={args.N:,}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet").filter(
        ~pl.col("converged"))
    n_nonconv = len(df)
    print(f"        non-convergent orbits: {n_nonconv:,}", flush=True)

    rng = np.random.default_rng(args.seed)
    sample_size = min(args.sample, n_nonconv)
    idx = rng.choice(n_nonconv, sample_size, replace=False)
    ns = df["n"].to_numpy()[idx]

    print(f"[detect] running Floyd on {sample_size} non-convergent orbits ...", flush=True)
    print(f"         max_steps={args.max_steps}, max_value={args.max_value:.1e}", flush=True)

    status_counter = Counter()
    cycle_counter = Counter()
    max_value = int(args.max_value)
    for nn in ns:
        status, info = floyd_detect(int(nn), args.q, args.max_steps, max_value)
        status_counter[status] += 1
        if status == "cycle":
            cyc_name = classify_cycle(info)
            cycle_counter[cyc_name] += 1

    print()
    print(f"=== Results ({sample_size} non-convergent orbits sampled) ===")
    print(f"{'status':>15}  count    fraction")
    for status, count in status_counter.most_common():
        frac = count / sample_size
        print(f"{status:>15}  {count:>5}  {frac:>9.4f}")

    print()
    print(f"=== Cycle breakdown (among orbits that entered a cycle) ===")
    n_in_cycle = status_counter.get("cycle", 0)
    if n_in_cycle == 0:
        print(f"  No orbits entered detectable cycles within {args.max_steps} steps.")
        print(f"  All non-convergent orbits either overflowed or ran out of steps.")
    else:
        print(f"{'cycle':>20}  count    fraction_of_cycles")
        for cyc, count in cycle_counter.most_common():
            frac = count / n_in_cycle
            print(f"{cyc:>20}  {count:>5}  {frac:>18.4f}")

    print()
    print(f"=== Overall fraction landing in non-trivial cycles ===")
    n_nontrivial = sum(c for cyc, c in cycle_counter.items() if not cyc.startswith("1"))
    n_overflow = status_counter.get("overflow", 0)
    n_maxsteps = status_counter.get("maxsteps", 0)
    print(f"  Non-trivial cycle landings: {n_nontrivial} of {sample_size} non-convergent ({100*n_nontrivial/sample_size:.2f}%)")
    print(f"  Overflowed (likely divergent): {n_overflow} ({100*n_overflow/sample_size:.2f}%)")
    print(f"  Hit max_steps without resolution: {n_maxsteps} ({100*n_maxsteps/sample_size:.2f}%)")
    print()
    print(f"Implication for paper:")
    if n_nontrivial / sample_size > 0.05:
        print(f"  Non-trivial cycle landings are >5% of non-convergent orbits.")
        print(f"  If these orbits have similar prefix-complexity structure to the (a)-class,")
        print(f"  excluding them might bias the conv_rate analysis. Worth checking j-distribution.")
    else:
        print(f"  Non-trivial cycle landings are <5% of non-convergent orbits.")
        print(f"  The (a)-only reach-1 analysis is methodologically clean. Other 'non-convergent'")
        print(f"  orbits are predominantly truly divergent (unbounded growth to overflow).")


if __name__ == "__main__":
    main()
