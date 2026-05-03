"""
Experiment 36 — q=5 fourth-cycle targeted search at m ≡ 33 mod 40.

Santos (2021, Table 1) lists three known cycles for the qx+1 problem at q=5
(smallest members {1, 13, 17}) and conjectures that at most one more cycle
exists, with smallest member ≡ 33 mod 40. Our exp 29 found no fourth cycle
in [1, 10^8]. This experiment extends the search to [10^8, 10^10] by
restricting to m ≡ 33 mod 40 and running Floyd's tortoise-hare cycle
detection on each.

m ≡ 33 mod 40 yields one start every 40 integers → ~2.5×10^8 starts in
[10^8, 10^10]. Most diverge fast at q=5 (orbit grows by ~5/2 per step).
Cycle landings classified by smallest member; any smallest ∉ {1, 13, 17}
is a new cycle.

Usage:
    python 36_q5_fourth_cycle_search.py --m_min 100000000 --m_max 10000000000
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

S_TRIVIAL  = np.int8(0)
S_NONTRIV  = np.int8(1)
S_DIVERGE  = np.int8(2)
S_TIMEOUT  = np.int8(3)


@njit(cache=True)
def step_q5(x):
    if x % 2 == 0:
        return x // 2
    else:
        return 5 * x + 1


@njit(cache=True)
def classify_one(start, max_value, max_steps):
    """Floyd's tortoise-hare for q=5. Returns (status, smallest_member, cycle_len)."""
    safe_cap = max_value // 5
    t = start
    h = start
    steps = 0
    while steps < max_steps:
        if (t % 2 == 1) and (t > safe_cap):
            return S_DIVERGE, np.int64(0), np.int64(0)
        t = step_q5(t)
        if (h % 2 == 1) and (h > safe_cap):
            return S_DIVERGE, np.int64(0), np.int64(0)
        h = step_q5(h)
        if h > max_value:
            return S_DIVERGE, np.int64(0), np.int64(0)
        if (h % 2 == 1) and (h > safe_cap):
            return S_DIVERGE, np.int64(0), np.int64(0)
        h = step_q5(h)
        if h > max_value:
            return S_DIVERGE, np.int64(0), np.int64(0)
        steps += 1
        if t == h:
            break
    if t != h:
        return S_TIMEOUT, np.int64(0), np.int64(0)

    # Locate cycle entry
    t2 = start
    while t2 != h:
        t2 = step_q5(t2)
        h  = step_q5(h)
    cycle_entry = t2

    # Walk one full cycle, find smallest
    smallest = cycle_entry
    cur = step_q5(cycle_entry)
    cycle_len = np.int64(1)
    while cur != cycle_entry:
        if cur < smallest:
            smallest = cur
        cur = step_q5(cur)
        cycle_len += 1
        if cycle_len > max_steps:
            return S_TIMEOUT, np.int64(0), np.int64(0)

    if smallest == 1:
        return S_TRIVIAL, np.int64(1), cycle_len
    return S_NONTRIV, smallest, cycle_len


@njit(parallel=True, cache=True)
def search_residue_class(m_min, m_max, residue, modulus, max_value, max_steps):
    """Search all m in [m_min, m_max] with m ≡ residue (mod modulus).
    Returns:
        n_processed, n_trivial, n_nontriv, n_diverge, n_timeout
        unique_smallest_seen (int64 array, padded with 0)
        max_cycle_length_seen
    """
    # Snap m_min to first m ≡ residue mod modulus, m >= m_min
    first_m = m_min + ((residue - m_min) % modulus)
    if first_m < m_min:
        first_m += modulus
    n_total = (m_max - first_m) // modulus + 1

    # Per-thread arrays
    n_threads = 32
    chunk = max(1, n_total // n_threads)
    n_chunks = (n_total + chunk - 1) // chunk

    cnt_processed = np.zeros(n_chunks, dtype=np.int64)
    cnt_trivial   = np.zeros(n_chunks, dtype=np.int64)
    cnt_nontriv   = np.zeros(n_chunks, dtype=np.int64)
    cnt_diverge   = np.zeros(n_chunks, dtype=np.int64)
    cnt_timeout   = np.zeros(n_chunks, dtype=np.int64)
    # Track first 100 unique smallest values per chunk (to detect new cycles)
    smallest_buf  = np.zeros((n_chunks, 100), dtype=np.int64)
    smallest_n    = np.zeros(n_chunks, dtype=np.int64)
    max_cycle_len = np.zeros(n_chunks, dtype=np.int64)
    new_cycle_starts = np.zeros((n_chunks, 100), dtype=np.int64)
    new_cycle_smallest = np.zeros((n_chunks, 100), dtype=np.int64)
    new_cycle_n = np.zeros(n_chunks, dtype=np.int64)

    for c in prange(n_chunks):
        i_start = c * chunk
        i_end = min(i_start + chunk, n_total)
        for i in range(i_start, i_end):
            m = first_m + i * modulus
            if m > m_max:
                break
            status, smallest, cycle_len = classify_one(m, max_value, max_steps)
            cnt_processed[c] += 1
            if status == S_TRIVIAL:
                cnt_trivial[c] += 1
            elif status == S_NONTRIV:
                cnt_nontriv[c] += 1
                # Check if smallest is "new" (not 1, 13, 17)
                if smallest != 1 and smallest != 13 and smallest != 17:
                    if new_cycle_n[c] < 100:
                        new_cycle_starts[c, new_cycle_n[c]] = m
                        new_cycle_smallest[c, new_cycle_n[c]] = smallest
                        new_cycle_n[c] += 1
                # Track unique smallest
                already = False
                for j in range(smallest_n[c]):
                    if smallest_buf[c, j] == smallest:
                        already = True
                        break
                if not already and smallest_n[c] < 100:
                    smallest_buf[c, smallest_n[c]] = smallest
                    smallest_n[c] += 1
                if cycle_len > max_cycle_len[c]:
                    max_cycle_len[c] = cycle_len
            elif status == S_DIVERGE:
                cnt_diverge[c] += 1
            else:
                cnt_timeout[c] += 1

    return (cnt_processed.sum(), cnt_trivial.sum(), cnt_nontriv.sum(),
            cnt_diverge.sum(), cnt_timeout.sum(),
            smallest_buf, smallest_n, max_cycle_len.max(),
            new_cycle_starts, new_cycle_smallest, new_cycle_n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m_min", type=int, default=100_000_000)
    ap.add_argument("--m_max", type=int, default=10_000_000_000)
    ap.add_argument("--residue", type=int, default=33)
    ap.add_argument("--modulus", type=int, default=40)
    ap.add_argument("--max_value", type=int, default=10**18)
    ap.add_argument("--max_steps", type=int, default=1_000_000)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[setup] q=5 cycle search on m ≡ {args.residue} mod {args.modulus}")
    print(f"        m range: [{args.m_min:,}, {args.m_max:,}]")
    print(f"        expected n_starts: ~{(args.m_max - args.m_min) // args.modulus:,}")
    print(f"        max_value: {args.max_value:.1e}")

    # Warmup
    _ = classify_one(np.int64(33), args.max_value, 1000)

    print(f"[run] Floyd cycle detection ...", flush=True)
    t0 = time.perf_counter()
    (cnt_proc, cnt_triv, cnt_nontriv, cnt_div, cnt_to,
     smallest_buf, smallest_n, max_cl,
     new_starts, new_smallest, new_n) = search_residue_class(
        args.m_min, args.m_max, args.residue, args.modulus,
        args.max_value, args.max_steps)
    t = time.perf_counter() - t0
    print(f"        done in {t:.1f}s")
    print()
    print(f"  total processed:  {cnt_proc:,}")
    print(f"  trivial cycle:    {cnt_triv:,}  ({100*cnt_triv/max(cnt_proc,1):.4f}%)")
    print(f"  non-trivial:      {cnt_nontriv:,}  ({100*cnt_nontriv/max(cnt_proc,1):.4f}%)")
    print(f"  divergent:        {cnt_div:,}  ({100*cnt_div/max(cnt_proc,1):.4f}%)")
    print(f"  timeout:          {cnt_to:,}")
    print(f"  max cycle length seen: {max_cl}")

    # Aggregate unique smallest values across chunks
    all_smallest = set()
    for c in range(smallest_buf.shape[0]):
        for j in range(smallest_n[c]):
            all_smallest.add(int(smallest_buf[c, j]))
    print(f"\n  unique cycle smallest members: {sorted(all_smallest)}")

    # Aggregate any "new" cycle landings (smallest not in {1, 13, 17})
    print()
    new_total = 0
    for c in range(new_n.shape[0]):
        new_total += new_n[c]
    if new_total == 0:
        print(f"  ===> NO new cycle landings found. The conjectured 4th cycle has")
        print(f"        smallest member > {args.m_max:,} = {np.log10(args.m_max):.2f} log10")
    else:
        print(f"  ===> NEW CYCLE LANDINGS: {new_total}")
        for c in range(new_n.shape[0]):
            for j in range(new_n[c]):
                print(f"      m = {new_starts[c, j]:,}  smallest = {new_smallest[c, j]:,}")
        # Save
        rows = []
        for c in range(new_n.shape[0]):
            for j in range(new_n[c]):
                rows.append({"m_start": int(new_starts[c, j]),
                             "smallest": int(new_smallest[c, j])})
        df = pl.DataFrame(rows)
        out = out_dir / f"36_q5_fourth_cycle_landings_m{args.m_min}_{args.m_max}.csv"
        df.write_csv(out)
        print(f"  [save] {out}")


if __name__ == "__main__":
    main()
