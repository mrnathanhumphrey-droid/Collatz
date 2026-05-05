"""
agent3_inverse_tree_3xm1_Dn_v2.py
================================
Agent 3 follow-up: per-basin D_n(k) tables with bumped N_MAX and pushed depth.

Differences from v1:
  - N_MAX = 10^10 (was 10^8)
  - MAX_DEPTH = 10 (was 6); early-stop when any layer exceeds 1M vertex cap
  - Reuses cached basin densities from v1 (loaded from agent3_basin_densities.csv)
  - Faster D_n(k) via numpy-vectorized residue histogram
  - Per-root CSV outputs include n=0..max_attained_depth, k=1..5

Same definitions as v1 (Plancherel mass at modulus 3^k of empirical depth-n
inverse-tree vertex residue distribution, summed over coprime xi).
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np
from numba import njit

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
N_MAX = 10**10
MAX_DEPTH = 10
DEPTH_VERTEX_CAP = 1_000_000
ROOTS = [1, 5, 17]


# ----------------------- Cycles -----------------------

@njit(cache=True)
def syracuse_3xm1(m):
    x = 3 * m - 1
    while (x & 1) == 0:
        x >>= 1
    return x


def build_cycle(root, max_len=100):
    seen = []
    m = root
    for _ in range(max_len * 2):
        seen.append(m)
        m = int(syracuse_3xm1(m))
        if m == root:
            return tuple(seen)
    raise RuntimeError(f"cycle of {root} did not close in {max_len} steps")


# ----------------------- BFS with early stop -----------------------

def inverse_tree_3xm1(root, cycle_members, max_depth=MAX_DEPTH,
                      n_max=N_MAX, vertex_cap=DEPTH_VERTEX_CAP):
    """BFS the inverse tree; stop expanding deeper if a layer exceeds vertex_cap.
    Returns layers + (early_stop_depth, last_layer_size) if early-stopped, or
    (None, ...) if completed normally.
    """
    cycle_set = set(cycle_members)
    layers = [[root]]
    early_stop = None
    for depth in range(1, max_depth + 1):
        next_layer = []
        for n in layers[-1]:
            mod3 = n % 3
            if mod3 == 0:
                continue
            v = 1 if mod3 == 1 else 2
            n_2v = n << v
            while True:
                num = n_2v + 1
                if num % 3 != 0:
                    break
                m_pred = num // 3
                if m_pred > n_max:
                    break
                if m_pred > 0 and (m_pred & 1) == 1 and m_pred not in cycle_set:
                    next_layer.append(m_pred)
                v += 2
                n_2v <<= 2
        layers.append(next_layer)
        if len(next_layer) > vertex_cap:
            early_stop = (depth, len(next_layer))
            break
    return layers, early_stop


# ----------------------- D_n(k) (numpy-vectorized) -----------------------

def D_n_k_vec(vertex_arr_int64, k):
    """Plancherel mass at modulus 3^k. vertex_arr_int64 is np.int64 array."""
    if vertex_arr_int64.size == 0:
        return float("nan")
    M = 3 ** k
    counts = np.bincount(vertex_arr_int64 % M, minlength=M).astype(np.float64)
    mu = counts / counts.sum()
    mu_hat = np.fft.fft(mu)
    abs_sq = mu_hat.real ** 2 + mu_hat.imag ** 2
    # Build coprime mask once (we recompute small cost; could cache).
    coprime_mask = np.array([(xi % 3) != 0 for xi in range(M)], dtype=bool)
    return float(abs_sq[coprime_mask].sum())


def D_n_k_table(layers, k_list):
    rows = []
    for n, layer in enumerate(layers):
        if len(layer) == 0:
            arr = np.empty(0, dtype=np.int64)
        else:
            arr = np.array(layer, dtype=np.int64)
        for k in k_list:
            d = D_n_k_vec(arr, k)
            rows.append({
                "n": n,
                "k": k,
                "vertex_count": len(layer),
                "D_n_k": d,
                "phi_3k": 2 * 3 ** (k - 1),
            })
    return rows


# ----------------------- main -----------------------

def main():
    print("=" * 78)
    print("Agent 3 v2: inverse-tree D_n(k) for 3x-1, deeper / wider")
    print("=" * 78)
    print(f"  N_MAX = {N_MAX:,}")
    print(f"  MAX_DEPTH = {MAX_DEPTH}")
    print(f"  vertex_cap (early stop trigger) = {DEPTH_VERTEX_CAP:,}")
    print()

    # Cycles
    cycles_info = []
    for r in ROOTS:
        cyc = build_cycle(r)
        print(f"  root {r}: cycle = {cyc}")
        cycles_info.append((r, cyc))
    print()

    # Try to load cached basin densities
    densities = None
    cache_path = os.path.join(OUT_DIR, "agent3_basin_densities.csv")
    if os.path.exists(cache_path):
        densities = {}
        with open(cache_path) as f:
            for row in csv.DictReader(f):
                key = row["root"]
                key = int(key) if key.isdigit() else key
                densities[key] = float(row["basin_density"])
        print(f"  loaded basin densities from cache: {densities}")
    print()

    k_list = [1, 2, 3, 4, 5]
    per_root_tables = {}

    for root, cyc in cycles_info:
        print(f"--- Root {root}, cycle = {cyc} ---")
        t0 = time.time()
        layers, early_stop = inverse_tree_3xm1(
            root, cyc, MAX_DEPTH, N_MAX, DEPTH_VERTEX_CAP)
        elapsed = time.time() - t0
        print(f"  BFS: {elapsed:.2f}s")
        for d, layer in enumerate(layers):
            warn = ("  [STOP > 1M cap]" if (early_stop and d == early_stop[0])
                    else "")
            print(f"    depth {d:>2}: {len(layer):>10,} vertices{warn}")
        attained_depth = len(layers) - 1
        print(f"  attained depth = {attained_depth}"
              + (f"  (early-stopped due to depth-{early_stop[0]} > 1M)"
                 if early_stop else "  (full)"))

        t1 = time.time()
        rows = D_n_k_table(layers, k_list)
        print(f"  D_n(k) computed in {time.time()-t1:.2f}s")
        per_root_tables[root] = rows

        # Print compact
        print(f"  D_n(k) for root {root}:")
        print(f"    {'n':>3} | "
              + " | ".join(f"k={k}".rjust(14) for k in k_list))
        for n in range(attained_depth + 1):
            row_strs = []
            for k in k_list:
                d = next(r["D_n_k"] for r in rows if r["n"] == n
                         and r["k"] == k)
                row_strs.append(f"{d:14.6e}" if not math.isnan(d)
                                else f"{'nan':>14}")
            print(f"    {n:>3} | " + " | ".join(row_strs))

        out = os.path.join(OUT_DIR, f"agent3_v2_Dn_root_{root}.csv")
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["n", "k", "vertex_count",
                                              "D_n_k", "phi_3k"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"    saved {out}")
        print()

    # Total D_n(k) weighted by basin density (use cache if available)
    if densities is None:
        print("  no basin densities cache; skipping total table")
    else:
        max_common_depth = min(
            max(r["n"] for r in per_root_tables[root])
            for root in ROOTS
        )
        weight_resolved = sum(densities[r] for r in ROOTS)
        print("=" * 78)
        print(f"Total D_n(k) (basin-weighted) for n = 0..{max_common_depth}")
        print("=" * 78)
        total_rows = []
        for n in range(max_common_depth + 1):
            for k in k_list:
                tot = 0.0
                for root in ROOTS:
                    d_lookup = next(
                        (r["D_n_k"] for r in per_root_tables[root]
                         if r["n"] == n and r["k"] == k), None)
                    if d_lookup is None or math.isnan(d_lookup):
                        continue
                    tot += densities[root] * d_lookup
                tot_norm = (tot / weight_resolved) if weight_resolved else float("nan")
                total_rows.append({
                    "n": n, "k": k,
                    "D_n_k_weighted": tot,
                    "D_n_k_renorm": tot_norm,
                    "weight_used": weight_resolved,
                })
        print(f"  {'n':>3} | "
              + " | ".join(f"k={k}".rjust(14) for k in k_list))
        for n in range(max_common_depth + 1):
            row_strs = []
            for k in k_list:
                t = next(r["D_n_k_renorm"] for r in total_rows
                         if r["n"] == n and r["k"] == k)
                row_strs.append(f"{t:14.6e}" if not math.isnan(t)
                                else f"{'nan':>14}")
            print(f"  {n:>3} | " + " | ".join(row_strs))
        out = os.path.join(OUT_DIR, "agent3_v2_Dn_total.csv")
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["n", "k", "D_n_k_weighted",
                                              "D_n_k_renorm", "weight_used"])
            w.writeheader()
            for r in total_rows:
                w.writerow(r)
        print(f"  saved {out}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
