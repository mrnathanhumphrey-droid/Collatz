"""
agent3_inverse_tree_3xm1_Dn.py
==============================
Agent 3 task: D_n(k) = depth-n inverse-tree Plancherel mass at modulus 3^k for
the (x+1)/3 inverse map of the 3x-1 forward dynamics.

Definitions (per the brief):
  Forward 3x-1 (Syracuse):    m -> (3m - 1) / 2^v   for odd m, v = v_2(3m - 1)
  Inverse:                    g_+(y) = (2^e · y + 1) / 3,
                              valid when 2^e · y ≡ -1 (mod 3) and result odd.

  Vertices: odd integers; vertices with y ≡ 0 (mod 3) have NO preimages.

Three roots, one per known 3x-1 cycle: {1}, {5, 7}, {17, 25, 37, 55, 41, 61, 91}.
For each root R, BFS the inverse tree to depth max_depth = 6, excluding
cycle members from preimage expansions (so the cycle does not loop) and
capping each preimage at m_pred <= N_MAX so vertex count stays bounded.

For each depth n in {0,1,...,6} and each k in {1,...,5}:
  N_n = #{depth-n vertices}
  mu_{n,k}(r) = #{depth-n vertices ≡ r mod 3^k} / N_n        (r in Z/3^k)
  mu_hat_{n,k}(xi) = sum_r mu(r) * exp(-2*pi*i * r*xi / 3^k)
  D_n(k) = sum_{xi in Z/3^k, gcd(xi, 3) = 1} |mu_hat_{n,k}(xi)|^2

Sanity (n = 0): single vertex (the root r_0). mu = delta_{r_0}, |mu_hat(xi)| = 1
for all xi. D_0(k) = phi(3^k) = 2 * 3^(k-1).

Outputs:
  agent3_Dn_root_1.csv,  agent3_Dn_root_5.csv,  agent3_Dn_root_17.csv
  agent3_Dn_total.csv  (basin-density-weighted sum)
  agent3_basin_densities.csv
  agent3_inverse_tree_3xm1_Dn_log.txt
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np
from numba import njit

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
N_MAX = 10**8          # m_pred cap during BFS (tunable)
MAX_DEPTH = 6
N_BASIN_SAMPLE = 1_000_000  # forward classification budget
FORWARD_STEP_CAP = 5000
FORWARD_VAL_CAP = 10**18
DEPTH_VERTEX_CAP = 1_000_000

ROOTS = [1, 5, 17]

# ----------------------- 3x-1 Syracuse step + cycle catalog -----------------

@njit(cache=True)
def syracuse_3xm1(m):
    x = 3 * m - 1
    while (x & 1) == 0:
        x >>= 1
    return x


def build_cycle(root, max_len=100):
    """Walk Syracuse-3x-1 from `root` until it returns; return tuple of members."""
    seen = []
    m = root
    for _ in range(max_len * 2):
        seen.append(m)
        m = int(syracuse_3xm1(m))
        if m == root:
            return tuple(seen)
    raise RuntimeError(f"cycle of {root} did not close in {max_len} steps")


# ----------------------- Inverse-tree BFS (exact int) -----------------------

def inverse_tree_3xm1(root, cycle_members, max_depth=MAX_DEPTH,
                      n_max=N_MAX, vertex_cap=DEPTH_VERTEX_CAP):
    """BFS the inverse tree of forward 3x-1 starting at `root`.

    At each vertex y (with y % 3 != 0), enumerate preimages
    g_+(y) = (2^e * y + 1) / 3 for e starting at v_min and stepping by 2,
    while m_pred <= n_max.
    Skip preimages that are in `cycle_members` (back-edges into cycle).

    Returns list of layers [layer_0, layer_1, ..., layer_max_depth].
    Reports if any depth exceeds vertex_cap (for log-only; we keep the layer
    intact and let the user decide).
    """
    cycle_set = set(cycle_members)
    layers = [[root]]
    capped_depths = []
    for depth in range(1, max_depth + 1):
        next_layer = []
        for n in layers[-1]:
            mod3 = n % 3
            if mod3 == 0:
                continue
            v = 1 if mod3 == 1 else 2  # v_min by parity-of-3 constraint
            n_2v = n << v
            while True:
                num = n_2v + 1
                if num % 3 != 0:
                    # Should not happen given v_min/step 2; safety bail.
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
            capped_depths.append((depth, len(next_layer)))
    return layers, capped_depths


# ----------------------- D_n(k) computation -----------------------

def D_n_k(vertices, k):
    """Plancherel mass at modulus 3^k of empirical residue distribution.

    D_n(k) = sum_{xi coprime to 3 in Z/3^k} |mu_hat(xi)|^2
      mu_hat(xi) = sum_r mu(r) * exp(-2 pi i r xi / 3^k)
    """
    if not vertices:
        return float("nan")
    M = 3 ** k
    N = len(vertices)
    counts = np.zeros(M, dtype=np.int64)
    for v in vertices:
        counts[v % M] += 1
    mu = counts.astype(np.float64) / N
    mu_hat = np.fft.fft(mu)
    abs_sq = mu_hat.real ** 2 + mu_hat.imag ** 2
    total = 0.0
    for xi in range(M):
        if xi % 3 != 0:
            total += abs_sq[xi]
    return float(total)


def D_n_k_table(layers, k_list):
    rows = []
    for n, layer in enumerate(layers):
        for k in k_list:
            d = D_n_k(layer, k)
            rows.append({
                "n": n,
                "k": k,
                "vertex_count": len(layer),
                "D_n_k": d,
                "phi_3k": 2 * 3 ** (k - 1),
            })
    return rows


# ----------------------- Basin densities -----------------------

@njit(cache=True)
def classify_with_cycle_set_max(start, cycle_min_arr, val_cap, step_cap):
    """Walk forward 3x-1; return cycle_min_arr index (1-based) of the cycle
    that 'start' lands in, or 0 if val_cap or step_cap exceeded.

    cycle_min_arr is a sorted np.int64 array of MIN members of each cycle.
    We declare landing when m equals any min; we trust that if m is in any
    cycle then iterating Syracuse from m will hit the min within
    cycle-length steps. To detect any-cycle-membership reliably, we walk
    Syracuse until m equals a min of any known cycle OR caps are hit.
    """
    m = start
    for _ in range(step_cap):
        for i in range(len(cycle_min_arr)):
            if m == cycle_min_arr[i]:
                return i + 1  # 1-based index
        if m > val_cap:
            return 0
        x = 3 * m - 1
        while (x & 1) == 0:
            x >>= 1
        m = x
    return 0


def basin_densities(cycles_info, n_sample):
    """cycles_info: list of (root, cycle_members_tuple). Returns dict
    {root: density, "unresolved": density}."""
    cycle_mins = np.array([min(cm) for (_, cm) in cycles_info],
                          dtype=np.int64)
    print(f"  Basin density: classifying {n_sample:,} odd starts...")
    print(f"  cycle minima used as detection set: {cycle_mins.tolist()}")
    # Sample odd starts in a wide range
    rng = np.random.default_rng(20260504)
    starts = (2 * rng.integers(1, 10**12, size=n_sample, dtype=np.int64)
              + 1)
    counts = np.zeros(len(cycles_info) + 1, dtype=np.int64)  # idx 0 = unresolved
    t0 = time.time()
    for s in starts:
        c = classify_with_cycle_set_max(int(s), cycle_mins,
                                        FORWARD_VAL_CAP, FORWARD_STEP_CAP)
        counts[c] += 1
    elapsed = time.time() - t0
    print(f"    classified in {elapsed:.1f}s")
    total = counts.sum()
    out = {}
    for i, (root, _) in enumerate(cycles_info):
        out[root] = counts[i + 1] / total
    out["unresolved"] = counts[0] / total
    return out


# ----------------------- main -----------------------

def main():
    print("=" * 78)
    print("Agent 3: inverse-tree Plancherel mass D_n(k) for 3x-1")
    print("=" * 78)
    print(f"  N_MAX (preimage cap) = {N_MAX:,}")
    print(f"  max_depth = {MAX_DEPTH}, k = 1..5")
    print(f"  vertex_cap (per-depth warn) = {DEPTH_VERTEX_CAP:,}")
    print()

    # Build cycles for each root
    print("Verifying 3x-1 Syracuse cycles for the three roots:")
    cycles_info = []
    for r in ROOTS:
        cyc = build_cycle(r)
        print(f"  root {r}: cycle = {cyc}")
        cycles_info.append((r, cyc))
    print()

    k_list = [1, 2, 3, 4, 5]

    # Per-root inverse trees
    per_root_tables = {}
    layer_counts = {}
    for root, cyc in cycles_info:
        print(f"--- Root {root}, cycle = {cyc} ---")
        t0 = time.time()
        layers, capped = inverse_tree_3xm1(root, cyc, MAX_DEPTH, N_MAX)
        elapsed = time.time() - t0
        print(f"  BFS depth={MAX_DEPTH}, N_MAX={N_MAX:,}: {elapsed:.2f}s")
        for d, layer in enumerate(layers):
            warn = "" if len(layer) <= DEPTH_VERTEX_CAP else "  [WARN > 1M cap]"
            print(f"    depth {d}: {len(layer):>10,} vertices{warn}")
        rows = D_n_k_table(layers, k_list)
        per_root_tables[root] = rows
        layer_counts[root] = [len(L) for L in layers]
        # Print compact
        print(f"  D_n(k) for root {root}:")
        print(f"    {'n':>3} | "
              + " | ".join(f"k={k}".rjust(14) for k in k_list))
        for n in range(MAX_DEPTH + 1):
            row_strs = []
            for k in k_list:
                d = next(r["D_n_k"] for r in rows if r["n"] == n
                         and r["k"] == k)
                row_strs.append(f"{d:14.6e}" if not math.isnan(d)
                                else f"{'nan':>14}")
            print(f"    {n:>3} | " + " | ".join(row_strs))
        print()
        out = os.path.join(OUT_DIR, f"agent3_Dn_root_{root}.csv")
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["n", "k", "vertex_count",
                                              "D_n_k", "phi_3k"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"    saved {out}")
        print()

    # Basin densities
    print("=" * 78)
    print("Basin density estimation")
    print("=" * 78)
    densities = basin_densities(cycles_info, N_BASIN_SAMPLE)
    print(f"  basin densities: {densities}")
    print()

    out = os.path.join(OUT_DIR, "agent3_basin_densities.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["root", "basin_density"])
        for r in ROOTS:
            w.writerow([r, densities[r]])
        w.writerow(["unresolved", densities["unresolved"]])
    print(f"  saved {out}")

    # Total D_n(k) = weighted sum over roots
    print()
    print("=" * 78)
    print("Total D_n(k) = sum over roots, weighted by basin density")
    print("(Renormalized by sum of resolved basin weights so unresolved doesn't")
    print(" contaminate; weight_used column reports the renormalization factor.)")
    print("=" * 78)
    total_rows = []
    weight_resolved = sum(densities[r] for r in ROOTS)
    for n in range(MAX_DEPTH + 1):
        for k in k_list:
            tot = 0.0
            for root in ROOTS:
                d = next(r["D_n_k"] for r in per_root_tables[root]
                         if r["n"] == n and r["k"] == k)
                w = densities[root]
                if not math.isnan(d):
                    tot += w * d
            tot_norm = (tot / weight_resolved) if weight_resolved > 0 \
                       else float("nan")
            total_rows.append({
                "n": n, "k": k,
                "D_n_k_weighted": tot,
                "D_n_k_renorm": tot_norm,
                "weight_used": weight_resolved,
            })
    print(f"  {'n':>3} | "
          + " | ".join(f"k={k}".rjust(14) for k in k_list))
    for n in range(MAX_DEPTH + 1):
        row_strs = []
        for k in k_list:
            t = next(r["D_n_k_renorm"] for r in total_rows
                     if r["n"] == n and r["k"] == k)
            row_strs.append(f"{t:14.6e}" if not math.isnan(t)
                            else f"{'nan':>14}")
        print(f"  {n:>3} | " + " | ".join(row_strs))
    print()
    out = os.path.join(OUT_DIR, "agent3_Dn_total.csv")
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
