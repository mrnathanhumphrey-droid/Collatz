"""
duality_followup_check.py
=========================
Two follow-up checks on the Agent 2 vs Agent 3 D_n asymmetry:

CHECK 1 — Matched-N: rebuild Agent 2's 3x+1 inverse tree from 1, and at each
  depth n compute D_n^matched(k) using only the first |V_n^Agent3| vertices
  (sorted by integer value, smallest first).

CHECK 2 — Per-vertex normalized: D~_n(k) := D_n(k) / |V_n|^2 for both Agents.

Sample-size-artifact ⇒ matched-N gives Agent 3-like values.
Structural fingerprint ⇒ matched-N still close to original Agent 2.
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8")

E_MAX = 30
K_VALUES = [1, 2, 3, 4, 5]
N_VALUES = [0, 1, 2, 3, 4, 5, 6]


# --------------------------------------------------------------------------- #
# Agent 2's inverse tree builder + D_n closed form                            #
# --------------------------------------------------------------------------- #

def predecessors_3xplus1(y, e_max=E_MAX):
    if y % 3 == 0:
        return []
    out = []
    e_start = 2 if (y % 3) == 1 else 1
    pw = (1 << e_start) * y
    e = e_start
    while e <= e_max:
        x = (pw - 1) // 3
        out.append((e, x))
        pw <<= 2
        e += 2
    return out


def build_inverse_tree(max_depth):
    levels = [[1]]
    for d in range(1, max_depth + 1):
        prev = levels[-1]
        cur = []
        for y in prev:
            if y % 3 == 0:
                continue
            for e, x in predecessors_3xplus1(y, E_MAX):
                if d == 1 and x == 1:
                    continue
                cur.append(x)
        levels.append(cur)
    return levels


def D_n_from_vertices(vertices, k):
    """Compute D_n(k) = 3^k * sum mu(r)^2 - 3^{k-1} * sum Q(s)^2
    where mu(r) = (#vertices ≡ r mod 3^k AND r coprime to 3) / |vertices|.
    Q(s) = sum_{r ≡ s mod 3^{k-1}, r coprime to 3} mu(r).

    Returns Fraction.
    """
    N = 3 ** k
    if not vertices:
        return Fraction(0)
    V = len(vertices)
    # Histogram of residues mod 3^k, only r coprime to 3 (i.e., r mod 3 != 0)
    cnt_full = Counter(v % N for v in vertices)
    cnt_coprime = {r: c for r, c in cnt_full.items() if r % 3 != 0}

    sum_mu_sq = Fraction(0)
    for r, c in cnt_coprime.items():
        sum_mu_sq += Fraction(c, V) ** 2

    # Q(s) = mass at residue s mod 3^{k-1} via lifts from coprime r
    if k == 1:
        # Q lives on Z/1; sum Q^2 = (sum mu)^2
        sum_Q_sq = (sum(Fraction(c, V) for c in cnt_coprime.values())) ** 2
    else:
        N_lower = 3 ** (k - 1)
        Q = {}
        for r, c in cnt_coprime.items():
            s = r % N_lower
            Q[s] = Q.get(s, Fraction(0)) + Fraction(c, V)
        sum_Q_sq = sum(q * q for q in Q.values())

    return Fraction(3 ** k) * sum_mu_sq - Fraction(3 ** (k - 1)) * sum_Q_sq


# --------------------------------------------------------------------------- #
# Load Agent 2 and Agent 3 data                                                #
# --------------------------------------------------------------------------- #

def load_agent2_full():
    """Full (n, k, |V_n|, D_n)."""
    out = {}
    with open(r"C:\Collatz\result_inverse_tree_residue.csv", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n = int(row["n"]); k = int(row["k"])
            V = int(row["total_vertices"])
            D = Fraction(int(row["D_n_k_num"]), int(row["D_n_k_den"]))
            out[(n, k)] = (V, D)
    return out


def load_agent3_root(root):
    out = {}
    with open(fr"C:\Collatz\agent3_Dn_root_{root}.csv", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n = int(row["n"]); k = int(row["k"])
            V = int(row["vertex_count"])
            D = float(row["D_n_k"])
            out[(n, k)] = (V, D)
    return out


def load_agent3_total():
    out = {}
    with open(r"C:\Collatz\agent3_Dn_total.csv", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n = int(row["n"]); k = int(row["k"])
            d_w = float(row["D_n_k_weighted"])
            out[(n, k)] = d_w
    return out


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def main():
    print("=" * 78)
    print("Follow-up: matched-N and per-vertex-normalized duality checks")
    print("=" * 78)
    print()

    A2 = load_agent2_full()
    A3_root1 = load_agent3_root(1)
    A3_root5 = load_agent3_root(5)
    A3_root17 = load_agent3_root(17)
    A3_total = load_agent3_total()

    # CHECK 1: rebuild Agent 2's tree, sort vertices at each depth, truncate to
    # Agent 3's per-root and total vertex counts, recompute D.
    print("CHECK 1: Matched-N — Agent 2's first M_n vertices (sorted by value)")
    print("-" * 78)
    print()

    print("  Building Agent 2's 3x+1 inverse tree to depth 6...")
    import time
    t0 = time.time()
    levels = build_inverse_tree(6)
    elapsed = time.time() - t0
    sizes_A2 = [len(L) for L in levels]
    print(f"  Built {sum(sizes_A2)} total vertices in {elapsed:.1f}s")
    print(f"  Per-depth sizes: {sizes_A2}")
    print()

    # Sort each level by integer value
    levels_sorted = [sorted(L) for L in levels]

    # Agent 3 root 1 vertex counts at each depth:
    A3_root1_sizes = [A3_root1[(n, 1)][0] for n in N_VALUES]
    A3_root5_sizes = [A3_root5[(n, 1)][0] for n in N_VALUES]
    A3_root17_sizes = [A3_root17[(n, 1)][0] for n in N_VALUES]

    print(f"  Agent 3 root 1 sizes:  {A3_root1_sizes}")
    print(f"  Agent 3 root 5 sizes:  {A3_root5_sizes}")
    print(f"  Agent 3 root 17 sizes: {A3_root17_sizes}")
    print()

    # Compute matched D on each truncation
    print("  Matched-N D_n (Agent 2 truncated to Agent 3 root-1 sizes), sorted-by-value")
    print("-" * 78)
    print(f"  {'(n,k)':>8}  {'A2 full':>14}  {'A2 matched':>14}  {'A3 root 1':>14}"
          f"  {'matched/full':>14}  {'matched/A3':>14}")
    rows_csv = []
    for n in N_VALUES:
        M_n = A3_root1_sizes[n]
        for k in K_VALUES:
            full_D = A2[(n, k)][1]
            verts_truncated = levels_sorted[n][:M_n]
            matched_D = D_n_from_vertices(verts_truncated, k)
            a3_D = A3_root1[(n, k)][1]
            r_full = float(matched_D / full_D) if full_D != 0 else float("inf")
            r_a3 = float(matched_D) / a3_D if a3_D > 0 else float("inf")
            print(f"  ({n},{k}):  {float(full_D):>14.4e}  "
                  f"{float(matched_D):>14.4e}  {a3_D:>14.4e}  "
                  f"{r_full:>14.4f}  {r_a3:>14.4f}")
            rows_csv.append({
                "n": n, "k": k,
                "A2_full": float(full_D),
                "A2_matched_to_A3root1": float(matched_D),
                "A3_root1": a3_D,
                "matched_to_full_ratio": r_full,
                "matched_to_A3_ratio": r_a3,
            })
        print()

    # CHECK 2: per-vertex normalized D~_n = D_n / |V_n|^2
    print()
    print("CHECK 2: Per-vertex normalized  D~_n(k) := D_n(k) / |V_n|^2")
    print("-" * 78)
    print(f"  {'(n,k)':>8}  {'A2 |V|':>10}  {'A2 D~':>14}  {'A3 r1 |V|':>10}  {'A3 r1 D~':>14}"
          f"  {'A3/A2 ratio':>14}")
    for n in N_VALUES:
        for k in K_VALUES:
            V_a2 = A2[(n, k)][0]
            D_a2 = float(A2[(n, k)][1])
            V_a3 = A3_root1[(n, k)][0]
            D_a3 = A3_root1[(n, k)][1]
            tilde_a2 = D_a2 / (V_a2 ** 2) if V_a2 > 0 else float("nan")
            tilde_a3 = D_a3 / (V_a3 ** 2) if V_a3 > 0 else float("nan")
            ratio = tilde_a3 / tilde_a2 if tilde_a2 > 0 else float("inf")
            print(f"  ({n},{k}):  {V_a2:>10}  {tilde_a2:>14.4e}  "
                  f"{V_a3:>10}  {tilde_a3:>14.4e}  {ratio:>14.4f}")
        print()

    # CSV
    out_csv = r"C:\Collatz\duality_followup_data.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        fields = ["n", "k", "A2_full", "A2_matched_to_A3root1", "A3_root1",
                   "matched_to_full_ratio", "matched_to_A3_ratio"]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows_csv)
    print(f"\n[csv: {out_csv}]")


if __name__ == "__main__":
    main()
