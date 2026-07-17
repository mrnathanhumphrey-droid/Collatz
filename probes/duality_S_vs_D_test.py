"""
duality_S_vs_D_test.py
======================
Forward-backward duality test between forward S_n^(3x+/-1) (= identical by my
proved chain symmetry) and inverse-tree D_n^((x-/+1)/3) tables from Agents 2/3.

Tests:
  1. D_n(k) * S_k = const?
  2. D_n(k) + S_k = const?
  3. D_n(k) / S_k stable across n?
  4. D_n(n) vs S_n: scaling relation along the diagonal?
  5. D_{n+1}(k)/D_n(k) "q/3 analogue" — what universal ratio (if any)?

Also tests the basin fingerprint: does Agent 3's three-basin total D_n
differ from Agent 2's single-basin D_n at matched (n, k)?
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")

# Forward S_k^(3x+1) values (identical to S_k^(3x-1) by my proved chain symmetry)
S_k = {
    1: Fraction(2, 3),
    2: Fraction(10, 21),
    3: Fraction(31370, 67963),
    4: Fraction(143195649659456490, 308468774477179141),
    # k=5: from R77.7 cache eps_5 = -1.1517469151e-3, S_5 = eps_5 + 7/15 ≈ 0.4655
    # rougher; use float for k=5
}


def load_agent2():
    """Load result_inverse_tree_residue.csv: 3x+1 single-basin D_n(k)."""
    out = {}
    with open(r"C:\Collatz\result_inverse_tree_residue.csv", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n = int(row["n"]); k = int(row["k"])
            num = int(row["D_n_k_num"]); den = int(row["D_n_k_den"])
            out[(n, k)] = Fraction(num, den)
    return out


def load_agent3_total():
    """Load agent3_Dn_total.csv: 3x-1 three-basin weighted D_n(k)."""
    out = {}
    with open(r"C:\Collatz\agent3_Dn_total.csv", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n = int(row["n"]); k = int(row["k"])
            d_w = float(row["D_n_k_weighted"])
            d_renorm = float(row["D_n_k_renorm"])
            out[(n, k)] = (d_w, d_renorm)
    return out


def load_agent3_per_root():
    """Load per-root tables for basin fingerprint test."""
    roots = {}
    for root in [1, 5, 17]:
        path = fr"C:\Collatz\agent3_Dn_root_{root}.csv"
        roots[root] = {}
        with open(path, encoding="utf-8") as fh:
            r = csv.DictReader(fh)
            for row in r:
                n = int(row["n"]); k = int(row["k"])
                d = float(row["D_n_k"])
                roots[root][(n, k)] = d
    return roots


def main():
    A2 = load_agent2()
    A3_tot = load_agent3_total()
    A3_per = load_agent3_per_root()

    print("=" * 78)
    print("Agent 2 (3x+1 inverse, single basin) vs Agent 3 (3x-1 inverse, three basins)")
    print("=" * 78)
    print()

    # 1. Direct numerical comparison D_2 (Agent 2) vs D_3_tot at matched (n, k)
    print("1. D_n^(3x+1)(k) vs D_n^(3x-1)(k) total — basin asymmetry direct check")
    print("-" * 78)
    print(f"  {'(n,k)':>8}  {'A2 (3x+1)':>15}  {'A3 (3x-1 total)':>17}  {'ratio A3/A2':>12}")
    for n in range(0, 7):
        for k in [2, 3, 4, 5]:
            if (n, k) in A2 and (n, k) in A3_tot:
                a2 = float(A2[(n, k)])
                a3 = A3_tot[(n, k)][0]  # weighted
                ratio = a3/a2 if a2 != 0 else float("inf")
                print(f"  ({n},{k}):  {a2:>15.6e}  {a3:>17.6e}  {ratio:>12.4f}")
    print()

    # 2. Duality candidates against forward S_k
    print("2. Duality test: forward S_k vs inverse D_n(k)")
    print("-" * 78)
    S_floats = {k: float(v) for k, v in S_k.items()}
    S_floats[5] = 0.4655315306  # approximate (k=5 not in exact cache via Fraction)
    print(f"  Forward S_k: {S_floats}")
    print()
    print(f"  {'(n,k)':>8}  {'D_n(k) [A2]':>15}  {'S_k':>12}  "
          f"{'D·S':>12}  {'D+S':>12}  {'D/S':>14}")
    for n in range(0, 7):
        for k in [1, 2, 3, 4, 5]:
            if (n, k) not in A2 or k not in S_floats: continue
            d = float(A2[(n, k)])
            s = S_floats[k]
            print(f"  ({n},{k}):  {d:>15.6e}  {s:>12.6f}  "
                  f"{d*s:>12.6e}  {d+s:>12.6f}  {d/s:>14.6e}")
    print()

    # 3. Diagonal D_n(n) vs S_n
    print("3. Diagonal D_n(n) vs S_n")
    print("-" * 78)
    print(f"  {'n':>3}  {'D_n(n) [A2]':>15}  {'S_n':>12}  {'ratio D/S':>14}")
    for n in [1, 2, 3, 4, 5]:
        if (n, n) in A2 and n in S_floats:
            d = float(A2[(n, n)])
            s = S_floats[n]
            print(f"  {n:>3}  {d:>15.6e}  {s:>12.6f}  {d/s:>14.6e}")
    print()

    # 4. q/3 analogue for inverse: ratio D_{n+1}(k) / D_n(k)
    print("4. Inverse 'q/3 analogue': D_{n+1}(k) / D_n(k) ratios across n")
    print("-" * 78)
    print(f"  {'k':>3}  ratios D_{{n+1}}(k)/D_n(k) for n=0..5")
    for k in [2, 3, 4, 5]:
        ratios = []
        for n in range(0, 6):
            if (n, k) in A2 and (n+1, k) in A2:
                d_n = float(A2[(n, k)])
                d_np1 = float(A2[(n+1, k)])
                if d_n > 0:
                    ratios.append(d_np1/d_n)
        print(f"  k={k}:  {[f'{r:.4f}' for r in ratios]}")
    print()
    print("  Compare to forward S_{k+1}/S_k = q/3 = 1 (q=3, S converges)")
    for k in [1, 2, 3, 4]:
        if (k+1) in S_floats and k in S_floats:
            r = S_floats[k+1] / S_floats[k]
            print(f"  S_{k+1}/S_{k} = {r:.4f}")
    print()

    # 5. Basin fingerprint: per-root D_n(k) at fixed (n, k)
    print("5. Basin fingerprint: per-root D_n^(root)(k) at (n=3, k=3) and (n=4, k=4)")
    print("-" * 78)
    for (n, k) in [(3, 3), (4, 4), (5, 5)]:
        print(f"  (n={n}, k={k}):")
        for root in [1, 5, 17]:
            d = A3_per[root].get((n, k), float("nan"))
            print(f"    root {root:>2}:  D = {d:.6e}")
        if (n, k) in A3_tot:
            print(f"    weighted total:  {A3_tot[(n,k)][0]:.6e}")
        if (n, k) in A2:
            print(f"    Agent 2 (3x+1):  {float(A2[(n,k)]):.6e}")
        print()

    # 6. Vertex growth comparison
    print("6. Vertex growth |V_n| comparison")
    print("-" * 78)
    print("  Agent 2 (3x+1 single basin): |V_n| = 1, 14, 135, 1350, 13500, 135000, 1350000")
    print("    growth ratio: ~14, ~9.6, 10, 10, 10, 10  (asymptotically 10)")
    print()
    print("  Agent 3 (3x-1 root 1): |V_n| = 1, 13, 56, 189, 459, 1061, 2247")
    print("    growth ratio: 13, 4.3, 3.4, 2.4, 2.3, 2.1  (DECAYING — basin/cap-bounded)")


if __name__ == "__main__":
    main()
