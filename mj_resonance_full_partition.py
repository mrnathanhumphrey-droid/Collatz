"""
mj_resonance_full_partition.py — extend Result 63 derivation to FULL-tree mod-3 partition.

Key finding from initial run: m_j atomic decomposition gives only 0.3-1.4%
of empirical |μ̂(1/3)|². The {m_j} chain is NOT the dominant atom set —
in fact m_j with j ≡ 0 mod 3 (m_j ≡ 0 mod 3) are LEAVES (no preds, w_j=1).

The CORRECT closed form is the full-population partition:

  |μ̂(ξ=1/3)|² = ½[(P_0-P_1)² + (P_1-P_2)² + (P_0-P_2)²] / Z²

where P_a = Σ_{m ≡ a mod 3, m in tree} w(m), summed over ALL integers (not just m_j).

In the limit k → ∞, the FFT at j/2^k closest to 1/3 converges to this.
"""
import csv
import math
import os
import sys
import time
import numpy as np
from collections import deque, defaultdict

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_inverse_tree(max_value):
    tree = {1: {'parent': None, 'depth': 0}}
    q = deque([1])
    while q:
        m = q.popleft()
        d = tree[m]['depth']
        if m % 3 == 0:
            continue
        v_start = 2 if (m % 3 == 1) else 1
        for v in range(v_start, 64, 2):
            num = m * (1 << v) - 1
            if num <= 0:
                continue
            pred = num // 3
            if pred > max_value:
                break
            if pred & 1 == 0 or pred == m:
                continue
            if pred not in tree:
                tree[pred] = {'parent': m, 'depth': d + 1}
                q.append(pred)
    return tree


def subtree_sizes(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    by_depth_desc = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    for m in by_depth_desc:
        for c in children[m]:
            size[m] += size[c]
    return size


def compute_M_k_and_fft(integers, weights, k):
    mod = 1 << k
    M = np.zeros(mod, dtype=np.float64)
    rs = integers % mod
    np.add.at(M, rs, weights)
    Z = M.sum()
    mu_hat = np.fft.ifft(M) * mod / Z
    return mu_hat


def main():
    print("# m_j FULL-TREE partition derivation (Result 63 extension)")
    N = 1 << 22
    t0 = time.perf_counter()
    tree = build_inverse_tree(N)
    sizes = subtree_sizes(tree)
    odd_ints = np.array([m for m in tree if m & 1], dtype=np.int64)
    weights = np.array([float(sizes[m]) for m in odd_ints])
    Z = weights.sum()
    print(f"# Tree: {len(odd_ints)} odd nodes, Z = {Z:.4e}")

    # Population partition by m mod q for various q
    def pop_partition(q):
        P = defaultdict(float)
        n_per = defaultdict(int)
        for m, w in zip(odd_ints, weights):
            r = int(m) % q
            P[r] += w
            n_per[r] += 1
        return dict(P), dict(n_per)

    print(f"\n# Population partition mod 3:")
    P3, N3 = pop_partition(3)
    for r in sorted(P3.keys()):
        print(f"  m ≡ {r} mod 3:  count = {N3[r]:>8}    Σw = {P3[r]:>14.0f}    Σw/Z = {P3[r]/Z:.6f}")

    # |μ̂(1/3)|² closed form from full population partition
    a, b, c = P3.get(0, 0), P3.get(1, 0), P3.get(2, 0)
    val_pop = 0.5 * ((a-b)**2 + (b-c)**2 + (c-a)**2) / (Z**2)
    print(f"\n  |μ̂(1/3)|² closed form (full pop): {val_pop:.6e}")

    # Compare to atomic-only (m_j chain)
    print(f"\n# Atomic-only partition (m_j chain, j s.t. m_j ≤ N):")
    m_j_vals = [(j, (4**j - 1)//3) for j in range(1, 30) if (4**j - 1)//3 <= N]
    P3_at = defaultdict(float)
    for j, m in m_j_vals:
        if m in sizes:
            P3_at[m % 3] += sizes[m]
    print(f"  S_0={P3_at.get(0,0):.0f}  S_1={P3_at.get(1,0):.0f}  S_2={P3_at.get(2,0):.0f}")
    a, b, c = P3_at.get(0, 0), P3_at.get(1, 0), P3_at.get(2, 0)
    val_at = 0.5 * ((a-b)**2 + (b-c)**2 + (c-a)**2) / (Z**2)
    print(f"  |μ̂_atomic(1/3)|² closed form: {val_at:.6e}")

    print(f"\n# Population/atomic ratio: {val_pop/val_at:.2f}x")
    print(f"  → atomic accounts for {val_at/val_pop*100:.2f}% of full-pop closed-form")

    # Verify: empirical at large k
    print(f"\n# Empirical FFT closest to ξ=1/3 across k:")
    print(f"  {'k':>3}  {'j':>8}  {'j/2^k':>10}  {'|μ̂(j/2^k)|²':>14}  {'closed form':>12}")
    for k in [10, 12, 14, 16, 18, 20]:
        mu = compute_M_k_and_fft(odd_ints, weights, k)
        # Closest to 1/3
        target = 1.0 / 3
        mod = 1 << k
        j = round(target * mod)
        if j % 2 == 0:
            j -= 1  # prefer odd j
        emp = abs(mu[j])**2
        print(f"  {k:>3}  {j:>8}  {j/mod:>10.7f}  {emp:>14.6e}  {val_pop:>12.6e}")

    # ============================================================
    # Same for ξ = 1/2: P_evens vs P_odds. m all odd in tree → all in P_1
    # |μ̂(1/2)|² = (Σ w · (-1))² / Z² = 1 (since all m odd)
    # ============================================================
    print(f"\n# ξ = 1/2: trivial — all m odd, μ̂(1/2) = -1 → |μ̂|² = 1")

    # ξ = 1/6: partition mod 6
    print(f"\n# ξ = 1/6: partition mod 6 (full pop)")
    P6, N6 = pop_partition(6)
    for r in sorted(P6.keys()):
        print(f"  m ≡ {r} mod 6:  count = {N6[r]:>8}    Σw = {P6[r]:>14.0f}")
    # m all odd → m mod 6 ∈ {1, 3, 5}
    # exp(2πi · m / 6): m=1→ζ, m=3→-1, m=5→ζ²=ζ̄ where ζ = exp(πi/3)
    ζ = complex(math.cos(math.pi/3), math.sin(math.pi/3))
    s = P6.get(1, 0) * ζ + P6.get(3, 0) * (-1) + P6.get(5, 0) * ζ.conjugate()
    val_16 = abs(s)**2 / (Z**2)
    print(f"  |μ̂(1/6)|² closed form: {val_16:.6e}")

    print(f"\n# ξ = 1/6 empirical:")
    for k in [10, 12, 14, 16]:
        mu = compute_M_k_and_fft(odd_ints, weights, k)
        target = 1.0 / 6
        mod = 1 << k
        j = round(target * mod)
        if j % 2 == 0:
            j -= 1
        emp = abs(mu[j])**2
        print(f"  k={k:>3}  j={j:>8}  ξ={j/mod:.7f}  |μ̂|²={emp:.6e}    closed form: {val_16:.6e}")

    # Also examine partition counts vs weight-sums (do residues "spread evenly" by count but unevenly by weight?)
    print(f"\n# Mass concentration: P_a / Z (mass per residue) vs N_a / N_tree (count per residue)")
    print(f"\n  mod 3:")
    print(f"    {'r':>3}  {'count':>8}  {'frac count':>10}  {'mass':>14}  {'frac mass':>10}")
    n_tree = sum(N3.values())
    for r in sorted(P3.keys()):
        print(f"    {r:>3}  {N3[r]:>8}  {N3[r]/n_tree:>10.4f}  {P3[r]:>14.0f}  {P3[r]/Z:>10.4f}")

    # Save
    out = os.path.join(OUTDIR, "mj_full_partition.csv")
    with open(out, 'w') as f:
        f.write("modulus,residue,count,frac_count,sum_weight,frac_mass\n")
        for q, (P, N_per) in [(3, (P3, N3)), (6, (P6, N6))]:
            n_total = sum(N_per.values())
            for r in sorted(P.keys()):
                f.write(f"{q},{r},{N_per[r]},{N_per[r]/n_total:.6f},{P[r]:.0f},{P[r]/Z:.6f}\n")
    print(f"\n[save] {out}")
    print(f"\nTotal time: {time.perf_counter()-t0:.2f}s")


if __name__ == "__main__":
    main()
