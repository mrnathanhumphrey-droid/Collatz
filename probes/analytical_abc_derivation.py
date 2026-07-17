"""
analytical_abc_derivation.py — Derive (a, b, c) limiting fractions from first principles.

Path-counting derivation:
- Random odd n ∈ [1, N] uniform → P(n ≡ r mod 3) = 1/3 for r ∈ {0, 1, 2}
- Forward Syracuse step T(m) = (3m+1)/2^v with v = v_2(3m+1) ~ Geom(1/2)
  - v even (P = 1/3) → T(m) ≡ 1 mod 3
  - v odd  (P = 2/3) → T(m) ≡ 2 mod 3
- Path n = m_0 → m_1 → ... → m_h = 1 has h+1 nodes.
- m_0 mod 3 ∈ {0, 1, 2} (1/3 each)
- m_h = 1 (always residue 1)
- m_i (1 ≤ i ≤ h-1): residue 1 with prob 1/3, residue 2 with prob 2/3

Average residue counts over random path:
  ⟨#r0⟩ = 1/3
  ⟨#r1⟩ = 1/3 + 1 + (h-1)·1/3 = (h+3)/3
  ⟨#r2⟩ = 1/3 + 0 + (h-1)·2/3 = (2h-1)/3
  Sum = h + 1 = D (path length in nodes) ✓

Mass fractions (D = avg path length in nodes):
  a = (1/3) / D = 1/(3D)
  b = (⟨D⟩+2) / (3D)   ← from substitution h = D-1
  c = (2⟨D⟩-3) / (3D)

Closed form |μ̂(1/3)|²:
  |μ̂(1/3)|² = (a² + b² + c² - ab - bc - ca)
            = ½[(a-b)² + (b-c)² + (a-c)²]
            = (D² - 4D + 7) / (3D²)

Asymptotic D → ∞: → 1/3.
"""
import math
import os
import sys
import time
import csv
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


def subtree_sizes_and_partition(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    by_depth_desc = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    for m in by_depth_desc:
        for c in children[m]:
            size[m] += size[c]
    # Partition by m mod 3
    P = defaultdict(float)
    counts = defaultdict(int)
    for m in tree:
        if m & 1:
            r = m % 3
            P[r] += size[m]
            counts[r] += 1
    Z = sum(P.values())
    return Z, dict(P), dict(counts), size


def closed_form_abc(D):
    """Analytical (a, b, c) given path length D (in nodes)."""
    a = 1.0 / (3 * D)
    b = (D + 2) / (3 * D)
    c = (2 * D - 3) / (3 * D)
    return a, b, c


def closed_form_mu_hat_third_sq(D):
    """Closed form |μ̂(1/3)|² = (D² − 4D + 7) / (3D²)."""
    return (D * D - 4 * D + 7) / (3 * D * D)


def main():
    print("# Analytical (a, b, c) derivation — first principles")
    print("# Random-path heuristic: v_i ~ Geom(1/2), P(v even)=1/3, P(v odd)=2/3")
    print("# Each path mod-3 residues: m_0 uniform {0,1,2}; m_h=1 (residue 1);")
    print("# intermediate m_i ∈ {1, 2} with prob 1/3, 2/3 respectively.")
    print()

    print("# Closed forms (path length D in nodes):")
    print("#   a = 1 / (3D)")
    print("#   b = (D + 2) / (3D)")
    print("#   c = (2D − 3) / (3D)")
    print("#   |μ̂(1/3)|² = (D² − 4D + 7) / (3D²)  →  1/3 as D → ∞")
    print()

    # Verify at multiple N
    Ns = [1<<k for k in [14, 16, 18, 20, 22]]
    rows = []
    print(f"{'N':>10}  {'#nodes':>10}  {'D=Z/N_t':>10}  {'a_emp':>8}  {'a_pred':>8}  {'b_emp':>8}  {'b_pred':>8}  {'c_emp':>8}  {'c_pred':>8}  {'|μ̂(1/3)|² emp':>14}  {'|μ̂|² pred':>10}  {'|μ̂|²(D→∞)':>10}")
    for N in Ns:
        t0 = time.perf_counter()
        tree = build_inverse_tree(N)
        Z, P, cnt, size = subtree_sizes_and_partition(tree)
        N_tree = sum(cnt.values())
        D = Z / N_tree

        # Empirical a, b, c
        a_e = P.get(0, 0) / Z
        b_e = P.get(1, 0) / Z
        c_e = P.get(2, 0) / Z

        # Predicted from closed form
        a_p, b_p, c_p = closed_form_abc(D)

        # Empirical |μ̂(1/3)|² via closed-form partition formula
        emp_sq = 0.5 * ((a_e - b_e)**2 + (b_e - c_e)**2 + (a_e - c_e)**2)

        # Predicted analytical (from D)
        pred_sq = closed_form_mu_hat_third_sq(D)

        # Asymptotic limit
        asymp = 1.0 / 3.0

        rows.append({'N': N, 'n_nodes': N_tree, 'D': D, 'a_emp': a_e, 'a_pred': a_p,
                     'b_emp': b_e, 'b_pred': b_p, 'c_emp': c_e, 'c_pred': c_p,
                     'emp_sq': emp_sq, 'pred_sq': pred_sq, 'asymp': asymp,
                     'time': time.perf_counter() - t0})

        print(f"{N:>10}  {N_tree:>10}  {D:>10.4f}  {a_e:>8.5f}  {a_p:>8.5f}  {b_e:>8.4f}  {b_p:>8.4f}  {c_e:>8.4f}  {c_p:>8.4f}  {emp_sq:>14.5f}  {pred_sq:>10.5f}  {asymp:>10.5f}")

    # Save
    out = os.path.join(OUTDIR, "analytical_abc.csv")
    with open(out, 'w') as f:
        f.write("N,n_nodes,D,a_emp,a_pred,b_emp,b_pred,c_emp,c_pred,mu_hat_third_sq_emp,mu_hat_third_sq_pred,mu_hat_third_sq_asymp,time_s\n")
        for r in rows:
            f.write(f"{r['N']},{r['n_nodes']},{r['D']:.4f},{r['a_emp']:.6f},{r['a_pred']:.6f},"
                    f"{r['b_emp']:.6f},{r['b_pred']:.6f},{r['c_emp']:.6f},{r['c_pred']:.6f},"
                    f"{r['emp_sq']:.6f},{r['pred_sq']:.6f},{r['asymp']:.6f},{r['time']:.2f}\n")
    print(f"\n[save] {out}")

    # Errors and convergence
    print(f"\n# Match quality:")
    print(f"  {'N':>10}  {'Δa abs':>10}  {'Δb abs':>10}  {'Δc abs':>10}  {'Δ|μ̂|² abs':>14}  {'Δ|μ̂|² → 1/3':>12}")
    for r in rows:
        da = abs(r['a_emp'] - r['a_pred'])
        db = abs(r['b_emp'] - r['b_pred'])
        dc = abs(r['c_emp'] - r['c_pred'])
        d_sq = abs(r['emp_sq'] - r['pred_sq'])
        d_asymp = abs(r['emp_sq'] - r['asymp'])
        print(f"  {r['N']:>10}  {da:>10.6f}  {db:>10.6f}  {dc:>10.6f}  {d_sq:>14.6f}  {d_asymp:>12.6f}")

    # Empirical D vs Lagarias-Sinai prediction
    print(f"\n# Lagarias-Sinai heuristic: ⟨h⟩ = ⟨log m⟩ / |log(3/4)| ≈ ⟨log m⟩ / 0.288")
    print(f"# For uniform m on [1, N]: ⟨log m⟩ ≈ log N − 1")
    print(f"# So ⟨D⟩ = ⟨h⟩ + 1 ≈ (log N − 1)/0.288 + 1 ≈ 3.47·log(N) − 2.47")
    print()
    print(f"  {'N':>10}  {'log N':>8}  {'D_pred (LS)':>12}  {'D_empirical':>12}  {'ratio':>8}")
    for r in rows:
        log_N = math.log(r['N'])
        D_LS = log_N / 0.288 - 2.47
        ratio = r['D'] / D_LS if D_LS > 0 else float('nan')
        print(f"  {r['N']:>10}  {log_N:>8.3f}  {D_LS:>12.3f}  {r['D']:>12.3f}  {ratio:>8.4f}")

    # Sub-residue analysis: verify v parity → next residue mapping
    print(f"\n# Verify v-parity → next residue mapping at N=2^22")
    N = 1 << 22
    tree = build_inverse_tree(N)
    # For each non-root, non-leaf-parent, check v parity matches expected residue
    # Reverse: for parent m, child c = pred. T(c) = m. So m mod 3 dictated by v=v_2(3c+1) parity.
    n_check = 0
    n_match = 0
    for c in tree:
        if not (c & 1):
            continue
        info = tree[c]
        if info['parent'] is None:
            continue
        m = info['parent']
        # v = v_2(3c+1)
        threec = 3 * c + 1
        v = 0
        while threec & 1 == 0:
            threec >>= 1
            v += 1
        # Predicted m mod 3: 1 if v even, 2 if v odd
        pred_mod = 1 if (v % 2 == 0) else 2
        actual_mod = m % 3
        n_check += 1
        if pred_mod == actual_mod:
            n_match += 1
        if n_check >= 1000000:
            break
    print(f"  Checked {n_check} (child, parent) pairs: {n_match} match v-parity rule = {100*n_match/n_check:.2f}%")

    # Distribution of v at random Syracuse step (under iid Geom(1/2) heuristic vs empirical)
    print(f"\n# Distribution of v = v_2(3c+1) for random c in tree (vs Geom(1/2) heuristic)")
    v_counts = defaultdict(int)
    for c in tree:
        if c & 1:
            threec = 3 * c + 1
            v = 0
            while threec & 1 == 0:
                threec >>= 1
                v += 1
            if v >= 1:
                v_counts[v] += 1
    total = sum(v_counts.values())
    print(f"  {'v':>3}  {'empirical':>12}  {'Geom(1/2)':>12}")
    for v in sorted(v_counts.keys())[:10]:
        emp_p = v_counts[v] / total
        geom_p = 0.5 ** v
        print(f"  {v:>3}  {emp_p:>12.6f}  {geom_p:>12.6f}")

    # P(v even) and P(v odd)
    pe = sum(v_counts[v] for v in v_counts if v % 2 == 0) / total
    po = sum(v_counts[v] for v in v_counts if v % 2 == 1) / total
    print(f"\n  P(v even) empirical: {pe:.6f}    heuristic: 1/3 = {1/3:.6f}")
    print(f"  P(v odd)  empirical: {po:.6f}    heuristic: 2/3 = {2/3:.6f}")


if __name__ == "__main__":
    main()
