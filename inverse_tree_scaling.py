"""inverse_tree_scaling.py — does variant (a) Pearson r vs D_emp survive at larger N?

If Pearson stays >0.85 as N → 2^20 / 2^22, outcome (α) potential.
If Pearson drops toward 0 (matching eigvec's r=-0.004), outcome (γ) finite-N artifact.
"""
import csv
import math
import os
import sys
import time
from collections import deque, defaultdict

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def load_eigvec_mod32():
    eig = {}
    with open(r"C:\Collatz\inverse_tree\inverse_tree_eigvec_mod32.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = int(row['residue_mod_32'])
            if r % 2 == 1:
                eig[r] = float(row['predicted_density'])
    m = sum(eig.values()) / len(eig)
    return {r: v / m for r, v in eig.items()}


def load_empirical_D(t=90):
    with open(r"C:\Collatz\experiments_output\chang_qsd_test.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row['t']) == t:
                return {int(k[3:]): float(v) for k, v in row.items() if k.startswith('D_r')}


def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n; my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return num / den if den > 0 else 0.0


def mae(xs, ys):
    return sum(abs(x - y) for x, y in zip(xs, ys)) / len(xs)


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


def subtree_sizes_and_metrics(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    by_depth_desc = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    sigma_sum = {m: float(tree[m]['depth']) for m in tree}  # Σ depths over subtree
    inv_sigma_sum = {m: (1.0 / max(1, tree[m]['depth'])) for m in tree}
    for m in by_depth_desc:
        for c in children[m]:
            size[m] += size[c]
            sigma_sum[m] += sigma_sum[c]
            inv_sigma_sum[m] += inv_sigma_sum[c]
    return size, sigma_sum, inv_sigma_sum


def predicted_D_mod32(weights):
    by_r = defaultdict(float)
    total = 0.0
    for m, w in weights.items():
        if m & 1:
            by_r[m % 32] += w
            total += w
    return {r: by_r[r] / total * 16 for r in range(1, 32, 2)}


def main():
    eig = load_eigvec_mod32()
    D_emp = load_empirical_D(t=90)
    odd_residues = sorted(eig.keys())
    eig_xs = [eig[r] for r in odd_residues]
    emp_xs = [D_emp[r] for r in odd_residues]

    Ns = [1 << k for k in (16, 18, 20)]
    results = []
    print(f"{'N':>10} {'#nodes':>10} {'time':>7}  {'(a) r':>8}  {'(a) MAE':>8}  {'(a) v eig':>10}  {'(c) r':>8}  {'(e) r':>8}")
    for N in Ns:
        t0 = time.perf_counter()
        tree = build_inverse_tree(N)
        size, sigma_sum, inv_sig = subtree_sizes_and_metrics(tree)
        t1 = time.perf_counter()

        D_a = predicted_D_mod32(size)
        D_c = predicted_D_mod32(inv_sig)
        D_e = predicted_D_mod32(sigma_sum)

        a_xs = [D_a[r] for r in odd_residues]
        c_xs = [D_c[r] for r in odd_residues]
        e_xs = [D_e[r] for r in odd_residues]

        ra = pearson(a_xs, emp_xs)
        rc = pearson(c_xs, emp_xs)
        re = pearson(e_xs, emp_xs)
        ma = mae(a_xs, emp_xs)
        ra_eig = pearson(a_xs, eig_xs)

        print(f"{N:>10}  {len(tree):>9}  {t1-t0:>6.2f}s  {ra:>+8.4f}  {ma:>8.4f}  {ra_eig:>+10.4f}  {rc:>+8.4f}  {re:>+8.4f}")
        results.append({'N': N, 'n_nodes': len(tree), 'time_s': t1-t0,
                        'D_a': D_a, 'D_c': D_c, 'D_e': D_e,
                        'r_a_emp': ra, 'r_c_emp': rc, 'r_e_emp': re,
                        'mae_a_emp': ma, 'r_a_eig': ra_eig})

    # Save full prediction table at largest N
    largest = results[-1]
    out_csv = os.path.join(OUTDIR, "inverse_tree_scaling.csv")
    with open(out_csv, 'w') as f:
        f.write("N,n_nodes,time_s,r_a_emp,mae_a_emp,r_c_emp,r_e_emp,r_a_eig\n")
        for r in results:
            f.write(f"{r['N']},{r['n_nodes']},{r['time_s']:.2f},"
                    f"{r['r_a_emp']:.6f},{r['mae_a_emp']:.6f},{r['r_c_emp']:.6f},"
                    f"{r['r_e_emp']:.6f},{r['r_a_eig']:.6f}\n")
    print(f"\n[save] {out_csv}")

    # Convergence trend
    if len(results) >= 2:
        trend_a = [r['r_a_emp'] for r in results]
        trend_e = [r['r_e_emp'] for r in results]
        print(f"\n# Trend Pearson(variant a, emp): {trend_a}")
        print(f"# Trend Pearson(variant e, emp): {trend_e}")
        if abs(trend_a[-1] - trend_a[0]) < 0.05:
            print("# variant (a) STABLE across N (does NOT decay to eigvec's near-zero)")
        elif trend_a[-1] < trend_a[0] - 0.1:
            print("# variant (a) DECAYING toward eigvec asymptote → outcome (γ) finite-N artifact")
        else:
            print("# variant (a) INCREASING — outcome (α) potential")

    # Print final D_a table at largest N alongside D_emp
    print("\n# At largest N:")
    print(f"  {'r':>3}  {'(a)':>7}  {'(e)':>7}  {'emp':>7}  {'(a)-emp':>8}")
    for r in odd_residues:
        a = largest['D_a'][r]; e = largest['D_e'][r]; em = D_emp[r]
        print(f"  {r:>3}  {a:>7.3f}  {e:>7.3f}  {em:>7.3f}  {a-em:>+8.3f}")


if __name__ == "__main__":
    main()
