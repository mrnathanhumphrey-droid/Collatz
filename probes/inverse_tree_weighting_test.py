"""
inverse_tree_weighting_test.py — test whether inverse tree weighting reproduces D_empirical(r).

Step 1 quick decisive test: does Result 23's leading eigvec of M_closed (k=5)
match D_empirical at t=90 from chang_qsd_test.csv?

If r ~ -0.20 (matches Result 23's prior comparison): outcome (β) REJECTED, then
we still test variants (a), (b), (c), (d) by building inverse tree explicitly.
"""
import csv
import math
import os
import sys
import time
from collections import deque, defaultdict

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# Load reference data
# ============================================================

def load_eigvec_mod32():
    """Load Result 23's eigvec of M_closed at k=5, restricted to 16 odd residues, normalized to mean 1."""
    path = r"C:\Collatz\inverse_tree\inverse_tree_eigvec_mod32.csv"
    eig = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = int(row['residue_mod_32'])
            if r % 2 == 1:
                eig[r] = float(row['predicted_density'])
    # Renormalize to mean = 1 over odd residues
    mean = sum(eig.values()) / len(eig)
    return {r: v / mean for r, v in eig.items()}


def load_empirical_D(t=90):
    """Load D_t(r) at given t from chang_qsd_test.csv. Returns dict r -> D."""
    path = r"C:\Collatz\experiments_output\chang_qsd_test.csv"
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row['t']) == t:
                D = {}
                for k, v in row.items():
                    if k.startswith('D_r'):
                        D[int(k[3:])] = float(v)
                return D
    return None


# ============================================================
# Statistics helpers
# ============================================================

def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (den_x * den_y) if den_x * den_y > 0 else 0.0


def spearman_rank(values):
    sorted_idx = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0] * len(values)
    for rank, i in enumerate(sorted_idx):
        ranks[i] = rank
    return ranks


def spearman(xs, ys):
    return pearson(spearman_rank(xs), spearman_rank(ys))


def mae(xs, ys):
    return sum(abs(x - y) for x, y in zip(xs, ys)) / len(xs)


# ============================================================
# Inverse Collatz tree builder
# ============================================================

def build_inverse_tree_from_one(max_value):
    """BFS-build the inverse Syracuse tree from m=1 up to max_value.

    Returns dict m -> {parent, depth, syracuse_v_used (for stepping back to parent)}
    Only odd integers.
    """
    tree = {1: {'parent': None, 'depth': 0, 'v': None}}
    q = deque([1])
    while q:
        m = q.popleft()
        d = tree[m]['depth']
        # Predecessors: pred = (m * 2^v - 1) / 3, must be positive odd integer ≤ max_value
        # 2^v · m ≡ 1 mod 3. Since 2 ≡ -1 mod 3, 2^v alternates 1, 2, 1, 2 (for v=0,1,2,3)
        # m ≡ 0 mod 3: never works (no predecessor of 1 from m ≡ 0)
        # Actually: 2^v · m ≡ 1 mod 3. If m ≡ 1 mod 3: need 2^v ≡ 1 → v even.
        # If m ≡ 2 mod 3: need 2^v ≡ 2 → v odd.
        # If m ≡ 0 mod 3: no solution.
        if m % 3 == 0:
            continue
        v_start = 2 if (m % 3 == 1) else 1
        for v in range(v_start, 64, 2):  # step by 2
            num = m * (1 << v) - 1
            if num <= 0:
                continue
            if num % 3 != 0:
                # shouldn't happen if v_start logic correct
                continue
            pred = num // 3
            if pred > max_value:
                break  # further v will only increase pred
            if pred & 1 == 0:
                continue
            if pred == m:
                continue  # avoid self-loop at m=1, v=2 case
            if pred not in tree:
                tree[pred] = {'parent': m, 'depth': d + 1, 'v': v}
                q.append(pred)
    return tree


def compute_subtree_sizes(tree):
    """For each m in tree, compute size of subtree rooted at m (including m itself)."""
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    # Topo sort: process by descending depth
    by_depth = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    for m in by_depth:
        for c in children[m]:
            size[m] += size[c]
    return size


def compute_orbit_length(m):
    """Number of Syracuse steps from m to 1. (Forward.)"""
    steps = 0
    n = m
    cap = 100000
    while n != 1 and steps < cap:
        if n & 1 == 0:
            n >>= 1
        else:
            n = 3 * n + 1
        steps += 1
    return steps


def compute_orbit_length_odd(m):
    """Steps in the ODD-only Syracuse iteration from m to 1."""
    steps = 0
    n = m
    cap = 100000
    while n != 1 and steps < cap:
        n = 3 * n + 1
        while n & 1 == 0:
            n >>= 1
        steps += 1
    return steps


# ============================================================
# Variants: predicted D(r) from per-integer weights
# ============================================================

def predicted_D_mod32(weights):
    """Marginalize weights to mod 32, restricted to odd residues, mean-1 normalized."""
    by_r = defaultdict(float)
    total = 0.0
    for m, w in weights.items():
        if m & 1:
            by_r[m % 32] += w
            total += w
    if total == 0:
        return {r: 0 for r in range(1, 32, 2)}
    n_residues = 16
    D = {}
    for r in range(1, 32, 2):
        D[r] = by_r[r] / total * n_residues
    return D


# ============================================================
# Main
# ============================================================

def main():
    print("# Step 1 — quick eigvec vs D_empirical comparison")
    eigvec = load_eigvec_mod32()
    D_emp = load_empirical_D(t=90)

    odd_residues = sorted(eigvec.keys())
    print(f"\n{'r':>3}  {'eigvec':>8}  {'empirical':>10}")
    eig_xs = []
    emp_xs = []
    for r in odd_residues:
        e = eigvec[r]
        d = D_emp[r]
        eig_xs.append(e)
        emp_xs.append(d)
        print(f"{r:>3}  {e:>8.4f}  {d:>10.4f}")

    r_pearson = pearson(eig_xs, emp_xs)
    r_spearman = spearman(eig_xs, emp_xs)
    err = mae(eig_xs, emp_xs)
    print(f"\n  Pearson r:     {r_pearson:+.4f}")
    print(f"  Spearman ρ:    {r_spearman:+.4f}")
    print(f"  Mean abs err:  {err:.4f}")

    if r_pearson < 0:
        print("\n  >>> outcome (β) REJECTED: eigvec & D_empirical are anticorrelated")
    elif r_pearson < 0.5:
        print("\n  >>> outcome (β) REJECTED: weak/no correlation")
    elif r_pearson < 0.85:
        print("\n  >>> outcome (γ) candidate: moderate correlation")
    else:
        print("\n  >>> outcome (β) potential: strong correlation")

    # Step 2 — build inverse tree at small N, compute variant (a) subtree-size, compare
    print("\n# Step 2 — inverse tree at N=2^16, variant (a) subtree-size")
    N = 1 << 16
    t0 = time.perf_counter()
    tree = build_inverse_tree_from_one(N)
    t1 = time.perf_counter()
    print(f"  Tree built: {len(tree)} odd nodes ≤ {N}, time {t1-t0:.2f}s")

    sizes = compute_subtree_sizes(tree)
    D_a = predicted_D_mod32(sizes)
    a_xs = [D_a[r] for r in odd_residues]
    r_a_eig = pearson(a_xs, eig_xs)
    r_a_emp = pearson(a_xs, emp_xs)
    err_a_emp = mae(a_xs, emp_xs)
    print(f"\n  Pearson(variant a, eigvec):     {r_a_eig:+.4f}  (sanity check, expect ~+1)")
    print(f"  Pearson(variant a, empirical):  {r_a_emp:+.4f}")
    print(f"  MAE(variant a, empirical):      {err_a_emp:.4f}")

    # Variant (d) — equivalent to (a) at finite N; verify numerically
    # variant(d)(m) = #{n ≤ N : m on path 1...n in inverse tree} = subtree size of m
    # → numerically identical to (a). Confirm.
    D_d = D_a
    d_xs = [D_d[r] for r in odd_residues]
    r_d_emp = pearson(d_xs, emp_xs)
    print(f"\n  variant (d) is numerically identical to (a) — same subtree count")
    print(f"  Pearson(variant d, empirical):  {r_d_emp:+.4f}")

    # Variant (b) — depth-weighted Σ 1/d
    print("\n# Variant (b) — depth-weighted Σ_{descendants} 1/depth(desc, m)")
    t0 = time.perf_counter()
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    # For each node m, compute Σ_{descendants d} 1/(depth(d) - depth(m))
    # BFS from each m. Could be O(N^2). For N=2^16 that's 4e9. Too slow.
    # Instead: post-order DP. Σ_{desc} 1/dist = sum over children c of (1 + Σ_{desc(c)} 1/(dist+1))
    # That recursion isn't decomposable simply. Skip for now and use approx:
    # depth-weighted (b) ≈ subtree size with depth normalization. Skip.
    print("  (deferred — non-decomposable DP, skip variant b for N>=2^14)")

    # Variant (c) — 1/σ_orbit, summed over descendants
    # σ_orbit(m) = forward odd Syracuse steps from m to 1
    # Each integer in tree has σ = depth in inverse tree (since inverse-tree depth = forward-orbit distance)
    print("\n# Variant (c) — Σ_{descendants} 1/σ_orbit(desc) (= 1/depth in inverse tree)")
    weights_c = defaultdict(float)
    # Σ 1/depth over descendants
    # Post-order: weight(m) = 1/d(m) + Σ weight(c) for c in children
    for m in sorted(tree.keys(), key=lambda x: -tree[x]['depth']):
        d_m = tree[m]['depth']
        own = 1.0 / d_m if d_m > 0 else 1.0  # m=1 has depth 0
        weights_c[m] = own + sum(weights_c[c] for c in children.get(m, []))
    D_c = predicted_D_mod32(weights_c)
    c_xs = [D_c[r] for r in odd_residues]
    r_c_emp = pearson(c_xs, emp_xs)
    err_c_emp = mae(c_xs, emp_xs)
    print(f"  Pearson(variant c, empirical):  {r_c_emp:+.4f}")
    print(f"  MAE(variant c, empirical):      {err_c_emp:.4f}")

    # Variant (e) — σ-weighted (NOT inverse-σ): orbits with long σ contribute more
    # This is the "survivor-conditioning" hypothesis: D_t(r) at large t favors long-orbit integers
    print("\n# Variant (e) NEW — σ-weighted Σ_{descendants} σ_orbit(desc)")
    weights_e = defaultdict(float)
    for m in sorted(tree.keys(), key=lambda x: -tree[x]['depth']):
        d_m = tree[m]['depth']
        own = float(d_m)
        weights_e[m] = own + sum(weights_e[c] for c in children.get(m, []))
    D_e = predicted_D_mod32(weights_e)
    e_xs = [D_e[r] for r in odd_residues]
    r_e_emp = pearson(e_xs, emp_xs)
    err_e_emp = mae(e_xs, emp_xs)
    print(f"  Pearson(variant e, empirical):  {r_e_emp:+.4f}")
    print(f"  MAE(variant e, empirical):      {err_e_emp:.4f}")

    # Variant (f) — SELF visit by residue (not descendant-summed): just count integers per residue
    # = uniform measure on odd integers in [1, N] mod 32. Should be ~uniform → all D=1
    counts_f = defaultdict(int)
    for m in tree:
        counts_f[m] = 1
    D_f = predicted_D_mod32(counts_f)
    f_xs = [D_f[r] for r in odd_residues]
    r_f_emp = pearson(f_xs, emp_xs)
    print(f"\n  Variant (f) uniform on tree: Pearson = {r_f_emp:+.4f}  (sanity, expect near 0)")

    # Print full table
    print("\n# Full D-prediction table (residues mod 32, mean-1 normalized)")
    print(f"  {'r':>3}  {'eig':>7}  {'(a)':>7}  {'(c)':>7}  {'(e)':>7}  {'emp':>8}")
    for r in odd_residues:
        print(f"  {r:>3}  {eigvec[r]:>7.3f}  {D_a[r]:>7.3f}  {D_c[r]:>7.3f}  {D_e[r]:>7.3f}  {D_emp[r]:>8.3f}")

    print("\n# Summary table")
    print(f"  {'variant':>10}  {'r vs emp':>10}  {'MAE vs emp':>10}")
    print(f"  {'eigvec':>10}  {r_pearson:>+10.4f}  {err:>10.4f}")
    print(f"  {'(a)':>10}  {r_a_emp:>+10.4f}  {err_a_emp:>10.4f}")
    print(f"  {'(c)':>10}  {r_c_emp:>+10.4f}  {err_c_emp:>10.4f}")
    print(f"  {'(e)':>10}  {r_e_emp:>+10.4f}  {err_e_emp:>10.4f}")

    # Save CSV
    out_csv = os.path.join(OUTDIR, "inverse_tree_predicted_D.csv")
    with open(out_csv, 'w') as f:
        f.write("r,eigvec,variant_a,variant_c,variant_e,empirical\n")
        for r in odd_residues:
            f.write(f"{r},{eigvec[r]:.6f},{D_a[r]:.6f},{D_c[r]:.6f},{D_e[r]:.6f},{D_emp[r]:.6f}\n")
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
