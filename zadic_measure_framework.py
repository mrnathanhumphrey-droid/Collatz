"""
zadic_measure_framework.py — Z_2 measure framework tests beyond Result 58.

Genuinely new tests:
1. Variant (b): branching-density weighted (number of predecessors per node)
2. Hausdorff dim estimate via box counting at modulus 2^k
3. σ-band conditional: partition tree by orbit-σ band, compute variant (a) per band
4. Conformality check: ratio μ(cylinder)/(2^-k)^δ across cylinders

Result 58 already established:
- variant (a) subtree-size: Pearson +0.86 with D_emp at t=90
- variant (c) 1/σ_orbit: Pearson +0.66
- variant (e) σ-weighted: Pearson +0.78
- M_closed eigvec: Pearson −0.004 (rejected)
- Pearson stable across N=2^16 to 2^22
"""
import csv
import math
import os
import sys
import time
from collections import deque, defaultdict

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


# ============================================================
# Reference data
# ============================================================
def load_empirical_D(t=90):
    with open(r"C:\Collatz\experiments_output\chang_qsd_test.csv") as f:
        for row in csv.DictReader(f):
            if int(row['t']) == t:
                return {int(k[3:]): float(v) for k, v in row.items() if k.startswith('D_r')}


def pearson(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    den = math.sqrt(sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys))
    return num/den if den > 0 else 0.0


def mae(xs, ys):
    return sum(abs(x-y) for x,y in zip(xs,ys)) / len(xs)


# ============================================================
# Inverse tree builder (with branching multiplicity tracked)
# ============================================================
def build_inverse_tree(max_value):
    """Returns:
      tree: dict m -> {parent, depth, n_branches}
        n_branches[m] = number of predecessors of m within bound
    """
    tree = {1: {'parent': None, 'depth': 0, 'n_branches': 0}}
    q = deque([1])
    while q:
        m = q.popleft()
        d = tree[m]['depth']
        if m % 3 == 0:
            continue
        v_start = 2 if (m % 3 == 1) else 1
        n_pred = 0
        for v in range(v_start, 64, 2):
            num = m * (1 << v) - 1
            if num <= 0:
                continue
            pred = num // 3
            if pred > max_value:
                break
            if pred & 1 == 0 or pred == m:
                continue
            n_pred += 1
            if pred not in tree:
                tree[pred] = {'parent': m, 'depth': d + 1, 'n_branches': 0}
                q.append(pred)
        tree[m]['n_branches'] = n_pred
    return tree


def all_subtree_metrics(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    by_depth_desc = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    sigma_sum = {m: float(tree[m]['depth']) for m in tree}
    branching_weighted = {m: float(tree[m]['n_branches']) for m in tree}
    for m in by_depth_desc:
        for c in children[m]:
            size[m] += size[c]
            sigma_sum[m] += sigma_sum[c]
            branching_weighted[m] += branching_weighted[c]
    return size, sigma_sum, branching_weighted, children


def predicted_D_mod32(weights):
    by_r = defaultdict(float); tot = 0.0
    for m, w in weights.items():
        if m & 1:
            by_r[m % 32] += w
            tot += w
    return {r: by_r[r] / tot * 16 for r in range(1, 32, 2)}


# ============================================================
# Step 1: variant (b) branching-density
# ============================================================
def step1_variant_b(tree, sizes, branching_weighted, D_emp):
    print("\n# Step 1: variant (b) branching-density weighted")
    odd_residues = sorted(D_emp.keys())
    emp = [D_emp[r] for r in odd_residues]

    # Definition: w(m) = sum over descendants of n_branches(d)
    # (already computed in branching_weighted accumulator)
    D_b = predicted_D_mod32(branching_weighted)
    b_xs = [D_b[r] for r in odd_residues]
    print(f"  Pearson(b, emp): {pearson(b_xs, emp):+.4f}")
    print(f"  MAE: {mae(b_xs, emp):.4f}")

    # Also try: w(m) = n_branches(m) directly (no descendant accumulation)
    w_direct = {m: float(tree[m]['n_branches']) for m in tree}
    D_b_direct = predicted_D_mod32(w_direct)
    bd_xs = [D_b_direct[r] for r in odd_residues]
    print(f"  Pearson(branching-only direct, emp): {pearson(bd_xs, emp):+.4f}")
    print(f"  MAE direct: {mae(bd_xs, emp):.4f}")

    # Variant (b'): subtree-size × branching-density (combined)
    w_combined = {m: sizes[m] * (1 + tree[m]['n_branches']) for m in tree}
    D_combined = predicted_D_mod32(w_combined)
    bc_xs = [D_combined[r] for r in odd_residues]
    print(f"  Pearson(size × branching, emp): {pearson(bc_xs, emp):+.4f}")
    print(f"  MAE combined: {mae(bc_xs, emp):.4f}")

    return D_b, D_b_direct, D_combined


# ============================================================
# Step 2: Hausdorff dim of measure support via box counting
# ============================================================
def step2_hausdorff(tree, sizes):
    print("\n# Step 2: Hausdorff dim via box counting")
    # For each modulus 2^k for k=5..15:
    # count distinct cylinders of mass > 0
    # also: compute "mass dimension" dim_μ via local dim at scale 2^-k
    # dim_μ(x) ≈ log(μ(B_k(x))) / log(2^-k)

    odd_integers = [m for m in tree if m & 1]
    total_mass = sum(sizes[m] for m in odd_integers)

    print(f"  {'k':>3}  {'2^k':>8}  {'#cyls':>8}  {'box-dim':>8}  {'mass-dim':>10}  {'avg-loc-dim':>12}")
    rows = []
    for k in range(5, 16):
        mod = 1 << k
        cyls = defaultdict(float)
        for m in odd_integers:
            cyls[m % mod] += sizes[m]

        n_active = len(cyls)
        # Box-counting dim: log(n_active) / log(2^k)
        box_dim = math.log(max(n_active, 1)) / (k * math.log(2))

        # Mass dim (Renyi q=0): same as box dim
        # Mass dim (q=2): -log(Σ p_i^2) / log(2^k) where p_i = mass/total
        sum_p_sq = sum((c / total_mass) ** 2 for c in cyls.values())
        mass_dim_q2 = -math.log(sum_p_sq) / (k * math.log(2))

        # Average local dim: avg over points of -log(p_i)/log(2^-k)
        avg_loc_dim = 0.0
        for c in cyls.values():
            p = c / total_mass
            if p > 0:
                avg_loc_dim += p * (-math.log(p)) / (k * math.log(2))

        print(f"  {k:>3}  {mod:>8}  {n_active:>8}  {box_dim:>8.4f}  {mass_dim_q2:>10.4f}  {avg_loc_dim:>12.4f}")
        rows.append((k, mod, n_active, box_dim, mass_dim_q2, avg_loc_dim))

    # Save
    out = os.path.join(OUTDIR, "zadic_hausdorff.csv")
    with open(out, 'w') as f:
        f.write("k,modulus,n_active_cyls,box_dim,mass_dim_q2,avg_local_dim\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.6f},{r[4]:.6f},{r[5]:.6f}\n")
    print(f"  [save] {out}")
    return rows


# ============================================================
# Step 3: σ-band conditional analog
# ============================================================
def step3_sigma_band_predicted(tree, sizes, D_emp):
    print("\n# Step 3: σ-band conditional D_predicted")
    # For each tree node m with depth d = σ_orbit(m), partition by quartile of σ
    odd_integers = sorted([m for m in tree if m & 1])
    depths = sorted([tree[m]['depth'] for m in odd_integers])
    n = len(odd_integers)
    quartiles = [
        depths[int(n * 0.25)],
        depths[int(n * 0.50)],
        depths[int(n * 0.75)],
    ]
    print(f"  σ quartile cutoffs: {quartiles}")
    print(f"  σ range: [{min(depths)}, {max(depths)}]")

    bands = {
        'q1 (0-25%)': (depths[0], quartiles[0]),
        'q2 (25-50%)': (quartiles[0], quartiles[1]),
        'q3 (50-75%)': (quartiles[1], quartiles[2]),
        'q4 (75-100%)': (quartiles[2], depths[-1]),
    }

    odd_residues = sorted(D_emp.keys())
    emp = [D_emp[r] for r in odd_residues]

    print(f"\n  {'band':>15}  {'n_nodes':>8}  {'D_pred (a) Pearson vs emp':>30}")
    band_predictions = {}
    for label, (lo, hi) in bands.items():
        # Restrict to integers with depth in [lo, hi]
        sub = {m: sizes[m] for m in odd_integers if lo <= tree[m]['depth'] <= hi}
        if not sub:
            continue
        D_band = predicted_D_mod32(sub)
        b_xs = [D_band[r] for r in odd_residues]
        r = pearson(b_xs, emp)
        print(f"  {label:>15}  {len(sub):>8}  {r:>+30.4f}")
        band_predictions[label] = D_band

    return band_predictions, bands


# ============================================================
# Step 4: σ-band empirical from chang_qsd_test
# ============================================================
def step4_sigma_band_empirical():
    """Cross-reference D_t(r) at multiple t values from chang_qsd_test.csv.
    Each t corresponds to "alive at iteration t" — these are σ-band proxies in time.
    Higher t = surviving longer = larger σ band."""
    print("\n# Step 4: σ-band empirical (using D_t at multiple t as σ proxies)")
    by_t = {}
    with open(r"C:\Collatz\experiments_output\chang_qsd_test.csv") as f:
        for row in csv.DictReader(f):
            t = int(row['t'])
            by_t[t] = {int(k[3:]): float(v) for k, v in row.items() if k.startswith('D_r')}
    return by_t


# ============================================================
# Step 5: Conformality check
# ============================================================
def step5_conformality(tree, sizes):
    """Test μ(T(A)) ∝ Σ|T'|^δ over a sample of cylinders.
    For 2-adic Z_2, T = Syracuse, so |T'|_2 = |3·m+1|_2 / |m|_2 / 2^v_2(3m+1)
    For odd m: |m|_2 = 1. |3m+1|_2 = 1/2^v_2(3m+1).
    So |T'|_2 = 2^{-2 v_2(3m+1)} (one factor for division by 2^v in the formula, another for 3m+1's 2-adic absolute value contribution? Need care).

    Actually for the OUTER inverse tree, we want how mass distributes under inverse map.
    Conformality of μ on an attractor of T means μ(T(A)) = ∫_A φ(x)^δ dμ for some "conformal Jacobian" φ.

    Practically test: for each cylinder mod 2^k, compute mass μ; for its image under one Syracuse step,
    compute mass; verify a power-law relationship across cylinders.
    """
    print("\n# Step 5: Conformality check (mass scaling under one Syracuse step)")
    k = 6
    mod = 1 << k
    odd_integers = [m for m in tree if m & 1]
    total = sum(sizes[m] for m in odd_integers)

    # Mass on each cylinder mod 64 (32 odd cylinders)
    cyl_mass = defaultdict(float)
    cyl_count = defaultdict(int)
    for m in odd_integers:
        r = m % mod
        cyl_mass[r] += sizes[m]
        cyl_count[r] += 1

    # For each cylinder r mod 64, compute Syracuse image residue (r' = Syracuse(r) mod 64).
    # Most r's are deterministic at mod 64; r=21 fans out.
    print(f"  {'r':>3}  {'mass(r)':>10}  {'count':>8}  {'Syracuse(r) mod 64':>18}  {'mass(image)':>12}  {'ratio':>8}")
    rows_to_print = sorted([r for r in cyl_mass if r in (5, 21, 7, 27, 31, 53, 37, 13)])
    for r in rows_to_print:
        if r % 2 == 0:
            continue
        # Sample lifts to find image
        # We'll just use the smallest valid lift
        m_sample = r
        threem = 3 * m_sample + 1
        while threem & 1 == 0:
            threem >>= 1
        image_r = threem % mod
        ratio = cyl_mass.get(image_r, 0) / cyl_mass[r] if cyl_mass[r] > 0 else 0
        print(f"  {r:>3}  {cyl_mass[r]:>10.0f}  {cyl_count[r]:>8}  {image_r:>18}  {cyl_mass.get(image_r, 0):>12.0f}  {ratio:>8.4f}")


# ============================================================
# Main
# ============================================================
def main():
    D_emp = load_empirical_D(t=90)

    N = 1 << 22
    print(f"# Building inverse tree at N=2^22 = {N}")
    t0 = time.perf_counter()
    tree = build_inverse_tree(N)
    sizes, sigma_sum, branching_w, children = all_subtree_metrics(tree)
    t1 = time.perf_counter()
    print(f"# Tree: {len(tree)} odd nodes, time {t1-t0:.2f}s")

    # Step 1
    D_b, D_b_direct, D_combined = step1_variant_b(tree, sizes, branching_w, D_emp)

    # Step 2
    haus_rows = step2_hausdorff(tree, sizes)

    # Step 3
    band_preds, bands = step3_sigma_band_predicted(tree, sizes, D_emp)

    # Step 4
    by_t = step4_sigma_band_empirical()

    # Compare predicted bands (from inverse tree depth) to empirical bands (from time slices)
    print("\n# Step 4b: predicted band vs empirical D_t — Pearson per pairing")
    odd_residues = sorted(D_emp.keys())
    print(f"  Pred band ↓ vs Empirical t →")
    print(f"  {'pred_band':>15}  {'t=10':>8}  {'t=30':>8}  {'t=50':>8}  {'t=70':>8}  {'t=90':>8}  {'t=110':>8}")
    for label, D_band in band_preds.items():
        b_xs = [D_band[r] for r in odd_residues]
        row = [label]
        for t in [10, 30, 50, 70, 90, 110]:
            if t in by_t:
                e_xs = [by_t[t][r] for r in odd_residues]
                row.append(f"{pearson(b_xs, e_xs):>+8.3f}")
            else:
                row.append(f"{'--':>8}")
        print(f"  {row[0]:>15}  {row[1]}  {row[2]}  {row[3]}  {row[4]}  {row[5]}  {row[6]}")

    # Step 5
    step5_conformality(tree, sizes)

    # Save predictions CSV
    out_csv = os.path.join(OUTDIR, "zadic_measure_predictions.csv")
    with open(out_csv, 'w') as f:
        cols = ['r', 'D_b_subtree_branching', 'D_b_branching_only', 'D_combined_size_x_branching']
        for label in band_preds:
            cols.append(f"D_pred_band_{label.split()[0]}")
        cols.append('D_emp_t90')
        f.write(",".join(cols) + "\n")
        for r in sorted(D_emp.keys()):
            row = [str(r), f"{D_b[r]:.6f}", f"{D_b_direct[r]:.6f}", f"{D_combined[r]:.6f}"]
            for label, D_band in band_preds.items():
                row.append(f"{D_band[r]:.6f}")
            row.append(f"{D_emp[r]:.6f}")
            f.write(",".join(row) + "\n")
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
