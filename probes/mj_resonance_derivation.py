"""
mj_resonance_derivation.py — atomic decomposition of the trajectory measure
on Z_2 at the {m_j = (4^j-1)/3} chain. Derives closed-form Fourier
coefficients at resonance frequencies ξ = 1/3, 1/2, 1/6, 2/3, 5/6 and
compares to empirical values from Result 62.

Result 62: empirical |μ̂(j/2^k)|² ≈ 0.034 at the 1/3-resonance peak (j=341, k=10).

Atomic decomposition:
  μ̂_atomic(ξ) = (1/Z) · Σ_j w_j · exp(2πi · ξ · m_j)

For ξ = 1/3, since m_j mod 3 cycles through {2, 0, 1, 2, 0, 1, ...} starting at j=2:

  μ̂_atomic(1/3) = (1/Z) · [S_0 + S_1·ω + S_2·ω²]

where ω = exp(2πi/3), S_a = sum of w_j over j with m_j ≡ a mod 3.

Magnitude squared identity (since ω + ω² = -1, ω·ω̄ = 1):

  |a + b·ω + c·ω²|² = a² + b² + c² - ab - bc - ca
                    = ½[(a-b)² + (b-c)² + (c-a)²]
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
    print("# m_j atomic decomposition derivation")
    N = 1 << 22
    t0 = time.perf_counter()
    tree = build_inverse_tree(N)
    sizes = subtree_sizes(tree)
    odd_ints = np.array([m for m in tree if m & 1], dtype=np.int64)
    weights = np.array([float(sizes[m]) for m in odd_ints])
    Z = weights.sum()
    print(f"# Tree built: {len(odd_ints)} odd nodes, Z = {Z:.4e}, time {time.perf_counter()-t0:.2f}s")

    # ============================================================
    # Step 1: Identify m_j and extract w_j
    # ============================================================
    print(f"\n# Step 1: m_j chain and atomic weights")
    m_j_vals = []
    j = 1
    while True:
        m = (4**j - 1) // 3
        if m > N:
            break
        m_j_vals.append((j, m))
        j += 1

    # Get w_j from sizes dict (m_j must be in tree as predecessor of 1)
    print(f"  {'j':>3}  {'m_j':>10}  {'m_j mod 3':>10}  {'subtree_size':>14}  {'w_j/Z':>10}")
    w_j = {}
    for j, m in m_j_vals:
        if m in sizes:
            w_j[j] = sizes[m]
            print(f"  {j:>3}  {m:>10}  {m % 3:>10}  {sizes[m]:>14}  {sizes[m]/Z:>10.6e}")

    # ============================================================
    # Step 2: Test geometric decay w_j ~ c · r^j
    # ============================================================
    print(f"\n# Step 2: Geometric decay test")
    js = sorted(w_j.keys())
    log_w = [math.log(w_j[j]) for j in js]
    if len(js) >= 3:
        slope, intercept = np.polyfit(js, log_w, 1)
        r_geom = math.exp(slope)
        c_geom = math.exp(intercept)
        ssr = sum((log_w[i] - (slope * js[i] + intercept)) ** 2 for i in range(len(js)))
        print(f"  log w_j = {slope:+.4f} · j + {intercept:.4f}")
        print(f"  → r = exp(slope) = {r_geom:.6f}    c = exp(intercept) = {c_geom:.4e}")
        print(f"  → w_j ≈ {c_geom:.4e} × {r_geom:.4f}^j")
        print(f"  Residual SS: {ssr:.4f} (low if pure geometric)")
        print()
        print(f"  {'j':>3}  {'empirical w_j':>14}  {'fit c·r^j':>14}  {'ratio':>8}")
        for i, j in enumerate(js):
            fit = c_geom * (r_geom ** j)
            print(f"  {j:>3}  {w_j[j]:>14}  {fit:>14.4e}  {w_j[j]/fit:>8.4f}")

    # ============================================================
    # Step 3: Closed-form |μ̂_atomic(1/3)|²
    # ============================================================
    print(f"\n# Step 3: Closed-form |μ̂_atomic(ξ)|² for various resonance ξ")

    def atomic_fourier(xi):
        """μ̂_atomic(ξ) = (1/Z) Σ w_j exp(2πi ξ m_j)"""
        s = 0.0 + 0j
        for j, m in m_j_vals:
            if j in w_j:
                s += w_j[j] * complex(math.cos(2 * math.pi * xi * m),
                                      math.sin(2 * math.pi * xi * m))
        return s / Z

    def partition_sum(divisor):
        """Group w_j by m_j mod divisor → return dict residue → sum w_j."""
        groups = defaultdict(float)
        for j, m in m_j_vals:
            if j in w_j:
                groups[m % divisor] += w_j[j]
        return dict(groups)

    # ξ = 1/3
    print(f"\n  ξ = 1/3:")
    P3 = partition_sum(3)
    print(f"    partition sums (m_j mod 3):  S_0={P3.get(0,0):.0f}  S_1={P3.get(1,0):.0f}  S_2={P3.get(2,0):.0f}")
    a, b, c = P3.get(0, 0.0), P3.get(1, 0.0), P3.get(2, 0.0)
    abs_sq_closed = (a*a + b*b + c*c - a*b - b*c - c*a) / (Z * Z)
    abs_sq_alt = 0.5 * ((a-b)**2 + (b-c)**2 + (c-a)**2) / (Z * Z)
    direct = abs(atomic_fourier(1.0/3))**2
    print(f"    closed form (a²+b²+c²-ab-bc-ca)/Z²: {abs_sq_closed:.6e}")
    print(f"    closed form ½[(a-b)²+(b-c)²+(c-a)²]/Z²: {abs_sq_alt:.6e}")
    print(f"    direct |Σ w_j ω^{{(m_j mod 3)}}|²/Z²: {direct:.6e}")

    # ξ = 1/2
    print(f"\n  ξ = 1/2:")
    # all m_j odd → exp(πi · m_j) = -1
    # μ̂_atomic(1/2) = -Σ w_j / Z
    sum_w = sum(w_j.values())
    abs_sq_half = (sum_w / Z) ** 2
    direct = abs(atomic_fourier(0.5))**2
    print(f"    Σ w_j = {sum_w:.0f}  Σw_j/Z = {sum_w/Z:.6e}")
    print(f"    closed form (Σw_j/Z)²: {abs_sq_half:.6e}")
    print(f"    direct: {direct:.6e}")

    # ξ = 1/6
    print(f"\n  ξ = 1/6:")
    P6 = partition_sum(6)
    # m_j is odd, so m_j mod 6 ∈ {1, 3, 5}
    print(f"    partition (m_j mod 6): {dict(sorted(P6.items()))}")
    # exp(2πi·m/6): m=1 → ζ; m=3 → -1; m=5 → ζ̄  where ζ=exp(πi/3)
    ζ = complex(math.cos(math.pi/3), math.sin(math.pi/3))
    s = P6.get(1, 0) * ζ + P6.get(3, 0) * (-1) + P6.get(5, 0) * ζ.conjugate()
    abs_sq_16 = abs(s)**2 / (Z**2)
    direct = abs(atomic_fourier(1.0/6))**2
    print(f"    closed form |P_1·ζ + P_3·(-1) + P_5·ζ̄|²/Z²: {abs_sq_16:.6e}")
    print(f"    direct: {direct:.6e}")

    # ξ = 2/3 (conjugate of 1/3, must give same magnitude)
    print(f"\n  ξ = 2/3:")
    direct = abs(atomic_fourier(2.0/3))**2
    print(f"    direct |μ̂_atomic(2/3)|²: {direct:.6e}    (should equal |μ̂_atomic(1/3)|² = {abs_sq_closed:.6e})")

    # ξ = 5/6 (conjugate of 1/6)
    print(f"\n  ξ = 5/6:")
    direct = abs(atomic_fourier(5.0/6))**2
    print(f"    direct |μ̂_atomic(5/6)|²: {direct:.6e}    (should equal |μ̂_atomic(1/6)|² = {abs_sq_16:.6e})")

    # ============================================================
    # Step 4: Compare to empirical |μ̂(j/2^k)|² at j/2^k → resonance ξ
    # ============================================================
    print(f"\n# Step 4: Empirical comparison via FFT at dyadic approximations")
    K_VALS = [8, 10, 12, 14, 16]
    mu_hat_per_k = {k: compute_M_k_and_fft(odd_ints, weights, k) for k in K_VALS}

    def closest_dyadic(target, k):
        """j s.t. j/2^k closest to target (j odd preferred)."""
        mod = 1 << k
        j_continuous = target * mod
        candidates = [int(math.floor(j_continuous)), int(math.ceil(j_continuous))]
        return min(candidates, key=lambda j: abs(j/mod - target))

    print(f"  {'ξ':>6}  {'k':>3}  {'j':>8}  {'j/2^k':>10}  {'|μ̂_emp|²':>12}  {'|μ̂_atomic|² (closed)':>22}")

    targets = [
        (1.0/3, '1/3', abs_sq_closed),
        (1.0/2, '1/2', abs_sq_half),
        (1.0/6, '1/6', abs_sq_16),
        (2.0/3, '2/3', abs_sq_closed),  # symmetric
        (5.0/6, '5/6', abs_sq_16),
    ]
    rows = []
    for target, label, atomic_pred in targets:
        for k in K_VALS:
            j = closest_dyadic(target, k)
            if 0 <= j < (1 << k):
                emp = abs(mu_hat_per_k[k][j])**2
                rows.append((label, target, k, j, j/(1<<k), emp, atomic_pred))
                print(f"  {label:>6}  {k:>3}  {j:>8}  {j/(1<<k):>10.6f}  {emp:>12.6e}  {atomic_pred:>22.6e}")

    # Save
    out = os.path.join(OUTDIR, "mj_resonance_predictions.csv")
    with open(out, 'w') as f:
        f.write("xi_label,xi_value,k,j,xi_dyadic,mu_hat_sq_empirical,mu_hat_sq_atomic_pred\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]:.6f},{r[2]},{r[3]},{r[4]:.6f},{r[5]:.6e},{r[6]:.6e}\n")
    print(f"\n[save] {out}")

    # Save w_j weights
    out_w = os.path.join(OUTDIR, "mj_atomic_weights.csv")
    with open(out_w, 'w') as f:
        f.write("j,m_j,m_j_mod_3,m_j_mod_6,subtree_size,w_over_Z\n")
        for j, m in m_j_vals:
            if j in w_j:
                f.write(f"{j},{m},{m%3},{m%6},{w_j[j]},{w_j[j]/Z:.6e}\n")
    print(f"[save] {out_w}")

    # ============================================================
    # Step 5: Decay fit summary
    # ============================================================
    if len(js) >= 3:
        out_fit = os.path.join(OUTDIR, "mj_decay_fit.csv")
        with open(out_fit, 'w') as f:
            f.write("param,value\n")
            f.write(f"slope,{slope:.6f}\n")
            f.write(f"intercept,{intercept:.6f}\n")
            f.write(f"r_geom,{r_geom:.6f}\n")
            f.write(f"c_geom,{c_geom:.6e}\n")
            f.write(f"ssr,{ssr:.6f}\n")
        print(f"[save] {out_fit}")

    print(f"\n# Total time: {time.perf_counter()-t0:.2f}s")

    # ============================================================
    # Step 6: Atomic / empirical ratio summary
    # ============================================================
    print(f"\n# Atomic / empirical magnitude ratio at high k (k=16):")
    print(f"  {'ξ':>6}  {'|μ̂_emp(k=16)|²':>16}  {'|μ̂_atomic|²':>14}  {'atomic/emp':>12}")
    for target, label, atomic_pred in targets:
        k = 16
        j = closest_dyadic(target, k)
        emp = abs(mu_hat_per_k[k][j])**2 if (0 <= j < (1<<k)) else float('nan')
        ratio = atomic_pred / emp if emp > 0 else float('nan')
        print(f"  {label:>6}  {emp:>16.6e}  {atomic_pred:>14.6e}  {ratio:>12.4f}")


if __name__ == "__main__":
    main()
