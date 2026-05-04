"""
decay_law_derivation.py — analytical derivation of |μ̂(a/3^k)|² for k=1, 2, 3.

Method: build Markov chain on coprime-to-3 residues mod 3^k. Transition kernel:
  K[r → s] = P(T(m) ≡ s mod 3^k | m ≡ r mod 3^k, v ~ Geom(1/2))

  T(m) = (3m+1) · 2^(-v) mod 3^k

For odd m, v = v_2(3m+1) ≥ 1 with P(v=j) = 2^(-j) (Lagarias-Sinai heuristic).

Two key identities for mod 3^k:
1. (3m+1) mod 3^k depends only on m mod 3^(k-1)
2. 2^(-v) mod 3^k cycles with period ord_{3^k}(2) = 2·3^(k-1)

Stationary π_r on coprime-to-3 residues. Then:
  μ̂(a/3^k) = Σ_r π_r · exp(2πi · a · r / 3^k)
  |μ̂(a/3^k)|² closed form

Compare to empirical from tree at N=2^22.
"""
import csv
import math
import os
import sys
import time
import cmath
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
    by_d = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    s = {m: 1 for m in tree}
    for m in by_d:
        for c in children[m]:
            s[m] += s[c]
    return s


def build_markov_chain(k):
    """Build transition matrix K on coprime-to-3 residues mod 3^k."""
    N = 3**k
    M = 2 * 3**(k-1)  # ord_{3^k}(2)
    inv2 = pow(2, -1, N)  # 2^(-1) mod N
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]  # 2^(-v) mod N for v=1..M

    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n_states = len(coprime_states)

    K = np.zeros((n_states, n_states))
    # P(v ≡ r_v mod M, v ≥ 1) = 2^(-r_v) / (1 - 2^(-M))
    Z_v = 1.0 - 2.0**(-M)  # normalization

    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = 2.0**(-r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r], state_idx[target]] += p

    return K, coprime_states


def stationary_distribution(K):
    """Solve π = πK. Returns left eigenvector with eigenvalue 1."""
    eigvals, eigvecs = np.linalg.eig(K.T)  # K.T for left eigvec
    idx = np.argmax(np.real(eigvals))
    pi = np.real(eigvecs[:, idx])
    pi = pi / pi.sum()
    return pi


def fourier_at(a, k, pi, states):
    """μ̂(a/3^k) given stationary distribution on coprime-to-3 states."""
    N = 3**k
    s = 0.0 + 0.0j
    for r, p in zip(states, pi):
        s += p * cmath.exp(2j * math.pi * a * r / N)
    return s


def empirical_partition(odd_ints, weights, k):
    """Empirical f_r = P_r / Z mod 3^k."""
    N = 3**k
    P = defaultdict(float)
    for m, w in zip(odd_ints, weights):
        P[int(m) % N] += w
    Z = sum(P.values())
    return {r: P[r] / Z for r in range(N)}


def main():
    print("# Analytical derivation of 3-adic decay law")
    N_max = 1 << 22
    t0 = time.perf_counter()
    tree = build_inverse_tree(N_max)
    sizes = subtree_sizes(tree)
    odd_ints = np.array([m for m in tree if m & 1], dtype=np.int64)
    weights = np.array([float(sizes[m]) for m in odd_ints])
    Z = weights.sum()
    print(f"# Tree built: {len(odd_ints)} nodes, Z={Z:.4e}, time {time.perf_counter()-t0:.2f}s")

    K_LEVELS = [1, 2, 3, 4]

    derived_results = {}
    empirical_results = {}

    for k in K_LEVELS:
        N_q = 3**k
        M = 2 * 3**(k-1)
        print(f"\n# === k = {k}, q = 3^k = {N_q}, ord_q(2) = {M} ===")

        # Build Markov chain
        K, coprime = build_markov_chain(k)
        print(f"  Markov chain: {len(coprime)} coprime states")

        # Stationary
        pi = stationary_distribution(K)
        print(f"  Stationary distribution π_r (sum = {pi.sum():.6f}):")
        if k == 1:
            for r, p in zip(coprime, pi):
                print(f"    r={r:>2}: π_r = {p:.6f}")
        elif k == 2:
            print(f"    {dict(zip(coprime, [f'{p:.4f}' for p in pi]))}")
            print(f"    As fractions of 63: {dict(zip(coprime, [f'{p*63:.2f}/63' for p in pi]))}")
        elif k == 3:
            print(f"    18 states, max π = {pi.max():.4f}, min = {pi.min():.4f}")
        else:
            print(f"    {len(coprime)} states, max π = {pi.max():.4f}, min = {pi.min():.4f}")

        # Compute analytical |μ̂(a/3^k)|² for each primitive a
        primitives = [a for a in range(1, N_q) if math.gcd(a, N_q) == 1]
        # For 3^k: gcd(a, 3^k) = 1 iff a not multiple of 3
        analytical_sq = {}
        for a in primitives:
            mh = fourier_at(a, k, pi, coprime)
            analytical_sq[a] = abs(mh)**2

        # Empirical f_r
        f_emp = empirical_partition(odd_ints, weights, k)

        # Compute empirical |μ̂(a/3^k)|² (closed form: π via empirical f)
        empirical_sq = {}
        for a in primitives:
            s = 0.0 + 0.0j
            for r in range(N_q):
                s += f_emp[r] * cmath.exp(2j * math.pi * a * r / N_q)
            empirical_sq[a] = abs(s)**2

        # Print comparison
        print(f"\n  |μ̂(a/3^{k})|² for each primitive a:")
        print(f"    {'a':>4}  {'analytical':>12}  {'empirical':>12}  {'ratio':>8}")
        for a in primitives:
            ratio = analytical_sq[a] / empirical_sq[a] if empirical_sq[a] > 1e-15 else float('nan')
            print(f"    {a:>4}  {analytical_sq[a]:>12.6f}  {empirical_sq[a]:>12.6f}  {ratio:>8.4f}")

        avg_anal = sum(analytical_sq.values()) / len(analytical_sq)
        avg_emp = sum(empirical_sq.values()) / len(empirical_sq)
        max_anal = max(analytical_sq.values())
        max_emp = max(empirical_sq.values())
        min_anal = min(analytical_sq.values())
        min_emp = min(empirical_sq.values())
        print(f"\n  Averages:  analytical = {avg_anal:.6f}, empirical = {avg_emp:.6f}")
        print(f"  Max:       analytical = {max_anal:.6f}, empirical = {max_emp:.6f}")
        print(f"  Min:       analytical = {min_anal:.6f}, empirical = {min_emp:.6f}")

        derived_results[k] = {'analytical': analytical_sq, 'empirical': empirical_sq,
                              'pi': pi, 'states': coprime, 'avg_anal': avg_anal,
                              'avg_emp': avg_emp, 'max_anal': max_anal, 'max_emp': max_emp,
                              'min_anal': min_anal, 'min_emp': min_emp}

    # ============================================================
    # Decay analysis
    # ============================================================
    print(f"\n\n# === Decay law analysis ===")
    print(f"\n  {'k':>3}  {'avg_anal':>10}  {'max_anal':>10}  {'min_anal':>10}  {'avg_emp':>10}  {'4·avg_prev':>10}  {'ratio':>8}")
    prev_avg = None
    for k in K_LEVELS:
        d = derived_results[k]
        ratio_4 = (prev_avg * 4 if prev_avg else float('nan'))
        ratio = d['avg_anal'] / prev_avg if prev_avg else float('nan')
        print(f"  {k:>3}  {d['avg_anal']:>10.6f}  {d['max_anal']:>10.6f}  {d['min_anal']:>10.6f}  {d['avg_emp']:>10.6f}  {ratio_4 if prev_avg else 'n/a':>10}  {1/ratio if prev_avg else float('nan'):>8.4f}")
        prev_avg = d['avg_anal']

    print(f"\n  Decay ratio analysis (analytical):")
    for i in range(1, len(K_LEVELS)):
        k1, k2 = K_LEVELS[i-1], K_LEVELS[i]
        a1 = derived_results[k1]['avg_anal']
        a2 = derived_results[k2]['avg_anal']
        ratio = a1 / a2  # higher k → smaller, so a1 > a2 if decay
        print(f"    k={k1} → k={k2}: ratio (avg_k1 / avg_k2) = {ratio:.4f}")

    # Conjecture: 0.31 × 4^(-(k-1))
    print(f"\n  Conjecture R65: |μ̂(a/3^k)|² ≈ 0.31 × 4^(-(k-1))")
    print(f"  {'k':>3}  {'conjecture':>10}  {'avg_anal':>10}  {'max_anal':>10}")
    for k in K_LEVELS:
        conj = 0.31 * 4**(-(k-1))
        d = derived_results[k]
        print(f"  {k:>3}  {conj:>10.6f}  {d['avg_anal']:>10.6f}  {d['max_anal']:>10.6f}")

    # Save results
    out = os.path.join(OUTDIR, "decay_law_derivation.csv")
    with open(out, 'w') as f:
        f.write("k,a,analytical_mu_hat_sq,empirical_mu_hat_sq\n")
        for k in K_LEVELS:
            d = derived_results[k]
            for a in sorted(d['analytical'].keys()):
                f.write(f"{k},{a},{d['analytical'][a]:.10f},{d['empirical'][a]:.10f}\n")
    print(f"\n  [save] {out}")

    out_pi = os.path.join(OUTDIR, "decay_law_stationary.csv")
    with open(out_pi, 'w') as f:
        f.write("k,r,pi_r\n")
        for k in K_LEVELS:
            d = derived_results[k]
            for r, p in zip(d['states'], d['pi']):
                f.write(f"{k},{r},{p:.10f}\n")
    print(f"  [save] {out_pi}")

    print(f"\nTotal time: {time.perf_counter()-t0:.2f}s")


if __name__ == "__main__":
    main()
