"""
subleading_eigenvalue_test.py — verify λ_2 = 1/2 conjecture for Markov chain on (Z/3^k Z)*.

If λ_2 = 1/2 exactly, this explains the empirical convergence rate ~1/2 of S_k → 7/15
and supports the rigorous S∞ = 7/15 proof.
"""
import sys
import math
import numpy as np
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")


def build_markov_chain(k):
    N = 3**k
    M = 2 * 3**(k-1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = np.zeros((n, n))
    Z_v = 1.0 - 2.0**(-M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = 2.0**(-r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r], state_idx[target]] += p
    return K, coprime_states


def main():
    print("# Subleading eigenvalue λ_2 of K_k")
    print()
    print(f"  {'k':>3}  {'N_k':>5}  {'#states':>8}  {'λ_1':>10}  {'|λ_2|':>10}  {'|λ_3|':>10}  {'top 5 |λ|':>40}")

    for k in range(1, 8):
        K, coprime = build_markov_chain(k)
        n = len(coprime)
        eigvals = np.linalg.eigvals(K)
        # Sort by absolute value descending
        abs_eigs = sorted(np.abs(eigvals), reverse=True)
        top_str = ", ".join(f"{x:.4f}" for x in abs_eigs[:5])
        l1 = abs_eigs[0]
        l2 = abs_eigs[1] if len(abs_eigs) > 1 else 0
        l3 = abs_eigs[2] if len(abs_eigs) > 2 else 0
        print(f"  {k:>3}  {3**k:>5}  {n:>8}  {l1:>10.6f}  {l2:>10.6f}  {l3:>10.6f}  {top_str:>40}")

    print()
    print("# Verify λ_2 = 1/2 conjecture")
    print(f"  Comparing |λ_2| to 1/2 = 0.5:")
    for k in range(2, 7):
        K, _ = build_markov_chain(k)
        eigvals = np.linalg.eigvals(K)
        abs_eigs = sorted(np.abs(eigvals), reverse=True)
        l2 = abs_eigs[1]
        diff = l2 - 0.5
        print(f"  k={k}: |λ_2| = {l2:.10f}, |λ_2| − 1/2 = {diff:+.2e}")


if __name__ == "__main__":
    main()
