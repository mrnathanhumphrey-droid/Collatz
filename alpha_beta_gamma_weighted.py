"""
alpha_beta_gamma_weighted.py — compute π²-weighted ψ and verify S_k formula.

S_k = 3^k · Σ_r π_{k-1}(r)² · (ψ_r - 1/3)

where ψ_r = α̃² + β̃² + γ̃² for the sub-cell shares at residue r mod 3^(k-1).

The π²-weighted ⟨ψ⟩ is what governs S_k convergence, NOT the unweighted ψ_avg.
"""
import sys
import math
import os
from fractions import Fraction
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_markov_rational(k):
    N = 3**k
    M = 2 * 3**(k-1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2**M - 1, 2**M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2**r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime_states


def stationary_rational(K):
    n = len(K)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            raise ValueError(f"Singular at col {col}")
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    return b


def stationary_numerical(k):
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
    pi = np.ones(n) / n
    for _ in range(200):
        pi = pi @ K
        pi = pi / pi.sum()
    return pi, coprime_states


def main():
    print("# π²-weighted ψ analysis")
    print()

    # Get stationary at all levels
    pis = {}
    for k in [1, 2, 3, 4, 5, 6]:
        if k <= 4:
            K, coprime = build_markov_rational(k)
            pi_q = stationary_rational(K)
            pi = [float(p) for p in pi_q]
            pis[k] = (pi, coprime, pi_q)
        else:
            pi, coprime = stationary_numerical(k)
            pis[k] = (list(pi), coprime, None)
        print(f"  k={k}: {len(pis[k][1])} states")

    # For each k → k+1 transition: compute ψ_r for each r mod 3^k, then π²-weighted average
    print()
    print(f"# Verify S_k formula: S_k = 3^k · Σ_r π_{{k-1}}²(r) · (ψ_r - 1/3)")
    print()
    print(f"  {'k':>3}  {'⟨ψ⟩_unweighted':>16}  {'⟨ψ⟩_π²-weighted':>17}  {'⟨ψ-1/3⟩_π²':>13}  {'X_{{k-1}}':>10}  {'S_k computed':>14}  {'S_k expected':>14}")

    # S_k known:
    expected = {1: 2/3, 2: 10/21, 3: 31370/67963, 4: 0.4642, 5: 0.4655}

    for k in [1, 2, 3, 4, 5]:
        pi_k, coprime_k, _ = pis[k]
        pi_kp1, coprime_kp1, _ = pis[k+1]
        N_k = 3**k
        state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}

        # Compute ψ_r for each r in coprime_k
        psi_list = []  # (r, π_k(r), ψ_r)
        for i, r in enumerate(coprime_k):
            lifts = [r, r + N_k, r + 2 * N_k]
            masses = [pi_kp1[state_idx_kp1[lift]] for lift in lifts if lift in state_idx_kp1]
            total = sum(masses)
            if total == 0:
                continue
            shares = [m / total for m in masses]
            psi = sum(s**2 for s in shares)
            psi_list.append((r, pi_k[i], psi))

        # Unweighted average ψ
        psi_uw = sum(p[2] for p in psi_list) / len(psi_list)

        # π²-weighted average ψ
        sum_pi2 = sum(p[1]**2 for p in psi_list)
        sum_pi2_psi = sum(p[1]**2 * p[2] for p in psi_list)
        psi_w = sum_pi2_psi / sum_pi2

        # ⟨ψ - 1/3⟩_π² = (sum_pi2_psi - sum_pi2/3) / sum_pi2
        psi_dev_w = psi_w - 1/3

        # X_{k-1} = X at level k = 3^k · Σ π_k² (this is the X for the level we're at, but in S_(k+1) formula it's X_k = previous level's X)
        # For computing S_(k+1) using k:
        # S_{k+1} = 3 · X_k · ⟨ψ - 1/3⟩_{π_k²-weighted}
        #        = 3 · 3^k · Σ π_k² · ⟨ψ - 1/3⟩_w
        #        = 3 · 3^k · sum_pi2_psi - 3^k · sum_pi2 — wait let me redo.

        # S_{k+1} = 3^{k+1} Σ π_{k+1}² - 3^k Σ π_k²
        # Σ π_{k+1}² = Σ_r π_k(r)² · ψ_r = sum_pi2_psi
        # 3^{k+1} sum_pi2_psi - 3^k sum_pi2 = 3^k (3·sum_pi2_psi - sum_pi2) = 3^k · sum_pi2 · (3·psi_w - 1)

        # X_k = 3^k Σ π_k² = 3^k · sum_pi2
        X_k = 3**k * sum_pi2
        S_kp1 = X_k * (3 * psi_w - 1)

        S_kp1_expected = expected.get(k+1, float('nan'))

        print(f"  {k:>3}  {psi_uw:>16.6f}  {psi_w:>17.6f}  {psi_dev_w:>13.6f}  {X_k:>10.6f}  {S_kp1:>14.6f}  {S_kp1_expected:>14.6f}")

    # Now: as k → ∞, what does π²-weighted ψ approach?
    print()
    print(f"# Convergence of π²-weighted ψ to 1/3:")
    print(f"  {'k':>3}  {'⟨ψ⟩_w':>10}  {'⟨ψ-1/3⟩_w':>12}  {'ratio_(k+1)/k':>14}")

    weighted_data = []
    for k in [1, 2, 3, 4, 5]:
        pi_k, coprime_k, _ = pis[k]
        pi_kp1, coprime_kp1, _ = pis[k+1]
        N_k = 3**k
        state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}

        psi_list = []
        for i, r in enumerate(coprime_k):
            lifts = [r, r + N_k, r + 2 * N_k]
            masses = [pi_kp1[state_idx_kp1[lift]] for lift in lifts if lift in state_idx_kp1]
            total = sum(masses)
            if total == 0:
                continue
            shares = [m / total for m in masses]
            psi = sum(s**2 for s in shares)
            psi_list.append((r, pi_k[i], psi))

        sum_pi2 = sum(p[1]**2 for p in psi_list)
        sum_pi2_psi = sum(p[1]**2 * p[2] for p in psi_list)
        psi_w = sum_pi2_psi / sum_pi2
        weighted_data.append((k, psi_w))

    prev = None
    for (k, psi_w) in weighted_data:
        dev = psi_w - 1/3
        ratio = (psi_w - 1/3) / prev if prev is not None and prev != 0 else float('nan')
        if not isinstance(ratio, str):
            print(f"  {k:>3}  {psi_w:>10.6f}  {dev:>12.6f}  {ratio:>14}")
        else:
            print(f"  {k:>3}  {psi_w:>10.6f}  {dev:>12.6f}  {ratio:>14}")
        prev = dev

    # The asymptotic behavior:
    # If S_k → 7/15 = 0.4667, then with X_{k-1} growing linearly as k·S∞:
    # ⟨ψ-1/3⟩_w → S_∞ / (3 · X_{k-1}) → 7/15 / (3 · k · 7/15) = 1/(3k)
    print()
    print(f"# Predicted ⟨ψ-1/3⟩_w from S∞ = 7/15:")
    print(f"  Asymptotic: ⟨ψ-1/3⟩_w · X_{{k}} → S_∞/3 = 7/45 = {7/45:.6f}")
    print(f"\n  {'k':>3}  {'⟨ψ-1/3⟩_w · X_k':>16}  {'7/45 = 0.155556':>14}  {'diff':>10}")
    for k in [1, 2, 3, 4, 5]:
        pi_k, coprime_k, _ = pis[k]
        pi_kp1, coprime_kp1, _ = pis[k+1]
        N_k = 3**k
        state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}

        psi_list = []
        for i, r in enumerate(coprime_k):
            lifts = [r, r + N_k, r + 2 * N_k]
            masses = [pi_kp1[state_idx_kp1[lift]] for lift in lifts if lift in state_idx_kp1]
            total = sum(masses)
            if total == 0:
                continue
            shares = [m / total for m in masses]
            psi = sum(s**2 for s in shares)
            psi_list.append((r, pi_k[i], psi))

        sum_pi2 = sum(p[1]**2 for p in psi_list)
        sum_pi2_psi = sum(p[1]**2 * p[2] for p in psi_list)
        psi_w = sum_pi2_psi / sum_pi2
        X_k = 3**k * sum_pi2
        product = (psi_w - 1/3) * X_k
        diff = product - 7/45
        print(f"  {k:>3}  {product:>16.6f}  {7/45:>14.6f}  {diff:>10.4e}")


if __name__ == "__main__":
    main()
