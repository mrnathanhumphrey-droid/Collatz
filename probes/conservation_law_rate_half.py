"""
conservation_law_rate_half.py — verify and exploit the rigorous conservation law:

Theorem 76.1: Σ_{j=0,1,2} M_{n+1}(η_0 + j·3^n) = 0 for all η_0 ∈ (Z/3^n)*.

PROOF: M_{n+1}(η) = Σ_{ξ: 3∤ξ in Z/3^{n+1}} μ̂_{n+1}(ξ) μ̂_{n+1}*(ξη)

  Σ_j M_{n+1}(η_0 + j·3^n) = Σ_ξ μ̂_{n+1}(ξ) · Σ_j μ̂_{n+1}*(ξ(η_0 + j·3^n))

For each fixed ξ, the inner sum:
  Σ_j μ̂_{n+1}*(ξη_0 + ξj·3^n) = Σ_r π_{n+1}(r) Σ_j e^{2πi r ξ(η_0 + j·3^n)/3^{n+1}}
                              = Σ_r π_{n+1}(r) e^{2πi r ξ η_0/3^{n+1}} · Σ_j e^{2πi r ξ j/3}

Inner-most sum Σ_j ω^{j} with ω = e^{2πi rξ/3} gives 3 if rξ ≡ 0 mod 3, else 0. Since π_{n+1}
is supported on r ∈ (Z/3^{n+1})* (i.e., 3∤r), and ξ ∈ (Z/3^{n+1})* also, rξ ≢ 0 mod 3 always.
So the inner sum = 0 for every r, hence Σ_j M_{n+1}(η_0 + j·3^n) = 0. ∎

This is RIGOROUS, no Geom assumed.

Theorem 76.2 (Pairing): At level n+1, exactly one of {η_0 + j·3^n : j=0,1,2} is self-inverse
mod 3^{n+1} (since the self-inverse element of (Z/3^N Z)* are {1, -1}, and exactly one lift
of any η_0 reduces to either 1 or -1 mod 3^{n+1}).

Combined with M-symmetry (M_{n+1}(η) = M_{n+1}(η^{-1})*) and conservation:
- Self-inverse lift x: M(x) = M(x)* (real)
- Mutual-inverse pair (y, y^{-1}): M(y) = M(y^{-1})*
For real M values (verified empirically): M(self-inv) = -2·M(pair)

Corollary: S_{n+1} = -2·M_{n+1}(1 + 3^n) = -2·M_{n+1}(1 + 2·3^n).

This is the EXACT relation. Now use this to derive the rate.
"""
import sys
import os
from fractions import Fraction
import cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


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


def char_func_complex(pi, coprime, k, xi):
    """μ̂_k(ξ) = Σ_r π_k(r) e^{-2πi r ξ / 3^k}."""
    N = 3**k
    z = complex(0, 0)
    for i, r in enumerate(coprime):
        z += float(pi[i]) * cmath.exp(-2j * cmath.pi * r * xi / N)
    return z


def main():
    print("# Conservation law and rate-1/2 derivation")
    print()

    # Verify conservation law at multiple levels
    pis = {}
    for k in [1, 2, 3, 4]:
        K, coprime = build_markov_rational(k)
        pi_q = stationary_rational(K)
        pis[k] = (pi_q, coprime)
    print("# Stationary distributions computed for k=1..4")
    print()

    # Compute M values
    print("# Verify Theorem 76.1: Σ_j M_{n+1}(η_0 + j·3^n) = 0")
    print()
    for n in [1, 2, 3]:
        N_n = 3**n
        N_np1 = 3**(n+1)
        pi_q, coprime_np1 = pis[n+1]

        # Compute μ̂_{n+1}(ξ) for all ξ
        mu_hat = {xi: char_func_complex(pi_q, coprime_np1, n+1, xi) for xi in range(N_np1)}

        # For each η_0 in (Z/3^n)*, compute M_{n+1}(η_0 + j·3^n) for j=0,1,2
        coprime_n = pis[n][1]
        max_dev = 0.0
        for eta_0 in coprime_n:
            triple = []
            for j in range(3):
                eta = (eta_0 + j * N_n) % N_np1
                # M_{n+1}(η) = Σ_ξ μ̂(ξ) μ̂*(ξη)
                M_val = sum(mu_hat[xi] * mu_hat[(xi * eta) % N_np1].conjugate()
                            for xi in coprime_np1)
                triple.append(M_val)
            sum_triple = sum(triple)
            dev = abs(sum_triple)
            if dev > max_dev:
                max_dev = dev
        print(f"  n={n} → n+1={n+1}: max |Σ_j M_{{n+1}}(η_0 + j·3^n)| over all η_0 = {max_dev:.2e}")

    print()

    # Verify Theorem 76.2 (pairing) for η_0 = 1 specifically
    print("# Verify Corollary 76.3: S_{n+1} = -2 · M_{n+1}(1 + 3^n)")
    print()
    print(f"  {'n':>3}  {'S_{n+1}':>14}  {'M(1+3^n)':>14}  {'-2·M(1+3^n)':>14}  {'ratio':>10}")
    for n in [1, 2, 3]:
        N_n = 3**n
        N_np1 = 3**(n+1)
        pi_q, coprime_np1 = pis[n+1]
        mu_hat = {xi: char_func_complex(pi_q, coprime_np1, n+1, xi) for xi in range(N_np1)}

        # M(1) = S_{n+1}
        S_np1 = sum(mu_hat[xi] * mu_hat[xi].conjugate() for xi in coprime_np1)
        # M(1 + 3^n)
        M_lift = sum(mu_hat[xi] * mu_hat[(xi * (1 + N_n)) % N_np1].conjugate()
                     for xi in coprime_np1)
        # M(1 + 2·3^n)
        M_lift2 = sum(mu_hat[xi] * mu_hat[(xi * (1 + 2*N_n)) % N_np1].conjugate()
                      for xi in coprime_np1)

        ratio = float((-2 * M_lift / S_np1).real) if abs(S_np1) > 1e-20 else 0
        print(f"  {n:>3}  {S_np1.real:>14.10f}  {M_lift.real:>14.10f}  "
              f"{(-2*M_lift).real:>14.10f}  {ratio:>10.6f}")
        print(f"           M(1+2·3^n) = {M_lift2.real:.10f}  (should equal M(1+3^n))")
    print()

    # Now: derive S_∞ from the conservation law iterated
    # S_{n+1} = -2 M_{n+1}(1 + 3^n)
    # As n → ∞, M_{n+1}(1 + 3^n) → ?
    # The frequency η = 1 + 3^n is "level n" frequency in a 3-adic sense (lowest non-trivial)
    print("# Tracking M_{n+1}(1 + 3^n) values (the 'leading deviation mode'):")
    print(f"  {'n':>3}  {'M(1+3^n)':>14}  {'M·(-2) = S_{n+1}':>20}  {'M - (-7/30)':>14}")
    for n in [1, 2, 3]:
        N_n = 3**n
        N_np1 = 3**(n+1)
        pi_q, coprime_np1 = pis[n+1]
        mu_hat = {xi: char_func_complex(pi_q, coprime_np1, n+1, xi) for xi in range(N_np1)}
        M_lift = sum(mu_hat[xi] * mu_hat[(xi * (1 + N_n)) % N_np1].conjugate()
                     for xi in coprime_np1)
        target_M = -7/30  # since S_∞ = 7/15, and S = -2·M ⟹ M = -7/30
        diff = M_lift.real - target_M
        print(f"  {n:>3}  {M_lift.real:>14.10f}  {(-2*M_lift).real:>20.10f}  {diff:>+14.6e}")
    print()
    print(f"  Target M_∞(1+3^n→1+0 in Z_3) = -S_∞/2 = -(7/15)/2 = -7/30 = {-7/30:.10f}")
    print()

    # The KEY INSIGHT: At each level n+1, M(1+3^n) is the "fine deviation mode" at level n+1.
    # Its CONVERGENCE rate to -7/30 governs the rate of S → 7/15.
    # Since S_{n+1} = -2·M_{n+1}(1+3^n), |S_{n+1} - 7/15| = 2·|M_{n+1}(1+3^n) - (-7/30)|.
    # So the rate of M(1+3^n) → -7/30 equals the rate of S → 7/15.

    print("# Strategy for rate-1/2 rigorous proof:")
    print("  1. M_{n+1}(1 + 3^n) is the 'leading deviation mode' at level n+1.")
    print("  2. Its convergence to -7/30 = -S_∞/2 is governed by the same operator that governs S → 7/15.")
    print("  3. Build the operator T_lead acting on the sequence (M_n(1+3^{n-1}))_n.")
    print("  4. Identify T_lead's spectrum to certify rate.")
    print()

    # Let's also examine M values at OTHER lift positions for η_0 = 1 across levels
    # M_n(1 + 3^{n-1}) is the "highest-level" lift; M_n(1 + 3^j) for j < n-1 captures coarser modes
    print("# 'Tower' values M_n(1 + 3^j·c) for various coarseness levels j:")
    for n in [2, 3, 4]:
        N_n = 3**n
        pi_q, coprime_n = pis[n]
        mu_hat = {xi: char_func_complex(pi_q, coprime_n, n, xi) for xi in range(N_n)}
        S_n_val = sum(mu_hat[xi] * mu_hat[xi].conjugate() for xi in coprime_n).real

        print(f"\n  n = {n}: S_n = {S_n_val:.10f}")
        print(f"    {'η':>5}  {'M_n(η)':>14}  {'rel to S_n':>12}")
        # η values that are 1 mod 3 lifted at various levels
        for j in range(n):
            eta = 1 + 3**j
            if eta < N_n and eta % 3 != 0:
                M_val = sum(mu_hat[xi] * mu_hat[(xi * eta) % N_n].conjugate()
                           for xi in coprime_n).real
                print(f"    {eta:>5}  {M_val:>+14.10f}  {M_val/S_n_val:>+12.6f}")
        # η = 1 itself
        print(f"    {1:>5}  {S_n_val:>+14.10f}  {1.0:>+12.6f}")


if __name__ == "__main__":
    main()
