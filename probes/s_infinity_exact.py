"""
s_infinity_exact.py — compute S_k exactly via rational arithmetic.

S_k = 3^k · Σ_r π_k(r)² − 3^(k-1) · Σ_r π_{k-1}(r)²

where π_k is the stationary distribution of the Markov chain on (Z/3^k Z)*
with transition K_k[r → s] = P(T(m) ≡ s mod 3^k | m ≡ r mod 3^k, v ~ Geom(1/2)).

If S_k → S∞ = 7/15 exactly, this should be visible in the rational sequence.
"""
import sys
import math
from fractions import Fraction
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


def build_markov_rational(k):
    """Build K_k as matrix of Fractions."""
    N = 3**k
    M = 2 * 3**(k-1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2**M - 1, 2**M)  # 1 - 2^(-M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2**r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime_states


def stationary_rational(K):
    """Solve π K = π exactly via Gauss elimination over Q.

    Equivalent to finding left null vector of K - I with sum = 1.
    """
    n = len(K)
    # Build augmented matrix: [K^T - I | 0] with last row replaced by sum constraint
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)

    # Gauss elimination
    for col in range(n):
        # Find pivot
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
        # Normalize row
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        # Eliminate
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    return b


def main():
    print("# Exact rational computation of S_k")
    print()
    print(f"  {'k':>3}  {'X_k = 3^k Σπ²':>20}  {'(decimal)':>12}  {'S_k':>15}  {'(decimal)':>12}")

    X_prev = Fraction(0)
    rows = []
    for k in [1, 2, 3, 4, 5]:
        K, coprime = build_markov_rational(k)
        pi = stationary_rational(K)
        sum_sq = sum(p * p for p in pi)
        X_k = Fraction(3**k) * sum_sq
        S_k = X_k - X_prev * Fraction(3**k, 3**(k-1)) if k > 1 else X_k - 1
        # Wait: S_k = X_k - X_{k-1}, where X_k = 3^k Σ π_k². Just direct sub.
        S_k = X_k - X_prev if k > 1 else X_k - 1  # at k=1: S_1 = X_1 - 1 (Parseval: X_1 = 1 + S_1)

        rows.append((k, X_k, sum_sq, S_k))
        print(f"  {k:>3}  {str(X_k):>20}  {float(X_k):>12.6f}  {str(S_k):>15}  {float(S_k):>12.6f}")
        X_prev = X_k

    print()
    print("# 7/15 = ", Fraction(7, 15), " = ", float(Fraction(7, 15)))
    print()
    print("# S_k - 7/15 (rational):")
    for k, X_k, sum_sq, S_k in rows:
        diff = S_k - Fraction(7, 15)
        print(f"  k={k}: S_k - 7/15 = {diff} = {float(diff):+.8f}")

    # Check: do S_k satisfy a recursion?
    print()
    print("# Increments S_{k+1} / S_k:")
    for i in range(1, len(rows)):
        S_prev = rows[i-1][3]
        S_curr = rows[i][3]
        ratio = S_curr / S_prev if S_prev != 0 else Fraction(0)
        print(f"  S_{rows[i][0]} / S_{rows[i-1][0]} = {ratio} = {float(ratio):.8f}")

    # Try fitting to 7/15 - c·r^k for some r
    print()
    print("# Candidates: assume S_k = 7/15 + (linear combo of geometric terms)")
    # If S_k - 7/15 → 0 geometrically, identify rate.
    # diffs: 1/105, -??, -??
    diffs = [rows[i][3] - Fraction(7, 15) for i in range(len(rows))]
    print(f"  Diffs: {[(rows[i][0], str(diffs[i]), float(diffs[i])) for i in range(len(rows))]}")

    # Compute S_∞ candidate from extrapolation
    if len(rows) >= 4:
        # Use last 3 to extrapolate
        # If S_k = A + B·r^k, fit A
        s3 = float(rows[-3][3])
        s4 = float(rows[-2][3])
        s5 = float(rows[-1][3])
        # (s5 - s4) / (s4 - s3) ≈ r
        if (s4 - s3) != 0:
            r = (s5 - s4) / (s4 - s3)
            # s5 - A = (s4 - A) · r → A = (s5 - r·s4) / (1 - r)
            if abs(1 - r) > 1e-10:
                A_extrap = (s5 - r * s4) / (1 - r)
                print(f"\n  Geometric extrapolation: r = {r:.6f}, A (= S∞) = {A_extrap:.8f}")
                print(f"  7/15 = {7/15:.8f}, |S∞ - 7/15| = {abs(A_extrap - 7/15):.2e}")


if __name__ == "__main__":
    main()
