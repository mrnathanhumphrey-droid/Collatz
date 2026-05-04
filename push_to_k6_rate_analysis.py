"""
push_to_k6_rate_analysis.py — push S_k exact through k=6 (486 states) and analyze
the rate-1/2 convergence pattern more thoroughly.

Goal: confirm |ε_n| · 2^n stabilizes at a clean constant, and look for closed-form
structure in (S_n - 7/15) sequence.
"""
import sys
import os
import time
from fractions import Fraction

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


def main():
    print("# Push S_k exact through k=6 + rate analysis")
    print()
    target = Fraction(7, 15)

    X = {0: Fraction(1)}
    S = {}
    for k in [1, 2, 3, 4, 5, 6]:
        t0 = time.time()
        K, coprime = build_markov_rational(k)
        pi_q = stationary_rational(K)
        sum_pi_sq = sum(p * p for p in pi_q)
        X[k] = Fraction(3**k) * sum_pi_sq
        elapsed = time.time() - t0
        if k > 0:
            S[k] = X[k] - X[k-1]
        print(f"  k={k} ({len(coprime)} states): solved in {elapsed:.1f}s, S_{k} = {float(S.get(k, 0)):.10f}")

    print()
    print("# Exact ε_n = S_n - 7/15:")
    eps = {}
    for k in [1, 2, 3, 4, 5, 6]:
        eps[k] = S[k] - target
        print(f"  ε_{k} = {float(eps[k]):+.15e}")
    print()

    # Ratios
    print("# Ratios ε_{n+1}/ε_n:")
    for n in range(1, 6):
        ratio = float(eps[n+1] / eps[n])
        print(f"  ε_{n+1}/ε_{n} = {ratio:+.10f}")
    print()

    # |ε_n| · 2^n
    print("# |ε_n| · 2^n (should approach a constant if rate is exactly 1/2):")
    g_seq = []
    for n in range(1, 7):
        g = float(abs(eps[n])) * (2**n)
        g_seq.append(g)
        print(f"  |ε_{n}| · 2^{n} = {g:.10f}")
    print()

    # Differences in g_n: if ε_n = C·(1/2)^n + D·(1/2)^{2n} or similar, check structure
    print("# Δg_n = g_{n+1} - g_n (would be 0 if pure 1/2 rate; nonzero indicates higher-order):")
    for n in range(1, 6):
        dg = g_seq[n] - g_seq[n-1]
        print(f"  Δg_{n} = {dg:+.10e}")
    print()

    # Try to fit ε_n = A · (1/2)^n + B · (1/4)^n + C · (1/8)^n
    # Use ε_3, ε_4, ε_5, ε_6 (4 equations, 3 unknowns) — overdetermined, can verify
    # ε_n ≈ A · (1/2)^n + B · (1/4)^n + C · (1/8)^n
    # ε_3 = A/8 + B/64 + C/512
    # ε_4 = A/16 + B/256 + C/4096
    # ε_5 = A/32 + B/1024 + C/32768
    # Solve linear system
    print("# Fit ε_n ≈ A·(1/2)^n + B·(1/4)^n + C·(1/8)^n  (post-transient n ≥ 3)")
    e3, e4, e5, e6 = float(eps[3]), float(eps[4]), float(eps[5]), float(eps[6])
    import numpy as np
    M_mat = np.array([
        [1/8, 1/64, 1/512],
        [1/16, 1/256, 1/4096],
        [1/32, 1/1024, 1/32768],
    ])
    b_vec = np.array([e3, e4, e5])
    sol = np.linalg.solve(M_mat, b_vec)
    A, B, C = sol
    print(f"  Fit (using n=3,4,5):  A = {A:+.10e}, B = {B:+.10e}, C = {C:+.10e}")
    e6_pred = A/64 + B/4096 + C/262144
    print(f"  Predicted ε_6 = {e6_pred:+.10e}, actual = {e6:+.10e}, diff = {e6 - e6_pred:+.10e}")
    print()

    # Look for clean rational A
    print(f"# Is A close to a clean rational?")
    print(f"  A·30 ≈ {A * 30:.6f}")
    print(f"  A·45 ≈ {A * 45:.6f}")
    print(f"  A·105 ≈ {A * 105:.6f}")
    print(f"  A·315 ≈ {A * 315:.6f}")
    print()

    # Save extended S_k table
    out = os.path.join(OUTDIR, "S_k_exact_through_6.csv")
    with open(out, 'w') as f:
        f.write("k,S_k_decimal,eps_k_decimal,abs_eps_times_2k\n")
        for k in [1, 2, 3, 4, 5, 6]:
            f.write(f"{k},{float(S[k]):.15f},{float(eps[k]):+.15e},{float(abs(eps[k])) * 2**k:.10f}\n")
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
