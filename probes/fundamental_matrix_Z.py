"""
fundamental_matrix_Z.py — compute Kemeny-Snell fundamental matrix Z = (I − P + 1·π^T)^(-1)
exactly over Q for our Markov chain K_k on (Z/3^k Z)*.

Goal: identify whether c = 7/45 emerges as a specific entry / row-sum / trace / quadratic form
involving Z. This is the rigorous machinery hint embedded in the reference index.

For an ergodic finite chain with stationary π:
  Z = (I − P + 1·π^T)^(-1)
  Z is the FUNDAMENTAL MATRIX. Entries Z_{ij} relate to mean first passage times:
    M_{ij} = (Z_{jj} − Z_{ij}) / π_j     (first passage time i → j)
  Kemeny constant: K = Σ_j π_j · M_{ij} = trace(Z) − 1   (independent of i)
  Variance of returns: σ²_j = ... involves Z_{jj}

Candidate identities to test:
  1. trace(Z) for k=1, 2, 3 — does its limit relate to 7/45 or 7/15?
  2. π^T · Z · π or similar quadratic forms
  3. Row sums / column sums
  4. Z restricted to deviation subspace
  5. ||d_{k+1}||² in terms of Z entries
"""
import sys
import os
from fractions import Fraction
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


def matmul_Q(A, B):
    nA = len(A)
    nB = len(B[0])
    nM = len(B)
    C = [[Fraction(0) for _ in range(nB)] for _ in range(nA)]
    for i in range(nA):
        for j in range(nB):
            s = Fraction(0)
            for k in range(nM):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C


def invert_rational(A):
    """Compute matrix inverse over Q via Gauss-Jordan."""
    n = len(A)
    M = [row[:] + [Fraction(1) if i == j else Fraction(0) for j in range(n)] for i, row in enumerate(A)]
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            raise ValueError(f"Singular at col {col}")
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
        piv = M[col][col]
        for j in range(2 * n):
            M[col][j] /= piv
        for row in range(n):
            if row != col and M[row][col] != 0:
                factor = M[row][col]
                for j in range(2 * n):
                    M[row][j] -= factor * M[col][j]
    return [row[n:] for row in M]


def compute_Z(K, pi):
    """Z = (I − K + 1·π^T)^(-1)."""
    n = len(K)
    A = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            ident = Fraction(1) if i == j else Fraction(0)
            A[i][j] = ident - K[i][j] + pi[j]
    return invert_rational(A)


def trace(M):
    return sum(M[i][i] for i in range(len(M)))


def fmt_frac(q, max_denom=10000):
    """Pretty print a fraction or its decimal."""
    if q.denominator > max_denom:
        return f"{float(q):.10f}"
    return f"{q} = {float(q):.10f}"


def main():
    print("# Kemeny-Snell fundamental matrix Z analysis for our K_k chain")
    print()

    target_745 = Fraction(7, 45)
    target_715 = Fraction(7, 15)

    for k in [1, 2, 3]:
        print(f"## k = {k}")
        print()
        K, coprime = build_markov_rational(k)
        n = len(K)
        pi = stationary_rational(K)
        print(f"  states: n = {n}, coprime residues mod 3^{k}")
        print(f"  π = {[str(p) for p in pi]}")
        print()

        # Build Z
        Z = compute_Z(K, pi)

        # Trace
        tr_Z = trace(Z)
        print(f"  trace(Z) = {fmt_frac(tr_Z)}")
        # Kemeny constant K = trace(Z) - 1
        kemeny = tr_Z - Fraction(1)
        print(f"  Kemeny constant K = trace(Z) − 1 = {fmt_frac(kemeny)}")
        print(f"     vs 7/45 = {float(target_745):.10f}")
        print(f"     vs 7/15 = {float(target_715):.10f}")
        print()

        # Row sums
        row_sums = [sum(Z[i][j] for j in range(n)) for i in range(n)]
        print(f"  row sums of Z: (theory says all = 1 for ergodic chain)")
        for i in range(n):
            print(f"    row {coprime[i]:>3}: {fmt_frac(row_sums[i])}")
        print()

        # Diagonal entries
        print(f"  diagonal Z_ii:")
        for i in range(n):
            print(f"    Z[{coprime[i]}][{coprime[i]}] = {fmt_frac(Z[i][i])}")
        print()

        # Look for 7/45 in any entry
        print(f"  scanning Z entries for 7/45 = {float(target_745):.6f}, 7/15 = {float(target_715):.6f}, 10/189...")
        target_10_189 = Fraction(10, 189)
        target_10_21 = Fraction(10, 21)
        target_2_3 = Fraction(2, 3)
        for tgt_name, tgt in [("7/45", target_745), ("7/15", target_715),
                               ("10/189", target_10_189), ("10/21", target_10_21),
                               ("2/3", target_2_3), ("1/3", Fraction(1, 3))]:
            hits = []
            for i in range(n):
                for j in range(n):
                    if Z[i][j] == tgt:
                        hits.append((coprime[i], coprime[j]))
            if hits:
                print(f"    Z entries equal to {tgt_name}: {hits}")

        # π^T Z π and π^T Z (deviation forms)
        # Sum_i π_i Z_ij gives a row vector — for ergodic chain, equals π_j (Z is "averaged")
        pi_Z = [sum(pi[i] * Z[i][j] for i in range(n)) for j in range(n)]
        print()
        print(f"  π^T · Z = {[str(x) for x in pi_Z]}")

        # Z · π — column averaging
        Z_pi = [sum(Z[i][j] * pi[j] for j in range(n)) for i in range(n)]
        print(f"  Z · π = {[str(x) for x in Z_pi]}  (theory: all equal Σ π_j = 1)")

        # π^T · Z · π
        piZpi = sum(pi[i] * Z[i][j] * pi[j] for i in range(n) for j in range(n))
        print(f"  π^T · Z · π = {fmt_frac(piZpi)}")

        # Quadratic form: Σ π_i² · Z_ii (sum on diagonal weighted by π²)
        diag_pi2 = sum(pi[i] ** 2 * Z[i][i] for i in range(n))
        print(f"  Σ π_i² · Z_ii = {fmt_frac(diag_pi2)}")

        # Σ π_i² (this is Σ π² which appears in S_k formula)
        sum_pi_sq = sum(p * p for p in pi)
        X_k = Fraction(3**k) * sum_pi_sq
        print(f"  Σ π_i² = {fmt_frac(sum_pi_sq)}")
        print(f"  X_k = 3^k · Σ π² = {fmt_frac(X_k)}")

        # ||d_{k+1}||² known exact for k=1: 10/189; k=2: ...
        if k == 1:
            d2_known = Fraction(10, 189)
            print(f"  ||d_2||² = 10/189 = {float(d2_known):.10f} (exact, from R74 verification)")
            print(f"  ||d_2||² · 3 = {fmt_frac(d2_known * 3)} (vs 7/45 = {float(target_745):.10f})")

        print()
        print()

    # Asymptotic test: as k grows, what does trace(Z) / something look like?
    print("# Asymptotic check: trace(Z) / n vs k")
    print()
    for k in [1, 2, 3]:
        K, coprime = build_markov_rational(k)
        n = len(K)
        pi = stationary_rational(K)
        Z = compute_Z(K, pi)
        tr = trace(Z)
        print(f"  k={k}: n={n}, trace(Z) = {float(tr):.10f},  trace(Z)/n = {float(tr) / n:.10f}, "
              f"(trace − 1)/n = {float(tr - 1) / n:.10f}")


if __name__ == "__main__":
    main()
