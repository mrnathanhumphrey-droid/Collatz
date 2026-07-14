"""
s6_exact_fraction.py — compute the EXACT rational S_6 and S_6 - 7/15 and write them
to disk as fractions (the original push_to_k6 script only saved decimals).

Reuses Nathan's exact build_markov_rational + stationary_rational verbatim.
Run: python s6_exact_fraction.py
Takes several minutes (486-state exact rational solve, ~400s on the original run).
Output: experiments_output/S6_exact_fraction.txt
"""
import os
import time
from fractions import Fraction

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
    target = Fraction(7, 15)
    X = {}
    for k in [5, 6]:
        t0 = time.time()
        K, coprime = build_markov_rational(k)
        pi_q = stationary_rational(K)
        X[k] = Fraction(3**k) * sum(p * p for p in pi_q)
        print(f"k={k} ({len(coprime)} states) solved in {time.time()-t0:.1f}s")

    S6 = X[6] - X[5]
    eps6 = S6 - target

    lines = []
    lines.append("# Exact S_6 and S_6 - 7/15 (rational)")
    lines.append("")
    lines.append("S_6 numerator:")
    lines.append(str(S6.numerator))
    lines.append("S_6 denominator:")
    lines.append(str(S6.denominator))
    lines.append("")
    lines.append("S_6 - 7/15 numerator:")
    lines.append(str(eps6.numerator))
    lines.append("S_6 - 7/15 denominator:")
    lines.append(str(eps6.denominator))
    lines.append("")
    lines.append(f"# decimal cross-check (must match log 0.4661687610): {float(S6):.10f}")
    lines.append(f"# eps_6 decimal (must match log -4.979056652203831e-04): {float(eps6):.15e}")
    lines.append("")
    lines.append("# Lean-ready forms:")
    lines.append(f"def S6 : Q := {S6.numerator}/{S6.denominator}")
    lines.append(f"-- S6 - 7/15 = {eps6.numerator}/{eps6.denominator}")

    out = os.path.join(OUTDIR, "S6_exact_fraction.txt")
    with open(out, "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    print(f"\n[saved] {out}")


if __name__ == "__main__":
    main()
