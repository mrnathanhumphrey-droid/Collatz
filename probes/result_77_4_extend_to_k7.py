"""result_77_4_extend_to_k7.py — exact-rational ε_7 = S_7 − 7/15 via Markov chain at k=6 and k=7.

Per push_to_k6_rate_analysis.py:
  X[k] := 3^k · Σ_{r coprime to 3} π_k[r]²   (π_k = stationary at level k)
  S[k] := X[k] − X[k-1]
  ε[k] := S[k] − 7/15

Need X[6] and X[7] both as exact Fractions to get exact-rational S[7] and ε[7].

Output: experiments_output/result_77_4_S_k_exact_through_7.csv (extending the k=6 file).
"""
from __future__ import annotations
import os, sys, csv, time
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_markov_rational(k):
    N = 3 ** k
    M = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime = [r for r in range(N) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(coprime)}
    n = len(coprime)
    K = [[Fraction(0)] * n for _ in range(n)]
    Z_v = Fraction(2 ** M - 1, 2 ** M)
    for r in coprime:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2 ** r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[idx[r]][idx[target]] += p
    return K, coprime


def stationary_rational(K, label=""):
    n = len(K)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    t0 = time.time()
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
        if col % 50 == 0 or col == n - 1:
            print(f"  [{label} gauss] col {col}/{n} elapsed={time.time()-t0:.0f}s", flush=True)
    return b


def compute_X(k):
    print(f"\n[k={k}] building Markov chain ({2*3**(k-1)} states)...", flush=True)
    t0 = time.time()
    K, coprime = build_markov_rational(k)
    print(f"  built in {time.time()-t0:.1f}s", flush=True)
    pi = stationary_rational(K, label=f"k={k}")
    sum_sq = sum(p * p for p in pi)
    X = Fraction(3 ** k) * sum_sq
    return X


def main():
    T0 = time.time()
    print("=== R77.4: exact-rational ε_7 ===\n")

    X = {0: Fraction(1)}
    for k in [1, 2, 3, 4, 5, 6, 7]:
        X[k] = compute_X(k)
        S_k = X[k] - X[k - 1]
        eps_k = S_k - Fraction(7, 15)
        print(f"  k={k}: X={float(X[k]):.10f}, S={float(S_k):+.10e}, ε={float(eps_k):+.6e}, "
              f"|ε|·2^k={float(abs(eps_k))*2**k:.6f}", flush=True)

    # Output CSV with exact rationals
    output_csv = os.path.join(OUTDIR, "result_77_4_S_k_exact_through_7.csv")
    with open(output_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "S_k_decimal", "eps_k_decimal", "abs_eps_times_2k",
                    "S_k_rat_num", "S_k_rat_den", "eps_k_rat_num", "eps_k_rat_den"])
        for k in [1, 2, 3, 4, 5, 6, 7]:
            S_k = X[k] - X[k - 1]
            eps_k = S_k - Fraction(7, 15)
            w.writerow([k,
                        f"{float(S_k):.16f}",
                        f"{float(eps_k):+.16e}",
                        f"{float(abs(eps_k) * 2 ** k):.10f}",
                        str(S_k.numerator), str(S_k.denominator),
                        str(eps_k.numerator), str(eps_k.denominator)])
        f.flush()
        os.fsync(f.fileno())

    print(f"\n[save] {output_csv}")
    print(f"[done] total elapsed {time.time() - T0:.1f}s")


if __name__ == "__main__":
    main()
