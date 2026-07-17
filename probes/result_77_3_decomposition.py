"""
result_77_3_decomposition.py
============================
R77.3 Stage A: Test whether the 3-mode model
    eps_n = A * (1/2)^n + B * (1/4)^n + C * (1/8)^n
is EXACT in exact rational arithmetic over Q.

Procedure:
  1. Run the project's exact-rational Markov chain at k = 1..6 to get
     exact S_k as Fractions (no floats anywhere).
  2. Form eps_k = S_k - 7/15 as Fractions.
  3. Solve the 3x3 linear system over Q from eps_1, eps_2, eps_3 for (A, B, C).
  4. Predict eps_4, eps_5, eps_6 from the closed form and compare to the
     actual exact-rational eps_k.
  5. If residuals at n = 4, 5, 6 are exactly zero in Q -> outcome (alpha)
     (3-mode model is exact, Nisoli bypass succeeds at order 3).
  6. If nonzero -> outcome (beta) (3-mode model is approximate; bypass fails).

The arithmetic is exact via fractions.Fraction. Stage B (T_4 spectrum) lives
in result_77_3_T_N_eigenvalue_extension.py.
"""
from __future__ import annotations

import os
import sys
import csv
import time
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


# --------------------------------------------------------------------------- #
# Exact-rational Markov chain (replicated from push_to_k6_rate_analysis.py)   #
# --------------------------------------------------------------------------- #

def build_markov_rational(k: int):
    N = 3 ** k
    M = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2 ** M - 1, 2 ** M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2 ** r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime_states


def stationary_rational(K):
    n = len(K)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0))
          for j in range(n)] for i in range(n)]
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


# --------------------------------------------------------------------------- #
# Linear-algebra utilities (exact over Q)                                     #
# --------------------------------------------------------------------------- #

def gauss_solve_Q(A, b):
    """Solve square linear system Ax = b over Fraction. A, b mutated copies used."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            return None
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
        piv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= piv
        for row in range(n):
            if row != col and M[row][col] != 0:
                factor = M[row][col]
                for j in range(col, n + 1):
                    M[row][j] -= factor * M[col][j]
    return [M[i][n] for i in range(n)]


# --------------------------------------------------------------------------- #
# Main R77.3 Stage A                                                          #
# --------------------------------------------------------------------------- #

def main():
    print("=" * 78)
    print("R77.3 Stage A: 3-mode exact-rational decomposition test")
    print("=" * 78)
    print()
    target = Fraction(7, 15)

    # 1. Compute S_k as exact Fractions through k = 6.
    X = {0: Fraction(1)}
    S = {}
    eps = {}
    for k in [1, 2, 3, 4, 5, 6]:
        t0 = time.time()
        K, coprime = build_markov_rational(k)
        pi_q = stationary_rational(K)
        sum_pi_sq = sum(p * p for p in pi_q)
        X[k] = Fraction(3 ** k) * sum_pi_sq
        S[k] = X[k] - X[k - 1]
        eps[k] = S[k] - target
        elapsed = time.time() - t0
        print(f"  k={k} ({len(coprime)} states): solved in {elapsed:.1f}s")
        print(f"    S_{k} = {S[k]}")
        print(f"    eps_{k} = {eps[k]}    (~ {float(eps[k]):+.6e})")
    print()

    # 2. Solve the 3x3 system for (A, B, C) using eps_1, eps_2, eps_3.
    #    eps_n = A * (1/2)^n + B * (1/4)^n + C * (1/8)^n
    #
    #    Row n: [(1/2)^n, (1/4)^n, (1/8)^n]
    print("-" * 78)
    print("Solving 3-mode system from eps_1, eps_2, eps_3 over Q:")
    print("-" * 78)
    M3 = [
        [Fraction(1, 2),  Fraction(1, 4),   Fraction(1, 8)],
        [Fraction(1, 4),  Fraction(1, 16),  Fraction(1, 64)],
        [Fraction(1, 8),  Fraction(1, 64),  Fraction(1, 512)],
    ]
    b3 = [eps[1], eps[2], eps[3]]
    sol = gauss_solve_Q(M3, b3)
    A, B, C = sol
    print(f"  A = {A}    (~ {float(A):+.6e})")
    print(f"  B = {B}    (~ {float(B):+.6e})")
    print(f"  C = {C}    (~ {float(C):+.6e})")
    print()

    # 3. Predict eps_4, eps_5, eps_6 from the closed form and compare.
    print("-" * 78)
    print("Verification: predict eps_n for n = 4, 5, 6 from closed form, compare")
    print("-" * 78)
    residuals = {}
    for n in [4, 5, 6]:
        pred = (A * Fraction(1, 2 ** n)
                + B * Fraction(1, 4 ** n)
                + C * Fraction(1, 8 ** n))
        actual = eps[n]
        r = actual - pred
        residuals[n] = (pred, actual, r)
        print(f"  n = {n}:")
        print(f"    predicted = {pred}")
        print(f"               (~ {float(pred):+.6e})")
        print(f"    actual    = {actual}")
        print(f"               (~ {float(actual):+.6e})")
        print(f"    residual  = {r}")
        print(f"               (~ {float(r):+.6e})")
        if r == 0:
            print("    *** EXACT MATCH IN Q ***")
        else:
            ratio = float(abs(r) / abs(pred)) if pred != 0 else float('inf')
            print(f"    |r/pred|  = {ratio:.4e}")
        print()

    # 4. Outcome classification.
    print("=" * 78)
    print("Outcome classification (Stage A)")
    print("=" * 78)
    all_zero = all(residuals[n][2] == 0 for n in [4, 5, 6])
    if all_zero:
        outcome = "alpha (3-mode model EXACT in Q; bypass succeeds at order 3)"
    else:
        max_rel = max(
            float(abs(residuals[n][2]) / abs(residuals[n][0]))
            for n in [4, 5, 6] if residuals[n][0] != 0
        )
        if max_rel < 1e-3:
            outcome = (
                f"beta (3-mode model APPROXIMATE; max |r/pred| = {max_rel:.3e})"
            )
        else:
            outcome = (
                f"beta-large (3-mode model fails; max |r/pred| = {max_rel:.3e})"
            )
    print(f"  Outcome: {outcome}")
    print()

    # ----------------------------------------------------------------------- #
    # Stage A.4: extend to 4-mode model if residuals nonzero                  #
    # ----------------------------------------------------------------------- #
    if not all_zero:
        print("-" * 78)
        print("Stage A.4: extend to 4-mode model")
        print("  eps_n = A * (1/2)^n + B * (1/4)^n + C * (1/8)^n + D * (1/16)^n")
        print("  Solve from eps_1, eps_2, eps_3, eps_4. Predict eps_5, eps_6.")
        print("-" * 78)
        M4 = [
            [Fraction(1, 2 ** n),
             Fraction(1, 4 ** n),
             Fraction(1, 8 ** n),
             Fraction(1, 16 ** n)]
            for n in [1, 2, 3, 4]
        ]
        b4 = [eps[1], eps[2], eps[3], eps[4]]
        sol4 = gauss_solve_Q(M4, b4)
        A4, B4, C4, D4 = sol4
        print(f"  A = {A4}    (~ {float(A4):+.6e})")
        print(f"  B = {B4}    (~ {float(B4):+.6e})")
        print(f"  C = {C4}    (~ {float(C4):+.6e})")
        print(f"  D = {D4}    (~ {float(D4):+.6e})")
        print()
        residuals4 = {}
        for n in [5, 6]:
            pred = (A4 * Fraction(1, 2 ** n)
                    + B4 * Fraction(1, 4 ** n)
                    + C4 * Fraction(1, 8 ** n)
                    + D4 * Fraction(1, 16 ** n))
            actual = eps[n]
            r = actual - pred
            residuals4[n] = (pred, actual, r)
            print(f"  n = {n}:")
            print(f"    predicted = {pred}")
            print(f"               (~ {float(pred):+.6e})")
            print(f"    actual    = {actual}")
            print(f"               (~ {float(actual):+.6e})")
            print(f"    residual  = {r}")
            print(f"               (~ {float(r):+.6e})")
            if r == 0:
                print("    *** EXACT MATCH IN Q ***")
            else:
                ratio = float(abs(r) / abs(pred)) if pred != 0 else float('inf')
                print(f"    |r/pred|  = {ratio:.4e}")
            print()
        all4_zero = all(residuals4[n][2] == 0 for n in [5, 6])
        if all4_zero:
            print("  4-mode model EXACT in Q at n = 5, 6.")
        else:
            print("  4-mode model NOT exact at n = 5, 6 (residuals nonzero).")
    else:
        sol4 = None
        residuals4 = {}
        A4 = B4 = C4 = D4 = None

    # ----------------------------------------------------------------------- #
    # CSV save                                                                #
    # ----------------------------------------------------------------------- #
    out_csv = os.path.join(OUTDIR, "result_77_3_eigenvalue_data.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["section", "key", "num", "den", "decimal"])
        for k in [1, 2, 3, 4, 5, 6]:
            w.writerow(["eps", f"n={k}",
                        eps[k].numerator, eps[k].denominator, float(eps[k])])
        w.writerow(["3mode_coef", "A",
                    A.numerator, A.denominator, float(A)])
        w.writerow(["3mode_coef", "B",
                    B.numerator, B.denominator, float(B)])
        w.writerow(["3mode_coef", "C",
                    C.numerator, C.denominator, float(C)])
        for n in [4, 5, 6]:
            pred, actual, r = residuals[n]
            w.writerow(["3mode_residual", f"n={n}",
                        r.numerator, r.denominator, float(r)])
        if sol4 is not None:
            w.writerow(["4mode_coef", "A",
                        A4.numerator, A4.denominator, float(A4)])
            w.writerow(["4mode_coef", "B",
                        B4.numerator, B4.denominator, float(B4)])
            w.writerow(["4mode_coef", "C",
                        C4.numerator, C4.denominator, float(C4)])
            w.writerow(["4mode_coef", "D",
                        D4.numerator, D4.denominator, float(D4)])
            for n in [5, 6]:
                pred, actual, r = residuals4[n]
                w.writerow(["4mode_residual", f"n={n}",
                            r.numerator, r.denominator, float(r)])
    print(f"[save] {out_csv}")

    # Also write a small Python-readable cache of eps for Stage B
    eps_cache = os.path.join(OUTDIR, "result_77_3_eps_exact.txt")
    with open(eps_cache, "w") as f:
        for k in [1, 2, 3, 4, 5, 6]:
            f.write(f"{k}\t{eps[k].numerator}\t{eps[k].denominator}\n")
    print(f"[save] {eps_cache}")


if __name__ == "__main__":
    main()
