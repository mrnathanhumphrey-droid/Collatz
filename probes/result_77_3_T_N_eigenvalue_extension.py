"""
result_77_3_T_N_eigenvalue_extension.py
=======================================
R77.3 Stage B: Build T_3 (and attempted T_4) companion matrices over Q from
the empirical eps_n sequence and read off their spectra.

The two structural questions Stage B answers:

(B.1) Direct recursion test.  Does the order-3 recursion
        eps_{n+3} = (7/8) eps_{n+2} - (7/32) eps_{n+1} + (1/64) eps_n
      (whose companion matrix has spectrum {1/2, 1/4, 1/8} over Q) hold
      EXACTLY for the eps_n produced by the project's Markov chain?
      We check at n = 1 (predicts eps_4), n = 2 (predicts eps_5), n = 3
      (predicts eps_6).  If equality is exact in Q, the order-3 model is
      exact; if not, the model is approximate.

(B.2) Empirical recursion fit.  Fit the order-3 recursion coefficients
      (alpha_2, alpha_1, alpha_0) from windows of eps_n itself (3 consecutive
      windows -> 3 equations -> unique solution).  Compare to the predicted
      (7/8, -7/32, 1/64).  If they match: T_3's predicted spec {1/2, 1/4, 1/8}
      is what the data actually says.  If not: read off the empirical
      companion matrix's spectrum.

(B.3) Order-4 attempt.  Fit a 4-mode closed form
        eps_n = A (1/2)^n + B (1/4)^n + C (1/8)^n + D (1/16)^n
      from eps_1..eps_4 (4 unknowns, 4 equations -> unique).  Predict eps_5,
      eps_6 and check residuals.  If residuals vanish: T_4 with spec
      {1/2, 1/4, 1/8, 1/16} is consistent with the data.  Otherwise: spectrum
      doesn't extend cleanly along (1/2)^k.

NOTE on execution: same harness as the rest of R77; if `python` is denied
this script documents what it would compute analytically (see the companion
markdown `result_77_3_nisoli_bypass.md`).
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
# Re-run Markov chain at k = 1..6 and extract eps_n exactly                   #
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


def get_eps_exact():
    target = Fraction(7, 15)
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
        print(f"  k = {k}: solved in {time.time()-t0:.1f}s")
    return eps


# --------------------------------------------------------------------------- #
# Linear-algebra helpers (exact over Q)                                       #
# --------------------------------------------------------------------------- #

def gauss_solve_Q(A, b):
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


def companion_from_recursion(coeffs):
    """coeffs = [alpha_{r-1}, alpha_{r-2}, ..., alpha_0]."""
    r = len(coeffs)
    C = [[Fraction(0)] * r for _ in range(r)]
    C[0] = list(coeffs)
    for i in range(1, r):
        C[i][i - 1] = Fraction(1)
    return C


def char_poly_rational_eigvals(coeffs):
    """For r=3 with coefficients (a2, a1, a0):
       char poly = lambda^3 - a2 lambda^2 + a1 lambda - a0  (sign convention
       for `companion_from_recursion`).

       Try rational roots of the form +/- p/q where p | constant_term,
       q | leading.  Caller passes the symbolic char poly coefficients
       [1, -a2, a1, -a0] for r = 3 ; we factor by trial.

       For the predicted (7/8, -7/32, 1/64), the roots are {1/2, 1/4, 1/8} per
       (lambda - 1/2)(lambda - 1/4)(lambda - 1/8).
    """
    raise NotImplementedError("Use sympy or hand-factor; see markdown writeup")


# --------------------------------------------------------------------------- #
# Main: Stages B.1, B.2, B.3                                                  #
# --------------------------------------------------------------------------- #

def main():
    print("=" * 78)
    print("R77.3 Stage B: companion-matrix spectrum tests")
    print("=" * 78)
    print()

    eps = get_eps_exact()

    # ----------------------------------------------------------------------- #
    # B.1: predicted recursion test                                           #
    #      Does eps_{n+3} = (7/8) eps_{n+2} - (7/32) eps_{n+1} + (1/64) eps_n #
    #      hold exactly at n = 1, 2, 3 (predicting eps_4, eps_5, eps_6)?      #
    # ----------------------------------------------------------------------- #
    print("-" * 78)
    print("B.1  Predicted recursion test (alpha = 7/8, -7/32, 1/64)")
    print("-" * 78)
    a2, a1, a0 = Fraction(7, 8), Fraction(-7, 32), Fraction(1, 64)
    for n in [1, 2, 3]:
        pred = a2 * eps[n + 2] + a1 * eps[n + 1] + a0 * eps[n]
        actual = eps[n + 3]
        r = actual - pred
        print(f"  n = {n} -> eps_{n+3}:")
        print(f"    pred   = {pred}    (~ {float(pred):+.6e})")
        print(f"    actual = {actual}    (~ {float(actual):+.6e})")
        print(f"    r      = {r}    (~ {float(r):+.6e})")
        if r == 0:
            print("    *** EXACT ***")
        else:
            print(f"    |r/pred| = {float(abs(r)/abs(pred)):.4e}")
        print()

    # ----------------------------------------------------------------------- #
    # B.2: empirical recursion fit                                            #
    #      Solve for (a2, a1, a0) from the 3 windows                          #
    #         eps_4 = a2 eps_3 + a1 eps_2 + a0 eps_1                          #
    #         eps_5 = a2 eps_4 + a1 eps_3 + a0 eps_2                          #
    #         eps_6 = a2 eps_5 + a1 eps_4 + a0 eps_3                          #
    # ----------------------------------------------------------------------- #
    print("-" * 78)
    print("B.2  Empirical order-3 recursion fit")
    print("-" * 78)
    M = [
        [eps[3], eps[2], eps[1]],
        [eps[4], eps[3], eps[2]],
        [eps[5], eps[4], eps[3]],
    ]
    b = [eps[4], eps[5], eps[6]]
    sol = gauss_solve_Q(M, b)
    a2_emp, a1_emp, a0_emp = sol
    print(f"  a2 = {a2_emp}    (~ {float(a2_emp):+.6f})    [predicted 7/8 = {float(Fraction(7,8))}]")
    print(f"  a1 = {a1_emp}    (~ {float(a1_emp):+.6f})    [predicted -7/32 = {float(Fraction(-7,32))}]")
    print(f"  a0 = {a0_emp}    (~ {float(a0_emp):+.6f})    [predicted 1/64 = {float(Fraction(1,64))}]")
    print()
    print(f"  Match predicted? {a2_emp == Fraction(7,8) and a1_emp == Fraction(-7,32) and a0_emp == Fraction(1,64)}")
    print()
    # The empirical companion T_3 (from data) and its char poly.
    T3_emp = companion_from_recursion([a2_emp, a1_emp, a0_emp])
    print("  Empirical T_3 (from data):")
    for row in T3_emp:
        print(f"    {[str(x) for x in row]}")
    print(f"  char poly: lambda^3 - ({a2_emp}) lambda^2 - ({a1_emp}) lambda - ({a0_emp})")
    print()

    # ----------------------------------------------------------------------- #
    # B.3: 4-mode test                                                        #
    #      eps_n = A (1/2)^n + B (1/4)^n + C (1/8)^n + D (1/16)^n             #
    #      Solve from eps_1..eps_4, predict eps_5, eps_6                      #
    # ----------------------------------------------------------------------- #
    print("-" * 78)
    print("B.3  4-mode closed-form test (T_4 with spec {1/2, 1/4, 1/8, 1/16})")
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
    for n in [5, 6]:
        pred = (A4 * Fraction(1, 2 ** n)
                + B4 * Fraction(1, 4 ** n)
                + C4 * Fraction(1, 8 ** n)
                + D4 * Fraction(1, 16 ** n))
        actual = eps[n]
        r = actual - pred
        print(f"  n = {n}:")
        print(f"    pred   = {pred}    (~ {float(pred):+.6e})")
        print(f"    actual = {actual}    (~ {float(actual):+.6e})")
        print(f"    r      = {r}    (~ {float(r):+.6e})")
        if r == 0:
            print("    *** EXACT ***")
        else:
            print(f"    |r/pred| = {float(abs(r)/abs(pred)):.4e}")
        print()

    # ----------------------------------------------------------------------- #
    # B.4 (optional): order-4 empirical recursion fit, T_4 spectrum from data #
    #      Need 4 windows of length 4 -> needs eps_1..eps_8.  We only have    #
    #      eps_1..eps_6 -> at most 2 windows -> under-determined.  Skip and   #
    #      document.                                                          #
    # ----------------------------------------------------------------------- #
    print("-" * 78)
    print("B.4  Order-4 empirical recursion fit -> SKIPPED")
    print("     (need eps_1..eps_8; only eps_1..eps_6 available)")
    print("-" * 78)
    print()


if __name__ == "__main__":
    main()
