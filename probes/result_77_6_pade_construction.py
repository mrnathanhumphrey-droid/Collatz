"""
result_77_6_pade_construction.py
================================
R77.6 — Generating function analysis of eps_n via Pade approximants.

Stage 1: rebuild eps_1..eps_6 as exact Fractions using the project's
         level-k Markov chain (push_to_k6_rate_analysis.py logic).

Stage 2: construct [m/n] Pade approximants of f_tilde(z) = E(z)/z^2 over Q
         for various (m, n) pairs. Since we have 5 nontrivial coefficients
         (eps_2..eps_6), Pade is well-defined for m + n <= 4.

Stage 3: report pole locations (roots of Q) for each [m/n], distance to
         z = 2 (the rate-1/2 endpoint), and stability across (m, n).

Output:
  - C:/Collatz/experiments_output/result_77_6_pade_poles.csv
  - stdout summary
"""
from __future__ import annotations

import csv
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
OUT_CSV = os.path.join(OUTDIR, "result_77_6_pade_poles.csv")


# --------------------------------------------------------------------------- #
# Stage 1: exact Markov chain                                                 #
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


def get_eps_exact(max_k: int = 6) -> dict:
    """Return {k: Fraction(eps_k)} for k = 1..max_k."""
    target = Fraction(7, 15)
    X = {0: Fraction(1)}
    S = {}
    eps = {}
    for k in range(1, max_k + 1):
        t0 = time.time()
        K, coprime = build_markov_rational(k)
        pi_q = stationary_rational(K)
        sum_pi_sq = sum(p * p for p in pi_q)
        X[k] = Fraction(3 ** k) * sum_pi_sq
        S[k] = X[k] - X[k - 1]
        eps[k] = S[k] - target
        print(f"  k = {k:>2} ({len(coprime):>4} states):  elapsed {time.time() - t0:.1f}s "
              f"  eps_{k} ~ {float(eps[k]):+.8e}")
    return eps


# --------------------------------------------------------------------------- #
# Stage 2: Pade construction over Q                                           #
# --------------------------------------------------------------------------- #

def gauss_solve_Q(A, b):
    """Solve A x = b over Q. A: n x n list of Fractions; b: list of Fractions."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            return None  # singular
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


def pade_approximant(coeffs, m: int, n: int):
    """Compute [m/n] Pade for series f(z) = sum_{k>=0} c_k z^k.

    Returns (P_coeffs, Q_coeffs) where P has degree <= m and Q has degree <= n,
    Q(0) = 1, and f(z) Q(z) - P(z) = O(z^{m+n+1}).

    Returns None if the linear system is singular.
    """
    if len(coeffs) < m + n + 1:
        raise ValueError(f"Need at least {m + n + 1} coefficients; got {len(coeffs)}")
    c = coeffs

    if n == 0:
        # Trivial: Q = 1, P is just the truncation of f to degree m
        P = [c[k] for k in range(m + 1)]
        return P, [Fraction(1)]

    # Set up n x n linear system for q_1..q_n:
    #   sum_{j=1}^{n} q_j * c_{m+i-j} = -c_{m+i},  for i = 1..n
    A = []
    b = []
    for i in range(1, n + 1):
        row = []
        for j in range(1, n + 1):
            idx = m + i - j
            if idx < 0:
                row.append(Fraction(0))
            else:
                row.append(c[idx])
        A.append(row)
        b.append(-c[m + i])

    sol = gauss_solve_Q(A, b)
    if sol is None:
        return None
    Q = [Fraction(1)] + sol  # Q(0) = 1

    # Compute P_k = sum_{j=0}^{min(k, n)} q_j * c_{k - j}, for k = 0..m
    P = []
    for k in range(m + 1):
        pk = Fraction(0)
        for j in range(0, min(k, n) + 1):
            pk += Q[j] * c[k - j]
        P.append(pk)
    return P, Q


def poly_eval_float(coeffs, z: complex) -> complex:
    """Evaluate sum_j c_j z^j at z (float coefficients)."""
    return sum(float(c) * (z ** j) for j, c in enumerate(coeffs))


def find_roots_numpy(coeffs):
    """Find roots of polynomial sum_j c_j z^j using numpy. Returns list of complex.
    coeffs: list of Fraction; numpy.roots takes coefficients in DECREASING order
    of degree (highest first)."""
    cf = [float(c) for c in coeffs[::-1]]
    if cf[0] == 0:
        # Strip leading zeros
        while cf and cf[0] == 0:
            cf.pop(0)
    if len(cf) <= 1:
        return []
    return list(np.roots(cf))


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def main():
    print("=" * 78)
    print("R77.6 — Generating function analysis of eps_n via Pade approximants")
    print("=" * 78)
    print()

    print("Stage 1: rebuild exact eps_1..eps_6 from project's Markov chain")
    print("-" * 78)
    eps = get_eps_exact(max_k=6)
    print()

    # E(z) = sum_{n>=1} eps_n z^n. We work with f_tilde(z) = (E(z) - eps_1 z) / z^2
    # to drop the n=1 transient (per R77 brief). Coefficients of f_tilde are
    # eps_2, eps_3, eps_4, eps_5, eps_6 at degrees 0, 1, 2, 3, 4.
    coeffs = [eps[2], eps[3], eps[4], eps[5], eps[6]]
    print("Stage 2: working series f_tilde(z) = eps_2 + eps_3 z + eps_4 z^2 + eps_5 z^3 + eps_6 z^4")
    print("-" * 78)
    print("(n=1 transient excluded; z=2 still corresponds to rate-1/2 since")
    print(" f_tilde shares E(z)'s singularities away from z=0.)")
    print()
    for j, c in enumerate(coeffs):
        print(f"  c_{j} = eps_{j + 2} = {c}")
        print(f"        ~ {float(c):+.10e}")
    print()

    # We have 5 coefficients, so [m/n] Pade is well-defined for m + n <= 4.
    # Brief asks for various (m, n) — all those with m+n=5 require eps_7
    # (only feasible by extending k=7; brief notes 5-point analysis is diagnostic only).
    targets = [
        # (m, n) — annotated by intended role
        (1, 1, "lowest-order diagonal"),
        (2, 1, "near-diagonal m+n=3"),
        (1, 2, "near-diagonal m+n=3"),
        (3, 1, "[3/1] m+n=4"),
        (2, 2, "diagonal [2/2] m+n=4"),
        (1, 3, "[1/3] m+n=4"),
        (4, 0, "Taylor truncation (sanity: no poles expected)"),
        (0, 4, "[0/4] all-pole approximant"),
    ]

    results = []
    print("Stage 3: Pade approximants and pole locations")
    print("-" * 78)

    for m, n, role in targets:
        if m + n + 1 > len(coeffs):
            print(f"  [m={m}/n={n}]: SKIP — requires {m + n + 1} coefficients, have {len(coeffs)}")
            continue
        out = pade_approximant(coeffs, m, n)
        if out is None:
            print(f"  [m={m}/n={n}]: SINGULAR linear system")
            continue
        P, Q = out
        # Print P, Q
        print(f"  [{m}/{n}]  ({role})")
        P_str = " + ".join(f"({p})z^{j}" for j, p in enumerate(P) if p != 0) or "0"
        Q_str = " + ".join(f"({q})z^{j}" for j, q in enumerate(Q) if q != 0) or "0"
        print(f"    P(z) = {P_str}")
        print(f"    Q(z) = {Q_str}")

        # Floats for human reading
        P_f = [float(p) for p in P]
        Q_f = [float(q) for q in Q]
        P_str_f = " + ".join(f"{p:+.6e} z^{j}" for j, p in enumerate(P_f))
        Q_str_f = " + ".join(f"{q:+.6e} z^{j}" for j, q in enumerate(Q_f))
        print(f"    P(z) (float) = {P_str_f}")
        print(f"    Q(z) (float) = {Q_str_f}")

        # Roots of Q (poles of [m/n])
        roots = find_roots_numpy(Q)
        if not roots:
            print(f"    poles: none (n=0 approximant)")
        else:
            print(f"    poles (roots of Q):")
            for r in roots:
                d2 = abs(r - 2.0)
                print(f"      z = {r.real:+.6e} {r.imag:+.6e}j   |z| = {abs(r):.6e}   |z - 2| = {d2:.6e}")

        # CSV rows
        if roots:
            for j_r, r in enumerate(roots):
                results.append({
                    "m": m, "n": n, "role": role,
                    "pole_idx": j_r,
                    "pole_real": r.real, "pole_imag": r.imag,
                    "pole_abs": abs(r),
                    "dist_to_2": abs(r - 2.0),
                    "dist_to_neg2": abs(r + 2.0),
                })
        else:
            results.append({
                "m": m, "n": n, "role": role,
                "pole_idx": -1,
                "pole_real": np.nan, "pole_imag": np.nan,
                "pole_abs": np.nan,
                "dist_to_2": np.nan,
                "dist_to_neg2": np.nan,
            })
        print()

    # CSV write
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        fieldnames = ["m", "n", "role", "pole_idx", "pole_real", "pole_imag",
                       "pole_abs", "dist_to_2", "dist_to_neg2"]
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    print(f"[csv: {OUT_CSV}]")
    print()

    # ---- Stage 4: diagnostic summary ---- #
    print("=" * 78)
    print("Stage 4: diagnostic summary")
    print("=" * 78)

    # For each (m, n), report all poles' distance to z = 2.
    # Diagnostic signatures:
    # - Pure geometric (1/2)^n: single pole exactly at z = 2 stable across (m, n)
    # - Power-law: poles converge to z = 2 along real ray [2, infty)
    # - Logarithmic: poles cluster near z = 2 with distinctive pattern
    print()
    print("Closest-pole-to-z=2 across (m, n):")
    print(f"  {'(m,n)':>8}  {'closest pole':>40}  {'|z - 2|':>14}")
    for row in results:
        if row["pole_idx"] < 0:
            continue
        # gather all rows for this (m, n)
        pass

    by_mn = {}
    for r in results:
        if r["pole_idx"] < 0:
            continue
        key = (r["m"], r["n"])
        by_mn.setdefault(key, []).append(r)

    for key in sorted(by_mn.keys()):
        rs = by_mn[key]
        rs.sort(key=lambda x: x["dist_to_2"])
        closest = rs[0]
        z = complex(closest["pole_real"], closest["pole_imag"])
        print(f"  ({key[0]},{key[1]}):  z = {z.real:+.4e} {z.imag:+.4e}j   "
              f"|z - 2| = {closest['dist_to_2']:.4e}")
    print()

    # All real poles in window |z - 2| < 1
    print("All poles within |z - 2| < 1 (potential rate-1/2 signal):")
    print(f"  {'(m,n)':>8}  {'pole z':>40}  {'|z - 2|':>14}")
    for row in results:
        if row["pole_idx"] < 0:
            continue
        if row["dist_to_2"] < 1.0:
            z = complex(row["pole_real"], row["pole_imag"])
            print(f"  ({row['m']},{row['n']}):  {z.real:+.4e} {z.imag:+.4e}j   "
                  f"{row['dist_to_2']:.4e}")
    print()

    # Diagonal Pade pole evolution: (1,1) and (2,2) — show poles
    print("Diagonal Pade pole evolution [n/n]:")
    for diag_n in [1, 2]:
        rs = by_mn.get((diag_n, diag_n), [])
        if not rs:
            continue
        print(f"  [{diag_n}/{diag_n}] poles:")
        for r in rs:
            z = complex(r["pole_real"], r["pole_imag"])
            print(f"    {z.real:+.6e}{z.imag:+.6e}j   |z - 2| = {r['dist_to_2']:.4e}")
    print()


if __name__ == "__main__":
    main()
