"""
result_77_7_extend_to_k7.py
===========================
R77.7 — Extend exact eps_n to k=7 and refresh the Pade analysis with 6
coefficients (m + n <= 5). Sharpens (G-power) vs (G-log) discrimination
beyond R77.6's N=5 budget.

Stage 1: Compute eps_7 via the level-7 Markov chain (N=1458 states).
         Cache all eps_1..eps_7 as exact Fractions (numer, denom) to JSON.

Stage 2: Construct [m/n] Pade approximants of f_tilde(z) = E(z)/z^2 over Q
         for all (m, n) with m + n <= 5.

Stage 3: Diagnostic — does the m+n=5 row sharpen the closest-pole-to-2
         distance? Does the convergence pattern firm up toward (G-power)
         or (G-log)?

This script is designed to be re-runnable: if cache exists, it skips
recomputation.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
EPS_CACHE = os.path.join(OUTDIR, "result_77_7_eps_exact_through_k7.json")
OUT_CSV = os.path.join(OUTDIR, "result_77_7_pade_poles.csv")


# --------------------------------------------------------------------------- #
# Markov chain (rational)                                                     #
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


def get_eps_exact(max_k: int = 7) -> dict:
    """Return {k: Fraction(eps_k)} for k = 1..max_k.

    Caches to JSON. Skips any level already in cache; only computes missing.
    """
    eps = {}
    if os.path.exists(EPS_CACHE):
        with open(EPS_CACHE, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        eps = {int(k): Fraction(int(d["num"]), int(d["den"])) for k, d in data.items()}
        if all(k in eps for k in range(1, max_k + 1)):
            print(f"  [cache hit] loaded eps_1..eps_{max_k} from {EPS_CACHE}")
            for k in range(1, max_k + 1):
                print(f"    eps_{k} ~ {float(eps[k]):+.10e}")
            return eps
        cached_levels = sorted(eps.keys())
        print(f"  [partial cache] have eps_{cached_levels} from {EPS_CACHE}")
        for k in cached_levels:
            print(f"    eps_{k} ~ {float(eps[k]):+.10e}  (cached)")

    target = Fraction(7, 15)
    # Reconstruct X_{k-1} = 1 + (k-1)*7/15 + sum(eps_1..eps_{k-1}) from cached eps
    for k in range(1, max_k + 1):
        if k in eps:
            continue
        t0 = time.time()
        K, coprime = build_markov_rational(k)
        t_build = time.time() - t0
        t0b = time.time()
        pi_q = stationary_rational(K)
        t_solve = time.time() - t0b
        sum_pi_sq = sum(p * p for p in pi_q)
        X_k = Fraction(3 ** k) * sum_pi_sq
        # X_{k-1} from cached eps_1..eps_{k-1}
        if k == 1:
            X_km1 = Fraction(1)
        else:
            sum_eps_prior = sum(eps[j] for j in range(1, k))
            X_km1 = Fraction(1) + Fraction(k - 1) * target + sum_eps_prior
        S_k = X_k - X_km1
        eps[k] = S_k - target
        print(f"  k = {k:>2} ({len(coprime):>5} states):  build {t_build:.1f}s   "
              f"solve {t_solve:.1f}s   eps_{k} ~ {float(eps[k]):+.10e}", flush=True)

        # Save partial cache after each level
        partial = {str(kk): {"num": str(eps[kk].numerator), "den": str(eps[kk].denominator)}
                   for kk in eps}
        with open(EPS_CACHE, "w", encoding="utf-8") as fh:
            json.dump(partial, fh, indent=2)

    return eps


# --------------------------------------------------------------------------- #
# Pade construction (over Q)                                                  #
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


def pade_approximant(coeffs, m: int, n: int):
    if len(coeffs) < m + n + 1:
        raise ValueError(f"Need at least {m + n + 1} coefficients; got {len(coeffs)}")
    c = coeffs

    if n == 0:
        P = [c[k] for k in range(m + 1)]
        return P, [Fraction(1)]

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
    Q = [Fraction(1)] + sol

    P = []
    for k in range(m + 1):
        pk = Fraction(0)
        for j in range(0, min(k, n) + 1):
            pk += Q[j] * c[k - j]
        P.append(pk)
    return P, Q


def find_roots_numpy(coeffs):
    cf = [float(c) for c in coeffs[::-1]]
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
    print("R77.7 — Extend exact eps_n to k=7; refresh Pade analysis (m+n <= 5)")
    print("=" * 78)
    print()

    print("Stage 1: rebuild exact eps_1..eps_7 (with cache)")
    print("-" * 78)
    eps = get_eps_exact(max_k=7)
    print()

    coeffs = [eps[2], eps[3], eps[4], eps[5], eps[6], eps[7]]
    print(f"Stage 2: f_tilde(z) coefficients (eps_2..eps_7), {len(coeffs)} terms")
    print("-" * 78)
    for j, c in enumerate(coeffs):
        ndig_n = len(str(abs(c.numerator)))
        ndig_d = len(str(abs(c.denominator)))
        print(f"  c_{j} = eps_{j + 2}:  num {ndig_n} digits, den {ndig_d} digits;  ~ {float(c):+.10e}")
    print()

    # All (m, n) with m + n <= 5 and m + n >= 2 (skip trivial m+n=0,1)
    targets = []
    for total in range(2, 6):
        for m in range(0, total + 1):
            n = total - m
            label = f"m+n={total}"
            if m == n:
                label += " diagonal"
            elif abs(m - n) == 1:
                label += " near-diag"
            targets.append((m, n, label))

    results = []
    print(f"Stage 3: {len(targets)} Pade approximants  (m+n in [2, 5])")
    print("-" * 78)

    for m, n, role in targets:
        if m + n + 1 > len(coeffs):
            print(f"  [{m}/{n}]: SKIP — needs {m + n + 1} coefficients, have {len(coeffs)}")
            continue
        out = pade_approximant(coeffs, m, n)
        if out is None:
            print(f"  [{m}/{n}] ({role}): SINGULAR linear system")
            continue
        P, Q = out
        roots = find_roots_numpy(Q)
        if roots:
            roots_sorted = sorted(roots, key=lambda r: abs(r - 2.0))
            closest = roots_sorted[0]
            d2 = abs(closest - 2.0)
            print(f"  [{m}/{n}] ({role}):  closest pole z = {closest.real:+.6f} {closest.imag:+.6f}j   "
                  f"|z - 2| = {d2:.4e}   ({len(roots)} pole(s) total)")
            for j_r, r in enumerate(roots_sorted):
                results.append({
                    "m": m, "n": n, "role": role,
                    "rank_by_dist_to_2": j_r,
                    "pole_real": r.real, "pole_imag": r.imag,
                    "pole_abs": abs(r),
                    "dist_to_2": abs(r - 2.0),
                })
        else:
            print(f"  [{m}/{n}] ({role}):  no poles (n=0 approximant)")

    print()
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        fieldnames = ["m", "n", "role", "rank_by_dist_to_2",
                       "pole_real", "pole_imag", "pole_abs", "dist_to_2"]
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    print(f"[csv: {OUT_CSV}]")
    print()

    # ---- Stage 4 diagnostic ---- #
    print("=" * 78)
    print("Stage 4: diagnostic")
    print("=" * 78)
    print()
    print("Closest-pole-to-z=2 by approximant order:")
    print(f"  {'(m,n)':>8}  {'closest pole':>40}  {'|z - 2|':>14}  {'role':>20}")
    by_mn = {}
    for r in results:
        if r["rank_by_dist_to_2"] != 0:
            continue
        by_mn[(r["m"], r["n"])] = r
    for key in sorted(by_mn.keys(), key=lambda x: (x[0] + x[1], x[0])):
        r = by_mn[key]
        z = complex(r["pole_real"], r["pole_imag"])
        print(f"  ({key[0]},{key[1]}):  z = {z.real:+.6e}{z.imag:+.6e}j   "
              f"{r['dist_to_2']:.4e}   {r['role']}")
    print()

    # Compare m+n=4 closest-to-2 vs m+n=5 closest-to-2: the SHARPENING
    print("Sharpening from m+n=4 to m+n=5:")
    closest_at_total = {}
    for r in results:
        if r["rank_by_dist_to_2"] != 0:
            continue
        total = r["m"] + r["n"]
        cur = closest_at_total.get(total)
        if cur is None or r["dist_to_2"] < cur:
            closest_at_total[total] = r["dist_to_2"]
    for total in sorted(closest_at_total.keys()):
        print(f"  m+n = {total}:  min |z - 2| over all (m,n) at this order = {closest_at_total[total]:.4e}")
    if 4 in closest_at_total and 5 in closest_at_total:
        ratio = closest_at_total[5] / closest_at_total[4]
        print(f"\n  Sharpening ratio (5 / 4) = {ratio:.4f}")
        if ratio < 0.4:
            print("  -> Strong tightening; consistent with power-law (G-power) with alpha large")
        elif ratio < 0.7:
            print("  -> Moderate tightening; consistent with branch-cut, type unresolved")
        elif ratio < 0.95:
            print("  -> Mild tightening; consistent with logarithmic (G-log) (slow convergence)")
        else:
            print("  -> No sharpening; suggests artifact or near-degenerate system at small N")


if __name__ == "__main__":
    main()
