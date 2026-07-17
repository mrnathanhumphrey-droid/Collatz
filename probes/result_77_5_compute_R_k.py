"""
result_77_5_compute_R_k.py — Stage 1 of R77.5

Compute the inter-level lift residual operator R_k := pi_{k+1} - T(pi_k)
explicitly over Q for k = 1..5 using exact rationals.

Lift map T_{k -> k+1}: for each coprime r' in Z/3^{k+1}, set
    T(pi_k)(r') = pi_k[r' mod 3^k] / 3.

This is well-defined because every coprime r' in Z/3^{k+1} has r' mod 3^k
coprime in Z/3^k (the 3 lifts of each coprime r in Z/3^k are r, r+3^k,
r+2*3^k, all coprime in Z/3^{k+1}).

Outputs:
  C:\\Collatz\\result_77_5_R_k_norms.csv

Reuses build_markov_rational(k) and stationary_rational(K) from
push_to_k6_rate_analysis.py.
"""
import sys
import os
import time
import csv
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
OUT_CSV = r"C:\Collatz\result_77_5_R_k_norms.csv"


def build_markov_rational(k):
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


def pi_dict(k):
    """Return (pi as dict r -> Fraction, coprime_list) at level k."""
    K, coprime = build_markov_rational(k)
    pi_vec = stationary_rational(K)
    pi = {coprime[i]: pi_vec[i] for i in range(len(coprime))}
    return pi, coprime


def lift_pi(pi_k, k):
    """Apply T_{k -> k+1}: each r' in coprime(Z/3^{k+1}) gets pi_k[r' mod 3^k] / 3."""
    Nk = 3 ** k
    Nk1 = 3 ** (k + 1)
    third = Fraction(1, 3)
    lift = {}
    for rp in range(Nk1):
        if rp % 3 == 0:
            continue
        r = rp % Nk
        # since rp is coprime to 3 in Z/3^{k+1} iff rp not divisible by 3,
        # and r = rp mod 3^k inherits the same coprime property.
        lift[rp] = pi_k[r] * third
    return lift


def squared_l2_norm(vec_dict):
    """Return Σ v^2 over all entries (Fraction)."""
    s = Fraction(0)
    for v in vec_dict.values():
        s += v * v
    return s


def sum_entries(vec_dict):
    s = Fraction(0)
    for v in vec_dict.values():
        s += v
    return s


def main():
    print("# R77.5 Stage 1: compute R_k := pi_{k+1} - T(pi_k)")
    print()

    # Compute pi_k for k = 1..6
    pi_at = {}
    coprime_at = {}
    for k in [1, 2, 3, 4, 5, 6]:
        t0 = time.time()
        pi, cop = pi_dict(k)
        pi_at[k] = pi
        coprime_at[k] = cop
        print(f"  pi_{k}: {len(cop)} states, solved in {time.time()-t0:.2f}s")

    print()

    # Save Fraction R_k for Stage 2
    R = {}
    norms_sq = {}

    for k in [1, 2, 3, 4, 5]:
        pi_k = pi_at[k]
        pi_k1 = pi_at[k + 1]
        T_pi_k = lift_pi(pi_k, k)

        # Sanity: T(pi_k) and pi_{k+1} should both sum to 1
        s_T = sum_entries(T_pi_k)
        s_pi_k1 = sum_entries(pi_k1)
        assert s_T == 1, f"lift sum != 1 at k={k}: {s_T}"
        assert s_pi_k1 == 1, f"pi_{k+1} sum != 1 at k={k}: {s_pi_k1}"

        # R_k = pi_{k+1} - T(pi_k); same coprime support in Z/3^{k+1}
        R_k = {}
        for rp in pi_k1:
            R_k[rp] = pi_k1[rp] - T_pi_k[rp]

        # Sum check (should be 0)
        s_R = sum_entries(R_k)
        assert s_R == 0, f"R_{k} sum != 0: {s_R}"

        norm_sq = squared_l2_norm(R_k)
        R[k] = R_k
        norms_sq[k] = norm_sq

        print(f"  k={k}: ‖R_{k}‖²_2 = {norm_sq.numerator}/{norm_sq.denominator}  ≈ {float(norm_sq):.10e}")

    print()
    print("# Ratios ‖R_k‖² / ‖R_{k-1}‖²:")
    ratios = {}
    for k in [2, 3, 4, 5]:
        r = norms_sq[k] / norms_sq[k - 1]
        ratios[k] = r
        print(f"  k={k}: ratio = {float(r):.6f}  (= {r.numerator}/{r.denominator})")

    print()
    print("# Stationary R_k norm scaling: ‖R_k‖² · 3^k (testing 1/2 rate ⇔ ratio 1/4):")
    for k in [1, 2, 3, 4, 5]:
        scaled = norms_sq[k] * (3 ** k)
        print(f"  k={k}: ‖R_{k}‖² * 3^{k} = {float(scaled):.10e}")

    print()
    # CSV output
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["k", "N_k", "N_kplus1", "norm_sq_num", "norm_sq_den",
                    "norm_sq_decimal", "norm_sq_times_3pow_k", "ratio_to_prev"])
        for k in [1, 2, 3, 4, 5]:
            n_sq = norms_sq[k]
            scaled = n_sq * (3 ** k)
            ratio = float(ratios[k]) if k in ratios else ""
            w.writerow([
                k,
                len(coprime_at[k]),
                len(coprime_at[k + 1]),
                str(n_sq.numerator),
                str(n_sq.denominator),
                f"{float(n_sq):.15e}",
                f"{float(scaled):.15e}",
                f"{ratio:.10f}" if ratio != "" else "",
            ])
    print(f"[save] {OUT_CSV}")

    # Persist R as a pickle-free txt for Stage 2 reuse
    # Actually, simpler: re-import this module from Stage 2.

    return R, norms_sq, pi_at, coprime_at


if __name__ == "__main__":
    main()
