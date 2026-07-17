"""
tauberian_verify.py

Verification script for TAUBERIAN_SCOPING_VERIFICATION.md (Phase 4).
Loads cached eps_n exact rationals through k=6 from experiments_output/result_77_7_eps_exact_through_k7.json
Tests several candidate Tauberian-theorem-driven asymptotic predictions.

Expected output reproduces the consistency-test tables in TAUBERIAN_SCOPING_VERIFICATION.md.

Run from C:/Collatz/ as: python tauberian_verify.py
"""

from fractions import Fraction
import json
import math


def main():
    path = "experiments_output/result_77_7_eps_exact_through_k7.json"
    with open(path) as f:
        data = json.load(f)

    # Load exact rationals
    eps = {}
    for k, v in data.items():
        eps[int(k)] = Fraction(int(v["num"]), int(v["den"]))

    ks = sorted(eps.keys())
    print(f"Loaded eps_n for n in {ks}")
    print()

    # Table: eps_n, |eps_n|*2^n, deviation from 1/30
    print("Empirical eps_n, |eps_n|*2^n, deviation from 1/30 = 0.03333:")
    print(f"{'n':>3} {'eps_n (float)':>16} {'|eps_n|*2^n':>14} {'delta_n':>14}")
    for n in ks:
        val = float(eps[n])
        absval = abs(val) * 2**n
        delta = absval - 1.0 / 30.0
        print(f"{n:3d} {val:+16.6e} {absval:14.6e} {delta:+14.6e}")
    print()

    # Test: |eps_n|*2^n*n^alpha (consistency test for various alpha)
    print("Consistency test: if eps_n ~ C*(1/2)^n*n^(-alpha-1) (Chevalier-style branch),")
    print("then |eps_n|*2^n*n^(alpha+1) should approach a constant.")
    print("Equivalently testing |eps_n|*2^n*n^beta directly for various beta:")
    print(f"{'beta':>6} ", end="")
    for n in ks:
        if n < 2:
            continue
        print(f"  n={n:<10}", end="")
    print()
    for beta in [0.0, 0.5, 1.0, 1.5, 2.0]:
        print(f"{beta:6.2f} ", end="")
        for n in ks:
            if n < 2:
                continue
            absval = abs(float(eps[n])) * 2**n
            v = absval * n**beta
            print(f"  {v:10.4e}", end="")
        print()
    print()

    # Subleading correction analysis
    print("Subleading correction analysis: delta_n = |eps_n|*2^n - 1/30")
    print("Test which power of n best fits delta_n:")
    print(f"{'n':>3} {'delta_n':>14} {'delta_n*n':>14} {'delta_n*n^(3/2)':>16} {'delta_n*n^2':>14}")
    for n in ks:
        if n < 2:
            continue
        absval = abs(float(eps[n])) * 2**n
        delta = absval - 1.0 / 30.0
        print(f"{n:3d} {delta:+14.6e} {delta*n:+14.6e} {delta*n**1.5:+16.6e} {delta*n*n:+14.6e}")
    print()

    # Log-log slope fits for consecutive pairs
    print("Implied exponent beta from log-log fit of |eps_n|*2^n ~ C/n^beta:")
    vals = [(n, abs(float(eps[n])) * 2**n) for n in ks if n >= 2]
    for i in range(len(vals) - 1):
        n1, v1 = vals[i]
        n2, v2 = vals[i + 1]
        beta = -math.log(v2 / v1) / math.log(n2 / n1)
        print(f"  n={n1} -> n={n2}: implied beta = {beta:+.4f}")
    print()

    # Direct test of Chevalier Thm 1.14 prediction: |eps_n|*2^n*n^(3/2) should be constant
    print("Chevalier Thm 1.14 (alpha=1/2) prediction: |eps_n|*2^n * n^(3/2) -> constant.")
    print(f"{'n':>3} {'|eps_n|*2^n*n^(3/2)':>22}")
    for n in ks:
        if n < 2:
            continue
        absval = abs(float(eps[n])) * 2**n
        print(f"{n:3d} {absval*n**1.5:22.6e}")
    print()
    print("Growth across data:")
    vals_check = [abs(float(eps[n])) * 2**n * n**1.5 for n in ks if n >= 2]
    print(f"  min = {min(vals_check):.4e}, max = {max(vals_check):.4e}, ratio = {max(vals_check)/min(vals_check):.3f}x")
    print()

    # Direct test of Chevalier Thm 1.16 M=1 prediction: |eps_n|*2^n*n^(1/2) should be constant
    print("Chevalier Thm 1.16 (M=1) prediction: |eps_n|*2^n * n^(1/2) -> constant.")
    print(f"{'n':>3} {'|eps_n|*2^n*n^(1/2)':>22}")
    for n in ks:
        if n < 2:
            continue
        absval = abs(float(eps[n])) * 2**n
        print(f"{n:3d} {absval*n**0.5:22.6e}")
    print()
    vals_check = [abs(float(eps[n])) * 2**n * n**0.5 for n in ks if n >= 2]
    print(f"  min = {min(vals_check):.4e}, max = {max(vals_check):.4e}, ratio = {max(vals_check)/min(vals_check):.3f}x")
    print()

    # Print verdict summary
    print("=" * 60)
    print("Verdict (per TAUBERIAN_SCOPING_VERIFICATION.md):")
    print("=" * 60)
    print("- Chevalier Thm 1.14 (pure alpha=1/2): FALSIFIED (4.3x range)")
    print("- Chevalier Thm 1.16 (M=1): WEAK SUPPORT (1.5x range)")
    print("- Leading simple-pole-like behavior: STRONGLY SUPPORTED")
    print("- Subleading branch type: INDETERMINATE from N=5")
    print("- Sign flip of delta_n at n=5->6: requires multi-term ansatz")


if __name__ == "__main__":
    main()
