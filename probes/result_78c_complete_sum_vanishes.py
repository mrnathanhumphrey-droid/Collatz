"""
result_78c_complete_sum_vanishes.py — KEY TEST: does the COMPLETE sum vanish?

Hypothesis (from Cochrane Theorem 2 + D = 0 analysis):
  S_complete(r, ℓ, ε, m) := Σ_{u=0}^{3^{r+1}-1} e_{3^{r+1}}(c_{ℓ,ε} · 4^u - 9mu) = 0 (or near 0)

If TRUE: the polynomial identification + Cochrane Theorem 2 predicts complete vanishing.
Our INCOMPLETE sum (length 3^{r-1}) is then a partial sum of a fully-cancelling sum.
Partial-sum bounds (Pólya-Vinogradov, Erdős-Turán) give O(√N · log N) saving.

This would CLOSE the rate-1/2 question via a different route: partial-sum analysis instead
of direct Cochrane Theorem 2.
"""
import sys
import os
import math
import cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def compute_complete_sum(r, ell_phase, eps, m):
    """S_complete = Σ_{u=0}^{3^{r+1}-1} e_{3^{r+1}}(c · 4^u - 9mu)."""
    N_phase = 3**(r+1)
    omega_r = 1 + 3**r
    if pow(omega_r, 3, N_phase) != 1:
        omega_r = 1 + 2 * (3**r)
    c_eps = pow(2, eps, N_phase)
    c_ell_eps = (c_eps * pow(omega_r, ell_phase, N_phase)) % N_phase

    total = complex(0, 0)
    pow_4 = 1
    for u in range(N_phase):
        phase_int = (c_ell_eps * pow_4 - 9 * m * u) % N_phase
        total += cmath.exp(2j * cmath.pi * phase_int / N_phase)
        pow_4 = (pow_4 * 4) % N_phase
    return total


def compute_partial_sum(r, ell_phase, eps, m, length):
    """S_partial = Σ_{u=0}^{length-1} e_{3^{r+1}}(c · 4^u - 9mu)."""
    N_phase = 3**(r+1)
    omega_r = 1 + 3**r
    if pow(omega_r, 3, N_phase) != 1:
        omega_r = 1 + 2 * (3**r)
    c_eps = pow(2, eps, N_phase)
    c_ell_eps = (c_eps * pow(omega_r, ell_phase, N_phase)) % N_phase

    total = complex(0, 0)
    pow_4 = 1
    for u in range(length):
        phase_int = (c_ell_eps * pow_4 - 9 * m * u) % N_phase
        total += cmath.exp(2j * cmath.pi * phase_int / N_phase)
        pow_4 = (pow_4 * 4) % N_phase
    return total


def main():
    print("# R78c: TEST — does the complete sum vanish?")
    print()

    print("# Test 1: |S_complete| at varying r, fixed (ℓ, ε, m)")
    print()
    print(f"  {'r':>3}  {'modulus':>10}  {'|S_complete|':>15}  {'|S_complete|/3^{r+1}':>22}  {'verdict':>8}")
    for r in [2, 3, 4]:
        N_phase = 3**(r+1)
        for ell, eps, m in [(0, 0, 0), (0, 0, 1), (1, 1, 5), (2, 0, 10)]:
            S = compute_complete_sum(r, ell, eps, m)
            verdict = "≈ 0" if abs(S) < 1e-6 else f"NOT 0"
            print(f"  {r:>3}  {N_phase:>10}  {abs(S):>15.6e}  {abs(S)/N_phase:>22.6e}  ({ell},{eps},{m}): {verdict}")
    print()

    # Test 2: comprehensive — does S_complete = 0 for ALL (ℓ, ε, m)?
    print("# Test 2: Comprehensive scan — is S_complete ≡ 0 over all (ℓ, ε, m)?")
    for r in [2, 3]:
        N_phase = 3**(r+1)
        N_r = 2 * 3**(r-1)
        max_abs = 0.0
        any_nonzero = False
        for ell in range(3):
            for eps in [0, 1]:
                for m in range(N_r):
                    S = compute_complete_sum(r, ell, eps, m)
                    if abs(S) > max_abs:
                        max_abs = abs(S)
                    if abs(S) > 1e-6:
                        any_nonzero = True
        print(f"  r={r}: max |S_complete| over all (ℓ, ε, m) = {max_abs:.6e}, "
              f"any non-zero? {any_nonzero}")
    print()

    # Test 3: partial sum sizes vs trivial bound
    print("# Test 3: |S_partial(length=3^{r-1})| vs |S_complete| (= 0?)")
    print()
    print(f"  Validates: Kalafatelis's incomplete sum is a PARTIAL SUM of a fully-cancelling complete sum")
    print()
    for r in [2, 3, 4, 5]:
        N_phase = 3**(r+1)
        N_partial = 3**(r-1)
        partial_sums = []
        for ell in range(3):
            for eps in [0, 1]:
                for m in range(2 * N_partial):
                    S_p = compute_partial_sum(r, ell, eps, m, N_partial)
                    partial_sums.append(abs(S_p))
        partial_arr = np.array(partial_sums)
        sqrt_N = math.sqrt(N_partial)
        print(f"  r={r}: max |S_partial| = {partial_arr.max():.4f}, "
              f"mean = {partial_arr.mean():.4f}, √N = {sqrt_N:.4f}")
        print(f"        max/√N = {partial_arr.max()/sqrt_N:.4f}, mean/√N = {partial_arr.mean()/sqrt_N:.4f}")
    print()

    # Save
    out = os.path.join(OUTDIR, "r78c_complete_partial.csv")
    with open(out, 'w', encoding='utf-8') as f:
        f.write("r,N_partial,max_S_partial,mean_S_partial,max_S_partial_over_sqrt_N\n")
        for r in [2, 3, 4, 5]:
            N_partial = 3**(r-1)
            partial_sums = []
            for ell in range(3):
                for eps in [0, 1]:
                    for m in range(2 * N_partial):
                        partial_sums.append(abs(compute_partial_sum(r, ell, eps, m, N_partial)))
            partial_arr = np.array(partial_sums)
            sqrt_N = math.sqrt(N_partial)
            f.write(f"{r},{N_partial},{partial_arr.max():.6f},{partial_arr.mean():.6f},"
                    f"{partial_arr.max()/sqrt_N:.6f}\n")
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
