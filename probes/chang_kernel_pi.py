"""
chang_kernel_pi.py — Build Chang/Quadrium 2603.11066v6 Definition C.5
cylinder-averaged mod-64 kernel P, solve for stationary pi exactly,
compare to our Result 40 empirical depletion/enhancement.

Door 2 of the Chang connection arc.
"""
import numpy as np
from fractions import Fraction
import sys

sys.stdout.reconfigure(encoding="utf-8")


def collatz_odd_step(n):
    """One odd-to-odd Syracuse step: n -> (3n+1)/2^v2(3n+1)."""
    m = 3 * n + 1
    while m % 2 == 0:
        m //= 2
    return m


def build_chang_kernel(use_fractions=True):
    """Build P per Definition C.5: 32x32, depth-13 lifts (2^7 lifts per residue)."""
    odd_residues = list(range(1, 64, 2))  # 32 residues
    idx = {r: i for i, r in enumerate(odd_residues)}
    N_LIFTS = 128  # 2^7

    if use_fractions:
        P = [[Fraction(0) for _ in range(32)] for _ in range(32)]
    else:
        P = np.zeros((32, 32))

    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(N_LIFTS):
            n = r + 64 * k
            next_n = collatz_odd_step(n)
            next_r = next_n % 64
            counts[idx[next_r]] += 1
        for j in range(32):
            if use_fractions:
                P[i][j] = Fraction(counts[j], N_LIFTS)
            else:
                P[i][j] = counts[j] / N_LIFTS

    return P, odd_residues, idx


def stationary_exact(P_frac):
    """Solve pi P = pi exactly using fractions. Returns list of 32 Fractions."""
    # We solve (P^T - I) pi = 0 with sum pi = 1
    # Using Gauss elimination with rationals
    n = 32
    A = [[P_frac[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    # Replace last row with sum constraint
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)

    # Gauss elimination
    for col in range(n):
        # Find pivot
        pivot = -1
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            raise ValueError("Singular matrix")
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        # Normalize
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] = A[col][j] / piv
        b[col] = b[col] / piv
        # Eliminate
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] = A[row][j] - factor * A[col][j]
                b[row] = b[row] - factor * b[col]

    return b


def main():
    print("# Chang/Quadrium kernel: building P (Definition C.5)...")
    P_frac, odd_residues, idx = build_chang_kernel(use_fractions=True)
    print("# Built. Verifying row sums = 1.")
    for i in range(32):
        s = sum(P_frac[i])
        assert s == 1, f"Row {i} sum = {s}"
    print("# Row sums OK. Solving for stationary pi (exact rationals).")

    pi = stationary_exact(P_frac)
    pi_float = [float(p) for p in pi]

    print(f"\n# Verification: pi sums to {sum(pi)} (should be 1)")
    print(f"# All non-negative: {all(p >= 0 for p in pi)}")

    # Mass on I_2 = {7, 27, 31, 59, 63}
    I2 = [7, 27, 31, 59, 63]
    mass_I2 = sum(pi[idx[r]] for r in I2)
    chang_reports = Fraction(10121, 65280)
    print(f"\n# pi(I_2) = {mass_I2} = {float(mass_I2):.7f}")
    print(f"# Chang reports 10121/65280 = {float(chang_reports):.7f}")
    print(f"# Exact match: {mass_I2 == chang_reports}")

    # Print full distribution
    print("\n# Stationary distribution mod 64 (32 odd residues)")
    print(f"# Uniform baseline = 1/32 = {1/32:.5f}")
    print(f"\n#   {'r':>3}  {'pi(r)':>10}  {'pi/unif':>9}  {'I_2':>5}")
    for i, r in enumerate(odd_residues):
        in_i2 = "yes" if r in I2 else ""
        ratio = pi_float[i] * 32
        print(f"#   {r:>3}  {pi_float[i]:>10.6f}  {ratio:>9.4f}  {in_i2:>5}")

    # Project to mod 32 to compare with Result 40
    # Each odd r mod 32 corresponds to {r, r+32} in mod 64
    print("\n# Projection to mod 32 (matches Result 40 observable)")
    print(f"# Uniform mod-32 baseline = 1/16 = {1/16:.5f}")
    print(f"\n#   {'r mod 32':>9}  {'pi_32':>10}  {'pi_32/unif':>11}")
    pi_mod32 = {}
    for r32 in range(1, 32, 2):
        # Lift to mod 64: r32 and r32+32
        val = pi[idx[r32]] + pi[idx[r32 + 32]]
        pi_mod32[r32] = val
        ratio = float(val) * 16
        print(f"#   {r32:>9}  {float(val):>10.6f}  {ratio:>11.4f}")

    # Result 40 comparison
    print("\n# === Result 40 comparison ===")
    print("# Result 40 (survivor-conditioned at large t):")
    print("#   r=21 mod 32: DEPLETED (argmin from t=30 to t=60)")
    print("#   r=5  mod 32: ENHANCED (argmax from t=40)")
    print()
    print("# Chang stationary pi (un-conditioned cylinder-averaged):")
    for r32 in [21, 5, 13, 1, 11, 31]:
        val = float(pi_mod32[r32])
        ratio = val * 16
        flag = ""
        if ratio > 1.05:
            flag = "ENHANCED"
        elif ratio < 0.95:
            flag = "DEPLETED"
        print(f"#   r={r32:>2}: pi = {val:.6f}, pi/unif = {ratio:.4f}  {flag}")

    # Save to CSV for record
    with open("experiments_output/chang_pi.csv", "w") as f:
        f.write("r,pi_exact,pi_float,pi_over_uniform,in_I2\n")
        for i, r in enumerate(odd_residues):
            in_i2 = 1 if r in I2 else 0
            f.write(f"{r},{pi[i]},{pi_float[i]:.10f},{pi_float[i]*32:.6f},{in_i2}\n")
    print("\n[save] experiments_output/chang_pi.csv")


if __name__ == "__main__":
    main()
