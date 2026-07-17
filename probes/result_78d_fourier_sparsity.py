"""
result_78d_fourier_sparsity.py — characterize the Fourier sparsity structure of f(u) = e_q(4^u - 9mu).

Key finding from r78b/c: complete sum vanishes (Cochrane Theorem 2 D=0 prediction holds)
AND F̂(ξ) is supported on 3·Z/q (multiples of 3), with 6/9 zero spikes and 3/9 non-zero spikes.

This is because 4 has order 3^r in (Z/3^{r+1})*, so f(u) has period 3^r in u, and its Fourier
transform on Z/3^{r+1} is supported on (3^{r+1}/3^r)·Z/3^{r+1} = 3·Z/3^{r+1}.

For r-level analysis: F̂(3a) for a ∈ Z/3^r is the 3^r-Fourier transform of f.

This script:
1. Computes F̂(ξ) for r = 2, 3 explicitly
2. Identifies the sparsity pattern
3. Bounds the magnitudes

Then derives the implication for the Pólya-Vinogradov-style partial sum bound.
"""
import sys
import os
import math
import cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def F_hat(r, c, m, xi):
    """F̂(ξ) = Σ_{u=0}^{q-1} e_q(c·4^u - 9mu - ξu) where q = 3^{r+1}."""
    q = 3**(r+1)
    total = complex(0, 0)
    pow4 = 1
    for u in range(q):
        phase = ((c * pow4 - 9 * m * u - xi * u) % q) / q
        total += cmath.exp(2j * cmath.pi * phase)
        pow4 = (pow4 * 4) % q
    return total


def main():
    print("# R78d: Fourier sparsity of f(u) = e_q(4^u - 9mu)")
    print()

    for r in [2, 3]:
        q = 3**(r+1)
        period = 3**r  # period of 4 in (Z/q)*
        print(f"## r = {r}: q = 3^{r+1} = {q}, period of 4 = {period}")
        print()

        # F̂ supported on (q/period)·Z/q = 3·Z/q
        # That's 3^r values: 0, 3, 6, ..., 3·(3^r - 1)
        # In Z/q, these are {3a : a ∈ Z/3^r}.
        c, m = 1, 0
        F_vals = {}
        for xi in range(q):
            val = F_hat(r, c, m, xi)
            F_vals[xi] = val

        # Check sparsity
        nonzero_xi = [xi for xi in range(q) if abs(F_vals[xi]) > 1e-10]
        print(f"  Non-zero F̂(ξ) at ξ ∈ {nonzero_xi}")
        print(f"  Count: {len(nonzero_xi)} out of {q}")
        for xi in nonzero_xi:
            print(f"    ξ = {xi}: F̂(ξ) = {F_vals[xi]:.4f}, |F̂| = {abs(F_vals[xi]):.4f}")
        print()

        # Sparsity ratio
        sparsity = len(nonzero_xi) / q
        print(f"  Sparsity: {len(nonzero_xi)}/{q} = {sparsity:.4f}")
        print()

        # Now check structure: F̂(3a) for a ∈ Z/period
        print(f"  F̂ on 3·Z/q (= multiples of 3, restricted to period {period}):")
        for a in range(period):
            xi = 3 * a
            if xi < q:
                val = F_vals[xi]
                marker = "  *" if abs(val) > 1e-10 else ""
                print(f"    a={a:>3} (ξ=3a={xi:>3}): F̂ = {abs(val):.4f}{marker}")
        print()

        # Identify which a values have non-zero F̂(3a)
        nonzero_a = [a for a in range(period) if abs(F_vals[3*a]) > 1e-10]
        print(f"  Non-zero a values: {nonzero_a}")
        print(f"  These match... let's see: 3a in nonzero_xi means 3a ∈ {nonzero_xi}, "
              f"so a ∈ {[xi//3 for xi in nonzero_xi]}")
        print()

        # max |F̂| relative to size measures
        max_F = max(abs(F_vals[xi]) for xi in nonzero_xi)
        print(f"  max |F̂| = {max_F:.4f}")
        print(f"  q = {q}, √q = {math.sqrt(q):.4f}")
        print(f"  max |F̂| / q = {max_F/q:.4f}")
        print(f"  max |F̂| / period = {max_F/period:.4f}")
        print(f"  max |F̂| / √q = {max_F/math.sqrt(q):.4f}")
        print(f"  max |F̂| / √period = {max_F/math.sqrt(period):.4f}")
        print()

        # Implied alpha: |F̂| = q^α
        if max_F > 0:
            alpha = math.log(max_F) / math.log(q)
            print(f"  Implied α (|F̂| = q^α): α = {alpha:.4f}")
        print()
        print()

    # Final analysis
    print("# Implication for Pólya-Vinogradov bound on partial sum")
    print()
    print("partial = (1/q) · Σ_{ξ ∈ 3·Z/q} 1̂(ξ) · F̂(ξ)")
    print(f"  Number of non-zero F̂ at level r: ~ 3 (out of q = 3^{{r+1}})")
    print(f"  |F̂| up to ~ q^α with α empirically ≈ 0.84 at r=2, varying")
    print(f"  |1̂| up to N = 3^{{r-1}}")
    print()
    print("Trivial bound: partial ≤ (1/q) · #non-zero · max|1̂| · max|F̂|")
    print("              = (1/q) · 3 · 3^{r-1} · q^α")
    print("              = 3^r / q · q^α")
    print("              = 3^r / 3^{r+1} · q^α")
    print("              = q^α / 3")
    print()
    print("For α = 0.84: bound = q^{0.84}/3 ≈ q · q^{-0.16}/3")
    print("For α = 0.5 (Gauss-sum strength): bound = √q/3 — would give square-root cancellation")
    print()


if __name__ == "__main__":
    main()
