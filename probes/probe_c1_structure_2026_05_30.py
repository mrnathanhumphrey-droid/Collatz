"""
probe_c1_structure_2026_05_30.py

Hunt for a V/T-style closed form at c(1). The exact c(1) is rational with ~10^34-digit
denominator because of the 2^{136} = 2^{ord_2(q^2)} period factor from one level of deepening.
Decompose c(1) = (2^{136} X + Y) / (2^{136} X' + Y') and inspect:
- X/X' = the leading-order rational (should match c(1) to O(2^{-72}) ~ 10^{-22}).
- The factorization of X, X', Y, Y' and their numerators/denominators after reduction.

If X/X' has small denominator with recognizable structure (like c(0)'s 19/127 = 19/(2^7-1)),
that's the V/T-analog at depth 1. If it's also huge, c(1) is genuinely non-clean.
"""
from fractions import Fraction
import math, sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8

def chi(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1

inv2_q = pow(2, q - 2, q)
pow_inv2 = [pow(inv2_q, a, q) for a in range(ord_2)]

def W(a_res):
    if a_res == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - a_res), 2**ord_2 - 1)

# N(s), T(s) for all s in Z/17
N = {}; T = {}
for s in range(q):
    Ns = Fraction(0); Ts = Fraction(0)
    for ar in range(ord_2):
        for br in range(ord_2):
            v = (s + pow_inv2[ar] - pow_inv2[br]) % q
            w = W(ar) * W(br)
            if v != 0:
                Ts += w
                Ns += w * chi(v, q)
    N[s] = Ns; T[s] = Ts

print("N(s), T(s) for s in Z/17:")
print(f"  {'s':>3} {'N(s)':>22} {'T(s)':>22}")
for s in range(q):
    print(f"  {s:>3} {str(N[s]):>22} {str(T[s]):>22}")

# Verify N(-s) = N(s) (i.e., N(17-s) = N(s)) since chi(-1) = +1 for q=17
print("\nSymmetry check N(s) == N(-s)?")
for s in range(1, 9):
    eq = N[s] == N[(q - s) % q]
    print(f"  N({s}) {'==' if eq else '!='} N({(q-s)%q}): {eq}")

# X = N(0) + 2*sum_{r=1}^{8} 2^{-8r} * N(2r mod 17)
# X' = T(0) + 2*sum_{r=1}^{8} 2^{-8r} * T(2r mod 17)
X = N[0]; Xp = T[0]
for r in range(1, 9):
    fac = Fraction(1, 2**(8 * r))
    X += 2 * fac * N[(2 * r) % q]
    Xp += 2 * fac * T[(2 * r) % q]

print(f"\nX  = {X}")
print(f"  reduced: num={X.numerator}, den={X.denominator}, factorization of num: {math.factorial(0)} ...")  # placeholder
print(f"X' = {Xp}")
print(f"\nX/X' = {X/Xp}")
print(f"  num digits: {len(str((X/Xp).numerator))}, den digits: {len(str((X/Xp).denominator))}")
print(f"  decimal: {float(X/Xp):.15f}")

# compare to c(1) measured from FFT
print(f"\nFFT c(1) measured = 0.153178230055")
print(f"X/X' - c(1)_fft   = {float(X/Xp) - 0.153178230055:.2e}")

# Y, Y' = the "subleading" combination with positive powers of 2^{8r}
Y = N[0]; Yp = T[0]
for r in range(1, 9):
    fac = 2**(8 * r)
    Y += 2 * fac * N[(2 * r) % q]
    Yp += 2 * fac * T[(2 * r) % q]
print(f"\nY  = {Y}")
print(f"Y' = {Yp}")
print(f"Y/Y' = {float(Y/Yp):.15f}")

# Exact c(1) = (2^136 X + Y)/(2^136 X' + Y')
two136 = 2**136
num = two136 * X + Y
den = two136 * Xp + Yp
c1 = num / den
print(f"\nc(1) exact = {float(c1):.15f}")
print(f"  num digits: {len(str(c1.numerator))}, den digits: {len(str(c1.denominator))}")
print(f"  match FFT to: {float(c1) - 0.153178230055:.2e}")

# Factor X's denominator and look for structure
print(f"\nStructure of X/X':")
xx = X/Xp
print(f"  X/X' = {xx.numerator} / {xx.denominator}")
# try to factor the denominator (small primes)
def factorize(n, lim=100000):
    facs = []
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 127, 257, 65537]:
        while n % p == 0:
            facs.append(p); n //= p
    if n > 1: facs.append(n)
    return facs

print(f"  num factorization: {factorize(xx.numerator)}")
print(f"  den factorization: {factorize(xx.denominator)}")
print(f"\nAlso c(0) = 19/127 for comparison:")
print(f"  19 prime; 127 = 2^7 - 1 (Mersenne)")
