"""
probe_step3_cofactor_2026_05_30.py

Step 3: 2675591024983 - 1 = 2·3·127·3511274311.
The 127 hit is suspicious. Now: factor 3511274311 fully, and check whether
its prime factors relate to ord_{q^k}(2) for q=17, or to other depth-related ordinals.
"""
from __future__ import annotations
import sys
from sympy import factorint, isprime
sys.stdout.reconfigure(encoding="utf-8")

q = 17
cof = 3511274311

print(f"Factoring 3511274311 = ", end="")
f = factorint(cof)
print(f)
print(f"  Prime? {isprime(cof)}")

# All factors
prime_check = sorted(f.keys())
print(f"\nPrime factors with structure check:")
print(f"{'p':>15} {'mod 17':>8} {'mod 8':>6} {'p-1 factored':>40} {'ord_p(2)':>15}")
for p in prime_check:
    pm1 = p - 1
    fpm1 = factorint(pm1)
    # ord_p(2) — try divisors of p-1
    ord2 = None
    if p > 2:
        from sympy import divisors
        for d in sorted(divisors(pm1)):
            if pow(2, d, p) == 1:
                ord2 = d; break
    print(f"{p:>15} {p % 17:>8} {p % 8:>6} {str(fpm1):>40} {str(ord2):>15}")

# Compare with ord_{q^k}(2) for k=1..5
print(f"\n=== ord_{{17^k}}(2) reference values ===")
for k in range(1, 6):
    mod = q ** k
    # ord_{q^k}(2) = ord_q(2) * q^{k-1} typically = 8 * 17^{k-1}
    o = 1; x = 2 % mod
    while x != 1 and o < 10**7:
        o += 1; x = (x * 2) % mod
    print(f"  ord_{{17^{k}}}(2) = {o} = {factorint(o)}")

# Does cofactor's factorization include any of these?
print(f"\n=== Cross-check: does 3511274311 divide 2^k - 1 for any structural k? ===")
import math
log2_p = math.log2(cof)
print(f"  log2(3511274311) = {log2_p:.4f}")
# ord_{cof}(2) is the smallest k with 2^k ≡ 1 mod cof.
# We already computed it above.

# Also: is cof a factor of 2^k - 1 for k = ord_{q^j}(2) for any j?
print(f"\n=== Does p | 2^{{8·17^j}} - 1 for j=0,1,2,3? ===")
for p in prime_check:
    if p == 2: continue
    for j in range(0, 4):
        k = 8 * (17**j)
        residue = pow(2, k, p)
        if residue == 1:
            print(f"  {p}: divides 2^{k} - 1 (i.e., 2^{{8·17^{j}}}-1)")
            break
    else:
        print(f"  {p}: does NOT divide 2^{{8·17^j}}-1 for j=0..3")

# And: is cof related to (2^k - 1)/(some divisor) for any k?
# Quick: does p ≡ 1 mod 17 (so 17 | p-1)?
print(f"\n=== Does p ≡ 1 mod 17 (i.e., 17 | p-1)? ===")
for p in prime_check:
    if (p - 1) % 17 == 0:
        print(f"  {p}: YES, p-1 contains 17 as factor")
    else:
        print(f"  {p}: NO ({p-1} mod 17 = {(p-1) % 17})")

# And the 13-digit prime 2675591024983 itself
print(f"\n=== Full context: 2675591024983 - 1 factorization ===")
p = 2675591024983
fpm1 = factorint(p - 1)
print(f"  p - 1 = {p-1} = {fpm1}")
# Compare to {2, 3, 127} pattern: 2·3·127 = 762. 762 * X = p - 1?
prod = 2 * 3 * 127
print(f"  2·3·127 = {prod}, (p-1)/{prod} = {(p-1)//prod} = {factorint((p-1)//prod)}")

# Check: do the 3 huge denominator factors have the SAME prime pattern in p-1?
print(f"\n=== Triple-check: 19841 - 1, 2675591024983 - 1, 1989596835661185763 - 1 ===")
big_dens = [
    ("19841", 19841),
    ("2675591024983", 2675591024983),
    ("1989596835661185763", 1989596835661185763),
]
for name, n in big_dens:
    f = factorint(n - 1)
    print(f"  {name} - 1: {f}")
    # 17 in factorization?
    if 17 in f:
        print(f"    ★ contains 17: 17^{f[17]}")
    if 127 in f:
        print(f"    ★ contains 127: 127^{f[127]}")
    # Mersenne-like?
    base = 1
    for p in sorted(f):
        base *= p**f[p]
    # number of distinct prime factors
    print(f"    #distinct primes: {len(f)};  largest prime: {max(f)}")

# The 29-digit num prime
print(f"\n=== huge1 - 1 ===")
p = 19777954930411325600042121089
f = factorint(p - 1)
print(f"  factorization: {f}")
if 17 in f:
    print(f"    ★ contains 17: 17^{f[17]}")
if 127 in f:
    print(f"    ★ contains 127: 127^{f[127]}")
