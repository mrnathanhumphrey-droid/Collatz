"""
probe_delta1_prime_structure_2026_05_30.py

Check whether the prime factors of Δ_1's numerator / denominator have
relationships to:
  - (Z/17)* multiplicative structure: residue mod 17, coset in (Z/17)*/⟨2⟩
  - 8-power residues mod 17 (=±1) and mod 256 (since 2^8 ≡ 1 mod 17)
  - residue mod ord_{q^n}(2) for n=1,2,3: 8, 136, 2312
  - Mersenne-like form: p = (2^k ± 1)/d for small d, or p = a·2^k + 1
  - 2-adic, 17-adic structure
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord2 = 8

# Δ_1 factorization
nums = {
    "num: 2":                              2,
    "num: 19841":                          19841,
    "num: huge1 (29-digit)":               19777954930411325600042121089,
    "den: 5":                              5,
    "den: 13":                             13,
    "den: 127":                            127,
    "den: 2675591024983 (13-digit)":       2675591024983,
    "den: 1989596835661185763 (19-digit)": 1989596835661185763,
}

def legendre(x, q):
    x %= q
    return 0 if x == 0 else (1 if pow(x, (q-1)//2, q) == 1 else -1)

def in_pow2_subgroup(x, q):
    """Is x in <2> = QR mod q (since ord_q(2)=8 and (Z/17)*= cyclic of order 16, <2>=QR)."""
    return legendre(x, q) == 1

# Discrete log of x in <2>: solve 2^k ≡ x mod 17 for k in 0..7. None if x ∉ <2>.
discrete_log_2 = {}
x = 1
for k in range(ord2):
    discrete_log_2[x] = k
    x = (x * 2) % q

print(f"<2> mod 17 = {sorted(discrete_log_2.keys())} (QR mod 17)")
print(f"Discrete-log table: {discrete_log_2}\n")

print(f"{'prime/integer':40} {'mod 17':>8} {'in <2>?':>8} {'dlog2':>6} {'mod 8':>6} {'mod 16':>7} {'mod 127':>8} {'mod 255':>8} {'mod 256':>8}")
print("-" * 120)
for name, n in nums.items():
    r17 = n % q
    r8 = n % 8
    r16 = n % 16
    r127 = n % 127
    r255 = n % 255
    r256 = n % 256
    in_pow2 = in_pow2_subgroup(n, q) if r17 != 0 else None
    dlog = discrete_log_2.get(r17, "-")
    print(f"{name:40} {r17:>8} {str(in_pow2):>8} {str(dlog):>6} {r8:>6} {r16:>7} {r127:>8} {r255:>8} {r256:>8}")

# Are the primes 8-power residues mod 17? (i.e., p^? ≡ 1 mod 17)
# 8-power residues mod 17 = {x^8 mod 17} = {1, -1} since (Z/17)*= cyclic of order 16.
# So an "8-power residue" is x with x^2 ≡ 1 mod 17, i.e., x ≡ ±1 mod 17.
print(f"\n=== 8-power residue mod 17 check ===")
print(f"  An 8-power residue mod 17 satisfies x^2 ≡ 1 mod 17, i.e., x ≡ ±1 mod 17.")
for name, n in nums.items():
    r17 = n % q
    is_8pow = r17 in (1, 16)
    print(f"  {name:40} mod 17 = {r17:>3}  is 8-power residue? {is_8pow}")

# Are the primes congruent to ±1 mod 8 (would be unusual)?
print(f"\n=== Quartic residue / mod-8 check ===")
for name, n in nums.items():
    r8 = n % 8
    qr_mark = "" if r8 not in (1, 7) else "  ★ ≡ ±1 mod 8"
    print(f"  {name:40} mod 8 = {r8}{qr_mark}")

# Check: are the huge primes related to 2^k for some k?
print(f"\n=== Mersenne / Proth check on huge primes ===")
def near_pow2(n, name):
    import math
    log2 = math.log2(n)
    k = int(round(log2))
    nearby = 2**k
    diff = n - nearby
    rel = diff / nearby
    print(f"  {name}: 2^{k:.0f} = {nearby:.6e}, diff = {diff:+}, rel = {rel:+.6e}")
    # check small prime divisibility of (n-1) and (n+1)
    if n > 0:
        m = n - 1
        bits = 0
        while m % 2 == 0:
            m //= 2; bits += 1
        print(f"    n-1 = 2^{bits} · {m}  (factored low: {m})")
        # check if m has small prime factor
        from sympy import factorint
        f = factorint(m)
        # only show small factors
        small = {p: e for p, e in f.items() if p < 10**6}
        big = {p: e for p, e in f.items() if p >= 10**6}
        print(f"    n-1 small part: {small};  big primes: {list(big.keys())[:3]}{'...' if len(big) > 3 else ''}")

for name, n in [("19841", 19841),
                ("huge1 (29-digit)", 19777954930411325600042121089),
                ("2675591024983 (13-digit)", 2675591024983),
                ("1989596835661185763 (19-digit)", 1989596835661185763)]:
    print(f"\n{name}:")
    near_pow2(n, name)

# Check 17-adic structure: does any of these primes split / ramify / have special form in Q(√17) or Q(ζ_17)?
# A prime p splits in Q(√17) iff Legendre(17, p) = +1 iff p is QR mod 17 (for p ≠ 17).
print(f"\n=== Q(√17) splitting behavior of huge primes ===")
for name, n in [("huge1 (29-digit)", 19777954930411325600042121089),
                ("2675591024983 (13-digit)", 2675591024983),
                ("1989596835661185763 (19-digit)", 1989596835661185763)]:
    leg = legendre(n, 17)
    behav = {1: "splits", -1: "inert", 0: "ramifies"}[leg]
    print(f"  {name}: Legendre(n, 17) = {leg} → {behav} in Q(√17)")

# Check related: are the primes related to factorizations of (2^k - 1)/d for k ~ ord computations?
# 2^8 - 1 = 255 = 3·5·17.
# 2^16 - 1 = 65535 = 3·5·17·257
# 2^136 - 1 = ? (ord_{q^2}(2) = 136)
print(f"\n=== Factor pattern in 2^k - 1 ===")
from sympy import factorint
for k in [8, 16, 32, 64, 128, 136]:
    n = 2**k - 1
    f = factorint(n)
    small = {p: e for p, e in f.items() if p < 10**8}
    big = [p for p in f if p >= 10**8]
    print(f"  2^{k} - 1: small={small}, big={[f'{p:.6e}' for p in big[:3]]}")

# Try: do 2675591024983, 1989596835661185763 divide 2^k - 1 for small k?
print(f"\n=== Multiplicative order checks ===")
big_primes = [
    ("19841", 19841),
    ("2675591024983", 2675591024983),
    ("1989596835661185763", 1989596835661185763),
    ("huge1 (29-digit)", 19777954930411325600042121089),
]
for name, p in big_primes:
    # ord of 2 mod p
    pm1 = p - 1
    # check small order
    o = 1; x = 2 % p
    while x != 1 and o < 1000:
        o += 1; x = (x * 2) % p
    if x == 1:
        print(f"  {name}: ord_p(2) = {o}  ({pm1}/{o} = {pm1//o if pm1 % o == 0 else 'not divisor'})")
    else:
        # try divisors of p-1
        from sympy import divisors
        divs = divisors(pm1)
        for d in sorted(divs):
            if pow(2, d, p) == 1:
                print(f"  {name}: ord_p(2) = {d} (from p-1 divisors); (p-1)/ord = {pm1//d}")
                break
        else:
            print(f"  {name}: ord_p(2) not found in p-1 divisors? bug.")
