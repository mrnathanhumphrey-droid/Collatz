"""
probe_c_exact_rationals_2026_05_29.py

c(m) is a finite weighted sum of character sums => RATIONAL at each finite m.
Compute c(0), c(1), c(2) EXACTLY via fractions, see the family of denominators
and how the rational sequence approaches c_inf.

c(0): N(0)/T(0) where N(s)/T(s) are "depth-0-like moments given shift s":
  N(s) = sum_{a,b>=1, s + 2^-a - 2^-b != 0 mod q} 2^{-(a+b)} chi(s + 2^-a - 2^-b)
  T(s) = same with chi -> 1
  c(0) = N(0)/T(0).

c(1): one level-0 deepening, then depth-0-like termination from shift eta_k.
  eta_k mod q = 2k mod q (verified for q=17 from 2^8 = 1 + 15*17 -> eta_1 = -15 ≡ 2 mod 17,
  generally eta_k ≡ 2k mod q via the binomial structure).
  weight on deepening event k: w_k = (1/3) * 2^{-ord_2 |k|} for k!=0, w_0 = 1/3.
  c(1) = sum_k w_k N(2k mod q) / sum_k w_k T(2k mod q).
"""
from fractions import Fraction
import sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8

def chi(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1

inv2_q = pow(2, q - 2, q)
pow_inv2 = [pow(inv2_q, a, q) for a in range(ord_2)]   # 2^-a mod q

def W(a_res):
    """Sum of 2^-a over a >= 1 with a ≡ a_res mod ord_2."""
    if a_res == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - a_res), 2**ord_2 - 1)

# N(s), T(s)
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

c0 = N[0] / T[0]
print(f"c(0) = {c0} = {float(c0):.15f}")
print(f"  19/127 = {Fraction(19,127)} = {19/127:.15f}   match: {c0 == Fraction(19,127)}")

# c(1)
K = 12
num = Fraction(0); den = Fraction(0)
for k in range(-K, K + 1):
    wk = Fraction(1, 3) if k == 0 else Fraction(1, 3 * 2**(ord_2 * abs(k)))
    sk = (2 * k) % q
    num += wk * N[sk]
    den += wk * T[sk]
c1 = num / den
print(f"\nc(1) (K={K} truncation) = {float(c1):.15f}")
print(f"  numerator   = {num}")
print(f"  denominator = {den}")
print(f"  as fraction = {c1.numerator}/{c1.denominator}")
print(f"  size:  num digits = {len(str(c1.numerator))}, den digits = {len(str(c1.denominator))}")
# compare to FFT-measured c(1) = 0.153178230055
print(f"  FFT measured c(1)  = 0.153178230055")
print(f"  diff = {float(c1) - 0.153178230055:.2e}")
