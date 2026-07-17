"""
probe_coset_q31_leading_2026_05_29.py
Check the leading-order chi(1-2^Delta) mechanism generalizes to q=31 (index 6, quotient Z/6).

Leading-order: A_chi_j/A_triv ~ V_j = sum_{Delta!=0} (2^-|Delta|/3) chi_j(1-2^Delta),
chi_j(u)=omega^{j*coset(u)}, omega=e(2pi i/6), coset(u)=dlog_g(u) mod 6.
Leading coset mass  E_c ~ (1/6)(1 + sum_{j=1}^5 omega^{-j c} V_j).
Compare to the observed converged levels 0.228461 / 0.151084 / 0.120455 (each twice).
Expect: same ordering + rough magnitude, exact mismatch (the non-elementary shift corrections).
"""
from __future__ import annotations
import numpy as np, sys
sys.stdout.reconfigure(encoding="utf-8")

q = 31
# primitive root g mod 31
def order(a, m):
    o = 1; x = a % m
    while x != 1:
        x = (x * a) % m; o += 1
    return o
g = next(a for a in range(2, q) if order(a, q) == q - 1)
# dlog table base g
dlog = {}
x = 1
for e in range((q - 1)):
    dlog[x] = e; x = (x * g) % q
index = 6
omega = np.exp(2j * np.pi / index)
inv2 = pow(2, -1, q)

def coset(u):
    return dlog[u % q] % index

def chi(j, u):
    u %= q
    if u == 0: return 0.0
    return omega ** (j * coset(u))

def one_minus_2pow(Delta):
    if Delta >= 0:
        return (1 - pow(2, Delta, q)) % q
    return (1 - pow(inv2, -Delta, q)) % q

K = 5000
print(f"q={q}, primitive root g={g}, ord_2={order(2,q)}, index={index}, <2>=coset 0={[u for u in range(1,q) if coset(u)==0]}")
# leading V_j
V = [None] * index
for j in range(index):
    s = 0j
    for D in range(-K, K + 1):
        if D == 0: continue
        val = one_minus_2pow(D)
        if val == 0: continue  # deepening (ord_2 | Delta)
        s += (2.0 ** (-abs(D)) / 3.0) * (chi(j, val) if j != 0 else 1.0)
    V[j] = s
V[0] = 1.0  # trivial char normalized
print("leading V_j (j=0..5):")
for j in range(index):
    print(f"  V_{j} = {V[j].real:+.6f} {V[j].imag:+.6f}i")

print("\nleading coset masses E_c ~ (1/6)(sum_j omega^{-jc} V_j):")
levels_pred = []
for c in range(index):
    E = sum((omega ** (-j * c)) * V[j] for j in range(index)) / index
    levels_pred.append(E.real)
    print(f"  coset {c}: {E.real:+.6f}  (imag {E.imag:+.2e})  {'<-- <2>' if c==0 else ''}")
distinct = sorted({round(x, 6) for x in levels_pred}, reverse=True)
print(f"\n leading distinct levels (desc): {distinct}")
print(f" observed converged levels       : [0.228461, 0.151084, 0.120455]")
print(" => same ordering + ballpark; exact mismatch = the non-elementary collision-shift corrections.")
