"""Verify Polya-Vinogradov decomposition for our partial sum.

partial = Σ_{u=0}^{N-1} f(u) = (1/q) Σ_ξ 1̂(ξ) · F̂(ξ)

where 1̂(ξ) = Σ_{v=0}^{N-1} e_q(ξv), F̂(ξ) = Σ_u f(u) e_q(-ξu).
"""
import cmath
import math

q = 27
N = 3

# f(u) = e_q(4^u)
def f(u):
    pow4 = pow(4, u, q)
    return cmath.exp(2j * cmath.pi * pow4 / q)

# Direct partial sum
partial = sum(f(u) for u in range(N))
print(f"Direct partial sum (N={N}, q={q}): {partial}")
print(f"  |partial| = {abs(partial):.6f}")

# Compute all F̂(ξ) for ξ ∈ Z/q
F_hat = [None] * q
for xi in range(q):
    F_hat[xi] = sum(f(u) * cmath.exp(-2j * cmath.pi * xi * u / q) for u in range(q))

# Compute all 1̂(ξ) for ξ ∈ Z/q
one_hat = [None] * q
for xi in range(q):
    one_hat[xi] = sum(cmath.exp(2j * cmath.pi * xi * v / q) for v in range(N))

# Polya-Vinogradov reconstruction
PV = sum(one_hat[xi] * F_hat[xi] for xi in range(q)) / q
print(f"PV reconstructed: {PV}")
print(f"  |PV| = {abs(PV):.6f}")
print(f"  diff = {abs(partial - PV):.6e}")

print()
print("F̂(ξ) values:")
for xi in range(q):
    print(f"  ξ={xi}: F̂(ξ) = {F_hat[xi]:.6f},  |F̂| = {abs(F_hat[xi]):.6f}")
