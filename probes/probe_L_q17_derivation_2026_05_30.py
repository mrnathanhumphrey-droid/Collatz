"""
probe_L_q17_derivation_2026_05_30.py

T-L-1: Construct L at q=17, diagonalize, identify dominant eigenvalues.

Basis ordering (locked):
  QR sorted by 2^k: [1, 2, 4, 8, 16, 15, 13, 9]
  NQR sorted by 3·2^k: [3, 6, 12, 7, 14, 11, 5, 10]
  Full basis: QR ++ NQR

Predicted eigenvalue from damped-osc fit: z ≈ 0.034 + 0.068i, |z| ≈ 0.076, arg = arctan(2).
If z = c·(1+2i) for clean c, we've located it analytically.
"""
from __future__ import annotations
import sys
import numpy as np
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8

def chi_2(x):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

# Build basis in locked order
QR = [pow(2, k, q) for k in range(ord_2)]  # [1, 2, 4, 8, 16, 15, 13, 9]
NQR_gen = 3
NQR = [(NQR_gen * pow(2, k, q)) % q for k in range(ord_2)]  # [3, 6, 12, 7, 14, 11, 5, 10]
basis = QR + NQR
print(f"Basis (QR then NQR): {basis}")
assert set(basis) == set(range(1, q)), "basis must cover (Z/17)*"

# Index lookup
idx = {s: i for i, s in enumerate(basis)}
n = q - 1

# Geom(4) marginal weights for a ≡ r mod 8
def P_geom4(r):
    if r == 0:
        return Fraction(3, 2**16 - 1)
    return Fraction(3 * 2**(16 - 2*r), 2**16 - 1)

# pow_inv2[a] = 2^(-a) mod q for a = 0..7
inv2_q = pow(2, -1, q)
pow_inv2 = [pow(inv2_q, a, q) for a in range(ord_2)]

# === L_dom: σ → 2^(-a)·σ ===
L_dom = np.zeros((n, n), dtype=np.complex128)
for i, sigma in enumerate(basis):
    for r in range(ord_2):
        sigma_new = (pow_inv2[r] * sigma) % q
        if sigma_new == 0: continue
        j = idx[sigma_new]
        L_dom[j, i] += float(P_geom4(r))

print(f"\nL_dom column sums: {[f'{x:.4f}' for x in L_dom.sum(axis=0)]}")
# Verify QR→QR and NQR→NQR (block diagonal)
QR_idx = [idx[s] for s in QR]; NQR_idx = [idx[s] for s in NQR]
QR_to_NQR = sum(L_dom[i, j].real for i in NQR_idx for j in QR_idx)
NQR_to_QR = sum(L_dom[i, j].real for i in QR_idx for j in NQR_idx)
print(f"L_dom QR→NQR mass: {QR_to_NQR:.2e}  (should be 0)")
print(f"L_dom NQR→QR mass: {NQR_to_QR:.2e}  (should be 0)")

# === L_subdom: σ → 2^(-a)·(σ - 2k) for k = ±1 ===
L_sub = np.zeros((n, n), dtype=np.complex128)
for i, sigma in enumerate(basis):
    for r in range(ord_2):
        for k in [-1, +1]:
            shifted = (sigma - 2*k) % q
            if shifted == 0:
                continue  # special case
            sigma_new = (pow_inv2[r] * shifted) % q
            if sigma_new == 0: continue
            j = idx[sigma_new]
            L_sub[j, i] += float(P_geom4(r))

# Cross-coset structure
QR_to_NQR_sub = sum(L_sub[i, j].real for i in NQR_idx for j in QR_idx)
NQR_to_QR_sub = sum(L_sub[i, j].real for i in QR_idx for j in NQR_idx)
print(f"\nL_sub QR→NQR mass: {QR_to_NQR_sub:.4f}  (should be nonzero — cross-coset)")
print(f"L_sub NQR→QR mass: {NQR_to_QR_sub:.4f}  (should be nonzero — cross-coset)")

# Full L_m=1: L_dom + ε · (L_sub - L_dom·2)/normalization
# More precisely: P(k=0) = (1 - 2ε), P(k=±1) = ε each. So
#   L_m = (1 - 2ε) L_dom + ε · L_sub (where L_sub already sums k = ±1)
eps = 1.0 / 256
L_m1 = (1 - 2*eps) * L_dom + eps * L_sub
print(f"\nL_m=1 col sums: {[f'{x.real:.4f}' for x in L_m1.sum(axis=0)]}")

# === Eigenvalues of L_dom ===
print(f"\n=== Eigenvalues of L_dom (asymptotic operator) ===")
evals_dom = np.linalg.eigvals(L_dom)
idx_sort = np.argsort(np.abs(evals_dom))[::-1]
for k, lam in enumerate(evals_dom[idx_sort]):
    print(f"  λ_{k}: {lam}  |λ| = {abs(lam):.6f}")

# === Eigenvalues of L_m=1 ===
print(f"\n=== Eigenvalues of L_m=1 (ε = 1/256 sub-dom included) ===")
evals_m1, evecs_m1 = np.linalg.eig(L_m1)
idx_sort = np.argsort(np.abs(evals_m1))[::-1]
for k, lam in enumerate(evals_m1[idx_sort]):
    print(f"  λ_{k}: {lam}  |λ| = {abs(lam):.6f}")

# === Check predicted eigenvalue z ≈ 0.034 + 0.068i ===
print(f"\n=== Searching for predicted eigenvalue z ≈ 0.034 + 0.068i ===")
predicted_z = 0.034 + 0.068j
print(f"Predicted |z| = {abs(predicted_z):.4f}, arg = {np.degrees(np.angle(predicted_z)):.2f}°")
print(f"Predicted arctan(2) = {np.degrees(np.arctan(2)):.4f}°")

print(f"\nClosest matches in L_m=1 spectrum:")
distances = [abs(lam - predicted_z) for lam in evals_m1]
top_matches = np.argsort(distances)[:5]
for i in top_matches:
    lam = evals_m1[i]
    print(f"  λ = {lam}, distance from z = {distances[i]:.6e}")

# Also check conjugate
distances_conj = [abs(lam - predicted_z.conjugate()) for lam in evals_m1]
top_matches_conj = np.argsort(distances_conj)[:5]
print(f"\nClosest to conjugate z̄ = 0.034 - 0.068i:")
for i in top_matches_conj:
    lam = evals_m1[i]
    print(f"  λ = {lam}, distance from z̄ = {distances_conj[i]:.6e}")

# === χ_2 in basis ===
chi_2_vec = np.array([chi_2(s) for s in basis], dtype=np.complex128)
print(f"\nχ_2 vector (in basis order): {[int(c.real) for c in chi_2_vec]}")
# = +1 for QR (first 8), -1 for NQR (last 8)

# Project χ_2 onto eigenvectors of L_m=1
chi_2_in_eigbasis = np.linalg.solve(evecs_m1, chi_2_vec)
print(f"\nχ_2 in eigenbasis of L_m=1 (showing |coef| > 0.01):")
for k in range(n):
    if abs(chi_2_in_eigbasis[k]) > 0.01:
        print(f"  λ = {evals_m1[k]:.4f}, coef = {chi_2_in_eigbasis[k]:.6f}")
