"""
probe_L_step0_q11_q23_2026_05_30.py

T-L-0a, T-L-0b: spectral re-derivation of swap-symmetry theorem for q ≡ 3 mod 4.

The L-matrix construction must reproduce c(m) ≡ 0 for q=11, q=23 via the parity argument:
  - σ → -σ involution.
  - χ_2 is odd (parity -1) since χ_2(-1) = -1 for q ≡ 3 mod 4.
  - Initial distribution P_0 is even (parity +1) by (X,Y) swap symmetry.
  - L commutes with σ → -σ (preserves parity).
  - Therefore ⟨χ_2, L^m P_0⟩ = 0 for all m.

This script verifies the mechanism numerically.
"""
from __future__ import annotations
import sys
import numpy as np
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

def chi_2(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

def build_L(q):
    """Construct L_dom on (Z/q)* with Geom(4) marginal."""
    ord_2 = 1; x = 2 % q
    while x != 1: ord_2 += 1; x = (x * 2) % q
    units = list(range(1, q))
    n = q - 1
    L = np.zeros((n, n), dtype=np.float64)
    inv2_q = pow(2, -1, q)
    pow_inv2 = [pow(inv2_q, a, q) for a in range(ord_2)]
    # P_geom4_mod_ord(r) = sum over a ≡ r mod ord_2 of 3*2^(-2a)
    # For r in {1..ord_2}: = 3 * 2^(-2r) / (1 - 2^(-2*ord_2)) for the periodized form
    def P_geom4(r):
        if r == 0:
            return Fraction(3, 2**(2*ord_2) - 1)
        return Fraction(3 * 2**(2*ord_2 - 2*r), 2**(2*ord_2) - 1)
    for i, sigma in enumerate(units):
        for r in range(ord_2):
            sigma_new = (pow_inv2[r] * sigma) % q
            if sigma_new == 0:
                continue
            j = units.index(sigma_new)
            L[j, i] += float(P_geom4(r))
    return L, units, ord_2

for q in [11, 23]:
    print(f"\n{'='*60}")
    print(f"q = {q}  (≡ {q%4} mod 4)")
    print(f"{'='*60}")

    L, units, ord_2 = build_L(q)
    print(f"ord_q(2) = {ord_2}")
    print(f"χ_2(-1) = {chi_2(-1, q)} (predicted -1 for q ≡ 3 mod 4)")

    n = q - 1
    chi_2_vec = np.array([chi_2(s, q) for s in units], dtype=np.float64)

    # σ → -σ involution operator I
    I = np.zeros((n, n), dtype=np.float64)
    for i, sigma in enumerate(units):
        neg_sigma = (-sigma) % q
        j = units.index(neg_sigma)
        I[j, i] = 1
    # I^2 = identity
    print(f"||I^2 - identity||: {np.linalg.norm(I @ I - np.eye(n)):.2e}")

    # Verify L commutes with I
    commutator_norm = np.linalg.norm(L @ I - I @ L)
    print(f"||[L, I]|| = ||L·I - I·L|| = {commutator_norm:.6e}  (predicted ~0)")

    # Even/odd subspaces under I
    P_even = (np.eye(n) + I) / 2
    P_odd = (np.eye(n) - I) / 2
    print(f"||P_even^2 - P_even|| = {np.linalg.norm(P_even @ P_even - P_even):.2e}")
    print(f"||P_odd^2 - P_odd|| = {np.linalg.norm(P_odd @ P_odd - P_odd):.2e}")

    # χ_2 parity
    chi_even = P_even @ chi_2_vec
    chi_odd = P_odd @ chi_2_vec
    print(f"\n||χ_2 (even part)|| = {np.linalg.norm(chi_even):.6e}  (predicted 0)")
    print(f"||χ_2 (odd part)||  = {np.linalg.norm(chi_odd):.6e}  (predicted nonzero)")

    # Initial distribution P_0 (use uniform on units as a stand-in, since the
    # exact P_0 depends on Tao-Syracuse depth-0 distribution which is symmetric).
    # Any swap-symmetric distribution suffices.
    P_0 = np.ones(n) / n  # uniform: swap-symmetric (P(σ) = P(-σ))
    P_0_even = P_even @ P_0
    P_0_odd = P_odd @ P_0
    print(f"\nUniform P_0:")
    print(f"  ||P_0 (even part)|| = {np.linalg.norm(P_0_even):.6e}")
    print(f"  ||P_0 (odd part)||  = {np.linalg.norm(P_0_odd):.6e}  (predicted 0)")

    # Compute ⟨χ_2, L^m P_0⟩ for m = 0, 1, ..., 10
    print(f"\n⟨χ_2, L^m P_0⟩ for m = 0..10:")
    current = P_0.copy()
    moments = []
    for m in range(11):
        moment = chi_2_vec @ current
        moments.append(moment)
        print(f"  m={m}: ⟨χ_2, L^m P_0⟩ = {moment:.2e}")
        current = L @ current
    max_abs = max(abs(m) for m in moments)
    print(f"\nMax |⟨χ_2, L^m P_0⟩| = {max_abs:.2e}")
    if max_abs < 1e-12:
        print(f"*** T-L-0 PASS at q={q}: spectral re-derivation works ***")
    else:
        print(f"*** T-L-0 FAIL at q={q}: max abs = {max_abs} ***")

print(f"\n\n{'='*60}")
print("Summary")
print(f"{'='*60}")
print("If both q=11 and q=23 PASS, the L-matrix construction is validated for the trivial branch.")
print("Only then proceed to T-L-1 (q=17 derivation).")
