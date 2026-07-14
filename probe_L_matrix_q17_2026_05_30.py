"""
probe_L_matrix_q17_2026_05_30.py

Construct the propagation operator L explicitly.

L acts on functions on (Z/q)* (16-dim for q=17). Two contributions:
  L_dominant: σ → 2^(-a)·σ mod q with a ~ Geom(4). Preserves cosets.
  L_subdom:   σ → 2^(-a)·(σ - 2k) mod q for k=±1, weight 2^(-8q^(m-1)) per depth m.

The c(m) → c(m+1) transition involves L_m = L_dominant + ε_m · L_subdom
where ε_m = 2^(-8q^(m-1)) (vanishes for m ≥ 2).

For m=1: ε_1 = 2^(-8) = 1/256 (visible cross-coset rate).
For m≥2: ε_m ≈ 0 (essentially L_dominant).

Goal: diagonalize L_dominant + ε_1 · L_subdom, identify dominant eigenvalues
in the χ_2-sector, see if they land along (1+2i) Gaussian direction.
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

units = list(range(1, q))  # (Z/17)* = {1, ..., 16}
chi_2_vec = np.array([chi_2(s) for s in units], dtype=np.float64)
print(f"χ_2 on units: {dict(zip(units, [int(c) for c in chi_2_vec]))}")
qr = [s for s in units if chi_2(s) == 1]
nqr = [s for s in units if chi_2(s) == -1]
print(f"QR coset:  {qr}")
print(f"NQR coset: {nqr}")

# Geom(4) marginal on a: P(a = j) = 3·2^(-2j) for j ≥ 1
# Periodize at mod ord_2 = 8: P(a ≡ r mod 8) = 3·2^(-2r)/(1-2^(-16)) for r ∈ {1..7},
#                              = 3·2^(-16)/(1-2^(-16)) for r = 0.
def P_a_mod_8(r):
    if r == 0:
        return Fraction(3, 2**16 - 1) * Fraction(2**0, 1)
        # = 3 / (2^16 - 1)
    return Fraction(3 * 2**(16 - 2*r), 2**16 - 1)

# Verify normalization
total = sum(P_a_mod_8(r) for r in range(ord_2))
print(f"\nGeom(4) mod-8 marginal: {[float(P_a_mod_8(r)) for r in range(ord_2)]}")
print(f"Sum = {float(total)} (expected 1)")

# pow_inv2 mod q for a ∈ {0..7}
inv2_q = pow(2, -1, q)
pow_inv2 = [pow(inv2_q, a, q) for a in range(ord_2)]
print(f"\n2^(-a) mod {q}: {pow_inv2}")

# === L_dominant: σ → 2^(-a)·σ ===
L_dom = np.zeros((q-1, q-1), dtype=np.float64)
for i, sigma in enumerate(units):
    for r in range(ord_2):
        sigma_new = (pow_inv2[r] * sigma) % q
        if sigma_new == 0:
            continue
        j = units.index(sigma_new)
        L_dom[j, i] += float(P_a_mod_8(r))

print(f"\nL_dominant row sums (should all = 1): {L_dom.sum(axis=0)}")
print(f"L_dominant col sums (should all = 1, since L: col → row distribution... wait)")
print(f"L_dom[:, 0] (col for σ=1): {L_dom[:, 0]}")
print(f"L_dom column sums: {L_dom.sum(axis=0)}")

# Check: L_dom should preserve cosets.
# For σ ∈ QR, image 2^(-a)σ should always be in QR.
qr_mask = np.array([chi_2(s) == 1 for s in units])
nqr_mask = ~qr_mask
print(f"L_dom: QR→QR sum = {L_dom[qr_mask][:, qr_mask].sum() / qr_mask.sum():.6f} (expected 1)")
print(f"L_dom: QR→NQR sum = {L_dom[nqr_mask][:, qr_mask].sum() / qr_mask.sum():.6f} (expected 0)")

# === L_subdom (sub-dominant cross-coset): σ → 2^(-a)·(σ - 2k) for k = ±1 ===
# Weight per step: 2^(-8) = 1/256 at depth m=1
# Within sub-dom, a still ~ Geom(4), independent of k.
# Sum over k = ±1 (2 contributions per per-step sub-dom event).
L_sub = np.zeros((q-1, q-1), dtype=np.float64)
for i, sigma in enumerate(units):
    for r in range(ord_2):
        for k in [-1, +1]:
            sigma_new = (pow_inv2[r] * ((sigma - 2 * k) % q)) % q
            if sigma_new == 0:
                continue
            j = units.index(sigma_new)
            L_sub[j, i] += float(P_a_mod_8(r))

print(f"\nL_sub matrix built.")
print(f"L_sub col sums (mass per col): {L_sub.sum(axis=0)}")
# This counts 2 per (a, k) → if we want the average over (a, k=±1), divide by 2.
# But we want the cross-coset contribution per step at depth m=1, which is 2^(-8) times the L_sub action.

# === L_m = L_dom + ε_m · (L_sub - L_dom_contribution_double_counted) ===
# Actually: the full per-step transition combines dominant and sub-dom. The full distribution at depth m+1:
#   p_{m+1}(σ') = (1 - 2ε_m) · L_dom[σ', σ] p_m(σ) + ε_m · L_sub[σ', σ] p_m(σ)
# where ε_m = 2^(-8q^{m-1}) and the (1 - 2ε_m) is normalization to keep total prob = 1.
# (Subtle: L_sub already sums over k=±1, so total subdom weight = 2 · ε_m on σ at this step.)

# Hmm actually let me think more carefully. The (a, k) joint at deepening event m:
#   P(a = j, k = κ) for j ≥ 1, κ ∈ Z:
#     κ = 0: 2^(-2j)/Z
#     κ = ±1: 2^(-2j - 8q^(m-1))/Z
#     ...
#   Z = (1 + 2·2^(-8q^(m-1)) + ...)/3
# So P(κ = 0) ≈ 1 - 2·2^(-8q^(m-1)), P(κ = ±1) ≈ 2^(-8q^(m-1)) each.

# For step m=1, ε_1 = 2^(-8q^0) = 2^(-8) = 1/256.
# L_m at depth m: column σ → row σ' has weights:
#   L_dom[σ', σ] · (1 - 2·ε_m) + L_sub_k=+1[σ', σ] · ε_m + L_sub_k=-1[σ', σ] · ε_m
# But L_sub above is the SUM over k = ±1, so we use L_sub[σ', σ] · ε_m.

eps_m1 = 1.0 / 256
L_m1 = (1 - 2*eps_m1) * L_dom + eps_m1 * L_sub
print(f"\nL_m=1 col sums: {L_m1.sum(axis=0)}  (should be 1)")

# === L_m for m ≥ 2: essentially L_dom (since ε_m ≈ 0) ===
L_asymp = L_dom

# === Diagonalize L_asymp ===
print(f"\n=== Eigenvalues of L_asymptotic (=L_dom) ===")
evals_dom, evecs_dom = np.linalg.eig(L_asymp)
# Sort by |λ|
idx = np.argsort(np.abs(evals_dom))[::-1]
evals_sorted = evals_dom[idx]
evecs_sorted = evecs_dom[:, idx]
for k, lam in enumerate(evals_sorted):
    print(f"  λ_{k} = {lam}  |λ| = {abs(lam):.6f}")

# === Diagonalize L_m=1 ===
print(f"\n=== Eigenvalues of L_m=1 (with ε_1 = 1/256 sub-dom) ===")
evals_m1, evecs_m1 = np.linalg.eig(L_m1)
idx = np.argsort(np.abs(evals_m1))[::-1]
evals_m1_sorted = evals_m1[idx]
evecs_m1_sorted = evecs_m1[:, idx]
for k, lam in enumerate(evals_m1_sorted):
    print(f"  λ_{k} = {lam}  |λ| = {abs(lam):.6f}")

# === Project χ_2 onto eigenvectors of L_dom ===
# χ_2 expressed in σ basis: chi_2_vec
# Eigenvector basis: evecs_dom columns
# Coefficients: solve evecs · c = chi_2_vec
print(f"\n=== χ_2 in eigenbasis of L_dom ===")
chi_2_in_eigbasis = np.linalg.solve(evecs_dom, chi_2_vec)
print(f"  Coefficients on each eigenvector:")
for k, (lam, coef) in enumerate(zip(evals_dom, chi_2_in_eigbasis)):
    if abs(coef) > 1e-10:
        print(f"    eigenvec_{k} (λ={lam:.4f}): coef = {coef:.6f}")

# === Damped-osc eigenvalue check ===
# From our fit, the damped-osc eigenvalue z ≈ 0.034 + 0.068i with |z| ≈ 0.076.
# Does L_m=1 have such an eigenvalue? Or L_dom?
print(f"\n=== Looking for damped-osc eigenvalue z ≈ 0.034 + 0.068i ===")
target_z = 0.034 + 0.068j
for k, lam in enumerate(evals_m1):
    if abs(lam - target_z) < 0.05 or abs(lam - target_z.conjugate()) < 0.05:
        print(f"  MATCH: λ_{k} = {lam}, diff from z = {lam - target_z}")

# === Initial vector v_0 (P at depth 0 — the c(0) configuration) ===
# Conditional on v_q(D) = 0 (well, ≥ 0 which is everything excluding D = 0).
# Or more precisely: conditional distribution of σ_0 = D mod q over (X-Y) iid Tao-Syracuse.
# We could compute this via the depth-0 character sum machinery.

# For now, just compute c_∞ predicted from the operator analysis.
# If the χ_2 component is preserved by L_dom (eigenvalue +1), c_∞ = χ_2 moment under the
# eigenvector of L_dom corresponding to λ = +1.

# But the χ_2 might decompose across multiple L_dom eigenvalues. The asymptotic c_∞ would be
# the coefficient on the λ=1 eigenvector.

# Let's check: which eigenvector(s) of L_dom does χ_2 project onto?
eigval1_idx = np.argmin(np.abs(evals_dom - 1.0))
print(f"\n  Eigenvalue closest to 1: λ_{eigval1_idx} = {evals_dom[eigval1_idx]}")
print(f"  Eigenvector: {np.real(evecs_dom[:, eigval1_idx])}")
print(f"  χ_2 coefficient on this eigenvec: {chi_2_in_eigbasis[eigval1_idx]:.6f}")
