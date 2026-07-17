"""
Path (a): Fourier-coefficient framing.

Extract leading eigenvector mu_inf of L (transfer op at q=17) to high precision.
Test each component against:
  - Q(i) algebraic numbers
  - Cyclotomic Q(zeta_17), Q(zeta_8)
  - Gauss-sum-like ratios involving sqrt(17), sqrt(5)
  - Simple rationals with small denominator

If mu_inf components have clean closed forms, c_inf = <chi_2, mu_inf> inherits one.

Note: L = L_dom + eps * (L_sub - 2 L_dom) with eps = 2^(-8), built on basis
of (Z/17)* with QR cosset first then NQR coset.
"""
from __future__ import annotations
import sys
from fractions import Fraction
from mpmath import mp, mpf, mpc, sqrt, pi, exp, log, pslq, eig
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 50

q = 17
ord_2 = 8

def chi_2(x):
    x = x % q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

# Basis: QR then NQR, sorted by powers
QR = [pow(2, k, q) for k in range(ord_2)]   # 2 is a QR mod 17 since 17 ≡ ±1 mod 8
NQR_gen = 3                                  # check: 3 mod 17 is NR? 3^8 mod 17 = ?
# Quick check
print(f"chi_2(3) = {chi_2(3)} (NR if -1)")
NQR = [(NQR_gen * pow(2, k, q)) % q for k in range(ord_2)]
basis_set = QR + NQR
print(f"Basis: {basis_set}")
assert set(basis_set) == set(range(1, q))
idx = {s: i for i, s in enumerate(basis_set)}
n = q - 1  # 16

# === Build L_dom (mpmath) ===
# Geom(4) marginal: a ~ {1..8}, P(a=r) = 3 * 2^(16-2r) / (2^16 - 1)
def P_geom4(r):
    if r < 1 or r > 8: return mpf(0)
    return mpf(3) * mpf(2)**(16 - 2*r) / mpf(2**16 - 1)

# Check normalization
total_p = sum(P_geom4(r) for r in range(1, 9))
print(f"Sum of Geom(4) probabilities for r=1..8: {total_p}")

# 2^(-a) mod q for a = 1..8
inv2 = pow(2, -1, q)  # 9 mod 17
pow_inv2 = [pow(inv2, a, q) for a in range(1, 9)]
print(f"2^-1 mod 17 = {inv2}")
print(f"2^-a mod 17 for a=1..8: {pow_inv2}")

# L_dom[j, i] = sum over a of P(a) if 2^-a * basis[i] = basis[j]
L_dom = [[mpf(0) for _ in range(n)] for _ in range(n)]
for i, sigma in enumerate(basis_set):
    for r in range(1, 9):
        sigma_new = (pow_inv2[r-1] * sigma) % q
        if sigma_new == 0: continue
        j = idx[sigma_new]
        L_dom[j][i] = L_dom[j][i] + P_geom4(r)

# Check column sums
print("L_dom column sums (should be 1):", [str(sum(L_dom[j][i] for j in range(n)))[:10] for i in range(n)])

# Verify block structure: QR -> QR, NQR -> NQR
qr_idx = [idx[s] for s in QR]
nqr_idx = [idx[s] for s in NQR]
qr_to_nqr = sum(L_dom[i][j] for i in nqr_idx for j in qr_idx)
nqr_to_qr = sum(L_dom[i][j] for i in qr_idx for j in nqr_idx)
print(f"L_dom QR->NQR mass = {qr_to_nqr} (should be 0)")
print(f"L_dom NQR->QR mass = {nqr_to_qr} (should be 0)")

# === Build L_sub (sigma -> 2^-a * (sigma - 2k) for k = pm 1) ===
L_sub = [[mpf(0) for _ in range(n)] for _ in range(n)]
for i, sigma in enumerate(basis_set):
    for r in range(1, 9):
        for k in [-1, +1]:
            shifted = (sigma - 2*k) % q
            if shifted == 0: continue
            sigma_new = (pow_inv2[r-1] * shifted) % q
            if sigma_new == 0: continue
            j = idx[sigma_new]
            L_sub[j][i] = L_sub[j][i] + P_geom4(r)

eps = mpf(1) / mpf(256)
L_full = [[(mpf(1) - 2*eps) * L_dom[j][i] + eps * L_sub[j][i] for i in range(n)] for j in range(n)]

# === Convert to mpmath matrix and eigendecompose ===
import mpmath
L_mp = mpmath.matrix(L_full)

print("\n=== Eigendecomposition ===")
eigvals, eigvecs = mpmath.eig(L_mp)
mods = [(abs(eigvals[k]), k) for k in range(n)]
mods.sort(reverse=True)
print("Top 6 eigenvalues by modulus:")
for r, (mod_v, k) in enumerate(mods[:6]):
    lam = eigvals[k]
    print(f"  #{r+1}: |lam| = {mod_v}, lam = {lam}")

# Find leading eigenvector (largest mod)
top_k = mods[0][1]
top_eval = eigvals[top_k]
print(f"\nTop eigenvalue: {top_eval} (should be 1 for stochastic L)")
print(f"  |1 - top| = {abs(top_eval - 1)}")

# Get the corresponding eigenvector (column)
mu_vec = [eigvecs[i, top_k] for i in range(n)]

# Normalize so sum = 1
s = sum(mu_vec)
mu_vec = [v / s for v in mu_vec]

print("\n=== Leading eigenvector mu_inf normalized (sum=1) ===")
print("idx | basis[idx] | chi_2 | mu_inf[idx]")
for i in range(n):
    print(f"  {i:2d} | {basis_set[i]:2d}    | {chi_2(basis_set[i]):+d}  | {mu_vec[i]}")

# Compute <chi_2, mu_inf>
chi_2_vec = [mpf(chi_2(s)) for s in basis_set]
c_inf_computed = sum(c * m for c, m in zip(chi_2_vec, mu_vec))
print(f"\n<chi_2, mu_inf> = {c_inf_computed}")
c_inf_ref = mpf("0.15298912060588517527891674877413229926086222622334")
print(f"reference c_inf = {c_inf_ref}")
print(f"diff           = {c_inf_computed.real - c_inf_ref}")

# Real / imag of mu components
mu_real = [v.real if hasattr(v, 'real') else mpf(v) for v in mu_vec]
mu_imag = [v.imag if hasattr(v, 'imag') else mpf(0) for v in mu_vec]

# === PSLQ each mu component against algebraic basis ===
print("\n=== PSLQ mu_inf components vs algebraic basis ===")
# Basis: small rationals + sqrt of small primes + zeta_8 = i*sqrt(2)/2 + 1/sqrt(2) + ...
basis_alg = [
    ("1", mpf(1)),
    ("1/16", mpf(1)/16), ("1/17", mpf(1)/17), ("1/8", mpf(1)/8), ("1/5", mpf(1)/5),
    ("1/sqrt(17)", 1/sqrt(mpf(17))),
    ("1/sqrt(5)", 1/sqrt(mpf(5))),
    ("1/sqrt(85)", 1/sqrt(mpf(85))),
    ("sqrt(17)/17", sqrt(mpf(17))/17),
    ("sqrt(5)/5",   sqrt(mpf(5))/5),
    ("(1+sqrt(17))/34", (1+sqrt(mpf(17)))/34),
    ("(sqrt(17)-1)/16", (sqrt(mpf(17))-1)/16),
    ("log(2)/4", log(mpf(2))/4),
    ("log(17)/16", log(mpf(17))/16),
]
names_alg = [n for n, v in basis_alg]
vals_alg = [v for n, v in basis_alg]

found_any = False
for i in range(n):
    v = mu_real[i]
    if abs(v) < mpf(10)**(-30):
        continue
    for tol_exp in [15, 25, 35]:
        tol = mpf(10) ** (-tol_exp)
        rel = pslq([v] + vals_alg, tol=tol, maxcoeff=500)
        if rel is not None and rel[0] != 0:
            terms = [f"({rel[0]:+d})*mu[{basis_set[i]}]"] + [f"({c:+d})*{n}" for c, n in zip(rel[1:], names_alg) if c != 0]
            if len(terms) <= 5:
                print(f"  basis={basis_set[i]:2d} mu={float(v):+.6e} : {' '.join(terms)} (tol=10^-{tol_exp})")
                found_any = True
                break

if not found_any:
    print("  No clean algebraic relation found at any tested precision/basis.")

# === Examine differences between QR and NQR sub-blocks ===
print("\n=== QR vs NQR block structure ===")
mu_qr = [mu_real[idx[s]] for s in QR]
mu_nqr = [mu_real[idx[s]] for s in NQR]
print(f"QR block sum  = {sum(mu_qr)}")
print(f"NQR block sum = {sum(mu_nqr)}")
print(f"Diff (QR-NQR) = {sum(mu_qr) - sum(mu_nqr)}")
print(f"Note: c_inf = <chi_2, mu> = (sum QR) - (sum NQR)")
print(f"Predicted c_inf = {sum(mu_qr) - sum(mu_nqr)}")

# === Examine QR block alone: is it cyclic? ===
print("\n=== QR block detail (8 components) ===")
print("k | 2^k mod 17 | mu_inf[2^k] | normalized to first")
mu_qr_normalized = [mu_qr[i] / mu_qr[0] for i in range(8)]
for i, s in enumerate(QR):
    print(f"{i} | {s:2d}        | {float(mu_qr[i]):+.6e} | {mu_qr_normalized[i]}")

# Test if these ratios are roots of unity / simple algebraic
print("\nPSLQ each mu_qr[k] / mu_qr[0] vs algebraic basis:")
basis_simple = [("1", mpf(1)), ("1/2", mpf(1)/2), ("1/3", mpf(1)/3), ("1/4", mpf(1)/4),
                ("1/sqrt(2)", 1/sqrt(mpf(2))), ("sqrt(2)", sqrt(mpf(2)))]
names_s = [n for n, v in basis_simple]
vals_s = [v for n, v in basis_simple]
for i in range(1, 8):
    r = mu_qr_normalized[i]
    for tol_exp in [15, 25]:
        tol = mpf(10) ** (-tol_exp)
        rel = pslq([r] + vals_s, tol=tol, maxcoeff=100)
        if rel is not None and rel[0] != 0 and abs(rel[0]) <= 10:
            terms = [f"({rel[0]:+d})*r"] + [f"({c:+d})*{n}" for c, n in zip(rel[1:], names_s) if c != 0]
            print(f"  k={i}: r = {float(r):+.6e}: {' '.join(terms)} (tol=10^-{tol_exp})")
            break

print("\n=== Done ===")
