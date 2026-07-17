"""
Lambert series approach at base 256.

The boundary integral involves weights 2^(-8|κ|) over κ ∈ Z\{0}, which sums
Σ_{κ≥1} 256^(-κ)·χ(κ) for various characters χ of (Z/17)*.

We computed:
  ⟨χ_L, μ_∞ shifted by s⟩ = (-χ_L(s) - c_∞)/16 for s ≠ 0  (uniform-within-coset μ_∞)
  χ_L(2κ·2^(-a)) = χ_L(κ)  (since 2 is QR mod 17)

So the boundary structural quantity is the Lambert sum at base 256:
  L_256(χ) := Σ_{κ≥1} 256^(-κ)·χ(κ)
             = T(χ) / (1 - 256^(-17))
  where T(χ) = Σ_{r=1}^{16} χ(r)·256^(-r)

Test: PSLQ c_∞ against {1, L_256(χ_L), L_256(ω), L_256(order-8), products, rationals}.
"""
from __future__ import annotations
import sys, math
from fractions import Fraction
from mpmath import mp, mpf, mpc, pslq, sqrt, log, pi
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 80
q = 17

# Character values at index r ∈ 1..16, indexed by k in chi_k = exp(2πi·k·dlog_3/16)
# Compute for all 16 characters of (Z/17)*

# Primitive root g = 3
g = 3
dlog = {}
x = 1
for k in range(16):
    dlog[x] = k
    x = (x * g) % q

# Character value table chi_k(r)
def chi_k_value(k, r):
    """chi_k(r) = exp(2pi i k dlog(r) / 16). Returns complex mpc."""
    if r % q == 0:
        return mpc(0)
    return mpc(mp.cos(2 * pi * k * dlog[r] / 16), mp.sin(2 * pi * k * dlog[r] / 16))

# T(chi_k) = sum_{r=1}^{16} chi_k(r) * 256^(-r) at EXACT mpf precision
def T_chi(k):
    """Exact Lambert numerator at base 256."""
    s = mpc(0)
    for r in range(1, q):
        s += chi_k_value(k, r) / mpf(256) ** r
    return s

# L_256(chi_k) = T(chi_k) / (1 - 256^(-17))
denom = mpf(1) - mpf(256) ** (-17)
def L_256(k):
    return T_chi(k) / denom

# Compute for all 16 characters
print(f"=== Lambert sums L_256(chi_k) at base 256 ===")
L_values = {}
for k in range(16):
    L_values[k] = L_256(k)
    order = 16 // math.gcd(k, 16) if k > 0 else 1
    print(f"  k={k:2d} (order {order:2d}): Re={float(L_values[k].real):+.20f}, Im={float(L_values[k].imag):+.20f}")

# Specifically: chi_L = chi_8 (Legendre), order 4 = chi_4 or chi_12
print(f"\nLegendre (k=8) Lambert sum: {L_values[8].real}")
print(f"Order-4 (k=4) Lambert sum: Re={L_values[4].real}, Im={L_values[4].imag}")

# c_inf reference
c_inf = mpf("0.15298912060588517527891674877413229926086222622334")
print(f"\nc_inf = {c_inf}")

# Try the naive prediction c_inf = -15 * L_256(chi_L)
naive_pred = -15 * L_values[8].real
print(f"\nNaive prediction c_inf = -15 * L_256(chi_L) = {naive_pred}")
print(f"Actual c_inf                                  = {float(c_inf)}")
print(f"Diff: {float(c_inf - naive_pred):+.3e}  → naive derivation is WRONG (as expected; case B cascade more complex)")

# Now PSLQ c_inf against Lambert basis + rationals
# Note: chi_k = chi_(-k) when conjugated, and Re/Im of L_256(chi_k) give independent reals

print(f"\n=== PSLQ basis ===")
basis_names = []
basis_vals = []

basis_names.append("1");                         basis_vals.append(mpf(1))
basis_names.append("1/17");                      basis_vals.append(mpf(1)/17)
basis_names.append("1/255");                     basis_vals.append(mpf(1)/255)
basis_names.append("1/256");                     basis_vals.append(mpf(1)/256)
basis_names.append("1/(256-1)");                 basis_vals.append(mpf(1)/255)
basis_names.append("256/(256^17 - 1)");          basis_vals.append(mpf(256)/(mpf(256)**17 - 1))
basis_names.append("log(2)");                    basis_vals.append(log(mpf(2)))
basis_names.append("log(17)");                   basis_vals.append(log(mpf(17)))
basis_names.append("1/sqrt(17)");                basis_vals.append(1/sqrt(mpf(17)))
basis_names.append("pi");                        basis_vals.append(pi)

# Lambert sums (real and imag parts) for ALL 16 chars
for k in range(16):
    order = 16 // math.gcd(k, 16) if k > 0 else 1
    val = L_values[k]
    if abs(val.real) > mpf("1e-60"):
        basis_names.append(f"Re L_256(chi_{k},ord{order})")
        basis_vals.append(val.real)
    if abs(val.imag) > mpf("1e-60"):
        basis_names.append(f"Im L_256(chi_{k},ord{order})")
        basis_vals.append(val.imag)

print(f"  Basis size: {len(basis_vals)}")

print(f"\n=== PSLQ c_inf vs Lambert basis ===")
for tol_exp in [10, 15, 20, 25, 30, 35, 40, 45]:
    tol = mpf(10) ** (-tol_exp)
    rel = pslq([c_inf] + basis_vals, tol=tol, maxcoeff=10**4)
    if rel is None:
        print(f"  tol=10^-{tol_exp}: NO relation")
        continue
    if rel[0] != 0:
        terms = [f"({rel[0]:+d})*c_inf"]
        for c, n in zip(rel[1:], basis_names):
            if c != 0:
                terms.append(f"({c:+d})*{n}")
        residual = sum(c*v for c, v in zip(rel, [c_inf]+basis_vals))
        if len(terms) <= 8:
            print(f"  tol=10^-{tol_exp}: RELATION → {' '.join(terms)}  [residual={float(residual):.2e}]")
        else:
            print(f"  tol=10^-{tol_exp}: RELATION ({len(terms)} terms, complex; residual={float(residual):.2e})")
            print(f"    Top 8 terms: {' '.join(terms[:8])}")
    else:
        nz = sum(1 for c in rel[1:] if c != 0)
        print(f"  tol=10^-{tol_exp}: trivial rel (c_inf coef 0, {nz} terms)")

# Also try: c_inf vs PRODUCTS of Lambert sums
print(f"\n=== PSLQ c_inf vs Lambert products ===")
extended_basis_names = list(basis_names)
extended_basis_vals = list(basis_vals)
# Add products of 2 Lambert sums (focus on chi_L Re and chi_4 Re/Im)
key_indices = [k for k, n in enumerate(basis_names) if "L_256" in n][:8]  # first few Lambert entries
for i in key_indices:
    for j in key_indices:
        if i <= j:
            prod_val = basis_vals[i] * basis_vals[j]
            if abs(prod_val) > mpf("1e-60"):
                extended_basis_names.append(f"({basis_names[i]})*({basis_names[j]})")
                extended_basis_vals.append(prod_val)

print(f"  Extended basis size: {len(extended_basis_vals)}")
for tol_exp in [25, 35, 45]:
    tol = mpf(10) ** (-tol_exp)
    rel = pslq([c_inf] + extended_basis_vals, tol=tol, maxcoeff=10**3)
    if rel is None:
        print(f"  tol=10^-{tol_exp}: NO relation")
        continue
    if rel[0] != 0:
        terms = [f"({rel[0]:+d})*c_inf"]
        for c, n in zip(rel[1:], extended_basis_names):
            if c != 0:
                terms.append(f"({c:+d})*{n}")
        residual = sum(c*v for c, v in zip(rel, [c_inf]+extended_basis_vals))
        if len(terms) <= 10:
            print(f"  tol=10^-{tol_exp}: RELATION → {' '.join(terms)}  [res={float(residual):.2e}]")
        else:
            print(f"  tol=10^-{tol_exp}: RELATION ({len(terms)} terms; res={float(residual):.2e})")
            print(f"    Top 10: {' '.join(terms[:10])}")
    else:
        nz = sum(1 for c in rel[1:] if c != 0)
        print(f"  tol=10^-{tol_exp}: trivial rel (c_inf coef 0, {nz} terms)")

# Single-ratio scan: c_inf / each Lambert sum (real part)
print(f"\n=== Single-ratio scan: c_inf / L_256(chi_k).real ===")
from fractions import Fraction
hits = []
for k in range(16):
    v = L_values[k].real
    if abs(v) < mpf("1e-30"):
        continue
    r = c_inf / v
    rf = float(r)
    if abs(rf) > 1e6 or abs(rf) < 1e-6:
        continue
    f = Fraction(rf).limit_denominator(1000)
    approx = mpf(f.numerator)/mpf(f.denominator)
    diff = r - approx
    if abs(diff) < mpf("1e-6"):
        order = 16 // math.gcd(k, 16) if k > 0 else 1
        hits.append((k, order, rf, f, float(diff)))

if hits:
    print(f"  {len(hits)} hits (denom ≤ 1000, diff < 1e-6):")
    for k, order, rf, frac, diff in hits:
        print(f"    c_inf / Re L_256(chi_{k},ord{order}) = {rf:.15f} ≈ {frac}, diff={diff:.2e}")
else:
    print("  No clean ratio hits")
