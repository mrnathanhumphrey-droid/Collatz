"""
Fantini-Rella Test 3 — sequence-level reframe.

Treat our c(m) sequence (m=0..5) as a Stokes-constant-like sequence S_m.
Test if any of these generate a Dirichlet L:
  L(s) = Σ_{m>=1} S_m / m^s

Also test:
  - Divisor-sum / twisted-divisor-sum structure: c(m) = Σ_{d|m} χ(d) (...)
  - Modular generating series y(q) = Σ c(m) q^m at small rational q
  - PSLQ generated L-shape values against known Dirichlet L-bank

Data:
  c(0) = 19/127       (exact)
  c(1) = 265011804960406635465672455997699 / 1730087916969634762193659498034425  (exact)
  c(2..5) numeric to 30-50 digits
  c_inf = 0.15298912060588517527891674877413229926086222622334 (50 digits)
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq, zeta
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 50

# === Inputs ===
c = {}
c[0] = mpf(19) / mpf(127)
c[1] = mpf(265011804960406635465672455997699) / mpf(1730087916969634762193659498034425)
c[2] = mpf("0.1532479207790874367716967397510852719979")
c[3] = mpf("0.1530053316915140267018082388619123587748")
c[4] = mpf("0.1529887090468211652268046601150049898549764588045")
c[5] = mpf("0.152988999413521905309166845365")
c_inf = mpf("0.15298912060588517527891674877413229926086222622334")

print("=== Data ===")
for m in range(6):
    print(f"  c({m}) = {c[m]}")
print(f"  c_inf  = {c_inf}")

# Delta_m
delta = {m: c[m] - c_inf for m in range(6)}
print("\n=== Δ_m = c(m) - c_inf ===")
for m in range(6):
    print(f"  Δ_{m} = {float(delta[m]):+.6e}")

# === Sequence-S candidate variants ===
print("\n=== Variant 1: S_m = c(m) ===")
print(f"  Dirichlet series at s=2: D(2) = Σ_{{m=1}}^5 c(m)/m^2 = {sum(c[m]/m**2 for m in range(1,6))}")
print(f"  Dirichlet series at s=3: D(3) = Σ c(m)/m^3                  = {sum(c[m]/m**3 for m in range(1,6))}")
print(f"  Dirichlet series at s=4: D(4)                                = {sum(c[m]/m**4 for m in range(1,6))}")

print("\n=== Variant 2: S_m = Δ_m ===")
# Δ_m alternates in sign and decays geometrically per damped-osc model
delta_series = lambda s: sum(delta[m]/m**s for m in range(1,6))
for s in [1, 1.5, 2, 3, 4]:
    print(f"  Σ Δ_m / m^{s} = {delta_series(mpf(s))}")

print("\n=== Variant 3: S_m = Δ_m · ρ^(-m), normalized ===")
# Extract ρ from earlier fits: from sign pattern of Δ and damping
# Δ_1 = -0.00... (sign), Δ_2 = +0.0003..., Δ_3 = +0.0001..., Δ_4 ~ -3e-7, Δ_5 ~ -1e-7
# Sign: Δ_1=-, Δ_2=+, Δ_3=+, Δ_4=-, Δ_5=-. The transfer-op extraction gave ρ ≈ 0.076.
rho_est = mpf("0.076")
S_normalized = {m: delta[m] / rho_est**m for m in range(1,6)}
print(f"  ρ_est = {rho_est}")
for m in range(1,6):
    print(f"  S_{m} = Δ_{m}/ρ^{m} = {float(S_normalized[m]):+.6e}")

# === Test divisor-sum structure ===
print("\n=== Variant 4: Twisted divisor sum test ===")
# If c(m) = Σ_{d|m} χ(d) · g(m/d), then for χ = Legendre mod q, small q,
# the sequence has specific structure. Test:
#   q=2: χ_2(1)=1, χ_2(d)=0 for even d
#   q=3: χ_3(d) = Legendre, period 3
#   q=5, 13, 17

def divisors(m):
    return [d for d in range(1, m+1) if m % d == 0]

# Test if c(m) - Σ_{d|m} (some character) has a clean pattern
# Simple test: does c(m) match a constant + Dirichlet series of small chars?
print("  c(m) values (numeric):")
for m in range(1, 6):
    print(f"    c({m}) = {float(c[m]):.10f}, divisors of {m}: {divisors(m)}")

# === Test 4-cycle structure ===
# The damped osc had θ near arctan(2) ≈ 1.107 rad. cos(m θ) cycles with period 2π/θ ≈ 5.68.
# So roughly 4-6 period cycle in sign of Δ_m. Test if c(m) - A·cos(mθ+φ) has clean form.
print("\n=== Variant 5: Extracted damped-osc fit residuals ===")
# From earlier work: c(m) ≈ A + B·ρ^m·cos(mθ+φ)
# With c_inf = A, fit B,ρ,θ,φ to (Δ_1,...,Δ_5)
# Use deltas to fit; we have ρ ≈ 0.076, θ near arctan(2) ≈ 1.107
# (Result from earlier: u ≈ 0.0688, v ≈ 0.00578, K ≈ -0.001168)
# Actually let me just compute residuals after subtracting fit

# fit B, theta, phi by least squares to Δ_1..Δ_5 with rho fixed
import numpy as np
m_arr = np.array([1,2,3,4,5])
delta_arr = np.array([float(delta[m]) for m in range(1,6)])
rho_f = 0.076
theta_f = np.arctan(2)  # 1.1071
# Δ_m = B·ρ^m·cos(mθ+φ); split: cos(mθ+φ) = cos(mθ)cos(φ) - sin(mθ)sin(φ)
A_cos = rho_f**m_arr * np.cos(m_arr * theta_f)
A_sin = rho_f**m_arr * np.sin(m_arr * theta_f)
M = np.column_stack([A_cos, -A_sin])
sol, res, _, _ = np.linalg.lstsq(M, delta_arr, rcond=None)
B_cos, B_sin = sol
B_amp = np.sqrt(B_cos**2 + B_sin**2)
phi_fit = np.arctan2(B_sin, B_cos)
print(f"  Fit: ρ={rho_f}, θ={theta_f:.4f}, B_amp={B_amp:.6e}, φ={phi_fit:.4f}")
fit_pred = B_amp * rho_f**m_arr * np.cos(m_arr * theta_f + phi_fit)
residuals = delta_arr - fit_pred
print("  Residuals (Δ_m - fit):")
for m_i, r in zip(m_arr, residuals):
    print(f"    m={m_i}: {r:+.3e}")

# === Variant 6: Look for known Dirichlet L identities ===
print("\n=== Variant 6: PSLQ Σ-objects against known L ===")
# Compute Σ S_m / m^s for several variants, PSLQ against:
# {1, L(1,χ_q), L(2,χ_q), log 2, log 17, pi, 1/sqrt(p), zeta(2), zeta(3), Catalan}

def build_chi(q_val, g, r, phi_q):
    table = [mpc(0)] * q_val
    x = 1
    for k in range(phi_q):
        table[x] = exp(2 * pi * mpc(0, 1) * r * k / phi_q)
        x = (x * g) % q_val
    return table

def L1_chi(chi_table, q_val):
    total = mpc(0)
    for a in range(1, q_val):
        total += chi_table[a] * digamma(mpf(a)/mpf(q_val))
    return -total / mpf(q_val)

def L_at_2(chi_table, q_val):
    # L(2, χ) via Hurwitz: L(2,χ) = (1/q^2) Σ_a χ(a) ζ(2, a/q)
    from mpmath import zeta as mp_zeta
    total = mpc(0)
    for a in range(1, q_val):
        if chi_table[a] != 0:
            total += chi_table[a] * mp_zeta(2, mpf(a)/mpf(q_val))
    return total / mpf(q_val)**2

# Dirichlet basis
L_4_17 = L1_chi(build_chi(17, 3, 4, 16), 17)
L_8_17 = L1_chi(build_chi(17, 3, 8, 16), 17)
L_2_17 = L1_chi(build_chi(17, 3, 2, 16), 17)
L_1_17 = L1_chi(build_chi(17, 3, 1, 16), 17)
L_4_5 = L1_chi(build_chi(5, 2, 1, 4), 5)
L_2_5 = L1_chi(build_chi(5, 2, 2, 4), 5)
L_4_13 = L1_chi(build_chi(13, 2, 3, 12), 13)
L2_4_17 = L_at_2(build_chi(17, 3, 4, 16), 17)
L2_2_17 = L_at_2(build_chi(17, 3, 8, 16), 17)
L2_4_5 = L_at_2(build_chi(5, 2, 1, 4), 5)
G_catalan = mpf("0.9159655941772190150546035149323841107741493742816721342664981196217630197762547694794")
zeta3 = zeta(3)

elem_basis = [
    ("1", mpf(1)),
    ("log2", log(mpf(2))), ("log5", log(mpf(5))), ("log17", log(mpf(17))),
    ("pi", pi), ("pi^2", pi**2), ("1/pi", 1/pi),
    ("1/sqrt(5)", 1/sqrt(mpf(5))), ("1/sqrt(17)", 1/sqrt(mpf(17))), ("1/sqrt(85)", 1/sqrt(mpf(85))),
    ("Catalan", G_catalan),
    ("zeta(2)", pi**2/6), ("zeta(3)", zeta3),
    ("ReL1_4_17", mpf(L_4_17.real)), ("ImL1_4_17", mpf(L_4_17.imag)),
    ("ReL1_8_17", mpf(L_8_17.real)), ("ImL1_8_17", mpf(L_8_17.imag)),
    ("ReL1_2_17", mpf(L_2_17.real)), ("ImL1_2_17", mpf(L_2_17.imag)),
    ("ReL1_4_5", mpf(L_4_5.real)), ("ImL1_4_5", mpf(L_4_5.imag)),
    ("ReL1_4_13", mpf(L_4_13.real)), ("ImL1_4_13", mpf(L_4_13.imag)),
    ("ReL2_4_17", mpf(L2_4_17.real)), ("ImL2_4_17", mpf(L2_4_17.imag)),
    ("ReL2_4_5", mpf(L2_4_5.real)), ("ImL2_4_5", mpf(L2_4_5.imag)),
]

# Test multiple "sequence summary" objects
summary_objects = [
    ("Σ c(m)/m^2", sum(c[m]/mpf(m)**2 for m in range(1,6))),
    ("Σ c(m)/m^3", sum(c[m]/mpf(m)**3 for m in range(1,6))),
    ("Σ Δ_m/m", sum(delta[m]/mpf(m) for m in range(1,6))),
    ("Σ Δ_m/m^2", sum(delta[m]/mpf(m)**2 for m in range(1,6))),
    ("Σ Δ_m/m^3", sum(delta[m]/mpf(m)**3 for m in range(1,6))),
    ("c(1) - 19/127", c[1] - c[0]),
    ("c_inf - 19/127 = Δ_∞", c_inf - mpf(19)/mpf(127)),
    ("(c_inf - 19/127) * 127", (c_inf - mpf(19)/mpf(127)) * mpf(127)),
    ("c_inf * 127", c_inf * mpf(127)),
    ("c_inf - c(5)", c_inf - c[5]),
    ("c_inf - c(4)", c_inf - c[4]),
]

print("\nSummary objects PSLQ vs elementary+Dirichlet basis:")
names = [n for n, v in elem_basis]
vals = [v for n, v in elem_basis]
for tname, tval in summary_objects:
    print(f"  Target: {tname} = {float(tval):.10e}")
    for tol_exp in [10, 15, 25, 35]:
        tol = mpf(10) ** (-tol_exp)
        rel = pslq([tval] + vals, tol=tol, maxcoeff=200)
        if rel is None:
            print(f"    tol=10^-{tol_exp}: nope")
        elif rel[0] != 0:
            terms = [f"({rel[0]:+d})*T"] + [f"({c:+d})*{n}" for c, n in zip(rel[1:], names) if c != 0]
            if len(terms) <= 6:
                print(f"    tol=10^-{tol_exp}: {' '.join(terms)}")
            else:
                print(f"    tol=10^-{tol_exp}: complex rel ({len(terms)} terms), skip")
            break

print("\n=== Done ===")
