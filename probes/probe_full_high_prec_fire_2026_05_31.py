"""
probe_full_high_prec_fire_2026_05_31.py

Fire everything with all high-precision data:
- c(0) exact = 19/127
- c(1) exact (33-digit num/den)
- c(2) at dps=40
- c(3) at dps=40
- c(4) at dps=50
- c(5) at dps=30 (overnight job)

1. Re-solve damped oscillation recurrence at high precision
2. Extract c_∞, ρ, θ at ~25-digit precision
3. Test eigenvalue z = ρ·e^(iθ) against (1+2i) direction at high precision
4. PSLQ c_∞, Δ_∞ at tol 10^-22 against full L-value basis
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq, identify
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 50

# === High-precision c(m) values ===
c0 = mpf(19) / mpf(127)
c1_num = 265011804960406635465672455997699
c1_den = 1730087916969634762193659498034425
c1 = mpf(c1_num) / mpf(c1_den)
c2 = mpf("0.1532479207790874367716967397510852719979")
c3 = mpf("0.1530053316915140267018082388619123587748")
c4 = mpf("0.1529887090468211652268046601150049898549764588045")
c5 = mpf("0.152988999413521905309166845365")

print(f"=== c(m) values ===")
for m, c in [(0, c0), (1, c1), (2, c2), (3, c3), (4, c4), (5, c5)]:
    print(f"  c({m}) = {c}")

# Δ_m
deltas = [c0 - c0, c1 - c0, c2 - c0, c3 - c0, c4 - c0, c5 - c0]
print(f"\n=== Δ_m ===")
for m, d in enumerate(deltas):
    print(f"  Δ_{m} = {d}")

# === Solve recurrence at high precision ===
import mpmath as mpmath
M = mpmath.matrix([
    [deltas[2], -deltas[1], mpf(1)],
    [deltas[3], -deltas[2], mpf(1)],
    [deltas[4], -deltas[3], mpf(1)],
])
rhs = mpmath.matrix([deltas[3], deltas[4], deltas[5]])
sol = mpmath.lu_solve(M, rhs)
u, v, K = sol[0], sol[1], sol[2]

print(f"\n=== Recurrence fit at full precision ===")
print(f"  u = {u}")
print(f"  v = {v}")
print(f"  K = {K}")

rho = sqrt(v) if v > 0 else mpf("NaN")
A_asymp = K / (1 - u + v)
print(f"  ρ = √v = {rho}")
print(f"  A = K/(1-u+v) = {A_asymp}")

c_inf = c0 + A_asymp
print(f"\n  c_∞ = c(0) + A = {c_inf}")
print(f"  Precision check: this should be ~25 digits at our input precision")

# Eigenvalue z extraction
if abs(u/(2*rho)) <= 1:
    theta = mpmath.acos(u / (2*rho))
    z = rho * (mpmath.cos(theta) + mpmath.mpc(0,1)*mpmath.sin(theta))
    print(f"\n  θ = {theta} rad = {mpmath.degrees(theta)}°")
    print(f"  z = ρ·e^(iθ) = {z}")
else:
    print(f"\n  u/(2ρ) = {u/(2*rho)} (real eigenvalues)")

# Test (1+2i) direction
print(f"\n=== Test z = c·(1+2i) hypothesis ===")
# z = c·(1+2i) implies z/(1+2i) is real
one_plus_2i = mpc(1, 2)
z_over = z / one_plus_2i
print(f"  z/(1+2i) = {z_over}")
print(f"  imag part / real part = {abs(z_over.imag)/abs(z_over.real) if abs(z_over.real)>0 else 'inf'}")
print(f"  → if very small, hypothesis confirmed; if O(1), refuted")

# Identify c via continued fraction
from fractions import Fraction
c_real = z_over.real
print(f"  c (real part) = {c_real}")
# Try to identify as rational
for max_denom in [100, 1000, 10000, 100000, 1000000]:
    f = Fraction(float(c_real)).limit_denominator(max_denom)
    approx = mpf(f.numerator)/mpf(f.denominator)
    diff = c_real - approx
    print(f"  c ≈ {f} (denom ≤ {max_denom}), diff = {float(diff):+.6e}")

# === L-value basis at full precision ===
print(f"\n=== Building L-value basis ===")
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

L_4_17 = L1_chi(build_chi(17, 3, 4, 16), 17)
L_2_17 = L1_chi(build_chi(17, 3, 8, 16), 17)
L_4_5  = L1_chi(build_chi(5, 2, 1, 4), 5)
L_2_5  = L1_chi(build_chi(5, 2, 2, 4), 5)
L_4_13 = L1_chi(build_chi(13, 2, 3, 12), 13)
L_8_17 = L1_chi(build_chi(17, 3, 2, 16), 17)
L_16_17 = L1_chi(build_chi(17, 3, 1, 16), 17)

print(f"  L(1, chi_4_17) = {L_4_17}")
print(f"  L(1, chi_2_17) = {L_2_17}")

# Targets at full precision
T1 = c_inf - c1   # depth-1 residual
T2 = c_inf - c0   # Δ_∞
T3 = c_inf        # c_∞

print(f"\n=== Targets at high precision ===")
print(f"  T1 = c_∞ - c(1)     = {T1}")
print(f"  T2 = Δ_∞ = c_∞ - c0 = {T2}")
print(f"  T3 = c_∞            = {T3}")

# === PSLQ at proper precision ===
basis = {
    '1': mpf(1),
    'Re L(chi_4_17)':  mpf(L_4_17.real),
    'Im L(chi_4_17)':  mpf(L_4_17.imag),
    'L(chi_2_17)':     mpf(L_2_17.real),
    'Re L(chi_4_5)':   mpf(L_4_5.real),
    'Im L(chi_4_5)':   mpf(L_4_5.imag),
    'L(chi_2_5)':      mpf(L_2_5.real),
    'Re L(chi_4_13)':  mpf(L_4_13.real),
    'Im L(chi_4_13)':  mpf(L_4_13.imag),
    'Re L(chi_8_17)':  mpf(L_8_17.real),
    'Im L(chi_8_17)':  mpf(L_8_17.imag),
    'Re L(chi_16_17)': mpf(L_16_17.real),
    'Im L(chi_16_17)': mpf(L_16_17.imag),
    '1/√17':           1/sqrt(mpf(17)),
    '1/√5':            1/sqrt(mpf(5)),
    'log(2)':          log(mpf(2)),
}
cand_names = list(basis.keys())
cand_vals = list(basis.values())

print(f"\n=== PSLQ at multiple tolerances ===")
for tol_exp in [10, 14, 18, 20, 22, 25]:
    tol = mpf(10) ** (-tol_exp)
    print(f"\n--- tol = 10^{-tol_exp}, maxcoeff = 10^6 ---")
    for tname, tval in [('T1', T1), ('T2', T2), ('T3', T3)]:
        rel = pslq([tval] + cand_vals, tol=tol, maxcoeff=10**6)
        if rel is None:
            print(f"  {tname}: no relation")
        else:
            terms = [f"({rel[0]})·{tname}"]
            for c, n in zip(rel[1:], cand_names):
                if c != 0:
                    terms.append(f"({c:+d})·{n}")
            residual = sum(c*v for c, v in zip(rel, [tval]+cand_vals))
            print(f"  {tname}: {' '.join(terms)}  [res={float(residual):.2e}]")
