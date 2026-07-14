"""
probe_pslq_with_c4_2026_05_30.py

Re-solve the damped osc recurrence with Δ_1 exact, Δ_2, Δ_3, Δ_4 at high precision,
Δ_5 still at float64 (12 digits). Now ONLY rhs[2]=Δ_5 contaminates the system.

Then PSLQ at improved tolerance against L-value basis.
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 50
q = 17

# Exact c(0), c(1)
c0_frac = Fraction(19, 127)
c1_frac = Fraction(265011804960406635465672455997699, 1730087916969634762193659498034425)
c0 = mpf(c0_frac.numerator) / mpf(c0_frac.denominator)
c1 = mpf(c1_frac.numerator) / mpf(c1_frac.denominator)

# High-precision c(2), c(3), c(4)
c2 = mpf("0.1532479207790874367716967397510852719979")
c3 = mpf("0.1530053316915140267018082388619123587748")
c4 = mpf("0.1529887090468211652268046601150049898549764588045")

# Float64 c(5)
c5 = mpf("0.1529889994135218")

print("=== c(m) at highest available precision ===")
print(f"c(0) = {c0}  (exact)")
print(f"c(1) = {c1}  (exact)")
print(f"c(2) = {c2}  (mpmath dps=40)")
print(f"c(3) = {c3}  (mpmath dps=40)")
print(f"c(4) = {c4}  (mpmath dps=50)")
print(f"c(5) = {c5}  (float64 ~12 digits)")

# Δ_m
deltas = {m: mpf(0) for m in range(6)}
deltas[1] = c1 - c0
deltas[2] = c2 - c0
deltas[3] = c3 - c0
deltas[4] = c4 - c0
deltas[5] = c5 - c0

print(f"\n=== Δ_m ===")
for m in range(1, 6):
    print(f"  Δ_{m} = {deltas[m]}")

# === Solve linear recurrence with mixed precision ===
import mpmath as mpmath
M = mpmath.matrix([
    [deltas[2], -deltas[1], mpf(1)],
    [deltas[3], -deltas[2], mpf(1)],
    [deltas[4], -deltas[3], mpf(1)],
])
rhs = mpmath.matrix([deltas[3], deltas[4], deltas[5]])
sol = mpmath.lu_solve(M, rhs)
u, v, K = sol[0], sol[1], sol[2]

print(f"\n=== Recurrence fit ===")
print(f"  u = {u}")
print(f"  v = {v}")
print(f"  K = {K}")

A = K / (1 - u + v)
print(f"  A = {A}")
c_inf_v1 = c0 + A
print(f"\nc_∞ (using Δ_5 float64) = {c_inf_v1}")

# === Alternative: solve at higher precision by NOT using Δ_5 ===
# We have 2 high-precision equations (eqs A, B from m=1, m=2):
#   eq A: u Δ_2 - v Δ_1 + K = Δ_3
#   eq B: u Δ_3 - v Δ_2 + K = Δ_4
# Both fully high-precision. 2 equations, 3 unknowns (u, v, K).
# Use the constraint A = c_∞ - c(0). We don't know A yet.
#
# Actually with 3 unknowns and 2 equations, we have a 1-parameter family of solutions.
# Parameterize by v. Then u and K are determined.
# For each v, compute A = K/(1-u+v) and check consistency with the EMPIRICAL TREND.
#
# Better: use Aitken extrapolation on c(2), c(3), c(4) (all high precision).
print(f"\n=== Aitken Δ² on c(2), c(3), c(4) ===")
# c_∞ ≈ c(4) - (c(4) - c(3))² / (c(4) - 2 c(3) + c(2))
num_aitken = (c4 - c3) ** 2
den_aitken = c4 - 2*c3 + c2
c_inf_aitken = c4 - num_aitken / den_aitken
print(f"  c_∞ (Aitken c(2),c(3),c(4)) = {c_inf_aitken}")

# Compare to v1
print(f"  c_∞ (recurrence) - c_∞ (Aitken) = {c_inf_v1 - c_inf_aitken}")

# === The TRUE precision after fix ===
# Δ_5 has ~12 digits. It enters rhs[2] only. The amplification factor through the inverse matrix is the
# (3,2) entry of M^{-1} multiplied by Δ_5's error.
# Compute M^{-1} explicitly
Minv = mpmath.matrix(M ** -1)
print(f"\n=== M^{-1} (for error propagation analysis) ===")
print(Minv)
# u_error ≈ |Minv[0,2]| × Δ_5_error
print(f"\n  Δ_5 float64 error ≈ 1e-13.")
print(f"  u error ≈ |M^-1[0,2]| × 1e-13 = {abs(Minv[0,2]) * mpf('1e-13')}")
print(f"  v error ≈ |M^-1[1,2]| × 1e-13 = {abs(Minv[1,2]) * mpf('1e-13')}")
print(f"  K error ≈ |M^-1[2,2]| × 1e-13 = {abs(Minv[2,2]) * mpf('1e-13')}")
# A = K / (1-u+v), so A_error ~ K_error / |1-u+v| + ... (dominant from K, u, v contributions cancel
# approximately).

# === Target c_∞ ===
print(f"\n=== Targets at the FIT precision ===")
T1 = c_inf_v1 - c1   # depth-1 residual
T2 = c_inf_v1 - c0   # Δ_∞
T3 = c_inf_v1        # c_∞
print(f"  T1 = c_∞ - c(1)   = {T1}")
print(f"  T2 = c_∞ - 19/127 = {T2}")
print(f"  T3 = c_∞          = {T3}")
print(f"  c_∞_aitken        = {c_inf_aitken}")
print(f"  diff (recurr - aitken) = {c_inf_v1 - c_inf_aitken}")

# === Improved L-value basis (drop redundant log(4+√17)/√17) ===
def build_chi(q, g, r, phi_q):
    table = [mpc(0)] * q
    x = 1
    for k in range(phi_q):
        table[x] = exp(2 * pi * mpc(0, 1) * r * k / phi_q)
        x = (x * g) % q
    return table

def L1_chi(chi_table, q):
    total = mpc(0)
    for a in range(1, q):
        total += chi_table[a] * digamma(mpf(a)/mpf(q))
    return -total / mpf(q)

L_4_17 = L1_chi(build_chi(17, 3, 4, 16), 17)
L_2_17 = L1_chi(build_chi(17, 3, 8, 16), 17)
L_4_5  = L1_chi(build_chi(5, 2, 1, 4), 5)
L_2_5  = L1_chi(build_chi(5, 2, 2, 4), 5)
L_4_13 = L1_chi(build_chi(13, 2, 3, 12), 13)
L_8_17 = L1_chi(build_chi(17, 3, 2, 16), 17)
L_16_17 = L1_chi(build_chi(17, 3, 1, 16), 17)

# Hecke L on Z[i] mod (1+2i)
print(f"\nComputing Hecke L at N_max=10000 for better precision...")
def hecke_at_1plus2i(a, b):
    r = (a + 2*b) % 5
    if r == 0: return None
    return [mpc(1,0), mpc(0,1), mpc(0,-1), mpc(-1,0)][r-1]
N_max = 10000
hecke_L = mpc(0)
for a_int in range(-N_max, N_max + 1):
    for b_int in range(1, N_max + 1):
        nrm = a_int*a_int + b_int*b_int
        if nrm == 0 or nrm > N_max: continue
        sym = hecke_at_1plus2i(a_int, b_int)
        if sym is None: continue
        hecke_L += sym / mpf(nrm)
for a_int in range(1, N_max + 1):
    sym = hecke_at_1plus2i(a_int, 0)
    if sym is None: continue
    hecke_L += sym / mpf(a_int * a_int)

basis = {
    '1': mpf(1),
    'Re L(1,chi_4_17)': mpf(L_4_17.real),
    'Im L(1,chi_4_17)': mpf(L_4_17.imag),
    'L(1,chi_2_17)':    mpf(L_2_17.real),
    'Re L(1,chi_4_5)':  mpf(L_4_5.real),
    'Im L(1,chi_4_5)':  mpf(L_4_5.imag),
    'L(1,chi_2_5)':     mpf(L_2_5.real),
    'Re L(1,chi_4_13)': mpf(L_4_13.real),
    'Im L(1,chi_4_13)': mpf(L_4_13.imag),
    'Re L(1,chi_8_17)': mpf(L_8_17.real),
    'Im L(1,chi_8_17)': mpf(L_8_17.imag),
    'Re L(1,chi_16_17)': mpf(L_16_17.real),
    'Im L(1,chi_16_17)': mpf(L_16_17.imag),
    'Re Hecke L(1, ψ_4/(1+2i))': mpf(hecke_L.real),
    'Im Hecke L(1, ψ_4/(1+2i))': mpf(hecke_L.imag),
    '1/√17': mpf(1)/sqrt(mpf(17)),
    '1/√5':  mpf(1)/sqrt(mpf(5)),
}

cand_names = list(basis.keys())
cand_vals = list(basis.values())

# Use Aitken-based c_∞ for PSLQ (cleaner, fewer assumptions about model)
print(f"\n=== PSLQ at multiple tolerances ===")
for cinf_label, cinf_val in [("recurrence", c_inf_v1), ("Aitken", c_inf_aitken)]:
    print(f"\n*** Using c_∞ from {cinf_label} = {cinf_val} ***")
    T1_x = cinf_val - c1
    T2_x = cinf_val - c0
    T3_x = cinf_val
    for tol_exp in [10, 12, 14, 16, 18, 20]:
        tol = mpf(10) ** (-tol_exp)
        print(f"  --- tol = 10^{-tol_exp} ---")
        for tname, tval in [('T1', T1_x), ('T2', T2_x), ('T3', T3_x)]:
            vec = [tval] + cand_vals
            rel = pslq(vec, tol=tol, maxcoeff=10**5)
            if rel is not None:
                terms = []
                if rel[0] != 0:
                    terms.append(f"({rel[0]})·{tname}")
                for c, name in zip(rel[1:], cand_names):
                    if c != 0:
                        terms.append(f"({c:+d})·{name}")
                residual = sum(c * v for c, v in zip(rel, vec))
                print(f"    {tname}: {' + '.join(terms)}  [residual {float(residual):.2e}]")
            else:
                print(f"    {tname}: no relation")
