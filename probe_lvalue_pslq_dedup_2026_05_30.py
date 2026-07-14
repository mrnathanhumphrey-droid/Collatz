"""
probe_lvalue_pslq_dedup_2026_05_30.py

Re-run PSLQ with:
1. Dedup basis: drop log(4+√17)/√17 since L(1,χ_2_17)/2 = ln(4+√17)/√17 exactly.
2. Realistic tolerance matching the c_∞ actual precision.

The c_∞ was extracted at HIGH RAW digits (40 printed) but the EFFECTIVE precision
is limited by Δ_4, Δ_5 (float64 ~ 12 digits) amplified by condition number 5e6.
So c_∞ trust ~ 7 digits.

Use tol = 10^-6, maxcoeff = 10^4 to find true relations only.
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40
q = 17

# Re-use the just-computed high-precision values
c_inf = mpf("0.1529891206058850808848604387217388211109")
A_asymp = mpf("0.003382821393286655688010045020951419536146")
c0 = mpf(19)/mpf(127)
T1 = mpf("-0.0001891094488005902410869256482049402539196")
T2 = A_asymp
T3 = c_inf

# === L-values ===
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

# Hecke L on Z[i] mod (1+2i) — recompute at higher N_max for ~5-digit precision
print("Computing Hecke L mod (1+2i) at N_max=5000...")
def hecke_at_1plus2i(a, b):
    r = (a + 2*b) % 5
    if r == 0: return None
    if r == 1: return mpc(1, 0)
    if r == 2: return mpc(0, 1)
    if r == 3: return mpc(0, -1)
    if r == 4: return mpc(-1, 0)

N_max = 5000
hecke_L = mpc(0)
for a_int in range(-N_max, N_max + 1):
    for b_int in range(1, N_max + 1):
        nrm = a_int*a_int + b_int*b_int
        if nrm == 0 or nrm > N_max:
            continue
        sym = hecke_at_1plus2i(a_int, b_int)
        if sym is None: continue
        hecke_L += sym / mpf(nrm)
for a_int in range(1, N_max + 1):
    sym = hecke_at_1plus2i(a_int, 0)
    if sym is None: continue
    hecke_L += sym / mpf(a_int * a_int)
print(f"  Hecke L = {hecke_L}")

# === DEDUPED basis ===
basis = {
    '1': mpf(1),
    'Re L(1,chi_4_17)':  mpf(L_4_17.real),
    'Im L(1,chi_4_17)':  mpf(L_4_17.imag),
    'L(1,chi_2_17)':     mpf(L_2_17.real),   # note: this is 2 × log(4+√17)/√17
    'Re L(1,chi_4_5)':   mpf(L_4_5.real),
    'Im L(1,chi_4_5)':   mpf(L_4_5.imag),
    'L(1,chi_2_5)':      mpf(L_2_5.real),
    'Re L(1,chi_4_13)':  mpf(L_4_13.real),
    'Im L(1,chi_4_13)':  mpf(L_4_13.imag),
    'Re L(1,chi_8_17)':  mpf(L_8_17.real),
    'Im L(1,chi_8_17)':  mpf(L_8_17.imag),
    'Re Hecke L(1, ψ_4/(1+2i))': mpf(hecke_L.real),
    'Im Hecke L(1, ψ_4/(1+2i))': mpf(hecke_L.imag),
    '1/√17':  mpf(1)/sqrt(mpf(17)),
    '1/√5':   mpf(1)/sqrt(mpf(5)),
    'π/√17':  pi/sqrt(mpf(17)),
    'π/√5':   pi/sqrt(mpf(5)),
    'log(2)': log(mpf(2)),
}

cand_names = list(basis.keys())
cand_vals = list(basis.values())

# === PSLQ with REALISTIC tolerance ===
print(f"\n=== PSLQ at realistic tolerance ===")
# c_∞'s true precision: limited by float64 Δ_4, Δ_5 amplified by condition number 5e6 → ~7 digits.
# So search for relations with residual < 10^-6, max integer coefficient < 10^5.
for tol_exp in [6, 8, 10, 12]:
    tol = mpf(10) ** (-tol_exp)
    print(f"\n--- tol = 10^{-tol_exp}, maxcoeff = 10^5 ---")
    for target_name, target_val in [('T1', T1), ('T2', T2), ('T3', T3)]:
        vec = [target_val] + cand_vals
        rel = pslq(vec, tol=tol, maxcoeff=10**5)
        if rel is None:
            print(f"  {target_name}: no relation")
        else:
            terms = []
            if rel[0] != 0:
                terms.append(f"({rel[0]})·{target_name}")
            for c, name in zip(rel[1:], cand_names):
                if c != 0:
                    terms.append(f"({c:+d})·{name}")
            residual = sum(c * v for c, v in zip(rel, vec))
            print(f"  {target_name}: {' + '.join(terms)}  [residual {float(residual):.2e}]")

# === Specifically test if T1 (or T2, T3) is rational ×L for individual L values ===
print(f"\n=== Single-candidate rational ratio test (target = q · candidate, q rational, denom ≤ 100) ===")
from fractions import Fraction
for target_name, target_val in [('T1', T1), ('T2', T2), ('T3', T3)]:
    print(f"\n--- {target_name} = {float(target_val):.10e} ---")
    for name, val in basis.items():
        if abs(val) < mpf(10)**(-30): continue
        ratio = target_val / val
        f = Fraction(float(ratio)).limit_denominator(100)
        approx = mpf(f.numerator) / mpf(f.denominator)
        diff = ratio - approx
        if abs(diff) < mpf(10)**(-5):
            print(f"  {name:30}: ratio = {float(ratio):+.8f} ≈ {f} (diff {float(diff):+.2e})")

# === Test: is c_∞ equal to a specific combination involving 19/127 and one L-value? ===
print(f"\n=== Test: c_∞ = 19/127 + q · L for small rational q ===")
for name, val in basis.items():
    if abs(val) < mpf(10)**(-30): continue
    delta_target = T3 - mpf(19)/mpf(127)  # = T2
    ratio = delta_target / val
    f = Fraction(float(ratio)).limit_denominator(1000)
    approx = mpf(f.numerator) / mpf(f.denominator)
    diff = ratio - approx
    if abs(diff) < mpf(10)**(-5):
        print(f"  19/127 + ({f})·{name}  diff {float(diff):+.2e}")
