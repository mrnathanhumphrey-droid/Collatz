"""
probe_hecke_highprec_2026_05_30.py

Compute Hecke L on Z[i] for quartic characters at:
  - π = 1+2i (norm 5)
  - π = 4+i  (norm 17)
  - π = 4-i  (norm 17, conjugate)
at high precision via direct sum at N_max = 50000.

Then PSLQ T1, T2, T3 against extended L-value basis.
"""
from __future__ import annotations
import sys, time
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40

q = 17

# Quartic residue symbol (α / π)_4 for various small Gaussian primes
def quartic_at_1_plus_2i(a, b):
    """π = 1 + 2i, N=5. iso Z[i]/π ≅ F_5 via i → 2.
    Lift: 1→1, 2→i, 3→-i, 4→-1."""
    r = (a + 2*b) % 5
    if r == 0: return None
    return [mpc(1,0), mpc(0,1), mpc(0,-1), mpc(-1,0)][r-1]

def quartic_at_4_plus_i(a, b):
    """π = 4 + i, N=17. iso Z[i]/π ≅ F_17 via i → -4 ≡ 13 (since 4+i ≡ 0 ⇒ i ≡ -4 mod π).
    Lift via discrete log: x = 3^k in F_17* → ψ(x) = i^k."""
    r = (a + 13*b) % 17
    if r == 0: return None
    # discrete log base 3 mod 17: table
    dlog3 = {1:0, 3:1, 9:2, 10:3, 13:4, 5:5, 15:6, 11:7, 16:8, 14:9, 8:10, 7:11, 4:12, 12:13, 2:14, 6:15}
    k = dlog3[r]
    # ψ(3^k) = i^k
    return [mpc(1,0), mpc(0,1), mpc(-1,0), mpc(0,-1)][k % 4]

def quartic_at_4_minus_i(a, b):
    """π = 4 - i, N=17. iso Z[i]/π ≅ F_17 via i → 4 mod 17."""
    r = (a + 4*b) % 17
    if r == 0: return None
    dlog3 = {1:0, 3:1, 9:2, 10:3, 13:4, 5:5, 15:6, 11:7, 16:8, 14:9, 8:10, 7:11, 4:12, 12:13, 2:14, 6:15}
    k = dlog3[r]
    return [mpc(1,0), mpc(0,1), mpc(-1,0), mpc(0,-1)][k % 4]

def hecke_L1(sym_fn, N_max):
    """Compute L(1, ψ) = sum over Gaussian integers α in upper half + positive real, ψ(α)/N(α)."""
    total = mpc(0)
    for a_int in range(-N_max, N_max + 1):
        for b_int in range(1, N_max + 1):
            nrm = a_int*a_int + b_int*b_int
            if nrm == 0 or nrm > N_max: continue
            s = sym_fn(a_int, b_int)
            if s is None: continue
            total += s / mpf(nrm)
    for a_int in range(1, N_max + 1):
        s = sym_fn(a_int, 0)
        if s is None: continue
        total += s / mpf(a_int*a_int)
    return total

print(f"Computing Hecke L at N_max = 30000 (expect ~25 digits convergence)...")
t0 = time.time()
N_max = 30000

print(f"  L(1, ψ_4 mod 1+2i)...")
t1 = time.time()
L_h_5 = hecke_L1(quartic_at_1_plus_2i, N_max)
print(f"    = {L_h_5}  ({time.time()-t1:.1f}s)")

print(f"  L(1, ψ_4 mod 4+i)...")
t1 = time.time()
L_h_17a = hecke_L1(quartic_at_4_plus_i, N_max)
print(f"    = {L_h_17a}  ({time.time()-t1:.1f}s)")

print(f"  L(1, ψ_4 mod 4-i)...")
t1 = time.time()
L_h_17b = hecke_L1(quartic_at_4_minus_i, N_max)
print(f"    = {L_h_17b}  ({time.time()-t1:.1f}s)")

print(f"\nTotal Hecke compute: {time.time()-t0:.1f}s")

# === Use the just-computed high-precision c_∞ ===
c0 = mpf(19)/mpf(127)
c1_num = 265011804960406635465672455997699
c1_den = 1730087916969634762193659498034425
c1 = mpf(c1_num)/mpf(c1_den)
c_inf = mpf("0.15298912060588505582178105855802487428330184973053")
T1 = c_inf - c1
T2 = c_inf - c0
T3 = c_inf

# === Dirichlet L values (for completeness) ===
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

L_d_4_17 = L1_chi(build_chi(17, 3, 4, 16), 17)
L_d_4_5  = L1_chi(build_chi(5, 2, 1, 4), 5)
L_d_8_17 = L1_chi(build_chi(17, 3, 2, 16), 17)

# Compare Hecke and Dirichlet
print(f"\n=== Hecke vs Dirichlet (for same character via iso) ===")
print(f"L(1, ψ_4 mod 1+2i)   = {L_h_5}")
print(f"L(1, χ_4 mod 5)      = {L_d_4_5}")
print(f"  ratio L_hecke / L_dirich = {L_h_5 / L_d_4_5}")
print(f"L(1, ψ_4 mod 4+i)    = {L_h_17a}")
print(f"L(1, ψ_4 mod 4-i)    = {L_h_17b}")
print(f"L(1, χ_4 mod 17)     = {L_d_4_17}")
print(f"  L_h(4+i) + L_h(4-i) = {L_h_17a + L_h_17b}")
print(f"  L_h(4+i) - L_h(4-i) = {L_h_17a - L_h_17b}")

# Specific test: is c_∞ related to L_h_5 by clean factor?
print(f"\n=== Test specific ratios ===")
for tname, tval in [('T1', T1), ('T2', T2), ('T3', T3)]:
    print(f"\n  {tname} = {tval}")
    for lname, lval in [('Re L_h(1+2i)', mpf(L_h_5.real)),
                         ('Im L_h(1+2i)', mpf(L_h_5.imag)),
                         ('|L_h(1+2i)|',  abs(L_h_5)),
                         ('Re L_h(4+i)',  mpf(L_h_17a.real)),
                         ('Im L_h(4+i)',  mpf(L_h_17a.imag)),
                         ('|L_h(4+i)|',   abs(L_h_17a))]:
        if abs(lval) > mpf(10)**(-30):
            ratio = tval / lval
            f = Fraction(float(ratio)).limit_denominator(10000)
            approx = mpf(f.numerator)/mpf(f.denominator)
            diff = ratio - approx
            if abs(diff) < mpf(10)**(-6):
                print(f"    {tname}/{lname:18} = {float(ratio):+.10e} ≈ {f}  (diff {float(diff):+.2e})")

# === Full PSLQ ===
basis = {
    '1': mpf(1),
    'Re L_d_4_17': mpf(L_d_4_17.real),
    'Im L_d_4_17': mpf(L_d_4_17.imag),
    'Re L_d_4_5':  mpf(L_d_4_5.real),
    'Im L_d_4_5':  mpf(L_d_4_5.imag),
    'L_d_8_17 (real part)': mpf(L_d_8_17.real),
    'L_d_8_17 (imag part)': mpf(L_d_8_17.imag),
    'Re L_h(1+2i)': mpf(L_h_5.real),
    'Im L_h(1+2i)': mpf(L_h_5.imag),
    'Re L_h(4+i)':  mpf(L_h_17a.real),
    'Im L_h(4+i)':  mpf(L_h_17a.imag),
    'Re L_h(4-i)':  mpf(L_h_17b.real),
    'Im L_h(4-i)':  mpf(L_h_17b.imag),
    'log(4+√17)/√17': log(mpf(4)+sqrt(mpf(17)))/sqrt(mpf(17)),
    '1/√17': mpf(1)/sqrt(mpf(17)),
    '1/√5': mpf(1)/sqrt(mpf(5)),
}
cand_names = list(basis.keys())
cand_vals = list(basis.values())

print(f"\n=== PSLQ at high precision (Hecke L's now at ~25 digits) ===")
for tol_exp in [10, 12, 14, 16, 18, 20, 22]:
    tol = mpf(10) ** (-tol_exp)
    print(f"\n--- tol = 10^{-tol_exp}, maxcoeff = 10^6 ---")
    for tname, tval in [('T1', T1), ('T2', T2), ('T3', T3)]:
        vec = [tval] + cand_vals
        rel = pslq(vec, tol=tol, maxcoeff=10**6)
        if rel is None:
            print(f"  {tname}: no relation")
        else:
            terms = []
            if rel[0] != 0:
                terms.append(f"({rel[0]})·{tname}")
            for c, name in zip(rel[1:], cand_names):
                if c != 0:
                    terms.append(f"({c:+d})·{name}")
            residual = sum(c * v for c, v in zip(rel, vec))
            print(f"  {tname}: {' + '.join(terms)}  [residual {float(residual):.2e}]")
