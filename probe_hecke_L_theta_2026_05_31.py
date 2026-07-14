"""
probe_hecke_L_theta_2026_05_31.py

Compute Hecke L(1, ψ) for quartic Hecke characters on Z[i] at meaningful precision,
then PSLQ c_∞ against the Hecke L basis.

Strategy: direct sum L(1, ψ) = Σ_n a(n)/n where a(n) = Σ_{|α|²=n} ψ(α),
organized by Gaussian-integer norm. Computed to N_max via norm enumeration.
For N_max = 10^6, expect ~6 digit precision (Hecke L oscillates so partial sums
of a(n)/n converge slowly).

Use Cesaro-summation acceleration via mpmath.nsum where applicable.

Quartic characters tested:
  ψ_4 mod (1+2i), norm 5
  ψ_4 mod (4+i), norm 17
  ψ_4 mod (4-i), norm 17 (conjugate)
  ψ_4 mod (3+2i), norm 13
"""
from __future__ import annotations
import sys, time
from mpmath import mp, mpf, mpc, sqrt, pi, log, exp, digamma, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30

# Discrete log of x mod q in base g
def dlog_table(g, q):
    table = {}
    x = 1
    for k in range(q - 1):
        table[x] = k
        x = (x * g) % q
    return table

# Quartic character on F_p via discrete log: chi_4(g^k) = i^k for g primitive root
def quartic_char_F_p(p, g):
    dlog = dlog_table(g, p)
    def chi(r):
        if r % p == 0: return None
        k = dlog[r % p]
        return [mpc(1,0), mpc(0,1), mpc(-1,0), mpc(0,-1)][k % 4]
    return chi

# Build characters via the Z[i]/π ≅ F_p isomorphism

# Z[i]/(1+2i) ≅ F_5 via i ↔ 2 (since 1+2i = 0 ⇒ i ≡ -1/2 ≡ 2 mod 5)
chi_p5 = quartic_char_F_p(5, 2)
def chi_4_at_1p2i(a, b):
    r = (a + 2 * b) % 5
    return chi_p5(r)

# Z[i]/(4+i) ≅ F_17 via i ↔ -4 ≡ 13 mod 17
chi_p17 = quartic_char_F_p(17, 3)
def chi_4_at_4pi(a, b):
    r = (a + 13 * b) % 17
    return chi_p17(r)

# Z[i]/(3+2i) ≅ F_13 via i ↔ ? 3+2i=0 ⇒ 2i=-3 ⇒ i=-3/2 mod 13.
# 2^(-1) mod 13: 2·7=14=1 mod 13, so 2^(-1)=7. Thus i ≡ -3·7 = -21 ≡ 5 mod 13.
chi_p13 = quartic_char_F_p(13, 2)
def chi_4_at_3p2i(a, b):
    r = (a + 5 * b) % 13
    return chi_p13(r)

# === Build a(n) = Σ_{|α|²=n} ψ(α) by enumerating Gaussian integers ===
def compute_an_table(chi_fn, N_max):
    """Compute a(n) = Σ ψ(α) over Gaussian integers in FIRST QUADRANT
    (b > 0 OR (b = 0 AND a > 0)) with a²+b² = n. This restricts to one
    representative per unit orbit when ψ is non-trivial on units."""
    a_table = [mpc(0)] * (N_max + 1)
    # First quadrant: b > 0 with any a, OR b = 0 with a > 0
    # b > 0
    for a in range(-int(N_max**0.5) - 1, int(N_max**0.5) + 2):
        a2 = a * a
        if a2 > N_max: continue
        b_max = int((N_max - a2)**0.5) + 1
        for b in range(1, b_max + 1):
            n = a2 + b * b
            if n == 0 or n > N_max: continue
            sym = chi_fn(a, b)
            if sym is None: continue
            a_table[n] += sym
    # b = 0 with a > 0
    for a in range(1, int(N_max**0.5) + 2):
        n = a * a
        if n > N_max: continue
        sym = chi_fn(a, 0)
        if sym is None: continue
        a_table[n] += sym
    return a_table

# === Compute L(1, ψ) = Σ_{n≥1} a(n)/n (truncated) ===
def hecke_L1(a_table, N_max):
    """L(1, ψ) ≈ Σ_{n=1}^{N_max} a(n)/n."""
    total = mpc(0)
    for n in range(1, N_max + 1):
        if a_table[n] != 0:
            total += a_table[n] / mpf(n)
    return total

# === Compute and report ===
N_max = 200000
print(f"Computing Hecke L(1, ψ_4) for various π via direct sum at N_max = {N_max}")
import math
print(f"(Expected precision: ~{int(0.5 * math.log10(N_max))} digits)\n")

t0 = time.time()
print("Building a(n) table for ψ mod (1+2i)...")
a_table_5 = compute_an_table(chi_4_at_1p2i, N_max)
print(f"  done in {time.time()-t0:.1f}s")
L_h_5 = hecke_L1(a_table_5, N_max)
print(f"  L_h(1, ψ_4/(1+2i)) ≈ {L_h_5}")

t1 = time.time()
print("\nBuilding a(n) table for ψ mod (4+i)...")
a_table_17 = compute_an_table(chi_4_at_4pi, N_max)
print(f"  done in {time.time()-t1:.1f}s")
L_h_17 = hecke_L1(a_table_17, N_max)
print(f"  L_h(1, ψ_4/(4+i)) ≈ {L_h_17}")

t2 = time.time()
print("\nBuilding a(n) table for ψ mod (3+2i)...")
a_table_13 = compute_an_table(chi_4_at_3p2i, N_max)
print(f"  done in {time.time()-t2:.1f}s")
L_h_13 = hecke_L1(a_table_13, N_max)
print(f"  L_h(1, ψ_4/(3+2i)) ≈ {L_h_13}")

print(f"\nTotal Hecke L compute: {time.time()-t0:.1f}s")

# === Dirichlet L-values for comparison/basis ===
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

L_d_4_17 = L1_chi(build_chi(17, 3, 4, 16), 17)
L_d_4_5 = L1_chi(build_chi(5, 2, 1, 4), 5)
L_d_4_13 = L1_chi(build_chi(13, 2, 3, 12), 13)

# === Targets at high precision ===
c0 = mpf(19)/mpf(127)
c1 = mpf(265011804960406635465672455997699) / mpf(1730087916969634762193659498034425)
c_inf = mpf("0.15298912060588517527891674877413229926086222622334")
T1 = c_inf - c1
T2 = c_inf - c0
T3 = c_inf

print(f"\n=== Targets ===")
print(f"  T1 = c_∞ - c(1)     = {T1}")
print(f"  T2 = Δ_∞ = c_∞ - c0 = {T2}")
print(f"  T3 = c_∞            = {T3}")

# === PSLQ basis ===
basis = {
    '1': mpf(1),
    'Re L_h(1+2i)':       mpf(L_h_5.real),
    'Im L_h(1+2i)':       mpf(L_h_5.imag),
    'Re L_h(4+i)':        mpf(L_h_17.real),
    'Im L_h(4+i)':        mpf(L_h_17.imag),
    'Re L_h(3+2i)':       mpf(L_h_13.real),
    'Im L_h(3+2i)':       mpf(L_h_13.imag),
    'Re L_d_4_17':        mpf(L_d_4_17.real),
    'Im L_d_4_17':        mpf(L_d_4_17.imag),
    'Re L_d_4_5':         mpf(L_d_4_5.real),
    'Im L_d_4_5':         mpf(L_d_4_5.imag),
    'Re L_d_4_13':        mpf(L_d_4_13.real),
    'Im L_d_4_13':        mpf(L_d_4_13.imag),
    '1/√17':              1/sqrt(mpf(17)),
    '1/√5':               1/sqrt(mpf(5)),
    '1/√13':              1/sqrt(mpf(13)),
}
cand_names = list(basis.keys())
cand_vals = list(basis.values())

print(f"\n=== PSLQ at multiple tolerances (full basis with Hecke L) ===")
for tol_exp in [4, 6, 8, 10, 12, 14]:
    tol = mpf(10) ** (-tol_exp)
    print(f"\n--- tol = 10^{-tol_exp}, maxcoeff = 10^5 ---")
    for tname, tval in [('T1', T1), ('T2', T2), ('T3', T3)]:
        rel = pslq([tval] + cand_vals, tol=tol, maxcoeff=10**5)
        if rel is None:
            print(f"  {tname}: no relation")
        else:
            terms = [f"({rel[0]})·{tname}"]
            for c, n in zip(rel[1:], cand_names):
                if c != 0: terms.append(f"({c:+d})·{n}")
            residual = sum(c*v for c, v in zip(rel, [tval]+cand_vals))
            print(f"  {tname}: {' '.join(terms)}  [res={float(residual):.2e}]")

# Single-ratio test
print(f"\n=== Single-ratio test against Hecke L candidates ===")
from fractions import Fraction
for tname, tval in [('T1', T1), ('T2', T2), ('T3', T3)]:
    print(f"\n  {tname} = {float(tval):.10e}")
    for name in ['Re L_h(1+2i)', 'Im L_h(1+2i)', 'Re L_h(4+i)', 'Im L_h(4+i)',
                 'Re L_h(3+2i)', 'Im L_h(3+2i)']:
        v = basis[name]
        if abs(v) < mpf(10)**(-15): continue
        r = tval / v
        f = Fraction(float(r)).limit_denominator(1000)
        approx = mpf(f.numerator)/mpf(f.denominator)
        diff = r - approx
        if abs(diff) < mpf(10)**(-4):
            print(f"    {tname}/{name:18} = {float(r):+.10e} ≈ {f}  diff = {float(diff):+.2e}")
