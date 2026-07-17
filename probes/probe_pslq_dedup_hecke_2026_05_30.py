"""
probe_pslq_dedup_hecke_2026_05_30.py

Re-run PSLQ after removing the redundant L_h(4-i) candidate (= L_h(4+i) by symmetry).
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40

# Pre-computed values
c_inf = mpf("0.15298912060588505582178105855802487428330184973053")
c1 = mpf(265011804960406635465672455997699) / mpf(1730087916969634762193659498034425)
c0 = mpf(19) / mpf(127)
T1 = c_inf - c1
T2 = c_inf - c0
T3 = c_inf

# Hecke L values from the previous run
L_h_5 = mpc("1.142154141932535297369304399351445366565", "0.5931873768458615026860096785063044015615")
L_h_17 = mpc("1.307710470774823380971794982090557298625", "0.5015833597976616484339861454768866896701")

# Dirichlet L values
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
L_d_4_13 = L1_chi(build_chi(13, 2, 3, 12), 13)
L_d_8_17 = L1_chi(build_chi(17, 3, 2, 16), 17)
L_d_16_17 = L1_chi(build_chi(17, 3, 1, 16), 17)

# DEDUPED basis (drop one of the redundant Hecke L's; drop log(4+√17)/√17 = L(1,χ_2_17)/2)
basis = {
    '1':                  mpf(1),
    'Re L_d(chi_4_17)':   mpf(L_d_4_17.real),
    'Im L_d(chi_4_17)':   mpf(L_d_4_17.imag),
    'Re L_d(chi_4_5)':    mpf(L_d_4_5.real),
    'Im L_d(chi_4_5)':    mpf(L_d_4_5.imag),
    'Re L_d(chi_4_13)':   mpf(L_d_4_13.real),
    'Im L_d(chi_4_13)':   mpf(L_d_4_13.imag),
    'Re L_d(chi_8_17)':   mpf(L_d_8_17.real),
    'Im L_d(chi_8_17)':   mpf(L_d_8_17.imag),
    'Re L_d(chi_16_17)':  mpf(L_d_16_17.real),
    'Im L_d(chi_16_17)':  mpf(L_d_16_17.imag),
    'Re L_h(1+2i)':       mpf(L_h_5.real),
    'Im L_h(1+2i)':       mpf(L_h_5.imag),
    'Re L_h(4+i)':        mpf(L_h_17.real),
    'Im L_h(4+i)':        mpf(L_h_17.imag),
    '1/√17':              mpf(1)/sqrt(mpf(17)),
    '1/√5':               mpf(1)/sqrt(mpf(5)),
    'log(4+√17)/√17':     log(mpf(4)+sqrt(mpf(17)))/sqrt(mpf(17)),
}

cand_names = list(basis.keys())
cand_vals = list(basis.values())

print(f"=== PSLQ deduped basis (size {len(basis)}) ===")
for tol_exp in [9, 10, 11, 12, 14, 16]:
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

# Also: separate small-basis PSLQ — just T against the four "primary suspect" Hecke / quartic L's
print(f"\n=== Small-basis PSLQ ===")
small_basis = [
    ('1',                mpf(1)),
    ('Re L_h(1+2i)',     mpf(L_h_5.real)),
    ('Im L_h(1+2i)',     mpf(L_h_5.imag)),
    ('Re L_h(4+i)',      mpf(L_h_17.real)),
    ('Im L_h(4+i)',      mpf(L_h_17.imag)),
    ('Re L_d(chi_4_17)', mpf(L_d_4_17.real)),
    ('Im L_d(chi_4_17)', mpf(L_d_4_17.imag)),
]
small_names = [n for n, v in small_basis]
small_vals = [v for n, v in small_basis]

for tol_exp in [9, 10, 11, 12]:
    tol = mpf(10) ** (-tol_exp)
    print(f"\n--- small basis tol = 10^{-tol_exp}, maxcoeff = 10^5 ---")
    for tname, tval in [('T1', T1), ('T2', T2), ('T3', T3)]:
        rel = pslq([tval] + small_vals, tol=tol, maxcoeff=10**5)
        if rel is None:
            print(f"  {tname}: no relation")
        else:
            terms = [f"({rel[0]})·{tname}"]
            for c, name in zip(rel[1:], small_names):
                if c != 0: terms.append(f"({c:+d})·{name}")
            residual = sum(c * v for c, v in zip(rel, [tval]+small_vals))
            print(f"  {tname}: {' + '.join(terms)}  [residual {float(residual):.2e}]")
