"""
PROBE CROSSING -- test Wilson's two-monotone-mode model against the long-period model.
Model: Lam_r ~ a*rho1^r - b*rho2^r (two monotone modes, opposite sign) + weak m=2,4 oscillatory correction.
A difference of two decaying exponentials crosses zero EXACTLY ONCE -> negative for all r beyond, never returns.

KEY: Lam_r = sum_m 4^{-m} A_r(m), and A_r(m)=gamma_r(tau_m)-gamma_{r-1}(tau_m) is computable from the nu measure
at ANY r (build_nu), NOT limited by the eps-data (which stopped at eps_12 => Lam_11). So we can get Lam_12,13,14.

DECISIVE: does Lam cross once near r=12 and STAY negative (two-mono model) or return positive (oscillation, KILL)?
MONO_r = 0.25 A_r(1) + (1/64) A_r(3);  OSC_r = 0.0625 A_r(2) + (1/256) A_r(4)   [Lam = MONO+OSC+tail]
"""
import os, sys, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_gapop_R28 import build_nu

R_MAX = 16
MMAX = 10                                       # sum m=1..10 (4^-m weights); tail m>10 ~ 1e-8, negligible
HERE = os.path.dirname(os.path.abspath(__file__))
EPSJSON = os.path.join(HERE, '..', 'experiments_output', 'result_77_7_eps_exact_through_k8_v2_vec_pool.json')
EPS_F = {1: 0.2, 2: 9.523809523809525e-3, 3: -5.091986325893010e-3, 4: -2.452258248318762e-3,
         5: -1.151746915130986e-3, 6: -4.979056652200001e-4, 7: -1.175236830400000e-3,
         8: -7.455463672900000e-4, 9: -7.520257156400000e-6, 10: 7.207509171100000e-4,
         11: 1.501967012082273e-3, 12: 2.274713720558208e-3}


def p_from_nu(dense, M, m):
    fac = pow(pow(4, -1, M), m, M)
    idx = (np.arange(M) * fac) % M
    return float(np.sum(dense * dense[idx]))


def main():
    t0 = time.time()
    print(f"# PROBE CROSSING -- two-monotone-mode vs long-period.  building nu to r={R_MAX}...\n")
    nus = build_nu(0.5, R_MAX)
    print(f"   nu built ({time.time()-t0:.1f}s)")
    gam = {m: {} for m in range(1, MMAX + 1)}
    for r in range(1, R_MAX + 1):
        M = 3 ** (r + 1)
        dense = np.zeros(M)
        for X, w in nus[r].items():
            dense[X] = w
        for m in range(1, MMAX + 1):
            gam[m][r] = 3 ** r * p_from_nu(dense, M, m)
        del dense
    A = {m: {r: gam[m][r] - gam[m][r - 1] for r in range(2, R_MAX + 1)} for m in range(1, MMAX + 1)}
    Lam = {r: sum(4.0 ** -m * A[m][r] for m in range(1, MMAX + 1)) for r in range(2, R_MAX + 1)}
    print(f"   Lam_r via nu built ({time.time()-t0:.1f}s)\n")

    # validate nu-Lam vs exact/eps-Lam
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in json.load(open(EPSJSON)).items()}
    LamEps = {r: (float(EPS[r + 1] - EPS[r]) / 2 if r + 1 <= 8 else (EPS_F[r + 1] - EPS_F[r]) / 2) for r in range(1, 12)}
    print("## VALIDATION  Lam_r(nu, sum_m 4^-m A) vs Lam_r(eps, telescoping)  r=3..11")
    print(f"   {'r':>2} {'Lam(nu)':>12} {'Lam(eps)':>12} {'rel diff':>10}")
    for r in range(3, 12):
        rd = abs(Lam[r] - LamEps[r]) / abs(LamEps[r]) if abs(LamEps[r]) > 1e-15 else float('nan')
        print(f"   {r:>2} {Lam[r]:>12.4e} {LamEps[r]:>12.4e} {rd:>10.2%}")
    print()

    # THE CROSSING: Lam_r r=8..14, extended PAST the eps data
    print("## THE CROSSING  Lam_r for r=8..14 (r>=12 is NEW, from nu, beyond the eps ladder)")
    print(f"   Wilson pred: crosses once near r=12, NEGATIVE for all r>=13, never returns. 2nd sign change > 12 = KILL.")
    print(f"   {'r':>2} {'Lam_r':>13} {'sign':>5} {'|Lam|/0.984^r':>13} {'source':>8}")
    prevsign = None; flips = []
    for r in range(8, R_MAX + 1):
        s = '+' if Lam[r] > 0 else '-'
        if prevsign and s != prevsign:
            flips.append(r)
        prevsign = s
        src = 'eps+nu' if r <= 11 else 'nu ONLY'
        print(f"   {r:>2} {Lam[r]:>13.4e} {s:>5} {abs(Lam[r])/0.984**r:>13.5f} {src:>8}")
    print(f"   sign flips in r=8..14: {flips}")
    ge12 = [Lam[r] for r in range(12, R_MAX + 1)]
    kill = any(ge12[i] > 0 for i in range(len(ge12))) and any(ge12[i] < 0 for i in range(len(ge12)))
    print(f"   Lam_12..14 = {['%.3e'%x for x in ge12]}")
    print(f"   => {'MODEL HOLDS so far: single crossing ~r=12, stays negative' if (all(x<0 for x in ge12)) else ('KILLED: sign change at r>12 => oscillation, not a difference of exponentials' if kill else 'inconclusive')}\n")

    # one-number check
    T = -EPS_F[12] / 2
    pred = T * (1 - 0.984)
    print(f"## ONE-NUMBER CHECK  Lam_12 vs prediction T*(1-rho) [if single-signed decay from r=12]")
    print(f"   T=Sum_{{r>=12}}Lam = -eps_12/2 = {T:.4e};  T*(1-0.984) = {pred:.3e};  Lam_12(nu) = {Lam[12]:.4e}")
    print(f"   Lam_11 = {Lam[11]:.4e} (+), Lam_12 = {Lam[12]:.4e} => crossing between r=11 and 12? {'YES' if Lam[11]>0>Lam[12] else 'no'}\n")

    # MONO / OSC split
    print("## MONO/OSC SPLIT  MONO=0.25 A(1)+(1/64)A(3);  OSC=0.0625 A(2)+(1/256)A(4)")
    print(f"   {'r':>2} {'MONO':>12} {'OSC':>12} {'|OSC/MONO|':>11} {'sign Lam':>9} {'sign MONO':>10}")
    for r in range(3, R_MAX + 1):
        MONO = 0.25 * A[1][r] + (1 / 64) * A[3][r]
        OSC = 0.0625 * A[2][r] + (1 / 256) * A[4][r]
        rr = abs(OSC) / abs(MONO) if abs(MONO) > 1e-15 else float('nan')
        print(f"   {r:>2} {MONO:>12.4e} {OSC:>12.4e} {rr:>11.3f} {'+' if Lam[r]>0 else '-':>9} {'+' if MONO>0 else '-':>10}")
    print("   [MONO dominates + MONO crosses ~r=12 => two-mono model; OSC>=MONO => KILL (long-period stands).]")


if __name__ == "__main__":
    main()
