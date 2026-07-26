"""
PROBE CHANNELDICH -- the enriched/depleted dichotomy across channels + gamma_inf(1) precision (Wilson, 2026-07-26).

gamma_r(k)=3^r C_r(k), C_r(k)=<rho_r, shift_k rho_r>.  WHITE BASELINE = 1 exactly (uniform rho: C=1/3^r => gamma=1).
So gamma_inf(k) > 1 = ENRICHED (positive autocorr at lag k vs white); < 1 = DEPLETED (anti-correlated).

WILSON'S DICHOTOMY (6/6 at k<=6): 3|k channels ENRICHED (fiber-periodic replica peaks), 3-nmid DEPLETED.
  PREDICTION (pre-registered): k=7,8 < 1 (depleted), k=9 > 1 (enriched), k=12 > 1.

STEP A  dichotomy table k=1..12: gamma_16(k), enriched/depleted, 3|k?, per-channel deparitied rate, gamma_inf(k).
STEP B  gamma_16(1) to 12 digits; gamma_inf(1) tight interval (deparitied rate + physical ~0.80 band); TEST 11/15=0.73333.
STEP C  relation check: 3*Sum 4^-k gamma_inf(k) = 3 S_inf/2 ; gamma_inf(2) predicted ~0.474 (weighted-mean residual) vs measured.
STEP D  denominator<=45 rational scan on gamma_inf(k), honestly scoped (resolves 1/45=0.022 spacing, nothing finer).

Reuses probe_channelfam rho-build (exact r<=5, float 1..11, cached 12..16). Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact
from probe_channelfam import rho_exact_norm, C_ex

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"
KMAX = 12


def nearest_rational(x, qmax=45, tol=0.005):
    best = None
    for q in range(1, qmax + 1):
        p = round(x * q)
        from math import gcd
        if gcd(p, q) != 1:
            continue
        v = p / q
        if abs(v - x) < tol and (best is None or abs(v - x) < best[2]):
            best = (p, q, abs(v - x))
    return best


def main():
    t0 = time.time()
    print("# PROBE CHANNELDICH -- enriched/depleted dichotomy (white baseline=1) + gamma_inf(1) precision\n")

    nex = build_nu_exact(5)
    rex = {r: rho_exact_norm(nex[r], r)[0] for r in range(1, 6)}
    nus = build_nu(0.5, 11)
    rho = {}
    for r in range(1, 12):
        N = 3 ** r
        mu = np.zeros(N)
        for X, w in nus[r].items():
            mu[(X - 1) // 3 % N] += float(w)
        d = R10.dlog_table(r)
        rr = np.zeros(N); rr[np.fromiter((d[a] for a in range(N)), np.int64, N)] = mu
        rho[r] = rr / rr.sum()
    del nus
    for r in range(12, 17):
        rr = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho[r] = rr / rr.sum()
    print(f"  (rho ready, {time.time()-t0:.1f}s)\n")

    gam = {r: {k: 3.0 ** r * float(np.dot(rho[r], np.roll(rho[r], -k))) for k in range(KMAX + 1)}
           for r in range(1, 17)}
    gam[0] = {k: 1.0 for k in range(KMAX + 1)}
    A = {r: {k: gam[r][k] - gam[r - 1][k] for k in range(KMAX + 1)} for r in range(1, 17)}

    def deparity_rate(k):
        Aseq = {r: A[r][k] for r in range(4, 17)}
        rr = [(abs(Aseq[x] / Aseq[x - 2])) ** 0.5 for x in (14, 15, 16) if abs(Aseq[x - 2]) > 1e-18]
        return float(np.median(rr)) if rr else float('nan')

    def gamma_inf(k, rate):
        if not (0 < rate < 1):
            return gam[16][k]        # oscillating/unstable -> use r=16 value
        return gam[16][k] + A[16][k] * rate / (1 - rate)

    # ---------- STEP A: dichotomy ----------
    print("## STEP A -- DICHOTOMY: gamma_inf(k) vs white baseline 1  (3|k ENRICHED, 3-nmid DEPLETED?)")
    print(f"   {'k':>2} {'3|k':>3} {'gamma_16(k)':>12} {'rate':>6} {'gamma_inf(k)':>12} {'vs 1':>10} {'predict':>9}")
    ginf = {}
    ok = True
    for k in range(1, KMAX + 1):
        rate = deparity_rate(k)
        gi = gamma_inf(k, rate)
        ginf[k] = gi
        enr = "ENRICHED" if gi > 1 else "depleted"
        three = (k % 3 == 0)
        pred = ("enriched" if three else "depleted")
        match = (gi > 1) == three
        ok = ok and match
        print(f"   {k:>2} {'yes' if three else '':>3} {gam[16][k]:>12.6f} {rate:>6.3f} {gi:>12.6f} "
              f"{enr:>10} {pred:>9} {'' if match else '  <-- MISS'}")
    print(f"   => dichotomy (3|k enriched / 3-nmid depleted) holds k=1..{KMAX}: {ok}")
    print(f"      [pre-registered: k=7,8 depleted, k=9 enriched, k=12 enriched]\n")

    # ---------- STEP B: gamma_inf(1) precision + 11/15 ----------
    print("## STEP B -- gamma_16(1) precision + gamma_inf(1) interval; TEST 11/15 = 0.7333333")
    print(f"   gamma_16(1) = {gam[16][1]:.12f}   (float rho_16; ~10-12 digit)")
    print(f"   A_16(1) = {A[16][1]:.9e} ; deparitied rate = {deparity_rate(1):.4f}")
    r_dep = deparity_rate(1)
    lo_rate, hi_rate = 0.78, min(0.95, r_dep + 0.03)     # physical band: deparity says falling toward ~0.80
    g_lo = gamma_inf(1, lo_rate); g_hi = gamma_inf(1, hi_rate)
    print(f"   gamma_inf(1) over rate [{lo_rate},{hi_rate}]: [{min(g_lo,g_hi):.6f}, {max(g_lo,g_hi):.6f}]")
    print(f"   11/15 = {11/15:.6f}  =>  {'INSIDE interval (survives)' if min(g_lo,g_hi) <= 11/15 <= max(g_lo,g_hi) else 'OUTSIDE (retired)'}")
    print(f"   (11/15 needs rate rho_1 = {(11/15 - gam[16][1])/((11/15 - gam[16][1]) + A[16][1]):.4f} "
          f"solving gamma_16 + A*rho/(1-rho) = 11/15)\n")

    # ---------- STEP C: relation + gamma_inf(2) prediction ----------
    print("## STEP C -- relation 3*Sum 4^-k gamma_inf(k) = 3 S_inf/2 ; gamma_inf(2) prediction")
    # use gamma_inf(1)=11/15 hypothesis to predict S_inf, then back out gamma_inf(2)
    S_from = lambda g1, g2: 2 * (4**-1 * g1 + 4**-2 * g2 + sum(4.0**-k * ginf[k] for k in range(3, KMAX + 1)))
    # weighted mean at S_inf=0.4737 (g1=11/15): known-channels residual gives g2
    S_hyp = 0.4737
    known = 4**-1 * (11/15) + sum(4.0**-k * ginf[k] for k in range(3, KMAX + 1))
    g2_pred = (S_hyp / 2 - known) / 4**-2
    print(f"   IF gamma_inf(1)=11/15 & S_inf=0.4737: predicted gamma_inf(2) = {g2_pred:.5f}")
    print(f"   measured gamma_16(2) = {gam[16][2]:.6f} (oscillating, rate {deparity_rate(2):.3f} -> ~this value)")
    print(f"   consistent: {abs(g2_pred - gam[16][2]) < 0.01}  (dichotomy: {g2_pred:.4f} < 1, 3-nmid 2 depleted {'OK' if g2_pred<1 else 'NO'})")
    S_mid = S_from(11/15, gam[16][2])
    print(f"   => S_inf (g1=11/15, g2=meas) = {S_mid:.5f}\n")

    # ---------- STEP D: denominator<=45 scan ----------
    print("## STEP D -- denominator<=45 scan on gamma_inf(k)  (resolves 1/45=0.022 spacing; nothing finer)")
    for k in range(1, 10):
        nr = nearest_rational(ginf[k], 45, 0.004)
        tag = f"{nr[0]}/{nr[1]} = {nr[0]/nr[1]:.5f} (|d|={nr[2]:.4f})" if nr else "no q<=45 within 0.004"
        star = " <=15|45?" if nr and nr[1] in (15, 45, 30, 45) else ""
        print(f"   k={k}: gamma_inf = {ginf[k]:.6f}  ->  {tag}{star}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
