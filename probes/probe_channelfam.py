"""
PROBE CHANNELFAM -- the channel-constant family gamma_inf(k) and the self-check Sum_k 4^-k gamma_inf(k) = S_inf/2.
(Wilson's probe order, 2026-07-26.)

gamma_r(k) = 3^r p_r(k),  p_r(k) = <rho_r, shift_k rho_r> (lag-k autocorr of the dlog profile); A_r(k)=gamma_r(k)-gamma_{r-1}(k)
(gamma_0(k)=1). Re dhat_r(n)=A_r(n)/A_r(0) for 3-nmid n (CHANNEL_ID). Lambda_r = Sum_{k>=1} 4^-k A_r(k) (banked =(eps_{r+1}-eps_r)/2).

STEP 1 -- NORMALIZATION GATE (derived against the record, sidesteps the Cauchy-transform factor Wilson flagged):
  two telescopes of Lambda_r  =>  Sum_{k>=1} 4^-k gamma_inf(k) = 1/3 + (S_inf - S_1)/2.
  So the boxed relation  Sum_k 4^-k gamma_inf(k) = S_inf/2  holds  IFF  S_1 = 2/3  (i.e. p_1(0)=5/9, gamma_1(0)=5/3).
  Gates: (a) S_1 == 2/3 EXACT; (b) bookkeeping Sum_{r<=R} Lambda_r == Sum_k 4^-k (gamma_R(k)-1) to machine precision;
         (c) Sum_{k>=1} 4^-k gamma_R(k)  vs  S_R/2  at R=16 (the finite-R headline).

STEP 2 -- full-precision gamma_r(k), r<=6 EXACT rationals, k=1..6; convention check Re dhat_r(k)=A_r(k)/A_r(0) vs banked d_k.
STEP 3 -- per-channel deparitied rate rho_k; gamma_inf(k) as an INTERVAL from the rate band (k=1..4 carry 99.7% of 4^-k weight).
STEP 4 -- HEADLINE: 2 * Sum_{k>=1} 4^-k gamma_inf(k)  vs 7/15=0.46667 and 0.473 (the tail-free determination of S_inf).
STEP 5 (PSLQ) deferred until step 4 lands.

Reuses probe_channel_audit rho-build: exact r<=5 (build_nu_exact), float r=1..11 (build_nu), cached rho_12..16.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"
KLAG = 24           # lags to carry for the Lambda sum + rates
NORM_S1 = F(2, 3)


def rho_exact_norm(nu, r):
    N = 3 ** r; d = R10.dlog_table(r)
    rho = {}
    for X, w in nu.items():
        s = d[(X - 1) // 3 % N]
        rho[s] = rho.get(s, F(0)) + w
    tot = sum(rho.values())
    return {s: w / tot for s, w in rho.items()}, N


def C_ex(rho, N, k):
    return sum(w * rho.get((s + k) % N, F(0)) for s, w in rho.items())


def main():
    t0 = time.time()
    print("# PROBE CHANNELFAM -- gamma_inf(k) family + self-check Sum 4^-k gamma_inf(k) = S_inf/2\n")

    # ---- exact rho r<=5 ----
    nex = build_nu_exact(5)
    rex = {r: rho_exact_norm(nex[r], r)[0] for r in range(1, 6)}
    # ---- float rho r=1..11 ----
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

    # p_r(k), gamma_r(k), A_r(k)
    p = {0: {k: 1.0 for k in range(KLAG + 1)}}
    for r in range(1, 17):
        p[r] = {k: float(np.dot(rho[r], np.roll(rho[r], -k))) for k in range(KLAG + 1)}
    gam = {r: {k: 3.0 ** r * p[r][k] for k in range(KLAG + 1)} for r in range(0, 17)}
    A = {r: {k: gam[r][k] - gam[r - 1][k] for k in range(KLAG + 1)} for r in range(1, 17)}
    Lam = {r: sum(4.0 ** -k * A[r][k] for k in range(1, KLAG + 1)) for r in range(1, 17)}

    # ================= STEP 1: normalization gate =================
    print("## STEP 1 -- NORMALIZATION GATE  (relation holds IFF S_1 = 2/3)")
    # (a) exact S_1
    N1 = 3
    p1_0 = C_ex(rex[1], N1, 0)
    g1_0 = 3 * p1_0
    S1 = g1_0 - 1          # A_1(0) = gamma_1(0) - gamma_0(0)
    print(f"   (a) EXACT: p_1(0)=||rho_1||^2 = {p1_0} ; gamma_1(0)=3p_1(0) = {g1_0} ; "
          f"S_1=A_1(0) = {S1}   ==2/3: {S1 == NORM_S1}  (=> relation = S_inf/2 EXACTLY)")
    # (b) bookkeeping telescope
    worst = 0.0
    for R in (8, 12, 16):
        lhs = sum(Lam[r] for r in range(1, R + 1))
        rhs = sum(4.0 ** -k * (gam[R][k] - 1.0) for k in range(1, KLAG + 1))
        worst = max(worst, abs(lhs - rhs))
    print(f"   (b) bookkeeping  Sum_{{r<=R}}Lambda_r == Sum_k 4^-k(gamma_R(k)-1): worst |diff| = {worst:.2e}  "
          f"[{'PASS' if worst < 1e-9 else 'FAIL'}]")
    # (c) finite-R headline: Sum 4^-k gamma_R(k) vs S_R/2
    print("   (c) Sum_{k>=1} 4^-k gamma_R(k)  vs  S_R/2   (S_R = A_R(0), banked ladder):")
    for R in (12, 14, 16):
        lhs = sum(4.0 ** -k * gam[R][k] for k in range(1, KLAG + 1))
        SR = A[R][0]
        print(f"       R={R}: Sum = {lhs:.6f} ; S_R/2 = {SR/2:.6f} ; 2*Sum = {2*lhs:.6f} vs S_R = {SR:.6f}  "
              f"rel {abs(2*lhs-SR)/SR:.1e}")
    print()

    # ================= STEP 2: exact gamma_r(k) + convention =================
    print("## STEP 2 -- EXACT gamma_r(k), r<=5, k=1..6  (rationals); convention Re dhat_r(k)=A_r(k)/A_r(0)")
    gex = {}
    for r in range(1, 6):
        N = 3 ** r
        gex[r] = {k: F(3) ** r * C_ex(rex[r], N, k) for k in range(0, 7)}
    for r in range(2, 6):
        Aex = {k: gex[r][k] - gex[r - 1][k] for k in range(0, 7)}
        d1 = Aex[1] / Aex[0]; d2 = Aex[2] / Aex[0]
        print(f"   r={r}: gamma_r(1..3) = {float(gex[r][1]):.5f},{float(gex[r][2]):.5f},{float(gex[r][3]):.5f} | "
              f"Re dhat(1)=A(1)/A(0) = {float(d1):+.6f}  Re dhat(2) = {float(d2):+.6f}")
    print(f"   [d1_2 should be banked 2/35={2/35:+.6f}]  gamma_1(0..3)={[str(gex[1][k]) for k in range(4)]}\n")

    # ================= STEP 3: per-channel rates + gamma_inf intervals =================
    print("## STEP 3 -- per-channel deparitied rate rho_k and gamma_inf(k) INTERVAL  (k=1..6)")
    print(f"   {'k':>2} {'3|k?':>4} {'gamma_16(k)':>12} {'A_16(k)':>11} {'rho_k(deparity)':>15} {'gamma_inf(k) band':>22}")
    ginf = {}
    for k in range(1, 7):
        # deparitied rate from |A_r(k)| pairs: s_r=(A_r+A_{r+1})/2, ratio at top; fallback two-step
        Aseq = {r: A[r][k] for r in range(4, 17)}
        # two-step magnitude rate around top
        def tstep(rr):
            return (abs(Aseq[rr] / Aseq[rr - 2])) ** 0.5 if abs(Aseq[rr - 2]) > 1e-18 else float('nan')
        rates = [tstep(rr) for rr in (14, 15, 16) if not np.isnan(tstep(rr))]
        rk = float(np.median(rates)) if rates else float('nan')
        band = [max(0.0, rk - 0.06), min(0.999, rk + 0.06)]
        lo = gam[16][k] + A[16][k] * band[0] / (1 - band[0])
        hi = gam[16][k] + A[16][k] * band[1] / (1 - band[1])
        g_lo, g_hi = min(lo, hi), max(lo, hi)
        ginf[k] = (g_lo, g_hi, 0.5 * (g_lo + g_hi))
        print(f"   {k:>2} {'yes' if k%3==0 else '':>4} {gam[16][k]:>12.6f} {A[16][k]:>+11.3e} {rk:>15.4f} "
              f"[{g_lo:.5f}, {g_hi:.5f}]")
    print()

    # ================= STEP 4: headline =================
    print("## STEP 4 -- HEADLINE: S_inf = 2 * Sum_{k>=1} 4^-k gamma_inf(k)  (tail-free)")
    wts = {k: 4.0 ** -k for k in range(1, 7)}
    wsum = sum(wts.values())
    print(f"   weights 4^-k (k=1..6): {[f'{wts[k]:.4f}' for k in range(1,7)]} ; cumulative {wsum:.4f} of 1/3 total "
          f"(= {wsum/(1/3)*100:.1f}%)")
    # use midpoint gamma_inf(k) for k=1..6; higher k negligible (weight < 4^-7 ~ 6e-5) -- approximate tail by gamma_16
    S_lo = 2 * (sum(4.0 ** -k * ginf[k][0] for k in range(1, 7))
                + sum(4.0 ** -k * gam[16][k] for k in range(7, KLAG + 1)))
    S_hi = 2 * (sum(4.0 ** -k * ginf[k][1] for k in range(1, 7))
                + sum(4.0 ** -k * gam[16][k] for k in range(7, KLAG + 1)))
    S_mid = 2 * (sum(4.0 ** -k * ginf[k][2] for k in range(1, 7))
                 + sum(4.0 ** -k * gam[16][k] for k in range(7, KLAG + 1)))
    print(f"   S_inf estimate = 2*Sum 4^-k gamma_inf(k) = {S_mid:.5f}   band [{min(S_lo,S_hi):.5f}, {max(S_lo,S_hi):.5f}]")
    print(f"   compare:  7/15 = {7/15:.5f}   |   0.473 (deparity/S2)   |   11/15 for gamma(1)? = {11/15:.5f}")
    sep = "separates" if (min(S_lo, S_hi) > 7/15 or max(S_lo, S_hi) < 7/15) else "does NOT separate 7/15"
    print(f"   => interval {sep} from 7/15.")
    # cross-check: does the finite-R Sum*2 already lean?
    print(f"   [finite-R anchor: 2*Sum 4^-k gamma_16(k) (no extrap) = "
          f"{2*sum(4.0**-k*gam[16][k] for k in range(1,KLAG+1)):.5f}]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
