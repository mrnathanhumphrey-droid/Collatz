"""
PROBE AC-LAGS -- the five-lag structure of the dominant mode Re dhat_r(1).

Pure measurement on the r<=16 data (validated build_nu/dlog/profile path). No new operator, no new depth.
Gives the pen the SHAPE of the inequality: which of the five autocorrelation lags does the work.

RATIO-2 closed form:  Re dhat_r(1) = [2 C(1) - (C(N/3-1)+C(N/3+1))] / [2 (C(0) - C(N/3))],  N=3^r,
   C(d)=sum_s rho(s) rho((s+d) mod N).  Wilson's decomposition (algebraically exact):
     G_r := [C(1) - C(N/3)] / [C(0) - C(N/3)]            (decay ratio: near-lag vs fiber-lag)
     K_r := Delta^2 C(N/3) / (2[C(0)-C(N/3)]),  Delta^2 C(N/3):=C(N/3-1)-2C(N/3)+C(N/3+1)   (curvature at fiber lag)
   ==> Re dhat_r(1) = G_r - K_r.   d1~3e-3 with G,K O(1) => near-cancellation; the 0.90 decay is it tightening.

AC-A gate (exact r<=7): 5-lag == banked d1 to 1e-12, and C(0)>C(N/3) strictly.
AC-B normalized table; AC-C cancellation depth; AC-D which term moves; AC-E fiber-lag symmetry; AC-F full shape r=12,16.
"""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_gapop_R28 import build_nu
from probe_ratio2 import rho_dense, d_hat_123, build_nu_exact, rho_exact_dict

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"

def five_lags_float(rho, N):
    n3 = N // 3
    C = lambda dl: float(np.dot(rho, np.roll(rho, -dl)))
    return {0: C(0), 1: C(1), n3 - 1: C(n3 - 1), n3: C(n3), n3 + 1: C(n3 + 1)}

def five_lags_exact(rho, N):
    n3 = N // 3
    def C(dl):
        s = F(0)
        for k, v in rho.items():
            w = rho.get((k + dl) % N)
            if w is not None:
                s += v * w
        return s
    return {0: C(0), 1: C(1), n3 - 1: C(n3 - 1), n3: C(n3), n3 + 1: C(n3 + 1)}

def GKd(C, N):
    n3 = N // 3
    C0, C1, Cm, Cn, Cp = C[0], C[1], C[n3 - 1], C[n3], C[n3 + 1]
    den = C0 - Cn
    G = (C1 - Cn) / den
    K = (Cm - 2 * Cn + Cp) / (2 * den)
    return G, K, G - K


def main():
    t0 = time.time()
    print("# PROBE AC-LAGS -- five-lag structure of Re dhat_r(1), r=2..16\n")
    print(f"building build_nu to r=16 ... (~9 min)")
    nus = build_nu(0.5, 16)
    print(f"  built ({time.time()-t0:.1f}s)\n")

    lags = {}; rho_keep = {}; d1_prof = {}
    for r in range(2, 17):
        rho, N = rho_dense(nus[r], r)
        lags[r] = five_lags_float(rho, N)
        d1_prof[r] = d_hat_123(rho, N)[0]                 # independent profile-path reference, all r
        if r in (12, 13, 14, 15, 16):
            np.save(os.path.join(SCRATCH, f"rho_r{r}.npy"), rho)
        if r in (12, 16):
            rho_keep[r] = rho                              # keep for AC-F
        else:
            del rho
    del nus
    print(f"  five lags + profile-path d1 computed, rho_12..16 dumped to scratchpad ({time.time()-t0:.1f}s)\n")

    # ---------- AC-A GATE ----------
    print("## AC-A  GATE: 5-lag closed form vs profile-path d1, r=2..16; exact r<=7; C(0)>C(N/3)")
    nex = build_nu_exact(7)
    worst = 0.0; posfail = []
    print(f"   {'r':>2} {'d1(5-lag)':>18} {'ref':>18} {'rel':>10}")
    d1 = {}
    for r in range(2, 17):
        N = 3 ** r
        _, _, d = GKd(lags[r], N); d1[r] = d
        if r <= 7:                                         # exact reference
            rho_e, Ne = rho_exact_dict(nex[r], r)
            _, _, de = GKd(five_lags_exact(rho_e, Ne), Ne)
            rel = abs(d - float(de)) / abs(float(de)); ref = f"{float(de):+.12e}"
        else:                                              # float profile-path reference
            rel = abs(d - d1_prof[r]) / abs(d1_prof[r]); ref = f"{d1_prof[r]:+.12e}"
        worst = max(worst, rel)
        if lags[r][0] <= lags[r][N // 3]:
            posfail.append(r)
        print(f"   {r:>2} {d:+.12e} {ref:>18} {rel:>10.1e}")
    print(f"\n   worst rel = {worst:.2e}  [{'GATE PASS' if worst < 1e-12 else 'GATE FAIL'}]")
    print(f"   C(0) > C(N/3) strictly at every r: {'YES' if not posfail else 'NO at '+str(posfail)}\n")
    if worst >= 1e-12 or posfail:
        print("   *** GATE issue -- stopping. ***"); return

    # ---------- AC-B  normalized table ----------
    print("## AC-B  normalized lags and the G-K decomposition (auditable: G-K should = d1)")
    print(f"   {'r':>2} {'C1/C0':>10} {'Cn3/C0':>10} {'Cm/C0':>10} {'Cp/C0':>10} {'G':>12} {'K':>12} {'G-K':>13} {'d1':>13}")
    G = {}; K = {}
    for r in range(2, 17):
        N = 3 ** r; L = lags[r]; C0 = L[0]; n3 = N // 3
        g, k, d = GKd(L, N); G[r] = g; K[r] = k
        print(f"   {r:>2} {L[1]/C0:>10.6f} {L[n3]/C0:>10.6f} {L[n3-1]/C0:>10.6f} {L[n3+1]/C0:>10.6f} "
              f"{g:>12.8f} {k:>12.8f} {d:>+13.6e} {d1[r]:>+13.6e}")
    print()

    # ---------- AC-C  cancellation depth ----------
    print("## AC-C  cancellation depth of G-K, and float floor on the difference")
    print(f"   {'r':>2} {'G':>12} {'K':>12} {'|G-K|':>12} {'digits coincide':>16} {'float floor rel(G-K)':>22}")
    for r in range(2, 17):
        g, k = G[r], K[r]; dd = abs(g - k)
        digits = -math.log10(dd / max(abs(g), abs(k))) if dd > 0 else float('inf')
        floor = 1e-16 * max(abs(g), abs(k)) / dd
        print(f"   {r:>2} {g:>12.8f} {k:>12.8f} {dd:>12.4e} {digits:>16.2f} {floor:>22.2e}")
    print("   [digits coincide = -log10(|G-K|/max(|G|,|K|)); float floor = eps*max(|G|,|K|)/|G-K| relative to d1.]\n")

    # ---------- AC-D  which term moves ----------
    print("## AC-D  are G,K converging? successive ratios/diffs alongside RATIO-2 rho_r=d1_{r+1}/d1_r")
    print(f"   {'r':>2} {'G':>11} {'dG':>11} {'K':>11} {'dK':>11} {'rho_r(d1)':>11}")
    for r in range(2, 17):
        dG = G[r] - G[r - 1] if r - 1 in G else float('nan')
        dK = K[r] - K[r - 1] if r - 1 in K else float('nan')
        rr = d1[r + 1] / d1[r] if r + 1 in d1 else float('nan')
        print(f"   {r:>2} {G[r]:>11.7f} {dG:>+11.2e} {K[r]:>11.7f} {dK:>+11.2e} {rr:>11.5f}")
    print("   readings: both->common L (diff decays 0.90)=approach question | one drifts=localizes | neither=wrong split\n")

    # ---------- AC-E  fiber-lag symmetry ----------
    print("## AC-E  fiber-lag antisymmetry  A_r = (C(N/3+1)-C(N/3-1))/C(0)")
    print(f"   {'r':>2} {'A_r':>14} {'|A_r|/(C1/C0)':>16}")
    for r in range(2, 17):
        N = 3 ** r; L = lags[r]; n3 = N // 3
        A = (L[n3 + 1] - L[n3 - 1]) / L[0]
        print(f"   {r:>2} {A:>+14.6e} {abs(A)/(L[1]/L[0]):>16.4e}")
    print("   [A_r==0 => C locally symmetric at N/3 => inequality shortens to C(1)>C(N/3+1). else size bounds the curvature simplification.]\n")

    # ---------- AC-F  full autocorrelation shape at r=12,16 ----------
    print("## AC-F  full C(k)/C(0) shape (r=12,16): fine windows +-5 at k=0,N/3,2N/3, plus coarse overview")
    for r in (12, 16):
        N = 3 ** r; n3 = N // 3; rho = rho_keep[r]; C0 = lags[r][0]
        def Cn(k):
            return float(np.dot(rho, np.roll(rho, -(k % N)))) / C0
        print(f"\n   --- r={r}, N={N} ---")
        for name, ctr in (("k=0", 0), ("k=N/3", n3), ("k=2N/3", 2 * n3)):
            vals = [(d, Cn(ctr + d)) for d in range(-5, 6)]
            print(f"   fine @ {name}: " + " ".join(f"{d:+d}:{v:.5f}" for d, v in vals))
        # coarse overview: 20 sampled points across [0,N)
        ks = [int(round(i * (N - 1) / 19)) for i in range(20)]
        cs = [Cn(k) for k in ks]
        print("   coarse C(k)/C0 (20 pts 0..N-1): " + " ".join(f"{v:+.4f}" for v in cs))
        # save full 200-grid for audit
        grid = [int(round(i * (N - 1) / 199)) for i in range(200)]
        gv = np.array([Cn(k) for k in grid])
        np.save(os.path.join(SCRATCH, f"ac_shape_r{r}.npy"), np.array([grid, gv]))
        print(f"   [full 200-pt grid saved -> scratchpad/ac_shape_r{r}.npy]  ({time.time()-t0:.1f}s)")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
