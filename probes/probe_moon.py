"""
PROBE MOON -- PHASE-DIVERSITY JOINT FIT. The m-channels A_r(m)=C_{r+1}(m)/3 = gamma_r(tau_m)-gamma_{r-1}(tau_m)
are different observables of the SAME bore operator, sharing (rho,theta) but with per-channel (A_m,phi_m).
Lam_r = sum_m 4^{-m} A_r(m). If the channel phases DIFFER, phase diversity substitutes for the r-window we can't get.

MODEL-ADEQUACY test, NOT a period measurement (R27-A: L(z) not rational => Lam_r not a finite sum of exponentials).
Exact A_r(m) via character ledger r<=7; float via build_nu r=8..12. M-A first; it can kill the probe.
"""
import os, sys, math, cmath, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_engine_R7 as R7
from probe_gamma_R9 import tau, v3
from probe_charledger_R10 import autocorr_dlog
from probe_gapop_R28 import build_nu

MS = [1, 2, 3, 4]
R_EXACT = 7
R_MAX = 12


def A_char(g, r, m):
    N = 3 ** r
    return sum(g[u] * R7.cram(r, (u - m) % N) for u in range(N))


def p_from_nu(dense, M, m):
    fac = pow(pow(4, -1, M), m, M)
    idx = (np.arange(M) * fac) % M
    return float(np.sum(dense * dense[idx]))


def main():
    t0 = time.time()
    print("# PROBE MOON -- PHASE-DIVERSITY JOINT FIT.\n")
    # exact A (r<=7)
    mu = {1: R7.mu1()}
    for k in range(2, R_EXACT + 1):
        mu[k] = R7.build_mu(mu[k - 1], k)
    gcache = {r: autocorr_dlog(mu[r], r) for r in range(1, R_EXACT + 1)}
    Aex = {m: {r: float(A_char(gcache[r], r, m)) for r in range(2, R_EXACT + 1)} for m in MS + [9, 27]}
    print(f"   exact A_r(m) built r<=7  ({time.time()-t0:.1f}s)")
    # float A (all r) via nu
    nus = build_nu(0.5, R_MAX)
    dense = {r: np.zeros(3 ** (r + 1)) for r in range(1, R_MAX + 1)}
    for r in range(1, R_MAX + 1):
        for X, w in nus[r].items():
            dense[r][X] = w
    gam = {m: {r: 3 ** r * p_from_nu(dense[r], 3 ** (r + 1), m) for r in range(1, R_MAX + 1)} for m in MS + [9, 27]}
    A = {m: {r: gam[m][r] - gam[m][r - 1] for r in range(2, R_MAX + 1)} for m in MS + [9, 27]}
    print(f"   float A_r(m) built r<=12  ({time.time()-t0:.1f}s)\n")

    # ================= M-A =================
    print("## M-A  THE DIVERSITY CHECK (kill condition).  A_r(m), r=2..12 [exact r<=7 cross-check in ()]")
    for m in MS:
        cells = []
        for r in range(2, R_MAX + 1):
            ex = f"({Aex[m][r]:+.4f})" if r in Aex[m] else ""
            cells.append(f"r{r}:{A[m][r]:+.5f}{ex}")
        print(f"   m={m}: " + "  ".join(cells))
    print("\n   sign sequences (r=2..12) and turnovers (sign flips):")
    turn = {}
    for m in MS:
        sg = ['+' if A[m][r] > 0 else '-' for r in range(2, R_MAX + 1)]
        flips = [r for r in range(3, R_MAX + 1) if (A[m][r] > 0) != (A[m][r - 1] > 0)]
        turn[m] = flips
        print(f"   m={m}: {' '.join(sg)}   sign-flips at r={flips}")
    print("\n   relative phase A_r(m)/A_r(1) (constant-real => in phase; varying/sign-flipping => diverse):")
    for m in (2, 3, 4):
        cells = [f"r{r}:{A[m][r]/A[1][r]:+.3f}" for r in range(6, R_MAX + 1)]
        print(f"   m={m}/1: " + "  ".join(cells))
    # kill decision: do turnovers coincide within <1 level?
    allflips = [f for m in MS for f in turn[m]]
    spread_ok = False
    if allflips:
        # cluster: are there flips at distinct r separated by >=1 across channels?
        distinct = sorted(set(allflips))
        spread_ok = (max(distinct) - min(distinct)) >= 2
    print(f"\n   turnover r-values across channels: {sorted(set(allflips))}")
    print(f"   => M-A: {'DIVERSITY PRESENT (turnovers spread >=2 levels) -> run M-B' if spread_ok else 'channels ~in-phase (turnovers coincide) -> KILL, redundant channels'}\n")
    if not spread_ok:
        print("   [M-A kill: the m-channels carry one shared phase; M-B/C void. Banked as an eigenvector-structure fact.]")
        return

    # ================= M-D (SNR, gates inclusion) =================
    print("## M-D  CHANNEL SNR (|A_r(m)| at r=10,11,12 vs transient residual ~ 0.0105*2^-r)")
    admit = list(MS)
    for m in MS + [9, 27]:
        trans = [0.0105 * 2 ** -r for r in (10, 11, 12)]
        sig = [abs(A[m][r]) for r in (10, 11, 12)]
        ok = all(s > t for s, t in zip(sig, trans))
        tag = 'ADMIT' if ok else 'reject'
        if m in (9, 27) and ok:
            admit.append(m)
        if m in (9, 27) or not ok:
            print(f"   m={m:>2}: |A| r10,11,12 = {sig[0]:.2e},{sig[1]:.2e},{sig[2]:.2e}  vs trans {trans[0]:.2e}.. [{tag}]")
    print(f"   admitted channels for M-B: {admit}\n")

    # ================= M-B  JOINT FIT =================
    print("## M-B  JOINT FIT  A_r(m)=A_m rho^r cos(r theta+phi_m); shared (rho,theta), r=6..12")
    rs = list(range(6, R_MAX + 1))
    def fit_channels(chans, rho, theta):
        tot = 0.0; params = {}
        for m in chans:
            B = np.array([[rho ** r * math.cos(theta * r), -rho ** r * math.sin(theta * r)] for r in rs])
            y = np.array([A[m][r] for r in rs])
            ab, *_ = np.linalg.lstsq(B, y, rcond=None)
            res = B @ ab - y
            tot += float(np.sum(res ** 2))
            params[m] = ab
        return tot, params
    best = None
    for ip in range(int((30 - 6) / 0.2) + 1):
        P = 6 + 0.2 * ip; th = 2 * math.pi / P
        for rho in np.linspace(0.94, 1.0, 31):
            tot, _ = fit_channels(admit, rho, th)
            if best is None or tot < best[0]:
                best = (tot, rho, th, P)
    _, rho, th, P = best
    _, params = fit_channels(admit, rho, th)
    phis = {m: math.atan2(params[m][1], params[m][0]) for m in admit}
    amps = {m: math.hypot(params[m][0], params[m][1]) for m in admit}
    phlist = sorted(phis[m] for m in admit)
    spread = phlist[-1] - phlist[0]
    print(f"   best joint fit: rho={rho:.4f}  period 2pi/theta = {P:.2f}  (window r=6..12 spans {7*th:.2f} rad = {math.degrees(7*th):.0f}deg)")
    for m in admit:
        print(f"     m={m}: amp={amps[m]:.4e}  phi={math.degrees(phis[m]):+.1f}deg")
    print(f"   phi spread = {spread:.2f} rad ({math.degrees(spread):.0f}deg)   [SUCCESS if >= 2 rad (~115deg)]")
    railed = (rho <= 0.941 or P >= 29.9)
    print(f"   => M-B free fit: {'RAILED to grid boundary (rho floor / P ceiling) -> no interior optimum, under-determined/misspecified' if railed else 'interior optimum found'}\n")
    # rho fixed at the established 0.984: do the phase-diverse channels agree on a PERIOD?
    print("   [rho fixed at 0.984 -- the phase-diversity payoff: do channels agree on P?]")
    def bestP(chans, rho_fix):
        b = None
        for ip in range(int((30 - 6) / 0.1) + 1):
            P2 = 6 + 0.1 * ip; th2 = 2 * math.pi / P2
            tot, _ = fit_channels(chans, rho_fix, th2)
            if b is None or tot < b[0]:
                b = (tot, P2)
        return b[1]
    print(f"     all 6 channels: best P = {bestP(admit,0.984):.1f}   (window spans only {math.degrees(7*2*math.pi/17):.0f}deg at P=17)")

    # ================= M-C  CROSS-VALIDATION =================
    print("\n## M-C  CROSS-VALIDATION (disjoint subsets; do NOT average) -- rho fixed 0.984, compare P")
    p13 = bestP([1, 3], 0.984); p24 = bestP([2, 4], 0.984)
    print(f"   subset {{1,3}} (the MONOTONE channels): best P = {p13:.1f}")
    print(f"   subset {{2,4}} (the OSCILLATORY channels): best P = {p24:.1f}")
    print(f"   separation |P13 - P24| = {abs(p13-p24):.1f}")
    print(f"   [agree => shared-mode adequate; disagree/rail => per-channel continuum contamination (R27-A) -- a finding.]")
    print(f"   also free (rho,theta) subsets:")

    # ================= M-C  CROSS-VALIDATION =================
    print("## M-C  CROSS-VALIDATION (disjoint channel subsets; do NOT average)")
    for label, sub in [("{1,3}", [1, 3]), ("{2,4}", [2, 4])]:
        b = None
        for ip in range(int((30 - 6) / 0.2) + 1):
            P2 = 6 + 0.2 * ip; th2 = 2 * math.pi / P2
            for rho2 in np.linspace(0.94, 1.0, 31):
                tot, _ = fit_channels(sub, rho2, th2)
                if b is None or tot < b[0]:
                    b = (tot, rho2, P2)
        print(f"   subset {label}: rho={b[1]:.4f}  period={b[2]:.2f}")
    print("   [agree => shared-mode model adequate in-range; disagree => per-channel continuum contamination (R27-A), a finding.]\n")

    print("## SCOPE: model-adequacy only. A consistent (rho,theta) = 'a two-mode model fits r=6..12 across channels',")
    print("   NOT 'the period is X' (R27-A: infinitely many modes). Period reported only with window+channels+phase-span.")


if __name__ == "__main__":
    main()
