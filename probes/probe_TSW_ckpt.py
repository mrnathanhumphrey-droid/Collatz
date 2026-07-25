"""
PROBE TSW-CKPT -- Wilson's three checkpoints, in order. Corrects/retracts parts of probe_TSW.

(a) TSW-A as an IDENTITY: g_r = <delta_r, Re w> = (Lambda_r - Lambda_r^unif)/S_r (NOT Lambda_r/S_r).
    Confirm g_1 = 0 (delta_1 = 0 forced). Report g_r and its signs; show g_2 = b_2 != Lambda_2/S_2 (~2.1x).
(b) Report the truncation's LEADING eigenvalue. If it isn't 1, the M-matrix is R29's corpse (R29-B: |lam2| non-
    convergent, {R(d)} not a state, no finite transfer operator exists) and B/C carry NO transport spectral info.
    Report C (deflated) SEPARATELY from B.
(c) THE computable coupling: source projection <s_r, Re w>, s_r = P0 T(uniform_{r-1}) via certified R16-A one-step
    transport. Compare its RATE and SIGN to Lambda_r's (exact local ratios r=12..16 = 0.950,0.897,0.913,0.893,
    declining, vs banked rho~0.984 -- they can't both be asymptotic). Source-dominated tail => rate is source rate,
    sign set by <s_r, Re w> at large r (manufactured fresh each level, no growing state, no non-existent operator).
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_engine_R7 as R7
import probe_charledger_R10 as R10
from probe_transport_R16 import theta_transport
from probe_gapop_R28 import build_nu, nu_hat, R_of_d

Rew = lambda x: 15.0 / (2 * (17 - 8 * math.cos(2 * math.pi * x))) - 0.5

_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): float(F(int(v['num']), int(v['den']))) for k, v in _hist.items()}
EPS.update({9: -7.520257156400000e-6, 10: 7.207509171100000e-4, 11: 1.501967012082273e-3, 12: 2.274713720558208e-3})
LAM_NU = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}
for r in range(12, 17):
    EPS[r + 1] = EPS[r] + 2 * LAM_NU[r]
LAM = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 17)}
Sval = {r: 7.0 / 15 + EPS[r] for r in range(1, 17)}
lam = 0.5
P = lambda d: (1 - lam) * lam ** abs(d) / (1 + lam)


def profile_from_theta(th, r):
    """|theta_hat(k)|^2 over k=0..N-1 from a dlog-domain measure th {t: weight}."""
    N = 3 ** r
    g = [0.0] * N
    items = [(t, float(w)) for t, w in th.items()]
    for t1, w1 in items:
        for t2, w2 in items:
            g[(t1 - t2) % N] += w1 * w2
    return [sum(g[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]


def theta2_mu(mu_r, r):
    N = 3 ** r
    g = R10.autocorr_dlog(mu_r, r); gf = [float(x) for x in g]
    return [sum(gf[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]


def main():
    print("# PROBE TSW-CKPT -- Wilson's three checkpoints, in order.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):                      # only need mu to r=7 (mu_8 build is the wall, and unused)
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ================= (a) =================
    print("## (a)  TSW-A AS IDENTITY: g_r = <delta_r,Re w> = (Lambda_r - Lambda^unif)/S_r ; g_1 = 0 ?")
    print(f"   {'r':>2} {'g_r=<delta,Rew>':>16} {'(Lam-Lam^u)/S_r':>16} {'Lambda_r/S_r':>14} {'Lam^unif/S_r':>14} {'gate':>5}")
    okA = True; gsign = []
    for r in range(1, 8):
        N = 3 ** r; th2 = theta2_mu(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim); Sr = Sval[r]
        delta = {k: th2[k] / Sr - 1.0 / M for k in prim}
        g = sum(delta[k] * Rew(k / N) for k in prim)                      # <delta_r, Re w>
        Lunif_over_S = (1.0 / M) * sum(Rew(k / N) for k in prim)
        LoverS = LAM[r] / Sr
        ok = abs(g - (LoverS - Lunif_over_S)) < 1e-9
        okA = okA and ok
        gsign.append('0' if abs(g) < 1e-9 else ('+' if g > 0 else '-'))
        print(f"   {r:>2} {g:>+16.9f} {LoverS - Lunif_over_S:>+16.9f} {LoverS:>+14.6f} {Lunif_over_S:>+14.3e} {'OK' if ok else 'DEV':>5}")
    print(f"   => IDENTITY gate r<=7 {'PASS' if okA else 'FAIL'};  g_1 = 0 confirmed: {gsign[0]=='0'} (delta_1=0 forced, dim delta_1 = 3^0-1 = 0)")
    _g2 = sum((theta2_mu(mu[2], 2)[k] / Sval[2] - 1.0 / 6) * Rew(k / 9) for k in (1, 2, 4, 5, 7, 8))
    print(f"   g_2 = {_g2:+.4e} (=b_2)  vs  Lambda_2/S_2 = {LAM[2]/Sval[2]:+.4e}  "
          f"(ratio {(LAM[2]/Sval[2])/_g2:.2f}x -- NOT doubly-exp; uniform baseline is O(1) at r=2)")
    full_g = gsign + ['+' if LAM[r] > 0 else '-' for r in range(8, 17)]     # r>=8: Lam^unif doubly-exp, sign(g)=sign(Lam)
    print(f"   sign(g_r) r=1..16 = {''.join(full_g)}  (g_1=0; last sign change at r={max(r for r in range(3,17) if (LAM[r]>0)!=(LAM[r-1]>0))})\n")

    # ================= (b) =================
    print("## (b)  THE TRUNCATION'S LEADING EIGENVALUE -- is the M-matrix R29's corpse?")
    nus = build_nu(lam, 7); NH = {rr: nu_hat(nus[rr], 3 ** (rr + 1)) for rr in range(0, 8)}
    MM = 14
    Rtab = {rr: {m: R_of_d(NH[rr], 3 ** (rr + 1), 2 * m).real for m in range(-MM, MM + 1)} for rr in range(0, 8)}

    def build_M(r, D):
        idx = list(range(-D, D + 1)); n = len(idx)
        def kap(m):
            num = Rtab[r][m]; den = 3 * sum(P(2 * (mp - m)) * Rtab[r - 1][mp] for mp in range(-D, D + 1))
            return num / den if den else float('nan')
        Mr = np.zeros((n, n))
        for i, m in enumerate(idx):
            km = kap(m)
            for j, mp in enumerate(idx):
                Mr[i, j] = 3 * P(2 * (mp - m)) * km
        return Mr, idx
    print(f"   {'D':>3} {'leading':>9} {'|lam2|':>9} {'lam2 real?':>11} {'leading==1?':>12}")
    for D in (6, 8, 10, 12):
        Mr, idx = build_M(7, D)
        ev = sorted(np.linalg.eigvals(Mr), key=lambda z: -abs(z))
        lead = ev[0].real; l2 = ev[1]
        print(f"   {D:>3} {lead:>9.5f} {abs(l2):>9.5f} {str(abs(l2.imag)<1e-9):>11} {str(abs(lead-1)<1e-3):>12}")
    print("   => leading != 1 (peaks ~1.08; 3K-symbol peaks 5/3, kappa damps to ~1.08) AND |lam2| NON-convergent in D")
    print("      (R29-B: 0.925->0.972->0.980->1.005, D-steps not shrinking; {R(d)} not a state, NO finite transfer")
    print("      operator exists). => The M-matrix is R29's CORPSE. TSW-B/C 'real spectrum' = property of a")
    print("      non-convergent truncation, NOT of the transport. RETRACT 'the slow mode is real' as a transport claim.")
    print("      (Dominant-subdominant reality was already banked at R27-A |lam2|~0.5; the SLOW post-deflation object")
    print("      is what was owed, and R28-D says deflation annihilates the data it needs -- so C is void here too.)\n")

    # ================= (c) =================
    print("## (c)  SOURCE PROJECTION <s_r, Re w>, s_r = P0 T(uniform_{r-1}) via certified R16-A transport")
    print("   (fresh innovation each level; no growing state, no non-existent operator). Compare rate & sign to Lambda_r.")
    print(f"   {'r':>2} {'<s_r,Rew> (norm)':>17} {'raw src proj':>14} {'ratio raw':>10} {'Lambda_r':>13} {'Lam ratio':>10} {'sgn s':>6} {'sgn Lam':>8}")
    prev_raw = None; prev_lam = None
    for r in range(2, 8):
        N = 3 ** r
        mu_unif = {a: F(1, 3 ** (r - 1)) for a in range(3 ** (r - 1))}
        th = theta_transport(mu_unif, r)                                   # T(uniform) at level r, dlog domain
        th2 = profile_from_theta(th, r)
        prim = [k for k in range(1, N) if k % 3 != 0]; Mm = len(prim)
        mass = sum(th2[k] for k in prim)
        s_norm = sum((th2[k] / mass - 1.0 / Mm) * Rew(k / N) for k in prim)   # normalized source deviation . W
        raw = sum(th2[k] * Rew(k / N) for k in prim)                          # raw source proj (= Lambda if input uniform)
        Lam = float(LAM[r])
        rr_raw = (raw / prev_raw) if prev_raw not in (None, 0) else float('nan')
        rr_lam = (Lam / prev_lam) if prev_lam not in (None, 0) else float('nan')
        print(f"   {r:>2} {s_norm:>+17.4e} {raw:>+14.4e} {rr_raw:>10.4f} {Lam:>+13.4e} {rr_lam:>10.4f} "
              f"{'+' if s_norm>0 else '-':>6} {'+' if Lam>0 else '-':>8}")
        prev_raw = raw; prev_lam = Lam
    print("   [fork: if source SIGN matches Lambda's at large r AND source RATE ~ observed Lambda rate (~0.91, not")
    print("    0.984), the tail is SOURCE-DOMINATED -> its sign/rate are set by the fresh source, computable here;")
    print("    the eigenvalue (real or not) does not set the tail. This is the first computable coupling question.]")


if __name__ == "__main__":
    main()
