"""
PROBE TSW -- TRANSPORT-ADJOINT DIRECTION CONVERGENCE.

Question (the whole ballgame for 7/15): is the SLOW eigenvalue of the transport-adjoint T*, restricted to the
mean-zero primitive subspace, REAL or a COMPLEX PAIR?
  REAL  -> g_r = <delta_r, Re w> eventually single-signed -> Sum Lambda_r converges -> 7/15 holds
           (the S2 r~36 rollover is the FINAL sign change).
  CPLX  -> g_r oscillates with lengthening period -> convergence not guaranteed -> S_inf ~ 0.477 live.

Certified operators only (guardrail): freq-domain identity from R14/R11 (Re w, profile, g_r); the linear transport
operator is the R28/R29 GAP MATRIX M (the one certified operator that advances the renewal r-1 -> r linearly on a
FIXED lattice {R(2m)}, leading -> rho, subdominant = lam2). No fresh construction.

TSW-A  freq-domain ground-truth gate: g_r = <delta_r, Re w> = Lambda_r/S_r (up to Lambda^unif, doubly-exp small);
       print the sign sequence the eigenvector story must explain.
TSW-B  direction convergence of (T*)^r W on M's mean-zero complement: overlap -> 1 (real) or oscillates (complex).
TSW-C  deflate the transient(s) first, re-run -> the slow mode's true character.
TSW-D  if REAL: slow eigenvector phi, sign of coupling, asymptotic sign of g_r, last-sign-change vs S2's r~36.
TSW-E  if COMPLEX: arg(lam_slow), phase increment per level (shrinking = lengthening period = beyond numerical reach).
"""
import os, sys, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_engine_R7 as R7
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu, nu_hat, R_of_d

# ---- Re w = W(x) (certified, R14/R11: Re[1/(4 e(x)-1)] = (4cos-1)/(17-8cos) = 15/(2(17-8cos)) - 1/2) ----
Rew = lambda x: 15.0 / (2 * (17 - 8 * math.cos(2 * math.pi * x))) - 0.5

# ---- banked exact/validated eps, Lambda, S (r=1..16) ----
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): float(F(int(v['num']), int(v['den']))) for k, v in _hist.items()}   # exact r<=8
EPS.update({9: -7.520257156400000e-6, 10: 7.207509171100000e-4, 11: 1.501967012082273e-3, 12: 2.274713720558208e-3})
LAM_NU = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}   # nu-validated (CROSSING)
for r in range(12, 17):
    EPS[r + 1] = EPS[r] + 2 * LAM_NU[r]
LAM = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 17)}
Sval = {r: 7.0 / 15 + EPS[r] for r in range(1, 17)}

lam = 0.5
P = lambda d: (1 - lam) * lam ** abs(d) / (1 + lam)


# ================= TSW-A =================
def theta2_profile(mu_r, r):
    """|theta_hat(k)|^2 on Z/3^r (real), from exact autocorrelation (R14 machinery)."""
    N = 3 ** r
    g = R10.autocorr_dlog(mu_r, r); gf = [float(x) for x in g]
    return [sum(gf[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]


def tsw_A():
    print("## TSW-A  GROUND-TRUTH GATE: g_r = <delta_r, Re w> = Lambda_r/S_r ; sign sequence")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)
    print(f"   {'r':>2} {'<delta_r,Re w>':>16} {'Lambda_r/S_r':>16} {'Lambda^unif/S_r':>16} {'gate':>6}")
    okA = True
    for r in range(1, 8):
        N = 3 ** r; th2 = theta2_profile(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim); Sr = Sval[r]
        delta = {k: th2[k] / Sr - 1.0 / M for k in prim}
        inner = sum(delta[k] * Rew(k / N) for k in prim)            # <delta_r, Re w>
        Lunif_over_S = (1.0 / M) * sum(Rew(k / N) for k in prim)     # = Lambda^unif / S_r
        LoverS = LAM[r] / Sr
        ok = abs(inner - (LoverS - Lunif_over_S)) < 1e-9
        okA = okA and ok
        print(f"   {r:>2} {inner:>+16.9f} {LoverS:>+16.9f} {Lunif_over_S:>+16.3e} {'OK' if ok else 'DEV':>6}")
    signs = ''.join('+' if LAM[r] > 0 else '-' for r in range(1, 17))
    print(f"   => identity gate r<=7 {'PASS' if okA else 'FAIL -- operator mis-built, STOP'} "
          f"(g_r = Lambda_r/S_r up to Lambda^unif, doubly-exp small)")
    print(f"   sign(g_r) r=1..16 = {signs}")
    print(f"   pre-registered:      --+++-+++++++++   [match: {signs == '--+++-+++++++++'}]")
    print(f"   LAST sign change at r={max(r for r in range(2,17) if (LAM[r]>0)!=(LAM[r-1]>0))} "
          f"(then single-signed {'+' if LAM[16]>0 else '-'} through r=16)\n")
    return okA


# ================= gap matrix M (certified R29) =================
def build_M(r, D, Rtab):
    idx = list(range(-D, D + 1)); n = len(idx)

    def kappa(m):
        num = Rtab[r][m]
        den = 3 * sum(P(2 * (mp - m)) * Rtab[r - 1][mp] for mp in range(-D, D + 1))
        return num / den if den else float('nan')
    Mr = np.zeros((n, n))
    for i, m in enumerate(idx):
        km = kappa(m)
        for j, mp in enumerate(idx):
            Mr[i, j] = 3 * P(2 * (mp - m)) * km
    return Mr, idx


def spectrum_report(Mr):
    ev, evec = np.linalg.eig(Mr)
    order = np.argsort(-np.abs(ev))
    ev = ev[order]; evec = evec[:, order]
    return ev, evec


def tsw_B_C_D_E(Rtab):
    print("## TSW-B/C  DIRECTION CONVERGENCE of (T*)^r W on M's mean-zero complement (certified R29 gap operator)")
    print("   T* = M^T; W_gap[m] = 4^-|m| (the certified Lambda-channel weights, = gap-domain Re w); P_0 removes leading.")
    r = 7
    for D in (8, 10, 12):
        Mr, idx = build_M(r, D, Rtab)
        n = len(idx)
        ev, evec = spectrum_report(Mr)
        # top spectrum real/complex
        top = ", ".join(
            f"{ev[i].real:+.4f}" + (f"{ev[i].imag:+.4f}i(|.|={abs(ev[i]):.4f},arg={math.degrees(cmath.phase(ev[i])):.1f}d)"
                                    if abs(ev[i].imag) > 1e-9 else "") for i in range(min(6, n)))
        ncomplex = int(np.sum(np.abs(ev.imag) > 1e-9))
        # leading right/left eigenvectors for Hotelling deflation
        uR = evec[:, 0].real; uR = uR / np.linalg.norm(uR)
        evL, evecL = np.linalg.eig(Mr.T)
        oL = np.argsort(-np.abs(evL)); uL = evecL[:, oL[0]].real; uL = uL / (uL @ uR)
        Mdef = Mr - ev[0].real * np.outer(uR, uL)                 # leading removed
        # seed: W_gap, projected off leading (mean-zero w.r.t. the stationary mode)
        W = np.array([4.0 ** -abs(m) for m in idx]); W = W - uR * (uL @ W)
        # TSW-B: power-iterate M^T (=Mdef^T after deflation) ; direction overlap + Rayleigh
        v = W / np.linalg.norm(v0 := W)
        print(f"   --- D={D} (r={r}, {n} states) --- top6 eig: {top}   [#complex in spectrum: {ncomplex}]")
        print(f"       {'it':>3} {'overlap<v_it,v_it-1>':>20} {'Rayleigh v.M^T v':>16} {'|proj on leading|':>17}")
        prev = None
        for it in range(1, 15):
            w = Mdef.T @ v
            nv = np.linalg.norm(w)
            vn = w / nv
            ov = float(vn @ v) if prev is not None else float('nan')
            rq = float(vn @ (Mdef.T @ vn))
            leadproj = abs(float(uR @ vn))                        # should stay ~0 (mean-zero guardrail)
            if it <= 8 or it == 14:
                print(f"       {it:>3} {ov:>20.6f} {rq:>16.6f} {leadproj:>17.2e}")
            prev = vn; v = vn
        # converged direction character
        rq_final = float(v @ (Mdef.T @ v))
        evd, evecd = np.linalg.eig(Mdef); od = np.argsort(-np.abs(evd))
        slow = evd[od[0]]
        slow_complex = abs(slow.imag) > 1e-9
        print(f"       => slow (deflated-dominant) eigenvalue = {slow.real:+.5f}"
              + (f"{slow.imag:+.5f}i  |.|={abs(slow):.5f} arg={math.degrees(cmath.phase(slow)):.2f}d  COMPLEX PAIR"
                 if slow_complex else f"  REAL")
              + f";  Rayleigh(converged dir)={rq_final:+.5f}  [#complex in FULL spectrum: {ncomplex}]")
        # TSW-D / TSW-E
        if not slow_complex:
            # correct coupling: g_r ~ lam^r <delta_0, phi_L>; phi_L = left slow eigvec (= v, the M^T power-iterate limit).
            # sign of g_r tail = sign of the observed data's projection onto phi_L.  Read it from the data-side tail:
            phi_L = v / np.linalg.norm(v)                        # converged left slow eigenvector
            phi_R = evecd[:, od[0]].real; phi_R /= np.linalg.norm(phi_R)
            # observed gap deviation at r=7: d[m] = R_7(2m) - leading reconstruction, projected onto phi_L
            dvec = np.array([Rtab[7][m] for m in idx], dtype=float)
            dvec = dvec - uR * (uL @ dvec)                       # remove leading (uniform) part
            proj = float(phi_L @ dvec)                           # data's slow-mode amplitude
            lam_local = LAM[16] / LAM[15]                         # observed local decay of Lambda near r=16
            print(f"       TSW-D: REAL slow mode lam={slow.real:+.4f} (obs local Lambda_16/Lambda_15={lam_local:.4f}). "
                  f"Real spectrum => FINITELY many sign changes => g_r EVENTUALLY single-signed.")
            print(f"              observed slow-mode amplitude (data proj on phi_L) sign = {'+' if proj>0 else '-'} "
                  f"({proj:+.3e}); Lambda_r POSITIVE + decaying at ~lam since r=7 (S2: through r~35).")
            print(f"              => IF this positive slow mode dominates: g_r stays + => eps rises => S_inf~0.477 (NOT 7/15).")
            print(f"              => 7/15 requires a SLOWER real mode (|lam|>{abs(slow.real):.3f}) with NEGATIVE coupling to")
            print(f"                 overtake and force ONE more (+ -> -) crossing (the r~36 rollover). Sign of that")
            print(f"                 coupling = the crux; NOT decided here (Wilson's derivation).")
        else:
            phinc = abs(math.degrees(cmath.phase(slow)))
            print(f"       TSW-E: COMPLEX slow mode. phase incr/level = {phinc:.3f} deg => period {360/phinc:.2f} levels.")
        print()
    print("   NOTE (truncation caveat, R29-B): |lam2| value drifts with D (0.92->1.00, D=4..10); its REALITY is")
    print("   robust across all reachable (D,r), but if lam2 is a continuous-spectrum EDGE the value is not a")
    print("   discrete eigenvalue. Direction convergence (overlap->1) is the robust real/complex read.\n")


def main():
    print("# PROBE TSW -- transport-adjoint direction convergence (is the slow mode REAL or COMPLEX?).\n")
    okA = tsw_A()
    if not okA:
        print("TSW-A gate FAILED -- operator identity mis-built. STOP."); return
    # build gap correlations R_r(2m) via certified build_nu (r<=7)
    nus = build_nu(lam, 7)
    NH = {rr: nu_hat(nus[rr], 3 ** (rr + 1)) for rr in range(0, 8)}
    MM = 14
    Rtab = {rr: {m: R_of_d(NH[rr], 3 ** (rr + 1), 2 * m).real for m in range(-MM, MM + 1)} for rr in range(0, 8)}
    tsw_B_C_D_E(Rtab)


if __name__ == "__main__":
    main()
