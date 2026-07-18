"""
PROBE 2c6 -- THE JOINT (g0,g1) CORRECTOR SEARCH + REVEAL. beta*=3/5 frozen. Direct/exact matvecs.
(Wilson's blind joint-optimum ruling is ON THE RECORD: decoupled = reference/direction only; the two
mod-27 couplings named as deviation owners. This is the SEARCH he asked me to reveal.)

Corrector:  h = h_beta * (1 + g0[tau] on v0 states) * (1 + g1[cls] on v>=1 states)
  g0 mean-zero on tau in {0,1,2} (v0 trit; tau via LADDER at beta*);  2 dof.
  g1 mean-zero on cls in {D9,U+,U-} (v1 trit);                        2 dof.
Bracket = [min,max] of r=(M^T h)/h over the tower. Minimize width.

Answers, in Wilson's order:
  Q1 real joint shrink (width vs rung-1 9/49 and vs v0-alone 1.14x);
  Q2 which structure carries the residual (cell/key spanning the surviving bracket);
  Q3 *** does anything in the post-joint residual finally MOVE with L *** (width, optimal g, residual key,
     and whether the residual now needs mod-27 = the coupling/cap-rung shell the shell-picture predicts).
"""
import numpy as np
from fractions import Fraction
from collections import defaultdict
from itertools import product

from probe_phase2a_q2b_q6 import build_M_gen, subgroup
from probe_phase2c5 import build, hbeta, cls3

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))

LADDER = {12: 0, 17: 1, 20: 2}
CLS = {"D9": 0, "U+": 1, "U-": 2}


def setup(L):
    Mt, gam, erho, a_arr, D, qL = build(3, L)
    nt = Mt.shape[0]
    hb = hbeta(erho, gam, a_arr, 3/5)
    r0 = Mt.T.dot(hb) / hb
    v0 = (gam % 3 != 0); v1 = ~v0
    # tau on v0 from the baseline ladder {12,17,20}/49
    tau = np.full(nt, -1)
    for i in np.where(v0)[0]:
        tau[i] = LADDER.get(int(round(r0[i] * 49)), -1)
    # cls on v1
    cl = np.full(nt, -1)
    for i in np.where(v1)[0]:
        lbl, b, eta = cls3(int(a_arr[i]), int(gam[i]))
        cl[i] = CLS[lbl]
    MtT = Mt.T.tocsr()
    par = erho % 2
    return dict(MtT=MtT, hb=hb, v0=v0, v1=v1, tau=tau, cl=cl, gam=gam, erho=erho,
                a_arr=a_arr, par=par, nt=nt, D=D, r0=r0)


def make_h(S, g0, g1):
    fac = np.ones(S["nt"])
    m0 = S["v0"]; fac[m0] *= 1 + np.array([g0[t] for t in S["tau"][m0]])
    m1 = S["v1"]; fac[m1] *= 1 + np.array([g1[c] for c in S["cl"][m1]])
    return S["hb"] * fac


def width_of(S, x):
    g0 = np.array([x[0], x[1], -x[0] - x[1]])
    g1 = np.array([x[2], x[3], -x[2] - x[3]])
    h = make_h(S, g0, g1)
    r = S["MtT"].dot(h) / h
    return r.max() - r.min()


def rvec_of(S, x):
    g0 = np.array([x[0], x[1], -x[0] - x[1]])
    g1 = np.array([x[2], x[3], -x[2] - x[3]])
    h = make_h(S, g0, g1)
    r = S["MtT"].dot(h) / h
    return r


def search(S, coarse=13, span=0.30):
    # coarse grid seed on 4 params (mean-zero param'd by 2 each)
    gr = np.linspace(-span, span, coarse)
    best = None
    for a in gr:
        for b in gr:
            for c in gr:
                for d in gr:
                    w = width_of(S, (a, b, c, d))
                    if best is None or w < best[0]:
                        best = (w, (a, b, c, d))
    # local refine (coordinate descent + shrink)
    x = list(best[1]); w = best[0]; step = span / (coarse - 1)
    for _ in range(60):
        improved = False
        for k in range(4):
            for dxt in (step, -step):
                y = x[:]; y[k] += dxt
                wy = width_of(S, y)
                if wy < w - 1e-15:
                    x, w = y, wy; improved = True
        if not improved:
            step *= 0.5
            if step < 1e-6:
                break
    return w, tuple(x)


def characterize(S, x, tag):
    r = rvec_of(S, x)
    lo, hi = r.min(), r.max()
    # residual carrier: cell (parity, v-class) min/max spans
    cells = {(0,0):"E,v0",(0,1):"E,v1",(1,0):"O,v0",(1,1):"O,v1"}
    log(f"      [{tag}] bracket=[{lo:.6f},{hi:.6f}] width={hi-lo:.6f}")
    span_cell = None; span_w = -1
    for (p, vk), nm in cells.items():
        mask = (S["par"] == p) & (S["v1"] if vk == 1 else S["v0"])
        if mask.any():
            cw = r[mask].max() - r[mask].min()
            log(f"         cell {nm}: [{r[mask].min():.5f},{r[mask].max():.5f}] width {cw:.5f}")
            if cw > span_w: span_w, span_cell = cw, nm
    log(f"         => residual carrier (widest cell) = {span_cell} (width {span_w:.5f})")
    # does the residual key through the trit resolution or need mod-27 ?
    for keyname, kf in [("trit key (e%6, v0:tau / v1:cls)",
                         lambda i: (int(S["erho"][i] % 6), 0, int(S["tau"][i])) if S["v0"][i]
                                    else (int(S["erho"][i] % 6), 1, int(S["cl"][i]))),
                        ("mod-9 key (a%9,g%9,e%6)",
                         lambda i: (int(S["a_arr"][i] % 9), int(S["gam"][i] % 9), int(S["erho"][i] % 6))),
                        ("mod-27 key (a%9,g%27,e%6)",
                         lambda i: (int(S["a_arr"][i] % 9), int(S["gam"][i] % 27), int(S["erho"][i] % 6)))]:
        km = defaultdict(set)
        for i in range(S["nt"]):
            km[kf(i)].add(round(float(r[i]), 8))
        wd = all(len(v) == 1 for v in km.values())
        nv = len({list(v)[0] for v in km.values() if len(v) == 1})
        log(f"         residual well-defined on {keyname}: {wd}" + (f" ({nv} vals)" if wd else ""))
    return hi - lo, span_cell, r


def main():
    log("# PROBE 2c6 -- JOINT (g0,g1) corrector search + reveal. beta*=3/5 frozen. Direct matvec.")
    res = {}
    for L in [2, 3]:
        log(f"\n{'='*74}\n## L={L}")
        S = setup(L)
        # reference points
        w_base = width_of(S, (0, 0, 0, 0))
        log(f"   baseline (no trit, beta*=3/5 only): width = {w_base:.6f}  (= rung-1 9/49 = {9/49:.6f})")
        # v0-alone (2c4): best over g0 only
        best0 = None
        gr = np.linspace(-0.30, 0.05, 141)
        for a in gr:
            for b in np.linspace(-0.10, 0.15, 101):
                w = width_of(S, (a, b, 0, 0))
                if best0 is None or w < best0[0]: best0 = (w, (a, b, 0, 0))
        log(f"   v0-trit ALONE best: width = {best0[0]:.6f} (shrink {w_base/best0[0]:.3f}x)  [matches 2c4 ~1.14x]")
        # JOINT
        wj, xj = search(S)
        g0 = (xj[0], xj[1], -xj[0]-xj[1]); g1 = (xj[2], xj[3], -xj[2]-xj[3])
        log(f"\n   *** JOINT (g0,g1) optimum ***")
        log(f"   g0(tau)  = ({g0[0]:+.4f},{g0[1]:+.4f},{g0[2]:+.4f})   [v0 trit dressing]")
        log(f"   g1(cls)  = ({g1[0]:+.4f},{g1[1]:+.4f},{g1[2]:+.4f})   [v1 trit: D9,U+,U-]")
        log(f"   Q1 JOINT WIDTH = {wj:.6f}   shrink vs rung-1 = {w_base/wj:.3f}x   vs v0-alone = {best0[0]/wj:.3f}x")
        wj2, carrier, rj = characterize(S, xj, "JOINT")
        res[L] = dict(w_base=w_base, w_v0=best0[0], wj=wj, xj=xj, g0=g0, g1=g1, carrier=carrier)
    # ------- Q3: L-motion -------
    log(f"\n{'='*74}\n## Q3 -- DOES ANYTHING MOVE WITH L?")
    a, b = res[2], res[3]
    log(f"   baseline width:  L=2 {a['w_base']:.6f}  L=3 {b['w_base']:.6f}   (delta {b['w_base']-a['w_base']:+.2e})")
    log(f"   v0-alone width:  L=2 {a['w_v0']:.6f}  L=3 {b['w_v0']:.6f}   (delta {b['w_v0']-a['w_v0']:+.2e})")
    log(f"   JOINT width:     L=2 {a['wj']:.6f}  L=3 {b['wj']:.6f}   (delta {b['wj']-a['wj']:+.2e})")
    log(f"   JOINT g0:        L=2 {tuple(round(v,4) for v in a['g0'])}  L=3 {tuple(round(v,4) for v in b['g0'])}")
    log(f"   JOINT g1:        L=2 {tuple(round(v,4) for v in a['g1'])}  L=3 {tuple(round(v,4) for v in b['g1'])}")
    log(f"   residual carrier: L=2 {a['carrier']}  L=3 {b['carrier']}")
    moved = abs(b['wj'] - a['wj']) > 1e-6
    log(f"   => post-joint residual {'MOVES with L' if moved else 'still L-INVARIANT'} "
        f"(|delta width| = {abs(b['wj']-a['wj']):.2e})")
    with open("logs/probe_phase2c6_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
