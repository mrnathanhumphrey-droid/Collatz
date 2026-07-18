"""
PROBE 2c2 -- first corrector (parity-sector dressing of the bad (odd,v3=0) cell).
Direct/entry algebra + a scalar line search. No eigensolves. Label ALGEBRAIC vs STATISTICAL.

G1: every (odd,v3=0) tower state has row survival in {1/9,4/9} EXACTLY; 4/9 <=> a+gamma==0 mod 3.
G2: every (even,v3=0) tower state has survival exactly 2/9.
CORRECTOR: h_beta(x) = h0(cell)*(1 + beta*sigma(x)*1_{(odd,v0)}),  sigma=+1 iff a+gamma==0 mod 3.
   Line search beta in (-1,1): Collatz-Wielandt bracket [min,max] of (M^T h_beta)/h_beta at L=2,3;
   optimal beta*, residual width, and WHICH cell/structure carries the residual (histogram).
   Pre-reg (structural only): width strictly shrinks; residual carrier no longer the (odd,v0) dichotomy.
"""
import numpy as np
from fractions import Fraction
from collections import defaultdict

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def build_tower(q, L, lam=0.5):
    qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    gam = np.array([s[2] for s in states])
    a_arr = np.array([s[0] for s in states])
    erho = np.array([dl[(s[1] * pow(s[0], -1, qL)) % qL] for s in states])
    tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsc()
    return Mt, gam[tw], erho[tw], a_arr[tw], D, qL


def cell_of(e, g):
    return (int(e % 2), 0 if g % 3 != 0 else 1)


def gates(L):
    Mt, gam, erho, a_arr, D, qL = build_tower(3, L)
    nt = Mt.shape[0]
    surv = np.asarray(Mt.sum(axis=0)).ravel()
    ov0 = (erho % 2 == 1) & (gam % 3 != 0)
    ev0 = (erho % 2 == 0) & (gam % 3 != 0)
    # G1
    s_ov0 = surv[ov0]; a_ov0 = a_arr[ov0]; g_ov0 = gam[ov0]
    in_set = np.all([np.isclose(v, 1/9) or np.isclose(v, 4/9) for v in s_ov0])
    sig_pred = ((a_ov0 + g_ov0) % 3 == 0)
    is49 = np.isclose(s_ov0, 4/9)
    dichotomy_ok = np.array_equal(is49, sig_pred)
    log(f"\n   G1 (odd,v3=0) L={L}: survival in {{1/9,4/9}} exactly: {in_set}  "
        f"[values {sorted(set(np.round(s_ov0,6)))}]; 4/9 <=> a+gamma==0 mod3: {dichotomy_ok}  "
        f"({'PASS' if in_set and dichotomy_ok else 'FAIL'})   (frac at 4/9 = {is49.mean():.4f})")
    # G2
    s_ev0 = surv[ev0]
    g2ok = np.allclose(s_ev0, 2/9)
    log(f"   G2 (even,v3=0) L={L}: survival == 2/9 exactly: {g2ok}  "
        f"[max|surv-2/9| = {np.max(np.abs(s_ev0 - 2/9)):.2e}]  ({'PASS' if g2ok else 'FAIL'})")
    return in_set and dichotomy_ok, g2ok


def corrector(L, partner):
    Mt, gam, erho, a_arr, D, qL = build_tower(3, L)
    nt = Mt.shape[0]
    par = erho % 2; v0 = (gam % 3 != 0)
    h0 = np.empty(nt)
    h0[(par == 0) & v0] = 2/3; h0[(par == 0) & ~v0] = 5/3
    h0[(par == 1) & v0] = 5/6; h0[(par == 1) & ~v0] = 4/3
    ov0 = (par == 1) & v0
    sig = np.zeros(nt)
    sig[ov0] = np.where((a_arr[ov0] + gam[ov0]) % 3 == 0, 1.0, -1.0)
    MtT = Mt.T.tocsr()

    def bracket(beta):
        hb = h0 * (1 + beta * sig)
        r = MtT.dot(hb) / hb
        return r.min(), r.max(), r

    betas = np.linspace(-0.98, 0.98, 1961)
    widths = np.array([(lambda mm: mm[1] - mm[0])(bracket(b)) for b in betas])
    i = int(np.argmin(widths))
    bstar = betas[i]
    # refine
    fine = np.linspace(bstar - 0.01, bstar + 0.01, 2001)
    fw = np.array([(lambda mm: mm[1] - mm[0])(bracket(b)) for b in fine])
    bstar = fine[int(np.argmin(fw))]
    m0, M0, r0 = bracket(0.0)
    mb, Mb, rb = bracket(bstar)
    log(f"\n   CORRECTOR L={L}:  baseline (beta=0) [{m0:.6f},{M0:.6f}] width {M0-m0:.6f}")
    log(f"      line search: beta* = {bstar:+.4f} (~{Fraction(bstar).limit_denominator(50)})  "
        f"bracket [{mb:.6f},{Mb:.6f}]  width {Mb-mb:.6f}   "
        f"shrink {(M0-m0)/(Mb-mb):.3f}x   1/3 inside: {mb <= 1/3 <= Mb}  partner inside: {mb <= partner <= Mb}")
    # a few beta samples
    log(f"      width(beta): " + "  ".join(f"{b:+.2f}:{bracket(b)[1]-bracket(b)[0]:.4f}"
                                           for b in [-0.5, -0.25, 0.0, 0.25, bstar, 0.5]))
    # residual carrier: histogram at beta*
    log(f"      residual carrier at beta* (per-cell ratio spread, and is (odd,v0) dichotomy gone?):")
    for c in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        mask = np.array([cell_of(erho[i2], gam[i2]) == c for i2 in range(nt)])
        rr = rb[mask]
        uniq = sorted(set(np.round(rr, 6)))
        uu = [str(Fraction(v).limit_denominator(200)) for v in uniq]
        log(f"        cell {c}: [{rr.min():.6f}, {rr.max():.6f}]  spread {rr.max()-rr.min():.6f}  "
            f"distinct ratios {uu if len(uu)<=6 else uu[:6]+['...('+str(len(uu))+')']}")
    # within (odd,v0), is the residual still split by sigma?
    rr = rb[ov0]; sg = sig[ov0]
    lo = rr[sg < 0]; hi = rr[sg > 0]
    log(f"        (odd,v0) by sigma: sigma=-1 -> [{lo.min():.6f},{lo.max():.6f}], "
        f"sigma=+1 -> [{hi.min():.6f},{hi.max():.6f}]  "
        f"(dichotomy {'GONE (overlap)' if hi.min() <= lo.max() and lo.min() <= hi.max() else 'still split'})")
    # which cell now carries the max residual
    spreads = {c: (lambda m: (rb[m].max()-rb[m].min()))(np.array([cell_of(erho[i2],gam[i2])==c for i2 in range(nt)]))
               for c in [(0,0),(0,1),(1,0),(1,1)]}
    carrier = max(spreads, key=spreads.get)
    log(f"        => residual max-spread carrier = cell {carrier} (spread {spreads[carrier]:.6f})")
    return M0 - m0, Mb - mb, bstar, carrier


def main():
    log("# PROBE 2c2 -- first corrector (parity dressing of the bad cell). Direct, no eigensolves.")
    log("\n## GATES")
    for L in [2, 3]:
        gates(L)
    log("\n## CORRECTOR (line search)")
    r2 = corrector(2, 0.346827)
    r3 = corrector(3, 0.333236)
    log(f"\n## SUMMARY: width {r2[0]:.4f}->{r2[1]:.4f} (L2, beta*={r2[2]:+.3f}, carrier {r2[3]}); "
        f"{r3[0]:.4f}->{r3[1]:.4f} (L3, beta*={r3[2]:+.3f}, carrier {r3[3]})")
    with open("logs/probe_phase2c2_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
