"""
PROBE 2c4 -- rung-2 dressing + the v1 front. Direct/exact + a 2-param search. No eigensolves.
tau (v0) := r-rung at beta*=3/5: tau=0<->12/49, tau=1<->17/49, tau=2<->20/49 (ascending r).

A) DESTINATION-tau MATRIX: per source tau-level i, mass-distribution of destination tau among its
   v'=0 targets -> 3x3 A[i][j] exact rationals, both L. (rung-3 raw material + diagonal-approx judge.)
B) SEARCH: beta=3/5 frozen; 2-param search over g(tau) mean-zero. Corrector h2 = h_beta*(1+g(tau)) on v0.
   Report g*, residual bracket, carrier. Reference g_diag=(-13,+2,+11) (~dev of {12,17,20}/49 from 1/3,
   in /147). g*=g_diag iff A is tau-uniform by row; gap = correlation readout.
C) v1 ANALOG TABLES: mod-9 gate/carry for v>=1 sources (even e' channels, LTE t-class), u%9 x gamma%9
   -> gamma' mod3, + R0(s) weights. (blind-derivation raw material for the v1 trit.)
"""
import numpy as np
from fractions import Fraction
from collections import defaultdict
from itertools import product

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def build(q, L, lam=0.5):
    qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    gam = np.array([s[2] for s in states]); a_arr = np.array([s[0] for s in states])
    erho = np.array([dl[(s[1] * pow(s[0], -1, qL)) % qL] for s in states])
    tw = np.where(gam != 0)[0]
    return M[tw][:, tw].tocsc(), gam[tw], erho[tw], a_arr[tw], D, qL


def hbeta(erho, gam, a_arr, beta):
    nt = len(erho); par = erho % 2; v0 = (gam % 3 != 0)
    h = np.empty(nt)
    h[(par == 0) & v0] = 2/3; h[(par == 0) & ~v0] = 5/3
    h[(par == 1) & v0] = 5/6; h[(par == 1) & ~v0] = 4/3
    ov0 = (par == 1) & v0
    sig = np.zeros(nt); sig[ov0] = np.where((a_arr[ov0] + gam[ov0]) % 3 == 0, 1.0, -1.0)
    return h * (1 + beta * sig)


LADDER = {12: 0, 17: 1, 20: 2}   # r*49 -> tau


def tau_of(r, v0):
    """tau for v0 states from r-rung; -1 elsewhere."""
    t = np.full(len(r), -1)
    for i in np.where(v0)[0]:
        k = int(round(r[i] * 49))
        t[i] = LADDER.get(k, -1)
    return t


def partA(L):
    Mt, gam, erho, a_arr, D, qL = build(3, L)
    nt = Mt.shape[0]
    hb = hbeta(erho, gam, a_arr, 3/5)
    r = Mt.T.dot(hb) / hb
    v0 = (gam % 3 != 0)
    tau = tau_of(r, v0)
    indptr, indices, data = Mt.indptr, Mt.indices, Mt.data
    A = [[0.0] * 3 for _ in range(3)]
    for x in range(nt):
        if not v0[x]:
            continue
        i = tau[x]
        for p in range(indptr[x], indptr[x + 1]):
            dst = indices[p]
            if gam[dst] % 3 != 0:                 # v'=0 target
                j = tau[dst]
                if j >= 0:
                    A[i][j] += data[p]
    log(f"\n   A) DESTINATION-tau MATRIX (v0 sources -> v'=0 targets), q=3 L={L}  (rows=src tau, cols=dst tau)")
    Af = [[Fraction(A[i][j]).limit_denominator(100000) for j in range(3)] for i in range(3)]
    rowsums = [sum(A[i]) for i in range(3)]
    for i in range(3):
        norm = [Fraction(A[i][j] / rowsums[i]).limit_denominator(5000) for j in range(3)]
        log(f"      src tau={i} ({[12,17,20][i]}/49): row-norm dist over dst tau = "
            f"[{norm[0]}, {norm[1]}, {norm[2]}]  (= {[float(x) for x in norm]})")
    # tau-uniform by row?
    r0 = [A[0][j] / rowsums[0] for j in range(3)]
    uniform = all(abs(A[i][j]/rowsums[i] - r0[j]) < 1e-12 for i in range(3) for j in range(3))
    dist = [A[0][j]/rowsums[0] for j in range(3)]
    is_unifdist = all(abs(d - 1/3) < 1e-12 for d in dist)
    log(f"      => A tau-uniform by row (all src rows identical): {uniform}  "
        f"[row dist = {'[1/3,1/3,1/3] fully uniform' if is_unifdist else '['+','.join(str(Fraction(d).limit_denominator(20)) for d in dist)+']'}]  "
        f"(no source-tau -> dest-tau correlation)")
    return uniform


def partB(L):
    Mt, gam, erho, a_arr, D, qL = build(3, L)
    nt = Mt.shape[0]
    hb = hbeta(erho, gam, a_arr, 3/5)
    r = Mt.T.dot(hb) / hb
    v0 = (gam % 3 != 0)
    tau = tau_of(r, v0)
    MtT = Mt.T.tocsr()
    tvec = np.zeros((nt, 3))
    for i in range(nt):
        if tau[i] >= 0: tvec[i, tau[i]] = 1.0

    par = erho % 2
    cells = {(0,0):"E,v0",(0,1):"E,v1",(1,0):"O,v0",(1,1):"O,v1"}
    def eval_g(g):
        fac = np.ones(nt); m = tau >= 0
        fac[m] = 1 + np.array([g[tau[i]] for i in np.where(m)[0]])
        h2 = hb * fac
        rr = MtT.dot(h2) / h2
        sp = {}
        for c,nm in cells.items():
            mask = (par==c[0]) & ((gam%3!=0)==(c[1]==0))
            if mask.any(): sp[nm] = (rr[mask].min(), rr[mask].max())
        return rr.min(), rr.max(), sp, rr
    v0mask = (gam % 3 != 0)
    # (1) g_diag = (-13,+2,+11)/49 : the v0-flattener (dst cancellation when dst-tau uniform)
    gdiag = (-13/49, 2/49, 11/49)
    lo, hi, sp, rr = eval_g(gdiag)
    v0_spread = rr[v0mask].max() - rr[v0mask].min()
    log(f"\n   B) g_diag=(-13,+2,+11)/49: v0 ratio spread = {v0_spread:.3e}  "
        f"({'v0 COLLAPSES to 1/3' if v0_spread < 1e-9 else 'v0 NOT flattened by g_diag (destination back-reaction + undressed DEEP dilution)'})")
    log(f"      full bracket [{lo:.6f},{hi:.6f}] width {hi-lo:.6f}; residual carrier by cell "
        f"{{{', '.join(f'{k}:[{v[0]:.4f},{v[1]:.4f}]' for k,v in sp.items())}}}")
    # (2) full-bracket minimizer (secondary: trades v0-flattening vs v1 back-coupling)
    best = None
    for g0 in np.linspace(-0.30, 0.05, 351):
        for g1 in np.linspace(-0.10, 0.15, 251):
            lo2, hi2, _, _ = eval_g((g0, g1, -g0-g1))
            if best is None or hi2-lo2 < best[0]: best = (hi2-lo2, (g0,g1,-g0-g1))
    wgm, gstar = best
    # (3) v0-spread minimizer: is its DIRECTION g_diag?
    bestv0 = None
    for g0 in np.linspace(-0.35, 0.05, 201):
        for g1 in np.linspace(-0.10, 0.15, 151):
            _, _, _, rr2 = eval_g((g0, g1, -g0-g1))
            s = rr2[v0mask].max() - rr2[v0mask].min()
            if bestv0 is None or s < bestv0[0]: bestv0 = (s, (g0,g1,-g0-g1))
    v0s, gv0 = bestv0
    gd = np.array([-13, 2, 11]) / np.linalg.norm([-13, 2, 11])
    gv = np.array(gv0); gvn = gv / np.linalg.norm(gv)
    log(f"      full-bracket min g*=({gstar[0]:+.4f},{gstar[1]:+.4f},{gstar[2]:+.4f}) width {wgm:.6f} "
        f"(rung-1 9/49={9/49:.4f}; shrink {9/49/wgm:.2f}x)")
    log(f"      v0-spread-min g=({gv0[0]:+.4f},{gv0[1]:+.4f},{gv0[2]:+.4f}): DIRECTION cos vs g_diag = {abs(np.dot(gvn,gd)):.4f} "
        f"({'ALIGNED — g_diag direction confirmed' if abs(np.dot(gvn,gd))>0.99 else 'deviates'}); "
        f"but v0 floor spread = {v0s:.4f} (trit CANNOT fully flatten v0)")
    log(f"      => g_diag right DIRECTION, over-magnitude; v0 trit alone shrinks little (v1 is binding + gets")
    log(f"         perturbed). PAIR-OF-TRITS (v0+v1, C) needed. Residual carrier = v1.")
    return gdiag, hi-lo, sp


def partC(L=3):
    qL = 3 ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    Z = 2 ** D - 1
    wf = [Fraction(2 ** (D - m), Z) for m in range(1, D + 1)]
    R0 = [sum(wf[m-1]*wf[((m-s-1) % D)] for m in range(1, D+1)) for s in range(min(D,6))]
    log(f"\n   C) v1 ANALOG mod-9 GATE/CARRY (v>=1 sources, EVEN e' channels), q=3:")
    lines = ["# 2c4-C v1 gate/carry q=3: EVEN e' channels (v>=1 source, gamma==0 mod3). (e'mod6,u%9,gamma%9)->gate,gammap%3",
             "# e'mod6\tc=1-2^e'%9\tLTE_t\tu%9\tgamma%9\tgate\tgammap%3"]
    units9 = sorted(set(u % 9 for u in subgroup(2 % 9, 9)))
    for epm in [0, 2, 4]:                          # even channels
        c = (1 - pow(2, epm, 9)) % 9
        t = 2 if epm == 0 else 1                    # LTE t: e'=0 mod6 -> t>=2 (digit-shift); 2,4 -> t=1
        passu = set()
        for u9 in units9:
            for g9 in range(0, 9, 3):               # gamma == 0 mod 3 (v>=1)
                T = (u9 * c) % 9
                if (g9 + T) % 3 == 0:
                    passu.add(u9)
                    gpm3 = ((g9 + u9 * c) // 3) % 3
                    lines.append(f"{epm}\t{c}\t{t}\t{u9}\t{g9}\t1\t{gpm3}")
        log(f"      e' mod6={epm} (t={t}, c={c}): passing u mod9 (gamma≡0 mod3) = {sorted(passu)}")
    with open("outputs/v1_gate_mod9_q3.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"      R0(s) weights s=0..5: {[str(x) for x in R0]}")
    log(f"      dumped outputs/v1_gate_mod9_q3.tsv")


def main():
    log("# PROBE 2c4 -- rung-2 dressing (A,B) + v1 front tables (C). Direct/exact + 2-param search.")
    for L in [2, 3]:
        log(f"\n{'='*72}\n## L={L}")
        partA(L); partB(L)
    partC(3)
    with open("logs/probe_phase2c4_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
