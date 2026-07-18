"""
L=4 JOINT (g0,g1) CORRECTOR SEARCH on GPU (CuPy SpMV). Cap-rung L-flow test.
Build M_tower(3,4) on CPU (scipy), move M^T to GPU, run the joint search there.
Mirrors probes/probe_phase2c6.py exactly; reports width, optimal g, residual carrier + resolution,
and the L=2/3/4 width sequence (the breathing-with-L verdict).
"""
import sys, time, json
import numpy as np
from fractions import Fraction
from collections import defaultdict
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

import cupy as cp
import cupyx.scipy.sparse as csp

LADDER = {12: 0, 17: 1, 20: 2}
CLS = {"D9": 0, "U+": 1, "U-": 2}

def cls3(a, gam):
    g3 = (gam % 9) // 3
    if g3 == 0: return 0            # D9
    a3 = a % 3
    return 1 if a3 == g3 else 2     # U+ / U-

def build(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    t0 = time.time()
    M, idx, n = build_M_gen(q, L, 2, raw)
    print(f"  build_M_gen L={L}: n={n} nnz={M.nnz} {time.time()-t0:.1f}s", flush=True)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    gam = np.array([s[2] for s in states]); a_arr = np.array([s[0] for s in states])
    erho = np.array([dl[(s[1] * pow(s[0], -1, qL)) % qL] for s in states])
    tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsr()
    return Mt, gam[tw], erho[tw], a_arr[tw], D, qL

def hbeta(erho, gam, a_arr, beta):
    nt = len(erho); par = erho % 2; v0 = (gam % 3 != 0)
    h = np.empty(nt)
    h[(par == 0) & v0] = 2/3; h[(par == 0) & ~v0] = 5/3
    h[(par == 1) & v0] = 5/6; h[(par == 1) & ~v0] = 4/3
    ov0 = (par == 1) & v0
    sig = np.zeros(nt); sig[ov0] = np.where((a_arr[ov0] + gam[ov0]) % 3 == 0, 1.0, -1.0)
    return h * (1 + beta * sig)

def main():
    L = 4
    import os, scipy.sparse as _sp
    CACHE = os.path.expanduser("~/l4_cache")
    if os.path.exists(CACHE + "_Mt.npz"):
        Mt = _sp.load_npz(CACHE + "_Mt.npz")
        z = np.load(CACHE + "_vec.npz")
        gam, erho, a_arr, D, qL = z["gam"], z["erho"], z["a_arr"], int(z["D"]), int(z["qL"])
        print(f"  loaded cache: tower={Mt.shape[0]} nnz={Mt.nnz}", flush=True)
    else:
        Mt, gam, erho, a_arr, D, qL = build(L)
        _sp.save_npz(CACHE + "_Mt.npz", Mt.tocsr())
        np.savez(CACHE + "_vec.npz", gam=gam, erho=erho, a_arr=a_arr, D=D, qL=qL)
        print(f"  cached build to {CACHE}", flush=True)
    nt = Mt.shape[0]
    hb = hbeta(erho, gam, a_arr, 3/5)
    # GPU operator = M^T (flow)
    MtT_gpu = csp.csr_matrix(Mt.T.tocsr().astype(np.float64))
    hb_gpu = cp.asarray(hb)
    r0 = (MtT_gpu.dot(hb_gpu) / hb_gpu).get()
    v0 = (gam % 3 != 0); v1 = ~v0
    tau = np.full(nt, -1)
    for i in np.where(v0)[0]:
        tau[i] = LADDER.get(int(round(r0[i] * 49)), -1)
    cl = np.full(nt, -1)
    for i in np.where(v1)[0]:
        cl[i] = cls3(int(a_arr[i]), int(gam[i]))
    print(f"  L=4 tower={nt} v0={v0.sum()} v1={v1.sum()}  tau-mapped={(tau>=0).sum()}/{v0.sum()} "
          f"cls-mapped={(cl>=0).sum()}/{v1.sum()}", flush=True)
    tau_g = cp.asarray(tau); cl_g = cp.asarray(cl)
    v0_g = cp.asarray(v0); v1_g = cp.asarray(v1)

    def width(x):
        g0 = cp.asarray([x[0], x[1], -x[0]-x[1], 0.0])   # index -1 -> 0 dressing
        g1 = cp.asarray([x[2], x[3], -x[2]-x[3], 0.0])
        fac = cp.ones(nt)
        fac = cp.where(v0_g, 1 + g0[tau_g], fac)
        fac = cp.where(v1_g, 1 + g1[cl_g], fac)
        h = hb_gpu * fac
        r = MtT_gpu.dot(h) / h
        return float(r.max() - r.min())

    def rvec(x):
        g0 = cp.asarray([x[0], x[1], -x[0]-x[1], 0.0])
        g1 = cp.asarray([x[2], x[3], -x[2]-x[3], 0.0])
        fac = cp.ones(nt)
        fac = cp.where(v0_g, 1 + g0[tau_g], fac)
        fac = cp.where(v1_g, 1 + g1[cl_g], fac)
        h = hb_gpu * fac
        return (MtT_gpu.dot(h) / h).get()

    w_base = width((0,0,0,0))
    print(f"  baseline width L=4 = {w_base:.6f}  (rung-1 9/49={9/49:.6f})", flush=True)

    # coarse grid (moderate) + seed from L=3 optimum g0~(-0.108,0,0.108), g1~0
    t0 = time.time()
    gr = np.linspace(-0.20, 0.20, 11)
    best = None
    for a in gr:
        for b in gr:
            for c in gr:
                for d in gr:
                    w = width((a,b,c,d))
                    if best is None or w < best[0]: best = (w, (a,b,c,d))
    # also test the L=3 seed explicitly
    for seed in [(-0.108,0.0,0.0,0.0), (-0.1084,0.0002,0.0,0.0)]:
        w = width(seed)
        if w < best[0]: best = (w, seed)
    # coordinate-descent refine
    x = list(best[1]); w = best[0]; step = 0.04
    for _ in range(80):
        improved = False
        for k in range(4):
            for dxt in (step, -step):
                y = x[:]; y[k] += dxt
                wy = width(y)
                if wy < w - 1e-15: x, w = y, wy; improved = True
        if not improved:
            step *= 0.5
            if step < 1e-6: break
    print(f"  search {time.time()-t0:.1f}s", flush=True)
    g0 = (x[0], x[1], -x[0]-x[1]); g1 = (x[2], x[3], -x[2]-x[3])
    # v0-alone reference
    best0 = None
    for a in np.linspace(-0.30,0.05,71):
        for b in np.linspace(-0.10,0.15,51):
            wv = width((a,b,0,0))
            if best0 is None or wv < best0[0]: best0 = (wv,(a,b))
    r = rvec(x)
    lo, hi = float(r.min()), float(r.max())
    # residual carrier by cell
    par = erho % 2
    cells = {(0,0):"E,v0",(0,1):"E,v1",(1,0):"O,v0",(1,1):"O,v1"}
    cellw = {}
    for (p,vk),nm in cells.items():
        m = (par==p) & (v1 if vk==1 else v0)
        if m.any(): cellw[nm] = float(r[m].max()-r[m].min())
    carrier = max(cellw, key=cellw.get)
    # residual resolution: gamma mod {9,27,81}, a mod {9,27}
    reso = {}
    for nm,kf in [("a9_g9_e6", lambda i:(int(a_arr[i]%9),int(gam[i]%9),int(erho[i]%6))),
                  ("a9_g27_e6", lambda i:(int(a_arr[i]%9),int(gam[i]%27),int(erho[i]%6))),
                  ("a9_g81_e6", lambda i:(int(a_arr[i]%9),int(gam[i]%81),int(erho[i]%6))),
                  ("a27_g81_e6", lambda i:(int(a_arr[i]%27),int(gam[i]%81),int(erho[i]%6)))]:
        km = defaultdict(set)
        for i in range(nt): km[kf(i)].add(round(float(r[i]),8))
        bad = sum(1 for v in km.values() if len(v)>1)
        reso[nm] = ("wd" if bad==0 else f"bad={bad}")

    out = dict(L=4, tower=int(nt), w_base=w_base, w_v0=best0[0],
               w_joint=w, g0=[round(v,5) for v in g0], g1=[round(v,5) for v in g1],
               shrink_rung1=round(w_base/w,4), shrink_v0=round(best0[0]/w,4),
               carrier=carrier, cellw={k:round(v,5) for k,v in cellw.items()},
               reso=reso, bracket=[round(lo,6),round(hi,6)])
    print("RESULT_JSON " + json.dumps(out), flush=True)
    print(f"\n  === L=4 JOINT ===", flush=True)
    print(f"  g0={out['g0']} g1={out['g1']}", flush=True)
    print(f"  width joint={w:.6f} (shrink rung1 {out['shrink_rung1']}x, v0-alone {out['w_v0']:.6f} -> {out['shrink_v0']}x)", flush=True)
    print(f"  carrier={carrier} cellw={out['cellw']}", flush=True)
    print(f"  residual resolution: {reso}", flush=True)

if __name__ == "__main__":
    main()
