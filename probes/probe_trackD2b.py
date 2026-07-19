"""
PROBE D2-b RIDER -- flux-weighted histogram of W's TOP SCALE. W = That + c = gamma' - floor(gamma/3).
Top-scale leading trit = floor(W / 3^{L-2}) in {0,1,2} (scale 3^{L-1}); + full W histogram (fine shape);
+ the W mod 3 marginal (Wilson pre-reg: exactly uniform). Flux-weighted along the mean-field stationary flow,
exact rationals, L=2,3.
PRE-REG (fork discriminator, SHAPE): {1/2,1/2,0} => factor two is a top-scale Bernoulli (bulk symbol);
uniform {1/3,1/3,1/3} => boundary-sourced (section-edge analysis).
"""
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
import sympy
from probe_phase2b_G1 import build_cell_transfer
from probe_phase2c0 import build_M_tower_and_coords

def cell_of(e, gam):
    return (e % 2) * 2 + (1 if gam % 3 == 0 else 0)

def perron_left(T):
    M = sympy.Matrix([[sympy.Rational(T[i][j].numerator, T[i][j].denominator) for j in range(4)] for i in range(4)])
    ns = (M.T - sympy.Rational(1, 3) * sympy.eye(4)).nullspace()
    v = ns[0]; v = v / sum(v)
    return [Fr(int(x.p), int(x.q)) for x in v]

def run(L):
    T, casc, cellcount, spread, D, qL = build_cell_transfer(L)
    v = perron_left(T)
    occ = {c: v[c] / cellcount[c] for c in range(4)}
    Mt, states, idx, tw, pos, twcoords, q, qL2, sub, D2, dl, w, two = build_M_tower_and_coords(L)
    wf = [Fr(2 ** (D - d), 2 ** D - 1) for d in range(1, D + 1)]
    scale = 3 ** (L - 2) if L >= 2 else 1
    toptrit = defaultdict(Fr)     # floor(W/3^{L-2}) -> flux
    Whist = defaultdict(Fr)       # full W -> flux
    Wm3 = defaultdict(Fr)         # W mod 3 -> flux (marginal)
    for src_t, (a, e, gam) in enumerate(twcoords):
        dla = dl[a]; os = occ[cell_of(e, gam)]
        for u in sub:
            dlu = dl[u]; da = (dla - dlu) % D; da = da if da != 0 else D
            for s in range(D):
                ep = (e + s) % D; bp = (u * two[ep]) % qL2; Tt = (u - bp) % qL2
                if (gam + Tt) % q != 0: continue
                gp = ((gam + Tt) // q) % qL2
                dst = (u, bp, gp)
                if dst not in idx or idx[dst] not in pos: continue
                db = (da - s) % D; db = db if db != 0 else D
                flux = os * wf[da - 1] * wf[db - 1]
                W = gp - (gam // 3)
                toptrit[W // scale] += flux
                Whist[W] += flux
                Wm3[W % 3] += flux
    def nd(d):
        t = sum(d.values()); return {k: d[k] / t for k in sorted(d)}
    tt = nd(toptrit); wh = nd(Whist); wm = nd(Wm3)
    print(f"\n## L={L}  (top scale 3^{{L-1}}={3**(L-1)}, leading trit = floor(W/{scale}))", flush=True)
    print(f"  TOP-SCALE LEADING TRIT split: " + ", ".join(f"t={k}:{tt[k]}({float(tt[k]):.4f})" for k in tt), flush=True)
    preds = {0: Fr(1,2), 1: Fr(1,2), 2: Fr(0)}
    bern = all(tt.get(k, Fr(0)) == preds[k] for k in (0,1,2)) and set(tt) <= {0,1,2}
    unif = all(tt.get(k, Fr(0)) == Fr(1,3) for k in (0,1,2)) and set(tt) <= {0,1,2}
    print(f"     => {'BERNOULLI {1/2,1/2,0} (bulk symbol)' if bern else ('UNIFORM {1/3,1/3,1/3} (boundary-sourced)' if unif else 'NEITHER pre-reg exactly -- reported as-is')}", flush=True)
    print(f"  W mod 3 marginal (pre-reg: uniform): " + ", ".join(f"{k}:{wm[k]}" for k in wm) +
          f"  ({'UNIFORM' if all(wm.get(k,Fr(0))==Fr(1,3) for k in (0,1,2)) else 'NOT uniform'})", flush=True)
    print(f"  full W histogram (fine shape): " + ", ".join(f"W={k}:{wh[k]}" for k in wh), flush=True)
    return tt, wm, wh

def main():
    print("# PROBE D2-b RIDER -- W top-scale trit histogram along the stationary flow. Exact rationals.")
    for L in [2, 3]:
        run(L)

if __name__ == "__main__":
    main()
