"""
PROBE D2-a -- (A) the 1/2-flux gate + (B) injection tables. CPU, exact rationals. Judge HELD (C).
A: on the mean-field 4-cell chain's stationary Perron flow, surviving flux splits EXACTLY 1/2 from
   unit-carry sources (v3=0) and 1/2 from divisible (v3>=1). population(2/3,1/3) x survival -> equality.
B: injection content W := That + c = gamma' - floor(gamma/3), tabulated ALONG the stationary flow:
   W mod 3 / mod 9 by (source v3-class, e' mod 6); (W mod 3) x (dest v'-class). Raw material for factor two.
Cells (G1 order): 0=(ev,v0) 1=(ev,v>=1) 2=(od,v0) 3=(od,v>=1); v3=0 cells {0,2}, v3>=1 cells {1,3}.
"""
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
import sympy
from probe_phase2b_G1 import build_cell_transfer
from probe_phase2c0 import build_M_tower_and_coords

def v3(n):
    if n == 0: return 99
    k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def cell_of(e, gam):
    return (e % 2) * 2 + (1 if gam % 3 == 0 else 0)

def perron_left(T):
    """exact LEFT Perron eigenvector (pi T = (1/3) pi) of the 4x4 = QSD occupation, normalized to sum 1.
    (build_cell_transfer is T[src][dst]: row-sums are survivals, so occupation is the left eigenvector.)"""
    M = sympy.Matrix([[sympy.Rational(T[i][j].numerator, T[i][j].denominator) for j in range(4)] for i in range(4)])
    ns = (M.T - sympy.Rational(1, 3) * sympy.eye(4)).nullspace()
    v = ns[0]; v = v / sum(v)
    return [Fr(int(x.p), int(x.q)) for x in v]

def fluxgate(L):
    T, casc, cellcount, spread, D, qL = build_cell_transfer(L)
    v = perron_left(T)                                   # cell occupation (QSD, left eigvec)
    surv = [sum(T[i][j] for j in range(4)) for i in range(4)]    # ROW sums = survival (T[src][dst])
    flux = [v[j] * surv[j] for j in range(4)]
    total = sum(flux)
    f_v0 = flux[0] + flux[2]                              # v3=0 sources
    f_v1 = flux[1] + flux[3]                              # v3>=1 sources
    print(f"  L={L}: cell occupation v={[str(x) for x in v]}  survival={[str(x) for x in surv]}", flush=True)
    print(f"        total surviving flux={total} (=1/3?); v3=0 flux={f_v0} share={f_v0/total}; "
          f"v3>=1 flux={f_v1} share={f_v1/total}", flush=True)
    ok = (f_v0 / total == Fr(1, 2) and f_v1 / total == Fr(1, 2))
    print(f"        => 1/2-FLUX GATE: {'PASS (exactly 1/2 each)' if ok else 'DEV'}", flush=True)
    return v, cellcount, ok

def injection_tables(L, v, cellcount):
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L)
    wf = [Fr(2 ** (D - d), 2 ** D - 1) for d in range(1, D + 1)]     # exact normalized halving weights
    # per-state stationary occupation = v[cell]/count[cell]
    occ = {c: v[c] / cellcount[c] for c in range(4)}
    # accumulators (all flux-weighted, exact Fraction)
    Wm3 = defaultdict(lambda: defaultdict(Fr))   # (v3class, e'%6) -> {W%3: flux}
    Wm9 = defaultdict(lambda: defaultdict(Fr))
    Wm3_dst = defaultdict(lambda: defaultdict(Fr))  # (W%3) -> {dst v'class: flux}
    for src_t, (a, e, gam) in enumerate(twcoords):
        dla = dl[a]; sc = cell_of(e, gam); vs = 0 if gam % 3 != 0 else 1
        os = occ[sc]
        for u in sub:
            dlu = dl[u]; da = (dla - dlu) % D; da = da if da != 0 else D
            for s in range(D):
                ep = (e + s) % D; bp = (u * two[ep]) % qL; T = (u - bp) % qL
                if (gam + T) % q != 0: continue
                gp = ((gam + T) // q) % qL
                dst = (u, bp, gp)
                if dst not in idx or idx[dst] not in pos: continue
                db = (da - s) % D; db = db if db != 0 else D
                flux = os * wf[da - 1] * wf[db - 1]
                W = gp - (gam // 3)                       # injected content = That + c
                vd = 0 if gp % 3 != 0 else 1
                Wm3[(vs, ep % 6)][W % 3] += flux
                Wm9[(vs, ep % 6)][W % 9] += flux
                Wm3_dst[W % 3][vd] += flux
    return Wm3, Wm9, Wm3_dst

def norm_dist(d):
    tot = sum(d.values())
    return {k: (d[k] / tot) for k in sorted(d)} if tot else {}

def main():
    print("# PROBE D2-a -- (A) 1/2-flux gate + (B) injection tables. Exact rationals. Judge HELD.")
    print("\n## A) THE 1/2-FLUX GATE")
    for L in [2, 3]:
        fluxgate(L)
    print("\n## B) INJECTION TABLES (W = That + c = gamma' - floor(gamma/3); flux-weighted along stationary flow)")
    L = 3
    v, cellcount, _ = fluxgate(L)
    Wm3, Wm9, Wm3_dst = injection_tables(L, v, cellcount)
    lines = [f"# D2-a injection tables q=3 L={L}. W=That+c=gamma'-floor(gamma/3). Flux-weighted (mean-field stationary).",
             "# === W mod 3 by (source v3-class, e' mod 6) ==="]
    print("\n  W mod 3 by (source v3-class, e' mod 6):", flush=True)
    for vs in (0, 1):
        for ep6 in range(6):
            dd = norm_dist(Wm3[(vs, ep6)])
            if dd:
                s = ", ".join(f"W%3={k}:{str(dd[k])}" for k in dd)
                tag = f"v3={'0' if vs==0 else '>=1'} e'%6={ep6}"
                print(f"    {tag}: {s}", flush=True)
                lines.append(f"{tag}\t" + "\t".join(f"{k}:{dd[k]}" for k in dd))
    lines.append("# === W mod 9 by (source v3-class, e' mod 6) ===")
    for vs in (0, 1):
        for ep6 in range(6):
            dd = norm_dist(Wm9[(vs, ep6)])
            if dd:
                lines.append(f"v3={'0' if vs==0 else '>=1'} e'%6={ep6}\t" + "\t".join(f"{k}:{dd[k]}" for k in dd))
    lines.append("# === (W mod 3) x (dest v'-class) joint (survival-conditioned injection) ===")
    print("\n  (W mod 3) x (dest v'-class) joint:", flush=True)
    for wm in sorted(Wm3_dst):
        dd = norm_dist(Wm3_dst[wm])
        s = ", ".join(f"v'={'0' if k==0 else '>=1'}:{str(dd[k])}" for k in dd)
        print(f"    W%3={wm}: {s}", flush=True)
        lines.append(f"W%3={wm}\t" + "\t".join(f"v'{k}:{dd[k]}" for k in dd))
    with open("outputs/injection_tables_q3.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n  dumped outputs/injection_tables_q3.tsv", flush=True)

if __name__ == "__main__":
    main()
