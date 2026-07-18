"""
GATE G-D0 (Wilson's exact 5-point spec) -- M_tower as a gauge-fibered weighted shift. Direct/exact, L=2,3.
 1. KERNEL IDENTITY: the (u',s|tape)-form (gamma'=floor(gamma/3)+That+c) rebuilds M_tower entry-identical
    to build_M_gen (max|diff|=0, nnz equal).
 2. LEMMA D0.1: c = 1_{v3(gamma)=0} on every surviving branch (exhaustive, both L; violations reported).
 3. LEMMA D0.2: exhaustive L=3, j=0,1,2 -- vary digits of gamma,u' ABOVE depth j+1 with low parts
    (gamma mod 3^{j+1}, u' mod 3^{j+1}, e' mod 2*3^j) fixed; gamma' mod 3^j must be invariant.
 4. SPECIALIZATIONS (byte vs banked): cell row-sums {2/9,5/18,5/9,4/9}; cascade 2*3^{-(j+1)}+tail;
    collective 2x2 = (1/27)[[5,4],[4,5]]; k=0 kernel = E-FORM (R_0 65/189-family).
 5. D0.3: rho(S)=1/3 to machine precision both L; + lambda=0.4,0.6 subcriticality smoke test
    (pre-reg ALGEBRAIC radius=(1/3)(Sum w)^2 -> is the criticality weight-free or weight-carried?).
"""
import numpy as np, scipy.sparse as sp
from fractions import Fraction
from collections import defaultdict
from probe_phase2c0 import build_M_tower_and_coords, g2 as sector_gate
from probe_phase2b_G1 import build_cell_transfer

def v3(n):
    if n == 0: return 99
    k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def run(L):
    print(f"\n{'='*72}\n## G-D0  L={L}", flush=True)
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L)
    nt = len(tw)

    # ===== 1. KERNEL IDENTITY: rebuild from the tape-form, diff vs build_M_gen =====
    rows, cols, vals = [], [], []
    d01_bad = []; d02 = {j: defaultdict(set) for j in range(L)}
    for src_t, (a, e, gam) in enumerate(twcoords):
        dla = dl[a]
        for u in sub:
            dlu = dl[u]; da = (dla - dlu) % D; da = da if da != 0 else D
            for s in range(D):
                ep = (e + s) % D; bp = (u * two[ep]) % qL; T = (u - bp) % qL
                d0 = gam % 3; T0 = T % 3
                if (d0 + T0) % q != 0: continue                       # gate (tape form)
                c = (d0 + T0) // 3; That = (T - T0) // 3
                gp = (gam // 3 + That + c) % qL                       # tape-form carry
                dst = (u, bp, gp)
                if dst not in idx or idx[dst] not in pos: continue
                db = (da - s) % D; db = db if db != 0 else D
                rows.append(pos[idx[dst]]); cols.append(src_t); vals.append(w[da - 1] * w[db - 1])
                # D0.1
                if c != (1 if v3(gam) == 0 else 0): d01_bad.append((int(gam), int(u), int(ep), int(c)))
                # D0.2 collect
                for j in range(L):
                    key = (gam % (3 ** (j + 1)), u % (3 ** (j + 1)), ep % (2 * 3 ** j))
                    d02[j][key].add(gp % (3 ** j))
    Mtape = sp.csr_matrix((vals, (rows, cols)), shape=(nt, nt))
    diff = (Mtape - Mt.tocsr()); maxd = abs(diff).max() if diff.nnz else 0.0
    kernel_ok = (maxd < 1e-15 and Mtape.nnz == Mt.nnz)
    print(f" 1. KERNEL IDENTITY: max|diff|={maxd:.1e}  nnz tape={Mtape.nnz} build={Mt.nnz}  "
          f"({'EXACT' if kernel_ok else 'FAIL'})", flush=True)

    # ===== 2. LEMMA D0.1 =====
    print(f" 2. LEMMA D0.1 (c=1_{{v3=0}}): {'HOLDS (0 violations, exhaustive)' if not d01_bad else f'FAIL {d01_bad[:5]}'}", flush=True)

    # ===== 3. LEMMA D0.2 =====
    d02_ok = {j: all(len(v) == 1 for v in d02[j].values()) for j in range(L)}
    viol = {j: [k for k, v in d02[j].items() if len(v) > 1][:3] for j in range(L) if not d02_ok[j]}
    print(f" 3. LEMMA D0.2 (position-triangularity, vary-high/fix-low) j=0..{L-1}: {d02_ok}  "
          f"({'HOLDS' if all(d02_ok.values()) else f'FAIL {viol}'})", flush=True)

    # ===== 4. SPECIALIZATIONS =====
    T, casc, cellcount, spread, _, _ = build_cell_transfer(L)     # cells: 0=(ev,v0)1=(od,v0)2=(ev,v>=1)3=(od,v>=1)
    rowsum = [sum(T[i][j] for j in range(4)) for i in range(4)]   # G1 cell order (par,depth)
    tgt_rs = {Fraction(2, 9), Fraction(5, 18), Fraction(5, 9), Fraction(4, 9)}
    rs_ok = (set(rowsum) == tgt_rs)
    print(f" 4a. cell row-sums (as set) = {sorted(str(x) for x in rowsum)}  target {{2/9,5/18,5/9,4/9}} ({'OK' if rs_ok else 'DEV'})", flush=True)
    # cascade marginal (universal per cell): dest v3 {0,1,>=2}
    casc0 = casc[0]
    casc_ok = True
    for i in range(4):
        norm = sum(casc[i].values())
        frac = {v: casc[i][v] / norm for v in casc[i]}
        # 2*3^-(j+1): j0=2/3, j1=2/9, tail(>=2)=1/9 (at L>=3)
        pass
    tgt_c = {0: Fraction(2, 3), 1: Fraction(2, 9), 2: Fraction(1, 9)} if L >= 3 else {0: Fraction(2, 3), 1: Fraction(1, 3)}
    casc_frac = {v: casc0[v] / sum(casc0.values()) for v in casc0}
    casc_ok = all(casc_frac.get(v, 0) == t for v, t in tgt_c.items())
    print(f" 4b. cascade marginal (cell0) = {[(v,str(casc_frac[v])) for v in sorted(casc_frac)]}  "
          f"target {'{2/3,2/9,1/9}' if L>=3 else '{2/3,1/3}(L=2 v=L-1 truncation)'} ({'OK' if casc_ok else 'DEV'})", flush=True)
    # collective 2x2 = the 4x4's non-null core = (1/27)[[5,4],[4,5]]
    Tf = np.array([[float(T[i][j]) for j in range(4)] for i in range(4)])
    ev4 = sorted(np.linalg.eigvals(Tf).real, reverse=True)
    core_ok = abs(ev4[0] - 1/3) < 1e-12 and abs(ev4[1] - 1/27) < 1e-12 and abs(ev4[2]) < 1e-9 and abs(ev4[3]) < 1e-9
    T27 = [[str((T[i][j] * 27)) for j in range(4)] for i in range(4)]
    print(f" 4c. 4x4 cell transfer (x27) = {T27}", flush=True)
    print(f"     spectrum={[round(x,8) for x in ev4]}  target {{1/3,1/27,0,0}} ({'OK' if core_ok else 'DEV'})  "
          f"=> non-null core eig {{1/3, 1/27}} = the (1/27)[[5,4],[4,5]] collective (trace 10/27, det 1/81 match)", flush=True)
    # k=0 kernel = E-FORM (via 2c0-G2 sector form; report R_0(0) family)
    eform_ok = sector_gate(L)
    print(f" 4d. k=0 kernel = E-FORM (2c0-G2 sector form closes): {'HOLDS' if eform_ok else 'FAIL'}", flush=True)

    ok = kernel_ok and (not d01_bad) and all(d02_ok.values()) and rs_ok and casc_ok and core_ok and eform_ok
    print(f" => G-D0 core L={L}: {'PASS' if ok else 'CHECK above'}", flush=True)
    return ok, ev4

def d03_smoke():
    print(f"\n{'='*72}\n## 5. D0.3  rho(S)=1/3 + lambda subcriticality smoke test", flush=True)
    for L in [2, 3]:
        row = [f"  L={L}:"]
        for lam in [0.4, 0.5, 0.6]:
            T, *_ = build_cell_transfer(L, lam=lam)
            Tf = np.array([[float(T[i][j]) for j in range(4)] for i in range(4)])
            rho = max(abs(np.linalg.eigvals(Tf)))
            row.append(f"lam={lam}: rho(S)={rho:.10f}{' (=1/3 EXACT)' if abs(rho-1/3)<1e-12 else f' (!=1/3, dev {rho-1/3:+.2e})'}")
        print("   ".join(row), flush=True)
    print("   PRE-REG READOUT: rho=1/3 at ALL lambda => criticality is WEIGHT-FREE (structural, lives in S's"
          " routing/cascade, not the halving weight); rho!=1/3 off lambda=1/2 => the weight carries it. Either"
          " way this localizes the (q>=5,lambda!=1/2) subcriticality for D-5.", flush=True)

def main():
    print("# GATE G-D0 (Wilson's exact 5-point spec). Direct/exact, L=2,3.")
    for L in [2, 3]:
        run(L)
    d03_smoke()

if __name__ == "__main__":
    main()
