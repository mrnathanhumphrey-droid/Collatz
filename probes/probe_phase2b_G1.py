"""
MICRO-PROBE G1 -- coarse tower transfer dump (the derivation's judge). Direct/exact at q=3.

Cells over the CLOSED carry tower (gamma!=0; no return to gamma=0 since gamma+T>=q => gamma'>=1):
   cell = (e_rho mod 2, v3(gamma) in {0, >=1}).  4 cells.
Deliverable (RAW transfer only -- no derived columns, no eigenvalues in the dump):
  (1) 4x4 cell-to-cell transfer, SOURCE-SIDE UNIFORM within cells, EXACT rationals, L=2 and L=3.
      T[c,c'] = (1/|c states|) * sum_{src in c} sum_{dst in c'} M[dst,src].
  (2) per-source-cell CASCADE split: destination v3(gamma') distribution {0, 1, >=2}, exact.
Judge (reported separately, NOT in the raw dump): row-sums vs {2/9,5/18,5/9,4/9} (pre-reg, gate-matched);
  Perron(4x4) and its L-trend vs rho(M_tower) = 0.346827 / 0.333236; within-cell spread = lumpability.

Exact arithmetic: w_delta = 2^{D-delta}/(2^D-1); w_da*w_db = 2^{2D-da-db}/(2^D-1)^2. Accumulate integer
numerators over common denom (2^D-1)^2; divide by source-cell state-count. e'_rho = (e_rho+da-db) mod D
(no modinv needed). Source-uniform => matches the compressed tower block lumped to cells.
"""
import numpy as np
from fractions import Fraction
from collections import defaultdict

from probe_phase2a_q2b_q6 import subgroup

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))

CELLS = [(0, 0), (0, 1), (1, 0), (1, 1)]        # (e parity, v3class: 0 => v3=0, 1 => v3>=1)
CLAB = {(0, 0): "even,v3=0", (0, 1): "even,v3>=1", (1, 0): "odd,v3=0", (1, 1): "odd,v3>=1"}
CI = {c: i for i, c in enumerate(CELLS)}


def v3(g):
    if g == 0: return 99
    v = 0
    while g % 3 == 0: g //= 3; v += 1
    return v


def build_cell_transfer(L, lam=0.5):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    inv2 = pow(2, -1, qL)
    pw = [pow(inv2, d, qL) for d in range(D + 1)]        # 2^{-d} mod qL
    # exact integer weight numerators wn[da][db] = 2^{2D-da-db}; common denom (2^D-1)^2
    twoD = 2 ** D
    wn = [[0] * (D + 1) for _ in range(D + 1)]
    wf = [[0.0] * (D + 1) for _ in range(D + 1)]
    Z = twoD - 1
    for da in range(1, D + 1):
        for db in range(1, D + 1):
            wn[da][db] = 2 ** (2 * D - da - db)
            wf[da][db] = (2 ** (D - da) / Z) * (2 ** (D - db) / Z)
    tcl = [0] * qL                                        # v3-class of gamma (0 => v3=0, 1 => v3>=1)
    v3c = [0] * qL                                        # min(v3,2) for cascade
    for g in range(1, qL):
        vv = v3(g); tcl[g] = 0 if vv == 0 else 1; v3c[g] = min(vv, 2)

    num = [[0] * 4 for _ in range(4)]                     # exact integer numerators, common denom
    cascade = [defaultdict(int) for _ in range(4)]        # per src cell: dest v3(0,1,>=2) -> numerator
    cellcount = [0] * 4
    # lumpability: per state dest-cell flow (float); track min/max per (src_cell, dst_cell)
    minv = np.full((4, 4), np.inf); maxv = np.full((4, 4), -np.inf)

    for a in sub:
        ap_all = [(a * pw[d]) % qL for d in range(D + 1)]
        for b in sub:
            er = dl[(b * (a ** 0) * pow(a, -1, qL)) % qL]  # e_rho = dlog2(b a^-1)
            par = er % 2
            bp_all = [(b * pw[d]) % qL for d in range(D + 1)]
            for g in range(1, qL):                        # tower source (gamma != 0)
                sc = CI[(par, tcl[g])]
                cellcount[sc] += 1
                statevec = [0.0, 0.0, 0.0, 0.0]
                for da in range(1, D + 1):
                    ap = ap_all[da]
                    for db in range(1, D + 1):
                        bp = bp_all[db]
                        T = (ap - bp) % qL
                        if (g + T) % q == 0:
                            gp = ((g + T) // q) % qL      # always >=1 (tower closed)
                            dpar = (par + da - db) % 2
                            dc = CI[(dpar, tcl[gp])]
                            num[sc][dc] += wn[da][db]
                            cascade[sc][v3c[gp]] += wn[da][db]
                            statevec[dc] += wf[da][db]
                for dc in range(4):
                    if statevec[dc] < minv[sc][dc]: minv[sc][dc] = statevec[dc]
                    if statevec[dc] > maxv[sc][dc]: maxv[sc][dc] = statevec[dc]
    denom = Z * Z
    T = [[Fraction(num[i][j], denom * cellcount[i]) for j in range(4)] for i in range(4)]
    casc = [{v: Fraction(cascade[i][v], denom * cellcount[i]) for v in sorted(cascade[i])} for i in range(4)]
    spread = maxv - minv
    return T, casc, cellcount, spread, D, qL


def run(L, rho_tower):
    log(f"\n{'='*78}\n## G1  q=3 L={L}  tower cell transfer (source-uniform, exact)   rho(M_tower)={rho_tower}")
    T, casc, cellcount, spread, D, qL = build_cell_transfer(L)
    # (1) 4x4 transfer
    log("\n   (1) 4x4 CELL TRANSFER  T[src -> dst]  (exact rationals)")
    log("       " + "src \\ dst".ljust(12) + "".join(CLAB[c].rjust(14) for c in CELLS) + "   row-sum")
    prereg = {(0, 0): Fraction(2, 9), (0, 1): Fraction(5, 9), (1, 0): Fraction(5, 18), (1, 1): Fraction(4, 9)}
    rowsum_ok = True
    for i, c in enumerate(CELLS):
        rs = sum(T[i])
        if rs != prereg[c]: rowsum_ok = False
        log("       " + CLAB[c].ljust(12) + "".join(str(T[i][j]).rjust(14) for j in range(4))
            + f"   {rs} {'OK' if rs==prereg[c] else '!= '+str(prereg[c])}")
    log(f"   row-sums vs pre-reg {{2/9,5/9,5/18,4/9}}: {'ALL MATCH' if rowsum_ok else 'MISMATCH'}")
    # (2) cascade split
    log("\n   (2) CASCADE SPLIT per source cell: destination v3(gamma') flow {0, 1, >=2} (exact; sums to survival)")
    for i, c in enumerate(CELLS):
        parts = "  ".join(f"v'={v}:{casc[i][v]}={float(casc[i][v]):.5f}" for v in sorted(casc[i]))
        tot = sum(casc[i].values())
        frac = "  ".join(f"v'={v}:{float(casc[i][v]/tot):.4f}" for v in sorted(casc[i]))
        log(f"       {CLAB[c].ljust(12)} flow: {parts}   | as frac of survival: {frac}")
    # (3) lumpability + Perron (the judge; reported, not in raw dump)
    log("\n   (3) JUDGE (reported, not in the raw dump):")
    log(f"       within-cell spread of per-state dest-cell flow: max = {spread.max():.3e}  "
        f"({'EXACT lumping (states in a cell are transfer-identical)' if spread.max()<1e-12 else 'approximate (uniform-averaged)'})")
    Tf = np.array([[float(T[i][j]) for j in range(4)] for i in range(4)])
    ev = np.linalg.eigvals(Tf)
    per = ev[np.argmax(np.abs(ev))]
    log(f"       Perron(4x4) = {per.real:.8f}{per.imag:+.2e}j   rho(M_tower) = {rho_tower}   "
        f"|Perron - rho| = {abs(per.real - rho_tower):.3e}")
    log(f"       4x4 eigenvalues: " + "  ".join(f"{z.real:+.6f}{z.imag:+.4f}j" for z in sorted(ev, key=lambda z:-abs(z))))
    # dump
    lines = [f"# G1 tower cell transfer q=3 L={L}. Cells (e_rho parity, v3(gamma):0 vs >=1). Source-uniform, EXACT.",
             "# (1) 4x4 transfer T[src->dst]",
             "src\\dst\t" + "\t".join(CLAB[c] for c in CELLS)]
    for i, c in enumerate(CELLS):
        lines.append(CLAB[c] + "\t" + "\t".join(str(T[i][j]) for j in range(4)))
    lines.append("# (2) cascade split: source cell -> dest v3(gamma') flow {0,1,>=2} (exact)")
    lines.append("src_cell\t" + "\t".join(f"v'={v}" for v in [0, 1, 2]))
    for i, c in enumerate(CELLS):
        lines.append(CLAB[c] + "\t" + "\t".join(str(casc[i].get(v, Fraction(0))) for v in [0, 1, 2]))
    with open(f"outputs/tower_cell_transfer_q3_L{L}.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"   dumped outputs/tower_cell_transfer_q3_L{L}.tsv")
    return per.real, rowsum_ok, spread.max()


def main():
    log("# MICRO-PROBE G1 -- coarse tower transfer dump (raw; Wilson derives blind, judged entry-by-entry).")
    p2, ok2, s2 = run(2, 0.346827)
    p3, ok3, s3 = run(3, 0.333236)
    log(f"\n{'='*78}")
    log(f"## SUMMARY: row-sums {'MATCH both L' if (ok2 and ok3) else 'MISMATCH'}; "
        f"lumping {'EXACT both L' if (s2<1e-12 and s3<1e-12) else 'approximate'}.")
    log(f"   Perron L-trend: L=2 {p2:.6f} (rho 0.346827, d {abs(p2-0.346827):.1e}), "
        f"L=3 {p3:.6f} (rho 0.333236, d {abs(p3-0.333236):.1e}).")
    with open("logs/probe_phase2b_G1_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
