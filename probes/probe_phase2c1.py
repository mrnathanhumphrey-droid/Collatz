"""
PROBE 2c1 -- zeroth bracket + parity-sector data + scale-ledger. Direct/entry algebra, no eigensolves.
INSTRUMENT LAW: direct at q=3 (and q=7 control). Claude gates; label ALGEBRAIC vs STATISTICAL.

A) Collatz-Wielandt baseline: h0 = cell-lift (2/3,5/3,5/6,4/3). Ratio (M^T h0)/h0 pointwise (M^T = flow
   convention, h0 is the coarse RIGHT eigvec of T[src,dst]=flow). Bracket [min,max] contains rho_L;
   histogram by cell. PRE-REG: 1/3 inside both L. width = defect scale.
B) Parity sector (k=D/2, chi=+-1): exact R_{D/2}(s) (s=0..D-1) L=2,3 (Fraction); N_{D/2}(e',g,g') L=2 exact int.
C) Scale ledger: ||D_k||^2 bucketed by a=v_q(k), q=3 (L=2,3) + q=7 L=2 control. PRE-REG: q=3 flat per scale
   (O(1) marginal); q=7 decaying (summable, gapped). Flat-vs-decaying = q3-vs-q7 marginality discriminator.
"""
import numpy as np
import scipy.sparse as sp
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
    ea = np.array([dl[s[0]] for s in states])
    erho = np.array([dl[(s[1] * pow(s[0], -1, qL)) % qL] for s in states])
    tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsc()
    return Mt, gam[tw], ea[tw], erho[tw], D, qL, np.array(raw) / sum(raw)


def cell_of(e, g):
    return (int(e % 2), 0 if g % 3 != 0 else 1)


# ---------------- A : bracket + histogram ----------------
def partA(L, partner):
    Mt, gam, ea, erho, D, qL, w = build_tower(3, L)
    nt = Mt.shape[0]
    par = erho % 2; v0 = (gam % 3 != 0)
    h0 = np.empty(nt)
    h0[(par == 0) & v0] = 2 / 3; h0[(par == 0) & ~v0] = 5 / 3
    h0[(par == 1) & v0] = 5 / 6; h0[(par == 1) & ~v0] = 4 / 3
    r = (Mt.T.dot(h0)) / h0                          # (M^T h0)/h0, flow convention
    m0, M0 = r.min(), r.max()
    log(f"\n   A) ZEROTH BRACKET L={L}: [min,max] = [{m0:.8f}, {M0:.8f}]  width={M0-m0:.6e}  "
        f"1/3 inside: {m0 <= 1/3 <= M0}  rho_L={partner} inside: {m0 <= partner <= M0}")
    log(f"      ratio histogram by cell:")
    for c in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        mask = np.array([cell_of(erho[i], gam[i]) == c for i in range(nt)])
        rr = r[mask]
        # distinct values (exact structure)
        uniq = sorted(set(np.round(rr, 8)))
        log(f"        cell {c}: n={mask.sum():5d}  [{rr.min():.6f}, {rr.max():.6f}]  mean {rr.mean():.6f}  "
            f"distinct ratios: {uniq if len(uniq) <= 6 else str(uniq[:6])+'...('+str(len(uniq))+')'}")
    return m0, M0


# ---------------- B : parity-sector tables ----------------
def partB(L):
    qL = 3 ** L; D = len(subgroup(2 % qL, qL)); k = D // 2
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    Z = 2 ** D - 1
    wf = [Fraction(2 ** (D - m), Z) for m in range(1, D + 1)]   # w_m exact
    # R_{D/2}(s) = sum_m w_m w_{m-s} (-1)^m   (chi_{D/2}(2^m)=(-1)^m)
    R = []
    for s in range(D):
        val = Fraction(0)
        for m in range(1, D + 1):
            mm = ((m - s - 1) % D) + 1
            val += wf[m - 1] * wf[mm - 1] * (-1) ** m
        R.append(val)
    log(f"\n   B) PARITY SECTOR k=D/2={k} (chi=(-1)^dlog), L={L}, D={D}:")
    log(f"      R_{{{k}}}(s), s=0..{D-1} (exact): " + "  ".join(f"{s}:{R[s]}" for s in range(min(D, 8)))
        + (" ..." if D > 8 else ""))
    if D > 8:
        log(f"        (full): " + "  ".join(f"{R[s]}" for s in range(D)))
    with open(f"outputs/parity_R_q3_L{L}.tsv", "w", encoding="utf-8") as f:
        f.write(f"# R_{{D/2={k}}}(s) exact, q=3 L={L}, D={D}. R_k(s)=sum_m w_m w_{{m-s}} (-1)^m\n# s\tR\n")
        for s in range(D): f.write(f"{s}\t{R[s]}\n")
    return R, k, D, qL, dl


def partB_N(L=2):
    """N_{D/2}(e',gamma,gamma') exact signed integer, L=2 full table."""
    qL = 3 ** L; sub = subgroup(2 % qL, qL); D = len(sub); k = D // 2
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    two = [pow(2, e, qL) for e in range(D)]
    N = defaultdict(int)
    for ep in range(D):
        c = (1 - two[ep]) % qL
        for u in sub:
            T = (u * c) % qL
            sign = (-1) ** dl[u]
            for gam in range(1, qL):
                if (gam + T) % 3 == 0:
                    gp = ((gam + T) // 3) % qL
                    N[(ep, gam, gp)] += sign
    nz = {kk: v for kk, v in N.items() if v != 0}
    log(f"      N_{{{k}}}(e',gamma,gamma') L={L}: {len(nz)} nonzero of {len(N)} entries; value range "
        f"[{min(nz.values())},{max(nz.values())}]; distinct values {sorted(set(nz.values()))}")
    # dump
    lines = [f"# N_{{D/2={k}}}(e',gamma,gamma') exact signed unit-count, q=3 L={L} (chi=(-1)^dlog u)",
             "# e'\tgamma\tgamma'\tN"]
    for (ep, gam, gp), v in sorted(nz.items()):
        lines.append(f"{ep}\t{gam}\t{gp}\t{v}")
    with open(f"outputs/parity_N_q3_L{L}.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return nz, k


# ---------------- C : scale ledger ----------------
def defect_amp(q, L):
    Mt, gam, ea, erho, D, qL, w = build_tower(q, L)
    nt = Mt.shape[0]
    orbits = defaultdict(dict)
    for i in range(nt): orbits[(int(erho[i]), int(gam[i]))][int(ea[i])] = i
    Vmat = np.exp(2j * np.pi * np.outer(np.arange(D), np.arange(D)) / D)
    amp = np.zeros(D)
    for key, dmap in orbits.items():
        cols = [dmap[m] for m in range(D)]
        S = Mt[:, cols].toarray()
        B = S @ Vmat
        amp += (np.abs(B) ** 2).sum(axis=0) / D
    return amp, D


def vq(k, q):
    if k == 0: return -1
    v = 0
    while k % q == 0: k //= q; v += 1
    return v


def partC(q, L):
    amp, D = defect_amp(q, L)
    buckets = defaultdict(list)
    for k in range(1, D):
        buckets[vq(k, q)].append(amp[k])
    log(f"\n   C) SCALE LEDGER q={q} L={L} (D={D}): defect mass ||D_k||^2 bucketed by a=v_{q}(k)")
    tot = sum(amp[1:])
    scales = sorted(buckets)
    for a in scales:
        arr = np.array(buckets[a])
        log(f"        a={a}: {len(arr):2d} sectors  total {arr.sum():9.5f} ({arr.sum()/tot*100:5.1f}%)  "
            f"mean/sector {arr.mean():.5f}")
    means = [np.mean(buckets[a]) for a in scales]
    ratios = [means[i+1]/means[i] for i in range(len(means)-1)]
    log(f"        mean/sector by scale: {[f'{m:.4f}' for m in means]}  ratios(a->a+1): {[f'{r:.3f}' for r in ratios]}")
    return {"means": means, "ratios": ratios, "share_deep": 1 - amp[[k for k in range(1,D) if vq(k,q)==0]].sum()/tot}


def main():
    log("# PROBE 2c1 -- zeroth bracket + parity-sector tables + scale ledger. Direct, no eigensolves.")
    log("\n## A) COLLATZ-WIELANDT ZEROTH BRACKET")
    partA(2, 0.346827); partA(3, 0.333236)
    log("\n## B) PARITY-SECTOR (k=D/2) CLOSED-FORM TABLES")
    partB(2); partB(3)
    partB_N(2)
    log("\n## C) SCALE LEDGER (marginality discriminator: q=3 flat vs q=7 decaying)")
    log("   -- q=3 (target):")
    c32 = partC(3, 2); c33 = partC(3, 3)
    log("   -- q=7 (CONTROL, pre-reg: NOT flat / decaying):")
    c72 = partC(7, 2)
    log("\n   VERDICT (honest, no fitting): compare count-normalized mean/sector ratios across scales:")
    log(f"      q=3 L=3 ratios {[f'{r:.3f}' for r in c33['ratios']]}  vs  q=7 L=2 ratio {[f'{r:.3f}' for r in c72['ratios']]}")
    log(f"      q=3 deep-scale mass share (a>=1) = {c33['share_deep']*100:.1f}% (L=3), {c32['share_deep']*100:.1f}% (L=2); "
        f"q=7 = {c72['share_deep']*100:.1f}% (L=2)")
    log("      => raw defect-mass ledger does NOT cleanly separate: q=3 mean/sector ratios (0.93 shallow, 0.64 deep)")
    log("         OVERLAP q=7 (0.72). q=3 spreads more mass into deep scales (26% vs 7%) BUT that is largely the")
    log("         3-adic sector COUNT (12/4/1) vs 7-adic (18/2), not per-sector intensity. The cheap raw-mass test")
    log("         is INCONCLUSIVE at this resolution; marginality (if real) lives in the tax-weighted / selection-")
    log("         graded combination (2c-3 pen), not the bare ||D_k||^2. NOT a clean q3-vs-q7 discriminator here.")
    with open("logs/probe_phase2c1_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
