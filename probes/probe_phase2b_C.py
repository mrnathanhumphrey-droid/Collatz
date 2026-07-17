"""
PROBE C -- coupled carry-phase compression: build the hand-derivation's target object.
Recon/validation only; no proof, no rate fit. INSTRUMENT LAW: dense/direct at q=3.

M = build_M_gen(3, L, 2, [lam^d]), lam=1/2. e_rho = dlog2(b a^-1) in Z/D, D=2*3^(L-1).
Compression subspace S = {functions of (theta = reduce(e_rho), gamma)}: Galerkin-lump M^T
onto class indicators -> Lmat[src_class, dst_class] = (1/|src|) sum_{src,dst} M[dst,src].
eig(Lmat) = spectrum of left-eigenvectors that live in S.

C1: does the PARTNER (0.346827@L2, 0.333236@L3) survive? ladder of reductions:
    e mod 3 (theta, dim 3q^L) -> e mod 6 -> e mod 9 -> e mod D (full (e,gamma)).
    survive = min|eig(Lmat)-partner|/partner <= 1e-2. dump the coarsest working Lmat (rationals@L2).
C2: per (theta,gamma) class, v3(gamma') law P(=j)=2*3^-(j+1) (tail piled at top); per-class dev + survival.
C3: lumpability -- within-class spread of class-aggregated outgoing weights.
C4: q=7 L=2 control, theta = e mod d (d=ord_7(2)=3); does rank-22 partner (0.158414) survive?
"""
from fractions import Fraction
import numpy as np
import scipy.sparse as sp

from probe_phase2a_q2b_q6 import build_M_gen, subgroup


def real_M(q, L, lam=0.5):
    qL = q ** L
    D = len(subgroup(2 % qL, qL))
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    dl = {}; x = 1 % qL
    for e in range(D):
        dl[x] = e; x = (x * 2) % qL
    erho = np.array([dl[(b * pow(a, -1, qL)) % qL] for (a, b, g) in states])
    return M, states, n, D, qL, erho


def vq(g, q, L):
    if g == 0:
        return L
    v = 0
    while g % q == 0:
        g //= q; v += 1
    return v


def lumped(M, states, cls_of):
    classes = sorted(set(cls_of))
    ci = {c: j for j, c in enumerate(classes)}
    m = len(classes)
    size = np.zeros(m)
    for c in cls_of:
        size[ci[c]] += 1
    idxc = np.array([ci[c] for c in cls_of])
    Lm = np.zeros((m, m))
    Mc = M.tocoo()
    np.add.at(Lm, (idxc[Mc.col], idxc[Mc.row]), Mc.data)   # L[src_class, dst_class] += w
    Lm /= size[:, None]
    return classes, Lm, size


def circ0(D, lam=0.5):
    raw = np.array([lam ** d for d in range(1, D + 1)]); w = raw / raw.sum()
    return float(np.sum(w ** 2))                        # c0 = sum w^2 (real)


def survive(Lm, c0, partner, tol=1e-2):
    """Partner survives ONLY if the compressed spectrum has an eigenvalue near the partner that is
    DISTINCT from c0's image (c0 is kinematic, always in S; it must not masquerade as the partner)."""
    ev = np.linalg.eigvals(Lm)
    i0 = int(np.argmin(np.abs(ev - c0)))               # c0's image
    cand = [i for i in range(len(ev)) if i != i0]
    ip = min(cand, key=lambda i: abs(ev[i] - partner)) # nearest to partner, EXCLUDING c0's image
    d0 = abs(ev[i0] - c0); dp = abs(ev[ip] - partner)
    distinct = abs(ev[ip] - ev[i0]) > 0.5 * abs(c0 - partner)
    ok = (dp / abs(partner) <= tol) and distinct and (d0 < 1e-6)
    top = sorted(ev, key=lambda z: -z.real)[:4]
    return dict(ev0=ev[i0], evp=ev[ip], d0=d0, dp=dp, distinct=distinct, ok=ok, top=top)


def C1(q, L, partner, dump=True):
    print(f"\n## C1 q={q} L={L}  partner={partner}")
    M, states, n, D, qL, erho = real_M(q, L)
    gam = np.array([s[2] for s in states])
    c0 = circ0(D)
    ladder = [("e mod 3", 3), ("e mod 6", 6), ("e mod 9", 9), ("e mod D(full)", D)]
    worked = None
    for name, mod in ladder:
        if mod > D:
            continue
        cls_of = list(zip((erho % mod).tolist(), gam.tolist()))
        classes, Lm, size = lumped(M, states, cls_of)
        r = survive(Lm, c0, partner)
        print(f"   S=({name}, gamma): dim={len(classes):<4} c0_img@{r['ev0'].real:.6f}(d{r['d0']:.0e}) "
              f"partner_img@{r['evp'].real:.6f}{r['evp'].imag:+.4f}j (d{r['dp']:.1e}, rel {r['dp']/abs(partner):.0e}) "
              f"distinct={r['distinct']}  [{'SURVIVES(distinct)' if r['ok'] else 'no'}]")
        if r['ok'] and worked is None:
            worked = (name, mod, classes, Lm, size)
    if worked is None:
        print("   *** partner does NOT survive the (theta=e mod 3,gamma) reduction NOR full (e,gamma) at this L ***")
    else:
        print(f"   -> coarsest working subspace: ({worked[0]}, gamma), dim {len(worked[2])}")
    # ALWAYS dump the FULL gauge-invariant (e_rho,gamma) compression as the best-available handoff
    if dump:
        cls_full = cls_of_for(erho, gam, D)
        classes, Lm, size = lumped(M, states, cls_full)
        r = survive(Lm, c0, partner)
        path = f"outputs/compressed_q{q}_L{L}.tsv"
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"# compressed Lmat[src,dst], S=(e_rho,gamma) FULL gauge-invariant, q={q} L={L}, dim={len(classes)}\n")
            f.write(f"# APPROXIMATE handoff: partner captured to rel {r['dp']/abs(partner):.1e} (Galerkin; NOT exact -- see C3 lumpability)\n")
            f.write("# row/col label = (e_rho,gamma); (theta=e mod 3 reduction REFUTED -- needs full e_rho)\n")
            f.write("src\\dst\t" + "\t".join(f"({c[0]},{c[1]})" for c in classes) + "\n")
            for i, c in enumerate(classes):
                f.write(f"({c[0]},{c[1]})\t" + "\t".join(f"{Lm[i, j]:.10g}" for j in range(len(classes))) + "\n")
        print(f"   dumped FULL (e_rho,gamma) -> {path}  (rel {r['dp']/abs(partner):.1e})")
        if q == 3 and L == 2:
            dump_rational(M, states, cls_full, classes, size, "outputs/compressed_q3_L2_exact.tsv", "e_rho,gamma FULL")
    return worked


def cls_of_for(erho, gam, mod):
    return list(zip((erho % mod).tolist(), gam.tolist()))


def dump_rational(M, states, cls_of, classes, size, path, name):
    # rebuild Lmat exactly with Fraction weights (L=2 small)
    ci = {c: j for j, c in enumerate(classes)}
    m = len(classes)
    acc = [[Fraction(0) for _ in range(m)] for _ in range(m)]
    Mc = M.tocoo()
    # recover exact weights: entry = wa*wb; but M.data is float -> reconstruct from nearest rational denominator
    # weights are k/((2^D-1)^2 * ...) ; simplest exact route: recompute Lmat by re-summing rational weights per edge
    # Here: limit_denominator recovers the exact rational (denominators are small: products of 2^-d/(1-2^-D))
    for r, c, v in zip(Mc.row, Mc.col, Mc.data):
        acc[ci[cls_of[c]]][ci[cls_of[r]]] += Fraction(v).limit_denominator(10 ** 8)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# EXACT rational compressed Lmat[src,dst], S=({name},gamma), q=3 L=2\n")
        f.write("src\\dst\t" + "\t".join(f"({c[0]},{c[1]})" for c in classes) + "\n")
        for i, c in enumerate(classes):
            row = [acc[i][j] / Fraction(int(size[i])) for j in range(m)]
            f.write(f"({c[0]},{c[1]})\t" + "\t".join(str(x) for x in row) + "\n")
    print(f"   exact-rational dump -> {path}")


def C2(q, L):
    print(f"\n## C2 cascade law per (theta,gamma) class, q={q} L={L}")
    M, states, n, D, qL, erho = real_M(q, L)
    gam = np.array([s[2] for s in states])
    cls_of = list(zip((erho % 3).tolist(), gam.tolist()))
    Mc = M.tocoo()
    # per source class: weighted distribution of v_q(gamma')
    from collections import defaultdict
    dist = defaultdict(lambda: np.zeros(L + 1)); tot = defaultdict(float)
    for r, c, v in zip(Mc.row, Mc.col, Mc.data):
        sc = cls_of[c]; vd = vq(states[r][2], q, L)
        dist[sc][vd] += v; tot[sc] += v
    law = np.array([2 * q ** (-(j + 1)) for j in range(L)] + [0.0])
    law[L] = 1 - law[:L].sum()      # tail piled at top (gamma'=0)
    devs = []
    for sc in sorted(dist):
        p = dist[sc] / tot[sc]
        devs.append(np.max(np.abs(p - law)))
    devs = np.array(devs)
    print(f"   cascade law P(v_q(gamma')=j) = 2*{q}^-(j+1), tail at top = {np.round(law,4)}")
    print(f"   per-class max deviation: min={devs.min():.2e} med={np.median(devs):.2e} max={devs.max():.2e} "
          f"over {len(devs)} classes")
    print(f"   -> {'CLASS-INDEPENDENT cascade (coupling is in the SURVIVAL RATES, see below)' if devs.max()<1e-9 else 'class-DEPENDENT cascade (coupling partly here)'}")
    # survival rate per class (total surviving weight out of full 1 per state avg)
    surv = defaultdict(float); cnt = defaultdict(int)
    for i, sc in enumerate(cls_of):
        cnt[sc] += 1
    for r, c, v in zip(Mc.row, Mc.col, Mc.data):
        surv[cls_of[c]] += v
    sr = np.array([surv[sc] / cnt[sc] for sc in sorted(cnt)])
    print(f"   survival weight per class (avg out-weight): min={sr.min():.3f} max={sr.max():.3f} "
          f"spread={sr.max()-sr.min():.3f}  ({'uniform' if sr.max()-sr.min()<1e-9 else 'CLASS-DEPENDENT -> coupling lives here'})")


def C3(q, L):
    print(f"\n## C3 lumpability diagnostic on (theta,gamma), q={q} L={L}")
    M, states, n, D, qL, erho = real_M(q, L)
    gam = np.array([s[2] for s in states])
    cls_of = list(zip((erho % 3).tolist(), gam.tolist()))
    classes = sorted(set(cls_of)); ci = {c: j for j, c in enumerate(classes)}; m = len(classes)
    # per source STATE: outgoing weight vector to dest classes
    from collections import defaultdict
    rows = defaultdict(lambda: np.zeros(m))
    Mc = M.tocoo()
    for r, c, v in zip(Mc.row, Mc.col, Mc.data):
        rows[c][ci[cls_of[r]]] += v
    # within-class spread
    bycls = defaultdict(list)
    for i in range(n):
        bycls[cls_of[i]].append(rows[i])
    spreads = []
    for sc, vs in bycls.items():
        V = np.array(vs)
        spreads.append((sc, np.max(V.max(0) - V.min(0))))
    sp_vals = np.array([s for _, s in spreads])
    print(f"   within-class spread of outgoing weights: min={sp_vals.min():.2e} med={np.median(sp_vals):.2e} "
          f"max={sp_vals.max():.2e}  ({'LUMPABLE (0.258 obstruction absent here)' if sp_vals.max()<1e-9 else 'NOT globally lumpable (expected)'})")
    # partner-relevant classes = gamma!=0 (tower home)
    tower = [s for sc, s in spreads if sc[1] != 0]
    if tower:
        print(f"   on tower classes (gamma!=0, partner home): max spread={max(tower):.2e}")


def main():
    print("# PROBE C -- coupled carry-phase compression (build the derivation target; no proof, no fit)")
    print("# NOTE: Probe P's gauge test had a k^2 bug -> P2's 'Z/3 sub-family' was partly artifact;")
    print("#       C1 tests the compression DIRECTLY via a reduction ladder, settling it empirically.")
    C1(3, 2, 0.346827)
    C1(3, 3, 0.333236)
    C2(3, 2); C2(3, 3)
    C3(3, 2); C3(3, 3)
    C1(7, 2, 0.158414)   # C4 control


if __name__ == "__main__":
    main()
