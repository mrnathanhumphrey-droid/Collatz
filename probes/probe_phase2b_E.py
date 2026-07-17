"""
PROBE E -- freeze the compressed (e_rho, gamma)-chain + L-series. Instrument only.
INSTRUMENT LAW: dense/direct at q=3. The compression is built DIRECTLY (never materialize
the full operator), so L=4 (dim 4374) is reachable where the 236k-state operator walls.

FROZEN CONVENTION (E1): Lmat[src_class, dst_class] = (1/|src_class|) * sum_{src in class}
sum_{dst in class} M[dst, src]  -- the SOURCE-SIDE UNIFORM average of outgoing weights
(each source state weighted 1/|class|; |class| = D). This is the Galerkin projection
P M^T P with P = orthogonal projection onto class-indicator functions under the UNIFORM
(counting) inner product on each class. The single free choice = UNIFORM source weighting
(not stationary-weighted, not right-eigenvector-weighted). Left-action; eigenvalues
approximate M-eigenvalues whose LEFT eigenvectors are (e_rho,gamma)-class functions.

Classes indexed idx = e_rho*qL + gamma, e_rho in Z/D, gamma in Z/q^L.
"""
from fractions import Fraction
from collections import defaultdict
import numpy as np

from probe_phase2a_q2b_q6 import subgroup


def setup(q, L, lam=0.5):
    qL = q ** L
    sub = subgroup(2 % qL, qL)
    D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D):
        dl[x] = e; x = (x * 2) % qL
    a_arr = np.array(sub, dtype=np.int64)              # the D representatives a
    e_a = np.array([dl[a] for a in sub])
    raw = np.array([lam ** d for d in range(1, D + 1)]); w = raw / raw.sum()
    pw = np.array([pow(pow(2, -1, qL), d, qL) for d in range(D + 1)], dtype=np.int64)  # 2^{-d} mod qL, d=0..D
    two = np.array([pow(2, j, qL) for j in range(D)], dtype=np.int64)                  # 2^j mod qL
    return qL, sub, D, a_arr, e_a, w, pw, two


def build_compressed(q, L, record=False):
    qL, sub, D, a_arr, e_a, w, pw, two = setup(q, L)
    m = D * qL
    Lmat = np.zeros((m, m))
    G = np.arange(qL)
    decomp = defaultdict(lambda: defaultdict(float)) if record else None
    for e_rho in range(D):
        b_arr = (a_arr * two[e_rho]) % qL
        base = e_rho * qL
        for da in range(1, D + 1):
            ap = (a_arr * pw[da]) % qL
            for db in range(1, D + 1):
                bp = (b_arr * pw[db]) % qL
                T = (ap - bp) % qL                          # over a (length D)
                wt = w[da - 1] * w[db - 1]
                ep = (e_rho + da - db) % D
                dbase = ep * qL
                S = G[:, None] + T[None, :]                 # qL x D
                mask = (S % q) == 0
                gp = (S // q) % qL
                if not mask.any():
                    continue
                rows = (base + G)[:, None] * np.ones((1, D), dtype=np.int64)
                cols = dbase + gp
                r = rows[mask]; c = cols[mask]
                np.add.at(Lmat, (r, c), wt)
                if record:
                    for rr, cc in zip(r.tolist(), c.tolist()):
                        decomp[(rr, cc)][(da, db)] += wt
    Lmat /= D                                              # source-side uniform average
    if record:
        for k in decomp:
            for kk in decomp[k]:
                decomp[k][kk] /= D
    return Lmat, D, qL, w, decomp


def circ_family(w, D):
    d = np.arange(1, D + 1)
    return np.array([np.sum(w ** 2 * np.exp(2j * np.pi * k * d / D)) for k in range(D)])


def analyze(Lmat, w, D, true_partner=None):
    ev = np.linalg.eigvals(Lmat)
    cks = circ_family(w, D)
    c0 = cks[0].real
    order = np.argsort(-ev.real)
    top = ev[order[:6]]
    # flag family membership
    def fam(z): return min(abs(z - cks)) < 1e-9
    # compressed Perron
    perron = ev[order[0]].real
    # compressed partner = top eigenvalue distinct from c0's image and not a family member
    i0 = int(np.argmin(np.abs(ev - c0)))
    cand = [i for i in range(len(ev)) if i != i0 and min(abs(ev[i] - cks)) > 1e-9]
    ip = max(cand, key=lambda i: ev[i].real) if cand else None
    partner_c = ev[ip] if ip is not None else None
    gap = abs(partner_c - ev[i0]) if partner_c is not None else None
    return dict(ev=ev, c0=c0, perron=perron, top=top, fam=[fam(z) for z in top],
                partner_c=partner_c, gap=gap,
                perr=(abs(partner_c - true_partner) / abs(true_partner)) if (true_partner and partner_c is not None) else None)


def emit_rational(Lmat, D, qL, path, label):
    m = D * qL
    lines = [f"# EXACT-rational compressed Lmat[src,dst], {label}, dim={m} (common-denom family)",
             "# label (e_rho,gamma); idx = e_rho*qL+gamma"]
    header = ["src\\dst"] + [f"({i // qL},{i % qL})" for i in range(m)]
    lines.append("\t".join(header))
    for i in range(m):
        row = [f"({i // qL},{i % qL})"]
        for j in range(m):
            if abs(Lmat[i, j]) < 1e-14:
                row.append("0")
            else:
                row.append(str(Fraction(Lmat[i, j]).limit_denominator(10 ** 9)))
        lines.append("\t".join(row))
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def emit_decomp(decomp, D, qL, w, path, label):
    lines = [f"# entry decomposition {label}: for each nonzero Lmat[src,dst], the (da,db) move-pairs and w-products",
             "# src(e_rho,gamma)\tdst(e_rho,gamma)\tentry\t(da,db):wprod ..."]
    for (r, c), dd in sorted(decomp.items()):
        entry = sum(dd.values())
        parts = " ".join(f"({da},{db}):{val:.6g}" for (da, db), val in sorted(dd.items()))
        lines.append(f"({r // qL},{r % qL})\t({c // qL},{c % qL})\t{entry:.8g}\t{parts}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main():
    print("# PROBE E -- compressed (e_rho,gamma)-chain L-series (frozen source-side uniform average)\n")

    # verify against Probe C's L=2 dump (65/189 self-loop etc.)
    Lm2, D2, qL2, w2, dec2 = build_compressed(3, 2, record=True)
    print(f"## VERIFY L=2: Lmat[(0,0)->(0,0)] = {Fraction(Lm2[0,0]).limit_denominator(10**7)} "
          f"(Probe C: 65/189={65/189:.6f}); c0=sum w^2={circ_family(w2,D2)[0].real:.6f}")

    print("\n## E2 -- q=3 L-series (a Perron->1/3, b top-6, c compressed gap, d error)")
    truep = {2: 0.346827, 3: 0.333236}
    for L in [2, 3, 4]:
        Lm, D, qL, w, _ = build_compressed(3, L)
        A = analyze(Lm, w, D, truep.get(L))
        print(f"   L={L} dim={D*qL}: Perron={A['perron']:.6f} (c0={A['c0']:.6f}, dist to 1/3 = {abs(A['perron']-1/3):.2e})")
        tops = " ".join(f"{z.real:+.5f}{'*' if fl else ''}" for z, fl in zip(A['top'], A['fam']))
        print(f"        top6 (*=family): {tops}")
        if A['partner_c'] is not None:
            print(f"        compressed-partner={A['partner_c'].real:+.6f}  compressed-gap(partner,c0)={A['gap']:.3e}"
                  + (f"  [true partner rel err {A['perr']:.1e}]" if A['perr'] else "  [L=4 true partner=G, not computed]"))
        if L in (2, 3):
            emit_rational(Lm, D, qL, f"outputs/compressed_q3_L{L}_frozen.tsv", f"q=3 L={L} frozen-uniform")

    print("\n## E3 -- entry decomposition (L=2,3) -> orbit lists")
    emit_decomp(dec2, D2, qL2, w2, "outputs/compressed_q3_L2_decomp.tsv", "q=3 L=2")   # correct L=2 qL
    _, D3, qL3, w3, dec3 = build_compressed(3, 3, record=True)
    emit_decomp(dec3, D3, qL3, w3, "outputs/compressed_q3_L3_decomp.tsv", "q=3 L=3")
    print("   dumped compressed_q3_L{2,3}_frozen.tsv + _decomp.tsv")

    print("\n## E4 -- q=7 control L-series (pre-reg: compressed gap stays OPEN, mirrors r7~0.38)")
    truep7 = {2: 0.158414}
    for L in [2, 3]:
        from probe_phase2b_E import setup as _setup
        _, _, Dchk, _, _, _, _, _ = _setup(7, L)
        dim = Dchk * 7 ** L
        if dim > 20000:                                    # dense eig walls (~20GB at dim 50421)
            print(f"   q=7 L={L}: compressed dim={dim} -> DENSE EIG WALLS (report the wall; L=2 stands).")
            continue
        Lm, D, qL, w, _ = build_compressed(7, L)
        A = analyze(Lm, w, D, truep7.get(L))
        print(f"   q=7 L={L} dim={D*qL}: Perron={A['perron']:.6f} (dist to 1/3={abs(A['perron']-1/3):.2e})  "
              f"compressed-partner={A['partner_c'].real if A['partner_c'] is not None else float('nan'):+.5f}  "
              f"gap(partner,c0)={A['gap']:.3e}"
              + (f"  [rel {A['perr']:.1e}]" if A['perr'] else ""))


if __name__ == "__main__":
    main()
