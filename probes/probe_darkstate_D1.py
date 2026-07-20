"""
PROBE D1 -- THE DARK STATE'S NAME. Discovery-first. Which involution S commutes with M+ and makes the
R4 dark doublet-member odd? Wilson's conjecture: P (pair swap). HONEST CAVEAT: Thm7 has P as the INTERTWINER
M- = P M+ P, so [M+,P] = (M+ - M-)P != 0 unless M+=M- (B1: they differ by 0.25) -- P likely FAILS. Discover the
real annihilator. (M real => [M, S.C] = [M,S].C, so C-composites don't change the commutator; census = permutations.)

D1-A commutator census; D1-B parity of banked eigenvectors under commuting S; D1-C symmetry-resolved product
census (even/odd); D1-D S-action in ladder coords. Dense L=3 banked vectors; no near-EP extraction.
"""
import numpy as np, scipy.sparse as sp
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

np.set_printoptions(linewidth=160, suppress=True)


def full_M(q, L, lam=0.5):
    qL = q ** L; sub = subgroup(2, qL); ordn = len(sub)
    raw = [lam ** d for d in range(1, ordn + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    return M.tocsr(), idx, states, n, qL


def involutions(states, idx, qL):
    """candidate involutions on (a,b,gamma) -> perm array (perm[i]=idx of S(state_i)); only valid ones kept."""
    inv = lambda x: pow(x, -1, qL)
    defs = {
        "P  (swap a<->b)":        lambda a, b, g: (b, a, g),
        "G  (a,b -> a^-1,b^-1)":  lambda a, b, g: (inv(a), inv(b), g),
        "PG (swap+invert)":       lambda a, b, g: (inv(b), inv(a), g),
        "Gm (invert, gamma->-g)": lambda a, b, g: (inv(a), inv(b), (-g) % qL),
        "PGm(swap+inv, g->-g)":   lambda a, b, g: (inv(b), inv(a), (-g) % qL),
        "Sig(-a,-b,-g) [B1]":     lambda a, b, g: ((-a) % qL, (-b) % qL, (-g) % qL),
        "Z  (gamma->-gamma)":     lambda a, b, g: (a, b, (-g) % qL),
    }
    perms = {}
    for name, f in defs.items():
        ok = True; perm = np.empty(len(states), dtype=np.int64)
        for i, (a, b, g) in enumerate(states):
            key = f(a, b, g)
            if key not in idx:
                ok = False; break
            perm[i] = idx[key]
        if ok and np.all(perm[perm] == np.arange(len(states))):    # involution check
            perms[name] = perm
    return perms


def commutator_rel(M, perm):
    """||M S - S M||_F / ||M||_F  for permutation S (involution). MS = M[:,perm], SM = M[perm,:]."""
    MS = M[:, perm]; SM = M[perm, :]
    D = MS - SM
    return np.sqrt((D.multiply(D)).sum()) / np.sqrt((M.multiply(M)).sum())


def parity(v, perm):
    Sv = v[perm]
    num = np.vdot(v, Sv); den = np.vdot(v, v)
    return (num / den).real, abs((num / den).imag)


def main():
    print("# PROBE D1 -- THE DARK STATE'S NAME. Discovery-first commutator census + parity. L=3 dense.")
    M, idx, states, n, qL = full_M(3, 3)
    perms = involutions(states, idx, qL)
    Mnorm = np.sqrt((M.multiply(M)).sum())

    # ---- D1-A commutator census ----
    print(f"\n## D1-A COMMUTATOR CENSUS (||[M,S]||_F/||M||_F, L=3 full op dim {n}):")
    commuting = {}
    for name, perm in perms.items():
        c = commutator_rel(M, perm)
        tag = "COMMUTES" if c < 1e-12 else ("anti/neither" )
        print(f"   {name:26s}: {c:.3e}   {tag}")
        if c < 1e-12:
            commuting[name] = perm
    print(f"   => commuting involutions: {list(commuting.keys()) or 'NONE among candidates'}")

    # ---- dense eig + identify modes ----
    ev, VR = np.linalg.eig(M.toarray())
    one = np.ones(n); v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
    # doublet members (nearest to seat k=1), k=2 pair, partner, top band modes
    up = [i for i in range(n) if ev[i].imag > 1e-9]
    def seat(k):
        th = 2 * np.pi * k / 9; return (1 / 3) * np.cos(th / 2) ** 2 * np.exp(1j * th)
    dbl = sorted(up, key=lambda i: abs(ev[i] - seat(1)))[:2]
    k2 = sorted(up, key=lambda i: abs(ev[i] - seat(2)))[:2]
    reals = [i for i in range(n) if abs(ev[i].imag) < 1e-9]
    partner = min(reals, key=lambda i: abs(ev[i] - 0.333236))
    topband = sorted(up, key=lambda i: -abs(ev[i]))[:20]

    # ---- D1-B parity ----
    if commuting:
        for name, perm in commuting.items():
            print(f"\n## D1-B PARITY under {name}:")
            # verify functional & init are S-symmetric
            pf = parity(one, perm); pv = parity(v0, perm)
            print(f"   readout <1| S-parity = {pf[0]:+.4f} (imag {pf[1]:.1e});  init v0 S-parity = {pv[0]:+.4f}")
            labels = [("doublet m0", dbl[0]), ("doublet m1", dbl[1]),
                      ("k2 pair m0", k2[0]), ("k2 pair m1", k2[1]), ("partner", partner)]
            for lab, j in labels:
                p, im = parity(VR[:, j], perm)
                print(f"   {lab:12s} lam={ev[j].real:+.5f}{ev[j].imag:+.5f}j  S-parity={p:+.4f} (imag {im:.1e})  "
                      f"{'EVEN' if p > 0.5 else ('ODD' if p < -0.5 else 'mixed')}")
    else:
        print("\n## D1-B: NO commuting involution -- Wilson's P-conjecture DIES (walk-back #29). "
              "P=1.374 (Thm7 intertwiner, [M,P]=(M-M_)P!=0), G/PG/Sig/Z ~ sqrt(2) (no relation).")
        # verify functional & init are P-symmetric (the mechanism's other half survives even though P is not a symmetry of M)
        P = perms.get("P  (swap a<->b)")
        if P is not None:
            print(f"   functional <1| P-symmetric: {np.allclose(one[P], one)};  init v0 P-fixed: {np.allclose(v0[P], v0)}  "
                  f"(both TRUE -- but P not a symmetry of M, so eigenvecs are not P-eigenstates)")
        # DIAGNOSIS: the dark state is a READOUT zero <1|r>=0 (DC-free mode), not a symmetry eigenstate.
        evL, VL = np.linalg.eig(M.toarray().T)
        def overlaps(j):
            r = VR[:, j]; jl = int(np.argmin(np.abs(evL - ev[j]))); l = VL[:, jl]
            r1 = np.dot(one, r); lv = np.dot(l, v0); lr = np.dot(l, r)
            return abs(r1), abs(lv), abs(r1 * lv / lr)
        print("\n## D1 DIAGNOSIS -- the dark state = READOUT (agreement-functional) annihilation, not a symmetry:")
        print(f"   {'mode':12s} {'|lam|':>7} {'phase':>7} {'|<1|r>|':>10} {'|<l|v0>|':>10} {'|A|=|prod|/..':>12}  channel")
        rows = [("doublet m0", dbl[0]), ("doublet m1", dbl[1]), ("k2 m0", k2[0]), ("k2 m1", k2[1]),
                ("partner", partner)] + [(f"band{t}", topband[t]) for t in range(8)]
        for lab, j in rows:
            r1, lv, A = overlaps(j)
            chan = "DARK: <1|r>~0 (DC-free)" if r1 < 1e-8 else ("dark: <l|v0>~0" if lv < 1e-8 else "visible")
            print(f"   {lab:12s} {abs(ev[j]):7.4f} {np.angle(ev[j]):+7.4f} {r1:10.3e} {lv:10.3e} {A:12.3e}  {chan}")
        print("\n   => SELECTION RULE (corrected): visibility = |<1|r_j>| (readout/DC coupling), NOT a global "
              "symmetry parity. The 'null column' = modes with <1|r_j> ~ 0 (mean-zero / DC-free). The all-ones "
              "readout is the DC projector; DC-free band modes are dark to the agreement observable.")

    # ---- D1-C symmetry-resolved census (raw dump) ----
    if commuting:
        name, perm = next(iter(commuting.items()))
        print(f"\n## D1-C SYMMETRY-RESOLVED PRODUCT CENSUS under {name} (raw; product=A_j(lam-1/3), R4-D conventions)")
        evL, VL = np.linalg.eig(M.toarray().T)
        def amp(j):
            r = VR[:, j]; jl = int(np.argmin(np.abs(evL - ev[j]))); l = VL[:, jl]
            return (np.dot(one, r)) * (np.dot(l, v0)) / np.dot(l, r)
        print(f"   {'mode':12s} {'|lam|':>7} {'phase':>7} {'parity':>7} {'|product|':>11}")
        for lab, j in [("doublet m0", dbl[0]), ("doublet m1", dbl[1]), ("k2 m0", k2[0]),
                       ("k2 m1", k2[1]), ("partner", partner)] + [(f"band{t}", topband[t]) for t in range(8)]:
            p, _ = parity(VR[:, j], perm)
            try:
                pr = amp(j) * (ev[j] - 1 / 3)
            except Exception:
                pr = complex('nan')
            print(f"   {lab:12s} {abs(ev[j]):7.4f} {np.angle(ev[j]):+7.4f} {p:+7.3f} {abs(pr):11.3e}  "
                  f"{'[EVEN/visible]' if p > 0.5 else ('[ODD/null]' if p < -0.5 else '[mixed]')}")


if __name__ == "__main__":
    main()
