"""
Verify the two Nathan-reported skeleton lemmas on the real q=3 operator
M = build_M_gen(3, L, 2, [lam^d]) (lam=1/2).

State classes:  Delta = {(a,a,0)},  C (carried-diagonal) = {(a,a,g): g!=0},
                O (off-diagonal) = {(a,b,g): a!=b}.
Edge src->dest iff M[dest,src] != 0 (build_M_gen stores M[dest,src]=weight).

L-A (no-return): NO path from carried-diagonal C back to Delta.
   one-step: M[Delta, C] == 0 ;  all-step: BFS from C never reaches Delta.
   (also report M[Delta, O] and whether Delta has ANY in-edge from outside Delta.)

L-B (gauge factorization, k=0 exact) -- verifying the CONCRETE CORE of the claim
   (my interpretation; the exact <=54/21 class reduction is Nathan's construction):
   gauge orbit = simultaneous phase shift (a,b,g)->(sa,sb,g), s in <2>; invariant
   label (rho=a*b^{-1} mod 3^L, g). Test: the k=0 (c0=sum w^2) eigenvector factors
   through the gauge quotient (constant on each (rho,g) orbit), while a k!=0 member
   does not. Also count total (rho,g) orbits and those reachable from Delta (cf. 54).
"""
import numpy as np
import scipy.sparse as sp

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

Q = 3


def real_M(L, lam=0.5):
    qL = Q ** L
    D = len(subgroup(2 % qL, qL))
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(Q, L, 2, raw)
    return M, idx, n, D, qL


def classify(idx, qL):
    Delta, C, O = set(), set(), set()
    for (a, b, g), i in idx.items():
        if a == b and g == 0:
            Delta.add(i)
        elif a == b:
            C.add(i)
        else:
            O.add(i)
    return Delta, C, O


def LA(L):
    M, idx, n, D, qL = real_M(L)
    Delta, C, O = classify(idx, qL)
    Mc = M.tocsr()
    Dl = np.array(sorted(Delta)); Cl = np.array(sorted(C)); Ol = np.array(sorted(O))
    # one-step weights into Delta from C and O
    w_from_C = Mc[Dl][:, Cl].sum()
    w_from_O = Mc[Dl][:, Ol].sum()
    non_delta = np.array(sorted(set(range(n)) - Delta))
    w_into_delta_ext = Mc[Dl][:, non_delta].sum()
    # BFS forward (src->dest) from all of C; can we reach Delta?
    # adjacency by columns: successors of src = rows with nonzero in column src
    Mcsc = M.tocsc()
    seen = set(C)
    frontier = list(C)
    reached_delta = set()
    while frontier:
        nxt = []
        for src in frontier:
            col = Mcsc.getcol(src)
            for dest in col.indices:
                if dest in Delta:
                    reached_delta.add(dest)
                if dest not in seen:
                    seen.add(dest); nxt.append(dest)
        frontier = nxt
    return dict(L=L, n=n, D=D, nDelta=len(Delta), nC=len(C), nO=len(O),
               w_from_C=float(w_from_C), w_from_O=float(w_from_O),
               w_into_delta_ext=float(w_into_delta_ext),
               bfs_reached_delta=len(reached_delta), reached_total=len(seen))


def gauge_label(a, b, g, qL):
    return (a * pow(b, -1, qL)) % qL, g


def LB(L, lam=0.5):
    M, idx, n, D, qL = real_M(L)
    Delta, C, O = classify(idx, qL)
    # circulant family
    raw = np.array([lam ** d for d in range(1, D + 1)]); w = raw / raw.sum()
    dd = np.arange(1, D + 1)
    c = np.array([np.sum(w ** 2 * np.exp(2j * np.pi * k * dd / D)) for k in range(D)])
    # eigenvectors of c_0 and c_1 via null space of (M - c_k I) (dense SVD; L=2 ok)
    Md = M.toarray()
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    # gauge orbits
    orb = {}
    for i, (a, b, g) in enumerate(states):
        orb.setdefault(gauge_label(a, b, g, qL), []).append(i)
    # reachable-from-Delta orbits
    Mcsc = M.tocsc(); seen = set(Delta); frontier = list(Delta)
    while frontier:
        nxt = []
        for src in frontier:
            for dest in Mcsc.getcol(src).indices:
                if dest not in seen:
                    seen.add(dest); nxt.append(dest)
        frontier = nxt
    reach_orbits = {gauge_label(*states[i], qL) for i in seen}

    def left_null(ck):
        # LEFT eigenvector (co-invariant functional) = right null of (M - ck I)^H
        U, S, Vh = np.linalg.svd((Md - ck * np.eye(n)).T.conj())
        return Vh.conj()[-1]

    def gauge_dev(vec):
        vec = vec / (np.abs(vec).max() + 1e-300)
        return max(np.max(np.abs(vec[m] - vec[m].mean())) for m in orb.values())

    # L-B tests the LEFT (co-invariant) eigenvector: k=0 factors through the gauge
    # quotient exactly; k!=0 does not.
    l0 = left_null(c[0]); dev0 = gauge_dev(l0)
    k1 = 1
    l1 = left_null(c[k1]); dev1 = gauge_dev(l1)
    l0n = l0 / (np.abs(l0).max() + 1e-300)
    delta_amp = (min(abs(l0n[i]) for i in Delta), max(abs(l0n[i]) for i in Delta))
    c_amp = max(abs(l0n[i]) for i in C)
    return dict(L=L, D=D, n=n, n_orbits=len(orb), n_reach_orbits=len(reach_orbits),
                c0=c[0].real, dev_c0=dev0, k1=k1, c1=c[k1], dev_c1=dev1,
                delta_amp=delta_amp, c_amp=c_amp)


def main():
    print("# Verify L-A (no-return) and L-B (gauge factorization k=0) -- real q=3 operator\n")
    print("## L-A -- no path from carried-diagonal C back to Delta")
    for L in [2, 3]:
        r = LA(L)
        ok_1 = r['w_from_C'] < 1e-14
        ok_ext = r['w_into_delta_ext'] < 1e-14
        ok_bfs = r['bfs_reached_delta'] == 0
        print(f"   L={L} dim={r['n']} (|D|={r['nDelta']} |C|={r['nC']} |O|={r['nO']}): "
              f"one-step w[Delta<-C]={r['w_from_C']:.1e} [{'0' if ok_1 else 'NONZERO'}], "
              f"w[Delta<-O]={r['w_from_O']:.1e}, w[Delta<-any-ext]={r['w_into_delta_ext']:.1e} "
              f"[{'Delta CLOSED from outside' if ok_ext else 'has ext in-edges'}]")
        print(f"        BFS from C ({r['nC']} states, reached {r['reached_total']} total): "
              f"Delta states reached = {r['bfs_reached_delta']}  "
              f"-> L-A {'HOLDS (no return)' if ok_bfs else 'FAILS'}")
    print()
    print("## L-B -- k=0 LEFT (co-invariant) eigenvector factors through gauge quotient (rho,g); k!=0 does not")
    for L in [2]:
        r = LB(L)
        print(f"   L={L} dim={r['n']}  gauge orbits (rho,g): total={r['n_orbits']}, "
              f"reachable-from-Delta={r['n_reach_orbits']}  (Nathan: object reduced to <=54)")
        print(f"      c0={r['c0']:.6f}: max within-orbit dev of LEFT eigvec = {r['dev_c0']:.2e}  "
              f"[{'GAUGE-INVARIANT (k=0 factors EXACTLY)' if r['dev_c0']<1e-8 else 'not gauge-invariant'}]")
        print(f"      ell_0 co-invariance: |val| on Delta in [{r['delta_amp'][0]:.3f},{r['delta_amp'][1]:.3f}] "
              f"(chi_0 constant), max|val| on carried-diagonal C = {r['c_amp']:.1e} (zero)")
        print(f"      c_{r['k1']}={r['c1'].real:+.4f}{r['c1'].imag:+.4f}j: max within-orbit dev = {r['dev_c1']:.2e}  "
              f"[{'gauge-invariant' if r['dev_c1']<1e-8 else 'NOT gauge-invariant (k=0 is special)'}]")


if __name__ == "__main__":
    main()
