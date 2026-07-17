"""
REQUEST F -- e=-1 sub-block spectral radius rho_minus (closes D1 maximality).

Frozen toy M(q,-1,lambda) via build_M_gen. The e=ab=-1 sector = states (a,b,g)
with b == -a mod q^L. Restrict transitions to STAY in that sector by taking the
PRINCIPAL SUBMATRIX M[S,S] (drops the us-exit branches to e=+1; NO renormalize
-> strictly sub-stochastic by construction, per guard). rho_minus = spectral
radius of that block (dense exact eig; no ESPRIT, no mass sequences).

PRE-REG: P1 rho_- < |lam2|=(1-l)/(1+l) STRICT everywhere.
         P2 rho_- at us=l/(1+l)^2 class or below, NOT u^2=1/(1+l)^2 class.
         P3 rho_- q-flat over {5,7,13}, L-stable over {1,2}.
GRID q in {5,7,13} x L in {1,2} x lam in {0.30,0.50,0.70}.
"""
import numpy as np

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

QS = [5, 7, 13]
LS = [1, 2]
LAMS = [0.30, 0.50, 0.70]


def e_minus_block(q, L, lam):
    qL = q ** L
    gen = (-1) % qL
    ordn = len(subgroup(gen, qL))                 # =2 for gen=-1, any L
    raw = [lam ** d for d in range(1, ordn + 1)]   # frozen weights lam^delta -> [u,s]
    M, idx, n = build_M_gen(q, L, gen, raw)
    # S = the e=-1 sector: b == -a mod qL
    S = [i for s, i in idx.items() if s[1] == (-s[0]) % qL]
    S = sorted(S)
    Msub = M.toarray()[np.ix_(S, S)]              # principal submatrix: drops exits, no renorm
    states_S = [s for s in idx if idx[s] in set(S)]
    return Msub, S, idx, states_S


def rho_and_vec(Msub):
    w, V = np.linalg.eig(Msub)
    k = int(np.argmax(np.abs(w)))
    return abs(w[k]), w[k], V[:, k]


def main():
    print("# REQUEST F -- e=-1 sub-block spectral radius rho_minus")
    print("# principal submatrix on b==-a (exits DROPPED, not renormalized); dense exact eig\n")
    rows = []
    fails = []
    for lam in LAMS:
        u = 1.0 / (1 + lam)
        s = lam / (1 + lam)
        lam2 = (1 - lam) / (1 + lam)      # |lam2|
        us = lam / (1 + lam) ** 2
        u2 = 1.0 / (1 + lam) ** 2
        print(f"## lambda={lam:.2f}   |lam2|={lam2:.4f}   us={us:.4f}   u^2={u2:.4f}")
        for L in LS:
            for q in QS:
                Msub, S, idx, states_S = e_minus_block(q, L, lam)
                rho, ev, vec = rho_and_vec(Msub)
                p1 = "PASS" if rho < lam2 - 1e-12 else "FAIL"
                # class placement: distance to us vs u2
                cls = "us-class-or-below" if rho <= us + 1e-9 else ("between" if rho < u2 - 1e-9 else "u^2-class")
                print(f"   q={q:<2} L={L}  dim={Msub.shape[0]:<4} rho_-={rho:.6f}  "
                      f"(rho/|lam2|={rho/lam2:.3f}, rho/us={rho/us:.3f})  P1={p1}  [{cls}]")
                rows.append((lam, L, q, Msub.shape[0], rho, lam2, us, u2, p1, cls))
                if p1 == "FAIL":
                    # localize offending eigenvector (adjudication branch)
                    vv = np.abs(vec) / (np.abs(vec).max() + 1e-300)
                    supp = sorted([(states_S[k], round(vv[k], 3)) for k in range(len(vv)) if vv[k] > 0.1],
                                  key=lambda t: -t[1])[:8]
                    fails.append((lam, L, q, rho, lam2, supp))
        print()

    # P3 checks: spread across q at each (lam,L), and L-stability
    print("## P3 -- q-flatness (max-min rho over q) and L-stability")
    import itertools
    bykey = {}
    for (lam, L, q, dim, rho, lam2, us, u2, p1, cls) in rows:
        bykey.setdefault((lam, L), {})[q] = rho
    for (lam, L) in sorted(bykey):
        vals = [bykey[(lam, L)][q] for q in QS]
        print(f"   lam={lam:.2f} L={L}: rho over q{QS} = "
              f"[{', '.join(f'{v:.6f}' for v in vals)}]  spread={max(vals)-min(vals):.2e}")
    for lam in LAMS:
        for q in QS:
            r1 = bykey[(lam, 1)][q]; r2 = bykey[(lam, 2)][q]
            print(f"   lam={lam:.2f} q={q}: L1={r1:.6f} L2={r2:.6f}  |L2-L1|={abs(r2-r1):.2e}")

    print("\n## VERDICT")
    npass = sum(1 for r in rows if r[8] == "PASS")
    print(f"   P1: {npass}/{len(rows)} grid points PASS (rho_- < |lam2| strict).")
    if fails:
        print("   P1 FAILURES (offending eigenvector localization):")
        for (lam, L, q, rho, lam2, supp) in fails:
            print(f"     lam={lam} L={L} q={q}: rho_-={rho:.6f} >= |lam2|={lam2:.4f}; support={supp}")
    else:
        print("   -> maximality confirmed numerically; D1's diagonal-ray subdominant is the true max.")


if __name__ == "__main__":
    main()
