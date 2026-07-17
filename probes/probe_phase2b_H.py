"""
Real q=3 operator, Delta-channel structure (Nathan's derivation). Fires H1+H2.

Real operator M = build_M_gen(3, L, 2, raw=[lam^delta]) -- phase group <2> mod 3^L,
order D=2*3^(L-1). Folded weights w_delta = lam^delta / sum (=2^-delta/(1-2^-D) at
lam=1/2). Circulant family (equal-move channel on Delta={(a,a,0)}):
    c_k = sum_delta w_delta^2 * exp(2i pi k delta / D),  k=0..D-1,  c_0 = sum w^2.

VALIDATION GATE (must pass or nothing counts): gap(L=2) ~ 2.9e-3, gap(L=3) ~ 1.0e-4;
c_0 = (1/3)(1-4^-D)/(1-2^-D)^2 = 0.343915 (L=2) / 0.333336 (L=3); L=1 family = {5/9,-1/3}.

H1  (amplitude, true v0=indicator(1,1,0)): EXACT modal amplitude A_i of each eigenvalue
    in the mass sequence s_k=1^T M^k v0.  Pre-reg: the k!=0 circulant members carry
    ZERO amplitude (= R32/R25's zero-amplitude tower cluster, now closed-form-identified);
    amplitude lives on c_0 and the dynamical partner.  [exact eigendecomp, NOT ESPRIT.]
H2  (family completeness at L=3): sparse shift-invert at sigma=c_k confirms each c_k is in
    the spectrum (residual ~0), incl. c_2..c_5 that were window-limited in the dense read.
"""
import numpy as np
import scipy.sparse.linalg as spla

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

Q = 3


def folded_weights(D, lam=0.5):
    raw = np.array([lam ** d for d in range(1, D + 1)], float)
    return raw / raw.sum()


def circ_family(w, D):
    d = np.arange(1, D + 1)
    return np.array([np.sum(w ** 2 * np.exp(2j * np.pi * k * d / D)) for k in range(D)])


def real_M(L, lam=0.5):
    qL = Q ** L
    D = len(subgroup(2 % qL, qL))
    raw = [lam ** dd for dd in range(1, D + 1)]
    M, idx, n = build_M_gen(Q, L, 2, raw)
    return M, idx, n, D


def modal_amplitudes(Md, v0):
    """A_i for s_k = 1^T M^k v0 = sum_i A_i mu_i^k, exact via left/right eigvecs."""
    mu, R = np.linalg.eig(Md)
    muL, L_ = np.linalg.eig(Md.T)
    # match left eigvecs to right by eigenvalue
    order = []
    used = set()
    for i in range(len(mu)):
        j = min((k for k in range(len(muL)) if k not in used),
                key=lambda k: abs(muL[k] - mu[i]))
        used.add(j); order.append(j)
    Lm = L_[:, order]                     # column i is left eigvec for mu_i
    ones = np.ones(Md.shape[0])
    A = np.zeros(len(mu), complex)
    for i in range(len(mu)):
        li = Lm[:, i]; ri = R[:, i]
        denom = li @ ri
        A[i] = (ones @ ri) * (li @ v0) / denom if abs(denom) > 1e-14 else 0.0
    return mu, A


def main():
    print("# Real q=3 Delta-channel operator -- validation + H1 + H2\n")

    for L in [1, 2, 3]:
        M, idx, n, D = real_M(L)
        w = folded_weights(D)
        c = circ_family(w, D)
        c0_closed = (1 / 3) * (1 - 4.0 ** (-D)) / (1 - 2.0 ** (-D)) ** 2
        print(f"## L={L}  D={D}  dim={n}   c0=sum w^2={c[0].real:.6f} (closed {c0_closed:.6f}, "
              f"diff {abs(c[0].real-c0_closed):.1e})")

        if L <= 2:
            Md = M.toarray()
            ev = np.linalg.eigvals(Md)
            mods = np.sort(np.abs(ev))[::-1]
            gap = mods[0] - mods[1]
            # each c_k an eigenvalue?
            dists = [min(abs(ev - ck)) for ck in c]
            print(f"   top|eig| = {mods[0]:.6f}, 2nd = {mods[1]:.6f}, gap = {gap:.3e}")
            print(f"   c_k in spectrum? max_k min|eig-c_k| = {max(dists):.2e}  "
                  f"(family: {', '.join(f'{ck.real:.4f}' for ck in c)})")
            # H1 amplitudes
            v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
            mu, A = modal_amplitudes(Md, v0)
            print("   H1 -- amplitudes A_i (|A|>1e-9 shown; others=zero-amplitude cluster):")
            items = sorted(zip(mu, A), key=lambda t: -abs(t[1]))
            for m_, a_ in items:
                if abs(a_) > 1e-9:
                    # is this a circulant family member?
                    kk = int(np.argmin(np.abs(c - m_)))
                    tag = f"= c_{kk}" if abs(c[kk] - m_) < 1e-6 else "(non-family / dynamical)"
                    print(f"      mu={m_.real:+.6f}{m_.imag:+.6f}j  |A|={abs(a_):.4e}  {tag}")
            # explicit: amplitude ON each c_k
            print("   H1 -- amplitude carried by each circulant member c_k:")
            for k, ck in enumerate(c):
                i = int(np.argmin(np.abs(mu - ck)))
                print(f"      c_{k}={ck.real:+.6f}{ck.imag:+.6f}j : matched mu, |A|={abs(A[i]):.3e}"
                      + ("   <-- ZERO-amplitude" if abs(A[i]) < 1e-9 else "   <-- carries amplitude"))
        else:
            # L=3: gap via sparse top pair; H2 via ROBUST LU-pivot singularity test.
            # NB: shift-invert eigs(sigma=c_k) is UNRELIABLE here -- the q=3 operator is
            # DEFECTIVE (R39 Jordan), which breaks ARPACK just as it broke ESPRIT/G0c.
            # min|diag(U)| of LU(M - sigma I) ~ 0 <=> sigma is an eigenvalue; use a
            # sigma+0.05 control to prove the test discriminates.
            import scipy.sparse as sp
            top = spla.eigs(M, k=3, which='LM', return_eigenvectors=False, maxiter=8000)
            tmods = np.sort(np.abs(top))[::-1]
            print(f"   top|eig| (LM) = {tmods[0]:.6f}, 2nd = {tmods[1]:.6f}, gap = {tmods[0]-tmods[1]:.3e}")
            Mc = M.tocsc(); Ieye = sp.identity(n, format='csc')

            def minpiv(sigma):
                lu = spla.splu((Mc - sigma * Ieye).astype(complex))
                return float(np.abs(lu.U.diagonal()).min())

            print("   H2 -- family completeness via LU-pivot singularity (control = c_k+0.05):")
            hits = 0
            for k in range(D):
                p = minpiv(c[k]); pc = minpiv(c[k] + 0.05)
                ok = p < 1e-9
                hits += ok
                par = " PARITY" if k in (0, D // 2) else ""
                if k < 6 or k in (0, D // 2) or not ok:
                    print(f"      c_{k:<2}={c[k].real:+.6f}{c[k].imag:+.6f}j  minpiv={p:.1e}  control={pc:.1e}  [{'EIG' if ok else 'absent'}]{par}")
            print(f"   H2 verdict: {hits}/{D} circulant members are exact eigenvalues of the L=3 operator.")
        print()


if __name__ == "__main__":
    main()
