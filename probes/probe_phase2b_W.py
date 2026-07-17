"""
PROBE W -- reweighted compression (Probe E's named next lever). Recon; no proof, no rate fit.
INSTRUMENT LAW: dense/direct at q=3.

Probe E: the UNIFORM source-average compression puts the partner on the WRONG side of c0
(compressed 0.334312 above c0; true 0.333236 below), error 3.2e-3 >> true gap 1e-4. Fix:
reweight the source measure. Candidates for the source weight mu(src) of the lumped operator
   Lmat_mu[c,c'] = ( sum_{src in c} mu(src) sum_{dst in c'} M[dst,src] ) / ( sum_{src in c} mu(src) ):
 (U)  uniform mu=1                          -- E's baseline
 (Q)  quasi-stationary mu=|right Perron|    -- M's dominant right eigenvector (carry-tower measure)
 (D)  c0-DEFLATION: compress M - c0*r0 l0^T / (l0.r0), l0 = Real-T1 ell_0 (exact), then uniform lump
Metric: compressed-partner vs TRUE partner (0.346827@L2, 0.333236@L3) -- distance AND correct side of c0.
"""
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
    gam = np.array([s[2] for s in states])
    w = np.array(raw) / sum(raw)
    return M, states, n, D, qL, erho, gam, w, dl


def circ(w, D, k):
    return complex(np.sum(w ** 2 * np.exp(2j * np.pi * k * (np.arange(D) + 1) / D)))


def subspace_pair(A, evals, vecs, c0, partner_val):
    """clean eigvecs of the near-degenerate {c0, partner} pair via 2x2 restriction (EP-robust)."""
    i0 = int(np.argmin(np.abs(evals - c0)))
    ip = int(np.argmin(np.abs(evals - partner_val)))
    U = np.column_stack([vecs[:, i0], vecs[:, ip]])
    Q, _ = np.linalg.qr(U)
    B = Q.conj().T @ (A @ Q)
    w2, V2 = np.linalg.eig(B)
    j0 = int(np.argmin(np.abs(w2 - c0))); jp = int(np.argmin(np.abs(w2 - partner_val)))
    return (Q @ V2[:, j0], w2[j0]), (Q @ V2[:, jp], w2[jp])


def ell0_realT1(states, qL, D, w, dl):
    """Real-T1 left eigenvector for k=0 (no twist): ell_0 = R_0(e_rho)/R_0(0) on gamma=0."""
    dd = np.arange(1, D + 1)
    R = np.array([np.sum(w * np.roll(w, -e) * (np.ones(D)) ) for e in range(D)])  # omega=1 at k=0
    R0 = R[0]
    l0 = np.zeros(len(states))
    for i, (a, b, g) in enumerate(states):
        if g == 0:
            l0[i] = R[dl[(b * pow(a, -1, qL)) % qL]] / R0
    return l0


def weighted_lump(M, cls_idx, m, mu):
    num = np.zeros((m, m)); den = np.zeros(m)
    Mc = M.tocoo()
    np.add.at(num, (cls_idx[Mc.col], cls_idx[Mc.row]), mu[Mc.col] * Mc.data)
    np.add.at(den, cls_idx, mu)
    den[den == 0] = 1.0
    return num / den[:, None]


def rank1_lump(cls_idx, m, mu, lvec, rvec, den):
    A = np.zeros(m); B = np.zeros(m)
    np.add.at(A, cls_idx, mu * lvec)     # per class: sum mu*l0
    np.add.at(B, cls_idx, rvec)          # per class: sum r0
    return np.outer(A / den, B)


def compressed_partner(Lmat, c0, cks, partner_val):
    ev = np.linalg.eigvals(Lmat)
    i0 = int(np.argmin(np.abs(ev - c0)))
    cand = [i for i in range(len(ev)) if i != i0 and min(abs(ev[i] - cks)) > 1e-9]
    ip = min(cand, key=lambda i: abs(ev[i] - partner_val)) if cand else None
    return ev[ip] if ip is not None else None, ev[i0]


def run(L, partner_val):
    print(f"\n## q=3 L={L}  true partner={partner_val}, c0=sum w^2")
    M, states, n, D, qL, erho, gam, w, dl = real_M(3, L)
    c0 = circ(w, D, 0).real
    cks = np.array([circ(w, D, k) for k in range(D)])
    cls = list(zip(erho.tolist(), gam.tolist()))
    classes = sorted(set(cls)); ci = {c: j for j, c in enumerate(classes)}; m = len(classes)
    cls_idx = np.array([ci[c] for c in cls])
    Md = M.toarray()
    # right + left eigenpairs for {c0, partner} (clean, subspace)
    evR, VR = np.linalg.eig(Md); evL, VL = np.linalg.eig(Md.T)
    (r0, _), (rp, _) = subspace_pair(Md, evR, VR, c0, partner_val)
    (l0s, _), (lp, _) = subspace_pair(Md.T, evL, VL, c0, partner_val)
    l0 = ell0_realT1(states, qL, D, w, dl)              # exact Real-T1 ell_0
    truep = partner_val; side_true = "below" if truep < c0 else "above"
    print(f"   c0={c0:.6f}; true partner {side_true} c0 by {abs(truep-c0):.2e}")

    # (U) uniform
    LU = weighted_lump(M, cls_idx, m, np.ones(n))
    pU, _ = compressed_partner(LU, c0, cks, partner_val)
    # (Q) quasi-stationary: mu = |M's dominant right eigenvector|
    dom = VR[:, int(np.argmax(np.abs(evR)))]
    muQ = np.abs(dom); muQ = muQ / muQ.max()
    LQ = weighted_lump(M, cls_idx, m, muQ)
    pQ, _ = compressed_partner(LQ, c0, cks, partner_val)
    # (D) c0-deflation with exact Real-T1 l0 and clean r0
    den = np.zeros(m); np.add.at(den, cls_idx, np.ones(n)); den[den == 0] = 1
    r0r = np.real(r0); l0r = l0                          # l0 exact real; r0 real part
    norm = float(l0r @ r0r)
    LD = LU - c0 * rank1_lump(cls_idx, m, np.ones(n), l0r, r0r, den) / (norm if abs(norm) > 1e-12 else 1.0)
    pD, _ = compressed_partner(LD, c0, cks, partner_val)

    def rep(tag, p):
        if p is None:
            print(f"   {tag}: no distinct partner"); return
        side = "below" if p.real < c0 else "above"
        ok = (side == side_true) and (abs(p - truep) / abs(truep) < abs(pU - truep) / abs(truep) if tag != "U" else True)
        print(f"   {tag}: compressed-partner={p.real:+.6f}{p.imag:+.4f}j  ({side} c0)  "
              f"rel err to true={abs(p-truep)/abs(truep):.2e}  side_correct={side==side_true}")
    rep("U(uniform)", pU)
    rep("Q(quasi-stat)", pQ)
    rep("D(c0-deflate)", pD)


def main():
    print("# PROBE W -- reweighted compression: get the partner on the correct side of c0")
    run(2, 0.346827)
    run(3, 0.333236)


if __name__ == "__main__":
    main()
