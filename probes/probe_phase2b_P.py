"""
PROBE P -- the dynamical PARTNER's home in the carry coordinate. Recon only.
INSTRUMENT LAW: direct/LU + dense eig only near the EP (q=3). Iterative used ONLY
for the gapped q=7 control, where the EP rationale does not apply (justified inline).

Object: build_M_gen(3, L, 2, [lam^d]), lam=1/2. Partner = eigenvalue nearest the known
partner value (0.346827 @L=2, 0.333236 @L=3), confirmed NOT equal to any family c_k.

P1 gamma-profile: m(v) = sum_{v_q(gamma)=v} |ell|^2, v=0(units)..L(gamma=0), normalized. LEFT and RIGHT.
P2 gauge: does ell_partner factor as omega^{-e_a} * f(e_rho, gamma)? best twist k + within-orbit residual.
P3 transfer: aggregate carry-level flow A[v->v'] over transitions (raw + partner-weighted).
P4 control: q=7 partner gamma-profile (pre-reg: same tower-graded shape).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from probe_phase2a_q2b_q6 import build_M_gen, subgroup


def real_M(q, L, lam=0.5):
    qL = q ** L
    D = len(subgroup(2 % qL, qL))
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    wf = np.array(raw) / sum(raw)
    return M, idx, n, D, qL, states, wf


def vq(g, q, L):
    if g == 0:
        return L
    v = 0
    while g % q == 0:
        g //= q; v += 1
    return v


def circ(wf, D, k):
    om = np.exp(2j * np.pi * k / D)
    return np.sum(wf ** 2 * om ** (np.arange(D) + 1))


def gamma_profile(vec, states, q, L):
    m = np.zeros(L + 1)
    p = np.abs(vec) ** 2
    for i, (a, b, g) in enumerate(states):
        m[vq(g, q, L)] += p[i]
    s = m.sum()
    return m / s if s > 0 else m


def find_partner(evals, wf, D, known):
    cks = np.array([circ(wf, D, k) for k in range(D)])
    nonfam = [i for i in range(len(evals)) if min(abs(evals[i] - cks)) > 1e-9]
    i = min(nonfam, key=lambda j: abs(evals[j] - known))
    dist_to_ck = min(abs(evals[i] - cks))
    dist_excl_perron = min(abs(evals[i] - cks[k]) for k in range(1, D))  # excl c0
    return i, evals[i], dist_to_ck, dist_excl_perron


def dlog_table(qL, D):
    dl = {}; x = 1 % qL
    for e in range(D):
        dl[x] = e; x = (x * 2) % qL
    return dl


def gauge_resid(vec, states, qL, D, dl, k):
    # CORRECTED 2026-07-16 (was `om**(dl[a]*k)` = k^2 twist bug; Probe C caught it): detwist by
    # omega^{+k e_a} with omega=exp(2i pi/D) -> exp(2i pi k e_a / D). P2's "Z/3 sub-family" was the
    # k^2 artifact; the real signal is k=0 (partner ~ gauge-invariant). C1 supersedes this test.
    om = np.exp(2j * np.pi * k / D)
    g = np.array([vec[i] * om ** dl[a] for i, (a, b, gam) in enumerate(states)])  # omega^{k e_a}
    g = g / (np.abs(g).max() + 1e-300)
    orb = {}
    for i, (a, b, gam) in enumerate(states):
        orb.setdefault((dl[(b * pow(a, -1, qL)) % qL], gam), []).append(i)
    return max(np.max(np.abs(g[m] - g[m].mean())) for m in orb.values())


def transfer_table(M, states, q, L, weight=None):
    Mc = M.tocoo()
    A = np.zeros((L + 1, L + 1))
    for r, c, v in zip(Mc.row, Mc.col, Mc.data):
        vs = vq(states[c][2], q, L); vd = vq(states[r][2], q, L)
        w = v if weight is None else v * weight[r] * weight[c]
        A[vs, vd] += w
    return A


def partner_via_subspace(A, evals, vecs, ic0_val, ipart_val):
    """Robust eigvec of the near-defective pair: restrict A to span{u_c0,u_part} (reliable
    even at the EP), diagonalize the 2x2, return the partner eigenvector in full space."""
    i0 = min(range(len(evals)), key=lambda j: abs(evals[j] - ic0_val))
    ip = min(range(len(evals)), key=lambda j: abs(evals[j] - ipart_val))
    U = np.column_stack([vecs[:, i0], vecs[:, ip]])
    Q, _ = np.linalg.qr(U)
    B = Q.conj().T @ (A @ Q)                      # 2x2 restriction
    w, V = np.linalg.eig(B)
    jp = int(np.argmin(np.abs(w - ipart_val)))
    return Q @ V[:, jp], w[jp]


def run_q3(L, known, do_right=True):
    print(f"\n## q=3 L={L}  (partner near {known})")
    M, idx, n, D, qL, states, wf = real_M(3, L)
    dl = dlog_table(qL, D)
    c0 = circ(wf, D, 0).real
    Md = M.toarray()
    print(f"   dim={n} D={D}. dense eig; partner via 2D invariant-subspace restriction (EP-robust) ...", flush=True)
    evalsL, L_ = np.linalg.eig(Md.T)                 # LEFT = eig(M^T)
    ip, cp, dck, dexcl = find_partner(evalsL, wf, D, known)
    print(f"   PARTNER eigenvalue = {cp.real:.6f}{cp.imag:+.6f}j  "
          f"(dist to nearest c_k = {dck:.2e}; nearest c_k EXCLUDING c0 = {dexcl:.2e}) -> distinct")
    # robust LEFT partner (deflate the coalescing pair against c0)
    ellp, _ = partner_via_subspace(Md.T, evalsL, L_, c0, cp)
    del L_, evalsL                                    # free 1.2GB before the RIGHT eig (L=3)
    prof_L = gamma_profile(ellp, states, 3, L)
    print(f"   P1 LEFT  (subspace) gamma-profile m(v), v=0..{L}: " + ' '.join(f'{x:.4f}' for x in prof_L))
    pure0 = prof_L[L] > 0.999
    print(f"        tower-graded? {'NO -- PURE gamma=0 (would contradict Real-T1)' if pure0 else 'YES (mass on gamma!=0 as Real-T1 requires)'}")
    # RIGHT partner (the home)
    evalsR, R_ = np.linalg.eig(Md)
    rp, _ = partner_via_subspace(Md, evalsR, R_, c0, cp)
    del R_, evalsR
    prof_R = gamma_profile(rp, states, 3, L)
    print(f"   P1 RIGHT (subspace) gamma-profile m(v), v=0..{L}: " + ' '.join(f'{x:.4f}' for x in prof_R)
          + "   <- the partner's HOME")
    # P2 gauge (on the robust LEFT partner)
    resids = [(k, gauge_resid(ellp, states, qL, D, dl, k)) for k in range(D)]
    bk, br = min(resids, key=lambda t: t[1])
    print(f"   P2 gauge: best twist k={bk}, within-(e_rho,gamma)-orbit residual={br:.2e}  "
          f"[{'FORK(a): factors in ONE sector' if br < 1e-9 else 'FORK(b): NO exact factorization -> mixes sectors'}]")
    print(f"        residual by k (min group): " + ' '.join(f'k{k}:{r:.1e}' for k, r in resids if r < 0.5) or "        (all ~1.0)")
    # P3 transfer
    wabs = np.abs(ellp)
    Araw = transfer_table(M, states, 3, L)
    Awt = transfer_table(M, states, 3, L, weight=wabs)
    np.savetxt(f"outputs/partner_transfer_raw_q3_L{L}.tsv", Araw, fmt='%.6e', delimiter='\t')
    np.savetxt(f"outputs/partner_transfer_wt_q3_L{L}.tsv", Awt, fmt='%.6e', delimiter='\t')
    np.savetxt(f"outputs/partner_gammaprofile_q3_L{L}.tsv",
               np.vstack([prof_L, prof_R if prof_R is not None else prof_L * np.nan]), fmt='%.6e', delimiter='\t')
    print(f"   P3 raw v->v' transfer (rows=src level 0..{L}, cols=dest):")
    for v in range(L + 1):
        rowsum = Araw[v].sum()
        norm = Araw[v] / rowsum if rowsum > 0 else Araw[v]
        print(f"        v={v}: " + ' '.join(f'{x:.3f}' for x in norm) + f"   (rowmass {rowsum:.3e})")
    return cp, prof_L


def run_q7_control(L=2):
    print(f"\n## P4 CONTROL q=7 L={L} (GAPPED r7~0.38, far from EP -> iterative safe; justified)")
    M, idx, n, D, qL, states, wf = real_M(7, L)
    cks = np.array([circ(wf, D, k) for k in range(D)])
    c0 = cks[0].real
    fam_sub = max(abs(cks[k]) for k in range(1, D)) / c0     # family (autocorrelation) subdominant ratio
    print(f"   dim={n} D={D}. top-30 eigenpairs via ARPACK (gapped, not near EP) ...", flush=True)
    vals, vecs = spla.eigs(M.T, k=30, which='LM', maxiter=10000)   # left eigvecs (M^T)
    order = np.argsort(-np.abs(vals))
    sub = vals[order[1]]                                       # actual subdominant eigenvalue
    print(f"   family subdominant ratio max|c_k,k!=0|/c0 = {fam_sub:.4f}  (r7 lit ~0.38)")
    print(f"   ACTUAL subdominant eigenvalue = {sub.real:.6f}{sub.imag:+.6f}j  |.|/c0={abs(sub)/c0:.4f}  "
          f"dist to nearest c_k = {min(abs(sub-cks)):.2e}  -> {'FAMILY member' if min(abs(sub-cks))<1e-6 else 'NON-family'}")
    # is there ANY non-family mode in the top-30, and where does it sit?
    nonfam = [j for j in order if min(abs(vals[j] - cks)) > 1e-6]
    if not nonfam:
        print("   FINDING: NO non-family mode in the top-30 -- at q=7 the near-top is ALL family (autocorrelations).")
        print("            => the q=7 subdominant is KINEMATIC (a c_k), NOT a tower-partner. Pre-reg P4 REFUTED:")
        print("            the q=3 'dynamical partner' (non-family, tower-home) has NO near-top analog at q=7.")
        return
    j = nonfam[0]
    prof = gamma_profile(vecs[:, j], states, 7, L)
    print(f"   top NON-family mode: eig={vals[j].real:.6f}{vals[j].imag:+.6f}j  |.|/c0={abs(vals[j])/c0:.4f} "
          f"(rank {list(order).index(j)+1} by modulus)")
    print(f"   P4 gamma-profile m(v), v=0..{L}: " + ' '.join(f'{x:.4f}' for x in prof))
    np.savetxt(f"outputs/partner_gammaprofile_q7_L{L}.tsv", prof.reshape(1, -1), fmt='%.6e', delimiter='\t')
    pure0 = prof[L] > 0.999
    print(f"        tower-graded? {'NO -- pure gamma=0 (family-like)' if pure0 else 'YES -- tower home like q=3 (but it is NOT near-top: see rank)'}")


def main():
    print("# PROBE P -- the partner's home in the carry coordinate (recon; no proof, no rate fit)")
    run_q3(2, 0.346827, do_right=True)
    try:
        run_q3(3, 0.333236, do_right=False)   # left only at L=3 to bound memory
    except (MemoryError, np.core._exceptions._ArrayMemoryError) as e:
        print(f"\n## q=3 L=3: DENSE EIG WALLED ({type(e).__name__}) -- reported, not extrapolated. L=2 + q=7 stand.")
    run_q7_control(2)


if __name__ == "__main__":
    main()
