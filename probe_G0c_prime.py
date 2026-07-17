"""
GATE G0c' (CARRIED -- blocks Phase 2, not Phase 1). Modal r_q = |lambda_2|/lambda_1 on build_M's
AMPLITUDE-CARRYING subspace, done right. q=7 primary (L=2 known bad 0.475), q=5 control (0.60-0.63).

FEASIBILITY: build_M at L=3 walls out (q7: dim 7.4M, ~1.6e11 nnz; q5: 1.25M, ~1e10 nnz). FIX (better
than L=3): by the mass identity Sum M^k v0 = ||pi_k||^2, the sequence ||pi_k||^2 = Sum_i A_i lambda_i^k
over build_M's AMPLITUDE-CARRYING eigenvalues (A_i=(1.r_i)(l_i.v0); tower/zero-amplitude modes are
absent from the sequence by construction). MATRIX-PENCIL/ESPRIT on the exact ||pi_k||^2 sequence
recovers {lambda_i, A_i} at effectively L=inf -- the R32 modal, no operator built.

THEORY (why L=2 raw modal missed): the amplitude-carrying eigenvalues are {1/3 (Perron), 1/q
(within-cell), r_q/3 (cross)}. |lambda|/lambda_1: within-> 3/q, cross-> r_q. At q=7, 3/7=0.4286 >
r_7=0.39, so the WITHIN mode is |lambda_2| and r_q=0.39 is |lambda_3|. The gate reads the CROSS
mode (the one ratio_within removes in probe_27), identified by ESPRIT among amplitude-carrying modes.

PRE-REG PASS BARS (committed; deviations reported AS deviations):
  q=7: cross-mode |lambda|/lambda_1 = 0.39 +/- 0.05.   q=5: 0.60-0.63.
  FALLBACK if ESPRIT ill-conditioned: L=1->2 build_M trend + probe_27 cross-rho as the anchor.
  ON FAIL (settled off by >0.05 or trend away): (M)-spectral framing miscalibrated -> report & STOP,
  Phase 2 restated in (c)-coordinates. NO patching.

SECONDARY (non-gating recon): q=3 build_M L=2,3 -- is the top pair exactly DEFECTIVE (no eigenbasis)
or distinct-and-merging (R39 vs R32)? Report eigenvector-matrix conditioning / pair-angle near lambda_1.

GUARDS: v-truncation identical all channels (vmax=64). q=3 alone never counts. Format per R26/R32.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from probe_27_high_k_rho_q5 import stationary_trunc
from probe_25_transfer_operator_Aprime import build_M
from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    try:
        print(m, flush=True)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def p2_sequence(q, kmax):
    s = {}
    for k in range(1, kmax + 1):
        pi, n = stationary_trunc(q, k, vmax=64)
        s[k] = float(np.dot(pi, pi))
        log(f"     ||pi_{k}||^2 = {s[k]:.14f}   (n={n})")
    return s


def esprit(seq, p):
    """ESPRIT: recover p geometric modes lambda_i of seq[k]=sum A_i lambda_i^k. seq is 1-indexed list."""
    s = np.array(seq, dtype=float)
    N = len(s)
    L = N // 2
    if L <= p:
        L = p + 1
    if N - L < p:
        return None, None
    # Hankel (N-L+1) x L
    H = np.array([s[i:i + L] for i in range(N - L + 1)])
    U, S, Vt = np.linalg.svd(H, full_matrices=False)
    p_eff = min(p, len(S))
    Us = U[:, :p_eff]
    U1 = Us[:-1, :]; U2 = Us[1:, :]
    Psi = np.linalg.pinv(U1) @ U2
    lam = np.linalg.eigvals(Psi)
    # amplitudes via Vandermonde least squares over k=1..N
    ks = np.arange(1, N + 1)
    Vd = np.vstack([lam ** k for k in ks])          # N x p
    A, *_ = np.linalg.lstsq(Vd, s, rcond=None)
    return lam, A


def main():
    log("# GATE G0c' -- amplitude-carrying modal r_q via ESPRIT on exact ||pi_k||^2 (L=inf). vmax=64.")
    log(f"#   PRE-REG: q7 cross-mode |l|/l1 = 0.39+/-0.05 ; q5 = 0.60-0.63.  (3/q ref: 3/7={3/7:.4f} 3/5={3/5:.3f})")
    log("")

    KMAX = {5: 8, 7: 7}
    for q in KMAX:
        log(f"## q={q}  (d={order_of_two(q)})  --- ESPRIT modal on ||pi_k||^2 ---")
        seq = p2_sequence(q, KMAX[q])
        s = [seq[k] for k in range(1, KMAX[q] + 1)]
        lam1 = 1.0 / 3.0
        for p in [2, 3, 4]:
            lam, A = esprit(s, p)
            if lam is None:
                log(f"   p={p}: too few points")
                continue
            order = np.argsort(-np.abs(lam))
            log(f"   p={p} modes (|lam|, |lam|/(1/3), A):")
            for i in order:
                log(f"        |lam|={abs(lam[i]):.5f}  ratio={abs(lam[i])/lam1:+.4f}  "
                    f"A={A[i].real:+.3e}{'' if abs(A[i].imag)<1e-6 else f'{A[i].imag:+.1e}j'}")
        # identify cross mode: among modes with |A| non-negligible, the one nearest r_q*(1/3),
        # i.e. NOT 1/3 (Perron) and NOT 1/q (within). Report all ratios for adjudication.
        log(f"   NB within-mode ratio 3/q={3/q:.4f} (=|lam|/l1 for lambda=1/q); cross-mode ratio = r_q.")
        log("")

    # ---- requested fallback: L=1->2 build_M amplitude trend (q5,q7); L=3 walls out ----
    log("## FALLBACK -- build_M raw+amplitude modal, L=1,2 (L=3 infeasible: q7 7.4M dim / ~1.6e11 nnz)")
    for q in [5, 7]:
        log(f"   q={q}:")
        for L in [1, 2]:
            M, idx, n = build_M(q, L)
            v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
            kk = min(30, n - 2)
            w = spla.eigs(M, k=kk, which='LM', return_eigenvectors=False)
            w = sorted(w, key=lambda z: -abs(z))
            l1 = abs(w[0])
            ratios = [abs(z) / l1 for z in w[:6]]
            log(f"     L={L} (dim={n}): raw |lam|/l1 top6 = {[f'{r:.4f}' for r in ratios]}")
    log("")

    # ---- SECONDARY (recon, non-gating): q=3 defectiveness at L=2,3 ----
    log("## SECONDARY (recon) -- q=3 build_M top-pair: defective (no eigenbasis) vs distinct-merging?")
    for L in [2, 3]:
        M, idx, n = build_M(3, L)
        A = M.toarray() if n <= 9000 else None
        if A is None:
            log(f"   q=3 L={L} dim={n} too large for dense -- skipped")
            continue
        w, R = np.linalg.eig(A)
        order = np.argsort(-np.abs(w))
        w = w[order]; R = R[:, order]
        # top-pair angle + eigenvector-matrix conditioning near top
        v1 = R[:, 0] / np.linalg.norm(R[:, 0])
        v2 = R[:, 1] / np.linalg.norm(R[:, 1])
        cos12 = abs(np.vdot(v1, v2))
        condR = np.linalg.cond(R)
        # conditioning of just the top-few eigenvectors
        cond_top = np.linalg.cond(R[:, :4])
        log(f"   q=3 L={L} (dim={n}): |lam1|={abs(w[0]):.5f} |lam2|={abs(w[1]):.5f} ratio={abs(w[1])/abs(w[0]):.5f}")
        log(f"        top-pair |<v1,v2>|={cos12:.3e} (->1 = merging/defective) ; cond(R)={condR:.2e} ; "
            f"cond(top4)={cond_top:.2e}")
    log("")
    log("## VERDICT (per pre-reg bars) written after inspecting ESPRIT cross-mode ratios above.")
    with open("result_G0c_prime_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("[DONE] result_G0c_prime_log.txt written")


if __name__ == "__main__":
    main()
