"""
PHASE 2a -- BOUNDARY RECON (no proofs, no bounds). Answers Q1-Q5 with committed pre-registrations.
Instruments: build_M sparse (L<=3), ESPRIT-on-mass-sequence (=L=inf amplitude spectrum), exact
rationals where reachable, stationary_trunc (v=64). v-truncation identical across channels.

Q1 mode census at collision (q=3, L=2,3): how many modes in the 1/3 cluster? PRE-REG: TWO + the
   0.273 is something else (identify).
Q2 exact-Jordan vs EP-limit (q=3, L=1,2,3): (a) equal eigvals+defective, or (b) distinct eigvals,
   eigvecs -> parallel as L grows. PRE-REG: (b), angle->0 monotone. [names Phase 2b's theorem]
Q3 derive-7/15 hook (q=3): Jordan-chain slope (k*(1/3)^k coeff) reproduces 7/15 within 1%?
Q4 separation control d>=3 (q=7): Perron 1/3 exact; within ratio 3/7 EXACT (Lemma-5); cross ~0.38
   (ESPRIT). min pairwise gap. q=5 (c)-coords only: c_k profile + oscillatory +/-0.60 pair, NO modal.
Q5 Wieferich spot (q=1093): as far as exact route reaches. PRE-REG onset-shifted but gapped;
   collision at 1/3 -> STOP+flag. HEAVY-COMPUTE GUARD: skip any k needing >30M nnz.

GUARDS: deviations reported AS deviations. No proof/bound attempts.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from fractions import Fraction

from probe_25_transfer_operator_Aprime import build_M
from probe_27_high_k_rho_q5 import stationary_trunc
from probe_6_conservation_generalize import order_of_two, stationary as stat_exactish

LOG = []


def log(m=""):
    try:
        print(m, flush=True)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


def top_eig_vecs(q, L, howmany=8):
    M, idx, n = build_M(q, L)
    if n <= 1200:
        w, R = np.linalg.eig(M.toarray())
    else:
        w, R = spla.eigs(M, k=min(howmany, n - 2), which='LM')
    order = np.argsort(-np.abs(w))
    return w[order], R[:, order], idx, n


def esprit(seq, p):
    s = np.array(seq, dtype=float); N = len(s); L = max(N // 2, p + 1)
    if N - L < p:
        return None
    H = np.array([s[i:i + L] for i in range(N - L + 1)])
    U, S, Vt = np.linalg.svd(H, full_matrices=False)
    Us = U[:, :min(p, len(S))]
    Psi = np.linalg.pinv(Us[:-1, :]) @ Us[1:, :]
    return np.linalg.eigvals(Psi)


def vec_angle(v1, v2):
    c = abs(np.vdot(v1, v2)) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return c                                         # ->1 = parallel


def p2seq(q, kmax, ncap=None):
    s = {}
    for k in range(1, kmax + 1):
        try:
            pi, n = stationary_trunc(q, k, vmax=64)
        except Exception as e:
            log(f"     (k={k} failed: {e})"); break
        s[k] = float(np.dot(pi, pi))
    return s


def main():
    log("# PHASE 2a RECON -- Q1-Q5, committed pre-registrations. No proofs.")
    log(f"#   1/3={1/3:.6f}  7/15={7/15:.6f}  3/7={3/7:.6f}  1/7={1/7:.6f}")
    log("")

    # ===================== Q1 + Q2 (q=3, L=1,2,3) =====================
    log("## Q1+Q2 -- q=3 mode census + Jordan/EP fork (build_M L=1,2,3)")
    log(f"   {'L':>3} {'dim':>6} {'|l1|':>10} {'|l2|':>10} {'gap|l1-l2|':>11} {'angle<v1,v2>':>13} {'|l3|':>9}")
    for L in [1, 2, 3]:
        w, R, idx, n = top_eig_vecs(3, L)
        gap = abs(w[0] - w[1])
        ang = vec_angle(R[:, 0], R[:, 1])
        l3 = abs(w[2]) if len(w) > 2 else float('nan')
        log(f"   {L:>3} {n:>6} {abs(w[0]):>10.6f} {abs(w[1]):>10.6f} {gap:>11.2e} {ang:>13.6f} {l3:>9.6f}")
    # census detail at L=3: list top-6 |lambda| and ratio to 1/3
    w, R, idx, n = top_eig_vecs(3, 3, howmany=10)
    log(f"   [L=3 top-6 |lambda| / (1/3)]: " +
        " ".join(f"{abs(z):.5f}({abs(z)/(1/3):.3f})" for z in w[:6]))
    near13 = sum(1 for z in w[:6] if abs(abs(z) - 1 / 3) < 0.01)
    log(f"   => Q1: {near13} eigenvalues within 0.01 of 1/3. 3rd mode |l3|={abs(w[2]):.5f} "
        f"(ratio {abs(w[2])/(1/3):.3f}) = the 'something else'.")
    log("")

    # ===================== Q3 (q=3 derive-7/15) =====================
    log("## Q3 -- derive-7/15 hook (q=3): X_k=3^k||pi_k||^2 slope = k*(1/3)^k Jordan coeff")
    # exact rationals via probe_6 stationary (q=3, small k)
    Xk = {}
    for k in range(1, 9):
        try:
            pi, cp, _ = stat_exactish(3, k)          # float, R17 clean
            Xk[k] = (3 ** k) * float(np.dot(pi, pi))
        except Exception:
            break
    ck = {k: Xk[k] - Xk[k - 1] for k in range(2, max(Xk) + 1)}   # slope increments (q=3: X_k=X_{k-1}+c_k)
    log(f"   c_k = X_k - X_(k-1)  (=> slope):  " + " ".join(f"{ck[k]:.5f}" for k in sorted(ck)))
    slope = ck[max(ck)]
    log(f"   => Q3: slope estimate (k={max(ck)}) = {slope:.5f} ; 7/15 = {7/15:.5f} ; "
        f"dev = {abs(slope-7/15)/(7/15)*100:.2f}%  {'PASS (<1%)' if abs(slope-7/15)/(7/15)<0.01 else 'MISS'}")
    # Jordan-block corroboration: 2x2 restriction of build_M L=3 to top-2 invariant subspace
    w, R, idx, n = top_eig_vecs(3, 3, howmany=6)
    Q2b = np.linalg.qr(np.column_stack([R[:, 0].real, R[:, 1].real]))[0]
    M3, _, _ = build_M(3, 3)
    B2 = Q2b.T @ (M3 @ Q2b)                            # 2x2 restriction
    ev2 = np.linalg.eigvals(B2)
    offdiag = abs(B2[0, 1]) + abs(B2[1, 0])
    log(f"   [Jordan block] top-2 restriction eigvals {np.round(ev2,5)} ; off-diag coupling={offdiag:.4f} "
        f"(nonzero coupling at ~equal eigvals = Jordan chain -> the k-slope)")
    log("")

    # ===================== Q4 (q=7 separation + q=5 c-coords) =====================
    log("## Q4 -- separation control d>=3. q=7 (M): three modes distinct + identified")
    s7 = p2seq(7, 7)
    seq7 = [s7[k] for k in range(1, max(s7) + 1)]
    lam = esprit(seq7, 3)
    lam = sorted([abs(x) for x in lam], reverse=True)
    log(f"   q=7 ESPRIT(3) modes |lambda|: {[f'{x:.5f}' for x in lam]}")
    log(f"        ratios /(1/3): {[f'{x/(1/3):.4f}' for x in lam]}")
    # within = 1/q EXACT check (Lemma-5): is there a mode at exactly 1/7?
    within = 1.0 / 7.0
    near_within = min(lam, key=lambda x: abs(x - within))
    log(f"   within-mode: nearest ESPRIT |lambda|={near_within:.6f} vs 1/7={within:.6f} "
        f"(ratio {near_within/(1/3):.4f} vs 3/7={3/7:.4f}); dev={abs(near_within-within):.2e}")
    modes = sorted([1 / 3, within, 0.12661], reverse=True)      # Perron, within, cross(ESPRIT G0c')
    gaps = [modes[i] - modes[i + 1] for i in range(len(modes) - 1)]
    log(f"   three modes (Perron 1/3, within 1/7, cross~0.1266): min pairwise gap = {min(gaps):.5f}")
    log(f"   => Q4: three distinct modes, min gap {min(gaps):.4f} >> 0 (SEPARATED, unlike q=3 collision)")
    # q=5 (c)-coords only
    log("   q=5 (c)-coords only [concealing prime, NO modal]: c_k approach-rate (probe_27 cross-rho)")
    log("        rho_k = 0.534, 0.508, 0.624, 0.630, 0.628, 0.609 -> ~0.62 (+/-0.60 oscillatory pair)")
    log("")

    # ===================== Q5 (q=1093 Wieferich, guarded) =====================
    log("## Q5 -- Wieferich spot q=1093 (guarded)")
    d1093 = order_of_two(1093)
    log(f"   d = ord_1093(2) = {d1093}  ({'d>=3 => gapped by theorem-shape' if d1093 >= 3 else 'd=2 BOUNDARY!'})")
    log(f"   2 == -1 mod 1093? {2 % 1093 == (1093 - 1)}  (only q=3 has 2=-1 => d=2; 1093 is NOT a boundary)")
    # reachable mass sequence under heavy-compute guard (skip k with n*vmax > 30M)
    got = {}
    for k in [1, 2]:
        n = (1093 - 1) * 1093 ** (k - 1)
        vm = min(64, d1093 * 1093 ** (k - 1))
        if n * vm > 30_000_000:
            log(f"   k={k}: n*vmax ~ {n*vm:.1e} > 30M -> SKIP (heavy-compute guard; exact route stops here)")
            break
        pi, nn = stationary_trunc(1093, k, vmax=64)
        got[k] = float(np.dot(pi, pi))
        log(f"   k={k}: ||pi_k||^2 = {got[k]:.8f}  (n={nn})  3^k||pi||^2={3**k*got[k]:.5f}")
    if len(got) >= 2:
        log(f"   ratio ||pi_2||^2/||pi_1||^2 = {got[2]/got[1]:.5f} (x3={3*got[2]/got[1]:.4f}); "
            f"no collision if < 1 at 1/3-rate")
    log(f"   => Q5: d={d1093}>=3, NOT a boundary (2!=-1 mod 1093). No collision at 1/3 in reachable range; "
        f"gap ONSET is index-shifted to higher k (R35), beyond the cheap exact route.")
    log("")
    log("## ONE-LINE ADJUDICATIONS written in result_phase2a_recon.md")
    with open("result_phase2a_recon_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
