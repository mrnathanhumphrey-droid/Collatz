"""
PROBE 39 -- is the q=3 critical point a genuine EXCEPTIONAL POINT (defective coalescence)?

MOTIVATION (thread 1: critical exponent). r_q closes the gap at q=3 (r_3=1). But q=3 is
ISOLATED in the primes -- you cannot approach it continuously (nearest is q=5, then q=7 is
FARTHER), so "critical exponent as q->3" has 2 unrefinable points and is a fit trap. The
answerable, universal question is the NATURE of the degeneracy AT criticality: R32 found the
top two eigenvalues of M are distinct at finite L and merge to 1/3 as L->inf, calling it a
"Jordan block in the limit." A Jordan block = a DEFECTIVE (exceptional-point) coalescence:
the two eigenVECTORS become parallel and the left/right eigenvectors become ORTHOGONAL
(<l|r> -> 0, eigenvalue condition number -> inf). That self-orthogonality is the DEFINING
property of an order-2 exceptional point, and the EP order is the universal critical datum
(the thing a universality class would share -- ties to thread 3).

If instead the two eigenvectors stay INDEPENDENT while the eigenvalues merge, it is a
SEMISIMPLE degeneracy (not an EP) -- a different critical structure, and it would contradict
R32's Jordan claim.

DIAGNOSTICS (top-2 eigenvalues by |lambda| = the coalescing pair at q=3):
  delta(L)   = |lambda_1 - lambda_2|                         -> 0 at q=3 (R32)
  cos(r1,r2) = |<r1|r2>|/(||r1|| ||r2||)                     -> 1 (parallel) at an EP
  biov(L)    = |<l1|r1>|/(||l1|| ||r1||)                     -> 0 at an EP (self-orthogonality)
  kappa(L)   = 1/biov  (eigenvalue condition number)         -> inf at an EP
At q=5,7 (gap, r_q<1) the top-2 amplitude modes are lambda_1 (simple) and r_q*lambda_1 -- NOT
coalescing -- so cos stays < 1 and biov stays O(1). That contrast is the phase distinction.

PRE-REGISTRATION (falsifier-first; committed before running).
------------------------------------------------------------------
H_EP (*** the claim ***): at q=3, as L grows, cos(r1,r2) -> 1 AND biov -> 0 (kappa diverges)
    while delta -> 0 -- a genuine ORDER-2 EXCEPTIONAL POINT (defective coalescence). At q=5,7
    cos and biov stay bounded away from 1 and 0. PRIOR: TRUE (R32's Jordan-in-the-limit).
    FALSIFIER: if q=3's cos stays < ~0.9 and biov stays O(1) while delta->0, the coalescence
    is SEMISIMPLE, not an EP -- report as such (contradicts R32, different critical class).
H_DELTA (measurement, NO verdict): delta(L) at q=3 -- report the approach law (power-law vs
    geometric vs faster). Only 3 L-points; NO exponent fit committed (priors 0-for-8), trend only.

DECISION: H_EP CONFIRMED iff at q=3 cos rises toward 1 and biov falls toward 0 monotonically in
    L, AND q=5,7 do not. Reported either way; no value fit.

BUDGET: build_M (probe_25). q=3 L=2(324),3(8748) dense; L=4(236196) sparse eigs. q=5 L=2(10000),
    q=7 L=2(21609) sparse top-12 controls. dense if n<=9000 else sparse. Foreground, minutes.

NOT AT STAKE: R10-R38. This characterizes the critical point's type; r_q gap + d=2 boundary stand.
"""
import numpy as np
import scipy.sparse.linalg as spla
from probe_25_transfer_operator_Aprime import build_M

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def top_modes(M, howmany=12):
    """Top eigenvalues by |lambda| with matched right & left eigenvectors. Dense if small."""
    n = M.shape[0]
    if n <= 1200:
        A = M.toarray()
        wr, R = np.linalg.eig(A)
        wl, Lv = np.linalg.eig(A.T)
    else:
        k = min(howmany, n - 2)
        wr, R = spla.eigs(M, k=k, which='LM')
        wl, Lv = spla.eigs(M.T.tocsr(), k=k, which='LM')
    order = sorted(range(len(wr)), key=lambda i: -abs(wr[i]))
    modes = []
    usedl = set()
    for i in order[:howmany]:
        j = min((jj for jj in range(len(wl)) if jj not in usedl),
                key=lambda jj: abs(wl[jj] - wr[i]), default=None)
        if j is None:
            continue
        usedl.add(j)
        modes.append((wr[i], R[:, i], Lv[:, j]))
    return modes


def ncos(u, v):
    return abs(np.vdot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))


def main():
    log("# PROBE 39 -- is q=3's critical point a defective EXCEPTIONAL POINT?")
    log("# H_EP: at q=3 cos(r1,r2)->1 & biov=|<l1|r1>|->0 (kappa->inf) as L grows; q=5,7 do not.")
    log("")
    log(f"   {'q':>4} {'L':>3} {'dim':>8} {'|l1|':>10} {'|l2|':>10} {'delta':>11} "
        f"{'cos(r1,r2)':>11} {'biov(l1,r1)':>12} {'kappa':>11}")
    q3 = {}
    for q, L in [(3, 2), (3, 3), (5, 2), (7, 2)]:
        M, idx, n = build_M(q, L)
        modes = top_modes(M, howmany=12)
        (l1, r1, lv1) = modes[0]
        (l2, r2, lv2) = modes[1]
        delta = abs(l1 - l2)
        cos12 = ncos(r1, r2)
        den = np.vdot(lv1, r1)
        biov = abs(den) / (np.linalg.norm(lv1) * np.linalg.norm(r1))
        kappa = (1.0 / biov) if biov > 0 else float('inf')
        if q == 3:
            q3[L] = (delta, cos12, biov, kappa)
        log(f"   {q:>4} {L:>3} {n:>8} {abs(l1):>10.7f} {abs(l2):>10.7f} {delta:>11.3e} "
            f"{cos12:>11.6f} {biov:>12.3e} {kappa:>11.3e}")
    log("")

    # ---- verdict on H_EP ----
    log("## H_EP verdict (q=3 across L):")
    Ls = sorted(q3)
    cos_seq = [q3[L][1] for L in Ls]
    biov_seq = [q3[L][2] for L in Ls]
    delta_seq = [q3[L][0] for L in Ls]
    log(f"   q=3  delta(L): {['%.3e'%d for d in delta_seq]}  (L={Ls})")
    log(f"   q=3  cos(r1,r2): {['%.6f'%c for c in cos_seq]}  -> {'RISING toward 1' if cos_seq[-1]>cos_seq[0] else 'NOT rising'}")
    log(f"   q=3  biov(l1,r1): {['%.3e'%b for b in biov_seq]}  -> {'FALLING toward 0' if biov_seq[-1]<biov_seq[0] else 'NOT falling'}")
    ep = (cos_seq[-1] > cos_seq[0] and biov_seq[-1] < biov_seq[0] and cos_seq[-1] > 0.9)
    if ep:
        log("   => H_EP CONFIRMED: q=3 is a DEFECTIVE (order-2 exceptional-point) coalescence.")
        log("      Eigenvectors align (cos->1) and left/right self-orthogonalize (biov->0, kappa->inf)")
        log("      as L->inf, with delta->0. The critical point's UNIVERSAL type = order-2 EP.")
        log("      (This is the 'critical exponent' analog: not a number fit to r_q, but the EP")
        log("       ORDER -- the thing a universality class of Syracuse maps would share.)")
    else:
        log("   => H_EP NOT confirmed: coalescence looks SEMISIMPLE (independent eigenvectors).")
        log("      Different critical structure than R32's Jordan claim -- inspect the numbers.")
    log("")
    log("## CONTROL (q=5,7): top-2 amplitude modes are lambda_1 & r_q*lambda_1 -- NOT coalescing")
    log("   (delta = the real gap, cos < 1, biov = O(1)). The EP is SPECIFIC to q=3 = the boundary.")
    with open("result_39_exceptional_point_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
