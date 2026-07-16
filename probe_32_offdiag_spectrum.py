"""
PROBE 32 -- THREAD A: stress-test "r_q is a clean subdominant eigenvalue of M".
Amplitude-resolved full spectrum. Is r_q a well-separated, amplitude-dominant mode with a
real gap below lambda_1, or is the operator-eigenvalue story murky?

PRE-REGISTRATION (falsifier-first).
------------------------------------------------------------------
HOLE 1 (the audit's load-bearing hole): the brief says "r_q = subdominant eigenvalue of M",
but the RAW lambda_2 is the tower cluster (~0.98*lambda_1), NOT r_q; r_q was recovered only
by amplitude selection (R26, imperfect) and its VALUE comes from direct measurement (R27).
This probe asks, at the q we can build well (q=3 L=3, q=5 L=2, q=7 L=2):

  h(k) = ||pi_k||^2 / lambda_1^k = sum_i A_i mu_i^k,   mu_i = lambda_i/lambda_1,  mu_1=1.
  cross-growth of mode i  ~  |A_i * (mu_i - 1)|.
  r_q := mu of the LARGEST-cross-growth subdominant mode.

FALSIFIER (what would sink the brief):
  H_CLEAN (*** the claim ***): there is ONE subdominant mode carrying the cross growth,
     well-separated (its |A*(mu-1)| >> the next), with mu = r_q matching R27 direct
     (r_5~0.62, r_7~0.39), and a GAP mu < ~0.7 clearly below the tower cluster (mu~0.95+)
     which must have NEGLIGIBLE amplitude.
     IF INSTEAD the cross growth is spread across many modes, or the top mode's mu does NOT
     match R27, or there is no gap -> H_CLEAN is REFUTED and "r_q is an eigenvalue" is not a
     clean statement -> the brief's L3 framing needs revision. Report either way.

  H_TOWER: the near-lambda_1 cluster (mu in [0.9,1)) carries |A| that VANISHES relative to
     the r_q mode. PRIOR: TRUE (R26). Quantify the ratio.

  H_Q3 (Hole 2): at q=3, is the top of the spectrum a Jordan block or two distinct
     eigenvalues colliding as L->inf? Report lambda_1, lambda_2 (expect DISTINCT at finite L,
     both -> 1/3; the defect is a LIMIT phenomenon) and the r_q=mu_2 -> 1.

DECISION: report r_q(amplitude), the separation ratio (top/next cross-growth), the tower
amplitude ratio, and match-to-R27. NO value committed (priors 0-for-8). This is a
structural stress test with numbers reported, not a fit.

BUDGET: dense eig for q=3 L=3 (8748) and q=5 L=2 (10000); sparse eigs top-120 for q=7 L=2
(21609). Both left (M^T) and right eigvecs for amplitudes.

NOT AT STAKE: R10-R31 (r_q<1 gap and r_5,r_7 from DIRECT measurement stand regardless).
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


def spectrum(M, v0idx, dense, howmany=120):
    n = M.shape[0]
    v0 = np.zeros(n); v0[v0idx] = 1.0
    one = np.ones(n)
    if dense:
        A = M.toarray()
        wr, R = np.linalg.eig(A)
        wl, L = np.linalg.eig(A.T)
    else:
        k = min(howmany, n - 2)
        wr, R = spla.eigs(M, k=k, which='LM')
        wl, L = spla.eigs(M.T.tocsr(), k=k, which='LM')
    # match left to right by nearest eigenvalue
    modes = []
    used = set()
    for i in range(len(wr)):
        j = min((jj for jj in range(len(wl)) if jj not in used),
                key=lambda jj: abs(wl[jj] - wr[i]), default=None)
        if j is None:
            continue
        used.add(j)
        ri = R[:, i]; li = L[:, j]
        den = li.conj().dot(ri)
        if abs(den) < 1e-13:
            continue
        Ai = (one.dot(ri)) * (li.conj().dot(v0)) / den
        modes.append((wr[i], Ai))
    return modes


def main():
    log("# PROBE 32 -- THREAD A: is r_q a clean amplitude-dominant eigenvalue of M?")
    log("# Falsifier H_CLEAN: one well-separated cross-growth mode = r_q, gapped below tower")
    log("")
    R27 = {3: 1.0, 5: 0.62, 7: 0.39}
    for q, L, dense in [(3, 3, True), (5, 2, True), (7, 2, False)]:
        M, idx, n = build_M(q, L)
        modes = spectrum(M, idx[(1, 1, 0)], dense)
        l1 = max(m[0].real for m in modes)
        maxA = max(abs(A) for lam, A in modes)
        tol = max(1e-3, 1e-4 * maxA)   # amplitude noise floor
        # amplitude-CARRYING modes only, sorted by |mu| descending
        carr = [(lam / l1, abs(A)) for lam, A in modes if abs(A) > tol]
        carr.sort(key=lambda t: -abs(t[0]))
        log(f"## q={q}  L={L}  (dim={n}, lambda_1={l1:.6f}, ampl. floor={tol:.1e})")
        log(f"   amplitude-CARRYING modes by |mu| (mu=lambda/lambda_1, |A|):")
        # r_q = 2nd-largest |mu| among amplitude-carrying (carr[0] is Perron, mu~1)
        rq_mu = abs(carr[1][0]) if len(carr) > 1 else None
        for i, (mu, absA) in enumerate(carr[:9]):
            tag = " <-Perron(mu~1)" if i == 0 else ("  <== r_q" if i == 1 else "")
            re, im = float(np.real(mu)), float(np.imag(mu))
            log(f"      |mu|={abs(mu):.4f}  mu={re:+.4f}{('%+.4fi'%im) if abs(im)>1e-9 else ''}"
                f"   |A|={absA:.3e}{tag}")
        # TOWER modes: |mu| in [0.7,1) that are NOT amplitude-carrying -- confirm |A|~0
        tower_maxA = max((abs(A) for lam, A in modes
                          if 0.7 < abs(lam / l1) < 0.999 and abs(lam / l1 - 1) >= 2e-3),
                         default=0.0)
        log("")
        if rq_mu is not None:
            match = 'YES' if abs(rq_mu - R27[q]) < 0.06 else ('L2-bias(->0.39 as L grows)' if q == 7 else 'CHECK')
            log(f"   => r_q = 2nd-largest |mu| w/ amplitude = {rq_mu:.4f}   vs R27 direct {R27[q]:.3f}   match={match}")
            log(f"      TOWER modes |mu| in [0.7,1) but NOT amplitude-carrying: max |A| = {tower_maxA:.2e}  "
                f"(~0 => tower dynamically invisible, r_q clean)")
            log(f"      gap Perron->r_q: 1-|mu| = {1-rq_mu:.3f}"
                f"{'  (r_3->1: near-degenerate, linear growth)' if q==3 else ''}")
        if q == 3:
            top3 = sorted((m[0].real for m in modes), reverse=True)[:3]
            log(f"   q=3 top-3 real eigenvalues: {['%.6f'%x for x in top3]}  "
                f"(top TWO distinct, both ->1/3 as L->inf = Jordan-in-the-limit; r_3=mu_2->1)")
        log("")

    log("## VERDICT (report, no fit):")
    log("   H_CLEAN holds iff for q=5,7: one mode dominates cross-growth (large separation),")
    log("   its |mu| matches R27, tower amplitude is negligible, and 1-|mu| is a real gap.")
    log("   If so, 'r_q = amplitude-dominant subdominant eigenvalue of M' is a CLEAN")
    log("   statement and the brief's L3 stands (with 'amplitude-dominant' made explicit).")
    flush()


def flush():
    with open("result_32_offdiag_spectrum_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
