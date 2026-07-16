"""
PROBE 33 -- close Holes 3 (B: M primitivity) and 4 (C: within(k) exactly k-flat?).

PRE-REGISTRATION.
------------------------------------------------------------------
B / HOLE 3: L2 (Perron: lambda_1 simple, strictly dominant) was OBSERVED (R25/R32), not
   guaranteed. Verify the sufficient condition: M is IRREDUCIBLE (single recurrent SCC on
   the support) and APERIODIC (lambda_1 the UNIQUE max-modulus eigenvalue). Then Perron-
   Frobenius GIVES simple strictly-dominant lambda_1. Expect: holds for q>=5; q=3 borderline
   (near-double top eigenvalue = the divergence).
   H_PRIM: for q>=5, M is irreducible + aperiodic. PRIOR: TRUE. If a complex eigenvalue sits
   at |lambda|=lambda_1, periodicity -> report it.

C / HOLE 4: we used "cross increments = total increments" assuming within(k) is EXACTLY
   k-flat for k>=3 (R15 H_FLAT, float only). R23 taught that "flat to float" can hide a
   doubly-exp correction (the 7/15 lesson). Compute within(k) in EXACT arithmetic and read
   the true k-dependence.
   H_WFLAT: within(k) is k-flat only to DOUBLY-EXP precision (~2^{-d q}), NOT exactly; the
   residual is far below r_q's scale so it does NOT affect r_q. PRIOR: TRUE (tower structure).
   If the k-dependence were O(1) or O(1/q) it WOULD contaminate r_q -> report the magnitude.

DECISION: report numbers. B: SCC count, |lambda_1|/|lambda_2|, any |lambda|=lambda_1.
C: exact within(k) and within(k)-within(k-1) for k=2..6, vs the 2^{-dq} scale.

NOT AT STAKE: R10-R32.
"""
from fractions import Fraction
from math import gcd
import numpy as np
import scipy.sparse.linalg as spla
import scipy.sparse.csgraph as csg

from probe_6_conservation_generalize import order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_25_transfer_operator_Aprime import build_M

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def reachable(M, s0):
    """states reachable from s0 following nonzero transitions (M[j,i]!=0 means i->j)."""
    n = M.shape[0]
    Mc = M.tocsc()
    seen = np.zeros(n, bool)
    stack = [s0]
    seen[s0] = True
    while stack:
        i = stack.pop()
        # successors j with M[j,i]!=0 -> column i of M
        col = M.getcol(i).nonzero()[0]
        for j in col:
            if not seen[j]:
                seen[j] = True
                stack.append(j)
    return np.where(seen)[0]


def within_exact(q, k):
    """ratio_within(q,k) as an exact Fraction."""
    d = order_of_two(q)
    M = order_of_two(q ** k)
    P2 = sum_p2_exact(M)                      # exact Fraction
    if not isinstance(P2, Fraction):
        P2 = Fraction(P2)
    prod = Fraction(1)
    for j in range(1, k):
        xj = Fraction(1, 2 ** (d * q ** (j - 1)))
        prod *= (1 + xj) / (3 * (1 - xj))
    return prod / (P2 ** (k - 1)) - 1


def main():
    log("# PROBE 33 -- Hole 3 (M primitivity) + Hole 4 (within exactly k-flat?)")
    log("")

    # ============ B / HOLE 3 ============
    log("## B / HOLE 3 -- is M irreducible + aperiodic (=> Perron simple, strictly dominant)?")
    log("")
    for q, L in [(3, 3), (5, 2), (7, 2)]:
        M, idx, n = build_M(q, L)
        s0 = idx[(1, 1, 0)]
        reach = reachable(M, s0)
        sub = M[reach][:, reach]
        # SCC on the reachable subgraph
        ncomp, labels = csg.connected_components(sub, directed=True, connection='strong')
        # dominant SCC = the one carrying the recurrent dynamics (largest)
        sizes = np.bincount(labels)
        big = sizes.max()
        # eigenvalues (top few) for aperiodicity: unique max modulus?
        kk = min(8, n - 2)
        w = spla.eigs(M, k=kk, which='LM', return_eigenvectors=False)
        mags = sorted((abs(z) for z in w), reverse=True)
        l1, l2 = mags[0], mags[1]
        # count eigenvalues within 1e-4 of l1 (periodicity signature)
        on_circle = sum(1 for z in w if abs(abs(z) - l1) < 1e-4 * l1)
        log(f"   q={q} L={L}: dim={n}, reachable from v0={len(reach)}, "
            f"SCCs(reachable)={ncomp}, largest SCC={big}")
        log(f"      |lambda_1|={l1:.6f}  |lambda_2|={l2:.6f}  ratio={l2/l1:.4f}  "
            f"#eigs at |lambda_1|={on_circle}")
        aperiodic = (on_circle == 1)
        verdict = ("PRIMITIVE (aperiodic; lambda_1 strictly dominant)" if aperiodic
                   else f"NOT strictly dominant ({on_circle} eigs on circle -- q=3 divergence" +
                        (" expected)" if q == 3 else " UNEXPECTED)"))
        log(f"      => {verdict}")
        log("")

    # ============ C / HOLE 4 ============
    log("## C / HOLE 4 -- within(k) exact: is it k-flat, or a doubly-exp tower correction?")
    log("")
    for q in [3, 5, 7, 11]:
        d = order_of_two(q)
        log(f"   q={q} (d={d}, tower scale 2^-dq = 2^-{d*q} = {float(2.0**(-d*q)):.2e})")
        prev = None
        for k in range(2, 7):
            w = within_exact(q, k)
            wf = float(w)
            dd = float(w - prev) if prev is not None else float('nan')
            note = ""
            if prev is not None and k >= 4:
                note = f"  |delta|/2^-dq = {abs(dd)/float(2.0**(-d*q)):.2f}" if float(2.0**(-d*q))>0 else ""
            log(f"      within({k}) = {wf:.16f}   delta = {dd:+.3e}{note}")
            prev = w
        log("")

    log("## READ:")
    log("   B: single recurrent SCC + unique max-modulus lambda_1 (q>=5) => M primitive =>")
    log("      L2's 'lambda_1 simple, strictly dominant' is GUARANTEED, not just observed.")
    log("   C: if within(k)-within(k-1) ~ 2^-dq (doubly-exp) for k>=4, within is k-flat to")
    log("      doubly-exp precision; residual << r_q scale => 'cross incr = total incr' is safe.")
    flush()


def flush():
    with open("result_33_holes_BC_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
