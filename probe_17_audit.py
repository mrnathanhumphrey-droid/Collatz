"""
PROBE 17 -- ADVERSARIAL AUDIT of the 2026-07-15 qx+1 session (R6-R16).

Not a discovery probe. The session produced an unusual run of confirmations; this tries
to BREAK them. Each check is designed so that a PASS is boring and a FAIL invalidates
named banked results.

THE FOUR REAL VULNERABILITIES (named before testing):

A1. SINGLE POINT OF FAILURE: every ratio in R6-R16 traces to probe_6.stationary(), a
    float power iteration with a fixed iteration cap and a 1e-15 L1 stopping rule.
    *** I NEVER CHECKED WHETHER IT CONVERGED. *** If it silently hits the cap, every
    number today is wrong. Test: measure the actual stationarity residual ||K^T pi - pi||_1
    and the probability-simplex error |sum(pi) - 1| at every (q,k) used.
    KILLS IF FAILS: R6, R7, R8, R9, R10, R11, R13, R15, R16 (all of them).

A2. UNTESTED EXTRAPOLATION: R14 verified the triangular grading ONLY at k=3. R15's
    ratio_within(k) USES it at k=4..8, and that is what produces the headline
    "within is frozen at 0.71958983896 while cross grows at 7/15".
    Test: verify the grading at k=4 (v1 mod d, v2 mod dq, v3 mod dq^2, v4 free).
    KILLS IF FAILS: R15's within/cross SPLIT (note: NOT the 7/15 slope -- see A5).

A3. NEW CODE PATH: R16 introduced a dense chain builder never cross-checked against the
    sparse one it replaced.
    Test: dense vs sparse ||pi||^2 at every (q,k) both can reach.
    KILLS IF FAILS: R16 (the dichotomy + the Aitken limits).

A4. THIN EXACT VERIFICATION: exact rational arithmetic was run at exactly ONE point
    (q=41, k=1, in R10). Everything else is float.
    Test: exact-rational ||pi_k||^2 at fresh (q,k), compared to float.
    KILLS IF FAILS: everything float-derived.

A5. ROBUSTNESS OF THE HEADLINE (the one claim I most want to be true, so test it hardest):
    R15's "cross grows at exactly 7/15" -- does it depend on the `within` formula at all?
    Since within is CONSTANT in k for k>=4, cross differences == total differences, so
    the slope should be independent of within. Test it by recomputing the slope from
    RAW TOTALS ONLY, never touching ratio_within.
    If the slope survives, the 7/15 result is immune to A2.

No new claims are made here. A clean audit changes nothing; a dirty one retracts.
"""
import sys
from math import gcd
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp

from probe_6_conservation_generalize import order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def build_sparse(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    cp = np.array([r for r in range(N) if gcd(r, q) == 1], dtype=np.int64)
    n = len(cp)
    Z = (2 ** M - 1) / 2 ** M
    inv_idx = np.full(N, -1, dtype=np.int64)
    inv_idx[cp] = np.arange(n)
    base_t = (q * cp + 1) % N
    rows, cols, vals = [], [], []
    src = np.arange(n)
    i2v = 1
    for v in range(1, M + 1):
        i2v = (i2v * inv2) % N
        rows.append(src); cols.append(inv_idx[(base_t * i2v) % N])
        vals.append(np.full(n, (0.5 ** v) / Z))
    K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n))
    return K, cp, N, M


def stationary_audited(q, k):
    """Same iteration as probe_6.stationary, but RETURNS THE RESIDUAL."""
    K, cp, N, M = build_sparse(q, k)
    Kt = K.T.tocsr()
    n = len(cp)
    pi = np.full(n, 1.0 / n)
    iters = 0
    for i in range(4000):
        nxt = Kt.dot(pi); s = nxt.sum()
        if s == 0:
            break
        nxt /= s
        iters = i + 1
        if np.abs(nxt - pi).sum() < 1e-15:
            pi = nxt; break
        pi = nxt
    resid = float(np.abs(Kt.dot(pi) - pi).sum())
    simplex = abs(float(pi.sum()) - 1.0)
    rowsum = float(np.abs(np.asarray(K.sum(axis=1)).ravel() - 1.0).max())
    return pi, resid, simplex, rowsum, iters, M, N


def main():
    log("# PROBE 17 -- ADVERSARIAL AUDIT of R6-R16. A PASS is boring; a FAIL retracts.")
    log("")

    # ================= A1 =================
    log("## A1 -- DID THE POWER ITERATION ACTUALLY CONVERGE? (never checked all session)")
    log("   resid = ||K^T pi - pi||_1 ; simplex = |sum(pi) - 1| ; rowsum = max|K.rowsum - 1|")
    log("   KILLS IF FAILS: R6, R7, R8, R9, R10, R11, R13, R15, R16 -- i.e. everything.")
    log("")
    log(f"{'q':>4} {'k':>3} {'n':>7} {'iters':>6} {'residual':>12} {'simplex err':>12} {'rowsum err':>12} {'':>6}")
    worst_resid = 0.0
    cases = [(3, 2), (3, 4), (3, 6), (3, 8), (5, 2), (5, 4), (7, 2), (7, 3),
             (11, 2), (11, 3), (13, 2), (17, 2), (31, 2), (41, 2), (47, 2), (59, 2)]
    for q, k in cases:
        try:
            pi, resid, simp, rs, iters, M, N = stationary_audited(q, k)
        except Exception as e:
            log(f"{q:>4} {k:>3}  ERROR {e}")
            continue
        worst_resid = max(worst_resid, resid)
        flag = "OK" if (resid < 1e-12 and simp < 1e-12 and rs < 1e-12) else "*** SUSPECT ***"
        log(f"{q:>4} {k:>3} {len(pi):>7} {iters:>6} {resid:>12.3e} {simp:>12.3e} {rs:>12.3e} {flag:>6}")
    log("")
    log(f"   worst residual over all audited cases: {worst_resid:.3e}")
    log(f"   A1: {'PASS -- the iteration converged everywhere' if worst_resid < 1e-12 else '*** FAIL -- RESULTS SUSPECT ***'}")
    log("")

    # ================= A2 =================
    log("## A2 -- R14's GRADING AT k=4 (R14 tested only k=3; R15 USES it at k=4..8)")
    log("   predict: v1 mod d | v2 mod d*q | v3 mod d*q^2 | v4 free")
    log("   KILLS IF FAILS: R15's within/cross SPLIT (the 7/15 SLOPE is tested separately in A5)")
    log("")
    for q in [3, 5, 7]:
        d = order_of_two(q)
        N = q ** 4
        M = order_of_two(N)
        inv2 = pow(2, -1, N)

        def val4(v1, v2, v3, v4):
            return (pow(inv2, v4, N)
                    + q * pow(inv2, v3 + v4, N)
                    + q * q * pow(inv2, v2 + v3 + v4, N)
                    + q ** 3 * pow(inv2, v1 + v2 + v3 + v4, N)) % N

        step = max(1, M // 7)
        others = [(a, b, c) for a in range(1, M + 1, step)
                  for b in range(1, M + 1, step) for c in range(1, M + 1, step)]
        bad = [0, 0, 0]
        tot = [0, 0, 0]
        for (b, c, e) in others:
            for v1 in range(1, M + 1 - d):
                tot[0] += 1
                if val4(v1 + d, b, c, e) != val4(v1, b, c, e):
                    bad[0] += 1
            for v2 in range(1, M + 1 - d * q):
                tot[1] += 1
                if val4(b, v2 + d * q, c, e) != val4(b, v2, c, e):
                    bad[1] += 1
            for v3 in range(1, M + 1 - d * q * q):
                tot[2] += 1
                if val4(b, c, v3 + d * q * q, e) != val4(b, c, v3, e):
                    bad[2] += 1
        log(f"   q={q:>3} d={d:>3} M=ord_{{q^4}}(2)={M:>6} (want d*q^3={d*q**3:>6}, "
            f"{'OK' if M == d*q**3 else 'WIEFERICH'})")
        log(f"      v1 mod d      : {tot[0]:>8} checks, {bad[0]:>3} bad -> {'EXACT' if not bad[0] else 'REFUTED'}")
        log(f"      v2 mod d*q    : {tot[1]:>8} checks, {bad[1]:>3} bad -> {'EXACT' if not bad[1] else 'REFUTED'}")
        log(f"      v3 mod d*q^2  : {tot[2]:>8} checks, {bad[2]:>3} bad -> {'EXACT' if not bad[2] else 'REFUTED'}")
    log("")

    # ================= A3 =================
    log("## A3 -- R16's NEW DENSE PATH vs the SPARSE path it replaced")
    log("   KILLS IF FAILS: R16 (the dichotomy + the Aitken limits)")
    log("")
    from probe_16_cross_convergence import total_ratio_dense
    log(f"{'q':>4} {'k':>3} {'sparse ratio':>18} {'dense ratio':>18} {'|rel diff|':>12}")
    a3_worst = 0.0
    for q, k in [(3, 4), (3, 6), (5, 3), (5, 4), (7, 3), (11, 3)]:
        pi, resid, simp, rs, iters, M, N = stationary_audited(q, k)
        nrm = float(np.dot(pi, pi))
        diag = float(sum_p2_exact(M) ** k)
        sparse_r = nrm / diag - 1
        dense_r, gb, _ = total_ratio_dense(q, k)
        rel = abs(sparse_r - dense_r) / abs(sparse_r)
        a3_worst = max(a3_worst, rel)
        log(f"{q:>4} {k:>3} {sparse_r:>18.12f} {dense_r:>18.12f} {rel:>12.2e}")
    log("")
    log(f"   A3: {'PASS -- dense == sparse' if a3_worst < 1e-9 else '*** FAIL ***'}  (worst {a3_worst:.2e})")
    log("")

    # ================= A4 =================
    log("## A4 -- EXACT RATIONAL ||pi_k||^2 vs float (exact ran at ONE point all session)")
    log("   KILLS IF FAILS: everything float-derived")
    log("")
    for q, k in [(5, 2), (7, 2), (3, 3)]:
        N = q ** k
        M = order_of_two(N)
        inv2 = pow(2, -1, N)
        cp = [r for r in range(N) if gcd(r, q) == 1]
        n = len(cp)
        idx = {r: i for i, r in enumerate(cp)}
        Z = F(2 ** M - 1, 2 ** M)
        K = [[F(0)] * n for _ in range(n)]
        for i, r in enumerate(cp):
            t0 = (q * r + 1) % N
            i2v = 1
            for v in range(1, M + 1):
                i2v = (i2v * inv2) % N
                K[i][idx[(t0 * i2v) % N]] += F(1, 2 ** v) / Z
        A = [[K[j][i] - (F(1) if i == j else F(0)) for j in range(n)] for i in range(n)]
        A[n - 1] = [F(1)] * n
        b = [F(0)] * n; b[n - 1] = F(1)
        for c in range(n):
            p = next((r for r in range(c, n) if A[r][c] != 0), None)
            if p is None:
                continue
            A[c], A[p] = A[p], A[c]; b[c], b[p] = b[p], b[c]
            pv = A[c][c]
            A[c] = [x / pv for x in A[c]]; b[c] = b[c] / pv
            for r in range(n):
                if r != c and A[r][c] != 0:
                    f = A[r][c]
                    A[r] = [A[r][j] - f * A[c][j] for j in range(n)]
                    b[r] = b[r] - f * b[c]
        exact_nrm = sum(x * x for x in b)
        pi, *_ = stationary_audited(q, k)
        fl = float(np.dot(pi, pi))
        rel = abs(float(exact_nrm) - fl) / fl
        log(f"   q={q:>3} k={k}: exact ||pi||^2 = {float(exact_nrm):.16f}   float = {fl:.16f}"
            f"   |rel| = {rel:.2e}  {'PASS' if rel < 1e-12 else '*** FAIL ***'}")
    log("")

    # ================= A5 =================
    log("## A5 -- IS THE 7/15 HEADLINE IMMUNE TO A2? Recompute the slope from RAW TOTALS,")
    log("   never touching ratio_within. (within is constant in k, so it should cancel.)")
    log("")
    tots = []
    for k in range(2, 9):
        pi, resid, simp, rs, iters, M, N = stationary_audited(3, k)
        nrm = float(np.dot(pi, pi))
        diag = float(sum_p2_exact(M) ** k)
        tots.append(nrm / diag - 1)
    D = [tots[i + 1] - tots[i] for i in range(len(tots) - 1)]
    tail = D[2:]
    mean = sum(tail) / len(tail)
    log(f"   raw total ratios k=2..8: {['%.8f' % t for t in tots]}")
    log(f"   differences            : {['%.6f' % d for d in D]}")
    log(f"   mean of last {len(tail)} diffs   : {mean:.6f}   vs 7/15 = {7/15:.6f}   "
            f"off {abs(mean-7/15)/(7/15):.2%}")
    log(f"   A5: {'PASS -- the 7/15 slope uses NO within formula and survives' if abs(mean-7/15)/(7/15) < 0.02 else '*** FAIL ***'}")
    flush()


def flush():
    with open("result_17_audit_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
