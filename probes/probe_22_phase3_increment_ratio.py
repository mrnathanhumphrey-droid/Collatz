"""
PROBE 22 -- PHASE 3 (the reconciliation): the WEIGHTED INCREMENT RATIO rho_k = c_{k+1}/c_k.

PRE-REGISTRATION (written before running; falsifier + pre-committed NUMBER).
------------------------------------------------------------------
THE OBJECT. R19 gave cross(k) = sum_{j<=k} c_j with c_k = cross(k) - cross(k-1) the
new-finest-coordinate mass, computable exactly (stationary route, R18 H_P0 gate:
   cross(k) = ||pi_k||^2 / P2^k  -  1  -  ratio_within(q,k)   ).
The Phase-3 target -- prove ||pi_k||^2 ~ C_q 3^{-k} for q>=5 -- is exactly: c_k -> 0
geometrically, i.e. the increment ratio rho_k := c_{k+1}/c_k obeys rho_k <= r_q < 1
uniformly in k. That is Result 1.

WHY THIS OBJECT AND NOT R20's 1/q. R20 measured an UNWEIGHTED per-level collision rate
~1/q. But 1/q = 1/3 at q=3 would give a CONVERGENT c_k -- and q=3 DIVERGES (R15/R16,
cross ~ (7/15)k linear). So the unweighted rate is NOT the decay that governs cross.
c_k is a WEIGHTED mass; its ratio rho_k is the honest Phase-3 quantity. The two are
different objects and this probe measures the right one.

THE TWO FACTS TO RECONCILE (both banked):
  q=3: cross(k) ~ (7/15)k linear  =>  c_k -> 7/15 constant  =>  rho_k -> 1.
  q>=5: cross(k) converges         =>  c_k -> 0 geometric     =>  rho_k -> r_q < 1.

PRE-COMMITTED PREDICTION (a NUMBER, stated to lose -- quantitative priors ~0-for-7 arc):
  *** rho_q := lim_k c_{k+1}/c_k = 3/q. ***
  This would reconcile everything at once: 3/q = 1 at q=3 (divergence), 3/q < 1 for q>=5
  (convergence), AND it repurposes R16's 3/q -- which FAILED as a raw per-level rate --
  as the correct object: the INCREMENT ratio, distinct from R20's unweighted 1/q.
  Predicted limits:  q=3 -> 1.0   q=5 -> 0.600000   q=7 -> 0.428571   q=11 -> 0.272727.
  NOTE UP FRONT: R16 already shows q=5 successive increments 0.0207,0.0106,0.0066,0.0042,
  i.e. rho ~ 0.512, 0.623, 0.636 -- RISING and already above 0.6. So 3/q may well be too
  low. That is fine: it is committed to be falsified if so.

HYPOTHESES / GATES (exact structure; the number is committed but is not a gate):
  G1 (*** THE REAL CLAIM ***): rho_k STABILIZES in k at each q -- the successive
      |rho_{k+1} - rho_k| shrink. If rho_k does NOT settle, there is no uniform r_q and
      the Phase-3 target statement itself is wrong (the important negative).
      Testable with a real sequence only at q=3 (>=4 ratios) and q=5 (>=2). SAID so.
  G2 (phase boundary): rho_k -> 1 at q=3; rho_k -> r_q < 1 at q=5,7,11.
      PRIOR: TRUE (seen 6x). If FALSE the whole picture is wrong.
  PRED (committed, NOT a gate): rho_q -> 3/q. Reported as hit/miss, no verdict rides on it.

PRECISION GUARD (stated before running): c_k = cross(k)-cross(k-1) is a difference of
floats converging to a limit; once c_k < 1e-9 (abs) the ratio is cancellation noise and
is FLAGGED unreliable, not reported as signal. In the k-ranges below c_k stays > 1e-3.

BUDGET (no silent truncation -- each q's max k is SAID): stationary cost ~ n*M,
n=phi(q^k)=(q-1)q^{k-1}, M=ord_2(q^k). Cap n*M <= 25e6 (all sizes below were run before
by probe_15/18 on this box; nothing heavier).
  q=3 (n*M=4*9^{k-1}): k<=8  -> rho_3..rho_7  (5 ratios)  <- the real G1 test
  q=5 (16*25^{k-1}):   k<=5  -> rho_3, rho_4  (2 ratios)
  q=7 (18*49^{k-1}):   k<=4  -> rho_3          (1 ratio)
  q=11(100*121^{k-1}): k<=3  -> rho_2          (1 ratio)
  q=7,11 give ONE ratio each: a point on rho_q vs 3/q, cannot test stabilization. SAID.

NOT AT STAKE: R10-R21, R5's rate, R6, R7, R12, THEOREM_C_745. This measures the Phase-3
object; a surprise refines the target, it does not unbank anything.
"""
import numpy as np

from probe_6_conservation_generalize import stationary, order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_15_tower_k_count import ratio_within

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def cross(q, k):
    """cross(k) via the stationary route (R18 H_P0 gate). Returns float."""
    N = q ** k
    M = order_of_two(N)
    pi, cp, _ = stationary(q, k)
    nrm = float(np.dot(pi, pi))
    diag = float(sum_p2_exact(M) ** k)
    tot = nrm / diag - 1.0
    return tot - ratio_within(q, k)


def kmax_for(q, cap=25_000_000):
    """Largest k with n*M <= cap. n=(q-1)q^{k-1}, M=ord_2(q^k)."""
    k = 1
    while True:
        n = (q - 1) * q ** (k - 1)
        M = order_of_two(q ** k)
        if n * M > cap:
            return k - 1
        k += 1


PRED = {3: 1.0, 5: 0.6, 7: 3 / 7, 11: 3 / 11}


def main():
    log("# PROBE 22 -- PHASE 3: weighted increment ratio rho_k = c_{k+1}/c_k")
    log("# Pre-reg: G1 (*** rho_k stabilizes ***) / G2 (phase boundary) / PRED committed: rho_q -> 3/q")
    log("")

    seqs = {}
    for q in [3, 5, 7, 11]:
        kmax = kmax_for(q)
        log(f"## q={q}  (d={order_of_two(q)},  k = 1..{kmax},  n*M budget 25e6)")
        cr = {}
        for k in range(1, kmax + 1):
            cr[k] = cross(q, k)
        # c_k = cross(k) - cross(k-1); cross(1) reported (expect ~0 -> c_2 = cross(2))
        cks = {k: cr[k] - cr.get(k - 1, 0.0) for k in range(2, kmax + 1)}
        log(f"   cross(1) = {cr.get(1, float('nan')):+.10f}   (expect ~0: no cross-cell at k=1)")
        log("")
        log(f"   {'k':>3} {'cross(k)':>16} {'c_k = cross(k)-cross(k-1)':>28}")
        for k in range(2, kmax + 1):
            log(f"   {k:>3} {cr[k]:>16.10f} {cks[k]:>28.10f}")
        log("")
        # rho_k = c_{k+1}/c_k
        rhos = {}
        log(f"   {'k':>3} {'rho_k = c_{k+1}/c_k':>22} {'|reliable?':>11}")
        for k in range(2, kmax):
            ck, ck1 = cks[k], cks[k + 1]
            reliable = abs(ck) > 1e-9 and abs(ck1) > 1e-9
            r = ck1 / ck if ck != 0 else float("nan")
            if reliable:
                rhos[k] = r
            log(f"   {k:>3} {r:>22.8f} {'yes' if reliable else 'NO (cancellation floor)':>11}")
        seqs[q] = rhos
        log(f"   PRED 3/q = {PRED[q]:.6f}   (committed limit for rho_q)")
        log("")

    # -------- G1 : does rho_k stabilize? --------
    log("## G1 (*** THE REAL CLAIM ***) -- does rho_k stabilize in k?")
    log("   (real test only where >=3 ratios exist: q=3. q=5 has 2; q=7,11 have 1 -- SAID.)")
    log("")
    for q in [3, 5, 7, 11]:
        ks = sorted(seqs[q])
        rs = [seqs[q][k] for k in ks]
        if len(rs) >= 3:
            deltas = [abs(rs[i + 1] - rs[i]) for i in range(len(rs) - 1)]
            shrinking = all(deltas[i + 1] <= deltas[i] + 1e-12 for i in range(len(deltas) - 1))
            log(f"   q={q}: rho = {['%.5f' % r for r in rs]}")
            log(f"          |delta rho| = {['%.5f' % d for d in deltas]}  "
                f"-> {'SHRINKING (stabilizing)' if shrinking else 'NOT monotone-shrinking'}")
        elif len(rs) >= 1:
            log(f"   q={q}: rho = {['%.5f' % r for r in rs]}  (too few to judge stabilization -- SAID)")
        else:
            log(f"   q={q}: no reliable ratio")
    log("")

    # -------- G2 + PRED : the limit vs 3/q, and vs 1 at q=3 --------
    log("## G2 (phase boundary) + PRED (rho_q -> 3/q?)")
    log("   last reliable rho_k = best available estimate of the limit r_q")
    log("")
    log(f"   {'q':>4} {'last rho_k':>14} {'3/q (pred)':>12} {'|miss|':>10} {'r_q<1?':>8}")
    for q in [3, 5, 7, 11]:
        ks = sorted(seqs[q])
        if not ks:
            log(f"   {q:>4}  no reliable ratio")
            continue
        last = seqs[q][ks[-1]]
        miss = abs(last - PRED[q])
        lt1 = "YES" if last < 1.0 - 1e-6 else ("~1 (q=3 divergence)" if q == 3 else "NO")
        log(f"   {q:>4} {last:>14.6f} {PRED[q]:>12.6f} {miss:>10.4f} {lt1:>8}")
    log("")
    log("   READING:")
    log("   - G2 CONFIRMED iff rho -> 1 at q=3 and rho < 1 at q=5,7,11.")
    log("   - PRED (3/q) is hit iff |miss| is small AND rho has settled; if rho is still")
    log("     climbing (q=5 R16 trend), the last value UNDER-reads the limit and 3/q is a")
    log("     lower bound at best. Reported honestly, no verdict rides on the number.")
    log("   - If G1 holds but PRED misses: r_q<1 EXISTS (Phase-3 target stands), value TBD.")
    log("   - If G1 FAILS at q=3 (rho not settling toward 1): the target statement is wrong.")
    flush()


def flush():
    with open("result_22_phase3_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
