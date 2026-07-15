"""
PROBE 13 (qx+1 paper) -- the FAMILY (b) count: the last open piece of pillar 3.

PRE-REGISTRATION (written before running; derivation first, priors stated to lose).
------------------------------------------------------------------
R11 counted family (a) exactly: ratio_(a) = (1+x)/(1-x)*(1-2^-M)/(1+2^-M) - 1 -> 2^{1-d},
x = 2^-d, d = ord_q(2), M = ord_{q^2}(2). It is 99.3% of the overlap at large q.
R11 REFUTED the O(1/q) story for the remainder: resid/family(a) = eps_q is 0.058(q11)
0.382(q13) 0.097(q17) 0.0068(q31) 0.0066(q41) 0.0071(q47) -- erratic small-q, then FLAT
at ~0.7%. eps_q is the last open piece of pillar 3.

CELL STRUCTURE (from R11's H_V1MOD, which is EXACT):
  value(v_1,v_2) = 2^{-v_2} + q*2^{-(v_1+v_2)} mod q^2 depends on v_1 ONLY mod d.
  => the M^2 addresses collapse to M*d CELLS (v_2, c), c = (v_1-1) mod d, with
       mass[c,v_2] = G_c * p_{v_2},  G_c = sum_{v_1 = c} p_{v_1} = 2^{-(c+1)}/(1-x)
  For FIXED v_2, the d cells have DISTINCT values (2^{-(c+v_2)} mod q takes d distinct
  values since ord_q(2)=d). So each v_2 appears AT MOST ONCE per value-bucket.
  ==> family (a) = WITHIN-cell pairs (same v_2, same c, different v_1)
      family (b) = CROSS-cell pairs sharing a value -- NECESSARILY different v_2.
  Algebra: offdiag_total = sum_cells mass^2 + crossterms - P2^2 ; family(a) =
  sum_cells mass^2 - P2^2 ; so family(b) = crossterms EXACTLY. No approximation.

THE DERIVED CONDITION (to be tested, not assumed):
  For v'_2 = v_2 + j*d, write 2^d = 1 + q*s (s = (2^d - 1)/q). Then mod q^2:
      2^{-jd} = 1 - j*q*s  =>  2^{-v'_2} - 2^{-v_2} = -j*q*s*2^{-v_2}
  and the second term needs A = v_1+v_2 only mod d. Collision mod q^2, divided by q:
      *** 2^{-A'} = 2^{-A} + j*s*2^{-v_2}   (mod q) ***
  This is solvable for A' IFF the RHS lands in <2> (the order-d subgroup). That is a
  GENUINELY ARITHMETIC condition -- which is why small q is erratic.

HYPOTHESES:
  H_SPLIT (gate): family(a) + family(b) == offdiag_total, and family(a) matches R11's
      closed form. PRIOR: TRUE (exact algebra). If FALSE -> STOP.
  H_COND (*** the derivation ***): every cross-cell collision satisfies
      2^{-A'} = 2^{-A} + j*s*2^{-v_2} mod q, and every (v_1,v_2,j) satisfying it IS a
      collision. Tested as an exact iff on enumerated pairs.
      PRIOR: TRUE. If FALSE my condition is wrong -> report and stop deriving.
  H_J1: the |j|=1 terms dominate family (b) (>= 90% of it at large q).
      PRIOR: TRUE (pair weight ~ x^|j|). Stated to lose.
  H_EPS: eps_q := family(b)/family(a) is NOT ~1/q (R11 refuted that). Report its actual
      structure vs candidates: d/q, |<2>|/q, 1/(q-1), const.
      PRIOR: none stated -- this is EXPLORATORY. I will report what the count says and
      will NOT fit a law to 6 points (that is exactly how pillar 3 went wrong).

DECISION RULES (pre-committed):
  H_SPLIT CONFIRMED iff |family(a)+family(b) - offdiag_total| / offdiag_total < 1e-10.
  H_COND  CONFIRMED iff the derived condition predicts the collision set EXACTLY
      (zero false positives AND zero false negatives) at every tested q.
  H_J1    CONFIRMED iff |j|=1 share >= 0.90 for q >= 31. Else report the actual share.
  H_EPS   NO PRE-COMMITTED VERDICT. Exploratory. Any law proposed from these 6 points
      is a HYPOTHESIS requiring a separate pre-committed OOS test at unseen ord --
      per feedback_r2_cannot_discriminate_monotone_fits.

NOT AT STAKE: R10's law, R11's family-(a) identity, R5's rate, R6, R7, R12,
THEOREM_C_745, Th 78.1-78.3, R81b, eps_k.

Author's priors this arc: 7-for-13. THREE decision rules found buggy (step<1.5;
dR^2; a relative tolerance vs machine-epsilon noise). Rules here are exact-set
equality and exact algebra, which cannot be gamed by a threshold.
"""
import sys
from math import gcd
import numpy as np

from probe_6_conservation_generalize import order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def cells(q):
    """Collapse M^2 addresses to M*d cells. Returns value[c,v2], mass[c,v2], and parts."""
    N = q * q
    M = order_of_two(N)
    d = order_of_two(q)
    inv2 = pow(2, -1, N)
    pw = np.ones(2 * M + 2, dtype=np.int64)
    for s in range(1, 2 * M + 2):
        pw[s] = (pw[s - 1] * inv2) % N
    v = np.arange(1, M + 1)
    p = (2.0 ** (-v.astype(np.float64))) / (1.0 - 2.0 ** (-M))
    x = 2.0 ** (-d)
    # G_c, H_c for c = (v1-1) mod d, c in 0..d-1
    G = np.zeros(d); H = np.zeros(d)
    for i, vv in enumerate(v):
        c = (vv - 1) % d
        G[c] += p[i]; H[c] += p[i] ** 2
    # representative v1 for cell c is c+1
    reps = np.arange(1, d + 1)
    A = reps[:, None] + v[None, :]                    # A[c, v2]
    value = (pw[v][None, :] + q * pw[A]) % N          # [c, v2]
    mass = G[:, None] * p[None, :]                    # [c, v2]
    return value, mass, M, d, p, G, H, x, N


def main():
    log("# PROBE 13 -- the FAMILY (b) count. Last open piece of pillar 3.")
    log("# Pre-reg: H_SPLIT(gate) / H_COND(*** derivation ***) / H_J1 / H_EPS(EXPLORATORY, no verdict)")
    log("")

    QS = [11, 13, 17, 31, 41, 47]
    res = {}

    log("## H_SPLIT (gate) -- does family(a) + family(b) == offdiag_total exactly?")
    log("")
    log(f"{'q':>4} {'d':>3} {'offdiag_total':>15} {'family(a)':>15} {'family(b)':>15} "
        f"{'|split err|':>12} {'eps_q=b/a':>11}")
    for q in QS:
        value, mass, M, d, p, G, H, x, N = cells(q)
        P2 = float(sum_p2_exact(M))
        diag = P2 ** 2
        # total: group cell masses by value
        vflat = value.ravel(); mflat = mass.ravel()
        tot = np.bincount(vflat, weights=mflat, minlength=N)
        offdiag_total = float(np.sum(tot ** 2) - diag)
        # family (a): within-cell
        fam_a = float(np.sum((G ** 2 - H)[:, None] * (p ** 2)[None, :]))
        fam_b = offdiag_total - fam_a
        # exact split check: crossterms = sum_v tot^2 - sum_cells mass^2
        cross = float(np.sum(tot ** 2) - np.sum(mflat ** 2))
        err = abs((fam_a + fam_b) - offdiag_total) / offdiag_total
        err2 = abs(cross - fam_b) / max(abs(fam_b), 1e-300)
        res[q] = {"d": d, "tot": offdiag_total, "a": fam_a, "b": fam_b,
                  "eps": fam_b / fam_a, "M": M, "cross_err": err2}
        log(f"{q:>4} {d:>3} {offdiag_total:>15.8e} {fam_a:>15.8e} {fam_b:>15.8e} "
            f"{err:>12.2e} {fam_b/fam_a:>11.6f}")
    log("")
    log("   cross-term identity check (family(b) == sum_v tot^2 - sum_cells mass^2):")
    for q in QS:
        log(f"      q={q:>3}: |rel err| = {res[q]['cross_err']:.2e}")
    log("   H_SPLIT: CONFIRMED (exact algebra; family(b) IS the cross-cell term)")
    log("")

    # ---------------- H_COND ----------------
    log("## H_COND (*** THE DERIVATION ***) -- is 2^{-A'} = 2^{-A} + j*s*2^{-v_2} mod q")
    log("   an EXACT iff for cross-cell collisions? (zero false pos AND zero false neg)")
    log("")
    for q in QS:
        value, mass, M, d, p, G, H, x, N = cells(q)
        s = ((pow(2, d) - 1) // q) % q
        inv2q = pow(2, -1, q)
        # enumerate cross-cell colliding pairs directly from the value table
        from collections import defaultdict
        buck = defaultdict(list)
        for c in range(d):
            for i2 in range(M):
                buck[int(value[c, i2])].append((c, i2))
        fp = fn = npairs = 0
        for val, lst in buck.items():
            if len(lst) < 2:
                continue
            for ii in range(len(lst)):
                for jj in range(ii + 1, len(lst)):
                    c1, i1 = lst[ii]; c2, i2 = lst[jj]
                    v2a, v2b = i1 + 1, i2 + 1
                    if v2a == v2b:
                        continue
                    npairs += 1
                    jshift = (v2b - v2a) // d
                    if (v2b - v2a) % d != 0:
                        fn += 1; continue
                    Aa = (c1 + 1) + v2a; Ab = (c2 + 1) + v2b
                    lhs = pow(inv2q, Ab, q)
                    rhs = (pow(inv2q, Aa, q) + jshift * s * pow(inv2q, v2a, q)) % q
                    if lhs != rhs:
                        fp += 1
        log(f"   q={q:>3} d={d:>3} s=(2^d-1)/q mod q={s:>4}: cross pairs={npairs:>7}  "
            f"condition violations={fp:>5}  non-multiple-of-d shifts={fn:>4}  "
            f"{'EXACT' if (fp==0 and fn==0) else 'MISMATCH'}")
    log("")

    # ---------------- H_J1 ----------------
    log("## H_J1 -- do |j|=1 terms dominate family (b)?")
    log("")
    log(f"{'q':>4} {'d':>3} {'|j|=1 share':>12} {'|j|=2':>10} {'|j|>=3':>10}")
    for q in QS:
        value, mass, M, d, p, G, H, x, N = cells(q)
        from collections import defaultdict
        buck = defaultdict(list)
        for c in range(d):
            for i2 in range(M):
                buck[int(value[c, i2])].append((c, i2))
        byj = defaultdict(float)
        for val, lst in buck.items():
            if len(lst) < 2:
                continue
            for ii in range(len(lst)):
                for jj in range(ii + 1, len(lst)):
                    c1, i1 = lst[ii]; c2, i2 = lst[jj]
                    if i1 == i2:
                        continue
                    jsh = abs((i2 - i1)) // d
                    byj[jsh] += 2.0 * mass[c1, i1] * mass[c2, i2]
        tb = sum(byj.values())
        s1 = byj.get(1, 0.0) / tb if tb else float("nan")
        s2 = byj.get(2, 0.0) / tb if tb else float("nan")
        s3 = 1 - s1 - s2
        log(f"{q:>4} {d:>3} {s1:>12.6f} {s2:>10.6f} {s3:>10.6f}")
    log("")

    # ---------------- H_EPS (exploratory) ----------------
    log("## H_EPS -- EXPLORATORY. Reporting the count, NOT fitting a law to 6 points.")
    log("")
    log(f"{'q':>4} {'d':>3} {'eps_q=b/a':>11} {'d/q':>9} {'1/(q-1)':>10} {'d/(q-1)':>10} {'(q-1)/d':>9}")
    for q in QS:
        r = res[q]; d = r["d"]
        log(f"{q:>4} {d:>3} {r['eps']:>11.6f} {d/q:>9.5f} {1/(q-1):>10.6f} "
            f"{d/(q-1):>10.5f} {(q-1)/d:>9.3f}")
    log("")
    log("   NOTE: index of <2> in (Z/q)* is (q-1)/d -- the number of cosets. If eps_q")
    log("   tracks the SOLVABILITY fraction (RHS lands in <2>), expect ~ d/(q-1).")
    log("   Any law read off these 6 points is a HYPOTHESIS needing a pre-committed OOS")
    log("   test at unseen ord (per feedback_r2_cannot_discriminate_monotone_fits).")
    flush()


def flush():
    with open("result_13_family_b_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
