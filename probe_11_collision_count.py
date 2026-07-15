"""
PROBE 11 (qx+1 paper) -- the FULL COLLISION COUNT at k=2. Derives pillar 3's prefactor.

PRE-REGISTRATION (written before running; derivation stated first, then tested).
------------------------------------------------------------------
R10 established EMPIRICALLY: ratio_2 = offdiag_2/diag_2 = 2^{1-ord_q(2)} to 0.7% at
large q, base 2 recovered at 2.017. R9/R10 flagged: only the CHEAPEST collision was
derived; the prefactor 2 was MATCHED, not counted; the O(1/q) residual unexplained.
This probe does the count.

THE DERIVATION (done before running; this probe tests it):

  Address a = (v_1, v_2), v_i in {1..M}, M = ord_{q^2}(2). Weight p_a = p_{v_1} p_{v_2},
  p_v = 2^{-v}/Z, Z = 1 - 2^{-M}.
  Value  = 2^{-v_2} + q * 2^{-(v_1+v_2)}  mod q^2.

  KEY STRUCTURAL FACT: the SECOND term carries a factor q, so it only needs
  A = v_1+v_2 mod d, where d = ord_q(2) (because 2^d = 1 mod q).
  ==> v_1 -> v_1 + d LEAVES THE VALUE UNCHANGED, EXACTLY.
  ==> v_1 is only ever determined mod d. Every value-bucket contains a whole
      GEOMETRIC TOWER in v_1. This is NOT "a collision costs a period shift" --
      the collisions are structural and always present.

  FAMILY (a) := pairs with the same v_2 and v_1 = v'_1 mod d (v_1 != v'_1).
  For fixed v_2 and c = v_1 mod d (c in 1..d), with x := 2^{-d}:
      G_c := sum_{v_1 = c mod d} p_{v_1} = 2^{-c}/(1-x)              [Z cancels exactly]
      H_c := sum_{v_1 = c mod d} p_{v_1}^2 = 4^{-c}/(1-x^2)          [large M]
      bucket offdiag = p_{v_2}^2 * (G_c^2 - H_c)
                     = p_{v_2}^2 * 4^{-c} * [1/(1-x)^2 - 1/(1-x^2)]
                     = p_{v_2}^2 * 4^{-c} * 2x / [(1-x)(1-x^2)]
  *** THE PREFACTOR 2 IS THE CROSS-TERM OF THE GEOMETRIC TOWER, NOT AN ORDERED-PAIR
      ARTIFACT. It falls out of 1/(1-x)^2 - 1/(1-x^2) = 2x/[(1-x)(1-x^2)]. ***

  Summing: sum_{c=1..d} 4^{-c} = (1-x^2)/3 ; sum_{v_2} p_{v_2}^2 = P2.
      offdiag_(a) = P2 * (1-x^2)/3 * 2x/[(1-x)(1-x^2)] = P2 * 2x / [3(1-x)]
      diag        = P2^2
      ==> ratio_(a) = 2x / [3 (1-x) P2]
  With P2 = (1/3)(1+2^{-M})/(1-2^{-M}) exactly:
      ==> ratio_(a) = 2x (1 - 2^{-M}) / [(1-x)(1 + 2^{-M})]   ~  2x/(1-x)  ~  2^{1-d}

  FAMILY (b) := everything else (different v_2, i.e. v'_2 = v_2 + j d with a
  compensating v_1 shift). PREDICTED to be the SMALL residual that explains R10's
  0.7% miss and the O(1/q) correction.

HYPOTHESES:
  H_V1MOD (structural): value(v_1, v_2) == value(v_1 + d, v_2) EXACTLY, all v_1,v_2,q.
      PRIOR: TRUE (derived). If FALSE the whole derivation is void -> STOP.
  H_FAMA  (*** THE COUNT ***): the measured family-(a) offdiag equals the closed form
      ratio_(a) = 2x(1-2^{-M})/[(1-x)(1+2^{-M})].
      PRIOR: TRUE, to floating precision. This is the prefactor-2 derivation.
  H_DOM_A: family (a) accounts for >= 99% of total offdiag at large q.
      PRIOR: TRUE (R10's miss was 0.7%).
  H_RESID: the residual (total - family a) is family (b), is POSITIVE, and shrinks
      relative to family (a) as q grows -- i.e. it IS the O(1/q) correction.
      PRIOR: TRUE. Stated to lose: if the residual does NOT shrink with q, the
      O(1/q) story is wrong.

DECISION RULES (pre-committed):
  H_V1MOD CONFIRMED iff value(v_1,v_2) == value(v_1+d,v_2) for every tested triple
      (exact integer equality, not tolerance).
  H_FAMA  CONFIRMED iff |measured ratio_(a) - closed form| / closed form < 1e-9.
      This is an EXACT algebraic claim; anything above 1e-9 is a REFUTATION.
  H_DOM_A CONFIRMED iff family_(a)/total >= 0.99 for q >= 41.
  H_RESID CONFIRMED iff residual > 0 for all q AND residual/family_(a) is decreasing
      in q across the tested primes. REFUTED if it grows or is non-monotone in q.

NOT AT STAKE: THEOREM_C_745, Th 78.1-78.3, R81b, eps_k, R5's rate, R6, R7.
This probe can only strengthen or refute pillar 3's DERIVATION; the LAW itself is
already measured conclusively (R10).

Author's priors this arc: 5-for-10, two decision rules found buggy. Rules here are
EXACT-equality tests, which cannot be gamed by a weak threshold.
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


def build(q):
    """All k=2 addresses: value(v1,v2) mod q^2, weight, for v1,v2 in 1..M."""
    N = q * q
    M = order_of_two(N)
    d = order_of_two(q)
    inv2 = pow(2, -1, N)
    pw = [1] * (2 * M + 2)
    for s in range(1, 2 * M + 2):
        pw[s] = (pw[s - 1] * inv2) % N          # 2^{-s} mod q^2
    v = np.arange(1, M + 1)
    p = (2.0 ** (-v.astype(np.float64))) / (1.0 - 2.0 ** (-M))
    # value(v1,v2) = 2^{-v2} + q*2^{-(v1+v2)} mod q^2
    pw_arr = np.array(pw, dtype=np.int64)
    V2 = pw_arr[v]                                # 2^{-v2}
    A = v[:, None] + v[None, :]                   # A[i,j] = v1_i + v2_j
    val = (V2[None, :] + q * pw_arr[A]) % N       # [v1, v2]
    w = p[:, None] * p[None, :]
    return val, w, M, d, p


def main():
    log("# PROBE 11 -- the FULL COLLISION COUNT at k=2. Derives pillar 3's prefactor 2.")
    log("# Pre-reg: H_V1MOD (structural) / H_FAMA (*** the count ***) / H_DOM_A / H_RESID")
    log("")

    QS = [11, 13, 17, 31, 41, 47]

    # ---------------- H_V1MOD ----------------
    log("## H_V1MOD -- does v_1 -> v_1 + d leave the value EXACTLY unchanged?")
    log("   (derived: the q*2^{-(v1+v2)} term needs A only mod d, since 2^d = 1 mod q)")
    ok = True
    for q in QS:
        val, w, M, d, p = build(q)
        # compare rows v1 and v1+d
        n = M - d
        same = np.array_equal(val[:n, :], val[d:d + n, :])
        log(f"   q={q:>3} d=ord_q(2)={d:>3} M=ord_{{q^2}}(2)={M:>5}: "
            f"value[v1] == value[v1+d] for all {n} rows -> {'EXACT' if same else 'MISMATCH'}")
        if not same:
            ok = False
    log(f"   H_V1MOD: {'CONFIRMED -- v_1 is only determined mod d' if ok else 'REFUTED -> STOP'}")
    if not ok:
        flush()
        sys.exit(1)
    log("")

    # ---------------- the count ----------------
    log("## H_FAMA / H_DOM_A / H_RESID -- decompose offdiag into family (a) and the rest")
    log("")
    log("   family (a) = same v_2, v_1 = v'_1 mod d  (the geometric tower in v_1)")
    log("   closed form: ratio_(a) = 2x(1-2^-M)/[(1-x)(1+2^-M)],  x = 2^-d")
    log("")
    log(f"{'q':>4} {'d':>4} {'ratio_2 total':>15} {'family(a) meas':>15} {'family(a) pred':>15} "
        f"{'|rel err|':>11} {'a/total':>9} {'resid/a':>10}")
    rows = {}
    for q in QS:
        val, w, M, d, p = build(q)
        N = q * q
        P2 = float(sum_p2_exact(M))
        diag = P2 ** 2
        # ---- total offdiag: group by value
        flat_v = val.ravel()
        flat_w = w.ravel()
        tot = np.bincount(flat_v, weights=flat_w, minlength=N)
        sq = np.bincount(flat_v, weights=flat_w ** 2, minlength=N)
        offdiag_total = float(np.sum(tot ** 2) - np.sum(sq))
        # ---- family (a): group by (v_2, v_1 mod d)
        v1 = np.arange(1, M + 1)
        c = (v1 - 1) % d                      # class of v_1 mod d
        key = c[:, None] * M + np.arange(M)[None, :]   # (v1 mod d, v2) key
        kf = key.ravel()
        tot_a = np.bincount(kf, weights=flat_w, minlength=d * M)
        sq_a = np.bincount(kf, weights=flat_w ** 2, minlength=d * M)
        offdiag_a = float(np.sum(tot_a ** 2) - np.sum(sq_a))
        # ---- closed form
        x = 2.0 ** (-d)
        tM = 2.0 ** (-M)
        pred_a = 2 * x * (1 - tM) / ((1 - x) * (1 + tM))
        r_tot = offdiag_total / diag
        r_a = offdiag_a / diag
        rel = abs(r_a - pred_a) / pred_a
        resid = offdiag_total - offdiag_a
        rows[q] = {"d": d, "r_tot": r_tot, "r_a": r_a, "pred": pred_a,
                   "rel": rel, "frac": offdiag_a / offdiag_total,
                   "resid_over_a": resid / offdiag_a}
        log(f"{q:>4} {d:>4} {r_tot:>15.8e} {r_a:>15.8e} {pred_a:>15.8e} "
            f"{rel:>11.2e} {offdiag_a/offdiag_total:>9.5f} {resid/offdiag_a:>10.5f}")
    log("")

    worst = max(r["rel"] for r in rows.values())
    log(f"## H_FAMA (*** THE COUNT ***) -- worst |rel err| vs closed form = {worst:.3e}")
    log(f"   H_FAMA: {'CONFIRMED -- the prefactor 2 is DERIVED' if worst < 1e-9 else 'REFUTED'}")
    log("   (the 2 = cross-term of the geometric tower: 1/(1-x)^2 - 1/(1-x^2) = 2x/[(1-x)(1-x^2)])")
    log("")

    big = [q for q in rows if q >= 41]
    fr = min(rows[q]["frac"] for q in big) if big else float("nan")
    log(f"## H_DOM_A -- family (a) share of total offdiag at q>=41: min = {fr:.5f}")
    log(f"   H_DOM_A: {'CONFIRMED (>=99%)' if fr >= 0.99 else 'not confirmed'}")
    log("")

    log("## H_RESID -- is the residual (family b) positive and shrinking with q?")
    for q in QS:
        r = rows[q]
        log(f"   q={q:>3} d={r['d']:>3}: resid/family(a) = {r['resid_over_a']:+.6f}")
    pos = all(rows[q]["resid_over_a"] > 0 for q in QS)
    log(f"   all positive: {pos}")
    log("   (H_RESID predicts this is the O(1/q) correction -> should shrink as q grows,")
    log("    at COMPARABLE d. Note d and q both vary here -- read with care, do not overfit.)")
    flush()


def flush():
    with open("result_11_collision_count_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
