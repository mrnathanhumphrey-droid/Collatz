"""
PROBE 14 (qx+1 paper) -- does R11's cell collapse GRADE to k=3? The route to Result 1.

PRE-REGISTRATION (written before running; derivation first; exact-equality tests only).
------------------------------------------------------------------
Result 1 (the paper's HEADLINE, M_k ~ (q/3)^k) is still NOT a theorem. R6 killed the
advertised route (R76 conservation: ports for free AND is insufficient; (q-1)/2=1 <=> q=3).
The live route is the OVERLAP COUNT, which worked at k=2 twice (R11 family (a), exact
identity; R13 family (b), exact iff). To reach Result 1 the count must run at every k.
This probe tests whether the k=2 cell structure GRADES.

THE DERIVATION (done before running):
  Iterating the chain and reducing mod q^k kills the r_0 term (R8), leaving
      value(v_1..v_k) = sum_{m=1}^{k} q^{m-1} * 2^{-S_m}  mod q^k,   S_m = v_{k-m+1}+...+v_k
  At k=3 explicitly (verified against the chain by hand):
      value(v1,v2,v3) = 2^{-v3} + q*2^{-(v2+v3)} + q^2*2^{-(v1+v2+v3)}  mod q^3
  The m-th term carries q^{m-1}, so it needs 2^{-S_m} only mod q^{k-m+1}; and 2^{-S} mod
  q^j depends on S only mod ord_{q^j}(2) = d*q^{j-1}  (d := ord_q(2); non-Wieferich).
  Therefore S_m matters only mod d*q^{k-m}, giving a TRIANGULAR grading:
      v_1 appears ONLY in S_k (needed mod q)    ==> v_1 matters only mod d
      v_2 appears in S_{k-1}, S_k               ==> v_2 matters only mod d*q
      ...
      v_k appears in all                        ==> v_k matters mod d*q^{k-1} = M
  At k=2 this reduces to R11's EXACT fact (v_1 mod d). k=3 is the first real test.

HYPOTHESES:
  H_M (gate): M := ord_{q^3}(2) == d*q^2, and ord_{q^2}(2) == d*q (non-Wieferich).
      PRIOR: TRUE. If FALSE for a tested q, that q is Wieferich-like -> report, skip it.
  H_G1 (*** the k=2 fact, must survive ***): value(v1+d, v2, v3) == value(v1,v2,v3) EXACTLY.
      PRIOR: TRUE (derived). If FALSE the grading is wrong -> STOP, route is dead.
  H_G2 (*** the NEW claim ***): value(v1, v2+d*q, v3) == value(v1,v2,v3) EXACTLY.
      PRIOR: TRUE (derived). This is what "grades" means. If FALSE the route needs rework.
  H_SHARP (the grading is TIGHT, not just an upper bound):
      (a) value(v1, v2+d, v3) != value(v1,v2,v3) for SOME (v1,v2,v3)  [d is NOT enough for v2]
      (b) value(v1, v2, v3+d*q) != value(v1,v2,v3) for SOME (v1,v2,v3) [d*q NOT enough for v3]
      PRIOR: TRUE both. If either is FALSE the collapse is STRONGER than derived -- which
      would be GOOD news (fewer cells), and must be reported as such, not buried.
  H_CELLS: the number of distinct cells is exactly d * (d*q) * M = d^3 * q^3.
      PRIOR: TRUE. NOTE (stated up front, not a surprise later): d^3*q^3 EXCEEDS the
      ~q^2*(q-1) ~ q^3 state count by a factor ~d^3, so the collapse ALONE does NOT give
      injectivity. The weights must do the work (effective 3^k). This probe tests the
      STRUCTURE, not the rate. It cannot prove Result 1.

DECISION RULES (pre-committed):
  H_G1/H_G2 CONFIRMED iff EXACT integer equality on EVERY tested triple (exhaustive where
      affordable, else full sweeps of the shifted coordinate over its whole range for a
      sample of the others). Any single mismatch = REFUTED.
  H_SHARP CONFIRMED iff at least one mismatch is exhibited for each of (a),(b).
  H_CELLS CONFIRMED iff the measured distinct-cell count equals d*(d*q)*M exactly.
  No tolerances anywhere -- these are exact integer identities.

NOT AT STAKE: R10's law, R11, R13, R5's rate, R6, R7, R12, THEOREM_C_745, Th 78.1-78.3.
A refutation here kills only MY proposed route to Result 1, not any banked result.

Author's priors this arc: 9-for-16, and FOUR decision rules found mis-specified
(step<1.5 vs linear growth; dR^2 vs free fits; a relative tolerance vs machine-eps noise;
a d-blind |j|=1 threshold). Hence: exact integer equality only, no thresholds to botch.
"""
import sys
from math import gcd
from itertools import product

from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def val3(q, N, inv2, v1, v2, v3):
    """value(v1,v2,v3) = 2^-v3 + q*2^-(v2+v3) + q^2*2^-(v1+v2+v3) mod q^3."""
    return (pow(inv2, v3, N)
            + q * pow(inv2, v2 + v3, N)
            + q * q * pow(inv2, v1 + v2 + v3, N)) % N


def main():
    log("# PROBE 14 -- does R11's cell collapse GRADE to k=3?  (the route to Result 1)")
    log("# Pre-reg: H_M(gate) / H_G1(k=2 fact must survive) / H_G2(*** NEW ***) / H_SHARP / H_CELLS")
    log("# Exact integer equality only. No tolerances.")
    log("")

    QS = [5, 7, 11, 13]

    # ---------------- H_M gate ----------------
    log("## H_M (gate) -- is ord_{q^j}(2) = d*q^{j-1} (non-Wieferich)?")
    good = []
    for q in QS:
        d = order_of_two(q)
        m2 = order_of_two(q * q)
        m3 = order_of_two(q ** 3)
        ok = (m2 == d * q) and (m3 == d * q * q)
        log(f"   q={q:>3}: d=ord_q(2)={d:>3}  ord_{{q^2}}={m2:>5} (want {d*q:>5})  "
            f"ord_{{q^3}}={m3:>6} (want {d*q*q:>6})  {'OK' if ok else 'WIEFERICH-LIKE -> skip'}")
        if ok:
            good.append(q)
    log("")

    # ---------------- H_G1 / H_G2 ----------------
    log("## H_G1 (v1 mod d, must survive from k=2) / H_G2 (v2 mod d*q, THE NEW CLAIM)")
    log("   exhaustive where affordable; else full sweeps of the shifted coord over 1..M")
    log("")
    g1_ok = g2_ok = True
    for q in good:
        d = order_of_two(q)
        N = q ** 3
        M = order_of_two(N)
        inv2 = pow(2, -1, N)
        exhaustive = M ** 3 <= 4_000_000
        n1 = n2 = 0
        bad1 = bad2 = 0
        if exhaustive:
            rng = range(1, M + 1)
            for v1, v2, v3 in product(rng, rng, rng):
                base = val3(q, N, inv2, v1, v2, v3)
                if v1 + d <= M:
                    n1 += 1
                    if val3(q, N, inv2, v1 + d, v2, v3) != base:
                        bad1 += 1
                if v2 + d * q <= M:
                    n2 += 1
                    if val3(q, N, inv2, v1, v2 + d * q, v3) != base:
                        bad2 += 1
            mode = f"EXHAUSTIVE ({M**3} triples)"
        else:
            # full sweeps of the shifted coordinate; sample the others deterministically
            import numpy as np
            others = [(a, b) for a in range(1, M + 1, max(1, M // 12))
                      for b in range(1, M + 1, max(1, M // 12))]
            for (v2, v3) in others:
                for v1 in range(1, M + 1 - d):
                    n1 += 1
                    if val3(q, N, inv2, v1 + d, v2, v3) != val3(q, N, inv2, v1, v2, v3):
                        bad1 += 1
            for (v1, v3) in others:
                for v2 in range(1, M + 1 - d * q):
                    n2 += 1
                    if val3(q, N, inv2, v1, v2 + d * q, v3) != val3(q, N, inv2, v1, v2, v3):
                        bad2 += 1
            mode = f"SWEEPS ({len(others)} slices x full range)"
        log(f"   q={q:>3} d={d:>3} M=ord_{{q^3}}(2)={M:>6}  [{mode}]")
        log(f"      H_G1  value(v1+d,  v2, v3) == value: {n1:>9} checks, {bad1:>4} mismatches "
            f"-> {'EXACT' if bad1 == 0 else 'REFUTED'}")
        log(f"      H_G2  value(v1, v2+d*q, v3) == value: {n2:>9} checks, {bad2:>4} mismatches "
            f"-> {'EXACT' if bad2 == 0 else 'REFUTED'}")
        if bad1:
            g1_ok = False
        if bad2:
            g2_ok = False
    log("")
    log(f"   H_G1: {'CONFIRMED -- the k=2 fact survives at k=3' if g1_ok else 'REFUTED -> route dead'}")
    log(f"   H_G2: {'CONFIRMED -- THE COLLAPSE GRADES' if g2_ok else 'REFUTED -> route needs rework'}")
    log("")

    # ---------------- H_SHARP ----------------
    log("## H_SHARP -- is the grading TIGHT? (d must NOT suffice for v2; d*q must NOT suffice for v3)")
    for q in good:
        d = order_of_two(q)
        N = q ** 3
        M = order_of_two(N)
        inv2 = pow(2, -1, N)
        a_bad = b_bad = 0
        a_tot = b_tot = 0
        for v1 in range(1, min(M, 6) + 1):
            for v2 in range(1, min(M, 6) + 1):
                for v3 in range(1, min(M, 6) + 1):
                    base = val3(q, N, inv2, v1, v2, v3)
                    if v2 + d <= M:
                        a_tot += 1
                        if val3(q, N, inv2, v1, v2 + d, v3) != base:
                            a_bad += 1
                    if v3 + d * q <= M:
                        b_tot += 1
                        if val3(q, N, inv2, v1, v2, v3 + d * q) != base:
                            b_bad += 1
        log(f"   q={q:>3}: (a) v2+d changes value in {a_bad}/{a_tot} cases "
            f"-> {'TIGHT (d insufficient for v2)' if a_bad else 'COLLAPSE IS STRONGER THAN DERIVED'}")
        log(f"          (b) v3+d*q changes value in {b_bad}/{b_tot} cases "
            f"-> {'TIGHT (d*q insufficient for v3)' if b_bad else 'COLLAPSE IS STRONGER THAN DERIVED'}")
    log("")

    # ---------------- H_CELLS ----------------
    log("## H_CELLS -- distinct cells == d * (d*q) * M = d^3*q^3 ?")
    log("   (stated up front: d^3*q^3 EXCEEDS the ~q^3 state count by ~d^3, so the collapse")
    log("    ALONE gives no injectivity. The weights must do the work. This tests STRUCTURE.)")
    log("")
    log(f"{'q':>4} {'d':>3} {'M':>7} {'cells d*(dq)*M':>16} {'d^3*q^3':>12} {'states ~q^2(q-1)':>17} {'cells/states':>13}")
    for q in good:
        d = order_of_two(q)
        M = order_of_two(q ** 3)
        cells = d * (d * q) * M
        states = q * q * (q - 1)
        log(f"{q:>4} {d:>3} {M:>7} {cells:>16} {d**3 * q**3:>12} {states:>17} {cells/states:>13.1f}")
    log("")
    log("   (cells == d^3*q^3 iff M == d*q^2, i.e. the H_M gate.)")
    flush()


def flush():
    with open("result_14_graded_collapse_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
