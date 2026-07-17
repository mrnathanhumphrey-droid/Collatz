"""
PROBE 15 (qx+1 paper) -- the k-fold TOWER count. Is the within-cell overlap k-INDEPENDENT,
and is the q=3 domination failure entirely CROSS-cell?

PRE-REGISTRATION (written before running; derivation first; predictions pre-committed).
------------------------------------------------------------------
R14 confirmed the triangular grading EXACTLY: coordinate j (j=1..k-1) matters only mod
m_j = d*q^{j-1}; coordinate k ranges over exactly its modulus M = d*q^{k-1} (no tower).
So each of v_1..v_{k-1} carries a GEOMETRIC TOWER with ratio x_j = 2^{-m_j} = 2^{-d*q^{j-1}}.

THE DERIVATION (done before running -- generalizes R11 coordinate-by-coordinate):
  For coordinate j with modulus m_j (M is a multiple of m_j), exactly as in R11:
      G^{(j)}_c = sum_{v = c mod m_j, 1..M} p_v = 2^{-c}/(1-x_j)        [Z cancels]
      sum_{c=1..m_j} (G^{(j)}_c)^2 = (1-x_j^2)/3 / (1-x_j)^2 = (1+x_j)/(3(1-x_j))
      sum_{c=1..m_j} H^{(j)}_c = P2                                      [all v]
  Cell = (c_1..c_{k-1}, v_k). mass = prod_j G^{(j)}_{c_j} * p_{v_k};
  within-cell sum of squares = prod_j H^{(j)}_{c_j} * p_{v_k}^2. Hence
      sum_cells mass^2      = [prod_j (1+x_j)/(3(1-x_j))] * P2
      sum_cells (within p^2) = P2^{k-1} * P2 = P2^k = diag
  ==> *** ratio_within(k) = [prod_{j=1}^{k-1} (1+x_j)/(3(1-x_j))] / P2^{k-1}  -  1 ***
  With P2 -> 1/3 this is  prod_{j=1}^{k-1} (1+x_j)/(1-x_j) - 1  ~  2*2^{-d}  for EVERY k>=2,
  because x_2 = 2^{-dq} is not merely small but NONEXISTENT at any available precision.
  At k=2 it reduces to R11's exact identity (1+x)/(1-x)*(1-2^-M)/(1+2^-M) - 1. [gate]

PRE-COMMITTED PREDICTIONS (computed BEFORE the run, from the formula above):
  q=3  (d=2, x_1=1/4):   ratio_within = 0.666667 at k=2;  ~0.71963 for ALL k>=3
        (x_2 = 2^-6 = 1/64 contributes; x_3 = 2^-18 and beyond do not)
        MEASURED TOTALS (R8): 1.013(k2) 1.604(k3) 2.069(k4) 2.534(k5) 3.000(k6)
        3.466(k7) 3.932(k8)  -- LINEAR, slope 0.4655 ~ 7/15.
        ==> PREDICT cross-cell = total - within carries ALL the growth:
            ~0.35(k2) 0.88(k3) 1.35(k4) 1.81(k5) 2.28(k6) 2.75(k7) 3.21(k8), slope ~0.4655.
  q=5  (d=4):  ratio_within ~ 0.133333 at every k>=2
  q=7  (d=3):  ratio_within ~ 0.285714 at every k>=2
  q=11 (d=10): ratio_within ~ 0.0019550 at every k>=2
  q=13 (d=12): ratio_within ~ 0.00048840 at every k>=2

HYPOTHESES:
  H_GATE: at k=2, ratio_within equals R11's exact identity to float precision, and equals
      R11's MEASURED family(a) (1.955034e-03 at q=11; 4.884005e-04 at q=13;
      7.843137e-03 at q=17; 6.451613e-02 at q=31). PRIOR: TRUE. If FALSE -> STOP.
  H_FLAT (*** THE CLAIM ***): ratio_within(k) is k-independent to <1e-9 relative for all
      k>=3 at fixed q (the x_2 term is the last one that moves it, and only at small d).
      PRIOR: TRUE.
  H_CROSS_GROWS (*** THE CONSEQUENCE ***): at q=3, cross(k) := total(k) - within(k) grows
      LINEARLY in k with slope ~0.4655 (~7/15), carrying ALL the k-growth.
      PRIOR: TRUE. STATED TO LOSE -- if cross is flat and within grows, my derivation is
      wrong and the route needs rework.
  H_CROSS_SMALL: at q>=5, cross(k)/within(k) stays small and does NOT grow with k.
      PRIOR: TRUE at q>=31 (R11 measured ~0.7% at k=2); UNKNOWN at q=5,7 (R11 measured
      eps = 0.058..0.382 at small q, erratic). Reported, not spun.

DECISION RULES (pre-committed):
  H_GATE CONFIRMED iff |ratio_within(k=2) - R11 measured family(a)|/measured < 1e-9.
  H_FLAT CONFIRMED iff max_k |ratio_within(k)/ratio_within(k=3) - 1| < 1e-9 over k=3..kmax.
  H_CROSS_GROWS CONFIRMED iff at q=3 the successive differences cross(k+1)-cross(k) are
      constant to within 2% over k=4..8 AND that constant is within 2% of 0.4655.
  H_CROSS_SMALL: EXPLORATORY at q=5,7 -- report, no verdict, no fit.

NOT AT STAKE: R10's law, R11, R13, R14's grading, R5's rate, R6, R7, R12, THEOREM_C_745.
This tests MY route to Result 1. A refutation kills the route, not a banked result.

Author's priors this arc: 13-for-20. FOUR mis-specified decision rules so far -- hence
the predictions above are NUMBERS committed before the run, per
feedback_r2_cannot_discriminate_monotone_fits.
"""
import sys
import numpy as np

from probe_6_conservation_generalize import stationary, order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def ratio_within(q, k):
    """[prod_{j=1..k-1} (1+x_j)/(3(1-x_j))] / P2^{k-1} - 1,  x_j = 2^{-d*q^{j-1}}."""
    d = order_of_two(q)
    M = order_of_two(q ** k)
    P2 = float(sum_p2_exact(M))
    prod = 1.0
    for j in range(1, k):
        mj = d * (q ** (j - 1))
        xj = 2.0 ** (-mj)
        prod *= (1 + xj) / (3 * (1 - xj))
    return prod / (P2 ** (k - 1)) - 1


def main():
    log("# PROBE 15 -- the k-fold TOWER count. Is within-cell overlap k-INDEPENDENT?")
    log("# Pre-reg: H_GATE / H_FLAT(*** the claim ***) / H_CROSS_GROWS(*** the consequence ***) / H_CROSS_SMALL")
    log("")

    # ---------------- H_GATE ----------------
    log("## H_GATE -- does ratio_within(k=2) reproduce R11's MEASURED family(a)?")
    R11 = {11: 1.955034213098729e-03, 13: 4.884004884004884e-04,
           17: 7.843137254901960e-03, 31: 6.451612903225806e-02}
    ok = True
    for q, meas in R11.items():
        pred = ratio_within(q, 2)
        rel = abs(pred - meas) / meas
        log(f"   q={q:>3}: ratio_within(2) = {pred:.15e}   R11 measured family(a) = {meas:.15e}"
            f"   |rel| = {rel:.2e}  {'OK' if rel < 1e-9 else 'FAIL'}")
        if rel >= 1e-9:
            ok = False
    log(f"   H_GATE: {'CONFIRMED' if ok else 'FAILED -> STOP'}")
    if not ok:
        flush(); sys.exit(1)
    log("")

    # ---------------- H_FLAT ----------------
    log("## H_FLAT (*** THE CLAIM ***) -- is ratio_within(k) k-INDEPENDENT for k>=3?")
    log("")
    log(f"{'q':>4} {'d':>3} " + "".join(f"{'k='+str(k):>17}" for k in range(2, 8)) + f"{'max dev k>=3':>14}")
    for q in [3, 5, 7, 11, 13]:
        d = order_of_two(q)
        vals = [ratio_within(q, k) for k in range(2, 8)]
        tail = vals[1:]
        dev = max(abs(v / tail[0] - 1) for v in tail)
        log(f"{q:>4} {d:>3} " + "".join(f"{v:>17.10e}" for v in vals) + f"{dev:>14.2e}")
    log("")
    log("   (k=2 differs from k>=3 by the x_2 = 2^-dq term; k>=3 should be flat because")
    log("    x_3 = 2^-dq^2 is nonexistent at any available precision.)")
    log("")

    # ---------------- H_CROSS_GROWS ----------------
    log("## H_CROSS_GROWS (*** THE CONSEQUENCE ***) -- at q=3, does cross = total - within")
    log("   carry ALL the k-growth, linearly, slope ~0.4655 (~7/15)?")
    log("")
    log(f"{'k':>3} {'total (measured)':>18} {'within (derived)':>18} {'cross = t - w':>15} {'d(cross)':>11}")
    prev = None
    diffs = []
    for k in range(2, 9):
        N = 3 ** k
        M = order_of_two(N)
        pi, cp, _ = stationary(3, k)
        nrm = float(np.dot(pi, pi))
        diag = float(sum_p2_exact(M) ** k)
        tot = nrm / diag - 1
        wit = ratio_within(3, k)
        cr = tot - wit
        dd = (cr - prev) if prev is not None else float("nan")
        if prev is not None and k >= 5:
            diffs.append(dd)
        log(f"{k:>3} {tot:>18.8f} {wit:>18.8f} {cr:>15.8f} {dd:>11.5f}")
        prev = cr
    log("")
    if diffs:
        mn, mx, av = min(diffs), max(diffs), sum(diffs) / len(diffs)
        spread = (mx - mn) / av
        # BUGFIX (R17 audit): this previously computed off vs my PREDICTION 0.4655 while
        # the label said 7/15, publishing 0.06% when the true deviation from 7/15 is 0.19%.
        # Report BOTH, each labelled with what it is actually compared against.
        off_pred = abs(av - 0.4655) / 0.4655
        off_715 = abs(av - 7 / 15) / (7 / 15)
        log(f"   cross increments k=5..8: mean {av:.5f}, spread {spread:.2%}")
        log(f"      vs my pre-committed prediction 0.4655 : off {off_pred:.2%}")
        log(f"      vs 7/15 = {7/15:.6f}                    : off {off_715:.2%}")
        off = off_715
        v = (spread < 0.02) and (off < 0.02)
        log(f"   H_CROSS_GROWS: {'CONFIRMED -- ALL k-growth is CROSS-cell' if v else 'NOT confirmed by the rule'}")
    log("")

    # ---------------- H_CROSS_SMALL ----------------
    log("## H_CROSS_SMALL -- at q>=5, does cross/within stay small and NOT grow with k?")
    log("   (EXPLORATORY at q=5,7 per pre-reg: report, no verdict, no fit)")
    log("")
    log(f"{'q':>4} {'k':>3} {'total':>16} {'within':>16} {'cross':>14} {'cross/within':>14}")
    for q in [5, 7, 11, 13]:
        for k in [2, 3, 4, 5]:
            N = q ** k
            try:
                M = order_of_two(N)
                ncp = (q - 1) * q ** (k - 1)
                if ncp * M > 40_000_000:
                    log(f"{q:>4} {k:>3}  skip (n*M={ncp*M:.1e} > budget)")
                    break
                pi, cp, _ = stationary(q, k)
            except Exception as e:
                log(f"{q:>4} {k:>3}  skip ({e})"); break
            nrm = float(np.dot(pi, pi))
            diag = float(sum_p2_exact(M) ** k)
            tot = nrm / diag - 1
            wit = ratio_within(q, k)
            log(f"{q:>4} {k:>3} {tot:>16.8e} {wit:>16.8e} {tot-wit:>14.6e} {(tot-wit)/wit:>14.6f}")
        log("")
    flush()


def flush():
    with open("result_15_tower_k_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
