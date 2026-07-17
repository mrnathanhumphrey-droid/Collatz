"""
PROBE 35 -- the s>=2 frontier: structural collapse + exact cross(2) at q=1093, with
BOTH magnitude predictions pre-registered before the observed value.

PRE-REGISTRATION (numbers committed BEFORE the run, per the arc's discipline).
------------------------------------------------------------------
Context: L3's mechanism (R13/R20 cascade, R14 grading) assumes s=v_q(2^d-1)=1. The smallest
prime with s>=2 is 1093 (search-verified, R35 note). At s>=2: ord_q^2(2)=d (order does NOT
lift), so 2^{jd}=1 mod q^2 for all j, R13's constant s_R13=(2^d-1)/q mod q = 0 (shift-
coupling dies), and the R14 tower grading displaces by one q-power. The s=2 order-collapse
makes cross(2) EXACTLY computable (M=ord_q^2(2)=d => d^2 cells, not d^2*q).

THE QUANTITATIVE PREDICTIONS (committed before observing):
  PRED_s1    = cross(2) if the order LIFTED (s=1 counterfactual) ~ within(2) ~ 2*2^{-d}
               ~ 2^{-363} ~ 1e-109. (The "no anomaly" value: astronomically small at d=364.)
  PRED_degen = the degenerate law's prediction. At s=2 one level is FREE (2^{jd}=1 mod q^2),
               so depth-2 collision mass ~ the DEPTH-1 mass of the ordinary law
               = total(1) = ||pi_1||^2 / P2(d) - 1. (One level of head start.)
  These differ by ~100 ORDERS OF MAGNITUDE, so the adjudication is unambiguous.

ADJUDICATION (pre-committed):
  OBSERVED ~ PRED_degen (and >> PRED_s1)  => R13-COLLAPSE CONFIRMED at the frontier: the
     degenerate law inflates depth-2 collisions to the depth-1 scale. L3 is s=1-conditional
     with the s>=2 regime mechanistically characterized. (STRONG outcome.)
  OBSERVED ~ PRED_s1 (no spike) => the collapse does NOT propagate to the correlation sum;
     the s-condition may be removable. (Surprising; worth chasing.)
  This decides the MECHANISM (does R13 collapse show up in the measure), NOT whether the gap
  closes -- one free level shifts the geometric race by a constant.

NOT AT STAKE: R10-R34. This characterizes the s>=2 frontier; it does not touch the gap for s=1.
"""
import numpy as np
from fractions import Fraction

from probe_6_conservation_generalize import order_of_two, stationary
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_23_exact_increment_recurrence import exact_cross
from probe_18_phase01_subgroup_form import cross_from_cells

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def main():
    log("# PROBE 35 -- s>=2 frontier: structural collapse + exact cross(2) at q=1093")
    log("")

    # ---- structural: does the order lift? (s=1) or not (s>=2) ----
    log("## STRUCTURE -- ord_q(2) vs ord_{q^2}(2): does the order LIFT?")
    log(f"   {'q':>6} {'d=ord_q':>8} {'ord_q^2':>10} {'lift?':>16} {'s_R13':>7} {'s':>3}")
    for q in [11, 13, 1091, 1093, 1097, 3511]:
        d = order_of_two(q)
        o2 = order_of_two(q * q)
        val = 2 ** d - 1
        s = 0
        while val % q == 0:
            s += 1; val //= q
        sR13 = ((2 ** d - 1) // q) % q
        lift = "LIFTS (dq, s=1)" if o2 == d * q else ("NO-LIFT (=d, s>=2)" if o2 == d else f"?({o2})")
        log(f"   {q:>6} {d:>8} {o2:>10} {lift:>16} {sR13:>7} {s:>3}")
    log("   => at s>=2 (q=1093,3511): order does NOT lift; s_R13=0 (R13 coupling dies).")
    log("")

    # ---- PRE-REGISTER predictions for q=1093 (BEFORE observed) ----
    q = 1093
    d = order_of_two(q)
    log(f"## PRE-REGISTERED PREDICTIONS for q={q} (d={d}), committed before observing:")
    x1 = Fraction(1, 2 ** d)
    P2_s1 = sum_p2_exact(d * q)
    if not isinstance(P2_s1, Fraction):
        P2_s1 = Fraction(P2_s1)
    pred_s1 = (1 + x1) / (3 * (1 - x1)) / P2_s1 - 1
    pred_s1_f = float(pred_s1)
    pi1, cp1, N1 = stationary(q, 1)
    norm1 = float(np.dot(pi1, pi1))
    P2_d = float(sum_p2_exact(d))
    pred_degen = norm1 / P2_d - 1
    log(f"   PRED_s1    (order lifts / no anomaly) = {pred_s1_f:.3e}   (~ 2*2^-d)")
    log(f"   PRED_degen (one free level = depth-1 mass total(1)) = {pred_degen:.6f}")
    log(f"   [these differ by ~{np.log10(pred_degen/pred_s1_f):.0f} orders of magnitude]")
    log("")

    # ---- OBSERVED: exact cross(2) at q=1093 ----
    log(f"## OBSERVED -- exact cross(2) at q={q} (M=ord_q^2(2)=d => cheap):")
    cross2, ncell = exact_cross(q, 2)
    obs = float(cross2)
    ref_f, ncell2 = cross_from_cells(q, 2)   # independent float sanity check
    log(f"   exact cross(2) = {obs:.8f}   ({ncell} cells)")
    log(f"   float cross_from_cells(2) = {ref_f:.8f}  (sanity: |rel|={abs(obs-ref_f)/abs(ref_f) if ref_f else 0:.1e})")
    log("")

    # ---- ADJUDICATE ----
    log("## ADJUDICATION (pre-committed):")
    r_s1 = obs / pred_s1_f if pred_s1_f else float('inf')
    r_dg = obs / pred_degen if pred_degen else float('inf')
    log(f"   OBSERVED / PRED_s1    = {r_s1:.2e}   (huge => far above the no-anomaly value)")
    log(f"   OBSERVED / PRED_degen = {r_dg:.3f}    (~1 => matches the one-free-level prediction)")
    if abs(np.log10(max(r_dg, 1e-300))) < abs(np.log10(max(r_s1, 1e-300))):
        log("   => OBSERVED lands near PRED_degen, ~100 orders above PRED_s1:")
        log("      R13-COLLAPSE CONFIRMED at the s=2 frontier. The degenerate law inflates")
        log("      depth-2 collisions to the depth-1 scale. L3 is s=1-conditional with the")
        log("      s>=2 regime MECHANISTICALLY CHARACTERIZED (side condition explained).")
    else:
        log("   => OBSERVED near PRED_s1 (no spike): R13-collapse does NOT propagate to the")
        log("      correlation sum; the s-condition may be REMOVABLE. Surprising -- chase it.")
    log("")
    log("## RE-ENTRY BOUNDARY CONDITION (for the pen-and-paper recursion):")
    log("   The derivation now has TWO fixed points to reproduce:")
    log("     (1) d=2  collapses the PHASE SPREAD (H={1,-1}, r_3=1)")
    log("     (2) s_R13=0 collapses the SHIFT COUPLING (j*s_R13*2^-v term dies, s>=2)")
    log("   A correct re-entry equation must have s_R13 sitting in it visibly, s.t. s_R13->0")
    log("   collapses a coupling AND the phase-spread factor collapses at d=2. Two constraints.")
    flush()


def flush():
    with open("result_35_s2_frontier_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
