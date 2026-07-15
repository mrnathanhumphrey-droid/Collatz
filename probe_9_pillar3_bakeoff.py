"""
PROBE 9 (qx+1 paper) -- pillar 3 bake-off: is delta_q ~ a/ord_q(2) or ~ a*2^{-ord_q(2)}?

PRE-REGISTRATION (written before running; falsifier first, prior stated to lose).
------------------------------------------------------------------
Pillar 3 (result_4_ctilde_ord2.md) claims delta_q ~ 0.82/ord_q(2), R^2=0.94, "OOS
validated at q=31,127,73". R8's overlap reframing predicts a DIFFERENT law and gives
it a MECHANISM.

MECHANISM (derived in R8, to be tested here):
  R8 gate proved pi_k IS the q-adic self-similar measure with address (v_1..v_k),
  value = sum_m q^{m-1} 2^{-S_m} mod q^k.
  At k=1: value = 2^{-v_1} mod q with v_1 in {1..ord_q(2)} = exactly one full period
    => the coding is a BIJECTION onto <2>  => ZERO collisions.
    (R8 measured offdiag = 0.00000e+00 at k=1 for q=5,7,11,13. Structural, not luck.)
  At k=2: a collision needs 2^{-v_2} = 2^{-v'_2} mod q, i.e. v'_2 = v_2 + j*ord_q(2).
    Pair weight ~ 2^{-2 v_2} * 2^{-j*ord}. Summing j>=1:
    ==> THE CHEAPEST COLLISION COSTS A FULL PERIOD SHIFT ==> overlap ~ 2^{-ord_q(2)}.
  1/ord_q(2) has NO mechanism -- nothing in the structure produces a reciprocal.

WHY 0.82/ord PLAUSIBLY PASSED (methodological hypothesis, tested here as H_R2TRAP):
  delta is CONVEX in 1/ord and spans 350x. A linear fit is dominated by the two
  largest points; the six small ones cluster near the origin where 1/ord is also
  small. So linear R^2 ~ 0.94 is attainable from a law that is wrong by 55x at q=11
  (0.82/10 = 0.082 predicted vs 0.0015 measured). Linear R^2 CANNOT discriminate
  monotone candidates across a 350x range -- the bake-off must be in LOG space.
  Also: the claimed "OOS" primes q=31,127,73 have ord = 5,7,9 -- INTERPOLATION between
  existing ord values (3,4,8,10,12). They test the fit, never the functional form.

TARGETS (three, all reported):
  ratio_2 := offdiag_2/diag_2  -- the FIRST collision; most direct probe of mechanism.
  ratio_kmax                   -- converged overlap ratio where affordable.
  delta_q(2)                   -- result_4's own target, recomputed from scratch.

HYPOTHESES:
  H_GATE : recomputing result_4's method (c~_q(2) = S_2/(q/3)^2, S_2 = X_2 - X_1,
           X_k = q^k*sum pi^2; delta = c~ - (q-3)/q) reproduces its published deltas
           {q5:+0.092, q7:+0.210, q11:+0.0015, q13:+0.0006, q17:+0.007, q31:+0.059,
            q73:+0.004, q127:+0.015}.  PRIOR: TRUE. If FALSE -> STOP, don't bake off
           against numbers I can't reproduce.
  H_EXP (*** THE TEST ***): log(target) is LINEAR in ord (exponential law), NOT linear
           in log(ord) (power law). Decided by log-space R^2 head-to-head.
           PRIOR: exponential WINS. STATED TO LOSE -- if power-law wins, my mechanism
           is wrong and R8's reframing loses a prediction.
  H_BASE : fitting target = a*c^{-ord} freely, the recovered base c ~ 2.
           PRIOR: c in [1.8, 2.2]. This is the sharp one: the mechanism names the base.
           A recovered c far from 2 = mechanism wrong even if exponential wins.
  H_R2TRAP: LINEAR R^2 of delta vs 1/ord is >= 0.9 (i.e. reproduces result_4's 0.94)
           WHILE the log-space test rejects the same law. PRIOR: TRUE. This is the
           methodological finding, and it is what makes the original claim forgivable.

DECISION RULES (pre-committed):
  H_GATE CONFIRMED iff every recomputed delta matches published to <= 2e-3 absolute
         (published values are given to 1-3 sig figs).
  H_EXP  EXPONENTIAL WINS iff log-space R^2(log t vs ord) - R^2(log t vs log ord)
         > 0.05 on the SAME points, for the ratio_2 target.
         POWER WINS iff the reverse. |diff| <= 0.05 -> INCONCLUSIVE, no verdict.
  H_BASE CONFIRMED iff recovered c in [1.8, 2.2]. Report c with its CI either way.
  H_R2TRAP CONFIRMED iff linear-R^2(delta vs 1/ord) >= 0.9 AND H_EXP says exponential.
  REAL OOS (not result_4's interpolation): fit on MID ord only (5..10), predict the
         EXTREMES ord=3 (q=7) and ord=12 (q=13). Report predicted/actual for both
         laws. This is the test result_4 never ran.

NOT AT STAKE: THEOREM_C_745, Th 78.1-78.3, R81b, eps_k, and R5's RATE (pillar 1).
Pillar 3 only. If H_EXP loses, R8's mechanism is wrong; R7's object identification
and R6's blocked-route verdict are untouched either way.

Author's structural priors this arc: 3-for-8 (H_EQUAL lost, H_ID_C's X_k branch lost,
H_DOM's decision rule was buggy). Stated to lose again.
"""
import sys
from math import gcd, log as ln
import numpy as np

from probe_6_conservation_generalize import stationary, order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


PUBLISHED = {5: 0.092, 7: 0.210, 11: 0.0015, 13: 0.0006,
             17: 0.007, 31: 0.059, 73: 0.004, 127: 0.015}
ORD = {}


def X_and_ratio(q, k):
    """X_k = q^k*||pi_k||^2 ; ratio_k = offdiag/diag with diag = (sum_v p_v^2)^k."""
    N = q ** k
    M = order_of_two(N)
    pi, cp, _ = stationary(q, k)
    nrm = float(np.dot(pi, pi))
    diag = float(sum_p2_exact(M) ** k)
    return (q ** k) * nrm, (nrm - diag) / diag, M


def linfit(x, y):
    """least squares y = a + b x ; returns (a, b, R^2)"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    b, a = np.polyfit(x, y, 1)
    yh = a + b * x
    ss_res = float(np.sum((y - yh) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return a, b, r2


def main():
    log("# PROBE 9 -- pillar 3 bake-off: a/ord  vs  a*2^{-ord}")
    log("# Pre-reg: H_GATE / H_EXP(*** the test ***) / H_BASE / H_R2TRAP. Prior: exponential wins.")
    log("")

    QS = [5, 7, 11, 13, 17, 31, 73, 127]
    for q in QS:
        ORD[q] = order_of_two(q)

    # ---------------- H_GATE ----------------
    log("## H_GATE -- recompute result_4's delta_q(2) from scratch")
    log(f"   method: c~_q(2)=S_2/(q/3)^2, S_2=X_2-X_1, X_k=q^k*sum(pi^2); delta=c~-(q-3)/q")
    log("")
    log(f"{'q':>4} {'ord':>4} {'X_1':>12} {'X_2':>12} {'c~_q(2)':>10} {'(q-3)/q':>9} "
        f"{'delta_mine':>11} {'published':>10} {'|diff|':>9}")
    data = {}
    gate_ok = True
    for q in QS:
        try:
            X1, r1, M1 = X_and_ratio(q, 1)
            X2, r2, M2 = X_and_ratio(q, 2)
        except Exception as e:
            log(f"{q:>4}  skip ({e})")
            continue
        ct = (X2 - X1) / ((q / 3) ** 2)
        base = (q - 3) / q
        dl = ct - base
        pub = PUBLISHED[q]
        diff = abs(dl - pub)
        if diff > 2e-3:
            gate_ok = False
        data[q] = {"delta": dl, "ratio2": r2, "ord": ORD[q], "X1": X1, "X2": X2}
        log(f"{q:>4} {ORD[q]:>4} {X1:>12.6f} {X2:>12.6f} {ct:>10.5f} {base:>9.5f} "
            f"{dl:>+11.5f} {pub:>+10.4f} {diff:>9.2e}")
    log("")
    log(f"   H_GATE: {'CONFIRMED -- reproduces result_4' if gate_ok else 'FAILED -> STOP'}")
    if not gate_ok:
        log("   (Cannot bake off against numbers I cannot reproduce.)")
        flush()
        sys.exit(1)
    log("")
    log("   NOTE: ratio_1 = 0 exactly at k=1 (bijection onto <2>) -- verify:")
    for q in QS:
        _, r1, _ = X_and_ratio(q, 1)
        log(f"      q={q:>3}: ratio_1 = {r1:+.3e}")
    log("")

    # ---------------- the three targets ----------------
    log("## Targets: ratio_2 (first collision), delta_q(2), vs ord")
    log("")
    log(f"{'q':>4} {'ord':>4} {'ratio_2':>13} {'2^-ord':>12} {'ratio_2/2^-ord':>15} "
        f"{'ratio_2*ord':>13} {'delta(2)':>11}")
    for q in QS:
        if q not in data:
            continue
        d = data[q]
        o = d["ord"]
        tp = 2.0 ** (-o)
        log(f"{q:>4} {o:>4} {d['ratio2']:>13.6e} {tp:>12.4e} {d['ratio2']/tp:>15.4f} "
            f"{d['ratio2']*o:>13.5f} {d['delta']:>+11.5f}")
    log("")
    log("   (If ratio_2/2^-ord is ~constant and ratio_2*ord is not => exponential.)")
    log("")

    # ---------------- H_EXP ----------------
    log("## H_EXP (*** THE TEST ***) -- log-space head-to-head, target = ratio_2")
    qs = [q for q in QS if q in data]
    o = np.array([data[q]["ord"] for q in qs], float)
    t = np.array([data[q]["ratio2"] for q in qs], float)
    lt = np.log(t)
    a_e, b_e, r2_e = linfit(o, lt)            # exponential: log t = a + b*ord
    a_p, b_p, r2_p = linfit(np.log(o), lt)    # power law  : log t = a + b*log(ord)
    log(f"   EXPONENTIAL  log(ratio_2) = {a_e:+.4f} {b_e:+.4f}*ord        log-space R^2 = {r2_e:.5f}")
    log(f"   POWER LAW    log(ratio_2) = {a_p:+.4f} {b_p:+.4f}*log(ord)   log-space R^2 = {r2_p:.5f}")
    diff = r2_e - r2_p
    verdict = ("EXPONENTIAL WINS" if diff > 0.05 else
               ("POWER LAW WINS -- mechanism REFUTED" if diff < -0.05 else "INCONCLUSIVE"))
    log(f"   R^2 difference = {diff:+.5f}  =>  {verdict}")
    log("")

    # ---------------- H_BASE ----------------
    log("## H_BASE -- fit ratio_2 = a*c^{-ord} freely; does the base c come out ~2?")
    c_rec = float(np.exp(-b_e))
    a_rec = float(np.exp(a_e))
    log(f"   recovered base c = exp(-slope) = {c_rec:.5f}    (mechanism predicts 2)")
    log(f"   recovered prefactor a = {a_rec:.5f}")
    log(f"   power-law exponent b = {b_p:.4f}   (the 1/ord law needs b = -1)")
    log(f"   H_BASE: {'CONFIRMED' if 1.8 <= c_rec <= 2.2 else 'REFUTED (base is not 2)'}")
    log("")

    # ---------------- H_R2TRAP ----------------
    log("## H_R2TRAP -- does LINEAR R^2 of delta vs 1/ord reproduce result_4's 0.94?")
    dl = np.array([data[q]["delta"] for q in qs], float)
    inv = 1.0 / o
    a_l, b_l, r2_lin = linfit(inv, dl)
    log(f"   LINEAR fit  delta = {a_l:+.5f} {b_l:+.5f}*(1/ord)   linear R^2 = {r2_lin:.5f}")
    log(f"   (result_4 reported ~0.94 and slope ~0.82)")
    log(f"   Worst absolute miss of the pure 0.82/ord law:")
    for q in qs:
        pred = 0.82 / data[q]["ord"]
        act = data[q]["delta"]
        log(f"      q={q:>3} ord={data[q]['ord']:>2}: 0.82/ord = {pred:.5f}   actual = {act:.5f}"
            f"   ratio = {pred/act if act else float('nan'):8.1f}x")
    trap = (r2_lin >= 0.9) and (diff > 0.05)
    log(f"   H_R2TRAP: {'CONFIRMED -- high linear R^2 coexists with a rejected law' if trap else 'not confirmed'}")
    log("")

    # ---------------- REAL OOS ----------------
    log("## REAL OOS -- fit on MID ord (5..10), predict the EXTREMES ord=3 (q=7), ord=12 (q=13)")
    log("   (result_4's 'OOS' q=31,127,73 have ord 5,7,9 = INTERPOLATION. This is the test it never ran.)")
    tr = [q for q in qs if 5 <= data[q]["ord"] <= 10]
    te = [q for q in qs if data[q]["ord"] in (3, 12)]
    log(f"   train q={tr} (ord {[data[q]['ord'] for q in tr]})")
    log(f"   test  q={te} (ord {[data[q]['ord'] for q in te]})")
    o_tr = np.array([data[q]["ord"] for q in tr], float)
    t_tr = np.array([data[q]["ratio2"] for q in tr], float)
    ae, be, _ = linfit(o_tr, np.log(t_tr))
    ap, bp, _ = linfit(np.log(o_tr), np.log(t_tr))
    log("")
    log(f"   {'q':>4} {'ord':>4} {'actual':>13} {'EXP pred':>13} {'x off':>8} "
        f"{'POWER pred':>13} {'x off':>8}")
    for q in te:
        oo = data[q]["ord"]; act = data[q]["ratio2"]
        pe = float(np.exp(ae + be * oo))
        pp = float(np.exp(ap + bp * ln(oo)))
        log(f"   {q:>4} {oo:>4} {act:>13.5e} {pe:>13.5e} {pe/act:>8.2f} "
            f"{pp:>13.5e} {pp/act:>8.2f}")
    log("")
    log("   (x off = predicted/actual; 1.00 is perfect. This extrapolates BOTH directions.)")
    flush()


def flush():
    with open("result_9_bakeoff_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
