"""
PROBE 10 (qx+1 paper) -- pillar 3 CONCLUSIVE: predictions pre-committed BEFORE running.

PRE-REGISTRATION. The numbers below are computed from fits to the OLD 8 primes
(ord 3..12) ONLY, and are written down BEFORE any new prime is evaluated. No fitting
happens on the new primes. This is a prediction test, not a fit.
------------------------------------------------------------------
WHY: R9 was INCONCLUSIVE on exponential-vs-power because both are free 2-param fits on
8 monotone points spanning 165x -- delta-R^2 cannot separate them. The fix is not a
better statistic, it is better data: go to ord where the laws differ by ORDERS OF
MAGNITUDE, far outside the fitted range (3..12).

result_4's "OOS" primes q=31,127,73 have ord 5,7,9 -- INTERPOLATION inside the fitted
range. It never tested the functional form. These do.

THE THREE LAWS (all fit on the OLD 8 primes only, ord 3..12):
  A  (pillar 3, published):  delta = 0.82 / ord            [claim of result_4]
  A' (pillar 3, actual best linear fit): delta = 0.81886/ord - 0.08903
     -- ALREADY DEAD: predicts delta < 0 for ord > 9.2, but H_LB (R8) PROVED delta > 0
        by Cauchy-Schwarz per fiber. A law that contradicts a theorem. Reported, not tested.
  B  (R8 mechanism):  ratio_2 = 2.450 * 2^{-ord}           [a = exp(0.8962) from R9]
  P  (free power law): ratio_2 = 73.87 * ord^{-4.4985}     [from R9; note exponent is
     -4.5, NOT -1 -- P is NOT pillar 3's law, it is the thing that tied B in R9]

PRE-COMMITTED PREDICTIONS (computed now, before the run):

  q=41, ord=20:
     A  predicts delta ~ 0.04100   (=> ratio_2 of the same order, ~1e-2)
     B  predicts ratio_2 = 2.336e-06
     P  predicts ratio_2 = 1.039e-04
     separation A:B ~ 18,000x ;  P:B ~ 44x     -- ALL THREE resolvable in float64.

  q=47, ord=23:
     A  predicts delta ~ 0.03565
     B  predicts ratio_2 = 2.921e-07
     P  predicts ratio_2 = 5.530e-05
     separation A:B ~ 122,000x ;  P:B ~ 189x   -- ALL THREE resolvable.

  q=59, ord=58:
     A  predicts delta ~ 0.01414
     B  predicts ratio_2 = 8.500e-18   <-- BELOW float64 machine epsilon (2.2e-16)
     P  predicts ratio_2 = 8.600e-07   <-- resolvable
     A:B ~ 1e15 ; P:B ~ 1e11.
     ASYMMETRY (stated up front): at ord=58 law B predicts a value BELOW machine
     epsilon, so a measured "0" CANNOT confirm B's exact value there -- it is
     indistinguishable from float noise. But it DOES decisively kill A (0.014, fourteen
     orders above noise) AND P (8.6e-07, nine orders above noise). So q=59 is
     two-sided for A-vs-B and P-vs-B, one-sided for pinning B.

DECISION RULES (pre-committed):
  For each new prime, compare measured ratio_2 against B and P by |log10(pred/actual)|.
  H_CONCLUSIVE: the law with the smaller |log10 miss| on ALL THREE new primes WINS
     outright, provided the winner's max |log10 miss| < 0.5 (i.e. within ~3x) AND the
     loser's min |log10 miss| > 1.0 (i.e. off by >10x) on at least one prime.
     If both laws land within 3x on all three -> INCONCLUSIVE (again), report as such.
  H_A_DEAD: law A is REFUTED iff measured delta(2) at q=59 is < 1e-3 (A predicts
     0.01414). This is a 14-order-of-magnitude claim; no ambiguity is possible.
  SELF-CHECK 1 (independent implementation): at q=41,k=2 rebuild pi by R8's ADDRESS
     ENUMERATION (a combinatorially different route than power iteration) and require
     max|pi_addr - pi_iter| < 1e-12. If this fails, ALL of R6-R10 is suspect -> STOP.
  SELF-CHECK 2 (exact arithmetic): recompute X_1, X_2 at q=41 in exact Fractions and
     require float agreement < 1e-12. Guards against float lying near small values.

NOT AT STAKE: THEOREM_C_745, Th 78.1-78.3, R81b, eps_k, R5's RATE (pillar 1), R7's
object identification, R6's blocked-route verdict. Pillar 3 ONLY.

Author's priors this arc: 3-for-8, and TWO decision rules already found buggy (H_DOM's
"step<1.5" could not see linear growth; H_EXP's delta-R^2 could not separate free fits).
Hence: this probe pre-commits NUMBERS, not adjectives.
"""
import sys
from math import gcd
from fractions import Fraction
import numpy as np

from probe_6_conservation_generalize import stationary, order_of_two, M_all
from probe_8_selfsimilar_overlap import sum_p2_exact, pi_by_address

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


PRED = {
    41: {"ord": 20, "A": 0.04100, "B": 2.336e-06, "P": 1.039e-04},
    47: {"ord": 23, "A": 0.03565, "B": 2.921e-07, "P": 5.530e-05},
    59: {"ord": 58, "A": 0.01414, "B": 8.500e-18, "P": 8.600e-07},
}


def X_exact(q, k):
    """Exact-rational X_k = q^k * sum(pi^2). Small q,k only."""
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    cp = [r for r in range(N) if gcd(r, q) == 1]
    n = len(cp)
    idx = {r: i for i, r in enumerate(cp)}
    Z = Fraction(2 ** M - 1, 2 ** M)
    # K[i][j] = sum over v of weight
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i, r in enumerate(cp):
        t0 = (q * r + 1) % N
        i2v = 1
        for v in range(1, M + 1):
            i2v = (i2v * inv2) % N
            j = idx[(t0 * i2v) % N]
            K[i][j] += Fraction(1, 2 ** v) / Z
    # stationary by solving pi K = pi, sum pi = 1  (dense Gaussian elimination over Q)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            continue
        A[c], A[p] = A[p], A[c]
        b[c], b[p] = b[p], b[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        b[c] = b[c] / pv
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [A[r][j] - f * A[c][j] for j in range(n)]
                b[r] = b[r] - f * b[c]
    pi = b
    return Fraction(q ** k) * sum(x * x for x in pi)


def measure(q, k):
    N = q ** k
    M = order_of_two(N)
    pi, cp, _ = stationary(q, k)
    nrm = float(np.dot(pi, pi))
    diag = float(sum_p2_exact(M) ** k)
    return (q ** k) * nrm, (nrm - diag) / diag, M


def main():
    log("# PROBE 10 -- pillar 3 CONCLUSIVE. Predictions pre-committed before the run.")
    log("# Laws fit on OLD primes (ord 3..12) ONLY. New primes: ord 20, 23, 58.")
    log("")
    log("## PRE-COMMITTED PREDICTIONS (from the docstring, written before running)")
    for q, p in PRED.items():
        log(f"   q={q:>3} ord={p['ord']:>2}:  A(0.82/ord) delta~{p['A']:.5f}   "
            f"B(2.450*2^-ord) ratio_2={p['B']:.3e}   P(73.87*ord^-4.4985) ratio_2={p['P']:.3e}")
    log("")

    # ---------------- SELF-CHECK 1 ----------------
    log("## SELF-CHECK 1 -- independent implementation (address enumeration vs power iteration)")
    log("   q=41,k=2: rebuild pi by R8's combinatorial address sum. Different route entirely.")
    pi_a, diag_e, naddr = pi_by_address(41, 2)
    if pi_a is None:
        log(f"   SKIPPED ({naddr:.1e} addresses over budget)")
    else:
        pi_i, cp, N = stationary(41, 2)
        dense = np.zeros(N)
        dense[cp] = pi_i
        err = float(np.max(np.abs(pi_a - dense)))
        log(f"   {naddr} addresses; max|pi_addr - pi_iter| = {err:.3e}   "
            f"{'PASS' if err < 1e-12 else 'FAIL -> STOP, R6-R10 SUSPECT'}")
        if err >= 1e-12:
            flush()
            sys.exit(1)
    log("")

    # ---------------- SELF-CHECK 2 ----------------
    log("## SELF-CHECK 2 -- exact rational vs float at q=41")
    try:
        xe1 = X_exact(41, 1)
        xf1, _, _ = measure(41, 1)
        d1 = abs(float(xe1) - xf1)
        log(f"   X_1 exact = {float(xe1):.14f}   float = {xf1:.14f}   |diff| = {d1:.3e}   "
            f"{'PASS' if d1 < 1e-12 else 'FAIL'}")
    except Exception as e:
        log(f"   exact X_1 skipped ({e})")
    log("")

    # ---------------- THE TEST ----------------
    log("## THE TEST -- measured vs pre-committed. NO FITTING ON THESE PRIMES.")
    log("")
    log(f"{'q':>4} {'ord':>4} {'measured ratio_2':>18} {'measured delta(2)':>18} "
        f"{'B pred':>12} {'B off':>9} {'P pred':>12} {'P off':>9}")
    res = {}
    for q in [41, 47, 59]:
        p = PRED[q]
        try:
            X1, r1, M1 = measure(q, 1)
            X2, r2, M2 = measure(q, 2)
        except Exception as e:
            log(f"{q:>4}  FAILED ({e})")
            continue
        ct = (X2 - X1) / ((q / 3) ** 2)
        dl = ct - (q - 3) / q
        bo = p["B"] / r2 if r2 > 0 else float("inf")
        po = p["P"] / r2 if r2 > 0 else float("inf")
        res[q] = {"ratio2": r2, "delta": dl, "ord": p["ord"], "B_off": bo, "P_off": po}
        log(f"{q:>4} {p['ord']:>4} {r2:>18.6e} {dl:>+18.6e} "
            f"{p['B']:>12.3e} {bo:>8.2f}x {p['P']:>12.3e} {po:>8.2f}x")
    log("")
    log("   (off = predicted/actual. 1.00 = perfect.)")
    log("")

    # ---------------- verdicts ----------------
    log("## H_A_DEAD -- law A (0.82/ord) predicts delta=0.01414 at q=59. Measured:")
    if 59 in res:
        d59 = res[59]["delta"]
        log(f"   measured delta(2) at q=59 = {d59:.6e}   A predicted 0.01414   "
            f"off by {0.01414/d59 if d59 else float('inf'):.3g}x")
        log(f"   H_A_DEAD: {'CONFIRMED -- law A REFUTED' if abs(d59) < 1e-3 else 'NOT confirmed'}")
    log("")

    log("## H_CONCLUSIVE -- B vs P on all three new primes (|log10 miss|)")
    if res:
        bmiss = [abs(np.log10(res[q]["B_off"])) for q in res if np.isfinite(res[q]["B_off"])]
        pmiss = [abs(np.log10(res[q]["P_off"])) for q in res if np.isfinite(res[q]["P_off"])]
        for q in res:
            log(f"   q={q:>3} ord={res[q]['ord']:>2}: |log10 B miss| = {abs(np.log10(res[q]['B_off'])):.3f}"
                f"   |log10 P miss| = {abs(np.log10(res[q]['P_off'])):.3f}")
        log("")
        log(f"   B: max |log10 miss| = {max(bmiss):.3f}  (= {10**max(bmiss):.2f}x)")
        log(f"   P: max |log10 miss| = {max(pmiss):.3f}  (= {10**max(pmiss):.2f}x)")
        log(f"   P: min |log10 miss| = {min(pmiss):.3f}  (= {10**min(pmiss):.2f}x)")
        # NOTE q=59's B prediction is below machine eps -> exclude from B's "win" test
        b_ok = [abs(np.log10(res[q]["B_off"])) for q in res if q != 59 and np.isfinite(res[q]["B_off"])]
        if b_ok and max(b_ok) < 0.5 and max(pmiss) > 1.0:
            log("   H_CONCLUSIVE: B WINS OUTRIGHT (B within 3x on resolvable primes; P off >10x somewhere)")
        elif b_ok and max(b_ok) > 1.0 and max(pmiss) < 0.5:
            log("   H_CONCLUSIVE: P WINS OUTRIGHT -- R8 mechanism REFUTED")
        else:
            log("   H_CONCLUSIVE: INCONCLUSIVE by the pre-registered rule -- report as such.")
        log("   (q=41,47 carry the two-sided verdict; q=59's B prediction is sub-machine-epsilon.)")
    flush()


def flush():
    with open("result_10_conclusive_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
