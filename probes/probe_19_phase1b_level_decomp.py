"""
PROBE 19 -- PHASE 1b: the GENERAL-k cross decomposition. Where does the k-dependence live?

PRE-REGISTRATION (written before running). EXPLORATORY -- no law will be fitted.
------------------------------------------------------------------
R18 established: cross(2) IS a weighted subgroup sum (exact). But R18 ALSO established
that Phase 1 as scoped cannot reach the target -- at k=2 there is no k to be bounded in.
THE TARGET is: for each fixed q >= 5, cross(k) is BOUNDED in k, i.e. the level-j
contribution c_j obeys c_j <= C_q * r_q^j with r_q < 1 (and r_3 = 1 -> linear).

R16 MEASURED c_j (as increments of cross(k)) decaying at ~0.6 at q=5:
    q=5 cross(k) = 0.0389, 0.0596, 0.0702, 0.0767, 0.0809  (k=2..6)
    increments   =   0.0207, 0.0106, 0.0066, 0.0042        ratios ~0.51, 0.62, 0.64
    q=3 increments -> 0.4655 constant (R15: slope 7/15, no decay)
*** THAT DECAY IS THE THING TO EXPLAIN. This probe asks WHERE it comes from. ***

WHAT I DO NOT KNOW (stated honestly): I do not have a derivation of the level structure
at general k. R13's condition is k=2-specific (it was derived from
value = 2^{-v_2} + q*2^{-(v_1+v_2)} mod q^2). At general k,
    value = sum_{m=1}^{k} q^{m-1} * 2^{-S_m},  S_m = v_{k-m+1}+...+v_k
and the collision condition is a system, not a single congruence. So this probe MEASURES
the structure rather than testing a derived formula. That is deliberate: my quantitative
priors have lost repeatedly this arc, so I am not guessing the answer first.

THE DECOMPOSITION TO MEASURE:
  Cells are (c_1..c_{k-1}, v_k). Two colliding cells differ in some subset of the k
  coordinates. Classify each colliding PAIR by:
    (a) DIFFSET -- the set of coordinates where the two cells differ
    (b) TOPDIFF -- the LOWEST-indexed coordinate where they differ (v_1 is "deepest"/
        coarsest: it is determined only mod d; v_k is finest)
    (c) j -- the shift (v'_k - v_k)/d in the last coordinate (R13's j at k=2)
  and report the MASS SHARE of each class.

HYPOTHESES (weak, exploratory -- all could lose):
  H_GATE: the classified masses sum to cross(k) from R18's exact cell expression.
      PRIOR: TRUE (bookkeeping). If FALSE the classifier is buggy -> STOP.
  H_NEWLEVEL: at level k, the mass carried by pairs whose DIFFSET includes the NEW
      coordinate (i.e. involves v_k) accounts for the INCREMENT c_k = cross(k)-cross(k-1),
      while pairs confined to the old coordinates reproduce cross(k-1).
      PRIOR: plausible, NOT derived. This is the natural "where does the new mass come
      from" split and if it holds it tells us what c_j IS.
  H_SHALLOW: the mass concentrates on pairs differing in FEW coordinates (|DIFFSET| small).
      PRIOR: plausible (weights are geometric), NOT derived.
  H_Q3: at q=3 the class structure does NOT decay with level, matching the linear growth;
      at q=5 it does. PRIOR: TRUE (this is the phase boundary, seen 6x).

DECISION RULES (pre-committed):
  H_GATE CONFIRMED iff |sum of class masses - cross(k)| / cross(k) < 1e-9.
  H_NEWLEVEL / H_SHALLOW / H_Q3: REPORTED AS MEASUREMENTS. NO VERDICT, NO FIT.
      Per feedback_r2_cannot_discriminate_monotone_fits: I will not fit a law to a
      handful of levels. Any structure seen here becomes a HYPOTHESIS for a separate
      pre-committed test, not a result.
  *** No threshold is used anywhere except the GATE, which is an exact identity. This is
      deliberate: six of my pre-registered thresholds have been mis-specified this arc,
      every one blind to a parameter of the mechanism (d, q, the noise floor, shell size).

BUDGET: cells(k) = d^k * q^{k(k-1)/2}; pairs are the real cost (cells^2 / #values).
  Affordable: q=3 k<=4 (11,664 cells / 54 values), q=5 k<=3 (8,000 / 100), q=7 k<=3.
  Anything larger is skipped and SAID so -- no silent truncation.

NOT AT STAKE: R10, R11, R13, R14, R15, R16, R17, R18's two gates. This is exploratory
structure-finding for Phase 1b.
"""
import sys
from math import gcd
from collections import defaultdict
from itertools import product

from probe_6_conservation_generalize import order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_18_phase01_subgroup_form import cross_from_cells

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def classify(q, k, pair_cap=40_000_000):
    """Enumerate cells, group by value, classify colliding pairs. Returns dicts of mass."""
    N = q ** k
    M = order_of_two(N)
    d = order_of_two(q)
    inv2 = pow(2, -1, N)
    P2 = float(sum_p2_exact(M))
    diag = P2 ** k
    mods = [d * (q ** (j - 1)) for j in range(1, k)]
    G = []
    for mj in mods:
        xj = 2.0 ** (-mj)
        G.append([(2.0 ** (-c)) / (1 - xj) for c in range(1, mj + 1)])
    pv = [(2.0 ** (-v)) / (1 - 2.0 ** (-M)) for v in range(1, M + 1)]
    ncell = 1
    for mj in mods:
        ncell *= mj
    ncell *= M

    buckets = defaultdict(list)   # value -> [(coords tuple, mass)]
    reps = [list(range(1, mj + 1)) for mj in mods]
    for combo in product(*reps):
        w1 = 1.0
        for j, c in enumerate(combo):
            w1 *= G[j][c - 1]
        for i, vk in enumerate(range(1, M + 1)):
            vs = list(combo) + [vk]
            val = 0
            s = 0
            for m in range(1, k + 1):
                s += vs[k - m]
                val = (val + (q ** (m - 1)) * pow(inv2, s, N)) % N
            buckets[val].append((tuple(vs), w1 * pv[i]))

    npairs = sum(len(v) * (len(v) - 1) // 2 for v in buckets.values())
    if npairs > pair_cap:
        return None, ncell, npairs, d, M

    by_diffset = defaultdict(float)
    by_topdiff = defaultdict(float)
    by_ndiff = defaultdict(float)
    by_j = defaultdict(float)
    total = 0.0
    for val, lst in buckets.items():
        n = len(lst)
        for i in range(n):
            for jj in range(i + 1, n):
                (a, wa), (b, wb) = lst[i], lst[jj]
                m2 = 2.0 * wa * wb
                diff = tuple(idx for idx in range(k) if a[idx] != b[idx])
                by_diffset[diff] += m2
                by_topdiff[min(diff) if diff else -1] += m2
                by_ndiff[len(diff)] += m2
                jsh = abs(a[k - 1] - b[k - 1]) // d
                by_j[jsh] += m2
                total += m2
    return {"diffset": by_diffset, "topdiff": by_topdiff, "ndiff": by_ndiff,
            "j": by_j, "total": total / diag, "diag": diag}, ncell, npairs, d, M


def main():
    log("# PROBE 19 -- PHASE 1b: the general-k cross decomposition. Where is the k-dependence?")
    log("# EXPLORATORY. Only the GATE has a threshold (an exact identity). No fits.")
    log("")

    CASES = [(3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2), (7, 3)]
    res = {}

    log("## H_GATE -- do the classified pair masses sum to R18's exact cross(k)?")
    log("")
    log(f"{'q':>4} {'k':>3} {'cells':>8} {'pairs':>12} {'classified':>16} {'cross (R18)':>16} {'|rel|':>10}")
    for q, k in CASES:
        out, ncell, npairs, d, M = classify(q, k)
        if out is None:
            log(f"{q:>4} {k:>3} {ncell:>8} {npairs:>12}   SKIPPED (over pair cap) -- SAID, not silent")
            continue
        cc, _ = cross_from_cells(q, k)
        rel = abs(out["total"] - cc) / abs(cc)
        res[(q, k)] = out
        log(f"{q:>4} {k:>3} {ncell:>8} {npairs:>12} {out['total']:>16.10f} {cc:>16.10f} {rel:>10.2e}")
    log("")
    ok = all(abs(res[c]["total"] - cross_from_cells(*c)[0]) / abs(cross_from_cells(*c)[0]) < 1e-9
             for c in res)
    log(f"   H_GATE: {'CONFIRMED -- classifier is exact' if ok else '*** FAILED -> STOP ***'}")
    if not ok:
        flush(); sys.exit(1)
    log("")

    # ---------- H_SHALLOW ----------
    log("## H_SHALLOW -- mass share by NUMBER of differing coordinates (|DIFFSET|)")
    log("   (coordinate 0 = v_1, the coarsest/deepest; coordinate k-1 = v_k, the finest)")
    log("")
    for q, k in sorted(res):
        o = res[(q, k)]
        tot = sum(o["ndiff"].values())
        parts = " ".join(f"|D|={n}: {o['ndiff'][n]/tot:6.4f}" for n in sorted(o["ndiff"]))
        log(f"   q={q} k={k}:  {parts}")
    log("")

    # ---------- topdiff ----------
    log("## MASS SHARE by LOWEST differing coordinate (which coordinate 'causes' the collision)")
    log("")
    for q, k in sorted(res):
        o = res[(q, k)]
        tot = sum(o["topdiff"].values())
        parts = " ".join(f"v_{c+1}: {o['topdiff'][c]/tot:6.4f}" for c in sorted(o["topdiff"]))
        log(f"   q={q} k={k}:  {parts}")
    log("")

    # ---------- H_NEWLEVEL ----------
    log("## H_NEWLEVEL -- does mass NOT involving the new finest coordinate reproduce cross(k-1)?")
    log("   split: pairs whose DIFFSET includes coordinate k-1 (the new one) vs those that don't")
    log("")
    log(f"{'q':>4} {'k':>3} {'cross(k)':>14} {'no-new-coord':>14} {'has-new-coord':>15} {'cross(k-1)':>14}")
    for q, k in sorted(res):
        if k < 3:
            continue
        o = res[(q, k)]
        newc = k - 1
        with_new = sum(m for ds, m in o["diffset"].items() if newc in ds)
        without = sum(m for ds, m in o["diffset"].items() if newc not in ds)
        prev, _ = cross_from_cells(q, k - 1)
        log(f"{q:>4} {k:>3} {o['total']:>14.8f} {without/o['diag']:>14.8f} "
            f"{with_new/o['diag']:>15.8f} {prev:>14.8f}")
    log("")
    log("   (H_NEWLEVEL would predict no-new-coord ~ cross(k-1). REPORTED, no verdict.)")
    log("")

    # ---------- j ----------
    log("## MASS SHARE by j = |v'_k - v_k|/d  (R13's shift index at k=2)")
    log("")
    for q, k in sorted(res):
        o = res[(q, k)]
        tot = sum(o["j"].values())
        parts = " ".join(f"j={n}: {o['j'][n]/tot:6.4f}" for n in sorted(o["j"])[:5])
        log(f"   q={q} k={k}:  {parts}")
    log("")
    log("## H_Q3 -- q=3 (linear, no decay) vs q=5 (converging). Compare the class structure.")
    log("   REPORTED as measurement. Any structure here is a HYPOTHESIS for a separate")
    log("   pre-committed test -- NOT a result. (feedback_r2_cannot_discriminate_monotone_fits)")
    flush()


def flush():
    with open("result_19_phase1b_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
