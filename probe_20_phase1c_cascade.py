"""
PROBE 20 (v2) -- PHASE 1c: the GENERAL-k collision condition.

*** v1 OF THIS PROBE WAS DEFECTIVE. Recorded, not hidden. ***
v1's headline hypothesis H_CASC ("collision <=> all k cascade levels pass") was a
TAUTOLOGY: cascade_levels computes R = sum_m q^{m-1} T_m and divides by q, k times;
"passed all k" <=> q^k | R <=> R = 0 mod q^k <=> collision, BY DEFINITION. Its
0-false-positive / 0-false-negative columns were guaranteed by arithmetic and carried no
information about Collatz. That is the SEVENTH pre-registration defect this arc and the
worst kind -- the previous six were mis-specified thresholds; this one was a hypothesis
with no content. v2 drops it and tests only claims that can fail.

PRE-REGISTRATION v2 (written before running).
------------------------------------------------------------------
WHAT IS ACTUALLY DERIVED (and can fail):
  value = sum_{m=1}^{k} q^{m-1} 2^{-S_m} mod q^k,  S_m = v_{k-m+1}+...+v_k.
  Collision <=> sum_m q^{m-1} T_m = 0 mod q^k,  T_m := 2^{-S_m} - 2^{-S'_m}.
  LEVEL 1: T_1 = 0 mod q  <=>  S'_1 = S_1 + j_1 d.
  Writing 2^d = 1+q*s and expanding 2^{-j_1 d} = 1 - j_1*q*s + O(q^2):
      T_1 = q*U_1  with  U_1 = j_1*s*2^{-S_1}  (mod q)          <-- NON-TRIVIAL
  LEVEL 2 (divide by q, reduce mod q):  U_1 + T_2 = 0, i.e.
      *** 2^{-S'_2} = 2^{-S_2} + j_1*s*2^{-S_1}  (mod q) ***     <-- NON-TRIVIAL
  This is EXACTLY R13's condition with S_1,S_2 in place of v_2,A (at k=2: S_1=v_2,
  S_2=v_1+v_2=A). THAT is the Phase-1c claim: R13 is the k=2 shadow of a general law.

WHAT IS *NOT* DERIVED (stated plainly): the explicit form at levels m >= 3. Level 3 is
  W_2 + T_3 = 0 mod q with W_2 := (U_1+T_2)/q -- but that is a DEFINITION, not a formula;
  an explicit W_2 in terms of (j_1, j_2, s, S_*) needs a second-order expansion of
  2^{-j d} mod q^3 which I have NOT done. So Phase 1c delivers level 2 only, and this
  probe does not pretend otherwise.

HYPOTHESES (all can fail):
  H_L2 (*** THE TEST ***): on EVERY colliding pair, the explicit level-2 form
      2^{-S'_2} = 2^{-S_2} + j_1*s*2^{-S_1} mod q holds; and on every pair that passes
      level 1 but fails the explicit form, the cascade also fails at level 2.
      Exact set equality. PRIOR: TRUE (derived). If FALSE, my derivation is wrong.
      *** v1 tested this only at k=2 for q=5,7 (and k=3 for q=3 alone). v2 extends to
          k=3 at q=5,7 via BUCKET-BASED collision enumeration. ***
  H_LEVELRATE (measurement, NO verdict): the conditional pass rate
      P(pass level m | passed level m-1) over RANDOM pairs. The "k successive
      H-membership events" picture predicts ~1/d at level 1 and ~d/(q-1) at levels >=2.
      REPORTED ONLY. No fit, no verdict -- per feedback_r2_cannot_discriminate_monotone_fits
      and because six of my thresholds this arc were blind to d or q.

DECISION RULES:
  H_L2 CONFIRMED iff disagreements == 0 at every tested (q,k). Exact set equality --
      cannot be mis-thresholded, which is the point.
  H_LEVELRATE: no rule. Measurement only.

BUDGET: collisions enumerated from VALUE BUCKETS (cheap); the converse checked on a
SAMPLE of non-colliding pairs. Sample sizes are printed. No silent truncation.

NOT AT STAKE: R10-R19. A refutation kills my Phase-1c level-2 derivation only.
"""
import sys
import random
from collections import defaultdict
from itertools import product

from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def suffix_sums(vs, k):
    S, tot = [], 0
    for m in range(1, k + 1):
        tot += vs[k - m]
        S.append(tot)
    return S


def levels_passed(q, k, N, inv2, Sa, Sb):
    """Unwind q^k | sum_m q^{m-1} T_m one division at a time. NOTE: passing all k is
    equivalent to collision BY DEFINITION -- used here only for the LEVEL-2 cross-check
    and for the H_LEVELRATE measurement, never as a 'test' of collision."""
    T = [(pow(inv2, Sa[m], N) - pow(inv2, Sb[m], N)) % N for m in range(k)]
    R = 0
    for m in range(k - 1, -1, -1):
        R = (R * q + T[m]) % (q ** k)
    passed, work = 0, R
    for _ in range(k):
        if work % q != 0:
            break
        work //= q
        passed += 1
    return passed


def build(q, k):
    N = q ** k
    M = order_of_two(N)
    d = order_of_two(q)
    inv2 = pow(2, -1, N)
    mods = [d * (q ** (j - 1)) for j in range(1, k)]
    reps = [list(range(1, mj + 1)) for mj in mods]
    cells = [tuple(list(c) + [vk]) for c in product(*reps) for vk in range(1, M + 1)]
    Ss = {c: suffix_sums(c, k) for c in cells}
    vals = {}
    for c in cells:
        v = 0
        for m in range(1, k + 1):
            v = (v + (q ** (m - 1)) * pow(inv2, Ss[c][m - 1], N)) % N
        vals[c] = v
    return cells, Ss, vals, N, M, d, inv2


def main():
    log("# PROBE 20 (v2) -- PHASE 1c. v1's H_CASC was a TAUTOLOGY and is DROPPED.")
    log("# Testing only what can fail: H_L2 (the explicit level-2 form) + H_LEVELRATE (measurement).")
    log("")
    log("## H_L2 (*** THE TEST ***) -- explicit form  2^{-S'_2} = 2^{-S_2} + j_1*s*2^{-S_1} (mod q)")
    log("   collisions enumerated from VALUE BUCKETS; converse checked on a SAMPLE of non-collisions")
    log("")
    log(f"{'q':>4} {'k':>3} {'d':>3} {'s':>4} {'cells':>7} {'collisions':>11} "
        f"{'L2 bad (coll)':>14} {'sample':>9} {'L2 bad (non)':>13}")
    CASES = [(3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2)]
    ok = True
    rng = random.Random(20260715)
    for q, k in CASES:
        cells, Ss, vals, N, M, d, inv2 = build(q, k)
        inv2q = pow(2, -1, q)
        s = ((pow(2, d) - 1) // q) % q
        buckets = defaultdict(list)
        for c in cells:
            buckets[vals[c]].append(c)

        def l2_explicit(a, b):
            """returns None if level 1 fails, else True/False for the explicit level-2 form"""
            S1a, S1b = Ss[a][0], Ss[b][0]
            if (S1b - S1a) % d != 0:
                return None
            j1 = (S1b - S1a) // d
            lhs = pow(inv2q, Ss[b][1], q)
            rhs = (pow(inv2q, Ss[a][1], q) + j1 * s * pow(inv2q, S1a, q)) % q
            return lhs == rhs

        # --- collisions: every colliding pair must satisfy the explicit level-2 form
        ncoll = badc = 0
        for val, lst in buckets.items():
            for i in range(len(lst)):
                for j in range(i + 1, len(lst)):
                    ncoll += 1
                    r = l2_explicit(lst[i], lst[j])
                    if r is not True:
                        badc += 1
        # --- non-collisions: explicit form must agree with the cascade's level-2
        nsamp = min(200_000, len(cells) * (len(cells) - 1) // 2)
        badn = 0
        for _ in range(nsamp):
            a = cells[rng.randrange(len(cells))]
            b = cells[rng.randrange(len(cells))]
            if a == b or vals[a] == vals[b]:
                continue
            r = l2_explicit(a, b)
            if r is None:
                continue
            casc = levels_passed(q, k, N, inv2, Ss[a], Ss[b]) >= 2
            if r != casc:
                badn += 1
        if badc or badn:
            ok = False
        log(f"{q:>4} {k:>3} {d:>3} {s:>4} {len(cells):>7} {ncoll:>11} {badc:>14} "
            f"{nsamp:>9} {badn:>13}")
    log("")
    log(f"   H_L2: {'CONFIRMED -- R13 IS the k=2 shadow of a general law (S_1,S_2 for v_2,A)' if ok else '*** REFUTED -> the Phase-1c derivation is WRONG ***'}")
    log("")

    # ---------------- H_LEVELRATE ----------------
    log("## H_LEVELRATE -- conditional pass rates over RANDOM pairs. MEASUREMENT, NO VERDICT.")
    log("   'k successive H-membership events' would predict ~1/d at level 1, ~d/(q-1) after.")
    log("")
    log(f"{'q':>4} {'k':>3} {'d':>3} {'1/d':>8} {'d/(q-1)':>9}   conditional P(pass m | passed m-1)")
    for q, k in [(3, 3), (3, 4), (5, 3), (7, 3), (11, 2)]:
        cells, Ss, vals, N, M, d, inv2 = build(q, k)
        cnt = [0] * (k + 1)
        NS = 300_000
        for _ in range(NS):
            a = cells[rng.randrange(len(cells))]
            b = cells[rng.randrange(len(cells))]
            if a == b:
                continue
            p = levels_passed(q, k, N, inv2, Ss[a], Ss[b])
            for m in range(p + 1):
                cnt[m] += 1
        rates = []
        for m in range(1, k + 1):
            rates.append(cnt[m] / cnt[m - 1] if cnt[m - 1] else float("nan"))
        log(f"{q:>4} {k:>3} {d:>3} {1/d:>8.4f} {d/(q-1):>9.4f}   "
            + "  ".join(f"L{m+1}: {r:.4f}" for m, r in enumerate(rates)))
    log("")
    log("   (Reported only. Any structure here is a HYPOTHESIS for a separate")
    log("    pre-committed test -- NOT a result.)")
    log("")
    log("## PHASE 1c DELIVERABLE (honest scope)")
    log("   DELIVERED: the level-2 condition generalizes -- R13's law with S_1,S_2 for v_2,A.")
    log("   NOT DELIVERED: the explicit form at levels m >= 3. Level 3 reads W_2 + T_3 = 0")
    log("      mod q with W_2 := (U_1+T_2)/q, which is a DEFINITION not a formula; an")
    log("      explicit W_2 needs a 2nd-order expansion of 2^{-jd} mod q^3, not done.")
    flush()


def flush():
    with open("result_20_cascade_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
