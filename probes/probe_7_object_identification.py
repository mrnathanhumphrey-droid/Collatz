"""
PROBE 7 (qx+1 paper) -- pin the object identification: WHICH S_k does the
(q/3)^k rate belong to, and is the paper's central object R76's object?

PRE-REGISTRATION (written before running; priors stated to lose).
------------------------------------------------------------------
R6 flagged: R76's S_n and the q-sweep's S_k^(q) share a name, and the
identification off q=3 was never verified. Reading R75/R74 makes the tension sharp:

  R74 identity : S_{k+1} = 3^{k+1} * ||d_{k+1}||^2
  R74 deviation: ||d_{k+1}||^2 := sum_{r'} (pi_{k+1}(r') - pi_k(parent(r'))/3)^2
  R75 Plancherel: S_k = sum_{xi: 3 not| xi} |mu_hat_k(xi)|^2
  => R76's S_k is built on the LEVEL-INCREMENTAL DEVIATION d_k.

  q-sweep (probe_5_universal_rate.X_gen): X_k := q^k * ||pi_k||^2
  => built on the RAW STATIONARY MASS pi_k.

These are different functionals. The memory/writeups assert BOTH
"S_k^(q) = q^k||pi_k||^2" AND "S_k^(3) -> 7/15". But R5 also established
||pi_k||^2 ~ 3^{-k} (participation ratio 3^k), which forces
    X_k(q=3) = 3^k * 3^{-k} -> 1,  NOT 7/15.
Both cannot be true of one object. Something is misnamed.

HYPOTHESES:
  H_ID_A: X_k == M_k(1).  PRIOR: FALSE (mu_hat(0)=1 alone puts >=1 into X_k
          but 0 into the coprime-restricted sum). Stated to lose.
  H_ID_B: M_k(1) == q^k * ||d_k||^2  (R74 identity generalized off q=3).
          PRIOR: TRUE. This is the load-bearing one -- if it holds, R76's
          machinery has a well-defined general-q object.
  H_ID_C: at q=3, M_k(1) -> 7/15 and X_k -> 1.  PRIOR: TRUE. If so, the
          "S_k^(3) -> 7/15" claim refers to M_k(1) (=R76's S_k), NOT to X_k,
          and the memory/writeup conflates them.
  H_ID_D: BOTH X_k and M_k(1) carry rate (q/3)^k.
          PRIOR: TRUE (R6's raw numbers already hint: M(1) at q=5 went
          1.367->2.267, ratio 1.658 vs q/3=1.667; q=7 4.256->9.942, ratio
          2.336 vs q/3=2.333). If TRUE, the rate result survives the
          misnaming and only the CONSTANT is object-dependent.

WHY THIS MATTERS: if H_ID_C holds, then R5's headline "S_k^(q) ~ (q/3)^k with
S_k^(3) -> 7/15" is describing two different functionals in one sentence. The
rate claim would still stand (per H_ID_D) but the paper's object needs naming
discipline before publication, and R6's "R76 route blocked" verdict needs to be
read against the RIGHT object.

DECISION RULES (pre-committed):
  H_ID_A CONFIRMED iff |X_k - M_k(1)| < 1e-10 for all tested. Else REFUTED.
  H_ID_B CONFIRMED iff |M_k(1) - q^k||d_k||^2| < 1e-9 for all tested. Else REFUTED.
  H_ID_C CONFIRMED iff |M_k(1) - 7/15| shrinks monotonically in k at q=3 AND
         |X_k - 1| < 1e-6 at q=3, k>=3.
  H_ID_D CONFIRMED iff ratio_k := (obj_k/obj_{k-1})/(q/3) -> 1 (within 5e-2 at
         the largest affordable k) for BOTH objects, all tested q.
  Any object whose ratio does NOT approach 1 -> report as such; no verdict spin.

Not at stake: THEOREM_C_745, Th 78.1-78.3, R81b, eps_k. Naming/identification
only -- this probe cannot falsify the rate mechanism, only relabel its subject.

Reuses probe_6_conservation_generalize.stationary / M_all (same chain as
probe_5_universal_rate.X_gen).
"""
import sys
from math import gcd
import numpy as np

from probe_6_conservation_generalize import stationary, M_all

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def X_of(pi, q, k):
    """q-sweep object: X_k = q^k * ||pi_k||^2 (probe_5 X_gen definition)."""
    return (q ** k) * float(np.dot(pi, pi))


def d_norm2(q, k, pi_k, cp_k, pi_km1, cp_km1):
    """R74 deviation generalized: ||d_k||^2 = sum_{r'} (pi_k(r') - pi_{k-1}(parent(r'))/q)^2
    with parent(r') = r' mod q^{k-1}."""
    Nkm1 = q ** (k - 1)
    idx_km1 = {int(r): i for i, r in enumerate(cp_km1)}
    tot = 0.0
    for i, rp in enumerate(cp_k):
        par = int(rp) % Nkm1
        ppar = float(pi_km1[idx_km1[par]]) if par in idx_km1 else 0.0
        tot += (float(pi_k[i]) - ppar / q) ** 2
    return tot


def main():
    log("# PROBE 7 -- object identification: X_k = q^k||pi||^2  vs  M_k(1) = R76's S_k")
    log("# Pre-registered H_ID_A(prior FALSE) / H_ID_B(TRUE) / H_ID_C(TRUE) / H_ID_D(TRUE)")
    log("")

    cache = {}

    def get(q, k):
        if (q, k) not in cache:
            pi, cp, N = stationary(q, k)
            Mv, _ = M_all(pi, cp, N)
            cache[(q, k)] = (pi, cp, N, Mv)
        return cache[(q, k)]

    log("## Table: both objects side by side")
    log("")
    log(f"{'q':>3} {'k':>2} {'X_k=q^k|pi|^2':>14} {'M_k(1)':>14} {'q^k||d_k||^2':>14} "
        f"{'X ratio/(q/3)':>14} {'M ratio/(q/3)':>14}")
    rows = {}
    for q in [3, 5, 7, 11]:
        kmax = 4 if q <= 5 else (3 if q == 7 else 2)
        prevX = prevM = None
        for k in range(1, kmax + 1):
            try:
                pi, cp, N, Mv = get(q, k)
            except Exception as e:
                log(f"{q:>3} {k:>2}  skip ({e})")
                continue
            X = X_of(pi, q, k)
            M1 = Mv[1].real
            if k >= 2:
                pim1, cpm1, _, _ = get(q, k - 1)
                dn = (q ** k) * d_norm2(q, k, pi, cp, pim1, cpm1)
            else:
                dn = float("nan")
            rX = (X / prevX) / (q / 3) if prevX else float("nan")
            rM = (M1 / prevM) / (q / 3) if prevM else float("nan")
            rows[(q, k)] = (X, M1, dn, rX, rM)
            log(f"{q:>3} {k:>2} {X:>14.8f} {M1:>14.8f} {dn:>14.8f} {rX:>14.5f} {rM:>14.5f}")
            prevX, prevM = X, M1
        log("")

    # ---- H_ID_A ----
    log("## H_ID_A -- is X_k == M_k(1)?  (prior: FALSE)")
    worst = max(abs(v[0] - v[1]) for v in rows.values())
    log(f"   max |X_k - M_k(1)| over all tested = {worst:.6e}")
    log(f"   H_ID_A: {'CONFIRMED' if worst < 1e-10 else 'REFUTED -- DIFFERENT OBJECTS'}")
    log("")

    # ---- H_ID_B ----
    log("## H_ID_B -- is M_k(1) == q^k*||d_k||^2 ?  (R74 identity generalized; prior: TRUE)")
    diffs = [(q, k, abs(v[1] - v[2])) for (q, k), v in rows.items() if not np.isnan(v[2])]
    for q, k, d in diffs:
        log(f"   q={q:2d} k={k}: |M_k(1) - q^k||d_k||^2| = {d:.6e}")
    wb = max(d for _, _, d in diffs) if diffs else float("nan")
    log(f"   worst = {wb:.6e}")
    log(f"   H_ID_B: {'CONFIRMED (R74 identity ports)' if wb < 1e-9 else 'REFUTED (R74 identity does NOT port)'}")
    log("")

    # ---- H_ID_C ----
    log("## H_ID_C -- at q=3: does M_k(1) -> 7/15 while X_k -> 1?  (prior: TRUE)")
    log(f"   7/15 = {7/15:.10f}")
    for k in range(1, 5):
        if (3, k) in rows:
            X, M1, _, _, _ = rows[(3, k)]
            log(f"   k={k}: M_k(1)={M1:.10f}  |M-7/15|={abs(M1-7/15):.3e}   "
                f"X_k={X:.10f}  |X-1|={abs(X-1):.3e}")
    log("")

    # ---- H_ID_D ----
    log("## H_ID_D -- do BOTH objects carry rate (q/3)^k?  (prior: TRUE)")
    log("   (ratio/(q/3) -> 1 means that object has the (q/3)^k rate)")
    for q in [3, 5, 7, 11]:
        ks = sorted(k for (qq, k) in rows if qq == q)
        if len(ks) < 2:
            continue
        kl = ks[-1]
        _, _, _, rX, rM = rows[(q, kl)]
        log(f"   q={q:2d} at k={kl}: X ratio/(q/3) = {rX:.5f}   M ratio/(q/3) = {rM:.5f}")
    log("")
    log("## Raw q=3 sanity: is X_k(q=3) flat at 1, i.e. is (q/3)^k trivial there?")
    for k in range(1, 5):
        if (3, k) in rows:
            log(f"   k={k}: X_k(q=3) = {rows[(3,k)][0]:.10f}")
    flush()


def flush():
    with open("result_7_identification_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
