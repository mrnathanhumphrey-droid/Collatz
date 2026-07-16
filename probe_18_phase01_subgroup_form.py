"""
PROBE 18 -- PHASE 0 + PHASE 1 of the Result-1 plan: put the cross-cell term in
weighted-subgroup form, and test whether it CONCENTRATES.

PRE-REGISTRATION (written before running; both phases are GATES, not discoveries).
------------------------------------------------------------------
TARGET (restated precisely, per the plan): for each fixed q >= 5, cross(k) is BOUNDED
in k. NOT "eps_q is small", NOT uniform in q. R16 measured the boundedness; we need a
proof, and a proof needs the sum in a standard form.

PHASE 0 -- the general-k cross expression.
  R14's grading (audited EXACT at k=3 and k=4 by R17-A2): coordinate j (j=1..k-1)
  matters only mod m_j = d*q^{j-1}; coordinate k ranges over M = d*q^{k-1}.
  So the M^k addresses collapse to CELLS (c_1..c_{k-1}, v_k), with
      mass[cell] = prod_{j<k} G^{(j)}_{c_j} * p_{v_k},   G^{(j)}_c = 2^{-c}/(1-x_j)
  and (R15, exact):
      cross(k) = [ sum_value (sum_{cells->value} mass)^2 - sum_cells mass^2 ] / diag
  This is EXACT -- no approximation. Gate: it must equal the measured total - within.

PHASE 1 -- reduce k=2 to a WEIGHTED SUBGROUP SUM.
  R13 proved (exact iff): cells (c_1,v_2) and (c'_1,v'_2) collide iff
      v'_2 = v_2 + j*d   AND   2^{-A'} = 2^{-A} + j*s*2^{-v_2} (mod q),
      A = c_1+v_2, A' = c'_1+v'_2, s = (2^d - 1)/q mod q.
  Write H = <2> in F_q^*, |H| = d;  a := 2^{-A} in H,  b := 2^{-v_2} in H. Then:
      *** the condition is  a + j*s*b  in  H  ***
  -- a MULTIPLICATIVE-SUBGROUP ADDITIVE-SHIFT incidence problem in F_q.
  Standard object (|{(a,b) in H^2 : a+sb in H}| ~ d^3/q + Weil error) EXCEPT for one
  twist: our weights are  w(h) ~ 2^{-ind(h)}  -- EXPONENTIAL IN THE DISCRETE LOG.
  So the sum is NOT equidistributed; it should CONCENTRATE on small ind.

HYPOTHESES:
  H_P0 (gate): cross_from_cells(k) == measured total(k) - within(k).
      PRIOR: TRUE (exact algebra). If FALSE -> the cell decomposition is wrong -> STOP.
  H_P1FORM (gate): the subgroup-form sum, built ONLY from R13's condition
      (a + j*s*b in H) and the cell weights, reproduces cross(2) exactly.
      PRIOR: TRUE. If FALSE, my reduction to a subgroup sum is wrong -> the whole
      Phase 2/3 literature plan is aimed at the wrong object -> STOP AND REPLAN.
  H_CONC (*** THE FALSIFIER ***): truncating the sum to v_2 <= T reproduces cross(2)
      with relative error decaying like ~2^{-T}, UNIFORMLY in q.
      PRIOR: TRUE. STATED TO LOSE. If the error does NOT decay geometrically in T, the
      concentration story is wrong, Phase 3a (truncation route) is DEAD before any math
      is invested, and only Phase 3b (Weil) remains -- which cannot work at small d
      (q=7, d=3), so the route would need replanning.
  H_JDOM: |j|=1 carries the sum (R13 measured 1.000000 at q=47). Re-checked here as a
      cross-check on the subgroup form, not a new claim.

DECISION RULES (pre-committed):
  H_P0 CONFIRMED iff |cross_cells - (total - within)| / |cross| < 1e-9 at every tested
      (q,k). Exact algebra -- anything larger is a real disagreement.
  H_P1FORM CONFIRMED iff |cross_subgroup - cross_cells| / |cross_cells| < 1e-9 at every
      tested q.
  H_CONC CONFIRMED iff for every tested q, err(T) := |cross_T - cross_full|/cross_full
      satisfies err(T+5)/err(T) <= 0.10 (i.e. at least ~2^-5 per 5 steps) across the
      range where err(T) > 1e-14 (above float noise).
      REFUTED iff err(T) plateaus above 1e-6 or decays slower than 2^{-T/4}.
      NOTE stated up front: the weights p_v ~ 2^-v are geometric BY CONSTRUCTION, so a
      trivial pass here proves only that the tail is geometric -- it does NOT prove the
      HEAD is boundable, which is the actual math. This gate is necessary, not
      sufficient, and is labelled as such so a pass is not oversold.

NOT AT STAKE: R10, R11, R13, R14, R15, R16, R17, R5's rate, R6, R7, R12,
THEOREM_C_745. A refutation kills the Phase-3a route, not a banked result.

Author's priors this arc: 17-for-26; FIVE mis-specified decision rules; and R17 caught a
published number that was wrong in the flattering direction. Hence: gates are exact
identities, and H_CONC's limitation is stated BEFORE the run.
"""
import sys
from math import gcd
from collections import defaultdict
import numpy as np

from probe_6_conservation_generalize import stationary, order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_15_tower_k_count import ratio_within

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def cell_data(q, k):
    """Cells (c_1..c_{k-1}, v_k). Returns list of (value, mass) plus diag."""
    N = q ** k
    M = order_of_two(N)
    d = order_of_two(q)
    inv2 = pow(2, -1, N)
    P2 = float(sum_p2_exact(M))
    diag = P2 ** k
    mods = [d * (q ** (j - 1)) for j in range(1, k)]      # m_1..m_{k-1}
    G = []
    for mj in mods:
        xj = 2.0 ** (-mj)
        G.append([(2.0 ** (-c)) / (1 - xj) for c in range(1, mj + 1)])
    pv = [(2.0 ** (-v)) / (1 - 2.0 ** (-M)) for v in range(1, M + 1)]
    ncell = 1
    for mj in mods:
        ncell *= mj
    ncell *= M
    if ncell > 12_000_000:
        return None, ncell, M, d, diag
    out = defaultdict(float)
    # iterate cells
    idx = [0] * (k - 1)
    reps = [list(range(1, mj + 1)) for mj in mods]
    from itertools import product
    for combo in product(*reps):
        w1 = 1.0
        for j, c in enumerate(combo):
            w1 *= G[j][c - 1]
        for i, vk in enumerate(range(1, M + 1)):
            # value = sum_{m=1}^{k} q^{m-1} * 2^{-S_m}, S_m = suffix sum of (combo..., vk)
            vs = list(combo) + [vk]
            val = 0
            s = 0
            for m in range(1, k + 1):
                s += vs[k - m]
                val = (val + (q ** (m - 1)) * pow(inv2, s, N)) % N
            out[val] += w1 * pv[i]
    return out, ncell, M, d, diag


def cross_from_cells(q, k):
    out, ncell, M, d, diag = cell_data(q, k)
    if out is None:
        return None, ncell
    tot_sq = sum(m * m for m in out.values())
    # sum_cells mass^2 -- recompute cheaply from the factorized structure
    mods = [d * (q ** (j - 1)) for j in range(1, k)]
    s_mass2 = 1.0
    for mj in mods:
        xj = 2.0 ** (-mj)
        s_mass2 *= sum(((2.0 ** (-c)) / (1 - xj)) ** 2 for c in range(1, mj + 1))
    s_mass2 *= sum(((2.0 ** (-v)) / (1 - 2.0 ** (-M))) ** 2 for v in range(1, M + 1))
    return (tot_sq - s_mass2) / diag, ncell


def main():
    log("# PROBE 18 -- PHASE 0 + PHASE 1: cross in weighted-SUBGROUP form; does it CONCENTRATE?")
    log("# Pre-reg: H_P0(gate) / H_P1FORM(gate) / H_CONC(*** falsifier ***) / H_JDOM")
    log("")

    # ================= PHASE 0 =================
    log("## PHASE 0 / H_P0 -- cross_from_cells(k) == measured total(k) - within(k) ?")
    log("")
    log(f"{'q':>4} {'k':>3} {'cells':>10} {'cross (cells)':>16} {'total - within':>16} {'|rel|':>10}")
    p0_ok = True
    for q, k in [(3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2)]:
        cc, ncell = cross_from_cells(q, k)
        if cc is None:
            log(f"{q:>4} {k:>3} {ncell:>10}  SKIP (over cell budget)")
            continue
        M = order_of_two(q ** k)
        pi, cp, _ = stationary(q, k)
        tot = float(np.dot(pi, pi)) / float(sum_p2_exact(M) ** k) - 1
        meas = tot - ratio_within(q, k)
        rel = abs(cc - meas) / abs(meas)
        if rel >= 1e-9:
            p0_ok = False
        log(f"{q:>4} {k:>3} {ncell:>10} {cc:>16.10f} {meas:>16.10f} {rel:>10.2e}")
    log("")
    log(f"   H_P0: {'CONFIRMED -- the cell decomposition is exact' if p0_ok else '*** FAILED -> STOP ***'}")
    if not p0_ok:
        flush(); sys.exit(1)
    log("")

    # ================= PHASE 1 =================
    log("## PHASE 1 / H_P1FORM -- rebuild cross(2) from R13's SUBGROUP condition alone")
    log("   H = <2> in F_q^*, |H| = d;  a = 2^{-A}, b = 2^{-v_2} in H;  condition: a + j*s*b in H")
    log("")

    def cross_subgroup(q, Tcap=None):
        d = order_of_two(q)
        N = q * q
        M = order_of_two(N)
        inv2N = pow(2, -1, N)
        inv2q = pow(2, -1, q)
        s = ((pow(2, d) - 1) // q) % q
        x = 2.0 ** (-d)
        P2 = float(sum_p2_exact(M))
        diag = P2 ** 2
        G = [(2.0 ** (-c)) / (1 - x) for c in range(1, d + 1)]
        pv = [(2.0 ** (-v)) / (1 - 2.0 ** (-M)) for v in range(1, M + 1)]
        Hset = {pow(2, i, q) for i in range(d)}
        # discrete-log table: for y in H, the c in 1..d with 2^{-c} = y
        dlog = {}
        for c in range(1, d + 1):
            dlog[pow(inv2q, c, q)] = c
        vmax = M if Tcap is None else min(M, Tcap)
        acc = 0.0
        for c1 in range(1, d + 1):
            for v2 in range(1, vmax + 1):
                A = c1 + v2
                a = pow(inv2q, A, q)
                b = pow(inv2q, v2, q)
                for j in range(1, (M - v2) // d + 1):
                    rhs = (a + j * s * b) % q
                    if rhs not in Hset:
                        continue
                    v2p = v2 + j * d
                    if v2p > vmax:
                        continue
                    cA = dlog[rhs]                       # A' = cA mod d
                    # c'_1 = A' - v'_2 mod d, mapped into 1..d
                    c1p = (cA - v2p - 1) % d + 1
                    acc += 2.0 * G[c1 - 1] * pv[v2 - 1] * G[c1p - 1] * pv[v2p - 1]
        return acc / diag

    log(f"{'q':>4} {'d':>3} {'s':>4} {'cross (subgroup)':>18} {'cross (cells)':>16} {'|rel|':>10}")
    p1_ok = True
    fulls = {}
    for q in [5, 7, 11, 13, 17, 31]:
        cs = cross_subgroup(q)
        cc, _ = cross_from_cells(q, 2)
        rel = abs(cs - cc) / abs(cc)
        fulls[q] = cs
        d = order_of_two(q)
        s = ((pow(2, d) - 1) // q) % q
        if rel >= 1e-9:
            p1_ok = False
        log(f"{q:>4} {d:>3} {s:>4} {cs:>18.12f} {cc:>16.12f} {rel:>10.2e}")
    log("")
    log(f"   H_P1FORM: {'CONFIRMED -- cross IS a weighted subgroup sum' if p1_ok else '*** FAILED -> REPLAN ***'}")
    if not p1_ok:
        log("   (the Phase 2/3 literature plan is aimed at the wrong object)")
        flush(); sys.exit(1)
    log("")

    # ================= H_CONC =================
    log("## H_CONC (*** THE FALSIFIER ***) -- does truncating v_2 <= T concentrate?")
    log("   err(T) = |cross_T - cross_full| / cross_full")
    log("   [!] STATED BEFORE THE RUN: p_v ~ 2^-v is geometric BY CONSTRUCTION, so a pass")
    log("      proves only the TAIL is geometric. It does NOT prove the HEAD is boundable,")
    log("      which is the actual math. NECESSARY, NOT SUFFICIENT.")
    log("")
    Ts = [5, 10, 15, 20, 25, 30, 40]
    log(f"{'q':>4} {'d':>3} " + "".join(f"{'T='+str(t):>12}" for t in Ts))
    conc_ok = True
    for q in [5, 7, 11, 13, 17, 31]:
        d = order_of_two(q)
        full = fulls[q]
        errs = []
        for T in Ts:
            ct = cross_subgroup(q, Tcap=T)
            errs.append(abs(ct - full) / abs(full))
        log(f"{q:>4} {d:>3} " + "".join(f"{e:>12.2e}" for e in errs))
        # check geometric decay while above noise
        for i in range(len(Ts) - 1):
            if errs[i] > 1e-14 and errs[i + 1] > 1e-14:
                if errs[i + 1] / errs[i] > 0.10:
                    conc_ok = False
    log("")
    log(f"   H_CONC: {'CONFIRMED -- the sum concentrates (tail geometric in T)' if conc_ok else 'NOT confirmed'}")
    log("   Reminder: this is the NECESSARY half. Phase 3a still needs the HEAD bounded.")
    flush()


def flush():
    with open("result_18_phase01_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
