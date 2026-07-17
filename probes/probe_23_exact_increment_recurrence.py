"""
PROBE 23 -- ROUTE 2: EXACT-RATIONAL cross(k) and increments c_k; hunt a linear recurrence.

PRE-REGISTRATION (written before running).
------------------------------------------------------------------
GOAL. R22 pinned rho_k = c_{k+1}/c_k as the Phase-3 object but its VALUE r_q (q>=5) was
float-wobbly (~0.3%) and pre-asymptotic. This probe computes cross(k) -- hence c_k --
EXACTLY (big-integer rationals, no float, no power iteration), to:
  (i)   CONFIRM q=3: is c_k = 7/15 EXACTLY (not just ~0.4655)?  [R15's slope, now exact]
  (ii)  rule the geometric (order-1) recurrence in/out EXACTLY at q>=5;
  (iii) fit a higher-order linear recurrence over Q if enough exact terms exist; its
        dominant root would be r_q as an algebraic number, and its ORDER previews the
        dimension of the route-1 transfer operator.

EXACT METHOD (no approximation anywhere):
  Cells (c_1..c_{k-1}, v_k), c_j in 1..m_j (m_j=d q^{j-1}), v_k in 1..M (M=ord_2(q^k)).
  mass(cell) = [prod_j 2^{-c_j}/(1-2^{-m_j})] * 2^{-v_k}/(1-2^{-M})
             = 2^{E} / D,  E = sum_j (m_j - c_j) + (M - v_k),
             D = prod_j (2^{m_j}-1) * (2^M-1).          [E integer, D fixed integer]
  value(cell) = sum_{m=1}^k q^{m-1} 2^{-S_m} mod q^k,  S_m = suffix sum. (only phi(q^k)
  distinct values). Accumulate N_r = sum_{cells->r} 2^E (exact big int), and
  SW = sum_cells 2^{2E}. Then, with P2 = (2^M+1)/(3(2^M-1)) exact (R8):
      cross(k) = (sum_r N_r^2 - SW) / D^2  /  P2^k        [ONE Fraction at the end]
  GATE: exact cross(k) must equal R18's float cross_from_cells to ~1e-12 rel.

DECISION RULES (pre-committed):
  G_EXACT (gate): |float(exact cross) - float cross_from_cells| / |..| < 1e-10 every (q,k).
      If FALSE the exact arithmetic is wrong -> STOP.
  Q3_SEVEN_FIFTEENTHS: c_k == 7/15 EXACTLY as Fractions for k>=4 at q=3.
      PRIOR: TRUE (R15 slope). If the exact Fractions are NOT 7/15 -> R15's slope was
      only asymptotic, a real correction -> report it.
  ORDER1 (geometric): c_{k+1}/c_k EQUAL as exact Fractions across k, per q.
      PRIOR: FALSE for q>=5 (R22 float ratios varied 0.51/0.62). Ruling it out EXACTLY
      (not "float differs") is the point.
  RECUR: report exact c_k and successive ratios; attempt an order-2 fit ONLY where >=5
      exact c-values exist. HONESTLY report "under-determined" otherwise -- NO fit to
      too-few points (feedback_r2_cannot_discriminate_monotone_fits).

BUDGET: cells = d^k q^{k(k-1)/2}. Cap 5e6 (integer ops, one-off ~1-2 min). No silent trunc.
  q=3: k<=5 (1.9M) -> c_2,c_3,c_4,c_5     q=5: k<=4 (4.0M) -> c_2,c_3,c_4
  q=7: k<=3 (9261) -> c_2,c_3             q=11: k<=3 (1.33M) -> c_2,c_3
  Under-determination at q>=5 is EXPECTED and is itself the finding that motivates route 1
  (build the operator; its dimension comes from STRUCTURE, not from fitting few points).

NOT AT STAKE: R10-R22, R5's rate, R6, R7, R12, THEOREM_C_745.
"""
import sys
from fractions import Fraction
from itertools import product

from probe_6_conservation_generalize import order_of_two
from probe_18_phase01_subgroup_form import cross_from_cells

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def exact_cross(q, k):
    """cross(k) as an exact Fraction via integer-numerator accumulation. No float."""
    N = q ** k
    M = order_of_two(N)
    d = order_of_two(q)
    inv2 = pow(2, -1, N)
    mods = [d * (q ** (j - 1)) for j in range(1, k)]      # m_1..m_{k-1}
    D = 1
    for mj in mods:
        D *= (2 ** mj - 1)
    D *= (2 ** M - 1)
    ncell = 1
    for mj in mods:
        ncell *= mj
    ncell *= M

    # precompute 2^{-c} mod N pieces are folded into the residue calc directly
    Nr = {}                 # residue -> integer numerator (sum of 2^E)
    SW = 0                  # sum_cells 2^{2E}   (the "within" piece, exact)
    reps = [range(1, mj + 1) for mj in mods]
    summods = sum(mods)
    for combo in product(*reps):
        base_E = summods - sum(combo)      # sum_j (m_j - c_j)
        for vk in range(1, M + 1):
            E = base_E + (M - vk)
            w = 1 << E                     # 2^E
            # value = sum_{m=1}^k q^{m-1} 2^{-S_m} mod N, S_m suffix sum of (combo..., vk)
            vs = combo + (vk,)
            val = 0
            s = 0
            for m in range(1, k + 1):
                s += vs[k - m]
                val = (val + (q ** (m - 1)) * pow(inv2, s, N)) % N
            Nr[val] = Nr.get(val, 0) + w
            SW += w * w
    sum_sq = sum(v * v for v in Nr.values())
    # cross = (sum_sq - SW)/D^2 / P2^k ;  P2 = (2^M+1)/(3(2^M-1))
    twoM = 2 ** M
    P2 = Fraction(twoM + 1, 3 * (twoM - 1))
    cross = Fraction(sum_sq - SW, D * D) / (P2 ** k)
    return cross, ncell


def main():
    log("# PROBE 23 -- ROUTE 2: EXACT cross(k), increments c_k, recurrence hunt")
    log("# Pre-reg: G_EXACT(gate) / Q3_SEVEN_FIFTEENTHS / ORDER1 / RECUR")
    log("")

    CASES = {3: [2, 3, 4, 5], 5: [2, 3, 4], 7: [2, 3], 11: [2, 3]}
    cross = {}

    log("## G_EXACT (gate) -- exact cross(k) vs R18 float cross_from_cells")
    log(f"   {'q':>4} {'k':>3} {'cells':>10} {'exact cross (float)':>22} {'R18 float':>18} {'|rel|':>10}")
    gate_ok = True
    for q, ks in CASES.items():
        for k in ks:
            ex, ncell = exact_cross(q, k)
            cross[(q, k)] = ex
            ref, _ = cross_from_cells(q, k)
            fe = float(ex)
            rel = abs(fe - ref) / abs(ref) if ref else 0.0
            if rel >= 1e-10:
                gate_ok = False
            log(f"   {q:>4} {k:>3} {ncell:>10} {fe:>22.15f} {ref:>18.12f} {rel:>10.2e}")
    log(f"   G_EXACT: {'CONFIRMED -- exact arithmetic matches' if gate_ok else '*** FAILED -> STOP ***'}")
    if not gate_ok:
        flush(); sys.exit(1)
    log("")

    # increments c_k
    log("## EXACT increments c_k = cross(k) - cross(k-1)  (cross(1)=0)")
    log("")
    ck = {}
    for q, ks in CASES.items():
        for k in ks:
            if k == 2:
                ck[(q, k)] = cross[(q, k)]           # cross(1)=0
            elif (q, k - 1) in cross:
                ck[(q, k)] = cross[(q, k)] - cross[(q, k - 1)]
        log(f"   q={q}:")
        for k in ks:
            if (q, k) in ck:
                c = ck[(q, k)]
                log(f"      c_{k} = {c}  = {float(c):.12f}")
        log("")

    # ---- Q3_SEVEN_FIFTEENTHS ----
    log("## Q3_SEVEN_FIFTEENTHS -- is c_k EXACTLY 7/15 at q=3 for k>=4?")
    target = Fraction(7, 15)
    for k in [3, 4, 5]:
        if (3, k) in ck:
            c = ck[(3, k)]
            eq = (c == target)
            log(f"   c_{k} = {c}   == 7/15 ? {eq}   (diff = {float(c - target):+.3e})")
    log("   (k=3 carries the x_2=2^-6 correction; k>=4 is where 7/15 should be exact.)")
    log("")

    # ---- ORDER1 (geometric) ----
    log("## ORDER1 -- is c_{k+1}/c_k EQUAL across k (geometric / order-1 recurrence)?")
    log("   EXACT ratios (not float 'differs'):")
    for q, ks in CASES.items():
        ratios = []
        for k in ks:
            if (q, k) in ck and (q, k + 1) in ck and ck[(q, k)] != 0:
                r = ck[(q, k + 1)] / ck[(q, k)]
                ratios.append((k, r))
        if not ratios:
            continue
        txt = "   ".join(f"c_{k+1}/c_{k} = {r} ({float(r):.6f})" for k, r in ratios)
        allsame = len({r for _, r in ratios}) == 1 if len(ratios) > 1 else None
        verdict = ("ALL EQUAL -> geometric" if allsame is True
                   else "NOT equal -> order >= 2" if allsame is False
                   else "only one ratio -- cannot judge")
        log(f"   q={q}: {txt}")
        log(f"          -> {verdict}")
    log("")

    # ---- RECUR ----
    log("## RECUR -- order-2 fit ONLY where >=5 exact c-values exist; else UNDER-DETERMINED")
    for q in CASES:
        cs = [ck[(q, k)] for k in CASES[q] if (q, k) in ck]
        n = len(cs)
        if n >= 5:
            # solve c_{k} = a c_{k-1} + b c_{k-2} from two eqns, verify on the rest
            import itertools as _it
            c0, c1, c2, c3, c4 = cs[:5]
            # [c2; c3] = [[c1,c0],[c2,c1]] [a;b]
            det = c1 * c1 - c0 * c2
            if det != 0:
                a = (c2 * c1 - c0 * c3) / det
                b = (c1 * c3 - c2 * c2) / det
                ok = all(cs[i] == a * cs[i - 1] + b * cs[i - 2] for i in range(2, n))
                log(f"   q={q}: order-2 fit a={a}, b={b}; verifies all terms? {ok}")
            else:
                log(f"   q={q}: singular -- cannot fit order-2")
        else:
            log(f"   q={q}: {n} exact c-values (<5) -> UNDER-DETERMINED for order-2. "
                f"Float (R22) shows order>=2, so this q needs route 1 (operator), NOT a fit.")
    log("")
    log("## READ: exact c_k confirm/deny 7/15 at q=3 and rule order-1 in/out at q>=5.")
    log("   Where under-determined, the ORDER is not fittable from data -> build the")
    log("   transfer operator (route 1); its dimension comes from the cell/subgroup")
    log("   STRUCTURE, and its 2nd eigenvalue IS r_q. That is the next probe.")
    flush()


def flush():
    with open("result_23_exact_increment_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
