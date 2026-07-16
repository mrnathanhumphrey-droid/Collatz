"""
PROBE 40 -- universality: is the phase boundary ALWAYS ord_q(w)=2, independent of the map?
(thread 2/3: a universality class of Syracuse-type maps.)

FRAMING. The qx+1 family (vary q, halve by 2) has ONE boundary: d=ord_q(2)=2 <=> q=3. "5x+1"
is NOT a different system -- it's q=5, already in-family. The genuine universality test is to
vary the HALVING BASE w: map (qx+1) w^{-v}, weights p_v=(w-1)w^{-v} (w-adic valuation model;
w=2 recovers Collatz). This is a well-defined q-adic self-similar measure for any w coprime to
q (T_v still contracts by 1/q q-adically; the "+1" factors out as a unit). The MECHANISM says
the gap closes when the halving phase w^{-v} mod q collapses to two values:
    gap absent  <=>  ord_q(w) = 2  <=>  w^2 ≡ 1, w ≢ 1  <=>  w ≡ -1 (mod q)  <=>  q | (w+1).
If TRUE independent of (q,w), the critical structure is a UNIVERSALITY CLASS keyed to d=2, and
it predicts a NEW critical point at (q=5, w=4) [4 ≡ -1 mod 5] AWAY from q=3.

GAP DIAGNOSTIC (direct, no operator). value(v_1..v_k) = sum_m q^{m-1} w^{-S_m} mod q^k,
S_m = suffix sums. Unnormalized cell mass m_cell = sum_{addresses->cell} prod_i w^{-v_i}.
With U2 = sum_{v<=V} w^{-2v} (= lambda_1 numerator, since lambda_1 = U2/U1^2 and the U1^{2k}
cancels), the normalized correlation is
    X_k = (sum_cell m_cell^2) / U2^k ,   Delta_k = X_k - X_{k-1}.
GAP  <=> Delta_k DECAYS geometrically (ratio Delta_{k+1}/Delta_k -> r < 1).
NO GAP <=> Delta_k -> const > 0 (linear growth in X_k), ratio -> 1.

PRE-REGISTRATION (falsifier-first; committed before running).
------------------------------------------------------------------
H_UNIV (*** the claim ***): gap absent IFF w ≡ -1 mod q (ord_q(w)=2), for ALL tested (q,w).
    Predicted NO-GAP (ratio ~1): (3,2)[known crit], (5,4), (7,6), (11,10).
    Predicted GAP  (ratio <1):   (5,2)[r5~.62], (5,3), (7,2)[r7~.39], (7,3), (11,2).
    PRIOR: TRUE -> universality class keyed to d=ord_q(w)=2; NEW critical point at (5,4)/(7,6).
    FALSIFIER: if (5,4) shows a GAP (ratio <1), the boundary is NOT ord_q(w)=2 -> w is
    RELEVANT, no universality class. Honest negative = a win (per the arc's ethos).
H_ARC (arc scar tissue): cross-domain "universality" has a body count this corpus (Collatz
    walls-not-threads; Solar/Cosmology killed). So the verdict is EXACT (ratio->1 vs <1 on the
    predicted split), not a vibe; if the split is not clean, report NO clean class.

DECISION RULE: for each (q,w), classify GAP if Delta_4/Delta_3 < 0.85, NO-GAP if > 0.95,
    AMBIG otherwise. H_UNIV CONFIRMED iff the GAP/NO-GAP split matches (w ≡ -1 mod q) exactly.

BUDGET: direct enumeration V^k. Large-w tails decay faster so V shrinks with w. Cap V^k <= 4M.
    Foreground, seconds each. No operator, no fit.

NOT AT STAKE: R10-R39 (qx+1 = w=2 column). A refutation kills the universality hypothesis only.
"""
from math import log10
from collections import defaultdict
from itertools import product

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def ord_mod(w, q):
    x, o = w % q, 1
    while x != 1:
        x = (x * w) % q
        o += 1
    return o


def X_of_k(q, w, k, V):
    """Normalized correlation X_k = sum_cell m_cell^2 / U2^k, m_cell unnormalized (prod w^{-v})."""
    N = q ** k
    invw = pow(w, -1, N)
    invw_pow = [1] * (k * V + 2)
    for i in range(1, len(invw_pow)):
        invw_pow[i] = (invw_pow[i - 1] * invw) % N
    qp = [q ** m for m in range(k)]
    # unnormalized address weight prod w^{-v_i} as float; group by value mod q^k
    cell = defaultdict(float)
    for addr in product(range(1, V + 1), repeat=k):
        s = 0
        val = 0
        wt = 1.0
        for m in range(1, k + 1):
            s += addr[k - m]            # suffix sum S_m
            val = (val + qp[m - 1] * invw_pow[s]) % N
            wt *= w ** (-addr[k - m])   # = prod w^{-v_i}
        cell[val] += wt
    U2 = sum(w ** (-2.0 * v) for v in range(1, V + 1))
    return sum(m * m for m in cell.values()) / (U2 ** k)


def analyze(q, w, kmax=4):
    V = max(10, int(13.0 / log10(w)) + 3)
    while V ** kmax > 4_000_000:
        V -= 1
    Xs = {}
    for k in range(1, kmax + 1):
        Xs[k] = X_of_k(q, w, k, V)
    d = {k: Xs[k] - Xs[k - 1] for k in range(2, kmax + 1)}
    ratio = (d[kmax] / d[kmax - 1]) if abs(d[kmax - 1]) > 1e-14 else float('nan')
    return V, Xs, d, ratio


def main():
    log("# PROBE 40 -- universality: gap absent IFF w = -1 mod q (ord_q(w)=2), independent of q?")
    log("# (w=2 = Collatz. NEW critical point predicted at (q=5,w=4).)")
    log("")
    CASES = [(3, 2), (5, 2), (5, 3), (5, 4), (7, 2), (7, 3), (7, 6), (11, 2), (11, 10)]
    log(f"   {'q':>4} {'w':>3} {'ord_q(w)':>9} {'w=-1?':>7} {'V':>4} "
        f"{'Delta_2':>11} {'Delta_3':>11} {'Delta_4':>11} {'D4/D3':>9} {'class':>8} {'predict':>8}")
    ok = True
    for q, w in CASES:
        o = ord_mod(w, q)
        wm1 = (w % q) == (q - 1)               # w == -1 mod q
        V, Xs, d, ratio = analyze(q, w, 4)
        if ratio != ratio:                      # nan
            cls = "FLAT0"
        elif ratio < 0.85:
            cls = "GAP"
        elif ratio > 0.95:
            cls = "NO-GAP"
        else:
            cls = "AMBIG"
        predict = "NO-GAP" if wm1 else "GAP"
        agree = (cls == predict) or (cls == "FLAT0" and predict == "NO-GAP")
        if not agree:
            ok = False
        log(f"   {q:>4} {w:>3} {o:>9} {str(wm1):>7} {V:>4} "
            f"{d[2]:>11.3e} {d[3]:>11.3e} {d[4]:>11.3e} {ratio:>9.4f} {cls:>8} {predict:>8}"
            f"{'' if agree else '  <-- MISMATCH'}")
    log("")
    if ok:
        log("## H_UNIV CONFIRMED -- the GAP/NO-GAP split matches (w = -1 mod q) at every (q,w).")
        log("   => the phase boundary is ord_q(w)=2 UNIVERSALLY, not q=3 specifically. A NEW")
        log("      critical point exists at (q=5,w=4), (q=7,w=6), ... -- a universality class of")
        log("      Syracuse-type maps keyed to d=2 (halving phase collapses to {1,-1}). w is")
        log("      IRRELEVANT to the boundary; only ord_q(w)=2 matters. (Ties thread 3 to thread 1:")
        log("      the SAME order-2 EP should sit at every ord_q(w)=2 boundary.)")
    else:
        log("## H_UNIV REFUTED (or ambiguous) -- the boundary does NOT track w = -1 mod q cleanly.")
        log("   Inspect mismatches: w may be RELEVANT, or the diagnostic needs higher k. Honest")
        log("   negative -- no clean universality class (the arc's cross-domain analogies have a")
        log("   body count; this would join them).")
    with open("result_40_universality_basew_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
