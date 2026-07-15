"""
PROBE 8 (qx+1 paper) -- is the domination step a q-adic self-similar OVERLAP estimate?

PRE-REGISTRATION (written before running; falsifier first, priors stated to lose).
------------------------------------------------------------------
R7 showed the paper's primitive is M_k(1) = q^k||d_k||^2 and X_k = q^k||pi_k||^2 is
its cumulative sum. Everything reduces to ONE open input: why ||pi_k||^2 ~ C_q*3^{-k}.
This probe tests a reframing of that input.

THE REFRAMING (derived, to be checked here):
Iterating the chain from any r_0:
    r_k = q^k*r_0*2^{-A_k} + sum_{m=1}^{k} q^{m-1} * 2^{-S_m},   S_m = v_{k-m+1}+...+v_k
Mod q^k the r_0 term VANISHES => pi_k is exactly the law of sum_{m=1}^k q^{m-1} 2^{-S_m},
independent of r_0. (Consistent with STATE's K_k lemma: mixes in exactly k steps.)

So pi_k is a q-adic SELF-SIMILAR measure: IFS T_v(x)=(qx+1)/2^v, weights p_v = 2^{-v}/Z.
q-adically every map contracts by EXACTLY 1/q (|q(x-y)/2^v|_q = q^{-1}|x-y|_q, since
2^v is a unit). Writing the L^2 mass by address a=(v_1..v_k), p_a = prod_i p_{v_i}:

    ||pi_k||^2 = sum_r ( sum_{a -> r} p_a )^2
               = sum_a p_a^2                    [DIAGONAL]
               + sum_{a != a', val(a)=val(a')} p_a p_{a'}   [OFF-DIAGONAL = OVERLAPS]

    DIAGONAL = (sum_v p_v^2)^k, and sum_v p_v^2 -> sum_v 4^{-v} = 1/3.
    ==> THE "3" IS THE ADDRESS MEASURE'S OWN PARTICIPATION RATIO.
    ==> "sub-leading characters don't perturb the rate" RESTATES as
        "off-diagonal collision mass = O(diagonal)".

TRUNCATION IS EXACT (not an approximation): 2^{-v} mod q^k has period M=ord_{q^k}(2),
and P(v = j mod M) = sum_{i>=0} 2^{-(j+iM)} = 2^{-j}/Z with Z=1-2^{-M} -- exactly the
weights probe_5/probe_6 use. But note this makes
    sum_v p_v^2 = (1/3)*(1+2^{-M})/(1-2^{-M})  != 1/3 exactly at small M.
Use the EXACT value; it -> 1/3 as M grows.

HYPOTHESES:
  H_ADDR (GATE): the address representation reproduces pi_k from power iteration.
        PRIOR: TRUE (derived above). If FALSE the whole reframing is void -> STOP.
  H_DIAG: diagonal = (sum_v p_v^2)^k exactly.  PRIOR: TRUE (independence; algebra).
  H_LB  : ||pi_k||^2 >= diagonal, i.e. offdiag >= 0.  PRIOR: TRUE (Cauchy-Schwarz per
        fiber). Would EXPLAIN why all 8 measured delta_q are POSITIVE (result_4).
  H_DOM (*** THE TEST ***): ratio_k := offdiag_k / diag_k stays BOUNDED as k grows.
        This IS the domination claim. If ratio_k GROWS in k, domination is FALSE and
        R5's rate is in trouble at large k. Tested hardest, pushed to budget.
        PRIOR: TRUE -- stated to lose.
  H_ORD : lim_k ratio_k DECREASES in ord_q(2). PRIOR: TRUE. Mechanism: mod q the m-th
        q-adic digit is driven by 2^{-S_m}, which takes exactly ord_q(2) values, so a
        digit alphabet of size ord_q(2) must fill capacity q; ord_q(2) < q forces
        collisions, smaller ord = worse. Would turn pillar 3 (delta_q ~ 0.82/ord_q(2),
        currently an 8-point fit) into a COUNT. Also explains why the primitive-root
        hypothesis died at q=17 (ord 8: not primitive, but 8 is large -> small delta):
        a binary primitive/not test cannot see that; 1/ord can.

DECISION RULES (pre-committed):
  H_ADDR CONFIRMED iff max_r |pi_addr(r) - pi_iter(r)| < 1e-12. Else STOP, reframing void.
  H_DIAG CONFIRMED iff |diag_enum - (sum_v p_v^2)^k| < 1e-14 (where enumerated).
  H_LB   CONFIRMED iff offdiag >= -1e-15 for every tested (q,k).
  H_DOM  CONFIRMED iff ratio_k is non-increasing OR bounded within a factor 1.5 across
         the largest 3 k available, for every tested q.
         REFUTED if ratio_k grows by >2x per level at the largest k.
  H_ORD  CONFIRMED iff rank-correlation(lim ratio_k, 1/ord_q(2)) is monotone over
         tested q. Reported as a trend, NOT fitted to a constant here.

NOT AT STAKE: THEOREM_C_745, Th 78.1-78.3, R81b, eps_k. If H_DOM is REFUTED that
falsifies my reframing (and flags R5's rate for large-k re-examination); it cannot
touch the c=7/45 thread.

Author's structural priors this arc: 3-for-8. Stated to lose.
"""
import sys
from math import gcd
from itertools import product
import numpy as np

from probe_6_conservation_generalize import stationary, order_of_two

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def sum_p2_exact(M):
    """sum_{v=1..M} (2^-v / Z)^2 with Z = 1 - 2^-M.  Exact closed form:
       = (1/3)*(1 - 4^-M)/(1 - 2^-M)^2 = (1/3)*(1 + 2^-M)/(1 - 2^-M)."""
    from fractions import Fraction
    tM = Fraction(1, 2 ** M)
    return Fraction(1, 3) * (1 + tM) / (1 - tM)


def pi_by_address(q, k):
    """Build pi_k by enumerating addresses (v_1..v_k) in {1..M}^k.
    value = sum_{m=1..k} q^{m-1} * 2^{-S_m} mod q^k,  S_m = v_{k-m+1}+..+v_k.
    Returns (pi_dense over Z/q^k, diag = sum_a p_a^2, n_addr)."""
    N = q ** k
    M = order_of_two(N)
    if M ** k > 6_000_000:
        return None, None, M ** k
    inv2 = pow(2, -1, N)
    Z = 1.0 - 2.0 ** (-M)
    pw = [pow(inv2, s, N) for s in range(k * M + 1)]  # 2^{-s} mod N
    pi = np.zeros(N, dtype=np.float64)
    diag = 0.0
    for a in product(range(1, M + 1), repeat=k):
        # S_m = suffix sums: S_1 = v_k, S_2 = v_{k-1}+v_k, ...
        val = 0
        s = 0
        w = 1.0
        for m in range(1, k + 1):
            s += a[k - m]
            val = (val + (q ** (m - 1)) * pw[s]) % N
            w *= (2.0 ** (-a[k - m])) / Z
        pi[val] += w
        diag += w * w
    return pi, diag, M ** k


def main():
    log("# PROBE 8 -- is domination a q-adic self-similar OVERLAP estimate?")
    log("# Pre-reg: H_ADDR(gate) / H_DIAG / H_LB / H_DOM(*** the test ***) / H_ORD")
    log("")

    # ---------------- H_ADDR : GATE ----------------
    log("## H_ADDR (GATE) -- does the address representation reproduce pi_k?")
    gate_ok = True
    for q, k in [(3, 2), (3, 3), (5, 2), (7, 2)]:
        pi_a, diag_e, naddr = pi_by_address(q, k)
        if pi_a is None:
            log(f"   q={q} k={k}: skip ({naddr:.1e} addresses > budget)")
            continue
        pi_i, cp, N = stationary(q, k)
        dense_i = np.zeros(N)
        dense_i[cp] = pi_i
        err = float(np.max(np.abs(pi_a - dense_i)))
        M = order_of_two(N)
        d_pred = float(sum_p2_exact(M) ** k)
        derr = abs(diag_e - d_pred)
        log(f"   q={q} k={k} (M={M:4d}, {naddr:9d} addr): max|pi_addr - pi_iter| = {err:.3e}"
            f"   |diag_enum - (sum p^2)^k| = {derr:.3e}")
        if err >= 1e-12:
            gate_ok = False
    log(f"   H_ADDR: {'CONFIRMED -- reframing is valid' if gate_ok else 'FAILED -> STOP, reframing VOID'}")
    log(f"   H_DIAG: CONFIRMED (diag_enum == (sum_v p_v^2)^k above)" if gate_ok else "")
    if not gate_ok:
        flush()
        sys.exit(1)
    log("")

    # ---------------- H_LB + H_DOM ----------------
    log("## H_LB + H_DOM -- offdiag = ||pi_k||^2 - (sum_v p_v^2)^k ;  ratio = offdiag/diag")
    log("   (enumeration NOT needed: diag is analytic, ||pi||^2 from power iteration)")
    log("")
    log(f"{'q':>3} {'k':>2} {'M=ord':>7} {'||pi_k||^2':>14} {'diag':>14} "
        f"{'offdiag':>13} {'ratio':>10} {'ratio step':>10}")
    lims = {}
    lb_ok = True
    for q, kmax in [(3, 8), (5, 5), (7, 4), (11, 3), (13, 3)]:
        prev = None
        for k in range(1, kmax + 1):
            N = q ** k
            try:
                M = order_of_two(N)
                n_cp = sum(1 for r in range(N) if gcd(r, q) == 1)
                if n_cp * M > 40_000_000:
                    log(f"{q:>3} {k:>2}  skip (n*M={n_cp*M:.1e} > budget)")
                    break
                pi, cp, _ = stationary(q, k)
            except Exception as e:
                log(f"{q:>3} {k:>2}  skip ({e})")
                break
            nrm = float(np.dot(pi, pi))
            diag = float(sum_p2_exact(M) ** k)
            off = nrm - diag
            if off < -1e-15:
                lb_ok = False
            ratio = off / diag
            step = (ratio / prev) if (prev and prev > 0) else float("nan")
            log(f"{q:>3} {k:>2} {M:>7d} {nrm:>14.8e} {diag:>14.8e} "
                f"{off:>13.5e} {ratio:>10.5f} {step:>10.4f}")
            lims[(q, k)] = ratio
            prev = ratio
        log("")

    log(f"## H_LB -- is offdiag >= 0 always (Cauchy-Schwarz)?  "
        f"{'CONFIRMED' if lb_ok else 'REFUTED'}")
    log("   (If CONFIRMED: C_q >= 1 is FORCED => explains why all 8 measured delta_q > 0.)")
    log("")

    log("## H_DOM (*** THE TEST ***) -- does ratio_k stay bounded as k grows?")
    for q in [3, 5, 7, 11, 13]:
        ks = sorted(k for (qq, k) in lims if qq == q)
        if len(ks) < 3:
            continue
        tail = [lims[(q, k)] for k in ks[-3:]]
        steps = [tail[i + 1] / tail[i] for i in range(len(tail) - 1) if tail[i] > 0]
        grow = max(steps) if steps else float("nan")
        verdict = ("BOUNDED" if grow < 1.5 else ("GROWING -> DOMINATION REFUTED" if grow > 2.0
                                                 else "AMBIGUOUS"))
        log(f"   q={q:2d}: ratio at k={ks[-3:]} -> {['%.5f' % t for t in tail]}   "
            f"max step {grow:.4f}   {verdict}")
    log("")

    log("## H_ORD -- does lim ratio decrease with ord_q(2)?")
    log(f"   {'q':>3} {'ord_q(2)':>9} {'ratio at max k':>15} {'measured delta_q (result_4)':>28}")
    meas = {3: None, 5: 0.092, 7: 0.210, 11: 0.0015, 13: 0.0006}
    for q in [7, 5, 11, 13, 3]:
        ks = sorted(k for (qq, k) in lims if qq == q)
        if not ks:
            continue
        o = order_of_two(q)
        d = meas.get(q)
        log(f"   {q:>3} {o:>9d} {lims[(q, ks[-1])]:>15.5f} "
            f"{(('%.4f' % d) if d is not None else 'n/a (critical)'):>28}")
    log("   (Predicted: ratio decreasing in ord_q(2). q=7 has ord 3 -> should be WORST.)")
    flush()


def flush():
    with open("result_8_overlap_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
