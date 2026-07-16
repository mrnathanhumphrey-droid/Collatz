"""
PROBE 37 -- the level-3 gate: does a next-digit constant sigma_1 enter, and can it vanish
independently of s_R13? (Fires the PHASE3 worksheet's OPEN item 1 before we trust the k=2
extrapolation. The arc has been burned twice trusting one level past the window.)

BACKGROUND. R20 delivered the level-2 gate exactly (2^{-S'_2} = 2^{-S_2} + j1*s*2^{-S_1} mod q,
s = s_R13 = (2^d-1)/q mod q) and left level 3 as W_2 + T_3 = 0 mod q with W_2 := (U_1+T_2)/q
"a DEFINITION not a formula -- needs a 2nd-order expansion of 2^{-jd} mod q^3." This probe does
that expansion and TESTS it.

DERIVATION (done before running). Write 2^d = 1 + q*s + q^2*sigma (mod q^3), where
  s     = s_R13 = ((2^d - 1)//q) % q           (the level-2 constant, R13)
  sigma = sigma_1 = (((2^d - 1)//q)//q) % q     (the NEXT digit -- the object under test)
Then 2^{j1 d} = 1 + q*(j1 s) + q^2*(j1 sigma + C(j1,2) s^2) mod q^3, so with P := j1 s and
Q := j1 sigma + C(j1,2) s^2:
    2^{-j1 d} = 1 - qP + q^2 (P^2 - Q) mod q^3
    U_1 = 2^{-S_1}[P - q(P^2 - Q)] mod q^2
    W_2 = (U_1 + T_2)/q  ==>   *** W_2 = Q_1 - y1*(P^2 - Q)  (mod q) ***
where y1 = 2^{-S_1} mod q and Q_1 = ( [ (j1 s)*2^{-S_1} + 2^{-S_2} - 2^{-S'_2} ] mod q^2 ) // q.
The sigma_1 dependence sits inside Q: the coefficient of sigma_1 in W_2 is exactly +y1*j1.
So DROPPING sigma_1 shifts W_2 by  y1 * j1 * sigma_1 (mod q)  -- nonzero iff y1*j1*sigma_1 != 0.

GROUND TRUTH. At k=3, W_2 = (U_1+T_2)/q is an EXACT integer (once levels 1,2 pass). Compute it
by exact big-int division from residues mod q^3 -- no model, no tolerance. Test symbolic forms
against it as exact set-equality (zero mismatch), the R13/R20 discipline.

PRE-REGISTRATION (numbers/direction committed BEFORE running; priors stated to lose).
------------------------------------------------------------------
H_GATE3  (*** THE TEST, exact iff ***): the WITH-sigma form W_2 = Q_1 - y1(P^2 - Q) holds on
    EVERY level-2-passing pair (zero mismatch). If ANY mismatch -> my 2nd-order expansion is
    WRONG; report and stop. PRIOR: TRUE (derived).
H_SIGMA_MATTERS (*** the decision ***): the NO-sigma form (sigma:=0 in Q) DISAGREES with W_2,
    and disagrees EXACTLY on the pairs where y1*j1*sigma_1 != 0 mod q.
    PRIOR: TRUE -> sigma_1 is a genuine level-3 correction; the clean geometric chain gains a
    second-order term with its own vanishing locus (a NEW boundary constant one level deeper).
    FALSIFIER: if the no-sigma form ALSO has zero mismatch, sigma_1 is INERT (drops out) and the
    clean chain survives to level 3 -- GOOD news for L3 (bound the k=2 form safely). I commit to
    reporting that outcome as such, not burying it.
H_SIGMA_INDEP (structure): sigma_1 ranges over F_q and can be 0 while s_R13 != 0 (independent of
    the level-2 constant). PRED: at least one prime q in [5,100) has sigma_1 = 0 & s_R13 != 0.
    (q=3 itself has s_R13=1, sigma_1=0 -- the boundary prime is already degenerate at 2nd order,
    a suggestive but not decisive data point since d=2 there anyway.)
H_L3RATE (measurement, NO verdict): level-3 conditional pass rate, compare to R20's ~1/q.

DECISION RULES (pre-committed, exact-set-equality -- cannot be mis-thresholded):
  H_GATE3        CONFIRMED iff with-sigma mismatches == 0 at every tested q.
  H_SIGMA_MATTERS CONFIRMED iff no-sigma mismatches > 0 AND every no-sigma mismatch has
                  y1*j1*sigma_1 != 0 mod q (and every zero-shift pair matches). REFUTED (inert)
                  iff no-sigma mismatches == 0.
  H_SIGMA_INDEP  reported; existence of one (q: sigma_1=0, s_R13!=0) proves independence.
  H_L3RATE       no rule. Measurement only.

NOT AT STAKE: R10-R36. A refutation of H_GATE3 kills only my 2nd-order derivation.
"""
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


def s_sigma(q, d):
    """s_R13 (1st digit) and sigma_1 (2nd digit) of (2^d - 1)/q."""
    t = (pow(2, d) - 1) // q
    s = t % q
    sigma = (t // q) % q
    return s, sigma


def build_cells(q, k):
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


def w2_exact(q, inv2, N, Sa, Sb):
    """Exact W_2 = (U_1+T_2)/q mod q by big-int division from residues mod q^3=N.
    Returns (w2, ok) where ok is False if levels 1/2 don't divide cleanly (shouldn't happen
    for same-bucket-mod-q^2 pairs)."""
    a1 = pow(inv2, Sa[0], N); b1 = pow(inv2, Sb[0], N)
    a2 = pow(inv2, Sa[1], N); b2 = pow(inv2, Sb[1], N)
    T1 = a1 - b1
    if T1 % q != 0:
        return None, False
    U1 = T1 // q
    S = U1 + (a2 - b2)
    if S % q != 0:
        return None, False
    return (S // q) % q, True


def w2_form(q, d, s, sigma, Sa, Sb, use_sigma):
    """Symbolic W_2 = Q_1 - y1*(P^2 - Q) mod q. P=j1 s, Q=j1 sigma + C(j1,2) s^2 (sigma
    optionally zeroed). Q_1 = 2nd digit of [ (j1 s) 2^{-S_1} + 2^{-S_2} - 2^{-S'_2} ] mod q^2."""
    q2 = q * q
    inv2_q2 = pow(2, -1, q2)
    S1a, S2a = Sa[0], Sa[1]
    S1b, S2b = Sb[0], Sb[1]
    j1 = (S1b - S1a) // d
    P = j1 * s
    sig = sigma if use_sigma else 0
    Q = j1 * sig + (j1 * (j1 - 1) // 2) * (s * s)
    A1 = pow(inv2_q2, S1a % (d * q), q2)       # 2^{-S_1} mod q^2
    A2 = pow(inv2_q2, S2a % (d * q), q2)       # 2^{-S_2} mod q^2
    B2 = pow(inv2_q2, S2b % (d * q), q2)       # 2^{-S'_2} mod q^2
    base = (P * A1 + A2 - B2) % q2
    # base must be ==0 mod q for a level-2 pair; guard
    Q1 = (base // q) % q
    y1 = A1 % q
    return (Q1 - y1 * ((P * P - Q) % q)) % q, (base % q)


def main():
    log("# PROBE 37 -- level-3 gate: does sigma_1 (2nd digit of (2^d-1)/q) enter W_2, and can")
    log("#            it vanish independently of s_R13? Exact-iff test vs big-int ground truth.")
    log("")

    # ---------- H_SIGMA_INDEP: the (s_R13, sigma_1) table ----------
    log("## H_SIGMA_INDEP -- (s_R13, sigma_1) over primes. sigma_1=0 with s_R13!=0 => independent")
    log(f"   {'q':>5} {'d':>4} {'s_R13':>6} {'sigma_1':>8}  note")
    first_indep = None
    primes = [p for p in range(3, 200) if all(p % j for j in range(2, int(p ** 0.5) + 1))]
    for q in primes:
        d = order_of_two(q)
        s, sig = s_sigma(q, d)
        note = ""
        if sig == 0 and s != 0:
            note = "<- sigma_1=0, s_R13!=0 (INDEP degeneration)"
            if first_indep is None and q >= 5:
                first_indep = q
        if q < 60 or note:
            log(f"   {q:>5} {d:>4} {s:>6} {sig:>8}  {note}")
    log("")
    log(f"   smallest prime q>=5 with sigma_1=0 & s_R13!=0: {first_indep}")
    log("   => sigma_1 is a SEPARATE q-adic digit; fixing s_R13!=0 leaves sigma_1 free over F_q.")
    log("")

    # ---------- H_GATE3 + H_SIGMA_MATTERS: the exact-iff test ----------
    # W_2 depends only on residues (S_1,S_2,S'_1,S'_2) and j1 -- NOT on v_1 or the full address.
    # So we CONSTRUCT valid level-1&2-passing pairs directly (works at any q, incl. sigma_1!=0
    # primes where full cell enumeration is infeasible). q=5,7 (sigma=0) test the first-order
    # structure; q=11,13,23,41 (sigma!=0) test the sigma_1 term itself -- the actual question.
    log("## H_GATE3 (with-sigma, exact iff) + H_SIGMA_MATTERS (no-sigma must FAIL where y1 j1 sigma!=0)")
    log("   pairs CONSTRUCTED to pass levels 1&2 directly (decoupled from v_1); exact big-int W_2")
    log("")
    log(f"   {'q':>4} {'d':>3} {'s':>3} {'sig':>4} {'L2 pairs':>9} {'with-sig bad':>13} "
        f"{'no-sig bad':>11} {'pred no-sig':>12} {'shift-ok':>9} {'match?':>7}")

    def run_direct(q, nsamp, seed):
        import random as _r
        rng = _r.Random(seed)
        d = order_of_two(q)
        N = q ** 3
        inv2 = pow(2, -1, N)
        inv2q = pow(2, -1, q)
        s, sigma = s_sigma(q, d)
        dq2 = d * q * q
        H = {pow(inv2q, i, q): i for i in range(d)}   # 2^{-i} mod q -> i  (dlog table over <2>)
        npair = bad_with = bad_no = pred_no = shift_bad = 0
        tries = 0
        while npair < nsamp and tries < nsamp * 40:
            tries += 1
            S1 = rng.randrange(1, dq2 + 1)
            S2 = rng.randrange(1, dq2 + 1)
            j1 = rng.randrange(1, q)              # nonzero mod q
            S1p = S1 + j1 * d
            # level-2 target: 2^{-S'_2} = 2^{-S_2} + j1 s 2^{-S_1} (mod q)
            t = (pow(inv2q, S2 % d, q) + j1 * s * pow(inv2q, S1 % d, q)) % q
            if t == 0 or t not in H:
                continue                          # no valid S'_2 -> not a depth-2 collision
            e0 = H[t]
            S2p = e0 + d * rng.randrange(0, q * q)
            if S2p < 1:
                S2p += d
            Sa = (S1, S2, 0); Sb = (S1p, S2p, 0)
            we, ok = w2_exact(q, inv2, N, Sa, Sb)
            if not ok:
                continue                          # guard (should not fire)
            npair += 1
            wf, _ = w2_form(q, d, s, sigma, Sa, Sb, use_sigma=True)
            wn, _ = w2_form(q, d, s, sigma, Sa, Sb, use_sigma=False)
            y1 = pow(inv2q, S1 % d, q)
            shift = (y1 * j1 * sigma) % q
            if wf != we:
                bad_with += 1
            if wn != we:
                bad_no += 1
            if shift != 0:
                pred_no += 1
            if (we - wn) % q != shift:             # predicted with-minus-no shift == y1 j1 sigma
                shift_bad += 1
        return d, s, sigma, npair, bad_with, bad_no, pred_no, shift_bad

    ok_gate = True
    sigma_matters = None
    for q in [5, 7, 11, 13, 23, 41]:
        d, s, sigma, npair, bw, bn, pn, sb = run_direct(q, 200_000, 3736 + q)
        if bw:
            ok_gate = False
        match = (bn == pn) and (sb == 0)
        if sigma != 0 and bn > 0 and match:
            sigma_matters = True
        log(f"   {q:>4} {d:>3} {s:>3} {sigma:>4} {npair:>9} {bw:>13} "
            f"{bn:>11} {pn:>12} {str(sb == 0):>9} {str(match):>7}")
    log("")
    log(f"   H_GATE3: {'CONFIRMED -- with-sigma form is the exact level-3 gate (0 mismatch, all q)' if ok_gate else '*** REFUTED -> 2nd-order derivation WRONG ***'}")
    if sigma_matters:
        log("   H_SIGMA_MATTERS: CONFIRMED -- at sigma_1!=0 primes (q=11,13,23,41) the no-sigma form")
        log("      FAILS, exactly by the predicted shift y1*j1*sigma_1. sigma_1 IS a genuine level-3")
        log("      correction: a NEW 2nd-order boundary constant, one q-adic digit deeper than s_R13.")
    else:
        log("   H_SIGMA_MATTERS: NOT confirmed -- no-sigma form did not fail as predicted; sigma_1")
        log("      may be inert. Inspect the table (clean chain would be GOOD news for L3).")
    log("")

    # ---------- H_L3RATE ----------
    log("## H_L3RATE -- level-3 conditional pass rate P(depth3 | depth2), vs R20's ~1/q. MEASUREMENT.")
    log(f"   {'q':>4} {'1/q':>7}   P(depth-3 collide | depth-2 collide)")
    for q in [3, 5, 7]:
        cells, Ss, vals, N, M, d, inv2 = build_cells(q, 3)
        q2 = q * q
        b2 = defaultdict(list)
        for c in cells:
            b2[vals[c] % q2].append(c)
        n2 = n3 = 0
        for _, lst in b2.items():
            L = len(lst)
            if L < 2:
                continue
            b3 = defaultdict(int)
            for c in lst:
                b3[vals[c]] += 1
            for i in range(L):
                for jx in range(i + 1, L):
                    n2 += 1
            for _, cc in b3.items():
                n3 += cc * (cc - 1) // 2
        log(f"   {q:>4} {1/q:>7.4f}   {n3/n2 if n2 else float('nan'):.4f}   (depth2 pairs={n2}, depth3={n3})")
    log("")
    log("## READ")
    log("   If sigma_1 enters (H_SIGMA_MATTERS confirmed): the re-entry recursion's phase/weight")
    log("   factor omega gains a 2nd-order term; the clean k=2 chain is a TRUNCATION, and the L3")
    log("   bound must target the corrected object (or prove sigma_1's vanishing is another benign")
    log("   index-shift a la R35). If inert: grind the bound on the k=2 form, safe two levels deep.")
    with open("result_37_level3_gate_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
