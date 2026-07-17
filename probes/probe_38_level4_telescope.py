"""
PROBE 38 -- level-4 gate: does a THIRD digit sigma_2 enter, and do the digit corrections
TELESCOPE into a generating function for omega? (Fires R37's follow-up: worksheet open item 3.)

CONTEXT. R37 derived the level-3 gate exactly and found the 2nd q-adic digit sigma_1 of
(2^d-1)/q enters W_2 with coefficient +y1*j1, but is NOT a new boundary (q=5,7 gap at sigma_1=0).
Two questions remain:
  (Q1) Does the 3rd digit sigma_2 enter the level-4 gate W_3? (structural)
  (Q2) Does it enter with the SAME coefficient y1*j1 as sigma_1 did? If every digit sigma_i of
       (2^d-1)/q enters W_{i+1} with the common coefficient y1*j1, the tower telescopes:
          digit-part of the cascade = y1*j1 * sum_i q^i sigma_i = y1*j1 * (2^d-1)/q
       i.e. the corrections REASSEMBLE into the exact tower constant -- a generating function
       for omega, and the closed form the L3 bound wants.

KEY METHOD (no hand-derivation of W_3). The digits (s, sigma_1, sigma_2, ...) enter the cascade
ONLY through 2^{-j1 d} (j1*d is the one guaranteed multiple of d = the level-1 shift). So we
SUBSTITUTE a truncated model for that single quantity and read the effect on the exact W_3:
  pow2_true   = 2^{-j1 d} mod q^4                 (exact; equals inverse of (1+qs+q^2 s1+q^3 s2)^j1)
  pow2_model  = inverse of (1 + q*s + q^2*sigma_1)^j1 mod q^4     (order-2: sigma_2 DROPPED)
Run the exact cascade (U1 -> W2 -> W3) with 2^{-S'_1} = 2^{-S_1}*pow2, once true once model.
  - pow2_model = pow2_true mod q^3 (differ only at the q^3 digit, exactly by j1*sigma_2), so the
    model still passes levels 1,2,3 cleanly (W_3 well-defined) and W_3_model differs from
    W_3_true only by the sigma_2 image. NO closed form for W_3 needed; ground truth is big-int.

GROUND TRUTH: W_3 = (W_2+T_3)/q, W_2=(U_1+T_2)/q, U_1=T_1/q -- exact integer divisions at k=4.
Pairs CONSTRUCTED to pass levels 1,2,3 (W_3 depends only on S_1,S_2,S_3,S'_1,S'_2,S'_3,j1).

BENIGN NOTE (settled, no new run): sigma_2=0 (indeed ALL higher digits =0) holds at q=5,7 since
(2^d-1)/q is a single digit there; both GAP (r5~0.62, r7~0.38, full-cascade). So "higher digit =0
still gaps" is already witnessed -- sigma_2 (like sigma_1) is not a new boundary. R38 tests only
whether sigma_2 ENTERS and whether it TELESCOPES.

PRE-REGISTRATION (numbers/direction before running; priors stated to lose).
------------------------------------------------------------------
H_GATE4_SANITY: the order-3 model (WITH sigma_2, base=(1+qs+q^2 s1+q^3 s2)) reproduces exact W_3
    on every pair (0 mismatch). [checks the plumbing; near-tautology by digit definition.]
H_SIGMA2_ENTERS (*** structural ***): the order-2 model (sigma_2 DROPPED) DISAGREES with exact
    W_3 on sigma_2!=0 primes. PRED: TRUE -> the cascade genuinely needs a 3rd digit at level 4.
    FALSIFIER: 0 disagreement -> the chain closes at 2 digits (omega needs only s, sigma_1);
    committed to report as such.
H_TELESCOPE (*** the payoff ***): the disagreement W_3_true - W_3_model == y1*j1*sigma_2 mod q --
    the SAME coefficient sigma_1 carried at level 3. PRED: TRUE -> digits enter linearly with the
    universal coefficient y1*j1 -> they reassemble into y1*j1*(2^d-1)/q = a GENERATING FUNCTION
    for omega. FALSIFIER: a DIFFERENT coefficient -> no clean telescoping; omega's digit
    dependence is genuinely higher-order and the closed form is harder.
H_SIGMA2_INDEP (structure): (s, sigma_1, sigma_2) over primes; sigma_2 ranges over F_q, vanishes
    independently. Reported.

DECISION RULES (exact-set-equality; cannot be mis-thresholded):
  H_GATE4_SANITY   CONFIRMED iff order-3 mismatches == 0 at every q.
  H_SIGMA2_ENTERS  CONFIRMED iff order-2 mismatches > 0 at sigma_2!=0 primes.
  H_TELESCOPE      CONFIRMED iff (W_3_true - W_3_model) == y1*j1*sigma_2 mod q on ALL pairs.
                   REFUTED iff the coefficient is anything else (report the actual coefficient).

NOT AT STAKE: R10-R37. A refutation of H_GATE4_SANITY kills only the level-4 plumbing.
"""
import random
from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def digits3(q, d):
    """(s, sigma_1, sigma_2) = first three q-adic digits of (2^d - 1)/q."""
    t = (pow(2, d) - 1) // q
    return t % q, (t // q) % q, (t // (q * q)) % q


def cascade_w3(q, a1, b1, a2, b2, a3, b3):
    """Exact W_3 mod q from residues a_m=2^{-S_m}, b_m=2^{-S'_m} (ints). Returns (W3, ok)."""
    T1 = a1 - b1
    if T1 % q != 0:
        return None, False
    U1 = T1 // q
    S = U1 + (a2 - b2)
    if S % q != 0:
        return None, False
    W2 = S // q
    S2 = W2 + (a3 - b3)
    if S2 % q != 0:
        return None, False
    return (S2 // q) % q, True


def run(q, nsamp, seed):
    rng = random.Random(seed)
    d = order_of_two(q)
    N = q ** 4
    inv2 = pow(2, -1, N)
    inv2q = pow(2, -1, q)
    s, s1, s2 = digits3(q, d)
    dq3 = d * q ** 3
    H = {pow(inv2q, i, q): i for i in range(d)}   # 2^{-i} mod q -> i
    # model powers of 2^{-j1 d}: order-2 (drop sigma_2) and order-3 (with sigma_2, == true)
    npair = bad_sanity = bad_order2 = pred_order2 = tele_bad = 0
    tries = 0
    while npair < nsamp and tries < nsamp * 60:
        tries += 1
        S1 = rng.randrange(1, dq3 + 1)
        S2 = rng.randrange(1, dq3 + 1)
        S3 = rng.randrange(1, dq3 + 1)
        j1 = rng.randrange(1, q)
        S1p = S1 + j1 * d
        a1 = pow(inv2, S1, N); b1 = pow(inv2, S1p, N)
        a2 = pow(inv2, S2, N)
        # level 2: 2^{-S'_2} = 2^{-S_2} + j1 s 2^{-S_1} (mod q)
        t2 = (pow(inv2q, S2 % d, q) + j1 * s * pow(inv2q, S1 % d, q)) % q
        if t2 == 0 or t2 not in H:
            continue
        S2p = H[t2] + d * rng.randrange(0, q ** 3)
        if S2p < 1:
            S2p += d
        b2 = pow(inv2, S2p, N)
        # W_2 (true) to build level-3 target
        w2, ok2 = cascade_w2(q, a1, b1, a2, b2)
        if not ok2:
            continue
        a3 = pow(inv2, S3, N)
        # level 3: W_2 + T_3 = 0 mod q -> 2^{-S'_3} = 2^{-S_3} + W_2 (mod q)
        t3 = (pow(inv2q, S3 % d, q) + w2) % q
        if t3 == 0 or t3 not in H:
            continue
        S3p = H[t3] + d * rng.randrange(0, q ** 3)
        if S3p < 1:
            S3p += d
        b3 = pow(inv2, S3p, N)

        w3_true, ok = cascade_w3(q, a1, b1, a2, b2, a3, b3)
        if not ok:
            continue
        npair += 1
        # model 2^{-j1 d}: order-2 (drop sigma_2), order-3 (with sigma_2)
        base2 = (1 + q * s + q * q * s1) % N
        base3 = (base2 + q ** 3 * s2) % N
        p2_o2 = pow(pow(base2, j1, N), -1, N)
        p2_o3 = pow(pow(base3, j1, N), -1, N)
        b1_o2 = (a1 * p2_o2) % N
        b1_o3 = (a1 * p2_o3) % N
        w3_o2, ok_a = cascade_w3(q, a1, b1_o2, a2, b2, a3, b3)
        w3_o3, ok_b = cascade_w3(q, a1, b1_o3, a2, b2, a3, b3)
        if not (ok_a and ok_b):
            # model broke a division (shouldn't for order-3; order-2 may at deeper digit)
            bad_sanity += 1 if not ok_b else 0
            continue
        if w3_o3 != w3_true:
            bad_sanity += 1
        y1 = a1 % q
        shift = (y1 * j1 * s2) % q
        if w3_o2 != w3_true:
            bad_order2 += 1
        if shift != 0:
            pred_order2 += 1
        if (w3_true - w3_o2) % q != shift:
            tele_bad += 1
    return d, s, s1, s2, npair, bad_sanity, bad_order2, pred_order2, tele_bad


def cascade_w2(q, a1, b1, a2, b2):
    T1 = a1 - b1
    if T1 % q != 0:
        return None, False
    U1 = T1 // q
    S = U1 + (a2 - b2)
    if S % q != 0:
        return None, False
    return (S // q) % q, True


def main():
    log("# PROBE 38 -- level-4 gate: does sigma_2 enter, and do the digit corrections TELESCOPE?")
    log("")

    log("## H_SIGMA2_INDEP -- (s, sigma_1, sigma_2) over primes (digits of (2^d-1)/q)")
    log(f"   {'q':>5} {'d':>4} {'s':>4} {'sig1':>5} {'sig2':>5}  note")
    primes = [p for p in range(3, 130) if all(p % j for j in range(2, int(p ** 0.5) + 1))]
    for q in primes:
        d = order_of_two(q)
        s, s1, s2 = digits3(q, d)
        note = "<- sigma_2 != 0 (3rd digit live)" if s2 != 0 else ""
        if q < 50 or note:
            log(f"   {q:>5} {d:>4} {s:>4} {s1:>5} {s2:>5}  {note}")
    log("")
    log("   (q=5,7 have s1=s2=...=0: single-digit (2^d-1)/q, and both GAP -> sigma_2=0 benign,")
    log("    same witness as sigma_1. R38 tests only ENTERS + TELESCOPES.)")
    log("")

    log("## H_GATE4_SANITY + H_SIGMA2_ENTERS + H_TELESCOPE (exact-iff vs big-int W_3)")
    log("   order-2 model = sigma_2 DROPPED from 2^{-j1 d}; order-3 = full (== true)")
    log("")
    log(f"   {'q':>4} {'d':>3} {'s':>3} {'sig1':>5} {'sig2':>5} {'pairs':>8} {'ord3 bad':>9} "
        f"{'ord2 bad':>9} {'pred ord2':>10} {'tele bad':>9} {'telescope?':>11}")
    ok_sanity = True
    enters = None
    tele = None
    for q in [11, 13, 23, 41]:
        d, s, s1, s2, n, bs, bo2, po2, tb = run(q, 120_000, 3800 + q)
        if bs:
            ok_sanity = False
        telescope = (tb == 0) and (bo2 == po2)
        if s2 != 0:
            if bo2 > 0 and telescope:
                enters = True
                tele = True if tele is None else (tele and True)
            elif bo2 > 0 and not telescope:
                enters = True
                tele = False
        log(f"   {q:>4} {d:>3} {s:>3} {s1:>5} {s2:>5} {n:>8} {bs:>9} {bo2:>9} "
            f"{po2:>10} {tb:>9} {str(telescope):>11}")
    log("")
    log(f"   H_GATE4_SANITY: {'CONFIRMED (order-3 model = exact W_3, 0 mismatch)' if ok_sanity else '*** BROKEN ***'}")
    if enters and tele:
        log("   H_SIGMA2_ENTERS: CONFIRMED -- order-2 (sigma_2-dropped) fails at sigma_2!=0 primes.")
        log("   H_TELESCOPE: *** CONFIRMED *** -- the failure == y1*j1*sigma_2, the SAME coefficient")
        log("      sigma_1 carried at level 3. => digits of (2^d-1)/q enter linearly with the common")
        log("      coefficient y1*j1 => they REASSEMBLE into y1*j1*(2^d-1)/q. A GENERATING FUNCTION")
        log("      for omega's digit-dependence: the whole correction tower is ONE term in the exact")
        log("      tower constant (2^d-1)/q. This is the closed form the L3 bound wants.")
    elif enters and tele is False:
        log("   H_SIGMA2_ENTERS: CONFIRMED, but H_TELESCOPE REFUTED -- sigma_2 enters with a")
        log("      DIFFERENT coefficient. No clean telescoping; omega's digit dependence is genuinely")
        log("      higher-order. Report the actual coefficient; the closed form is harder.")
    elif enters is None:
        log("   H_SIGMA2_ENTERS: NOT confirmed -- order-2 model did not fail; sigma_2 may be inert")
        log("      (chain closes at 2 digits). Inspect table -- would mean omega needs only s, sigma_1.")
    log("")
    log("## READ")
    log("   If telescoping holds: omega's dependence on the (2^d-1)/q tower is LINEAR with a single")
    log("   coefficient -> write it as y1*j1*(2^d-1)/q inside the character sum, and the level-3")
    log("   truncation (R37) is the k=2 shadow of a closed generating function. The L3 bound then")
    log("   targets a clean object -- no infinite tower of digit corrections to track.")
    with open("result_38_level4_telescope_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
