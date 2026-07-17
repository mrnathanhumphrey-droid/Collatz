"""
PROBE 44 -- THE RIG: does the order-reciprocity spine transmit force, or is (2,3) a coincidence?
(Platform brief "Order-Reciprocity Spine". Tests Claim A [see the object] vs Claim B [seeing=landing].)

THE SPINE (observed at (2,3), flagged as observation not theorem in R43):
   ord_8(3)=2 (Chang 2-adic face)  vs  ord_3(2)=2 (Nathan q-adic face), both order-2 of the
   reciprocal prime (2^2-1=3, 3^2-1=8=2^3, the Catalan pair).
H_RIG: there is ONE fixed relation R with (d_q=ord_q(2)) reciprocated by the 2-adic return-order,
   for all odd q, the SAME R that gives 1/4 at q=3. If so, r_q becomes derivable through R (rig
   pulls). If (2,3) is a smallest-prime coincidence, the rig is decorative (tools yank independently).

PRE-REGISTERED (brief's own, record 0-for-9 on 'and therefore'):
   H_RIG MOST LIKELY FAILS ((2,3) special via Catalan). PASS bar (brief S3): the persistence
   constants for q=3,5,7,11,13 lie on ONE monotone curve in d_q, no free params, hitting 1/4 at q=3.
   Scatter / flat / non-existence = FAIL.

STRUCTURAL PRE-CHECK (this probe, cheap arithmetic):
   (A) (Z/2^s)^* = Z/2 x Z/2^{s-2}: ALL element orders are POWERS OF 2. So ord_{2^s}(q) is always
       a power of 2, for every odd q. But d_q=ord_q(2) is generically NOT a power of 2.
       => order-matching reciprocity ord_{2^{s_q}}(q)=d_q is IMPOSSIBLE when d_q is not a power of 2.
   (B) ord_8(q): since q^2 ≡ 1 (mod 8) for EVERY odd q (all odd squares are 1 mod 8), ord_8(q) in
       {1,2} for all odd q (=2 unless q≡1 mod 8). So 'ord_8(3)=2' is GENERIC, not q=3-special.
   (C) the mod-8 persistence constant Pr[persistent]_q = Sum_k 2^{-k} * frac(persistent mu | k):
       since q^k mod 8 is period<=2 for all odd q, this is 1/4 for EVERY odd q (q not 1 mod 8) --
       FLAT in d_q, does not track it.

TEST: compute (A),(B),(C) for q=3,5,7,11,13,17,19; adjudicate H_RIG on the brief's bar.

NOT AT STAKE: R1-R43. This tests whether the order-reciprocity is a LAW or a (2,3) coincidence.
"""
from fractions import Fraction

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def ordmod(a, n):
    if n == 1:
        return 1
    x, o = a % n, 1
    while x != 1:
        x = (x * a) % n; o += 1
        if o > 4 * n:
            return -1
    return o


def is_pow2(x):
    return x >= 1 and (x & (x - 1)) == 0


def main():
    log("# PROBE 44 -- the order-reciprocity SPINE: law or (2,3) coincidence?")
    log("")

    QS = [3, 5, 7, 11, 13, 17, 19]

    # (A)+(reciprocity existence): does exists s_q with ord_{2^{s_q}}(q) = d_q = ord_q(2)?
    log("## (A) reciprocity ord_{2^{s_q}}(q) = d_q = ord_q(2)?  (2-adic orders are ALWAYS powers of 2)")
    log(f"   {'q':>4} {'d_q=ord_q(2)':>13} {'d_q pow2?':>10} {'ord_{2^s}(q) for s=2..6':>28} {'s_q exists?':>12}")
    recip_ok = {}
    for q in QS:
        dq = ordmod(2, q)
        orders2 = [ordmod(q, 2 ** s) for s in range(2, 7)]
        sq = None
        for s in range(2, 12):
            if ordmod(q, 2 ** s) == dq:
                sq = s; break
        recip_ok[q] = (sq is not None)
        log(f"   {q:>4} {dq:>13} {str(is_pow2(dq)):>10} {str(orders2):>28} "
            f"{('s='+str(sq)) if sq else 'NONE (impossible)':>12}")
    log("   => 2-adic orders are powers of 2; d_q reciprocable ONLY when d_q is a power of 2.")
    log(f"      d_q power-of-2 at q={[q for q in QS if is_pow2(ordmod(2,q))]}; "
        f"IMPOSSIBLE at q={[q for q in QS if not is_pow2(ordmod(2,q))]}.")
    log("")

    # (B) ord_8(q) is generic
    log("## (B) is 'ord_8(3)=2' special to q=3?  (q^2 mod 8 for odd q)")
    log(f"   {'q':>4} {'q mod 8':>8} {'q^2 mod 8':>10} {'ord_8(q)':>9}")
    for q in QS:
        log(f"   {q:>4} {q % 8:>8} {(q * q) % 8:>10} {ordmod(q, 8):>9}")
    log("   => q^2 ≡ 1 (mod 8) for EVERY odd q (all odd squares are 1 mod 8). ord_8(q)=2 unless")
    log("      q≡1 mod 8. So ord_8(3)=2 is GENERIC (every odd q≢1 mod8), NOT a q=3 reciprocal fact.")
    log("")

    # (C) the mod-8 persistence constant per q -- flat?
    log("## (C) mod-8 persistence constant Pr[persistent]_q = Sum_{k>=1} 2^{-k} * frac(mu: q^k mu ≡ 7 mod 8)")
    log(f"   {'q':>4} {'d_q':>4} {'Pr[persistent]_q':>17}  (tracks d_q? or flat 1/4?)")
    consts = {}
    for q in QS:
        dq = ordmod(2, q)
        # sum over k>=1 of 2^-k * (#{mu in {1,3,5,7}: q^k mu ≡7 mod8}/4).  q^k mod8 period<=2.
        # closed: for k with q^k invertible mod8 (always), exactly one mu works => frac=1/4 each k.
        # so constant = (sum_k 2^-k)*(1/4) = 1/4 -- unless q≡1 mod8 (then q^k≡1, need mu≡7 => still one).
        tot = Fraction(0)
        for k in range(1, 40):
            a = pow(q, k, 8)
            npers = sum(1 for mu in (1, 3, 5, 7) if (a * mu) % 8 == 7)
            tot += Fraction(1, 2 ** k) * Fraction(npers, 4)
        consts[q] = tot
        log(f"   {q:>4} {dq:>4} {str(tot):>17}")
    vals = set(consts.values())
    log(f"   => constants = {sorted(set(str(v) for v in consts.values()))}  "
        f"{'ALL EQUAL 1/4 -> FLAT in d_q (does NOT track it)' if vals=={Fraction(1,4)} else 'vary'}")
    log("")

    # ---- VERDICT ----
    log("## VERDICT on H_RIG:")
    all_recip = all(recip_ok.values())
    flat = (vals == {Fraction(1, 4)})
    log(f"   (A) order-reciprocity exists for ALL q? {all_recip}  "
        f"(fails at non-power-of-2 d_q: q=7,11,13,19)")
    log(f"   (B) ord_8(q)=2 generic (all odd q≢1 mod8)? YES -- not a q=3 reciprocal fact.")
    log(f"   (C) persistence constant tracks d_q? {'NO -- FLAT at 1/4 for all q' if flat else 'maybe'}")
    log("")
    log("   *** H_RIG FAILS (as pre-registered), with a SHARPER reason than 'small-prime scatter': ***")
    log("   (Z/2^s)^* has ONLY power-of-2 element orders, so the 2-adic side cannot reciprocate")
    log("   d_q=ord_q(2) whenever d_q is not a power of 2 -- provably impossible at q=7 (d=3),")
    log("   q=11 (d=10), q=13 (d=12), q=19 (d=18). And the mod-8 persistence constant is a FLAT 1/4")
    log("   for every odd q (the generic odd^2≡1 mod8 fact), NOT a function of d_q.")
    log("   The (2,3) [and (2,5),(2,17)] matches are because d_q=2,4,8 happen to be POWERS OF 2")
    log("   (reinforced at (2,3) by the Catalan pair 2^2-1=3, 3^2-1=8). NOT a reciprocity law.")
    log("")
    log("   => Claim A TRUE (one object, four faces, Tao's measure -- real, in the files).")
    log("      Claim B FALSE (the spine does NOT transmit force; the order-reciprocity is a")
    log("      coincidence of powers-of-2, not a coupling). The four tools are aimed at one")
    log("      animal; their forces do NOT sum through this spine. The mammoth keeps the sky.")
    log("      0-for-9 -> 0-for-10 on 'and therefore'. The machine working.")
    with open("result_44_reciprocity_spine_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
