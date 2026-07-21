"""
PROBE R7 -- THE CHANNEL ENGINE (the <4>-orbit law). Exact rationals; the closed-form off-diagonal engine.

  OffDiag_k = (2/3) sum_{m>=1} 4^{-m} C_k(m),      C_k periodic in m mod 3^{k-1} = ord_{3^k}(4).
  C_k(m) = sum_{a,a'} mu_{k-1}(a) mu_{k-1}(a') * c_{3^k}( 4^{-m}*(1+3a) - (1+3a') )
         = the twisted <4^m>-orbit correlation character sum of the level-(k-1) measure over the
           primitive shell (c = Ramanujan sum mod 3^k).
  Class assembly:  OffDiag_k = (2/3) sum_{r=1}^{P} [ 4^{-r}/(1-4^{-P}) ] C_k(r),   P = 3^{k-1}.

Derivation of the channel weight (2/3)4^{-m} (even gap g=2m, both orderings):
  sum_{v'>=1} 2 * 2^{-(v'+2m)} 2^{-v'} = 2 * 2^{-2m} * (1/3) = (2/3) 4^{-m}.
The outer valuation-difference v-v' is the gap; the inner unit factor 2^{-v'} is a unit mod 3^k, so
c_{3^k}(2^{-v'} delta) = c_{3^k}(delta) (Ramanujan sum depends only on gcd) -> C_k depends on m alone.

GATES (pre-registered, forced):
  R7-A  C_2 = (-1,-1,+2) at m=(1,2,0 mod3);  engine OffDiag_2 = -4/21.
        engine OffDiag_3 = -2980/203889 (FROZEN increment, from the engine not the S-table);
        single m=1 channel (2/3)4^{-1} C_3(1) = +2/147  (R6-C's gap-2 sign-flip, DERIVED from C_3(1)).
  R7-B  C-tables C_k(m) exact, k=2..5.
  R7-C  running ledger sum OffDiag_k vs -1/5, engine vs frozen S_k - S_{k-1}.
  (guard) odd-gap channels vanish identically (conjugate-kill) -- verified, not assumed.
"""
from fractions import Fraction as F

# frozen shells S_k (CollatzVerify/Basic.lean); OffDiag_k = S_k - S_{k-1}
S = {1: F(2, 3), 2: F(10, 21), 3: F(31370, 67963),
     4: F(143195649659456490, 308468774477179141),
     5: F(2490699741144069281815149277465294323722281911695597204406570,
          5350418720142111510029542161258891403960563152740894082816203)}


def mu1():
    return {1: F(1, 3), 2: F(2, 3)}                 # Judge One base, mod 3


def build_mu(mu_prev, k):
    """mu_k on Z/3^k via the renewal X = 2^{-v}(1+3 a), a ~ mu_{k-1}, v~Geom(1/2). Exact (folded valuations)."""
    M = 3 ** k
    inv2 = pow(2, -1, M)
    ordv, x = 1, 2 % M
    while x != 1:
        x = (x * 2) % M; ordv += 1                  # multiplicative order of 2 mod 3^k
    denom = 1 - F(1, 2 ** ordv)
    mu = {}
    for j in range(1, ordv + 1):
        wv = F(1, 2 ** j) / denom                   # sum_{v==j mod ord} 2^{-v}
        u = pow(inv2, j, M)                         # 2^{-v} mod 3^k for class j
        for a, pa in mu_prev.items():
            r = (u * (1 + 3 * a)) % M
            mu[r] = mu.get(r, F(0)) + wv * pa
    return mu


def cram(k, n):
    """Ramanujan sum c_{3^k}(n): 3^k|n -> 2*3^{k-1};  3^{k-1}||n -> -3^{k-1};  else 0."""
    M = 3 ** k
    n %= M
    if n == 0:
        return 3 ** k - 3 ** (k - 1)
    if n % (3 ** (k - 1)) == 0:
        return -3 ** (k - 1)
    return 0


def Ck(k, m, mu_prev, factor):
    """twisted <4^m>-orbit correlation character sum. factor = 4^{-m} mod 3^k (rotation), or
       for an ODD-gap channel g=2m+1 the caller passes 2^{-g} mod 3^k instead."""
    M = 3 ** k
    tot = F(0)
    lift = [((1 + 3 * a) % M, pa) for a, pa in mu_prev.items()]
    for A, pa in lift:
        for Ap, pap in lift:
            cr = cram(k, factor * A - Ap)
            if cr:
                tot += pa * pap * cr
    return tot


def engine_offdiag(k, mu_prev, want_C=False):
    """OffDiag_k assembled from the engine; returns (value, C-table dict m->C_k(m))."""
    M = 3 ** k
    P = 3 ** (k - 1)                                # ord_{3^k}(4)
    inv4 = pow(4, -1, M)
    s = F(0); C = {}
    for r in range(1, P + 1):
        fac = pow(inv4, r, M)                       # 4^{-r} mod 3^k
        C[r] = Ck(k, r, mu_prev, fac)
        classw = F(1, 4 ** r) / (1 - F(1, 4 ** P))  # sum_{m==r mod P} 4^{-m}
        s += classw * C[r]
    val = F(2, 3) * s
    return (val, C) if want_C else val


def odd_gap_check(k, mu_prev, mmax=4):
    """verify odd-gap channels g=2m+1 have vanishing character sum (conjugate-kill)."""
    M = 3 ** k
    inv2 = pow(2, -1, M)
    vals = []
    for m in range(0, mmax):
        g = 2 * m + 1
        fac = pow(inv2, g, M)                       # 2^{-g} mod 3^k
        vals.append(Ck(k, m, mu_prev, fac))
    return vals


def main():
    print("# PROBE R7 -- THE CHANNEL ENGINE (the <4>-orbit law). Exact rationals.\n")

    # build the measure hierarchy mu_1..mu_4 (mu_{k-1} feeds level k)
    mu = {1: mu1()}
    for k in range(2, 5):
        mu[k] = build_mu(mu[k - 1], k)
    for k in range(1, 5):
        print(f"   |support(mu_{k})| = {len(mu[k])}  (mod 3^{k}={3**k})")
    # weld: mu_2 must reproduce R6's pi2 total mass 1 and S-consistency is downstream
    print(f"   sum mu_4 = {sum(mu[4].values())} (=1 exact)\n")

    print("## R7-A  THE ENGINE AT k=2 and k=3 (gate)")
    v2, C2 = engine_offdiag(2, mu[1], want_C=True)
    print(f"   C_2(m), m=1..3 : {[C2[r] for r in (1,2,3)]}   (PRE-REG (-1,-1,+2); "
          f"{'OK' if [C2[1],C2[2],C2[3]]==[F(-1),F(-1),F(2)] else 'DEV'})")
    print(f"   engine OffDiag_2 = {v2}  (PRE-REG -4/21; {'GATE PASS' if v2==F(-4,21) else 'DEV'})")
    print(f"   frozen  OffDiag_2 = {S[2]-S[1]}\n")

    v3, C3 = engine_offdiag(3, mu[2], want_C=True)
    ch1_k3 = F(2, 3) * F(1, 4) * C3[1]              # single m=1 channel
    print(f"   C_3(1) = {C3[1]} = {float(C3[1]):+.6f}")
    print(f"   single m=1 channel (2/3)4^-1 C_3(1) = {ch1_k3}  (PRE-REG R6-C gap-2 = +2/147; "
          f"{'GATE PASS -> sign-flip DERIVED' if ch1_k3==F(2,147) else 'DEV'})")
    print(f"   engine OffDiag_3 = {v3} = {float(v3):+.8f}")
    print(f"   frozen OffDiag_3 = {S[3]-S[2]} = {float(S[3]-S[2]):+.8f}   "
          f"({'GATE PASS' if v3==S[3]-S[2] else 'DEV'})\n")

    print("## R7-B  THE C-TABLES  C_k(m) exact, k=2..5  (period 3^{k-1})")
    Cfull = {2: C2, 3: C3}
    off = {2: v2, 3: v3}
    for k in (4, 5):
        off[k], Cfull[k] = engine_offdiag(k, mu[k - 1], want_C=True)
    for k in (2, 3, 4, 5):
        Ck_tab = Cfull[k]; P = 3 ** (k - 1)
        show = P if P <= 9 else 9
        row = "  ".join(f"C({r})={Ck_tab[r]}" for r in range(1, show + 1))
        tail = "" if P <= 9 else f"  ... (P={P} classes)"
        print(f"   k={k} (P={P}): {row}{tail}")
    print()

    print("## R7-C  RUNNING LEDGER  (engine vs frozen; target sum = -1/5)")
    print(f"   {'k':>2} {'engine OffDiag_k':>22} {'float':>12} {'frozen':>12} {'match':>6} {'running':>12} {'vs -1/5':>10}")
    run = F(0)
    for k in (2, 3, 4, 5):
        fr = S[k] - S[k - 1]; run += off[k]
        match = 'OK' if off[k] == fr else 'DEV'
        print(f"   {k:2d} {str(off[k])[:22]:>22} {float(off[k]):>+12.6f} {float(fr):>+12.6f} {match:>6} "
              f"{float(run):>+12.7f} {float(run-F(-1,5)):>+10.6f}")
    print(f"   sum_(k>=2) OffDiag target = S_inf - S_1 = 7/15 - 2/3 = {F(7,15)-F(2,3)} = -1/5")
    print(f"   sum_(k>=3) = 7/15 - 10/21 = {F(7,15)-F(10,21)} = -1/105\n")

    print("## GUARD  odd-gap channels vanish (conjugate-kill; NOT assumed)")
    for k in (2, 3):
        ov = odd_gap_check(k, mu[k - 1])
        allz = all(x == 0 for x in ov)
        print(f"   k={k}: odd-gap C (g=1,3,5,7) = {ov}   ({'ALL ZERO' if allz else 'NONZERO!'})")


if __name__ == "__main__":
    main()
