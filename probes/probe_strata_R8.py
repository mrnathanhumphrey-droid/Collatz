"""
PROBE R8 -- UNIFORM KILL, LEDGER WELD, SIGN FRONT, STRATA. Exact rationals.
Reuses R7's engine verbatim (probe_engine_R7): only a measure swap + bookkeeping.

Engine (R7): OffDiag_k = (2/3) sum_{m>=1} 4^{-m} C_k(m),
  C_k(m) = sum_{a,a'} mu_{k-1}(a) mu_{k-1}(a') c_{3^k}(4^{-m}(1+3a) - (1+3a')),   c = Ramanujan sum.
R7's C uses RAW c (c_9 in {6,-3,0}); C_2(1)=-1, C_3(1)=4/49 are the RAW anchors (used in R8-C/D).
R8-A's "=0" gate is normalization-invariant (zero is zero).

Ordering load-bearing: R8-A gates R8-D; run A first, STOP if it fails.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7

# frozen shells S_k (CollatzVerify/Basic.lean)
S = {1: F(2, 3), 2: F(10, 21), 3: F(31370, 67963),
     4: F(143195649659456490, 308468774477179141),
     5: F(2490699741144069281815149277465294323722281911695597204406570,
          5350418720142111510029542161258891403960563152740894082816203)}


def v3(n):
    n = abs(int(n))
    if n == 0:
        return 10 ** 9
    j = 0
    while n % 3 == 0:
        n //= 3; j += 1
    return j


def C_table(k, mu_prev):
    """full C_k(m), m=1..P (P=3^{k-1}), from R7's raw-c engine."""
    M = 3 ** k
    P = 3 ** (k - 1)
    inv4 = pow(4, -1, M)
    return {m: R7.Ck(k, m, mu_prev, pow(inv4, m, M)) for m in range(1, P + 1)}


# ---------------------------------------------------------------- R8-A
def R8A(mu):
    print("## R8-A  THE UNIFORM KILL  (hard gate; PRE-REG: every C_k(m) == 0 exact)")
    ok = True
    for k in range(2, 6):
        Pm1 = 3 ** (k - 1)
        unif = {a: F(1, Pm1) for a in range(Pm1)}          # uniform on Z/3^{k-1}
        T = C_table(k, unif)
        nz = [(m, c) for m, c in T.items() if c != 0]
        status = "ALL ZERO" if not nz else f"NONZERO x{len(nz)}"
        print(f"   k={k} (P={Pm1}):  {status}")
        if nz:
            ok = False
            for m, c in nz[:5]:
                # raw pair counts for the failing cell
                M = 3 ** k; inv4 = pow(4, -1, M); fac = pow(inv4, m, M)
                ge_k = sum(1 for a in range(Pm1) for ap in range(Pm1)
                           if (fac * (1 + 3 * a) - (1 + 3 * ap)) % M == 0)
                eq_km1 = sum(1 for a in range(Pm1) for ap in range(Pm1)
                             if v3(fac * (1 + 3 * a) - (1 + 3 * ap)) == k - 1)
                print(f"      FAIL (k={k},m={m}): C={c}  #v3>=k={ge_k}  #v3=k-1={eq_km1}")
    print(f"   => R8-A {'PASS' if ok else 'FAIL -- STOP (walk-back #31); do NOT trust R8-C/D'}")
    return ok


# ---------------------------------------------------------------- R8-B
def R8B(off):
    print("\n## R8-B  LEDGER-DEVIATION WELD  (PRE-REG: sum_{k=2}^K OffDiag = d_K - 1/5, exact)")
    d = {k: S[k] - F(7, 15) for k in range(1, 6)}
    print(f"   anchors: d_1={d[1]} (=1/5? {d[1]==F(1,5)})  d_2={d[2]} (=1/105? {d[2]==F(1,105)})  "
          f"d_3={d[3]} (=-5191/1019445? {d[3]==F(-5191,1019445)})")
    run = F(0); ok = True
    print(f"   {'K':>2} {'sum OffDiag':>16} {'d_K - 1/5':>16} {'OffDiag_K':>16} {'d_K-d_(K-1)':>16} {'weld':>6}")
    for K in range(2, 6):
        run += off[K]
        lhs = run; rhs = d[K] - F(1, 5)
        perlevel_l = off[K]; perlevel_r = d[K] - d[K - 1]
        good = (lhs == rhs) and (perlevel_l == perlevel_r)
        ok = ok and good
        print(f"   {K:2d} {str(lhs)[:16]:>16} {str(rhs)[:16]:>16} {str(perlevel_l)[:16]:>16} "
              f"{str(perlevel_r)[:16]:>16} {'OK' if good else 'DEV':>6}")
    print(f"   => R8-B {'WELD HOLDS' if ok else 'MISMATCH -- table-provenance error'}")
    return ok


# ---------------------------------------------------------------- R8-C
def R8C(Ctabs):
    print("\n## R8-C  SIGN-FLIP FRONT  (measurement; signed tables verbatim, NO fit, NO threshold detector)")
    sgn = lambda x: '+' if x > 0 else ('-' if x < 0 else '0')
    for k in range(2, 6):
        T = Ctabs[k]; P = 3 ** (k - 1)
        # (i) by raw m
        raw = "".join(sgn(T[m]) for m in range(1, P + 1))
        print(f"\n   k={k} (P={P})  sign by raw m (m=1..P): {raw}")
        if k <= 3:
            print(f"      values: {[str(T[m]) for m in range(1, P+1)]}")
        # (ii) grouped by stratum j = v3(m)
        strat = {}
        for m in range(1, P + 1):
            j = v3(m) if m < P else (k - 1)                # m=P is DC (v3=k-1)
            strat.setdefault(j, []).append((m, T[m]))
        for j in sorted(strat):
            tag = "  [DC]" if j == k - 1 else ""
            signs = "".join(sgn(c) for _, c in strat[j])
            print(f"      stratum j={j}{tag}: signs={signs}  members m={[m for m,_ in strat[j]]}")
    print("\n   anchors: C_2(1)=-1, C_3(1)=+4/49, DC entries positive -> verify above.")


# ---------------------------------------------------------------- R8-D
def R8D1():
    print("\n## R8-D1  BAND COUNTS  (PRE-REG: #v3>=k = 3^{k-1}, #v3=k-1 = 2*3^{k-1}, m-independent)")
    ok = True
    for k in (3, 4, 5):
        M = 3 ** k; Pm1 = 3 ** (k - 1); inv4 = pow(4, -1, M)
        exp_hi, exp_lo = Pm1, 2 * Pm1
        cells = []
        for m in range(1, Pm1 + 1):
            fac = pow(inv4, m, M)
            hi = sum(1 for a in range(Pm1) for ap in range(Pm1)
                     if (fac * (1 + 3 * a) - (1 + 3 * ap)) % M == 0)
            lo = sum(1 for a in range(Pm1) for ap in range(Pm1)
                     if v3(fac * (1 + 3 * a) - (1 + 3 * ap)) == k - 1)
            cells.append((hi, lo))
        hi_ok = all(h == exp_hi for h, _ in cells)
        lo_ok = all(l == exp_lo for _, l in cells)
        uniform = len(set(cells)) == 1
        ok = ok and hi_ok and lo_ok and uniform
        print(f"   k={k}: expected (#v3>=k, #v3=k-1)=({exp_hi},{exp_lo}); "
              f"observed set={sorted(set(cells))}  "
              f"[{'m-INDEP + MATCH' if (hi_ok and lo_ok and uniform) else 'DEV'}]")
    print(f"   => R8-D1 {'PASS' if ok else 'DEV'}")
    return ok


def R8D2(Ctabs):
    print("\n## R8-D2  PALINDROME  C_k(r)=C_k(3^{k-1}-r) for r!=0  (k=4,5; k=2,3 prior)")
    ok = True
    for k in (2, 3, 4, 5):
        T = Ctabs[k]; P = 3 ** (k - 1)
        bad = [r for r in range(1, P) if T[r] != T[P - r]]
        tag = "(prior)" if k <= 3 else ""
        print(f"   k={k} (P={P}): palindrome {'HOLDS' if not bad else f'FAILS at r={bad[:5]}'} {tag}")
        ok = ok and not bad
    print(f"   => R8-D2 {'PASS' if ok else 'DEV'}")
    return ok


def W_closed(j):
    """closed form W_j = sum_{m>=1, v3(m)=j} 4^{-m} = x/(1-x) - x^3/(1-x^3), x=4^{-3^j}."""
    x = F(1, 4 ** (3 ** j))
    return x / (1 - x) - x ** 3 / (1 - x ** 3)


def R8D3(Ctabs):
    print("\n## R8-D3  STRATUM WEIGHTS + AVERAGED CORRELATIONS  (measurement; NO fit)")
    for k in range(2, 6):
        T = Ctabs[k]; P = 3 ** (k - 1)
        geom = 1 - F(1, 4 ** P)                              # 1 - 4^{-P}
        # residue-grouped route: W_j = sum_{r mod P, v3(r)=j} 4^{-r}/(1-4^{-P})
        Wgrp = {}; num = {}                                  # num[j] = sum 4^{-r} C_k(r)
        for r in range(1, P + 1):
            j = v3(r) if r < P else (k - 1)                  # r=P -> DC (j=k-1)
            w = F(1, 4 ** r) / geom
            Wgrp[j] = Wgrp.get(j, F(0)) + w
            num[j] = num.get(j, F(0)) + w * T[r]
        print(f"\n   k={k} (P={P}):  strata j=0..{k-2} bulk, j={k-1}=DC")
        print(f"      {'j':>3} {'W_j (grouped)':>18} {'W_j (closed)':>18} {'match':>6} {'Cbar_k(j)':>18} {'W_j*Cbar':>14}")
        offchk = F(0)
        for j in sorted(Wgrp):
            Wg = Wgrp[j]
            if j <= k - 2:
                Wc = W_closed(j); match = 'OK' if Wg == Wc else 'DEV'; wcs = str(Wc)[:18]
            else:                                             # DC aggregate: 4^{-P}/(1-4^{-P})
                Wc = F(1, 4 ** P) / geom; match = 'OK' if Wg == Wc else 'DEV'; wcs = str(Wc)[:18] + " [DC]"
            Cbar = num[j] / Wg
            contrib = Wg * Cbar
            offchk += contrib
            print(f"      {j:>3} {str(Wg)[:18]:>18} {wcs:>18} {match:>6} {str(Cbar)[:18]:>18} {float(contrib):>+14.3e}")
        print(f"      (2/3)*sum W_j Cbar = {F(2,3)*offchk} = OffDiag_k ? {F(2,3)*offchk == S[k]-S[k-1]}")
        # DC deadness (k>=3): DC weight x DC value
        dcw = F(1, 4 ** P) / geom; dcv = T[P]
        print(f"      DC: weight={float(dcw):.3e} x value={dcv}={float(dcv):.4f} -> {float(dcw*dcv):+.3e}  "
              f"(exact {dcw*dcv})")


def main():
    print("# PROBE R8 -- UNIFORM KILL / LEDGER WELD / SIGN FRONT / STRATA. Exact rationals.\n")
    # build real measure hierarchy (R7)
    mu = {1: R7.mu1()}
    for k in range(2, 5):
        mu[k] = R7.build_mu(mu[k - 1], k)

    a_ok = R8A(mu)                                            # RUN FIRST
    if not a_ok:
        print("\n*** R8-A FAILED -- halting per spec (do not trust downstream). ***")
        return

    # real C-tables (raw c) + OffDiag for B/C/D
    Ctabs = {}; off = {}
    for k in range(2, 6):
        off[k], Ctabs[k] = R7.engine_offdiag(k, mu[k - 1], want_C=True)

    R8B(off)
    R8C(Ctabs)
    R8D1()
    R8D2(Ctabs)
    R8D3(Ctabs)


if __name__ == "__main__":
    main()
