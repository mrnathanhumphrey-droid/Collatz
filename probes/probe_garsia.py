"""
PROBE GARSIA -- the arithmetic type of nu via support-count + Garsia ENTROPY dimension (2026-07-27)

Wilson's two-thought spec. nu_r = law of X_r under the certified renewal X' = 1 + 3*2^{-v} X
on Z/3^{r+1} (build_nu, from probe_gapop_R28). Decides singular-vs-a.c. TYPE of nu -> whether
S_inf can be rational at all (7/15 "wrong in kind" iff nu singular).

The p-adic Garsia dictionary (Wilson thought 2): local dimension = log(reach per level)/log 3.
  - SUPPORT count  s_r = |supp(nu_r)|.  s_r ~ beta^r  => box-dim = log beta / log 3.
  - ENTROPY        H_r = -sum w log w.   h = lim H_r/(r log 3) = Garsia entropy dimension.
    Streck: a.c. with Fourier decay  <=>  entropy MAXIMAL (h=1).  h<1 => singular.

*** Two DIFFERENT dimensions. Support fatness is necessary, not sufficient, for a.c. ***
The MEASURE decides via ENTROPY (mass scaling), not support (point count).

PART 1  SUPPORT -- exact, algebraic (2 is a primitive root mod 3^r):
        supp(nu_r) = {1+3t : t in (Z/3^r)^*}  => s_r = 2*3^{r-1}, beta=3 EXACTLY, box-dim 1.
        Verified: ord_{3^r}(2)=2*3^{r-1} (r=1..15); full-group closure = full unit coset (r=1..8).
        WARNING gated: len(build_nu, tol=1e-18) < s_r for r>=5 (tol walks only ~60 multipliers) --
        so naive len() FAKES a singular answer. Do NOT use len(build_nu) as the support.

PART 2  ENTROPY -- the decisive instrument. H_r from build_nu weights, r=1..RMAX.
        h_r = H_r/(r log3); per-level dH_r/log3 = (H_r-H_{r-1})/log3 (Garsia entropy dim, incremental).
        tol robustness at r=10,12 (mass-weighted => tol-safe even though support is not).
        Heuristic: per level adds H(v)=entropy of Geom(1/2)=2log2=1.386 nats vs scale log3=1.099
        => H(step)/log3=1.26>1 => predicts SATURATION h=1 (a.c.-type). Measurement adjudicates.

PART 3  MAHLER cross-check (Breuillard-Varju), light: base lambda=2, Mahler M(2)=2.
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_gapop_R28 import build_nu

L3 = math.log(3); L2 = math.log(2)


def ord_(a, m):
    x = a % m; k = 1
    while x != 1:
        x = (x * a) % m; k += 1
        if k > 4 * m: return -1
    return k


def full_group_support(r):
    """EXACT support of nu_r: level-by-level pushforward with the FULL multiplier group <2^{-1}>,
       no tol, no weights. Returns the support set on Z/3^{r+1}."""
    supp = {1}                                            # nu_0 = delta_1 on Z/3
    for lev in range(1, r + 1):
        Ml = 3 ** lev; Mlp = 3 ** (lev + 1); i2 = pow(2, -1, Ml)
        grp = set(); u = 1
        for _ in range(2 * 3 ** (lev - 1)):               # one full cycle of <i2> (order 2*3^{lev-1})
            u = (u * i2) % Ml; grp.add(u)
        supp = {(1 + 3 * ((u * Xp) % Ml)) % Mlp for Xp in supp for u in grp}
    return supp


def entropy(nu):
    w = np.fromiter(nu.values(), float)
    w = w / w.sum()
    return float(-(w * np.log(w)).sum())


def main():
    print("# PROBE GARSIA -- arithmetic type of nu via support + entropy dimension\n")

    # =============================== PART 1: SUPPORT (exact) ===============================
    print("## PART 1  SUPPORT -- exact / algebraic")
    print("  (a) is 2 a primitive root mod 3^r?  ord_{3^r}(2) should = phi(3^r) = 2*3^{r-1}")
    okpr = True
    for r in range(1, 16):
        o = ord_(2, 3 ** r); want = 2 * 3 ** (r - 1); ok = (o == want)
        okpr = okpr and ok
        if r <= 6 or not ok:
            print(f"     r={r:2d}: ord_{{3^{r}}}(2)={o}  want 2*3^{r-1}={want}  {'ok' if ok else 'FAIL'}")
    print(f"     => 2 is a primitive root mod 3^r for r=1..15: {okpr}")
    print("  (b) full-group closure support  vs  the full unit coset {1+3t: gcd(t,3)=1} (size 2*3^{r-1})")
    okcoset = True
    for r in range(1, 9):
        S = full_group_support(r)
        want_set = {(1 + 3 * t) % 3 ** (r + 1) for t in range(3 ** r) if t % 3 != 0}
        ok = (S == want_set); okcoset = okcoset and ok
        print(f"     r={r}: |supp|={len(S):6d}  2*3^{r-1}={2*3**(r-1):6d}  == full unit coset: {ok}")
    print(f"     => supp(nu_r) = full unit coset, s_r = 2*3^{{r-1}} EXACTLY, beta=3, box-dim=1: {okcoset}")
    print("  (c) TOL-TRAP: len(build_nu, tol=1e-18) vs true s_r=2*3^{r-1} -- naive len UNDERCOUNTS for r>=5")
    nus18 = build_nu(0.5, 12, tol=1e-18)
    for r in range(1, 13):
        s_true = 2 * 3 ** (r - 1); s_len = len(nus18[r])
        flag = "" if s_len == s_true else "  <-- tol undercount (would fake singular!)"
        print(f"     r={r:2d}: len(build_nu)={s_len:8d}   true s_r={s_true:8d}{flag}")
    print()

    # =============================== PART 2: ENTROPY (decisive) ===============================
    print("## PART 2  ENTROPY -- the Garsia entropy dimension (THE instrument)")
    RMAX = 14
    print(f"  building nu_r, r=1..{RMAX}, tol=1e-18 (entropy is mass-weighted => tol-safe) ...")
    nus = build_nu(0.5, RMAX, tol=1e-18)
    H = {r: entropy(nus[r]) for r in range(1, RMAX + 1)}
    Hmax = {r: math.log(2 * 3 ** (r - 1)) for r in range(1, RMAX + 1)}   # log|supp| = max possible H
    print(f"  {'r':>2} {'H_r(nats)':>11} {'h_r=H/(r*L3)':>13} {'dH_r/L3':>9} {'H_r/log|supp|':>14} {'|supp|':>9}")
    prev = 0.0
    for r in range(1, RMAX + 1):
        h_r = H[r] / (r * L3)
        dH = (H[r] - prev) / L3
        ratio_full = H[r] / Hmax[r]
        print(f"  {r:>2} {H[r]:>11.5f} {h_r:>13.5f} {dH:>9.5f} {ratio_full:>14.5f} {2*3**(r-1):>9d}")
        prev = H[r]
    # incremental dimension: robust late-window slope of H_r vs r, in units of log3
    rr = np.arange(8, RMAX + 1)
    slope = np.polyfit(rr, [H[r] for r in rr], 1)[0] / L3
    print(f"  => late-window (r=8..{RMAX}) entropy slope dH/dr / log3 = {slope:.5f}"
          f"   (=1 -> maximal/a.c.-type [Streck];  <1 -> singular, dim<1)")

    print("\n  tol robustness of ENTROPY (r=10,12): tol=1e-12,1e-18,1e-24 should agree to ~1e-6")
    for r in (10, 12):
        row = []
        for tol in (1e-12, 1e-18, 1e-24):
            nr = build_nu(0.5, r, tol=tol)[r]
            row.append(f"tol{tol:.0e}:H={entropy(nr):.6f}(|s|={len(nr)})")
        print(f"     r={r}: " + "   ".join(row))

    # =============================== PART 3: MAHLER cross-check ===============================
    print("\n## PART 3  Mahler / Breuillard-Varju cross-check (light)")
    Hv = 2 * L2                      # entropy of the step multiplier v ~ Geom(1/2): H = 2 log 2
    print(f"  per-level step entropy H(v)=H(Geom 1/2)=2 log2 = {Hv:.5f} nats;  scale step = log3 = {L3:.5f}")
    print(f"  H(step)/log3 = {Hv/L3:.5f}  (>1 => entropy production outruns the scale => saturates at h=1)")
    print(f"  base lambda=2 (algebraic integer, min poly x-2): Mahler M(2)=2, no conjugates =>")
    print(f"    Breuillard-Varju: no conjugate on/inside forcing entropy defect; entropy defect (if any)")
    print(f"    comes from OVERLAP of the 3-adic digit maps, not from the base. Read h_r vs 1 above.\n")


if __name__ == "__main__":
    main()
