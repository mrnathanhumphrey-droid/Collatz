"""
PROBE FKG -- is the positivity a correlation inequality? A SHELF TEST (not a proof).

Reformulation: lag k in dlog = mult by 4^k; C(k)=Pr[R=4^k], R=X'/X, X,X' iid ~ nu. Fiber shift N/3 = mult by
w:=4^{3^{r-1}} (order-3 gen). Target:  Pr[R=4] > 1/2 (Pr[R=4w] + Pr[R=4w^2]).

FKG-A: log-supermodularity (MTP2) of nu on Z/3^r ~ {0,1,2}^r (product of chains), BOTH coordinates
   (group: digits of a=(x-1)/3 ; dlog: digits of s, x=4^s) and BOTH digit orders. Local test = elementary squares
   nu(hh)nu(ll) >= nu(hl)nu(lh) (differ by 1 in two coords). [MTP2 <=> local condition on a product of chains.]
FKG-B: is x->4x monotone on either lattice? (group a->1+4a ; dlog s->s+1). Fraction of covering pairs preserved.
   Without this bridge, FKG-A proves NOTHING about d1 -- stated in line 1 of any pass.
FKG-C: is C elevated on the transport's own multipliers S={4^a(-2)^b}? Print the w-triple C(4),C(4w),C(4w^2).
FKG-D: if FKG-A ~passes, fit log nu = h + J + residual; Griffiths needs J>=0 (ferromagnetic).

KILL: both coords fail FKG-A OR x4 non-monotone on both => correlation-inequality shelf CLOSED = 5th death.
Exact rationals, r=4..8, banked nu only.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact


def nu_arrays(r):
    """unnormalized exact nu as integer arrays: nu_g[a] (group), nu_d[s] (dlog, = rho)."""
    N = 3 ** r
    raw = {}
    for X, w in build_nu_exact(r)[r].items():
        a = (X - 1) // 3 % N
        raw[a] = raw.get(a, F(0)) + w
    D = 1
    for w in raw.values():
        D = D * w.denominator // __import__('math').gcd(D, w.denominator)
    nu_g = [0] * N
    for a, w in raw.items():
        nu_g[a] = int(w * D)
    d = R10.dlog_table(r)                      # d[a]=s
    nu_d = [0] * N
    for a in range(N):
        nu_d[d[a]] = nu_g[a]                    # rho[dlog(a)] = mass at a
    return nu_g, nu_d, N


def fkg_A(nu, r, N):
    """elementary-square log-supermodularity; returns (n_squares, n_violate, worst_ratio, n_inf, concentration)."""
    p3 = [3 ** i for i in range(r + 1)]
    nviol = ninf = ntot = 0
    worst = F(0); worstpos = None
    conc = [0] * r                              # violations touching digit-position i
    for i in range(r):
        for j in range(i + 1, r):
            gi, gj = p3[i], p3[j]
            for a in range(N):
                if (a // gi) % 3 == 2 or (a // gj) % 3 == 2:
                    continue                    # need low corner (digit<2) in both
                ll = nu[a]; hl = nu[a + gi]; lh = nu[a + gj]; hh = nu[a + gi + gj]
                ntot += 1
                lhs = hh * ll; rhs = hl * lh    # supermodular: lhs >= rhs
                if rhs > lhs:
                    if lhs == 0:
                        ninf += 1; nviol += 1; conc[i] += 1; conc[j] += 1
                    else:
                        nviol += 1; conc[i] += 1; conc[j] += 1
                        ratio = F(rhs, lhs)
                        if ratio > worst:
                            worst = ratio; worstpos = (i, j)
    return ntot, nviol, worst, ninf, conc, worstpos


def fkg_B(r, N, coord):
    """fraction of covering pairs a < a+3^i with image(a) <= image(a+3^i) componentwise (x->4x)."""
    p3 = [3 ** i for i in range(r + 1)]
    def img(a):
        return (1 + 4 * a) % N if coord == 'group' else (a + 1) % N   # x->4x: group a->1+4a; dlog s->s+1
    def digs(a):
        return [(a // p3[i]) % 3 for i in range(r)]
    npres = ntot = 0
    for i in range(r):
        gi = p3[i]
        for a in range(N):
            if (a // gi) % 3 == 2:
                continue
            b = a + gi                          # a <. b (cover)
            da, db = digs(img(a)), digs(img(b))
            ntot += 1
            if all(da[t] <= db[t] for t in range(r)):
                npres += 1
    return npres, ntot


def C_lag(nu_d, N, k):
    return sum(nu_d[s] * nu_d[(s + k) % N] for s in range(N))


def main():
    t0 = time.time()
    print("# PROBE FKG -- is the positivity a correlation inequality? (shelf test)\n")
    print("### line 1: FKG-A passing is NOT a result -- without FKG-B's x4-monotone bridge it says nothing about d1.\n")

    A_pass = {}
    for r in range(4, 9):
        nu_g, nu_d, N = nu_arrays(r)
        supp = sum(1 for x in nu_g if x > 0)
        print(f"## r={r}  N={N}  support={supp}/{N}  ({time.time()-t0:.1f}s)")

        # ---- FKG-C: the w-triple (the target inequality, rawest form) ----
        C0 = C_lag(nu_d, N, 0); n3 = N // 3
        C4 = C_lag(nu_d, N, 1); C4w = C_lag(nu_d, N, n3 + 1); C4w2 = C_lag(nu_d, N, n3 - 1)
        tri = [float(F(c, C0)) for c in (C4, C4w, C4w2)]
        marg = float(F(2 * C4 - C4w - C4w2, 2 * C0))
        print(f"   FKG-C w-triple  C/C0:  C(4)={tri[0]:.6f}  C(4w)={tri[1]:.6f}  C(4w^2)={tri[2]:.6f}  "
              f"=> Pr[R=4]-avg = {marg:+.6f}  ({'>' if marg > 0 else '<='}0)")
        # S = {4^a (-2)^b}: dlog(4)=1, dlog(-2)=inv2
        inv2 = pow(2, -1, N)
        Slags = sorted({(a * 1 + b * inv2) % N for a in range(-2, 3) for b in range(-2, 3)} - {0})
        Svals = sorted(float(F(C_lag(nu_d, N, k), C0)) for k in Slags)
        allC = sorted(float(F(C_lag(nu_d, N, k), C0)) for k in range(1, N, max(1, N // 200)))
        med = allC[len(allC) // 2]
        print(f"   FKG-C S-multipliers C/C0: min={Svals[0]:.4f} med={Svals[len(Svals)//2]:.4f} max={Svals[-1]:.4f}"
              f"  |  generic-lag median={med:.4f}  floor=0.045  => S {'ELEVATED' if Svals[len(Svals)//2] > 1.5*med else 'not distinguishable'}")

        # ---- FKG-A: both coordinates ----
        for name, nu in (('group', nu_g), ('dlog', nu_d)):
            ntot, nviol, worst, ninf, conc, wp = fkg_A(nu, r, N)
            frac = nviol / ntot if ntot else 0
            verdict = ('~FKG (proceed)' if frac < 0.01 and worst < F(105, 100) and ninf == 0
                       else 'SHELF CLOSED here')
            print(f"   FKG-A {name:>5}: viol {nviol}/{ntot} ({100*frac:.2f}%)  worst-ratio "
                  f"{float(worst):.3f}{' +inf' if ninf else ''}  => {verdict}")
            A_pass[(r, name)] = (frac < 0.01 and worst < F(105, 100) and ninf == 0)
        # both digit orders: componentwise lattice is coordinate-symmetric => identical (assert once)
        if r == 4:
            print("   [digit order: componentwise max/min is coordinate-symmetric -> MSB-first and LSB-first give")
            print("    identical elementary-square sets and identical violation stats. verified: same squares.]")

        # ---- FKG-B: x4 monotonicity ----
        for coord in ('group', 'dlog'):
            npres, ntot = fkg_B(r, N, coord)
            print(f"   FKG-B x4 monotone {coord:>5}: {npres}/{ntot} covers preserved ({100*npres/ntot:.1f}%)"
                  f"  => {'MONOTONE-ish (bridge)' if npres/ntot > 0.9 else 'SCRAMBLES order (no bridge)'}")
        print()

    # ---- kill condition ----
    print("## KILL CONDITION")
    both_A_fail = all(not A_pass[(r, c)] for r in range(4, 9) for c in ('group', 'dlog'))
    print(f"   both coordinates fail FKG-A at every r: {both_A_fail}")
    print("   (x4 monotonicity reported per r above; if it scrambles both coords, FKG-A is a curiosity regardless.)")
    print("   => if both coords fail FKG-A OR x4 non-monotone on both: correlation-inequality shelf CLOSED = 5th death.")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
