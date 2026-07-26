"""
PROBE VALPROFILE (Wilson's norm-route program) -- v_p(N) at every p; product formula -> geometric mean (2026-07-26).

N = Prod_{3-nmid a} pi-hat_k(a) = N_{Q(zeta_{3^k})/Q}(pi-hat(1)), rational (result_NORMCHECK).
Program:
 1. v_p(N) at every relevant p. Expect: 3 (ramified, v_3=1-phi confirmed), divisors of 2^{ord}-1 (ord=ord_2(3^k)),
    and small numerator primes. (LTE: v_3(2^ord-1)=k since ord=2*3^{k-1}.)
 2. Product formula: |N|_inf = Prod_p p^{v_p(N)} (auto for rational N) -> gives |N|_inf from the valuations.
 3. Geometric mean of spectrum = |N|_inf^{1/phi} = exp(Sum_a ln|pi-hat(a)| / phi). Compare typical sqrt(k)*3^{-k/2}.
    If arithmetic reproduces it => a spectral moment DERIVED not measured (technique aimable at S_inf).

Reuses probe_normcheck.exact_stationary. Factors 2^ord-1 (k<=5).
"""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(2000000)
import numpy as np
import sympy as sp
from probe_normcheck import exact_stationary, v3
from probe_singlerec import fwd_hat
from probe_6_conservation_generalize import order_of_two


def vp(x, p):
    """p-adic valuation of a sympy Rational."""
    x = sp.Rational(x)
    if x == 0:
        return None
    num, den, v = x.p, x.q, 0
    while num % p == 0:
        num //= p; v += 1
    while den % p == 0:
        den //= p; v -= 1
    return v


def main():
    t0 = time.time()
    print("# PROBE VALPROFILE -- v_p(N) profile + product formula -> geometric mean of the spectrum\n")
    T = sp.symbols('T')
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for k in range(2, KMAX + 1):
        N, U, pi = exact_stationary(k)
        phiN = 2 * 3 ** (k - 1)
        D = sp.ilcm(*[sp.Rational(p).q for p in pi])
        Q = sum(int(D * p) * T ** x for p, x in zip(pi, U))
        Norm = sp.Rational(sp.resultant(sp.cyclotomic_poly(N, T), Q, T), D ** phiN)

        ordv = order_of_two(N)
        Mers = 2 ** ordv - 1
        mfac = sp.factorint(Mers)                      # 2^ord - 1 factorization
        primes = sorted(set([2, 3]) | set(mfac.keys()))
        prof = {p: vp(Norm, p) for p in primes}
        # residual unit after stripping listed primes
        R = sp.Rational(Norm)
        for p in primes:
            R = R / sp.Rational(p) ** prof[p]
        Rfac = sp.factorint(R.p) if abs(R.p) > 1 else {}
        Rden = sp.factorint(R.q) if R.q > 1 else {}

        lnN = math.log(abs(Norm.p)) - math.log(abs(Norm.q))
        gm = lnN / phiN                                 # ln(geometric mean)
        typ = math.log(math.sqrt(k)) - (k / 2) * math.log(3)   # ln(sqrt(k)*3^-k/2)

        print(f"## k={k}  phi={phiN}  ord_2(3^k)={ordv}  2^ord-1 factors: {dict(mfac)}")
        print(f"   v_p(N): " + "  ".join(f"v_{p}={prof[p]}" for p in primes)
              + f"   [v_3 pred {1-phiN}: {'OK' if prof[3]==1-phiN else 'NO'}; "
              f"v_3(2^ord-1)={mfac.get(3,0)} vs k={k}]")
        resid = "±1" if abs(R.p) == 1 and R.q == 1 else f"num{dict(Rfac)}/den{dict(Rden)}"
        print(f"   residual after stripping {{2,3,div(2^ord-1)}}: {resid}")
        print(f"   ln(geo mean)=ln|N|/phi = {gm:.5f}   vs typical ln(sqrt(k)3^-k/2)={typ:.5f}   "
              f"(diff {gm-typ:+.4f})")
        print()

    print("   [residual '±1' => N is supported on {2,3,div(2^ord-1)} = a fully closed arithmetic form for |N|_inf.")
    print("    geo mean vs typical: if they track, the arithmetic REPRODUCES a spectral moment (derived, not measured).]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
