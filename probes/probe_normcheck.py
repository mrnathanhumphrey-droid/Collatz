"""
PROBE NORMCHECK (free check) -- is Prod_{3-nmid a} pi-hat_k(a) rational? value + 3-adic valuation (2026-07-26).

Galois/orbit claim (Wilson): the set {3-nmid a} = units (Z/3^k)* is Galois-stable, so
    Prod_{a in (Z/3^k)*} pi-hat(a) = N_{Q(zeta_{3^k})/Q}(pi-hat(1))  is RATIONAL  (field norm).
|N(pi-hat(1))| = geometric mean of the spectrum (AGGREGATE-side functional, where the channels live).
If rational => the Mahler/house/Dwork rationality shelf opens (Siegel diss p.92-93 = Dwork's Thm) + the Siegel pairing.
If NOT rational => Wilson's Galois observation is wrong.

EXACT: forward-Syracuse stationary pi on units mod 3^k (Fractions, FULL v-period M=ord_2(3^k), untruncated).
  P_poly(T)=Sum_{x in U} pi(x) T^x ; D=lcm denominators ; Q=D*P_poly in Z[T].
  N(P(zeta))=Res_T(Phi_{3^k}, Q)/D^{phi(3^k)}   (Phi monic => Res=Prod_{units}Q(zeta^a)).
  v_3(N) = v_3(Res) - phi(3^k)*v_3(D).
Float cross-check vs fwd_hat (sign convention): Prod|pi-hat(a)| = |N|.
"""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(2000000)
from math import gcd
from fractions import Fraction
import numpy as np
import sympy as sp
from probe_singlerec import fwd_hat
from probe_6_conservation_generalize import order_of_two


def v3(x):
    """3-adic valuation of a nonzero rational (sympy Rational or Fraction)."""
    x = sp.nsimplify(x) if not isinstance(x, sp.Rational) else x
    x = sp.Rational(x)
    if x == 0:
        return None
    num, den = x.p, x.q
    v = 0
    while num % 3 == 0:
        num //= 3; v += 1
    while den % 3 == 0:
        den //= 3; v -= 1
    return v


def exact_stationary(k):
    """Exact forward-Syracuse stationary measure on units mod 3^k (Fractions, full period)."""
    N = 3 ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    U = [a for a in range(1, N) if a % 3 != 0]
    idx = {a: i for i, a in enumerate(U)}
    n = len(U)
    # transition weights (unnormalized): from r, for v=1..M, target=(3r+1)*inv2^v mod N, weight 2^-v
    P = sp.zeros(n, n)
    for r in U:
        p = 1
        row = {}
        for v in range(1, M + 1):
            p = (p * inv2) % N
            t = ((3 * r + 1) * p) % N
            row[t] = row.get(t, 0) + sp.Rational(1, 2 ** v)
        s = sum(row.values())
        for t, w in row.items():
            P[idx[r], idx[t]] = w / s           # normalized row (divide by 1-2^-M)
    # stationary: left null vector of (P - I): pi (P - I) = 0  => (P^T - I) pi^T = 0
    A = P.T - sp.eye(n)
    ns = A.nullspace()
    assert len(ns) == 1, f"k={k}: nullspace dim {len(ns)}"
    pi = ns[0]
    tot = sum(pi)
    pi = [pi[i] / tot for i in range(n)]
    return N, U, pi


def main():
    t0 = time.time()
    print("# PROBE NORMCHECK -- Prod_{3-nmid a} pi-hat_k(a) = N(pi-hat(1)) rational? value + v_3\n")
    T = sp.symbols('T')
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print(f"   {'k':>2} {'phi(3^k)':>9} {'v_3(Norm)':>10} {'pred':>10} {'Norm size':>16} {'ln|N| vs Sum ln|pi|':>24}")
    for k in range(2, KMAX + 1):
        N, U, pi = exact_stationary(k)
        phiN = 2 * 3 ** (k - 1)
        D = sp.ilcm(*[sp.Rational(p).q for p in pi])
        Q = sum(int(D * p) * T ** x for p, x in zip(pi, U))
        PhiN = sp.cyclotomic_poly(N, T)
        Res = sp.resultant(PhiN, Q, T)          # = Prod_{units} Q(zeta^a) (Phi monic)
        Norm = sp.Rational(Res, D ** phiN)
        vN = v3(Norm)
        # float cross-check (math.log on big ints, no huge-str)
        ph, _ = fwd_hat(k)
        prodabs = float(np.sum([np.log(abs(ph[a])) for a in U]))   # log|prod|
        logNorm = math.log(abs(Norm.p)) - math.log(abs(Norm.q))
        chk = f"{logNorm:.3f} vs {prodabs:.3f}"
        pred = -(phiN - 1)
        ndig = f"num{len(str(abs(Norm.p)))}d/den{len(str(abs(Norm.q)))}d"
        print(f"   {k:>2} {phiN:>9} {str(vN):>10} {str(pred):>10} {ndig:>16} {chk:>24} "
              f"[{'OK' if abs(logNorm-prodabs)<1e-3 else 'MISMATCH'}]")
    print("\n   [v_3 pred = -(phi(3^k)-1) = 1-2*3^{k-1}. Norm rational => Galois/orbit claim holds;")
    print("    float check: ln|Norm| must equal Sum_a ln|pi-hat(a)| (sign-convention-free).]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
