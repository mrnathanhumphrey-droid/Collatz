"""
PROBE R9 -- THE COLLISION IDENTITY GATE + gamma-TABLES. Exact rationals.
Gates Wilson's session:  S_K = 2 * sum_{m>=1} 4^{-m} gamma_{K-1}(tau_m),  exact at every finite K.

  tau_m := (4^{-m} - 1)/3  in Z_3   (v_3(tau_m) = v_3(m), LTE); computed mod 3^n from 4^{-m} mod 3^{n+1}.
  gamma_n(tau) := 3^n * Pr_{a,a' ~ mu_n iid}[ (a-a') + tau(1+3a) == 0 mod 3^n ]
               =  3^n * sum_a mu_n(a) * mu_n( (a + tau(1+3a)) mod 3^n )   [unique a' per a; O(supp)].
  gamma_n(0) = 3^n ||mu_n||^2 = X_n  (the qx+1 corpus accumulation).

gamma computed DIRECTLY from the mu tables (R7 build), NOT routed through the C-tables (that tests nothing).
Reuses R7 build_mu / mu1 / cram.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7

S = {1: F(2, 3), 2: F(10, 21), 3: F(31370, 67963),
     4: F(143195649659456490, 308468774477179141),
     5: F(2490699741144069281815149277465294323722281911695597204406570,
          5350418720142111510029542161258891403960563152740894082816203),
     6: F(3073594394508858158063860398646918469118621544451223655237787736044140305600491108232683716375149129809572776883214501465626550917830849901208399259893949668774475847028354296093382145075270,
          6593308371642093718618660948586623421072140639443126614173839021148571641533385540268510160579160359417192756329060861771589941737564918045220171980607356372262183045471725888175409486651187)}


def v3(n):
    n = abs(int(n))
    if n == 0:
        return 10 ** 9
    j = 0
    while n % 3 == 0:
        n //= 3; j += 1
    return j


def tau(m, n):
    """tau_m mod 3^n, from 4^{-m} mod 3^{n+1}."""
    Mp = 3 ** (n + 1)
    r = pow(pow(4, -1, Mp), m, Mp)         # 4^{-m} mod 3^{n+1}
    return ((r - 1) // 3) % (3 ** n)        # (4^{-m}-1)/3 mod 3^n


def gamma(mu_n, n, t):
    """gamma_n(tau=t) exact, via unique-collision partner a' = a + t(1+3a) mod 3^n."""
    M = 3 ** n
    s = F(0)
    for a, pa in mu_n.items():
        ap = (a + t * (1 + 3 * a)) % M
        pb = mu_n.get(ap)
        if pb:
            s += pa * pb
    return (3 ** n) * s


def S_from_gamma(K, mu):
    """S_K = 2 sum_{r=1}^{P} [4^{-r}/(1-4^{-P})] gamma_{K-1}(tau_r),  P=3^{K-1}."""
    n = K - 1
    P = 3 ** n
    geom = 1 - F(1, 4 ** P)
    tot = F(0)
    for r in range(1, P + 1):
        t = tau(r, n)
        w = F(1, 4 ** r) / geom
        tot += w * gamma(mu[n], n, t)
    return 2 * tot


def prim_plancherel(mu_n, n):
    """R66 object: sum_{xi primitive mod 3^n} |mu_hat_n(xi)|^2 = sum_{a,a'} mu mu' c_{3^n}(a-a')  (exact)."""
    s = F(0)
    for a, pa in mu_n.items():
        for ap, pap in mu_n.items():
            s += pa * pap * R7.cram(n, a - ap)
    return s


def main():
    print("# PROBE R9 -- THE COLLISION IDENTITY GATE + gamma-TABLES. Exact rationals.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 6):
        mu[k] = R7.build_mu(mu[k - 1], k)
    for n in range(1, 6):
        print(f"   |supp mu_{n}| = {len(mu[n])} (mod 3^{n})")
    print()

    # ---- R9-A ----
    print("## R9-A  COLLISION IDENTITY GATE  (S_K = 2 sum 4^-m gamma_{K-1}(tau_m); gamma from mu tables)")
    okA = True
    for K in range(2, 7):
        val = S_from_gamma(K, mu)
        good = (val == S[K])
        okA = okA and good
        print(f"   K={K}: engine S_K = {str(val)[:34]}{'...' if len(str(val))>34 else ''}  "
              f"vs frozen {str(S[K])[:20]}{'...' if len(str(S[K]))>20 else ''}  [{'PASS' if good else 'FAIL'}]")
        if not good:
            print(f"      FAIL K={K}: engine={val}\n              frozen={S[K]}")
    # hand anchor
    n = 1
    g = {r: gamma(mu[1], 1, tau(r, 1)) for r in (1, 2, 3)}
    print(f"   hand anchor K=2: gamma_1(tau_1..3)={[str(g[r]) for r in (1,2,3)]} (=2/3,2/3,5/3); "
          f"2*(16*{g[1]}+4*{g[2]}+1*{g[3]})/63... = {S_from_gamma(2,mu)} (=10/21)")
    print(f"   => R9-A {'GATE PASS' if okA else 'FAIL -- gamma route dies at normalization (#31)'}\n")

    # ---- R9-B ----
    print("## R9-B  gamma-TABLES + X_n WELD  (n=1..5; gamma_n(0)=X_n vs corpus X=1+sum_{j<=n}S_j)")
    Xrun = F(1)
    for n in range(1, 6):
        Xrun_prev = Xrun
        Xn = F(1) + sum(S[j] for j in range(1, n + 1))
        g0 = gamma(mu[n], n, 0)
        print(f"   n={n}: gamma_n(0)={g0} = X_n?  corpus X_n={Xn}  [{'WELD' if g0==Xn else 'DEV'}]  (float {float(g0):.6f})")
    print("   gamma_n(tau_m) by stratum j=v3(m) (distinct values per stratum; DC=highest j is tau=0 line):")
    for n in range(1, 6):
        P = 3 ** n
        strat = {}
        for r in range(1, P + 1):
            j = v3(r) if r < P else n
            strat.setdefault(j, set()).add(gamma(mu[n], n, tau(r, n)))
        parts = []
        for j in sorted(strat):
            vals = sorted(strat[j])
            tag = "DC" if j == n else f"j={j}"
            vs = ",".join(str(v) for v in vals) if len(vals) <= 4 else f"{len(vals)} vals ~{float(min(vals)):.3f}..{float(max(vals)):.3f}"
            parts.append(f"[{tag}: {vs}]")
        print(f"   n={n}: " + " ".join(parts))
    print()

    # ---- R9-C ----
    print("## R9-C  OFF-DC CONVERGENCE  (fixed m, gamma_n(tau_m) n=1..5; measurement, NO fit)")
    print(f"   {'m':>3} {'v3(m)':>5} | " + " ".join(f"{'n='+str(n):>22}" for n in range(1, 6)))
    for m in (1, 2, 3, 4, 9):
        cells = []
        for n in range(1, 6):
            gv = gamma(mu[n], n, tau(m, n))
            cells.append(f"{str(gv)[:14]:>14}({float(gv):.3f})")
        print(f"   {m:>3} {v3(m):>5} | " + " ".join(f"{c:>22}" for c in cells))
    print("   (Wilson's pre-reg: each fixed-m column CONVERGES; divergence only at tau=0. Verbatim above.)\n")

    # ---- R9-D ----
    print("## R9-D  NUMERAL WELDS  (measurement, table only, NO conclusion)")
    print("   (i) R66 primitive object P_n=sum_prim|mu_hat|^2  vs  <gamma_n>=sum 4^-m gamma_n(tau_m)/(1/3)")
    print(f"   {'n':>2} {'P_n=sum_prim|muhat|^2':>26} {'S_n(frozen)':>14} {'<gamma_n>=3*sum4^-m gamma':>28} {'3*S_{n+1}/2':>14}")
    for n in range(2, 6):
        Pn = prim_plancherel(mu[n], n)
        geom = 1 - F(1, 4 ** (3 ** n))
        gbar = 3 * sum(F(1, 4 ** r) / geom * gamma(mu[n], n, tau(r, n)) for r in range(1, 3 ** n + 1))
        print(f"   {n:>2} {str(Pn)[:26]:>26} {str(S[n])[:14]:>14} {str(gbar)[:28]:>28} {str(3*S[n+1]/2)[:14]:>14}  "
              f"[P_n==S_n:{Pn==S[n]}; <g>==3S_{{n+1}}/2:{gbar==3*S[n+1]/2}]")
    print("   (ii) #27 chain amplitude: NOT RECOVERABLE from archive (no frozen definition found; grep empty).")
    print("        Reported as NOT RECOVERABLE per spec -- not reconstructed.")


if __name__ == "__main__":
    main()
