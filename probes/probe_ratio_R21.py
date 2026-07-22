"""
PROBE R21 -- THE RATIO LAW. Reuses R7/R9/R10. The plainest reading of gamma (no new coordinate; R12-B settled that).

Defining congruence: (1+3a)(1+3tau) == (1+3a') mod 3^{r+1}, u=1+3tau_m=4^{-m}  =>  X'/X == 4^{-m}.
  rho_r(u) := Pr_{a,a'~mu_r}[ X'/X == u mod 3^{r+1} ],  X=1+3a, X'=1+3a'.   gamma_r(tau_m) = 3^r * rho_r(4^{-m}).
  f(u) := 3^r*rho_r(u) -> Haar density on 1+3Z_3 of the ratio of two iid Syracuse values.
  THEOREM:  Sum_{m>=1} 4^{-m} f(4^{-m}) = 7/30   (<=>  Sum 4^{-m}[f-1] = -1/10;  1/3-7/30=1/10=-Sum Lambda_r).

R21-A group-division gate (independent route).  R21-B full density.  R21-C argmax r=7 (derived prediction 2^7).
R21-D orbit vs bulk (geometric weight on atypical or representative part).
"""
import os, sys, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from math import lcm
import probe_engine_R7 as R7
import probe_gamma_R9 as R9
import probe_charledger_R10 as R10

v3 = R9.v3
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
S = {k: F(7, 15) + EPS[k] for k in EPS}


def ratio_law(mu_r, r):
    """rho_r(u)=Pr[X'/X==u mod 3^{r+1}] by DIRECT group division (integer weights). Independent of tau/C/engine."""
    Mp = 3 ** (r + 1)
    D = 1
    for pa in mu_r.values():
        D = lcm(D, pa.denominator)
    n = [((1 + 3 * a) % Mp, pa.numerator * (D // pa.denominator)) for a, pa in mu_r.items()]
    rho = {}
    for X, nX in n:
        iX = pow(X, -1, Mp)
        for Xp, nXp in n:
            R = (Xp * iX) % Mp
            rho[R] = rho.get(R, 0) + nX * nXp
    return {u: F(c, D * D) for u, c in rho.items()}


def mu_hat(mu_r, r, xi):
    N = 3 ** r
    return sum(complex(p) * cmath.exp(2j * math.pi * (xi * a % N) / N) for a, p in mu_r.items())


def main():
    print("# PROBE R21 -- THE RATIO LAW.  gamma_r(tau_m)=3^r rho_r(4^-m); f=Haar ratio density.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    RHO = {}
    for r in range(2, 7):
        RHO[r] = ratio_law(mu[r], r)

    # ================= R21-A =================
    print("## R21-A  DIRECT RATIO LAW (forced, a real gate): group-division rho_r vs banked collision-gamma, all m")
    okA = True
    for r in range(2, 7):
        Mp = 3 ** (r + 1); i4 = pow(4, -1, Mp); bad = 0
        for m in range(1, 3 ** r + 1):
            u = pow(i4, m, Mp)
            if 3 ** r * RHO[r].get(u, F(0)) != R9.gamma(mu[r], r, R9.tau(m, r)):
                bad += 1
        good = (bad == 0) and (sum(RHO[r].values()) == 1)
        okA = okA and good
        print(f"   r={r}: gamma_r(tau_m)==3^r rho_r(4^-m) all {3**r} m? {bad==0}   sum rho=1? {sum(RHO[r].values())==1}  [{'PASS' if good else 'FAIL'}]")
    print(f"   => R21-A {'GATE PASS -- ratio-density route independent of tau/C-tables/engine, reproduces gamma' if okA else 'FAIL (#41)'}\n")

    # weld: Sum_m 4^-m f_r(4^-m) == S_{r+1}/2 (-> 7/30)
    print("## R21 weld  Sum_{m=1..3^r} 4^-m f_r(4^-m) == S_{r+1}/2  (ratio-density reproduces 7/30)")
    for r in range(2, 7):
        Mp = 3 ** (r + 1); i4 = pow(4, -1, Mp); P = 3 ** r; geom = 1 - F(1, 4 ** P)
        wsum = sum(F(1, 4 ** m) / geom * 3 ** r * RHO[r].get(pow(i4, m, Mp), F(0)) for m in range(1, P + 1))
        tgt = S[r + 1] / 2 if r + 1 in S else None
        print(f"   r={r}: Sum 4^-m f_r = {float(wsum):.8f}  vs S_{{r+1}}/2={float(tgt):.8f}  [{'OK' if wsum==tgt else 'DEV'}]  (7/30={7/30:.8f})")
    print()

    # ================= R21-B =================
    print("## R21-B  THE FULL DENSITY (measurement, NO fit): f_r(u)=3^r rho_r(u), all u in 1+3Z/3^{r+1}")
    for r in (4, 5, 6):
        Mp = 3 ** (r + 1)
        dlog = {}                                   # 4^t -> t
        el = 1
        for t in range(3 ** r):
            dlog[el] = t; el = (el * 4) % Mp
        f = {u: float(3 ** r * RHO[r].get(u, F(0))) for u in dlog}   # all group elements
        vals = list(f.values())
        mean = sum(vals) / len(vals)
        var = sum((x - 1) ** 2 for x in vals) / len(vals)            # ||f-1||^2 / |G|
        # concentration of |f-1| by v3(u-1)
        strat = {}
        for u, fu in f.items():
            j = v3(u - 1) if u != 1 else r + 99
            strat.setdefault(min(j, r), []).append(abs(fu - 1))
        conc = "  ".join(f"v3(u-1)={j}:{sum(strat[j])/len(strat[j]):.3f}" for j in sorted(strat) if j <= r)
        srt = sorted(f.items(), key=lambda kv: kv[1])
        print(f"   r={r} (|G|=3^r={3**r}): min={min(vals):.4f} max={max(vals):.4f} mean={mean:.4f} "
              f"||f-1||^2/|G|={var:.4f}")
        print(f"        mean|f-1| by stratum: {conc}")
        top = "  ".join(f"m{(-dlog[u])%(3**r)}:{fu:.3f}" for u, fu in srt[-5:])
        bot = "  ".join(f"m{(-dlog[u])%(3**r)}:{fu:.3f}" for u, fu in srt[:5])
        print(f"        top5 f (u=4^-m): {top}")
        print(f"        bot5 f (u=4^-m): {bot}")
    print("   [Q: bounded? max grow w/ r? where does |f-1| concentrate? (orbit=whole group, so this is the full density)]\n")

    # ================= R21-C =================
    print("## R21-C  ARGMAX AT r=7 (tests derived prediction argmax=2^7=128 or conj 2059)")
    N = 3 ** 7
    xs = {xi: abs(mu_hat(mu[7], 7, xi)) for xi in range(1, N) if xi % 3 != 0}
    xstar = max(xs, key=xs.get)
    pred, conj = 128, N - 128
    d0 = min(xstar / N, 1 - xstar / N)
    print(f"   argmax xi*={xstar}  |mu_hat|={xs[xstar]:.6f}   x=xi*/N={xstar/N:.5f}  dist-to-0={d0:.5f}  (2/3)^7={(2/3)**7:.5f}")
    print(f"   prediction 2^7=128: |mu_hat(128)|={xs.get(128,float('nan')):.6f}  (conj 2059: {xs.get(2059,float('nan')):.6f}); "
          f"argmax==128 or 2059? {xstar in (128, 2059)}")
    top5 = sorted(xs.items(), key=lambda kv: -kv[1])[:5]
    print(f"   top5: " + "  ".join(f"{xi}(2^{round(math.log2(xi)) if xi&(xi-1)==0 else '?'}):{val:.5f}" for xi, val in top5))
    print("   [pre-registered miss is informative: near-trivial region is near-degenerate under superpoly decay]\n")

    # ================= R21-D =================
    print("## R21-D  ORBIT vs BULK (measurement, NO fit): geometric-weighted |f-1| vs unweighted mean |f-1|")
    for r in (4, 5, 6):
        Mp = 3 ** (r + 1); i4 = pow(4, -1, Mp); P = 3 ** r; geom = 1 - F(1, 4 ** P)
        f = {u: float(3 ** r * RHO[r].get(u, F(0))) for u in RHO[r]}
        # full group for unweighted:
        allu = {}
        el = 1
        for t in range(P):
            allu[el] = float(3 ** r * RHO[r].get(el, F(0))); el = (el * 4) % Mp
        unw = sum(abs(x - 1) for x in allu.values()) / P
        wsum = sum(float(F(1, 4 ** m) / geom) * abs(3 ** r * float(RHO[r].get(pow(i4, m, Mp), F(0))) - 1) for m in range(1, P + 1))
        print(f"   r={r}: geometric-weighted mean|f-1|={wsum:.4f}   unweighted mean|f-1|={unw:.4f}   "
              f"ratio={wsum/unw:.3f}  [{'representative' if 0.5<wsum/unw<2 else 'ATYPICAL'}]")
    print("   [Q: does the geometric weight 4^-m sit on an atypical part of the density, or a representative one?]")


if __name__ == "__main__":
    main()
