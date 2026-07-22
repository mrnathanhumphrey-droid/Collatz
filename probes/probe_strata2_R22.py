"""
PROBE R22 -- IS f STRATUM-ONLY. Reuses R7/R9/R10/R21. Tests THE assumption under the stratum-reduction framework:
does f_r(u)=3^r rho_r(u) depend only on j=v3(u-1)-1, or is it genuinely u-dependent within a stratum?

If stratum-only, the theorem Sum_m 4^-m f(4^-m)=7/30 collapses to Sum_j W_j F(j)=7/30 (R8 strata weights), with a
second (Haar) constraint Sum_j (2/3)3^-j F(j)=1. R22-A is the gate; B/C/D/E are the reduction + diagnostics.
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_gamma_R9 as R9
import probe_ratio_R21 as R21

v3 = R9.v3


def Wj(j):
    """W_j = Sum_{m: v3(m)=j} 4^-m = x/(1-x) - x^3/(1-x^3), x=4^{-3^j}  (R8, exact)."""
    x = F(1, 4 ** (3 ** j))
    return x / (1 - x) - x ** 3 / (1 - x ** 3)


def stats(vals):
    n = len(vals); mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / n
    return min(vals), max(vals), mean, math.sqrt(max(0.0, var))


def main():
    print("# PROBE R22 -- IS f STRATUM-ONLY.  f_r(u)=3^r rho_r(u); stratum j=v3(u-1)-1.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)
    RHO = {r: R21.ratio_law(mu[r], r) for r in (4, 5, 6)}

    # densities grouped by stratum (exclude u=1 DC)
    Ftab = {}                                   # Ftab[r][j] = list of f values in stratum j
    for r in (4, 5, 6):
        Mp = 3 ** (r + 1); grp = {}
        el = 1
        for t in range(3 ** r):                  # walk the whole group via <4>
            u = el; el = (el * 4) % Mp
            if u == 1:
                continue
            j = v3(u - 1) - 1
            grp.setdefault(j, []).append(float(3 ** r * RHO[r].get(u, F(0))))
        Ftab[r] = grp

    # ================= R22-A =================
    print("## R22-A  STRATUM-ONLY GATE (measurement; the assumption under everything): within-stratum spread of f")
    for r in (4, 5, 6):
        print(f"   r={r}:")
        print(f"     {'j':>2} {'#':>5} {'min':>8} {'max':>8} {'mean':>8} {'std':>8} {'std/mean':>9} {'max/mean':>9}")
        for j in sorted(Ftab[r]):
            vals = Ftab[r][j]
            mn, mx, mean, sd = stats(vals)
            print(f"     {j:>2} {len(vals):>5} {mn:>8.4f} {mx:>8.4f} {mean:>8.4f} {sd:>8.4f} {sd/mean:>9.4f} {mx/mean:>9.4f}")
    print("   [Q: does within-stratum std/mean -> 0 with r (stratum-only in limit), stay O(1) (u-dependent), or ~const?]\n")

    # ================= R22-B =================
    print("## R22-B  THE F(j) TABLE (measurement, NO fit): stratum means + successive diffs vs 7/15 and 0.6632")
    for r in (4, 5, 6):
        Fj = {j: sum(Ftab[r][j]) / len(Ftab[r][j]) for j in sorted(Ftab[r])}
        js = sorted(Fj)
        row = "  ".join(f"F({j})={Fj[j]:.4f}" for j in js)
        diffs = "  ".join(f"dF({j})={Fj[j+1]-Fj[j]:+.4f}" for j in js if j + 1 in Fj)
        print(f"   r={r}: {row}")
        print(f"        diffs: {diffs}   [vs 7/15=0.4667, 0.6632]")
    print("   [Q: diff sequence flat, or decreasing toward 7/15 from above (convexity)?]\n")

    # ================= R22-C =================
    print("## R22-C  HAAR CHECK (forced): Sum_j (2/3)3^-j F_r(j) + f(1)/3^r == 1 (normalization audit; fail=binning bug)")
    for r in (4, 5, 6):
        Fj = {j: sum(Ftab[r][j]) / len(Ftab[r][j]) for j in sorted(Ftab[r])}
        haarsum = sum((2 / 3) * 3 ** (-j) * Fj[j] for j in Fj)
        fDC = float(3 ** r * RHO[r].get(1, F(0)))
        total = haarsum + fDC / 3 ** r
        print(f"   r={r}: Sum_j (2/3)3^-j F(j) = {haarsum:.6f}  + f(1)/3^r={fDC/3**r:.6f}  = {total:.6f}  "
              f"[{'OK' if abs(total-1)<1e-9 else 'BINNING BUG'}]  (f(1)=X_r={fDC:.4f})")
    print()

    # ================= R22-D =================
    print("## R22-D  F(0) CONVERGENCE (measurement, NO fit): gamma_r(tau_1,2,4) [all j=0] r=1..7, exact+float, spread")
    print(f"   {'r':>2} {'gamma_r(tau_1)':>16} {'gamma_r(tau_2)':>16} {'gamma_r(tau_4)':>16} {'spread':>9} {'vs 0.66841':>11}")
    for r in range(1, 8):
        g = [R9.gamma(mu[r], r, R9.tau(m, r)) for m in (1, 2, 4)]
        gf = [float(x) for x in g]
        spread = max(gf) - min(gf)
        print(f"   {r:>2} {gf[0]:>16.6f} {gf[1]:>16.6f} {gf[2]:>16.6f} {spread:>9.2e} {sum(gf)/3-0.66841:>+11.5f}")
    print("   [Q: do the j=0 orbit points agree at each r (stratum-only ON the orbit)? where is the common value heading?]\n")

    # ================= R22-E =================
    print("## R22-E  WEIGHT/HAAR MISMATCH (measurement, cheap): W_j vs (2/3)3^-j, ratio (leverage of the 2nd equation)")
    print(f"   {'j':>2} {'W_j (geometric)':>18} {'(2/3)3^-j (Haar)':>18} {'ratio W_j/Haar':>15}")
    for j in range(0, 5):
        w = Wj(j); haar = F(2, 3) * F(1, 3 ** j)
        print(f"   {j:>2} {float(w):>18.8f} {float(haar):>18.8f} {float(w/haar):>15.6f}")
    print(f"   Sum W_j = {float(sum(Wj(j) for j in range(0,40))):.6f} (=1/3=Lambda_0)   "
          f"Sum (2/3)3^-j = {float(sum(F(2,3)*F(1,3**j) for j in range(0,40))):.6f} (=1)")
    print("   [the two constraints (7/30 geometric, 1 Haar) differ only through this ratio.]")


if __name__ == "__main__":
    main()
