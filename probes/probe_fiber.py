"""
PROBE FIBER -- the fiber-mean / fiber-fluctuation split (replaces the retracted interference ledger).

Level-lift is 3-to-1: primitive chi_k at level r -> psi_m at level r-1, m = k mod 3^{r-1},
fiber(m) = {m, m+3^{r-1}, m+2*3^{r-1}}. Split delta_r:
  fiber-mean       dbar_r(k) = mean of delta_r over its fiber   (constant on fibers; level r-1 content)
  fiber-fluctuation dpr_r(k) = delta_r(k) - dbar_r(k)            (mean-zero within each fiber; the new digit)

Fiber-average of Re w: M(Re w) = Re w with 4 -> 4^3 = 64 (kills modes not div by 3, reindexes). Doubly-exp.
Claim: Lambda_r^unif = S_r * (primitive mean of Re w) = -1/(2(4^{3^{r-1}}-1))-ish, closed form, doubly-exp.
And sign(g_r) = sign(<fiber-fluctuation, Re w>) with the fiber-mean channel doubly-exp dead (no propagating tail).

GATES (in order, all falsifying):
 1. M(Re w) == Re w with 4->64 to machine precision (lift-indexing plumbing).
 2. primitive mean of Re w == Lambda_r^unif/S_r closed form == -1/7, -0.00793078, -8.8e-7 (the split is RIGHT).
 3. <fiber-fluctuation, Re w> for r=2..12: sign and rate (should be the ladder 0.89-0.91, the only channel now).
"""
import os, sys, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

Rew = lambda x: 15.0 / (2 * (17 - 8 * np.cos(2 * np.pi * x))) - 0.5
Rew_q = lambda y, q: (q * np.cos(2 * np.pi * y) - 1) / (q * q - 2 * q * np.cos(2 * np.pi * y) + 1)  # Re[1/(q e(y)-1)]
JMAX = 12

_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): float(F(int(v['num']), int(v['den']))) for k, v in _hist.items()}
EPS.update({9: -7.520257156400000e-6, 10: 7.207509171100000e-4, 11: 1.501967012082273e-3, 12: 2.274713720558208e-3})
LAM_NU = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}
for r in range(12, 17):
    EPS[r + 1] = EPS[r] + 2 * LAM_NU[r]
LAM = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 17)}
Sval = {r: 7.0 / 15 + EPS[r] for r in range(1, 17)}


def delta_from_nu(nuj, j):
    N = 3 ** j
    mu = np.zeros(N)
    for X, w in nuj.items():
        mu[(X - 1) // 3 % N] += float(w)
    d = R10.dlog_table(j)
    rho = np.zeros(N)
    np.add.at(rho, np.array([d[a] for a in range(N)]), mu)
    prof = np.abs(np.fft.fft(rho)) ** 2
    prim = np.array([k for k in range(1, N) if k % 3 != 0]); M = len(prim)
    S = float(prof[prim].sum())
    dd = np.zeros(N); dd[prim] = prof[prim] / S - 1.0 / M
    return dd, S


def lam_unif_closed(r):
    """exact primitive-mean of Re w over order-3^r characters = (1/(2*3^{r-1}))[3^r/(4^{3^r}-1)-3^{r-1}/(4^{3^{r-1}}-1)]."""
    Nm = 3 ** (r - 1); N = 3 ** r
    from fractions import Fraction as FF
    return float(FF(1, 2 * Nm) * (FF(N, 4 ** N - 1) - FF(Nm, 4 ** Nm - 1)))


def main():
    t0 = time.time()
    print(f"# PROBE FIBER -- fiber-mean/fluctuation split, r=2..{JMAX}.\n")
    nus = build_nu(0.5, JMAX)
    dlt = {};
    for j in range(2, JMAX + 1):
        dlt[j], _ = delta_from_nu(nus[j], j)
    print(f"  delta_j built ({time.time()-t0:.1f}s)\n")

    # ---- GATE 1: M(Re w) == Re w with 4->64 ----
    print("## GATE 1  fiber-average M(Re w) == Re w with 4->64 (lift-indexing plumbing)")
    for r in range(2, 8):
        N = 3 ** r; Nm = 3 ** (r - 1)
        rew = Rew(np.arange(N) / N)
        Mrew = (rew[0:Nm] + rew[Nm:2 * Nm] + rew[2 * Nm:3 * Nm]) / 3       # fiber-average, fn of m
        rew64 = Rew_q(np.arange(Nm) / Nm, 64.0)
        err = float(np.max(np.abs(Mrew - rew64)))
        print(f"   r={r}: max|M(Re w) - Re w(4->64)| = {err:.2e}  [{'OK' if err < 1e-12 else 'FAIL'}]")
    print()

    # ---- GATE 2: primitive mean == Lambda^unif/S closed form ----
    print("## GATE 2  primitive-mean of Re w == Lambda_r^unif/S_r (closed form, doubly-exp)")
    print(f"   {'r':>2} {'prim-mean(numeric)':>18} {'closed form':>14} {'Wilson':>12}")
    wil = {1: -1.0 / 7, 2: -0.00793078, 3: -8.8e-7}
    for r in range(1, 7):
        N = 3 ** r
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim)
        pm = float(np.mean([Rew(k / N) for k in prim]))
        cf = lam_unif_closed(r)
        w = wil.get(r, float('nan'))
        print(f"   {r:>2} {pm:>18.9e} {cf:>14.6e} {w if not math.isnan(w) else '':>12}")
    print("   [prim-mean = closed form = -1/(2(4^{3^{r-1}}-1))-ish; explains banked doubly-exp death of Lambda^unif.]\n")

    # ---- fiber decomposition + GATE 3 ----
    print("## GATE 3  fiber split of delta_r: <fiber-mean,Rew>, <fluctuation,Rew>, g_r; sign & rate of fluctuation")
    print(f"   {'r':>2} {'<dbar,Rew>':>13} {'<dfluct,Rew>':>14} {'g_r':>13} {'Lam^unif/S':>13} {'sgn fluct':>10}")
    gfl = {}
    for r in range(2, JMAX + 1):
        N = 3 ** r; Nm = 3 ** (r - 1)
        d = dlt[r]; rew = Rew(np.arange(N) / N)
        dbar_m = (d[0:Nm] + d[Nm:2 * Nm] + d[2 * Nm:3 * Nm]) / 3          # fiber-mean, fn of m
        dbar = np.concatenate([dbar_m, dbar_m, dbar_m])                    # constant on fibers, level r
        dfluct = d - dbar
        g_mean = float(np.sum(dbar * rew))
        g_fluct = float(np.sum(dfluct * rew))
        g_r = float(np.sum(d * rew))
        gfl[r] = g_fluct
        print(f"   {r:>2} {g_mean:>+13.4e} {g_fluct:>+14.5e} {g_r:>+13.5e} {lam_unif_closed(r):>+13.3e} "
              f"{'+' if g_fluct > 0 else '-':>10}")
    print("   rate <fluct,Rew>_r / <fluct,Rew>_{r-1}:")
    rr = [gfl[r] / gfl[r - 1] for r in range(4, JMAX + 1) if abs(gfl[r - 1]) > 1e-15]
    print("   " + " ".join(f"{x:+.3f}" for x in rr))
    print(f"   [fluctuation is the ONLY channel (fiber-mean doubly-exp dead) => sign(g_r)=sign<fluct,Rew>; rate must be ~0.89-0.91.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
