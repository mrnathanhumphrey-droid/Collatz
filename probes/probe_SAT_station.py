"""
PROBE SAT -- THE STATION AT INFINITY. Use the global constraint Sum_{r>=1} Lam_r = -1/10 as a long-baseline datum
on the unreachable tail, and TEST candidate periods against it (period is input, not output => immune to D-B's
under-determination).

KEY EXACT FACT (must be stated): 2 Lam_r = eps_{r+1}-eps_r telescopes, so
    T_M := Sum_{r>=M} Lam_r = (eps_inf - eps_M)/2 = -eps_M/2   (given eps_inf=0, i.e. S->7/15).
So the 'station' is NOT independent info beyond the boundary datum eps_M; the -1/10 constraint is tautologically
-eps_1/2 (eps_1=1/5) + convergence. BUT the tail AMPLITUDE relative to a single term is a genuine PERIOD estimator
for a slowly-damped oscillation: |T|/|term| ~ 1/|1 - rho e^{i theta}|, large for small theta (long period).
"""
import os, sys, json, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EPSJSON = os.path.join(HERE, '..', 'experiments_output', 'result_77_7_eps_exact_through_k8_v2_vec_pool.json')
EPS_F = {1: 0.2, 2: 9.523809523809525e-3, 3: -5.091986325893010e-3, 4: -2.452258248318762e-3,
         5: -1.151746915130986e-3, 6: -4.979056652200001e-4, 7: -1.175236830400000e-3,
         8: -7.455463672900000e-4, 9: -7.520257156400000e-6, 10: 7.207509171100000e-4,
         11: 1.501967012082273e-3, 12: 2.274713720558208e-3}


def fit_ab(Lam, rs, rho, theta):
    """fit Lam_r = rho^r (a cos(theta r) - b sin(theta r)) to rs; return a,b,residual."""
    Aeq = np.array([[rho ** r * math.cos(theta * r), -rho ** r * math.sin(theta * r)] for r in rs])
    y = np.array([Lam[r] for r in rs])
    ab, res, *_ = np.linalg.lstsq(Aeq, y, rcond=None)
    pred = Aeq @ ab
    return ab[0], ab[1], float(np.sqrt(np.mean((pred - y) ** 2)))


def tail(a, b, rho, theta, M=12, N=4000):
    return sum(rho ** r * (a * math.cos(theta * r) - b * math.sin(theta * r)) for r in range(M, N))


def main():
    print("# PROBE SAT -- THE STATION AT INFINITY.\n")
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in json.load(open(EPSJSON)).items()}   # exact k<=8
    Lam = {r: (float((EPS[r + 1] - EPS[r]) / 2) if r + 1 <= 8 else (EPS_F[r + 1] - EPS_F[r]) / 2) for r in range(1, 12)}

    # ============ SAT-A ============
    print("## SAT-A  THE STATION (exact via telescoping: T_M = Sum_{r>=M} Lam_r = -eps_M/2)")
    print(f"   {'M':>2} {'T_M = -eps_M/2':>16} {'|T_M|':>11} {'|T_M|/|Lam_{M-1}|':>17} {'sign':>5} {'eps_M exact?':>12}")
    for M in (8, 9, 10, 11, 12):
        eM = float(EPS[M]) if M <= 8 else EPS_F[M]
        T = -eM / 2
        ratio = abs(T) / abs(Lam[M - 1]) if (M - 1) in Lam else float('nan')
        print(f"   {M:>2} {T:>16.4e} {abs(T):>11.4e} {ratio:>17.4f} {'+' if T > 0 else '-':>5} {'yes' if M <= 8 else 'float':>12}")
    T12 = -EPS_F[12] / 2
    print(f"   => STATION T := T_12 = Sum_{{r>=12}} Lam_r = -eps_12/2 = {T12:.4e}  (|T|/|Lam_11| = {abs(T12)/abs(Lam[11]):.3f})")
    print("   HONEST: T_M = -eps_M/2 is the TELESCOPING identity -- the 'station' is the boundary datum eps_M, NOT")
    print("   independent info beyond k=12. The -1/10 constraint is tautologically -eps_1/2 (eps_1=1/5)+convergence.")
    print("   What IS informative: |T|/|term| ~ 3 -- a large tail from a rho~0.984 oscillation needs a LONG period.\n")

    # ============ SAT-B ============
    print("## SAT-B  THE PERIOD SCAN (fit a,b to r=8..11 with rho=0.984 fixed; predict T_12; compare to exact T)")
    print(f"   target exact tail T = {T12:.4e}.  A period P is CONSISTENT if T_pred(P) matches T.")
    rho = 0.984
    rows = []
    for i in range(int((30 - 6) / 0.25) + 1):
        P = 6 + 0.25 * i
        th = 2 * math.pi / P
        a, b, resid = fit_ab(Lam, [8, 9, 10, 11], rho, th)
        Tp = tail(a, b, rho, th, 12)
        rel = abs(Tp - T12) / abs(T12)
        rows.append((P, Tp, rel, resid))
    cons = [r for r in rows if r[2] < 0.15]
    print("   P grid (showing P where T_pred within 15% of T, plus reference anchors):")
    for P in (9, 12, 18, 22, 26):
        r = min(rows, key=lambda x: abs(x[0] - P))
        print(f"     P={r[0]:>5.1f}: T_pred={r[1]:>+.3e}  rel_err_vs_T={r[2]:>7.1%}  fit_resid={r[3]:.1e}")
    if cons:
        Ps = [c[0] for c in cons]
        print(f"   CONSISTENT-P set (|T_pred-T|/|T|<15%): P in [{min(Ps):.1f}, {max(Ps):.1f}]  ({len(Ps)} grid pts)")
    else:
        best = min(rows, key=lambda x: x[2])
        print(f"   NO P within 15%; best P={best[0]:.1f} at rel {best[2]:.1%}")
    print()

    # ============ SAT-C ============
    print("## SAT-C  ROBUSTNESS (windows r=7..11 and r=9..11; rho free in [0.95,1.0])")
    for win in ([7, 8, 9, 10, 11], [9, 10, 11]):
        best_by_P = []
        for i in range(int((30 - 6) / 0.5) + 1):
            P = 6 + 0.5 * i; th = 2 * math.pi / P
            bestrho = None
            for rr in np.linspace(0.95, 1.0, 51):
                a, b, resid = fit_ab(Lam, win, rr, th)
                Tp = tail(a, b, rr, th, 12)
                rel = abs(Tp - T12) / abs(T12)
                if bestrho is None or resid < bestrho[3]:
                    bestrho = (P, Tp, rel, resid, rr)
            best_by_P.append(bestrho)
        cons = [c for c in best_by_P if c[2] < 0.20]
        Ps = [c[0] for c in cons]
        rng = f"[{min(Ps):.1f}, {max(Ps):.1f}]" if Ps else "empty"
        print(f"   window r={win[0]}..{win[-1]}: consistent-P (rel<20%) = {rng}  ({len(Ps)} pts)")
    print()

    # ============ SAT-D ============
    print("## SAT-D  NEARLY MODEL-FREE half-cycle length from the tail")
    # tail of a damped alternating half-cycle sum ~ -(2/pi) h |Lam_13| / (1+rho^h); |Lam_13| ~ rho*|Lam_11|
    L13 = rho * abs(Lam[11])
    print(f"   |Lam_13| ~ rho*|Lam_11| = {L13:.3e};  T = {T12:.3e}.  solve T ~ -(2/pi) h L13/(1+rho^h) for h:")
    best = None
    for h10 in range(20, 300):
        h = h10 / 10
        pred = -(2 / math.pi) * h * L13 / (1 + rho ** h)
        d = abs(pred - T12)
        if best is None or d < best[2]:
            best = (h, pred, d)
    print(f"   h = {best[0]:.1f}  => predicted tail {best[1]:.3e} vs T {T12:.3e}  => PERIOD P = 2h = {2*best[0]:.1f}")
    print("   [Wilson pre-reg: h~11, P~22.  A result near P=9 refutes his arithmetic; near 20 refutes period-9.]\n")

    # ============ SAT-E ============
    print("## SAT-E  THE CAVEAT (stated): this test ASSUMES Sum Lam = -1/10 (i.e. S->7/15). It tests 7/15 and")
    print("   the period JOINTLY. Output = 'if 7/15, then P in [range]'. eps_12 is measured RELATIVE to 7/15, so")
    print("   the tail -eps_12/2 is conditional on 7/15. If the consistent-P range is implausible/empty, that is")
    print("   evidence against the single-mode model OR against 7/15 -- either is a finding.")


if __name__ == "__main__":
    main()
