"""
PROBE R27 -- THE RECURRENCE FIT. Reuses R25 builder for R27-C. Exact-rational for A/B/D.

Wilson's amendment: fit the two-term recurrence on Lambda_r, NOT on S_r or p_r=S_r/rho^r. S_r/p_r carry the LEADING
mode (constant 7/15 at criticality, eigenvalue 1) which contaminates a 2-term fit -- the plausible source of R26's
spurious 1.12. Lambda_r -> 0 by construction (difference kills the constant), so it carries only SUBDOMINANT modes.
And Lambda_1..Lambda_7 are EXACT rationals => solve the 2x2 exactly, no drift.

  Lambda_r = (S_{r+1}-S_r)/2 = (eps_{r+1}-eps_r)/2  (eps_k=S_k-7/15, exact through k=8). Signs -,-,+,+,+,-,+.
 A: 2-term Lambda_{r+1}=a Lambda_r + b Lambda_{r-1}, four exact solves (pairs {2,3}/{3,4}/{4,5}/{5,6}); roots of
    lam^2-a lam-b=0; |lam2|=sqrt|b| (pre-reg 0.475-0.5); 2pi/arg (pre-reg ~9, EXPECTED TO FAIL vs R26 >16).
 B: 3-term order test (residual structure): exact recurrence? else |lam3|=residual rate.
 C: kappa=(|lam2|/rho)/(2lam^2) at eps=0.05,0.075,0.08,0.1 -> monotone to 1?  (2-term fit on subcritical Lambda_r).
 D: sign-pattern cross-check: propagate the recurrence from Lambda_1,Lambda_2, compare signs to banked -,-,+,+,+,-,+.
"""
import os, sys, math, json, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np


def solve2(L, i):
    """Exact 2-term solve from equations at r=i and r=i+1: L[i+1]=a L[i]+b L[i-1], L[i+2]=a L[i+1]+b L[i]."""
    m00, m01 = L[i], L[i - 1]
    m10, m11 = L[i + 1], L[i]
    det = m00 * m11 - m01 * m10
    a = (L[i + 1] * m11 - L[i + 2] * m01) / det
    b = (m00 * L[i + 2] - m10 * L[i + 1]) / det
    return a, b


def roots_of(a, b):
    """roots of lam^2 - a lam - b = 0 (a,b float)."""
    disc = a * a + 4 * b
    if disc >= 0:
        return [(a + math.sqrt(disc)) / 2, (a - math.sqrt(disc)) / 2]
    return [complex(a / 2, math.sqrt(-disc) / 2), complex(a / 2, -math.sqrt(-disc) / 2)]


def build_mu_qf(arr, k, q, lam, tol=1e-18):
    M = q ** k; inv2 = pow(2, -1, M)
    a_idx = np.nonzero(arr)[0]; a_val = arr[a_idx]
    base = (1 + q * a_idx) % M
    mu = np.zeros(M); u = inv2; v = 1
    while (1 - lam) * lam ** (v - 1) > tol:
        wv = (1 - lam) * lam ** (v - 1)
        mu += np.bincount((u * base) % M, weights=wv * a_val, minlength=M)
        u = (u * inv2) % M; v += 1
    return mu


def shells(q, lam, RMAX):
    arr = np.array([1.0]); Y = {0: 1.0}
    for k in range(1, RMAX + 1):
        arr = build_mu_qf(arr, k, q, lam)
        Y[k] = q ** k * float(np.sum(arr * arr))
    return {r: Y[r] - Y[r - 1] for r in range(1, RMAX + 1)}


def main():
    print("# PROBE R27 -- THE RECURRENCE FIT (on Lambda_r, exact).\n")
    hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                        'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in hist.items()}      # eps_k exact through k=8
    Lam = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}                     # Lambda_1..7 exact
    print("## Lambda_r (exact, = (eps_{r+1}-eps_r)/2), signs should be -,-,+,+,+,-,+")
    for r in range(1, 8):
        print(f"   Lambda_{r} = {str(Lam[r])[:34]:>34} = {float(Lam[r]):+.6e}  ({'+' if Lam[r]>0 else '-'})")
    print()

    # ================= R27-A =================
    print("## R27-A  TWO-TERM RECURRENCE, EXACT: Lambda_{r+1}=a Lambda_r + b Lambda_{r-1}")
    L = {r: Lam[r] for r in range(1, 8)}
    sols = {}
    for i in (2, 3, 4, 5):                       # pairs {2,3},{3,4},{4,5},{5,6}
        a, b = solve2(L, i)
        sols[i] = (a, b)
        rts = roots_of(float(a), float(b))
        lam2 = rts[0]
        mag = abs(lam2)
        if isinstance(lam2, complex) and abs(lam2.imag) > 1e-12:
            per = 2 * math.pi / abs(cmath.phase(lam2))
            argd = math.degrees(cmath.phase(lam2))
        else:
            per, argd = float('inf'), 0.0
        print(f"   solve{{{i},{i+1}}}: a={float(a):+.6f} b={float(b):+.6f}  |lam2|=sqrt|b|={math.sqrt(abs(float(b))):.5f} "
              f"root|.|={mag:.5f}  arg={argd:+.2f}deg period={per:.3f}")
    # spread
    aa = [float(sols[i][0]) for i in (2, 3, 4, 5)]
    bb = [float(sols[i][1]) for i in (2, 3, 4, 5)]
    alleq = all(sols[i] == sols[2] for i in (3, 4, 5))
    print(f"   spread: a in [{min(aa):.5f},{max(aa):.5f}] (dev {max(aa)-min(aa):.2e}); "
          f"b in [{min(bb):.5f},{max(bb):.5f}] (dev {max(bb)-min(bb):.2e})")
    print(f"   ALL FOUR SOLVES EXACTLY EQUAL? {alleq}  {'<<< JACKPOT: exact finite 2-term recurrence' if alleq else '(no: 2-term approximate, spread = 3rd-mode leakage)'}")
    print(f"   [PRE-REG: |lam2|~0.475-0.5; period ~9 EXPECTED TO FAIL vs R26>16 -- a miss kills period-9 corpus-wide]\n")

    # ================= R27-B =================
    print("## R27-B  ORDER TEST: three-term Lambda_{r+1}=a Lambda_r+b Lambda_{r-1}+c Lambda_{r-2}")
    # solve exactly from equations at r=3,4,5 (uses Lambda_1..6), check residual at r=6 (Lambda_7)
    import itertools
    M3 = [[L[3], L[2], L[1]], [L[4], L[3], L[2]], [L[5], L[4], L[3]]]
    rhs3 = [L[4], L[5], L[6]]
    # solve 3x3 exactly (Cramer)
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D = det3(M3)
    abc = []
    for col in range(3):
        Mc = [row[:] for row in M3]
        for row in range(3):
            Mc[row][col] = rhs3[row]
        abc.append(det3(Mc) / D)
    a3, b3, c3 = abc
    resid = L[7] - (a3 * L[6] + b3 * L[5] + c3 * L[4])         # predict Lambda_7, exact residual
    print(f"   3-term (a,b,c) from r=3,4,5: a={float(a3):+.5f} b={float(b3):+.5f} c={float(c3):+.5f}")
    print(f"   residual predicting Lambda_7: {float(resid):+.3e}  (exact {'ZERO -> exact 3-term recurrence' if resid==0 else 'nonzero -> >=3 modes'})")
    # 2-term residual for comparison: use solve{2,3}, predict Lambda_5,6,7
    a2, b2 = sols[2]
    for rr in (5, 6, 7):
        pred = a2 * L[rr - 1] + b2 * L[rr - 2]
        print(f"   2-term (solve{{2,3}}) predict Lambda_{rr}: resid={float(L[rr]-pred):+.3e}")
    print()

    # ================= R27-D =================
    print("## R27-D  SIGN-PATTERN CROSS-CHECK: propagate recurrence from Lambda_1,Lambda_2 vs banked -,-,+,+,+,-,+")
    for tag, (a, b) in [("solve{2,3}", sols[2]), ("solve{4,5}", sols[4])]:
        seq = [L[1], L[2]]
        for r in range(3, 8):
            seq.append(a * seq[-1] + b * seq[-2])
        signs = "".join('+' if float(x) > 0 else '-' for x in seq)
        actual = "".join('+' if Lam[r] > 0 else '-' for r in range(1, 8))
        print(f"   {tag}: propagated signs {signs}  vs actual {actual}  [{'MATCH' if signs==actual else 'MISMATCH -> 3rd mode matters'}]")
    print()

    # ================= R27-C =================
    print("## R27-C  kappa=(|lam2|/rho)/(2lam^2) at eps=0.05,0.075,0.08,0.1  (2-term fit on subcritical Lambda_r)")
    print(f"   {'eps':>6} {'lam':>6} {'rho_pred':>9} {'rho(fit)':>9} {'|lam2|':>8} {'|lam2|/rho':>10} {'2lam^2':>8} {'kappa':>8}")
    for eps in (0.05, 0.075, 0.08, 0.1):
        lam = 0.5 + eps; rho = 3 * (1 - lam) / (1 + lam)
        S = shells(3, lam, 14)
        Lf = {r: (S[r + 1] - S[r]) / 2 for r in range(1, 14)}     # subcritical Lambda_r (float)
        # exact-style 2-term solve at deep r (avg over a few), roots -> rho(leading), lam2
        res = []
        for i in (8, 9, 10):
            m00, m01, m10, m11 = Lf[i], Lf[i - 1], Lf[i + 1], Lf[i]
            det = m00 * m11 - m01 * m10
            a = (Lf[i + 1] * m11 - Lf[i + 2] * m01) / det
            b = (m00 * Lf[i + 2] - m10 * Lf[i + 1]) / det
            res.append((a, b))
        a = np.mean([x[0] for x in res]); b = np.mean([x[1] for x in res])
        rts = roots_of(a, b)
        rts = sorted(rts, key=lambda z: -abs(z))
        rho_fit = rts[0].real if not isinstance(rts[0], complex) or abs(rts[0].imag) < 1e-9 else abs(rts[0])
        lam2 = rts[1]
        ratio = abs(lam2) / rho_fit
        kappa = ratio / (2 * lam * lam)
        print(f"   {eps:>6} {lam:>6.3f} {rho:>9.5f} {rho_fit:>9.5f} {abs(lam2):>8.5f} {ratio:>10.5f} {2*lam*lam:>8.5f} {kappa:>8.5f}")
    print("   [Q: is kappa monotone toward 1? three points beat two.]")


if __name__ == "__main__":
    main()
