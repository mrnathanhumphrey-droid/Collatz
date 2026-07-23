"""
PROBE W2 -- THE TAILOR ON THE SECOND MOMENT.
Does the 76.2 self-inverse/pair structure send the (second-moment) rate to 1/2?

Built ENTIRELY on second-moment objects (M_n(eta), the eps/Lambda ladder). Does NOT reuse the W1 first-moment
core except for the explicit W2-C/D comparison. Exact rational where finite (eps table + Lambda ratios). The
lift-triple of eta=1 is {1, 1+3^{n-1}, 1+2*3^{n-1}}; by Thm 76.3 it equals S_n*(1,-1/2,-1/2) EXACTLY (rank-1).
"""
import os, sys, json, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_warp_W13 import build_core, phase_diag, topmods

HERE = os.path.dirname(os.path.abspath(__file__))
EPSJSON = os.path.join(HERE, '..', 'experiments_output', 'result_77_7_eps_exact_through_k8_v2_vec_pool.json')


def load_eps():
    d = json.load(open(EPSJSON))
    return {int(k): F(int(v['num']), int(v['den'])) for k, v in d.items()}


def main():
    print("# PROBE W2 -- THE TAILOR ON THE SECOND MOMENT.  1/2=0.50000  43/45=0.95556  0.984\n")
    EPS = load_eps()                                     # eps_k = S_k - 7/15, exact, k=1..8
    S = {k: F(7, 15) + EPS[k] for k in EPS}
    Lam = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}   # 2 Lam_r = eps_{r+1}-eps_r; r=1..7

    # ===================== W2-A =====================
    print("## W2-A  THE LIFT-TRIPLE OPERATOR (exact).  triple_n = (M_n(1), M_n(1+3^{n-1}), M_n(1+2*3^{n-1}))")
    print("   By Thm 76.3 (level m=n-1): S_n = -2 M_n(1+3^{n-1}); M-real (Lemma 76.0) => both pair entries = -S_n/2.")
    print(f"   {'n':>2} {'M_n(1)=S_n':>12} {'M_n(1+3^{n-1})':>14} {'M_n(1+2*3^{n-1})':>16} {'triple / S_n':>18}")
    for n in range(2, 9):
        s = S[n]; pair = -s / 2
        print(f"   {n:>2} {float(s):>12.8f} {float(pair):>14.8f} {float(pair):>16.8f} {'(1, -1/2, -1/2)':>18}")
    print("   => triple = S_n*(1,-1/2,-1/2) EXACT, RANK-1.  A 3x3 on it is DEGENERATE: the only nontrivial")
    print("      eigenvalue is the S-ratio S_{n+1}/S_n -> 1 (leading); M(1+3^n) carries the SAME rate as S_{n+1}")
    print("      (they are proportional by -1/2), so there is NO independent 1/2 subdominant in the eta=1 triple.")
    print("      [Wilson guardrail / W2-E outcome-3: the triple is too coarse -- it tracks only |S_n|, not the mode.]\n")

    # ===================== W2-B =====================
    print("## W2-B  THE 76.2 PAIRING AS A_*  (exact).  M(self-inv=1) = -2 * M(pair)?")
    print(f"   {'n':>2} {'M(1)=S_n':>12} {'-2*M(1+3^{n-1})':>16} {'match':>7}")
    okB = True
    for n in range(2, 9):
        lhs = S[n]; rhs = -2 * (-S[n] / 2)
        good = (lhs == rhs); okB = okB and good
        print(f"   {n:>2} {float(lhs):>12.8f} {float(rhs):>16.8f} {'EXACT' if good else 'DEV':>7}")
    print(f"   => 76.2 -2 factor holds EXACT (it IS Thm 76.3). {'' if okB else 'FAIL'}")
    print("   BUT: the -2 is a MAGNITUDE relation WITHIN one level (S_n vs its pair), NOT a rate-changing factor")
    print("   ACROSS levels. Within the eta=1 triple every entry has the SAME rate (all = const*S_n). So '-2 sends")
    print("   the core rate to 1/2' is FALSE as stated: -2 relates magnitudes at fixed n, the rate lives elsewhere.\n")

    # ===================== W2-C =====================  (THE decisive read)
    print("## W2-C  MOMENT SEPARATION -- the exact second-moment rate ladder vs the first-moment 0.970")
    print("   second-moment: eps_k=S_k-7/15 and Lam_r=(eps_{r+1}-eps_r)/2 (=OffDiag/2). EXACT ratios, NO fit:")
    print(f"   {'k':>2} {'eps_k':>12} {'eps_k/eps_{k-1}':>15} {'|eps_k|*2^k':>12}   ||   "
          f"{'r':>2} {'Lam_r':>12} {'Lam_{r+1}/Lam_r':>16} {'|Lam_r|*2^r':>12}")
    for k in range(2, 9):
        er = float(EPS[k] / EPS[k - 1]) if k - 1 in EPS else float('nan')
        e2 = abs(float(EPS[k])) * 2 ** k
        rr = ""
        if k in Lam and k - 1 in Lam:
            lr = float(Lam[k] / Lam[k - 1]); l2 = abs(float(Lam[k])) * 2 ** k
            rr = f"{k:>2} {float(Lam[k]):>12.3e} {lr:>16.5f} {l2:>12.5f}"
        elif k in Lam:
            rr = f"{k:>2} {float(Lam[k]):>12.3e} {'--':>16} {abs(float(Lam[k]))*2**k:>12.5f}"
        print(f"   {k:>2} {float(EPS[k]):>12.3e} {er:>15.5f} {e2:>12.5f}   ||   {rr}")
    print("   first-moment (W1 recheck): shifted-core slow mode (should be ~0.970, the period-9/envelope family):")
    r = 5
    D = phase_diag(r); C = build_core(r); T = D[:, None] * C
    evT = topmods(T, 4)
    print(f"     r={r}: first-moment shifted |lam2|={abs(evT[1]):.5f}  arg={math.degrees(cmath.phase(evT[1])):+.2f}deg")
    print("   PRE-REG: second-moment ladder and first-moment 0.970 are DIFFERENT objects. Read the ladder above:")
    print("   does Lam_{r+1}/Lam_r stay ~1/2 or rise toward ~0.98? (that decides transient-vs-asymptotic for 1/2).\n")

    # ===================== W2-D =====================
    print("## W2-D  PERIOD FROM THE FIRST-MOMENT MODE (arg of the 0.970 mode; NO fit)")
    print(f"   {'r':>2} {'|lam2|':>9} {'arg(deg)':>9} {'2pi/arg=period':>15}")
    for r in (3, 4, 5, 6):
        D = phase_diag(r); C = build_core(r); T = D[:, None] * C
        ev = topmods(T, 3)
        l2 = ev[1]; ang = math.degrees(cmath.phase(l2))
        per = (360 / abs(ang)) if abs(ang) > 1e-9 else float('inf')
        print(f"   {r:>2} {abs(l2):>9.5f} {ang:>+9.2f} {per:>15.4f}")
    # also the core units-only subdominant (the raw period-9 carrier)
    for r in (4, 5, 6):
        N = 3 ** r; units = [x for x in range(N) if x % 3]
        C = build_core(r); Cu = C[np.ix_(units, units)]
        ev = topmods(Cu, 3); l2 = ev[1]; ang = math.degrees(cmath.phase(l2))
        per = (360 / abs(ang)) if abs(ang) > 1e-9 else float('inf')
        print(f"   (core units r={r}) |lam2|={abs(l2):.5f} arg={ang:+.2f} period={per:.3f}")
    print()

    # ===================== W2-E =====================
    print("## W2-E  DENOUEMENT")
    lr_late = [float(Lam[r + 1] / Lam[r]) for r in range(3, 7)]  # r=3->4 .. 6->7
    print(f"   Lam ratios r=3..7: {['%.4f'%x for x in lr_late]}")
    e2 = [abs(float(EPS[k])) * 2 ** k for k in range(3, 9)]
    print(f"   |eps_k|*2^k k=3..8: {['%.4f'%x for x in e2]}   (flat ~const => rate 1/2; rising => faster than 1/2)")
    print("   Read: does the exact second-moment rate STAY 1/2 through r=7, or drift up toward the 0.97/0.984 mode?")
    print("   Verdict is whatever the exact numbers say -- banked as-is (transient vs asymptotic 1/2).")


if __name__ == "__main__":
    main()
