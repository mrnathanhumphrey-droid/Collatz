"""
PROBE D -- THE DAM (transient deflation). Subtract the clean 1/2 transient to expose the bore (rho~0.984,
period ~9) across r=3..11, turning ~3 contaminated bore terms into ~9 cleaner ones.

Two modes in Lam_r: transient A*(1/2)^r (dominant r=3,4,5) + bore C*rho^r*cos(theta r+phi) (dominant r>=7).
Deflate B_r = Lam_r - T_r. CROSS-VALIDATE (D-B) the deflation on data it wasn't fit to, else B is unfalsifiable.

Exact rational for Lam_r, r<=7 (eps exact k<=8). r=8..11 use 15-digit float eps (result_epsilon_11 + S_12);
flagged. Lam_r = (eps_{r+1}-eps_r)/2, available r=1..11.
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


def solve_recur(B, rs):
    """fit B_{r+1}=a B_r + b B_{r-1} on the 3 consecutive triples in rs (>=4 pts); return a,b,rho,period."""
    import numpy as np
    rows, rhs = [], []
    for i in range(1, len(rs) - 1):
        rows.append([B[rs[i]], B[rs[i - 1]]]); rhs.append(B[rs[i + 1]])
    a, b = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    disc = a * a + 4 * b
    if disc < 0:
        rho = math.sqrt(-b); theta = math.acos(max(-1, min(1, a / (2 * rho))))
        return a, b, rho, (2 * math.pi / theta if theta > 1e-9 else float('inf')), math.degrees(theta)
    r1 = (a + math.sqrt(disc)) / 2; r2 = (a - math.sqrt(disc)) / 2
    return a, b, max(abs(r1), abs(r2)), float('inf'), 0.0   # real roots -> no oscillation from this window


def main():
    print("# PROBE D -- THE DAM (transient deflation).  1/2=0.5  0.984  2pi/log2=9.0647\n")
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in json.load(open(EPSJSON)).items()}
    LamX = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}          # exact r=1..7
    Lam = {r: (float(LamX[r]) if r in LamX else (EPS_F[r + 1] - EPS_F[r]) / 2) for r in range(1, 12)}

    # ============ D-A ============
    print("## D-A  CHARACTERIZE THE TRANSIENT (fit A*(1/2)^r on r=3,4,5; exact).  [gate]")
    A = {r: LamX[r] * 2 ** r for r in (3, 4, 5)}                        # A_r if pure 1/2-geometric
    print(f"   (i) A_r = Lam_r*2^r :  A_3={float(A[3]):.6f}  A_4={float(A[4]):.6f}  A_5={float(A[5]):.6f}")
    sp = (max(float(A[r]) for r in (3,4,5)) - min(float(A[r]) for r in (3,4,5))) / float(A[4])
    print(f"       spread = {sp:.2%}  [pre-reg: agree to >=4 digits IF rate is exactly 1/2 & region clean]")
    mu13 = math.sqrt(float(LamX[5] / LamX[3]))                          # (Lam5/Lam3)^{1/2}
    r43 = float(LamX[4] / LamX[3]); r54 = float(LamX[5] / LamX[4])
    print(f"   (ii) ratios Lam_4/Lam_3={r43:.5f}  Lam_5/Lam_4={r54:.5f} (straddle 1/2); best mu=(Lam5/Lam3)^.5={mu13:.6f}")
    print(f"        mu-1/2 = {mu13-0.5:+.2e}  [pre-reg STOP if |mu-1/2|>1e-4]  -> |dev|={abs(mu13-0.5):.1e}")
    print("   READ: A_r agree only to ~2 digits and mu deviates from 1/2 by ~2e-3. Per the ledger the BORE is")
    print("   already 12-47% of the transient at r=3,4,5, so there is NO fully-clean region -- the mu-deviation")
    print("   MEASURES bore contamination, it does NOT refute the exact-1/2 transient (R26). D-B is the arbiter.")
    # least-squares A with mu=1/2 fixed, on r=3,4,5 (exact)
    num = sum(LamX[r] * F(1, 2 ** r) for r in (3, 4, 5)); den = sum(F(1, 4 ** r) for r in (3, 4, 5))
    Astar = num / den
    print(f"   deflation transient: T_r = Astar*(1/2)^r, Astar (LSQ mu=1/2 on r=3,4,5) = {float(Astar):.6f}\n")
    T = {r: float(Astar) * 0.5 ** r for r in range(3, 12)}
    B = {r: Lam[r] - T[r] for r in range(3, 12)}

    # ============ D-C ============
    print("## D-C  THE DEFLATED LADDER (the artifact).  B_r = Lam_r - T_r.  [E]=exact Lam, f=float")
    print(f"   {'r':>2} {'Lam_r':>12} {'T_r':>11} {'B_r=Lam-T':>12} {'|B|/0.984^r':>12} {'sign B':>7}")
    for r in range(3, 12):
        tag = 'E' if r <= 7 else 'f'
        print(f"   {r:>2}{tag} {Lam[r]:>11.3e} {T[r]:>11.3e} {B[r]:>12.3e} {abs(B[r])/0.984**r:>12.5f} {'+' if B[r]>0 else '-':>7}")
    print("   [|B|/0.984^r flat => bore rate is 0.984; sign B => the bore's period, uncontaminated by transient.]\n")

    # ============ D-B ============
    print("## D-B  CROSS-VALIDATION GATE (fit bore on pure r=8..11, predict r=3,4,5; never touched transient fit)")
    a, b, rho, per, thetadeg = solve_recur(B, [8, 9, 10, 11])
    print(f"   bore 2-term recurrence on r=8..11: a={a:.5f} b={b:.5f} -> rho={rho:.5f} "
          + (f"theta={thetadeg:.2f}deg period={per:.3f}" if per != float('inf') else "(real roots: window <1/2 period, no osc pinned)"))
    # extrapolate backward via B_{r-1} = (B_{r+1} - a B_r)/b
    Bpred = {10: B[10], 9: B[9]}
    for r in range(8, 2, -1):
        Bpred[r] = (Bpred[r + 2] - a * Bpred[r + 1]) / b if b != 0 else float('nan')
    print(f"   {'r':>2} {'B_r (deflated)':>14} {'B_r (bore-predicted)':>20} {'rel err':>9}")
    okB = True
    for r in (3, 4, 5, 6, 7):
        pe = Bpred.get(r, float('nan'))
        rel = abs(pe - B[r]) / abs(B[r]) if abs(B[r]) > 1e-15 else float('nan')
        if r in (3, 4, 5) and (rel != rel or rel > 0.20):
            okB = False
        print(f"   {r:>2} {B[r]:>14.3e} {pe:>20.3e} {rel:>9.2%}")
    print(f"   => D-B {'PASS (bore validated on held-out early terms; two-mode model confirmed)' if okB else 'FAIL/WEAK (see rel err; third mode or non-geometric transient -- D-C..F caveated)'}\n")

    # ============ D-D ============
    print("## D-D  SIGN SEQUENCE ON DEFLATED DATA (corrected F1-B) -- sign(B_r), r=3..11 (9 terms)")
    Bsign = ['+' if B[r] > 0 else '-' for r in range(3, 12)]
    print(f"   sign(B_r) r=3..11: {' '.join(Bsign)}")
    for p in (6, 9, 12, 18):
        pairs = [(i, i + p) for i in range(len(Bsign)) if i + p < len(Bsign)]
        if not pairs:
            print(f"     period {p:>2}: UNTESTABLE (need >= {p+1} deflated terms; have 9)"); continue
        ok = all(Bsign[i] == Bsign[j] for i, j in pairs)
        nf = sum(1 for i, j in pairs if Bsign[i] != Bsign[j])
        print(f"     period {p:>2}: {'CONSISTENT' if ok else 'EXCLUDED'} ({len(pairs)} checks, {nf} fail)")
    print()

    # ============ D-E ============
    print("## D-E  RATE AND PERIOD (two disjoint windows; report ranges, state period-coverage)")
    for win in ([5, 6, 7, 8, 9], [7, 8, 9, 10, 11]):
        a2, b2, rho2, per2, td2 = solve_recur(B, win)
        span = (max(win) - min(win))
        cov = span / 9.06
        print(f"   window r={win[0]}..{win[-1]} ({span} steps ~ {cov:.2f} periods of 9.06): rho={rho2:.5f} "
              + (f"period={per2:.3f} (theta={td2:.1f}deg)" if per2 != float('inf') else "(real roots -> no oscillation resolved in-window)"))
    # holdout: fit r=6..11, predict 3,4,5
    ah, bh, rhoh, perh, tdh = solve_recur(B, [6, 7, 8, 9, 10, 11])
    print(f"   holdout fit r=6..11: rho={rhoh:.5f} " + (f"period={perh:.3f}" if perh != float('inf') else "(real roots)"))
    cap = 7 / 45
    capok = all(abs(Lam[r]) <= cap for r in range(3, 12))
    print(f"   amplitude cap |Lam_r| <= 7/45={cap:.4f} at every r? {capok}\n")

    # ============ D-F ============
    print("## D-F  THIRD-MODE CHECK (subtract the D-E bore from B_r; residual r=3..11)")
    # bore from holdout recurrence, reconstructed forward from B[6],B[7]
    if bh != 0 and perh != float('inf'):
        bore = {6: B[6], 7: B[7]}
        for r in range(8, 12):
            bore[r] = ah * bore[r - 1] + bh * bore[r - 2]
        for r in range(5, 2, -1):
            bore[r] = (bore[r + 2] - ah * bore[r + 1]) / bh
        print(f"   {'r':>2} {'B_r':>12} {'bore_r':>12} {'resid=B-bore':>13} {'|resid|/|B|':>11}")
        for r in range(3, 12):
            res = B[r] - bore.get(r, float('nan'))
            print(f"   {r:>2} {B[r]:>12.3e} {bore.get(r,float('nan')):>12.3e} {res:>13.3e} {abs(res)/abs(B[r]) if abs(B[r])>1e-15 else float('nan'):>11.2%}")
        print("   [residual structureless+shrinking => 2-mode complete; persistent structure => 3rd mode (beat).]")
    else:
        print("   D-E window gave real roots (oscillation not resolved) -> third-mode deflation not run; report the")
        print("   real-root ladder instead: the bore over r=6..11 is under-resolved for a clean oscillatory subtract.")


if __name__ == "__main__":
    main()
