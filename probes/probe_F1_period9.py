"""
PROBE F1 -- THE PERIOD-9 FRONTIER.
The asymptotic 2nd-moment rate of S_n->7/15 is the critical-only period-9 mode (~0.984), not 1/2 (R81).
This probe pins its structure: integer-vs-irrational period (structural, not decimal), where the mode lives,
and the lattice check that (per Lapidus-Hung) FORCES the oscillation in the p-adic/renewal setting.

Data: exact eps_k through k=8 (JSON); float eps_k k=9..12 (result_epsilon_11.csv + S_12). Lam_r=(eps_{r+1}-eps_r)/2.
"""
import os, sys, json, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from gate_M_reality_760 import build_K_float, stationary_float

HERE = os.path.dirname(os.path.abspath(__file__))
EPSJSON = os.path.join(HERE, '..', 'experiments_output', 'result_77_7_eps_exact_through_k8_v2_vec_pool.json')

# float eps_k, k=1..12 (k<=8 also exact; k=9..12 float from result_epsilon_11.csv + S_12)
EPS_F = {1: 0.2, 2: 9.523809523809525e-3, 3: -5.091986325893010e-3, 4: -2.452258248318762e-3,
         5: -1.151746915130986e-3, 6: -4.979056652200001e-4, 7: -1.175236830400000e-3,
         8: -7.455463672900000e-4, 9: -7.520257156400000e-6, 10: 7.207509171100000e-4,
         11: 1.501967012082273e-3, 12: 2.274713720558208e-3}


def main():
    print("# PROBE F1 -- THE PERIOD-9 FRONTIER.\n")
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in json.load(open(EPSJSON)).items()}  # exact k<=8
    LamX = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}          # exact Lam r=1..7
    LamF = {r: (EPS_F[r + 1] - EPS_F[r]) / 2 for r in range(1, 12)}     # float Lam r=1..11
    RHO = 0.984
    ROOT = math.sqrt(0.97008)                                          # F1-D candidate

    # ===================== F1-A =====================
    print("## F1-A  THE EXACT LADDER (single reference artifact; exact k<=8, float k=9..12).  [E]=exact")
    print(f"   {'k':>2} {'eps_k':>13} {'|eps|*2^k':>10} {'|eps|/.984^k':>12} {'|eps|/.98493^k':>14} {'sign':>5} || "
          f"{'r':>2} {'Lam_r':>12} {'|Lam|*2^r':>10} {'Lam_{r+1}/Lam_r':>15} {'sign':>5}")
    for k in range(2, 13):
        e = EPS_F[k]; tag = 'E' if k <= 8 else 'f'
        e2 = abs(e) * 2 ** k; er = abs(e) / RHO ** k; er2 = abs(e) / ROOT ** k
        rr = ""
        if k in LamF:
            lam = LamF[k]; l2 = abs(lam) * 2 ** k
            lr = (LamF[k + 1] / LamF[k]) if (k + 1) in LamF else float('nan')
            rr = f"{k:>2} {lam:>12.3e} {l2:>10.5f} {lr:>15.5f} {'+' if lam > 0 else '-':>5}"
        print(f"   {k:>2}{tag} {e:>12.3e} {e2:>10.5f} {er:>12.5f} {er2:>14.5f} {'+' if e > 0 else '-':>5} || {rr}")
    print("   [|eps|*2^k flat => rate 1/2; |eps|/.984^k flat => rate .984; .98493=sqrt(0.97008) is F1-D.]\n")

    # ===================== F1-B =====================
    print("## F1-B  INTEGER-vs-IRRATIONAL PERIOD (structural: sign-sequence periodicity, exact-robust)")
    Lsign = [('+' if LamF[r] > 0 else '-') for r in range(1, 12)]      # Lam_1..11 signs
    Esign = [('+' if EPS_F[k] > 0 else '-') for k in range(1, 13)]     # eps_1..12 signs
    print(f"   Lam sign seq r=1..11: {' '.join(Lsign)}   (Lam_4/Lam_3, Lam_5/Lam_4 are the ~1/2 R29-D read)")
    print(f"   eps sign seq k=1..12: {' '.join(Esign)}   (eps flips at k=3 and k=10)")
    print("   test integer period p: sign(Lam_{r+p}) == sign(Lam_r) for all valid r:")
    for p in (2, 3, 4, 5, 6, 7, 8, 9):
        pairs = [(r, r + p) for r in range(1, 12) if r + p <= 11]
        if not pairs:
            print(f"     p={p:>2}: UNTESTABLE (need >= {p+1} terms)"); continue
        ok = all(Lsign[r - 1] == Lsign[rp - 1] for r, rp in pairs)
        nfail = sum(1 for r, rp in pairs if Lsign[r - 1] != Lsign[rp - 1])
        print(f"     p={p:>2}: {'CONSISTENT' if ok else 'EXCLUDED'} ({len(pairs)} checks, {nfail} fail)")
    for p in (12, 18):
        print(f"     p={p:>2}: UNTESTABLE (need >= {p+1} terms; have 11)")
    print(f"   lattice prediction: 2*pi/log2 = {2*math.pi/math.log(2):.4f} (IRRATIONAL) -- an irrational period")
    print("   produces NO exact integer sign-period, consistent with all integers 2..9 excluded above.\n")

    # ===================== F1-C =====================
    print("## F1-C  WHERE THE MODE LIVES -- the lift tower M_n(1+3^j)/S_n (float pi), n=2..7")
    pis = {}
    for k in range(1, 8):
        K, cop, N = build_K_float(k); pis[k] = (stationary_float(K), cop, N)
    def M_of(n, eta):
        pi, cop, N = pis[n]
        pf = np.zeros(N)
        for i, r in enumerate(cop):
            pf[r] = pi[i]
        mh = np.fft.fft(pf); g = mh.copy(); g[0::3] = 0.0
        idx = (np.arange(N) * (eta % N)) % N
        return float(np.sum(g * np.conj(mh[idx])).real)
    towers = {}
    for n in range(2, 8):
        Sn = M_of(n, 1)
        row = []
        for j in range(n):
            eta = 1 + 3 ** j
            row.append(M_of(n, eta) / Sn)
        towers[n] = row
        cells = "  ".join(f"j{j}:{row[j]:+.5f}" for j in range(n))
        print(f"   n={n} S_n={Sn:.6f}  M_n(1+3^j)/S_n:  {cells}")
    print("   (highest j=n-1 is -1/2 EXACT = the rank-1 eta=1 part; lower j = candidate independent modes.)")
    print("   fixed-j convergence (does a lower-j column carry an independent ~0.984 mode?):")
    for j in (0, 1, 2):
        seq = [towers[n][j] for n in range(max(2, j + 2), 8)]  # M_n(1+3^j)/S_n across n
        difs = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
        rat = [difs[i + 1] / difs[i] if abs(difs[i]) > 1e-15 else float('nan') for i in range(len(difs) - 1)]
        print(f"     j={j}: ratio/S_n seq={['%.5f'%x for x in seq]}  successive-diff ratios={['%.3f'%x for x in rat]}")
    print()

    # ===================== F1-D =====================
    print("## F1-D  THE 0.970/0.984 RELATION (numerology cell -- NOT a finding, needs a mechanism)")
    print(f"   first-moment shifted mode (W1/W3) = 0.97008 (stable).  sqrt(0.97008) = {ROOT:.6f}")
    print(f"   second-moment envelope rate ~ 0.984 (corpus).  |sqrt(0.97008) - 0.984| = {abs(ROOT-0.984):.5f}")
    print("   => sqrt(1st-moment mode) ~ 2nd-moment rate to ~1e-3; can't confirm/deny at current precision.")
    print("   FLAGGED: this is the numerology class that killed 5 leads. Probe cell, not a finding.\n")

    # ===================== F1-E =====================
    print("## F1-E  LATTICE CHECK (Feller archimedean vs Lapidus-Hung p-adic; R81 claim, corrected)")
    l2 = math.log(2); l3 = math.log(3)
    # Feller-lattice test: exists d>0 with ALL (log3 - v log2) in d*Z (through origin)?
    #   <=> log3 in d*Z AND log2 in d*Z  <=>  log3/log2 in Q.  It is not:
    r32 = l3 / l2
    print(f"   ARCHIMEDEAN step: log M = log3 - v*log2.  Feller-lattice needs some d>0 with all values in d*Z,")
    print(f"   i.e. log3 and log2 both multiples of d  <=>  log3/log2 in Q.  log3/log2 = {r32:.12f}.")
    print(f"   log3/log2 is IRRATIONAL (elementary: 3^q=2^p impossible by unique factorization) => log3,log2")
    print(f"   INCOMMENSURATE => NO common d => the archimedean renewal step is NON-LATTICE. [Wilson's premise CORRECT]")
    print(f"   Feller/renewal theorem: a NON-lattice step gives NO archimedean log-periodic oscillation.")
    print(f"   (The 'coset log3+log2*Z' reading is NOT Feller-lattice -- it needs the origin, which incommensurability blocks.)")
    print(f"   YET period-~9 IS measured => the oscillation is NOT archimedean in origin. It comes from the 3-adic")
    print(f"   LEVEL structure -- integer levels r, trivially lattice. Lapidus-Hung: every p-adic self-similar string")
    print(f"   is LATTICE (zeta rational, poles periodic) => the oscillation is STRUCTURALLY FORCED by the p-adic")
    print(f"   setting and FORBIDDEN in the archimedean analogue. Its period is set by the level-transfer operator's")
    print(f"   arg(lambda_2), NOT by 2*pi/log2.  [2*pi/log2 = {2*math.pi/l2:.4f}, 2*pi/log3 = {2*math.pi/l3:.4f} -- suggestive vs")
    print(f"   measured ~9.2-9.5, but UNEARNED numerology like F1-D until a mechanism ties arg(lambda_2) to either.]")


if __name__ == "__main__":
    main()
