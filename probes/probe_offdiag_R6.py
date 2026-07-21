"""
PROBE R6 -- THE INTERFERENCE LEDGER. R6-A gates the off-diagonal formalism at k=2 DIRECTLY from the
v != v' pair-character (Ramanujan) sums -- independent of the S-table.
Shell = primitive-frequency Plancherel: S_k = sum_{xi primitive mod 3^k} |mu_hat_k(xi)|^2
      = E_{X,X'}[ c_{3^k}(X - X') ],  c = Ramanujan sum (sum over xi coprime to 3).
At k=2: X = 2^{-v}(1+3Y) mod 9, v ~ Geom(1/2) (2^{-v} mod 9 has period 6), Y ~ pi_1 (mod 3).
DIAGONAL v=v': weight P(v=v') = sum 4^{-v} = 1/3 (infinite valuations, EXACT); phase collapses -> replicates S_1.
OFF-DIAGONAL v!=v': OffDiag_2, PRE-REGISTERED = -4/21 exactly.
R6-B ledger table (from frozen S). R6-C channel anatomy by |v-v'|.
"""
from fractions import Fraction as F

# 2^{-v} mod 9 for v in class j (v ≡ j mod 6, j=1..6; j=6 means v≡0):  [5,7,8,4,2,1]
INV9 = {1: 5, 2: 7, 3: 8, 4: 4, 5: 2, 6: 1}
# folded prob of class j:  W_j = sum_{v≡j mod6, v>=1} 2^{-v} = 2^{-j} * 64/63
W = {j: F(1, 2 ** j) * F(64, 63) for j in range(1, 7)}
# exact same-value (v=v') weight for class j:  D_j = sum_{v≡j} 4^{-v} = 4^{-j} * 4096/4095
D = {j: F(1, 4 ** j) * F(4096, 4095) for j in range(1, 7)}
PI1 = {1: F(1, 3), 2: F(2, 3)}                 # pi_1: residue 1 -> 1/3, residue 2 -> 2/3 (Judge One)


def c9(delta):
    d = delta % 9
    if d == 0: return 6                        # Ramanujan c_9: 9|delta
    if d % 3 == 0: return -3                   # 3|delta, 9 not | delta
    return 0                                    # coprime


def Xval(j, Y):
    return (INV9[j] * (1 + 3 * Y)) % 9         # X = 2^{-v}(1+3Y) mod 9


def R6A():
    print("## R6-A  GATE THE FORMALISM: OffDiag_2 direct from v!=v' Ramanujan sums (PRE-REG -4/21)")
    # DIAGONAL v=v': same class j (same 2^{-v}), independent Y,Y'; weight D_j
    diag = F(0)
    for j in range(1, 7):
        for Y in (1, 2):
            for Yp in (1, 2):
                diag += D[j] * PI1[Y] * PI1[Yp] * c9(Xval(j, Y) - Xval(j, Yp))
    # OFF-DIAGONAL v!=v'
    off = F(0)
    for j in range(1, 7):
        for jp in range(1, 7):
            w = (W[j] * W[jp] - D[j]) if j == jp else (W[j] * W[jp])   # v!=v' weight
            for Y in (1, 2):
                for Yp in (1, 2):
                    off += w * PI1[Y] * PI1[Yp] * c9(Xval(j, Y) - Xval(jp, Yp))
    print(f"   diagonal channel   = {diag} = {float(diag):+.8f}   (PRE-REG S_1 = 2/3;  {'REPLICATES' if diag==F(2,3) else 'DEV'})")
    print(f"   off-diagonal (v!=v')= {off} = {float(off):+.8f}   (PRE-REG -4/21;  {'GATE PASS' if off==F(-4,21) else 'DEV'})")
    print(f"   sum                = {diag+off} = {float(diag+off):+.8f}   (S_2 = 10/21;  {'OK' if diag+off==F(10,21) else 'DEV'})")
    print(f"   => diagonal weight = sum D_j = {sum(D.values())} (=1/3 exact); 3-to-1 lift x 1/3 = replication.")
    return off


def R6B():
    print("\n## R6-B  THE LEDGER TABLE (OffDiag_k = S_k - S_(k-1); frozen S from Basic.lean)")
    S = {1: F(2, 3), 2: F(10, 21), 3: F(31370, 67963),
         4: F(143195649659456490, 308468774477179141),
         5: F(2490699741144069281815149277465294323722281911695597204406570,
              5350418720142111510029542161258891403960563152740894082816203)}
    print(f"   {'k':>2} {'OffDiag_k':>14} {'float':>12} {'running total':>16} {'vs -1/5':>10}")
    run = F(0)
    for k in range(2, 6):
        od = S[k] - S[k - 1]; run += od
        print(f"   {k:2d} {str(od)[:14]:>14} {float(od):>+12.6f} {float(run):>+16.8f} {float(run - F(-1,5)):>+10.6f}")
    print(f"   Sum_(k>=2) OffDiag -> S_inf - S_1 = {F(7,15)-S[1]} = -1/5 (the constant's target);  "
          f"Sum_(k>=3) = {F(7,15)-S[2]} = -1/105")
    print(f"   signs: {['-' if S[k]-S[k-1] < 0 else '+' for k in range(2,6)]}  (overshoot: --++, k>=4 tail = "
          f"{F(7,15)-S[3]} = +{float(F(7,15)-S[3]):.5f} NET POSITIVE)")


def R6C(cap=54):
    print(f"\n## R6-C  CHANNEL ANATOMY: off-diagonal by |v-v'| (exact rationals, v,v' <= {cap})")
    # class of exact valuation v: ((v-1) mod 6)+1 ; 2^{-v} mod 9 = INV9[class]
    def cls(v): return ((v - 1) % 6) + 1
    bins = {1: F(0), 2: F(0), "3+": F(0)}
    for v in range(1, cap + 1):
        for vp in range(1, cap + 1):
            if v == vp: continue
            d = abs(v - vp)
            w = F(1, 2 ** v) * F(1, 2 ** vp)             # P(v)P(v') unnormalized (Geom(1/2))
            s = F(0)
            for Y in (1, 2):
                for Yp in (1, 2):
                    s += PI1[Y] * PI1[Yp] * c9(Xval(cls(v), Y) - Xval(cls(vp), Yp))
            key = d if d <= 2 else "3+"
            bins[key] += w * s
    tot = sum(bins.values())
    print(f"   |v-v'|=1 : {float(bins[1]):+.8f}")
    print(f"   |v-v'|=2 : {float(bins[2]):+.8f}")
    print(f"   |v-v'|>=3: {float(bins['3+']):+.8f}")
    print(f"   total off-diagonal (v,v'<= {cap}) = {float(tot):+.8f}   (-> -4/21 = {float(F(-4,21)):+.8f})")
    print(f"   [the two-sign tail structure: which |v-v'| bins carry which sign -- raw material]")


def pi2():
    """level-2 measure mod 9: pi_2(r) = sum_{j1,Y} W_{j1} PI1[Y] 1[Xval(j1,Y)=r]."""
    p = {}
    for j in range(1, 7):
        for Y in (1, 2):
            r = Xval(j, Y)
            p[r] = p.get(r, F(0)) + W[j] * PI1[Y]
    return p


def R6C_k3(cap=54):
    print(f"\n## R6-C (k=3)  channel anatomy of OffDiag_3 by |v_1-v'_1| (exact, v<= {cap}; c_27 Ramanujan)")
    inv27 = {}
    x = 1
    inv2_27 = pow(2, -1, 27)
    for v in range(1, 19):
        x = (x * inv2_27) % 27
        inv27[v] = x                            # 2^{-v} mod 27, class v mod 18 (v=18 -> class 18)
    P2 = pi2()

    def c27(delta):
        d = delta % 27
        if d == 0: return 18
        if d % 9 == 0: return -9
        if d % 3 == 0: return 0                 # (c_27 = 0 unless 9|delta)
        return 0

    def X3(cls18, x2):
        return (inv27[cls18] * (1 + 3 * x2)) % 27
    def cl18(v): return ((v - 1) % 18) + 1

    # diagonal v_1=v'_1 (should replicate S_2 = 10/21)
    D18 = {j: F(1, 4 ** j) * F(2 ** 36, 2 ** 36 - 1) for j in range(1, 19)}   # sum_{v≡j mod18} 4^{-v}
    diag = F(0)
    for j in range(1, 19):
        for x2, w2 in P2.items():
            for x2p, w2p in P2.items():
                diag += D18[j] * w2 * w2p * c27(X3(j, x2) - X3(j, x2p))
    bins = {1: F(0), 2: F(0), "3+": F(0)}
    for v in range(1, cap + 1):
        for vp in range(1, cap + 1):
            if v == vp: continue
            wv = F(1, 2 ** v) * F(1, 2 ** vp)
            s = F(0)
            for x2, w2 in P2.items():
                for x2p, w2p in P2.items():
                    s += w2 * w2p * c27(X3(cl18(v), x2) - X3(cl18(vp), x2p))
            key = abs(v - vp) if abs(v - vp) <= 2 else "3+"
            bins[key] += wv * s
    off = sum(bins.values())
    print(f"   diagonal channel = {diag} = {float(diag):+.8f}   (PRE-REG S_2 = 10/21; {'REPLICATES' if diag==F(10,21) else f'DEV {float(diag-F(10,21)):+.1e}'})")
    print(f"   off-diag by |v_1-v'_1|:  =1: {float(bins[1]):+.8f}   =2: {float(bins[2]):+.8f}   >=3: {float(bins['3+']):+.8f}")
    print(f"   total off-diag (v<= {cap}) = {float(off):+.8f}   (-> OffDiag_3 = -2980/203889 = {float(F(-2980,203889)):+.8f})")
    print(f"   [two-sign raw material: sign per |v-v'| bin at k=3]")


def main():
    print("# PROBE R6 -- THE INTERFERENCE LEDGER. Exact rationals.")
    R6A()
    R6B()
    R6C()
    R6C_k3()


if __name__ == "__main__":
    main()
