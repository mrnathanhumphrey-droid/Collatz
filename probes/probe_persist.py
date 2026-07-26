"""
PROBE PERSIST -- gate Wilson's level-1 dichotomy proof + measure the persistence-bound targets (2026-07-26).

LEVEL-1 PROOF (Wilson): 4 has mult order 3 mod 9  =>  3|k <=> 4^k = 1 mod 9 (trivial at 2nd tower level). In the
level-1 dlog coord 4 acts as a cyclic shift of order 3; nu_1=(0,1/3,2/3). Level-1 collision factor Sum_y nu_1(y)nu_1(4^k y):
  3|k -> trivial shift -> Sum nu_1^2 = 5/9    (diagonal);   3-nmid -> Sum rho(y)rho(y+-1) = 2/9   (off-diagonal).
=> gamma_1(k) = 3*factor = 5/3 (3|k) / 2/3 (else). Enrichment = Cauchy-Schwarz (diag >= off-diag), strict since nu_1
not shift-invariant. GATE: ord_9(4)=3; the two exact factors; matches banked gamma_1=[5/3,2/3,2/3,5/3].

PERSISTENCE as a BOUND not a sign: gamma_r(k) = gamma_1(k) * Prod_{j=2}^r (3 q_j(k)),  q_j(k)=p_j(k)/p_{j-1}(k).
Dichotomy survives to inf IFF  Prod_j 3q_j(1) < 3/2  (depleted stays <1)  AND  Prod_j 3q_j(3) > 3/5  (enriched stays >1).
log form: Sum_j log(3 q_j(k)) < log(3/2)=0.4055 (dep) / > log(3/5)=-0.5108 (enr). Since 3q_j=1+3(q_j-1/3),
the target is a two-sided bound on the RELAXATION Sum_j (q_j(k)-1/3) (signed) and Sum_j |q_j(k)-1/3| (total).
MEASURE these for k=1 (dep, rises to white), k=2 (dep, moves AWAY), k=3 (enr, falls to white) -- the targets Hank matches.

Reuses probe_channelfam rho-build. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from math import gcd, log
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact
from probe_channelfam import rho_exact_norm, C_ex

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def ord_mod(a, m):
    x = a % m; k = 1
    while x != 1:
        x = (x * a) % m; k += 1
        if k > m:
            return -1
    return k


def main():
    t0 = time.time()
    print("# PROBE PERSIST -- level-1 dichotomy proof + persistence-bound targets\n")

    # ---- level-1 gate ----
    print("## LEVEL-1 PROOF GATE")
    print(f"   ord_9(4) = {ord_mod(4, 9)}  (Wilson: 3)   [4,16=7,64=1 mod 9]")
    nex = build_nu_exact(5)
    r1 = rho_exact_norm(nex[1], 1)[0]
    v = [r1.get(s, F(0)) for s in range(3)]
    diag = sum(x * x for x in v)
    off = sum(v[y] * v[(y + 1) % 3] for y in range(3))
    print(f"   nu_1 (dlog) = {[str(x) for x in v]} ; Sum nu_1^2 = {diag} (=5/9? {diag==F(5,9)}) ; "
          f"cross = {off} (=2/9? {off==F(2,9)})")
    print(f"   gamma_1(k) = 3*factor: 3|k -> {3*diag} (=5/3? {3*diag==F(5,3)}) ; else -> {3*off} (=2/3? {3*off==F(2,3)})")
    print(f"   enrichment 5/9 >= 2/9 = Cauchy-Schwarz (diag>=off), strict since nu_1 not shift-invariant.\n")

    # ---- build rho ----
    nus = build_nu(0.5, 11)
    rho = {}
    for r in range(1, 12):
        N = 3 ** r
        mu = np.zeros(N)
        for X, w in nus[r].items():
            mu[(X - 1) // 3 % N] += float(w)
        d = R10.dlog_table(r)
        rr = np.zeros(N); rr[np.fromiter((d[a] for a in range(N)), np.int64, N)] = mu
        rho[r] = rr / rr.sum()
    del nus
    for r in range(12, 17):
        rr = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho[r] = rr / rr.sum()

    def p(r, k):
        return float(np.dot(rho[r], np.roll(rho[r], -k)))
    gam1 = {1: 2.0 / 3, 2: 2.0 / 3, 3: 5.0 / 3}

    # ---- persistence ladders ----
    print("## PERSISTENCE: q_j(k)=p_j(k)/p_{j-1}(k), relaxation (q_j-1/3), products vs thresholds")
    for k in (1, 2, 3):
        q = {j: p(j, k) / p(j - 1, k) for j in range(2, 17)}
        dev = {j: q[j] - 1.0 / 3 for j in q}
        prod = 1.0
        for j in range(2, 17):
            prod *= 3 * q[j]
        gam_r = gam1[k] * prod
        sig = sum(dev.values()); tot = sum(abs(d) for d in dev.values())
        thr = "< 3/2 (=1.5)" if k % 3 else "> 3/5 (=0.6)"
        prod_ratio = gam_r / gam1[k]
        inside = (prod_ratio < 1.5) if k % 3 else (prod_ratio > 0.6)
        print(f"   k={k} ({'ENR' if k%3==0 else 'dep'}): gamma_1={gam1[k]:.4f} -> gamma_16={gam_r:.5f} ; "
              f"Prod 3q_j = {prod_ratio:.5f}  needs {thr}: {inside}")
        print(f"        signed Sum(q_j-1/3) = {sig:+.5f} ; total Sum|q_j-1/3| = {tot:.5f} ; "
              f"log Prod 3q = {log(prod_ratio):+.5f} (thr {'log1.5=+0.405' if k%3 else 'log0.6=-0.511'})")
        # per-level deviation head (shows direction + decay)
        head = " ".join(f"{dev[j]:+.4f}" for j in range(2, 9))
        print(f"        (q_j-1/3) j=2..8: {head} ...")
    print()

    # ---- the concrete two-sided target for the pen ----
    print("## PEN TARGET (two-sided persistence inequality, measured room to threshold)")
    for k in (1, 3):
        q = {j: p(j, k) / p(j - 1, k) for j in range(2, 17)}
        prod = 1.0
        for j in range(2, 17):
            prod *= 3 * q[j]
        lp = log(prod)
        if k % 3:
            print(f"   k=1 depleted: log Prod 3q_j = {lp:+.5f}  vs ceiling log(3/2)=+0.40546  "
                  f"room = {0.40546 - lp:+.5f}  (tail beyond r=16 must not add > this)")
        else:
            print(f"   k=3 enriched: log Prod 3q_j = {lp:+.5f}  vs floor log(3/5)=-0.51083  "
                  f"room = {lp - (-0.51083):+.5f}  (tail must not subtract > this)")
    print("   => Hank's re-task: bound Sum_j |q_j(k)-1/3| (upper-bound / relaxation type, NOT sign) --")
    print("      Tao 1.14/1.17, BGK, Heilbronn, Bourgain are uniform |nu-hat| UPPER bounds = right type at last.")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
