"""
PROBE HIERARCHY -- Wilson's v_3(k) hierarchy: gamma_j(k)=X_j for j<=v_3(k), then departs & relaxes (2026-07-26).

MECHANISM (one line of group theory): G_j = <4> has order 3^j, so 4^k=id in G_j  <=>  3^j | k  <=>  v_3(k)>=j.
For every level j<=v_3(k) the channel-k collision IS the m=0 event => gamma_j(k)=gamma_j(0)=X_j (EXACT, base case at
every valuation). At the autocorr level this is automatic: lag k on the 3^j-point profile wraps to 0 when 3^j|k.
Then for j>v_3(k) the channel departs X_j and relaxes; dichotomy = whether it stays on the correct side of white=1.

CHECKS:
 (1) EXACT: gamma_j(k)=X_j for j<=v_3(k), and gamma_{v_3+1}(k) != X (departs). k=3 (v=1), k=9 (v=2), k=27 (v=3).
     X_1=5/3, X_2=9*p_2(0)=9*(5/21)=15/7=2.1429, X_3 (exact).
 (2) PREDICTION: gamma_inf(27), v_3=3, tracks X_3 then relaxes -> land ~2.3-2.9 (X_3 minus a little).
 (3) GROUP by v_3: v=0 depleted, v=1 {3,6,12}, v=2 {9}, v=3 {27}. gamma_inf increases with v_3 (tracks X_{v_3}).
 (4) BUDGET (Wilson): relaxation Pi_{j>v_3} 3q_j from base X_{v_3} (v>=1) or gamma_1=2/3 (v=0); budget used =
     |log Pi_relax| / |threshold|, threshold = log(3/2) depleted (ceiling) / log(3^{v+1}? ...)=-log X_{v_3}... floor.
     Identify the BINDING channels (Wilson: k=4 63%, k=3 58%; k=1 loose 23%).

Reuses probe_channelfam rho-build; k up to 28. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from math import log
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact
from probe_channelfam import rho_exact_norm, C_ex

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"
KMAX = 28


def v3(k):
    v = 0
    while k % 3 == 0:
        k //= 3; v += 1
    return v


def main():
    t0 = time.time()
    print("# PROBE HIERARCHY -- gamma_j(k)=X_j for j<=v_3(k), then relaxes; k=9 exact check + k=27 prediction\n")

    # ---- exact X_j and hierarchy anchors (r<=5) ----
    nex = build_nu_exact(5)
    rex = {r: rho_exact_norm(nex[r], r)[0] for r in range(1, 6)}
    print("## (1) EXACT hierarchy: gamma_j(k) == X_j for j <= v_3(k)  (Fractions)")
    Xex = {}
    for j in range(1, 6):
        N = 3 ** j
        Xex[j] = F(3) ** j * C_ex(rex[j], N, 0)
    print(f"   X_1={Xex[1]}={float(Xex[1]):.4f}  X_2={Xex[2]}={float(Xex[2]):.4f}  X_3={Xex[3]}={float(Xex[3]):.4f}  "
          f"X_4={float(Xex[4]):.4f}  X_5={float(Xex[5]):.4f}")
    for k in (3, 9, 27):
        vk = v3(k)
        print(f"   k={k} (v_3={vk}):")
        for j in range(1, min(vk + 2, 6)):
            N = 3 ** j
            g = F(3) ** j * C_ex(rex[j], N, k % N)      # lag k mod 3^j
            eq = (g == Xex[j])
            tag = "== X_j (tracks m=0)" if j <= vk else "!= X_j (DEPARTS)"
            print(f"       gamma_{j}({k}) = {g} = {float(g):.4f}   X_{j}={float(Xex[j]):.4f}   {'EQ' if eq else 'NE'}  [{tag}]")
    print()

    # ---- float rho for r=16 relaxed values + k=27 ----
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

    def gam(r, k):
        return 3.0 ** r * float(np.dot(rho[r], np.roll(rho[r], -k)))
    Xf = {j: gam(j, 0) for j in range(1, 17)}

    # ---- (2) k=27 prediction ----
    print("## (2) PREDICTION k=27 (v_3=3): tracks X_3 then relaxes -> land ~2.3-2.9")
    g27_j = [gam(j, 27 % 3 ** j if 3 ** j <= 27 else 27) for j in (1, 2, 3)]
    print(f"   gamma_1..3(27) = {g27_j[0]:.4f},{g27_j[1]:.4f},{g27_j[2]:.4f}  vs X_1..3 = "
          f"{Xf[1]:.4f},{Xf[2]:.4f},{Xf[3]:.4f}  (should track)")
    print(f"   gamma_16(27) [relaxed] = {gam(16,27):.5f}   X_3={Xf[3]:.4f}  => "
          f"{'IN [2.3,2.9] PREDICTION HOLDS' if 2.3<=gam(16,27)<=2.9 else 'OUT OF PREDICTED BAND'}\n")

    # ---- (3) group by v_3 ----
    print("## (3) gamma_16(k) grouped by v_3(k)  (increases with v_3, tracks X_{v_3})")
    byv = {}
    for k in range(1, KMAX + 1):
        byv.setdefault(v3(k), []).append(k)
    for v in sorted(byv):
        ks = [k for k in byv[v] if k <= KMAX]
        vals = ", ".join(f"{k}:{gam(16,k):.3f}" for k in ks[:8])
        base = f"X_{v}={Xf[v]:.3f}" if v >= 1 else "gamma_1=0.667"
        print(f"   v_3={v} (base {base}, {'ENRICHED>1' if v>=1 else 'depleted<1'}): {vals}")
    print()

    # ---- (4) budget ----
    print("## (4) BUDGET: relaxation |log Pi_relax| / |threshold|, per channel  (binding = highest %)")
    print(f"   {'k':>3} {'v3':>3} {'base':>7} {'gamma_16':>9} {'Pi_relax':>9} {'log':>8} {'thresh':>8} {'budget%':>8}")
    rows = []
    for k in range(1, 13):
        v = v3(k)
        base = Xf[v] if v >= 1 else (2.0 / 3)
        g16 = gam(16, k)
        pi = g16 / base
        lp = log(pi)
        if v == 0:                                   # depleted: ceiling log(3/2)
            thr = log(1.5); budget = max(0.0, lp) / thr
        else:                                        # enriched: floor -log(X_v) (need gamma>1 => pi>1/X_v)
            thr = log(base); budget = max(0.0, -lp) / thr
        rows.append((k, budget))
        print(f"   {k:>3} {v:>3} {base:>7.4f} {g16:>9.5f} {pi:>9.5f} {lp:>+8.4f} {thr:>8.4f} {budget*100:>7.1f}%")
    rows.sort(key=lambda x: -x[1])
    print(f"   BINDING (top budget): " + ", ".join(f"k={k}({b*100:.0f}%)" for k, b in rows[:4]))
    print(f"   => proof reduces to: X_v>1 for v>=1 (CS, m=0 divergence, ALREADY PROVED) + bounded relaxation on the "
          f"binding channels only.")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
