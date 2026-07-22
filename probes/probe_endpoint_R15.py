"""
PROBE R15 -- THE ENDPOINT SPLIT AND THE BULK CORRELATION. Reuses R7/R9/R10; exact where marked.
Corrected retarget (R14): Lambda_r = S_r*<delta_r, Re w>_bulk + Lambda_r^unif,
  Lambda_r^unif = S_r * [3^r/(4^{3^r}-1) - 3^{r-1}/(4^{3^{r-1}}-1)] / (2*3^{r-1})   (R10-C trace, MEASURE-FREE),
  bulk S_r*b_r := Lambda_r - Lambda_r^unif.
R15-A gate anchors (r=2): Lambda_2^unif=-110/29127, S_2*b_2=-240/67963, sum=-1490/203889=Lambda_2.
R15-B: is x in {0,1/2} a primitive angle? (3^r odd). R15-C: bulk ledger. R15-D: two-band split. R15-E: housekeeping.
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

S = dict(R10.S)
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
for k in (7, 8):
    S[k] = F(7, 15) + EPS[k]

Rew = lambda x: 15 / (2 * (17 - 8 * math.cos(2 * math.pi * x))) - 0.5


def Lunif_exact(r):
    """Lambda_r^unif from R10-C trace, measure-free exact."""
    tr = F(3 ** r, 4 ** (3 ** r) - 1) - F(3 ** (r - 1), 4 ** (3 ** (r - 1)) - 1)
    return S[r] * tr / (2 * 3 ** (r - 1))


def Lambda_exact(mu, r):
    if r <= 5:
        return R10.Lambda_r(mu, r)[0]
    return (EPS[r + 1] - EPS[r]) / 2          # certified R12-F (r=6,7)


def theta2(mu_r, r):
    N = 3 ** r
    g = R10.autocorr_dlog(mu_r, r); gf = [float(x) for x in g]
    return [sum(gf[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]


def main():
    print("# PROBE R15 -- ENDPOINT SPLIT + BULK CORRELATION. Exact gates + labeled measurement.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 6):
        mu[k] = R7.build_mu(mu[k - 1], k)

    Lam = {r: Lambda_exact(mu, r) for r in range(2, 8)}
    Lun = {r: Lunif_exact(r) for r in range(2, 8)}
    bulk = {r: Lam[r] - Lun[r] for r in range(2, 8)}

    # ---- R15-A ----
    print("## R15-A  ENDPOINT IDENTITY (forced, exact): Lambda_r = S_r*b_r + Lambda_r^unif")
    okA = True
    for r in range(2, 8):
        good = (bulk[r] + Lun[r] == Lam[r])
        print(f"   r={r}: Lambda_r^unif={str(Lun[r])[:24]:>24}  S_r*b_r={str(bulk[r])[:24]:>24}  "
              f"sum==Lambda_r? {bulk[r]+Lun[r]==Lam[r]}")
        okA = okA and good
    a_ok = (Lun[2] == F(-110, 29127) and bulk[2] == F(-240, 67963) and Lam[2] == F(-1490, 203889))
    print(f"   ANCHOR r=2: Lambda_2^unif={Lun[2]} (=-110/29127? {Lun[2]==F(-110,29127)}), "
          f"S_2*b_2={bulk[2]} (=-240/67963? {bulk[2]==F(-240,67963)}), Lambda_2={Lam[2]} (=-1490/203889? {Lam[2]==F(-1490,203889)})")
    print(f"   => R15-A {'GATE PASS' if okA and a_ok else 'FAIL (#34)'}\n")

    # ---- R15-B ----
    print("## R15-B  FIXED-POINT CONCENTRATION (forced): which of {0,1/2} is a primitive angle? (3^r odd)")
    for r in range(2, 6):
        N = 3 ** r
        prim = [k for k in range(1, N) if k % 3 != 0]
        has0 = (0 in prim)                       # x=0 <-> k=0
        half_k = N / 2
        has_half = (half_k == int(half_k)) and (int(half_k) in prim)
        # Sum_{k prim} Re w = R10-C trace (real, exact); certify Lambda^unif uses it
        num_sumRew = sum(Rew(k / N) for k in prim)
        trace = float(F(3 ** r, 4 ** (3 ** r) - 1) - F(3 ** (r - 1), 4 ** (3 ** (r - 1)) - 1))
        print(f"   r={r} (N={N}): x=0 primitive? {has0}  x=1/2 primitive? {has_half} (k=N/2={half_k} "
              f"{'non-integer' if half_k!=int(half_k) else 'integer'})  | "
              f"Sum_prim Re w={num_sumRew:.3e} vs R10-C trace={trace:.3e} [{'MATCH' if abs(num_sumRew-trace)<1e-9 else 'DEV'}]")
    print("   => NEITHER x=0 (trivial char, non-primitive) NOR x=1/2 (non-lattice, 3^r odd) is a primitive angle.")
    print("      No self-conjugate primitive angle => Lambda^unif is the primitive-SAMPLING RESIDUAL (=R10-C trace),")
    print("      NOT a literal endpoint atom; the 'both ends' picture needs the near-{0,1/2} NEIGHBORHOODS.\n")

    # ---- R15-C ----
    print("## R15-C  BULK CORRELATION LEDGER (measurement): S_r*b_r = Lambda_r - Lambda_r^unif")
    # Sum_{r>=1} Lambda_r = -1/10; Sum S_r b_r = -1/10 - Sum Lambda_r^unif  (bulk_1 = 0)
    Lun1 = Lunif_exact(1)
    sumLun = Lun1 + sum(Lun[r] for r in range(2, 8))
    print(f"   Lambda_1^unif={Lun1} (=Lambda_1=-2/21? {Lun1==F(-2,21)} => bulk_1=0, layer 1 all-uniform)")
    print(f"   {'r':>2} {'S_r*b_r (=bulk)':>22} {'float':>13} {'running Sum S_r b_r':>20}")
    run = F(0)
    for r in range(2, 8):
        run += bulk[r]
        print(f"   {r:>2} {str(bulk[r])[:22]:>22} {float(bulk[r]):>+13.6e} {float(run):>+20.8f}")
    tail_target = F(-1, 10) - sumLun    # Sum_{r>=2} S_r b_r target (bulk_1=0)
    print(f"   target Sum_(r>=2) S_r b_r = -1/10 - Sum_(r>=1)Lambda^unif = {float(tail_target):.8f} "
          f"(Sum Lambda^unif={float(sumLun):.3e}, ~all at r=1,2)")
    print(f"   signs S_r*b_r: {['+' if bulk[r]>0 else '-' for r in range(2,8)]}\n")

    # ---- R15-D ----
    print("## R15-D  TWO-BAND DECOMPOSITION (measurement, NO fit): b_r split at Re w=0 (x=arccos(1/4)/2pi=0.2088)")
    print(f"   {'r':>2} {'a_r (Re w>0 band)':>18} {'c_r (Re w<0 band)':>18} {'b_r=a_r-c_r':>14} {'a_r/c_r':>10} {'vs 5/3':>8}")
    for r in range(2, 7):
        N = 3 ** r
        th2 = theta2(mu[r], r) if r <= 5 else None
        if r == 6:
            mu[6] = R7.build_mu(mu[5], 6); th2 = theta2(mu[6], 6)
        prim = [k for k in range(1, N) if k % 3 != 0]; M = len(prim); Sr = float(S[r])
        delta = {k: th2[k] / Sr - 1.0 / M for k in prim}
        a = sum(Sr * delta[k] * Rew(k / N) for k in prim if Rew(k / N) > 0)
        c = -sum(Sr * delta[k] * Rew(k / N) for k in prim if Rew(k / N) < 0)
        br = a - c
        print(f"   {r:>2} {a:>+18.6e} {c:>+18.6e} {br:>+14.6e} {a/c if c else float('nan'):>10.4f} "
              f"{'>' if (a/c if c else 0)>5/3 else '<':>8}")
    print("   (does a_r/c_r cross 5/3 where b_r flips sign; monotone or oscillate? table only, no law fitted.)\n")

    # ---- R15-E ----
    print("## R15-E  CONVENTION + FEASIBILITY (housekeeping)")
    L2_from_mu = R10.Lambda_r(mu, 2)[0]
    print(f"   (i) r=2 from mu_2 table: Lambda_2={L2_from_mu} (==-1490/203889? {L2_from_mu==F(-1490,203889)}); "
          f"Lambda_2^unif(trace)={Lun[2]} -> bulk={L2_from_mu-Lun[2]} (==-240/67963? {L2_from_mu-Lun[2]==F(-240,67963)}) [2nd gate OK]")
    print(f"   (ii) r=8 wall: exact b_8 needs Lambda_8=(eps_9-eps_8)/2 but eps_9 is FLOAT-only (k>=9 wall); "
          f"D-band(r=8) needs mu_8 profile (supp 4374, ~19M-Fraction autocorr). Bulk-sequence exact ceiling = r=7. DEFER.")


if __name__ == "__main__":
    main()
