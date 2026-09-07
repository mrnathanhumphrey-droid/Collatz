"""
p-ADIC-DICH-DEPTH (Wilson #2, 2026-09-06) -- do the FORWARD failures recover with exact depth, or stay dead?
For the failing primes, tabulate gamma_r(k) at the MULTIPLES of p across r, and read the trend of the core
claim (k=p enriched; fiber gamma(p)<gamma(2p)).  Reuses probe_p_adic_dich_sweep.build_chain_float.
Targets: p=11 -> r=6, p=13 -> r=5, p=23 -> r=4 (r=5 ~ 6.4M support ~ 1GB, memory-blocked; read the r=2..4 trend).
"""
import os, sys, time, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_p_adic_dich as P
from probe_p_adic_dich_sweep import build_chain_float


def gamma_at(rr, ks):
    N = len(rr)
    return {k: N * float(np.dot(rr, np.roll(rr, -k))) for k in ks}


def depth_study(p, q, RMAX, mults):
    print(f"\n## p={p}, q={q}: gamma_r(k) at multiples of p, r=2..{RMAX}   (mults {mults})")
    hdr = "   " + "r".rjust(3) + "".join(f"g({k})".rjust(11) for k in mults) + "   forward?"
    print(hdr)
    trend = {k: [] for k in mults}
    for r in range(2, RMAX + 1):
        nus = build_chain_float(p, q, r)
        rr, _ = P.rho_dense(nus[r], p, r)
        del nus
        g = gamma_at(rr, mults)
        del rr; gc.collect()
        for k in mults:
            trend[k].append(g[k])
        all_enr = all(g[k] > 1 for k in mults)
        fiber = (g[mults[0]] < g[mults[1]]) if len(mults) >= 2 else True
        fwd = all_enr and fiber
        print("   " + f"{r:>3}" + "".join(f"{g[k]:>11.4f}" for k in mults)
              + f"   {'OK' if fwd else 'FAIL'} (k=p {'enr' if g[mults[0]]>1 else 'DEPL'}, fiber {'ok' if fiber else 'no'})")
    # trend read on the core pair
    gp = trend[mults[0]]; g2p = trend[mults[1]] if len(mults) >= 2 else None
    def direction(seq):
        if len(seq) < 2:
            return "n/a"
        d = seq[-1] - seq[-2]
        return f"{'rising' if d>0 else 'falling'} (last step {d:+.4f})"
    print(f"   TREND g({mults[0]}=p): {[f'{x:.3f}' for x in gp]}  -> {direction(gp)} ; "
          f"crosses 1? {'below all r' if all(x<1 for x in gp) else ('above all r' if all(x>1 for x in gp) else 'CROSSES')}")
    if g2p:
        print(f"   TREND g({mults[1]}=2p): {[f'{x:.3f}' for x in g2p]}  -> {direction(g2p)} ; "
              f"{'below all r' if all(x<1 for x in g2p) else ('above all r' if all(x>1 for x in g2p) else 'CROSSES')}")


def main():
    t0 = time.time()
    print("# p-ADIC-DICH-DEPTH -- do FORWARD failures recover with depth? (q=2)")
    depth_study(11, 2, 6, [11, 22, 33, 44, 55])
    print(f"   ({time.time()-t0:.1f}s)")
    depth_study(13, 2, 5, [13, 26, 39, 52])
    print(f"   ({time.time()-t0:.1f}s)")
    depth_study(23, 2, 4, [23, 46])   # r=5 memory-blocked (~6.4M support); read r=2..4 trend
    print(f"   ({time.time()-t0:.1f}s)")
    print(f"\n# DONE ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
