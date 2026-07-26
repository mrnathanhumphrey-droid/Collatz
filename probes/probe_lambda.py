"""
PROBE LAMBDA -- the aligned table + the power-saving delta against the 2^m ceiling (2026-07-26).

Wilson's ceiling: at i=k, 3^k=0 => pi-hat=1 => G(k,j)=2^j; phases=1 preserves it downward => |G(0,m)| <= 2^m
UNCONDITIONALLY. So sup|pi-hat| = 2^-m|G(0,m)| <= 1 (trivial), and the decay IS the saving of G below its ceiling.

Per-UNIT-m the ceiling growth is exactly 2. Define lambda_m = per-unit-m growth of |G(0,m)| = (|G|_k/|G|_{k'})^{1/dm}.
Saving delta = 1 - lambda_m/2. Qualitative delta>0 <=> geometric decay of sup|pi-hat|.

CONVENTIONS (explicit, since they crossed this session):
  sup|pi-hat|   = |pi-hat(xi*)|      (LINEAR modulus)          -- e.g. 0.25224 at k=3
  R66 = sup|pi-hat|^2                (SQUARED)                 -- e.g. 0.06362 at k=3
  |G(0,m)| = 2^m * sup|pi-hat|       (linear; ceiling 2^m)

One row per k. Fit lambda_m; check convergence/drift. Reuses probe_singlerec.fwd_hat. fwd_hat(k) cost grows ~3^k.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_singlerec import fwd_hat


def v2(n):
    n = int(n); v = 0
    while n and n % 2 == 0:
        n //= 2; v += 1
    return v


def main():
    t0 = time.time()
    print("# PROBE LAMBDA -- aligned table + power-saving delta vs the 2^m ceiling\n")
    print("## conventions: sup=|pi-hat(xi*)| LINEAR ; R66=sup^2 ; |G(0,m)|=2^m*sup (ceiling 2^m)\n")

    KMAX = 15
    print(f"   {'k':>2} {'xi*':>10} {'m':>3} {'dm':>3} {'sup=|pi|':>10} {'R66=|pi|^2':>11} "
          f"{'|G(0,m)|':>11} {'|G|/2^m':>8} {'Srate':>7} {'lam_m':>7} {'delta':>7}")
    prev = None
    lam_list = []
    for k in range(3, KMAX + 1):
        N = 3 ** k
        try:
            ph, _ = fwd_hat(k)
        except MemoryError:
            print(f"   k={k}: MemoryError -- stop.")
            break
        idx = np.arange(N)
        absph = np.abs(ph)
        absph[idx % 3 == 0] = -1.0        # exclude 3|xi (and DC)
        xistar = int(np.argmax(absph)); sup = float(absph[xistar])
        even = xistar if xistar % 2 == 0 else N - xistar
        m = v2(even)
        ispow = (even == 2 ** m)
        G = (2.0 ** m) * sup
        r66 = sup * sup
        if prev is None:
            dm = ''; srate = ''; lam = ''; delta = ''
        else:
            pk, pm, psup, pG = prev
            d = m - pm
            dm = f"{d:+d}"
            srate = f"{sup/psup:.4f}"
            lam_m = (G / pG) ** (1.0 / d) if d != 0 else float('nan')
            lam = f"{lam_m:.4f}"
            delta = f"{1 - lam_m/2:.4f}"
            lam_list.append((k, lam_m))
        tag = '' if ispow else ' !NOTPOW'
        print(f"   {k:>2} {xistar:>10} {m:>3} {dm:>3} {sup:>10.6f} {r66:>11.6f} "
              f"{G:>11.4f} {sup:>8.4f} {srate:>7} {lam:>7} {delta:>7}{tag}")
        prev = (k, m, sup, G)

    print("\n   [|G|/2^m == sup exactly (identity check). lam_m = per-unit-m growth of |G|; ceiling 2; delta=1-lam_m/2.]")
    if lam_list:
        print("\n## lambda_m drift (does the saving converge?)")
        print("   " + "  ".join(f"k{k}:{lm:.3f}" for k, lm in lam_list))
        # tail estimate over last few
        tail = [lm for _, lm in lam_list[-4:]]
        print(f"   tail mean lambda_m (last {len(tail)}) = {np.mean(tail):.4f}  =>  delta = {1-np.mean(tail)/2:.4f} "
              f"=>  per-unit-m sup-saving = {np.mean(tail)/2:.4f}")
        print(f"   [ceiling per k = 2^(dm/dk); if lambda_m -> const < 2, delta>0 is a genuine power saving.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
