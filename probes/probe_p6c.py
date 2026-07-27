"""
PROBE P6C (Wilson) -- the support-separation check at levels 2,3 (2026-07-26).

Wilson's mechanism candidate: a mod-3 support separation. If supp(nu_e) subset {s == alpha mod 3} and
supp(nu_o) subset {s == beta mod 3} with alpha != beta, then same-parity correlations need s and s+m in the
same class (=> 3|m) and cross-parity needs m == beta-alpha != 0 (3-nmid m) -- EXACTLY the P6B table.

Wilson checked LEVEL 1 by hand and it FAILS: x==1 mod3 => 3x+1==4 mod9 always, so x' = 4*2^-a mod 9; even a
runs x' over {1,7,4} = base-4 dlogs {0,2,1} = ALL THREE classes. No naive mod-3 separation at level 1.

QUESTION (this probe): does a separation appear at a FINER modulus / different grading at levels 2,3(,4,5)?
Grade the CERTIFIED, GATED profiles rho_e, rho_o (P6B: rho_from o partial_pihat, rho_e+rho_o=rho_full) by
base-4 dlog s -- Wilson's exact coordinate -- mod 3 and mod 9. Report signed mass and L1 mass per class, and the
concentration (max class share). If one class carries ~all the (parity-specific) mass at some level => mechanism
settled and the collapse becomes a corollary. If it stays spread => the separation is not in this grading.

Zero new transport: reuses probe_p1.build_level/bridge + probe_p6.partial_pihat + P6B rho_from. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p1 import build_level
from probe_p6 import partial_pihat
from probe_p6b import rho_from


def classmass(rho, s, mod):
    """signed sum and L1 mass of rho grouped by (s % mod)."""
    r = s % mod
    sig = np.array([rho[r == c].sum() for c in range(mod)])
    l1 = np.array([np.abs(rho[r == c]).sum() for c in range(mod)])
    return sig, l1


def main():
    t0 = time.time()
    print("# PROBE P6C -- support-separation check (base-4 dlog s), levels 2..5\n")
    print("Wilson: mod-3 separation of supp(nu_e)/supp(nu_o) would GIVE the P6B table (3-nmid m => cross).")
    print("Level-1 negative (by hand): even-a support = {s=0,1,2} = all classes. Does a finer modulus separate?\n")

    for j in (2, 3, 4, 5):
        L = build_level(j); N = 3 ** j; W = L['What']; s = np.arange(N)
        pe = partial_pihat(W, N, lambda a: a % 2 == 0)
        po = partial_pihat(W, N, lambda a: a % 2 == 1)
        rho_e = rho_from(L, pe); rho_o = rho_from(L, po); rho_full = rho_from(L, W)
        gate = np.max(np.abs(rho_e + rho_o - rho_full))
        print(f"## j={j}  (N=3^{j}={N})   [gate rho_e+rho_o==rho_full: {gate:.1e}]")

        for mod in (3, 9):
            print(f"   -- grading: s mod {mod} --")
            for name, rho in (("rho_e", rho_e), ("rho_o", rho_o)):
                sig, l1 = classmass(rho, s, mod)
                tot = l1.sum() + 1e-30
                share = l1 / tot                       # fraction of L1 mass in each class
                conc = share.max()                     # concentration: 1 => single-class support
                sig_str = " ".join(f"{v:+.4f}" for v in sig)
                shr_str = " ".join(f"{v:.3f}" for v in share)
                print(f"      {name}: signed[{sig_str}]  L1-share[{shr_str}]  maxshare={conc:.3f}"
                      + ("  <== SEPARATED" if conc > 0.90 else ""))
        print()

    # Direct forward-measure cross-check at level 1 (Wilson's hand computation, reproduced from machinery)
    print("## level-1 sanity (base-4 dlog of even-branch images): expect {0,1,2} spanned")
    j = 1; L = build_level(j); N = 3 ** j; W = L['What']
    pe = partial_pihat(W, N, lambda a: a % 2 == 0)
    rho_e = rho_from(L, pe)
    supp = np.where(np.abs(rho_e) > 1e-9)[0]
    print(f"   j=1 rho_e nonzero at s in {sorted(supp.tolist())}  (level-1: all of {{0,1,2}} if no separation)")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
