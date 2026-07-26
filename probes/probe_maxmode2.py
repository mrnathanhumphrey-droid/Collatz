"""
PROBE MAXMODE2 -- the decisive saturating-vs-climbing check on the max Fourier mode + U^2 identity + argmax migration.
(Wilson's hard flag, 2026-07-26.) "Slowing, 4% margin" is the 0.89 configuration -- deparity+two-step before ANY lemma.

|mu-hat_k(a)|^2 = |fft(rho_k)(a)|^2 (power spectrum; ps[0]=1). Channel identity (Parseval over k):
    mean_k |gamma_k(m)-1|^2 = Sum_{a!=0} |mu-hat_k(a)|^4  = Gowers U^2 norm of rho_k (minus trivial mode).
=> aggregate (Sum|mu-hat|^2 = X-1 -> inf) is the WRONG norm; the native object is ell^4 = U^2 (additive combinatorics).

CHECKS:
 (1) DECISIVE: deparitied + two-step rate of P_max(k)=max_{3-nmid a}|mu-hat_k(a)|^2, k=2..16. SATURATING (flat trend)
     => a uniform bound may exist; CLIMBING (rate rising with k) => NO uniform bound, 4% margin is an artifact, lemma FALSE.
 (2) argmax MIGRATION: which a achieves the max at each k? fixed or moving (as fraction a/3^k)?
 (3) U^2 identity: mean_k|gamma_k(m)-1|^2 == Sum_{a!=0} ps[a]^2, machine precision; report U^2 value + trend.
 (4) R66 answer (b): is ~2^-k FITTED or DERIVED? -> R66 doc sec 9: "empirically follows (1/2)^(k-c)... not yet clean
     closed form" = FITTED. (avg 3^-(k-1) IS derived via Parseval+invariant; max is NOT.)

Reuses channelfam rho-build. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def main():
    t0 = time.time()
    print("# PROBE MAXMODE2 -- saturating-vs-climbing (DECISIVE) + U^2 identity + argmax migration\n")

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

    Pmax, amax, U2 = {}, {}, {}
    for k in range(2, 17):
        N = 3 ** k
        ps = np.abs(np.fft.fft(rho[k])) ** 2
        m = np.arange(N)
        prim = (m % 3 != 0)
        idx = np.where(prim)[0]
        j = idx[np.argmax(ps[idx])]
        Pmax[k] = float(ps[j]); amax[k] = int(j)
        U2[k] = float((ps[1:] ** 2).sum())          # Sum_{a!=0} |mu-hat|^4 = U^2 (minus trivial)

    # ---- (1) DECISIVE saturating check ----
    print("## (1) DECISIVE: max-mode rate, raw / two-step / deparitied  (CLIMBING => lemma false)")
    print(f"   {'k':>2} {'P_max':>11} {'raw pow-ratio':>13} {'2step amp':>10} {'deparity amp':>12}")
    dep = {}
    for k in range(2, 17):
        dep[k] = 0.5 * (Pmax[k] + Pmax[k + 1]) if (k + 1) in Pmax else None
    for k in range(4, 17):
        raw = Pmax[k] / Pmax[k - 1]
        two = (Pmax[k] / Pmax[k - 2]) ** 0.25              # per-level AMPLITUDE (power^(1/2) per level, ^(1/2) again)
        dpr = ((dep[k] / dep[k - 2]) ** 0.25) if (dep.get(k) and dep.get(k - 2)) else float('nan')
        print(f"   {k:>2} {Pmax[k]:>11.4e} {raw:>13.4f} {two:>10.4f} {dpr:>12.4f}")
    # trend: is the two-step amplitude rate rising over k=8..16?
    amps = [(Pmax[k] / Pmax[k - 2]) ** 0.25 for k in range(8, 17)]
    early = np.mean(amps[:3]); late = np.mean(amps[-3:])
    print(f"   two-step amp: early(k8-10) mean {early:.4f} ; late(k14-16) mean {late:.4f} ; "
          f"drift {late-early:+.4f} => {'CLIMBING (no uniform bound!)' if late-early > 0.02 else 'SATURATING/flat (bound may exist)'}")
    # linear fit of rate vs k
    ks = np.arange(8, 17); slope = np.polyfit(ks, amps, 1)[0]
    print(f"   linear slope of two-step amp vs k (k=8..16): {slope:+.5f}/level "
          f"({'rising' if slope > 0.005 else ('falling' if slope < -0.005 else 'flat')})\n")

    # ---- (2) argmax migration ----
    print("## (2) argmax migration: which primitive a achieves P_max, as fraction a/3^k")
    for k in range(3, 17, 2):
        print(f"   k={k:>2}: a_max = {amax[k]:>10}  a/3^k = {amax[k]/3**k:.5f}  (3-a/3^k = {1-amax[k]/3**k:.5f})")
    print()

    # ---- (3) U^2 identity ----
    print("## (3) U^2 identity: mean_k|gamma_k(m)-1|^2 == Sum_{a!=0}|mu-hat|^4  (= Gowers U^2 of rho)")
    for k in (12, 14, 16):
        N = 3 ** k
        ps = np.abs(np.fft.fft(rho[k])) ** 2
        C = np.fft.ifft(ps).real
        gam = 3.0 ** k * C
        lhs = float(((gam - 1) ** 2).mean())
        rhs = float((ps[1:] ** 2).sum())
        print(f"   k={k}: mean_k(gamma-1)^2 = {lhs:.6e} ; Sum_{{a!=0}}|mu-hat|^4 = {rhs:.6e} ; rel {abs(lhs-rhs)/rhs:.1e}")
    u2r = np.median([U2[k] / U2[k - 1] for k in range(13, 17)])
    print(f"   U^2 per-level ratio (late) = {u2r:.4f} (amplitude {u2r**0.5:.4f}); U^2(16) = {U2[16]:.4e}")
    print()

    print("## (4) R66 (b): is ~2^-k max decay FITTED or DERIVED?")
    print("   R66 sec 9 (on disk): 'Max over primitive a: empirically follows (1/2)^(k-c)... NOT yet clean closed form.'")
    print("   => the ~2^-k max law is FITTED/empirical. Only the AVERAGE 3^-(k-1) is DERIVED (Parseval + invariant S_inf).")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
