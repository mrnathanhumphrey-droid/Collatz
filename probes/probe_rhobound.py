"""
PROBE RHOBOUND -- gate Wilson's rho-hat recursion via its UNCONDITIONAL bound + argmax parity (2026-07-26).

Wilson's corrected recursion (the DLOG/multiplicative object, the one that binds the channels):
    rho-hat_n(a) = D~(a) * E_{X~nu_{n-1}}[chi_a(3X+1)],   D~(a) = 1/(2 - (-1)^a e(-a/(2*3^n))).
Since |E[chi_a(3X+1)]| <= 1, this gives an UNCONDITIONAL, checkable bound:
    |rho-hat_n(a)|^2 <= |D~(a)|^2 = 1/(5 - 4(-1)^a cos(pi a/3^n)).
Near a=0: EVEN a -> denom ~5-4=1 -> bound ~1 (NO control); ODD a -> denom ~5+4=9 -> bound ~1/9 (ninefold suppression).
=> the sup must live at a EVEN with a/3^n small (or its conjugate). PREDICTION: argmax_a|rho-hat_k(a)| is EVEN, a/3^k small.

DECISIVE GATES (data-only, fft(rho) only, no recursion implementation):
 (1) bound NEVER violated: max_a |rho-hat_k(a)|^2 / |D~(a)|^2 <= 1 (else the (-1)^a parity is mis-indexed vs codebase a).
 (2) argmax parity: is argmax even with a/3^k small? (conjugate N-a has opposite parity, a/3^k near 1.)
If either fails, Wilson's recursion is still mis-indexed and nothing built on it counts.

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
    print("# PROBE RHOBOUND -- Wilson's |rho-hat|^2 <= 1/(5-4(-1)^a cos(pi a/3^k)) + argmax parity\n")

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

    print("## (1) BOUND CHECK: max_a |rho-hat|^2 / |D~(a)|^2  (must be <= 1; >1 => (-1)^a parity mis-indexed)")
    print(f"   {'k':>2} {'max ratio':>10} {'at a':>10} {'a parity':>8} {'a/3^k':>8}  {'bound-viol?':>11}")
    args = {}
    for k in range(2, 17):
        N = 3 ** k
        a = np.arange(N)
        ps = np.abs(np.fft.fft(rho[k])) ** 2
        with np.errstate(divide='ignore'):
            Dbound = 1.0 / (5 - 4 * ((-1.0) ** a) * np.cos(np.pi * a / N))
        prim = (a % 3 != 0)
        ratio = np.where(prim, ps / Dbound, 0.0)
        j = int(np.argmax(ratio))
        argmax_ps = int(a[prim][np.argmax(ps[prim])])
        args[k] = argmax_ps
        viol = ratio[j] > 1.0 + 1e-9
        print(f"   {k:>2} {ratio[j]:>10.4f} {j:>10} {'even' if j%2==0 else 'odd':>8} {j/N:>8.4f}  "
              f"{'VIOLATED' if viol else 'ok':>11}")
    print()

    print("## (2) ARGMAX PARITY: argmax_a|rho-hat_k(a)| -- Wilson predicts EVEN with a/3^k small")
    print(f"   {'k':>2} {'a_max':>11} {'parity':>7} {'a/3^k':>8} | {'conj N-a':>11} {'parity':>7} {'(N-a)/3^k':>10}")
    even_small = 0; tot = 0
    for k in range(2, 17):
        N = 3 ** k; am = args[k]; conj = N - am
        # pick representative in (0, N/2) to define "small"
        small = am if am < N / 2 else conj
        small_par = 'even' if small % 2 == 0 else 'odd'
        tot += 1
        if small % 2 == 0:
            even_small += 1
        print(f"   {k:>2} {am:>11} {'even' if am%2==0 else 'odd':>7} {am/N:>8.4f} | "
              f"{conj:>11} {'even' if conj%2==0 else 'odd':>7} {conj/N:>10.4f}   small-side={small_par}")
    print(f"\n   small-side EVEN in {even_small}/{tot} levels "
          f"({'CONFIRMS Wilson (even<->small)' if even_small >= tot-1 else 'FLIPPED or mixed => parity convention mis-indexed'})")

    print("\n## (3) small-a suppression: EVEN a bound ~1 (sup lives), ODD a bound ~1/9 (suppressed) -- check |rho-hat|^2")
    for k in (12, 14, 16):
        N = 3 ** k
        ps = np.abs(np.fft.fft(rho[k])) ** 2
        # smallest few even and odd a
        ev = [a for a in range(2, 40) if a % 2 == 0 and a % 3 != 0][:4]
        od = [a for a in range(1, 40) if a % 2 == 1 and a % 3 != 0][:4]
        print(f"   k={k}: even-a |rho-hat|^2 = {[f'{ps[a]:.2e}' for a in ev]} ; "
              f"odd-a = {[f'{ps[a]:.2e}' for a in od]}  (even >> odd near 0?)")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
