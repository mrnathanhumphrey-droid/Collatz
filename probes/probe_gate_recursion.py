"""
PROBE GATE_RECURSION -- gate Wilson's mu-hat recursion against R66 BEFORE building the (delta,A) escape scan (2026-07-26).

Wilson's recursion (from Syrac_n = 2^-a (3 Syrac_{n-1}+1), a~Geom(1/2)):
    mu-hat_n(xi) = Sum_{a>=1} 2^-a e(xi 2^-a / 3^n) mu-hat_{n-1}(xi 2^-a mod 3^{n-1}).
Wilson FLAGGED: he derived this from standard Syracuse structure, NOT R66's code. Gate: reproduce R66's mu-hat_k (k<=6).
Index conventions to test: a>=1 vs a>=0; 2^-a = modular inverse; arg reduced mod 3^{n-1}; phase mod 3^n.

CRITICAL 2nd question (mine): the CHANNEL-binding object (MAXMODE2, verified via U^2) is fft(rho) = the DLOG profile
(multiplicative Fourier of nu). R66's mu-hat and Wilson's recursion are the RAW/additive Fourier fft(nu). If
fft(nu) != fft(rho), the escape lemma bounds the WRONG object for the channels. Compare max|fft(nu)|^2 vs max|fft(rho)|^2
vs R66's table (k2 max=0.1428, k3=0.0636, k4=0.0313).

Reuses build_nu + R10.dlog_table. build_nu(0.5,6). Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

R66_MAX = {2: 0.1428, 3: 0.0636, 4: 0.0313, 5: 0.0167, 6: 0.00924}
R66_K2 = {1: 0.04935, 2: 0.04592, 4: 0.14283, 5: 0.14283, 7: 0.04592, 8: 0.04935}


def dense_nu(nu_k, N):
    a = np.zeros(N)
    for X, w in nu_k.items():
        a[X % N] += float(w)
    return a / a.sum()


def dense_rho(nu_k, N, r):
    d = R10.dlog_table(r)
    mu = np.zeros(N)
    for X, w in nu_k.items():
        mu[(X - 1) // 3 % N] += float(w)
    rr = np.zeros(N); rr[np.fromiter((d[a] for a in range(N)), np.int64, N)] = mu
    return rr / rr.sum()


def main():
    t0 = time.time()
    print("# PROBE GATE_RECURSION -- validate Wilson's recursion vs R66; identify the object\n")
    nus = build_nu(0.5, 6)

    nu_d = {k: dense_nu(nus[k], 3 ** k) for k in range(1, 7)}
    rho_d = {k: dense_rho(nus[k], 3 ** k, k) for k in range(1, 7)}
    mnu = {k: np.fft.fft(nu_d[k]) for k in range(1, 7)}      # additive Fourier of raw nu
    mrho = {k: np.fft.fft(rho_d[k]) for k in range(1, 7)}    # channel object (dlog / multiplicative)

    # ---- WHICH object is R66's mu-hat? ----
    print("## (0) which Fourier object is R66's mu-hat?  (R66 k=2 max=0.1428, k=3=0.0636, k=4=0.0313)")
    print(f"   {'k':>2} {'max|fft(nu)|^2':>15} {'max|fft(rho)|^2':>16} {'R66 max':>9}")
    for k in range(2, 7):
        N = 3 ** k; m = np.arange(N); prim = m % 3 != 0
        vnu = (np.abs(mnu[k]) ** 2)[prim].max()
        vrho = (np.abs(mrho[k]) ** 2)[prim].max()
        print(f"   {k:>2} {vnu:>15.5f} {vrho:>16.5f} {R66_MAX[k]:>9.5f}  "
              f"-> R66 = {'fft(nu) RAW' if abs(vnu-R66_MAX[k])<abs(vrho-R66_MAX[k]) else 'fft(rho) DLOG'}")
    print("   k=2 per-a |fft(nu)|^2 vs R66 (a=1,2,4):", end=" ")
    for a in (1, 2, 4):
        print(f"a{a}:{abs(mnu[2][a])**2:.5f}(R66 {R66_K2[a]:.5f})", end="  ")
    print("\n   k=2 per-a |fft(rho)|^2 (a=1,2,4):", end=" ")
    for a in (1, 2, 4):
        print(f"a{a}:{abs(mrho[2][a])**2:.5f}", end="  ")
    print("\n")

    # ---- test Wilson's recursion on BOTH objects, several conventions ----
    def test_recursion(mh, label, amax=64):
        print(f"## recursion test on {label}:  mu_n(xi) =?= Sum_a 2^-a e(xi 2^-a/3^n) mu_{{n-1}}(xi 2^-a mod 3^{{n-1}})")
        for n in range(3, 7):
            N = 3 ** n; Nm = 3 ** (n - 1)
            inv2 = pow(2, -1, N)
            prim = [xi for xi in range(1, N) if xi % 3 != 0]
            # variant weights: raw 2^-a  and normalized 2^-a/(1-2^-Mm)
            errs = []
            for xi in prim[:2000]:
                rhs = 0j; p = 1
                for a in range(1, amax + 1):
                    p = (p * inv2) % N
                    wa = (xi * p) % N
                    rhs += (0.5 ** a) * np.exp(2j * np.pi * wa / N) * mh[n - 1][wa % Nm]
                lhs = mh[n][xi]
                errs.append(abs(rhs - lhs))
            errs = np.array(errs)
            rel = errs.mean() / (np.abs([mh[n][xi] for xi in prim[:2000]]).mean())
            print(f"   n={n}: mean|rhs-lhs| = {errs.mean():.3e}  rel = {rel:.2e}  "
                  f"[{'REPRODUCES' if rel < 1e-3 else ('close' if rel<0.05 else 'NO')}]")
        print()

    test_recursion(mnu, "fft(nu) [RAW/additive]")
    test_recursion(mrho, "fft(rho) [DLOG/channel object]")

    print("## VERDICT: (see above) which object R66/Wilson's recursion is, and whether it = the channel object fft(rho).")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
