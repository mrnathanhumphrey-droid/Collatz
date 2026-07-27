"""
PROBE P6J -- ENRICHED: the same-parity hemisphere (3|k channels, M+ = 5/3) (2026-07-26).

Our T_i, Lambda_i sum ALL k channels (= 1/2 S). Split by parity of the base-4 lag k:
  T^-_i (cross / depleted, 3-nmid k, M- = 2/3) -- where the cascade/d1 lives
  T^+_i (same  / enriched,  3|k,     M+ = 5/3)
via gamma_i(k) = 3^i <rho_i, shift_k rho_i> (certified numerator-profile channel).

Cheap: build the base-4 numerator profile directly (no dense bridge matrix), gate vs build_level's rho.

E-A gate: (1/3)M+ + (2/3)M- = 1 every level (avg gamma over all k == 1 exactly), and 3|k shells reproduce
          certified gamma_inf(3)=1.2372, gamma_inf(6)=1.3717, gamma_inf(9)=2.112.
E-B: T^+_i, Lambda^+_i i=1..15; sign of Lambda^+; deparitied rate alongside ours.
E-C: is T^+ nearer its bound than T is to 7/30?
E-D: Lambda^- + Lambda^+ == full Lambda (= 1/2 dS) -- free redundancy check.

Reuses probe_27.stationary_trunc + probe_p1.build_level (gate only, low i). No new transport.
"""
import os, sys, time
from math import gcd
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_27_high_k_rho_q5 import stationary_trunc
from probe_p1 import build_level


def build_rho4(n):
    """certified base-4 numerator profile rho_i (= build_level's rho), no dense matrix."""
    q = 3 ** (n + 1); Nn = 3 ** n
    DL = np.full(q, -1, dtype=np.int64)
    g = 1
    for s in range(Nn):
        DL[g] = s; g = (g * 4) % q
    piW, _ = stationary_trunc(3, n)
    r = np.arange(Nn); cp = r[r % 3 != 0]
    nu = np.asarray(piW, float); nu = nu / nu.sum()
    Y = (3 * cp + 1) % q                      # numerator, always ==1 mod 3 (coset-1 = <4>), no fold
    return np.bincount(DL[Y], weights=nu, minlength=Nn)


def autocorr(f):
    F = np.fft.fft(f); return np.fft.ifft(F * np.conj(F)).real


def weighted_sum(C, Nn, which, KMAX=64):
    """Sum_{k>=1, which(k)} 4^-k C(k mod Nn); 4^-k negligible past k~30, so sum k=1..KMAX (C periodic mod Nn)."""
    s = 0.0
    for k in range(1, KMAX + 1):
        if which(k):
            s += (4.0 ** -k) * C[k % Nn]
    return s


def main():
    t0 = time.time()
    print("# PROBE P6J -- ENRICHED hemisphere (3|k, M+=5/3)\n")

    # ---- GATE build_rho4 == build_level rho ----
    print("## GATE build_rho4 == build_level['rho'] (base-4 numerator profile, i=2..6)")
    for i in range(2, 7):
        r4 = build_rho4(i); L = build_level(i)
        print(f"   i={i}: max|build_rho4 - build_level rho| = {np.max(np.abs(r4 - L['rho'])):.1e}")
    print()

    IMAX = 15
    C = {}; Nn = {}
    for i in range(1, IMAX + 1):
        rho = build_rho4(i); Nn[i] = 3 ** i
        C[i] = autocorr(rho)

    # ---- E-A: pinned mean + certified channel constants ----
    print("## E-A: pinned mean (1/3)M+ + (2/3)M- = 1 ; gamma_inf(3,6,9) = 1.2372, 1.3717, 2.112")
    for i in (3, 6, 9, 12, 15):
        g = 3 ** i * C[i]                                  # gamma_i(k) = 3^i C(k)
        k3 = np.arange(Nn[i]) % 3 == 0
        Mp = g[k3].mean(); Mm = g[~k3].mean()
        pinned = Mp / 3 + 2 * Mm / 3
        print(f"   i={i}: M+={Mp:.5f} M-={Mm:.5f}  (1/3)M+ +(2/3)M- = {pinned:.8f}  "
              f"| gamma_i(3)={g[3]:.4f} gamma_i(6)={g[6]:.4f} gamma_i(9)={g[9]:.4f}")
    print("   [gamma_i(k) -> gamma_inf(k) as i grows; check k=3->1.2372, k=6->1.3717, k=9->2.112]\n")

    # ---- E-B / E-D: T^+, T^-, Lambda^+, Lambda^-, and their sum ----
    Tp = {0: 0.0}; Tm = {0: 1.0 / 3}          # T^-_0 = 1/3 anchor (full T_0); split T^+_0=0 (no 3|k below k=Nn)
    for i in range(1, IMAX + 1):
        base = 3 ** i
        Tp[i] = base * weighted_sum(C[i], Nn[i], lambda k: k % 3 == 0)
        Tm[i] = base * weighted_sum(C[i], Nn[i], lambda k: k % 3 != 0)
    print("## E-B/E-D: T^+ (enriched 3|k), T^- (depleted 3-nmid k), Lambda each, sum vs full")
    print(f"   {'i':>2} {'T^+':>11} {'T^-':>11} {'T^++T^-':>11} {'Lam^+':>12} {'Lam^-':>12}")
    Lp = {}; Lm = {}
    for i in range(1, IMAX + 1):
        Lp[i] = Tp[i] - Tp[i - 1]; Lm[i] = Tm[i] - Tm[i - 1]
        print(f"   {i:>2} {Tp[i]:>11.7f} {Tm[i]:>11.7f} {Tp[i]+Tm[i]:>11.7f} {Lp[i]:>+12.8f} {Lm[i]:>+12.8f}")
    print(f"   [T^++T^- should = full T_i=1/2 S_(i+1); T_15 full ~ 0.235676]\n")

    # ---- E-B: deparitied rates of Lambda^+ and Lambda^- ----
    print("## E-B: deparitied two-step rate (Lam_i/Lam_(i-2))^(1/2), enriched vs depleted")
    print(f"   {'i':>2} {'Lam^+ rate':>12} {'Lam^- rate':>12}")
    for i in range(8, IMAX + 1):
        rp = (Lp[i] / Lp[i - 2]) ** 0.5 if Lp[i] * Lp[i - 2] > 0 else float('nan')
        rm = (Lm[i] / Lm[i - 2]) ** 0.5 if Lm[i] * Lm[i - 2] > 0 else float('nan')
        print(f"   {i:>2} {rp:>12.5f} {rm:>12.5f}")

    # ---- E-C: T^+ vs its limit/bound; T^- vs 7/30 ----
    print("\n## E-C: T^+ trajectory vs T^- (ours). T^-->7/30=0.23333? T^+->?")
    print(f"   T^+_15={Tp[15]:.6f} (Lam^+_15={Lp[15]:+.2e})   T^-_15={Tm[15]:.6f} (Lam^-_15={Lm[15]:+.2e}, vs 7/30={7/30:.6f})")
    print(f"   T^+ sign of Lambda: {'POSITIVE' if all(Lp[i]>0 for i in range(3,IMAX+1)) else 'MIXED/NEG'}; "
          f"T^- sign: {'POSITIVE' if all(Lm[i]>0 for i in range(8,IMAX+1)) else 'MIXED/NEG'}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
