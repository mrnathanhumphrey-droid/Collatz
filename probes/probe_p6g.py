"""
PROBE P6G (Wilson) -- the three-point kernel: gate K = cos2theta - 1/4 and the telescoped Lambda_i (2026-07-26).

Wilson assembled K in Fourier (c := cos theta):
  flanking 1/2[shift+shift^-1] = c ; full C_nu-hat = (5/4 + c)R_e-hat + B-hat ; deconv (5-4c) ; channel weight
  Re w -> (15/2)/(25-16c^2) - 1/2. Compose: K-hat = 1/4(25-16c^2)[(15/2)/(25-16c^2) - 1/2] = 2c^2 - 5/4 = cos2theta - 1/4.
  => K supported on 3 lags: 1/2 at n=+-2, -1/4 at n=0.  <C_rho, Re w> = R_e(2) - 1/4 R_e(0) + boundary.
  Per-level pairing = TWO numbers: ratio-4 autocorrelation R_e(2) and self-collision R_e(0).
  Positivity: <C_rho,Re w> > 0  <=>  R_e(2) > 1/4 R_e(0)  (Chebyshev/covariance family = m=0 proof's family).

The shell is a level-difference: A_i(k) = 3^i C_i(k) - 3^{i-1}C_{i-1}(k), so
  Lambda_i = 3^i[R_e^(i)(2) - 1/4 R_e^(i)(0)] - 3^{i-1}[R_e^(i-1)(2) - 1/4 R_e^(i-1)(0)] + boundary.

GATE: R_e^(i)(2), R_e^(i)(0) DIRECT at i=2..6; form 3^i-weighted level-differences; check vs certified Lambda_i
(shellA). Miss (if any) is in the 3^i normalization or the boundary -- both bounded/explicit -- NOT the kernel.

Reuses probe_p6b.shellA (certified channel) + probe_p6d.build_base2 (certified nu_e) + probe_p1.build_level. No new transport.
"""
import os, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p1 import build_level
from probe_p6b import shellA
from probe_p6d import build_base2


def autocorr(f):
    F = np.fft.fft(f)
    return np.fft.ifft(F * np.conj(F)).real


def kernel_symbolic_check():
    """verify K-hat = 1/4(25-16c^2)[(15/2)/(25-16c^2)-1/2] == 2c^2-5/4 on a theta grid."""
    th = np.linspace(0, np.pi, 41)[1:-1]; c = np.cos(th)
    Rew = (15.0 / 2) / (25 - 16 * c ** 2) - 0.5
    Khat = 0.25 * (25 - 16 * c ** 2) * Rew         # = (5-4c)(5/4+c) * Rew ; (5-4c)(5/4+c)=1/4(25-16c^2)
    target = 2 * c ** 2 - 1.25
    return float(np.max(np.abs(Khat - target)))


def lambda_certified(i, Kmax=60):
    """Lambda_i = Sum_{k>=1} 4^-k A_i(k mod Nn), A_i = shellA (certified primitive shell at level i)."""
    L = build_level(i); Nn = L['Nn']; W = L['What']
    A = np.array([shellA(L, W, r) for r in range(Nn)])          # A_i(r), r=0..Nn-1 (periodic)
    # exact geometric weight per residue r: Sum_{k>=1, k==r mod Nn} 4^-k
    denom = 1 - Fr(1, 4 ** Nn)
    lam = Fr(0)
    wr = []
    for r in range(Nn):
        kmin = r if r >= 1 else Nn
        w = Fr(1, 4 ** kmin) / denom
        wr.append(float(w))
    lam = float(np.dot(A, np.array(wr)))
    return lam, A


def main():
    t0 = time.time()
    print("# PROBE P6G -- three-point kernel K=cos2theta-1/4 + telescoped Lambda_i\n")

    print(f"## kernel symbolic check: max|K-hat - (2c^2-5/4)| on theta-grid = {kernel_symbolic_check():.2e}")
    print("   => K = {n=+-2: 1/2, n=0: -1/4};  <C_rho,Re w> = R_e(2) - 1/4 R_e(0) + boundary\n")

    # ---------- certified Lambda_i (pin convention against Lambda_1 = -2/21) ----------
    print("## certified Lambda_i via shellA (pin: Lambda_1 should be -2/21 = %.6f)" % (-2 / 21))
    lam = {}
    for i in range(1, 7):
        lam[i], _ = lambda_certified(i)
        tag = "  <== should be -2/21" if i == 1 else ""
        print(f"   Lambda_{i} = {lam[i]:+.6f}{tag}")
    tail = sum(lam[i] for i in range(2, 7))
    print(f"   Sum_{{i=2..6}} Lambda_i = {tail:+.6f}   (target 7/15 <=> -1/210 = {-1/210:+.6f}; measured ~ -0.00125)\n")

    # ---------- R_e^(i)(2), R_e^(i)(0) and the telescoped reconstruction ----------
    print("## R_e^(i)(lag) from nu_e (base-2), and 3^i-weighted level-difference vs certified Lambda_i")
    P = {}   # per-level pairing R_e(2) - 1/4 R_e(0)
    Re2 = {}; Re0 = {}
    for i in range(1, 7):
        S = build_base2(i)
        nu_e = S['R_e']
        Re = autocorr(nu_e)
        Re0[i] = float(Re[0]); Re2[i] = float(Re[2])
        P[i] = Re2[i] - 0.25 * Re0[i]
    print(f"   {'i':>2} {'R_e(2)':>12} {'R_e(0)':>12} {'P_i=R_e(2)-R_e(0)/4':>22} {'3^i P_i':>12}")
    for i in range(1, 7):
        print(f"   {i:>2} {Re2[i]:>12.6f} {Re0[i]:>12.6f} {P[i]:>22.6f} {3**i * P[i]:>12.6f}")
    print()
    print(f"   {'i':>2} {'Lambda_i cert':>14} {'3^i P_i - 3^(i-1) P_(i-1)':>26} {'boundary = diff':>16}")
    for i in range(2, 7):
        recon = 3 ** i * P[i] - 3 ** (i - 1) * P[i - 1]
        print(f"   {i:>2} {lam[i]:>14.6f} {recon:>26.6f} {lam[i] - recon:>16.6f}")

    # also: the positivity condition per level
    print("\n## positivity per level: R_e(2) > 1/4 R_e(0) ?  (Chebyshev/covariance family)")
    for i in range(1, 7):
        print(f"   i={i}: R_e(2)={Re2[i]:+.6f}  1/4 R_e(0)={0.25*Re0[i]:+.6f}  "
              f"{'>' if Re2[i] > 0.25*Re0[i] else '<='}  (P_i={P[i]:+.6f})")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
