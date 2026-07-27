"""
PROBE P6H (Wilson) -- the telescope collapse: S_{i+1} = 2*3^i[4 R_e(2) - R_e(0)], one limit (2026-07-26).

Wilson's telescope: T_i := 3^i[4 R_e^(i)(2) - R_e^(i)(0)] = Sum_k 4^-k gamma_i(k) = 1/2 S^(i) = 1/2 S_{i+1}.
So Lambda_i = T_i - T_{i-1} telescopes; the WHOLE S-ladder is two autocorrelation values of nu_e:
    S_{i+1} = 2 * 3^i [4 R_e^(i)(2) - R_e^(i)(0)] ;  S_inf = 2 lim T_i ; 7/15 <=> lim T_i = 7/30.
Anchors: T_0 = 1/3 = S_1/2, T_1 = 5/21 = S_2/2.
Asymptotics: 3^i R_e(0) and 3^i R_e(2) are COLLISION probabilities (like X_i=3^i Sum nu^2, LINEARLY divergent);
T_i is the RESIDUE of their cancellation -> S_inf/2. 4*slope(R2)==slope(R0) would be the cancellation identity.

Identity gate (0) uses shellA (dense matrix) -> LOW i only (i<=6). Asymptotics (1)(2)(3) use build_base2 only -> high i.
Reuses probe_p6d.build_base2 + probe_p6b.shellA + probe_p1.build_level. No new transport.
"""
import os, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p6d import build_base2
from probe_p6b import shellA
from probe_p1 import build_level


def autocorr(f):
    F = np.fft.fft(f)
    return np.fft.ifft(F * np.conj(F)).real


def lambda_cert(i):
    L = build_level(i); Nn = L['Nn']; W = L['What']
    A = np.array([shellA(L, W, r) for r in range(Nn)])
    denom = 1 - Fr(1, 4 ** Nn)
    w = np.array([float(Fr(1, 4 ** (k if k >= 1 else Nn)) / denom) for k in range(Nn)])
    return float(np.dot(A, w))


def main():
    t0 = time.time()
    IMAX = 12
    R0 = {}; R2 = {}; Xfull = {}
    for i in range(1, IMAX + 1):
        S = build_base2(i)
        nu_e = S['R_e']; nu = S['R_e'] + S['R_o']
        Re = autocorr(nu_e)
        R0[i] = float(Re[0]); R2[i] = float(Re[2])
        Xfull[i] = float(3 ** i * np.sum(nu ** 2))
    print(f"# PROBE P6H -- telescope collapse, asymptotics (i=1..{IMAX})  [build {time.time()-t0:.1f}s]\n")

    # ----- (0) identity gate: S_{i+1} = 2 T_i^Re  vs certified S-ladder (i<=6) -----
    print("## (0) identity:  S_{i+1} = 2*3^i[4R_e(2)-R_e(0)]  vs  certified 2/3 + 2 Sum_{j<=i} Lambda_j  (i<=6)")
    lam = {i: lambda_cert(i) for i in range(1, 7)}
    print(f"   {'i':>2} {'T_i':>12} {'S_(i+1)=2T_i':>13} {'S cert(shellA)':>15} {'diff':>10}")
    for i in range(1, 7):
        Ti = 3 ** i * (4 * R2[i] - R0[i]); Sre = 2 * Ti
        Scert = 2.0 / 3 + 2 * sum(lam[j] for j in range(1, i + 1))
        tag = "  = 10/21" if i == 1 else ""
        print(f"   {i:>2} {Ti:>12.8f} {Sre:>13.8f} {Scert:>15.8f} {Sre-Scert:>10.1e}{tag}")
    print(f"   [S_1=2/3={2/3:.6f}, T_1=5/21={5/21:.6f}; S_2=10/21={10/21:.6f}]\n")

    # ----- (1)(2) the two divergent sequences + slope cancellation -----
    print("## (1)(2) 3^i R_e(0), 3^i R_e(2) diverge linearly; is 4*slope(R2)==slope(R0)? (cancellation)")
    A0 = {i: 3 ** i * R0[i] for i in range(1, IMAX + 1)}
    A2 = {i: 3 ** i * R2[i] for i in range(1, IMAX + 1)}
    print(f"   {'i':>2} {'3^i R_e(0)':>12} {'4*3^i R_e(2)':>13} {'X_i':>10} {'d[3^iR0]':>10} {'d[4*3^iR2]':>11}")
    for i in range(1, IMAX + 1):
        d0 = A0[i] - A0[i - 1] if i > 1 else float('nan')
        d2 = 4 * A2[i] - 4 * A2[i - 1] if i > 1 else float('nan')
        print(f"   {i:>2} {A0[i]:>12.6f} {4*A2[i]:>13.6f} {Xfull[i]:>10.5f} {d0:>10.6f} {d2:>11.6f}")
    ii = np.arange(7, IMAX + 1)
    s0 = np.polyfit(ii, [A0[i] for i in ii], 1)[0]
    s2 = np.polyfit(ii, [4 * A2[i] for i in ii], 1)[0]
    sx = np.polyfit(ii, [Xfull[i] for i in ii], 1)[0]
    print(f"   slope[3^iR_e(0)]={s0:.6f}   slope[4*3^iR_e(2)]={s2:.6f}   residual={s2-s0:+.6f}   "
          f"slope[X_i]={sx:.6f}  (fit i=7..{IMAX})\n")

    # ----- (3) T_i -> S_inf/2, monotone? which side of 7/30 -----
    print(f"## (3) T_i = 1/2 S_{{i+1}} -> S_inf/2 : trajectory vs 7/30={7/30:.6f} (7/15 target); T_1=5/21={5/21:.6f}")
    print(f"   {'i':>2} {'T_i':>12} {'Lambda_i=T_i-T_(i-1)':>20} {'T_i-7/30':>12}")
    Tprev = 1.0 / 3
    Tvals = {}
    for i in range(1, IMAX + 1):
        Ti = 3 ** i * (4 * R2[i] - R0[i]); Tvals[i] = Ti
        print(f"   {i:>2} {Ti:>12.8f} {Ti-Tprev:>+20.8f} {Ti-7/30:>+12.8f}")
        Tprev = Ti
    partial = sum(Tvals[i] - Tvals[i - 1] for i in range(2, IMAX + 1))  # Sum_{2..IMAX} Lambda_i
    print(f"\n   Sum_{{i=2..{IMAX}}} Lambda_i = {partial:+.6f}   (7/15 needs T_inf-T_1 = 7/30-5/21 = -1/210 = {-1/210:+.6f})")
    print(f"   T_{IMAX} = {Tvals[IMAX]:.6f}; T_{IMAX}-7/30 = {Tvals[IMAX]-7/30:+.6f}  "
          f"(>0 => S_inf>7/15 unless turnover)")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
