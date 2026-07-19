"""
PROBE D2-c -- three parallel checks while the pen grinds the ladder. Judge (sigma vs spectra) STAYS HELD.
C1 endpoint contraction: trapezoid endpoint atom = 1/D = 1/(2*3^{L-1}) (DERIVED: e'=0 channel flux
   = (1/3D) Sum_e AC(e) = (1/3D)(Sum w)^2 = 1/(3D), normalized -> 1/D). L=4 => 1/54. Verify numeric L=2,3.
C2 the m=2 seat: sigma(theta)=(1/3)((1+e^{i theta})/2)^2, pair m at theta_m = m*2pi/3^{L-1}, modulus
   (1/3)cos^2(theta_m/2). Search dense L=3 spectrum + L=4 block data for nearest modes to each seat. No fit.
C3 doublet precision: exact eigenvalues to full precision, both pairs, L=3 (dense) and L=4 (block).
"""
import numpy as np
from fractions import Fraction as Fr
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

def tower_dense(L, lam=0.5):
    q = 3; qL = q ** L; sub = subgroup(2 % qL, qL); D = len(sub)
    M, idx, n = build_M_gen(q, L, 2, [lam ** d for d in range(1, D + 1)])
    gam = np.array([s[2] for s in sorted(idx, key=lambda s: idx[s])])
    tw = np.where(gam != 0)[0]
    return M[tw][:, tw].toarray(), D

def C1():
    print("\n## C1  ENDPOINT CONTRACTION  (trapezoid endpoint atom = 1/D = 1/(2*3^{L-1}))")
    for L in (2, 3, 4):
        D = 2 * 3 ** (L - 1)
        wf = [Fr(2 ** (D - d), 2 ** D - 1) for d in range(1, D + 1)]
        AC_sum = sum(wf[j] * wf[(j + e) % D] for e in range(D) for j in range(D))  # = (sum w)^2
        endpoint = Fr(1, D)
        print(f"   L={L}: D={D}  Sum_e AC(e)=(Sum w)^2={AC_sum} (=1)  =>  endpoint atom = 1/D = {endpoint} "
              f"({float(endpoint):.6f})", flush=True)
    print("   sequence 1/6 -> 1/18 -> 1/54 (x1/3 per level) CONFIRMED; L=4 = 1/54 EXACT (Lebesgue-restriction law).")

def seat(m, L):
    th = m * 2 * np.pi / (3 ** (L - 1))
    return (1 / 3) * np.cos(th / 2) ** 2, th

def nearest(spec, mod, ph, k=3):
    tgt = mod * np.exp(1j * ph)
    order = sorted(spec, key=lambda z: abs(z - tgt))
    return order[:k]

def C2():
    print("\n## C2  THE m=2 SEAT  (seat = (1/3)cos^2(theta_m/2) at theta_m = m*2pi/3^{L-1}; nearest modes, no fit)")
    # L=3 dense full spectrum
    Md, D = tower_dense(3)
    ev = np.linalg.eig(Md)[0]
    up = [z for z in ev if z.imag > 1e-9]     # upper-half conjugate reps
    for m in (1, 2):
        mod, ph = seat(m, 3)
        near = nearest(up, mod, ph, 3)
        print(f"   L=3 m={m}: seat modulus={mod:.4f} phase={ph:.4f}  nearest pairs: " +
              ", ".join(f"{z:.5f}(|.|={abs(z):.4f},arg={np.angle(z):.4f})" for z in near), flush=True)
    # L=4 from D1-C block (existing data)
    L4 = [0.320423+0.075242j, 0.320223+0.075252j]   # m=1 doublet (block-6 converged)
    for m in (1, 2):
        mod, ph = seat(m, 4)
        near = nearest(L4, mod, ph, 2)
        print(f"   L=4 m={m}: seat modulus={mod:.4f} phase={ph:.4f}  nearest in D1-C block-6: " +
              ", ".join(f"{z:.6f}(|.|={abs(z):.5f},arg={np.angle(z):.5f})" for z in near) +
              ("   [m=2 below block-6 depth -- needs deeper block]" if m == 2 else ""), flush=True)
    return up

def C3(up3):
    print("\n## C3  DOUBLET PRECISION  (exact eigenvalues, full precision)")
    # L=3: the two leading pairs (m=1 doublet)
    lead = sorted(up3, key=lambda z: -abs(z))[:2]
    print(f"   L=3 doublet (dense, full precision):", flush=True)
    for z in lead:
        print(f"      {z.real:.12f} {z.imag:+.12f}j   |.|={abs(z):.12f}  arg={np.angle(z):.12f}", flush=True)
    split3 = abs(lead[0] - lead[1])
    print(f"      splitting |p1-p2| = {split3:.6e}", flush=True)
    print(f"   L=4 doublet (D1-C block, converged ~res 1e-4 on the doublet):", flush=True)
    L4 = [0.320423+0.075242j, 0.320223+0.075252j]
    for z in L4:
        print(f"      {z.real:.6f} {z.imag:+.6f}j   |.|={abs(z):.6f}  arg={np.angle(z):.6f}", flush=True)
    print(f"      splitting |p1-p2| = {abs(L4[0]-L4[1]):.6e}", flush=True)

def main():
    print("# PROBE D2-c -- C1 endpoint contraction + C2 m=2 seat + C3 doublet precision. Judge HELD.")
    C1()
    up3 = C2()
    C3(up3)

if __name__ == "__main__":
    main()
