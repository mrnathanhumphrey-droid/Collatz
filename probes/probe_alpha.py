"""
PROBE ALPHA -- Wilson's final cheap measurement: the four alpha_c and their common argument (2026-07-26).

The isotropy of the direction distribution (CARRYMAT) kills the anisotropic |beta| part outright, so the target is
    sign( Sum_c p_c Re(omega^c alpha_c) ),   alpha_c = regression coeff of z_{4x} on z_x within cell c.
Cell-skew Sum_c p_c omega^c = 0.2425+0.0163i (modulus 0.243 = Wilson's 1/4) -- already confirmed from banked p_c.
Open datum: arg(alpha_c) -- the typical PHASE ROTATION of the digit profile under x4. Positivity needs it inside a
quarter turn. Wilson's rough estimate ~80deg (= 4pi/9, ninths natural for 3-adic digits). Recognizable => closed-form
structural constant (first in the MECHANISM not the bookkeeping); nondescript 78.6 => a number, done looking.

z_x = complex coord of P(dpi_x) (aX), z_{4x} = complex coord of P(dpi_{4x}) (aL) in CARRYMAT's orthonormal basis.
alpha_reg = E_mu[z_{4x} conj(z_x)] / E_mu[|z_x|^2]  (the sign-carrying coeff, since E[z^2]~0 & resid _|_).
alpha_ratio = E_mu[z_{4x}/z_x]  (Wilson's literal E[z4/z], reported for comparison).
Reuses probe_carrymat.cell_data. One build_nu(0.5,16) ~6 min.
"""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_gapop_R28 import build_nu
from probe_carrymat import cell_data, dense

RTOP = 16


def main():
    t0 = time.time()
    print("# PROBE ALPHA -- the four alpha_c, their common argument, phase rotation of the profile under x4\n")
    print(f"building build_nu to {RTOP} ... (~6 min)")
    nus = build_nu(0.5, RTOP)
    print(f"  built ({time.time()-t0:.1f}s)\n")

    r = RTOP
    hi = dense(nus[r], 3 ** (r + 1))
    c, cm, wgt, aX, bC, aL, _ = cell_data(hi, r)
    del hi
    z = aX[:, 0] + 1j * aX[:, 1]
    z4 = aL[:, 0] + 1j * aL[:, 1]
    Z = wgt.sum()
    p = np.array([wgt[c == j].sum() / Z for j in range(4)])
    w2 = np.exp(2j * np.pi / 3)
    skew = sum(p[j] * w2 ** j for j in range(4))
    print(f"## cell-skew  Sum_c p_c omega^c = {skew:.4f}   |.| = {abs(skew):.4f}  arg = {np.degrees(np.angle(skew)):+.1f}deg"
          f"   (Wilson's 1/4)\n")

    print("## alpha_c per cell  (regression = sign-carrying; ratio = Wilson's literal E[z4/z])")
    print(f"   {'c':>1} {'rot':>3} {'p_c':>6} | {'|a_reg|':>8} {'arg a_reg':>9} | {'|a_ratio|':>9} {'arg a_ratio':>11} | {'E[z^2]/E|z|^2':>13}")
    a_reg = {}
    for j in range(4):
        sel = (c == j); w = wgt[sel]
        zc, z4c = z[sel], z4[sel]
        sw = w.sum()
        Ezz = float((w * np.abs(zc) ** 2).sum() / sw)
        areg = complex((w * z4c * np.conj(zc)).sum() / sw / Ezz)
        # ratio estimator, guard tiny |z|
        good = np.abs(zc) > 1e-9
        aratio = complex((w[good] * (z4c[good] / zc[good])).sum() / w[good].sum())
        Ez2 = complex((w * zc ** 2).sum() / sw) / Ezz          # phase-uniformity check (should be ~0)
        a_reg[j] = (areg, Ezz)
        print(f"   {j:>1} {j%3:>3} {p[j]:>6.4f} | {abs(areg):>8.4f} {np.degrees(np.angle(areg)):>+8.1f} | "
              f"{abs(aratio):>9.4f} {np.degrees(np.angle(aratio)):>+10.1f} | {abs(Ez2):>13.4f}")
    print()

    # reconstruction: U_c ?= Re(omega^{+/-c} alpha_reg) * E|z|^2 ; pick sign that matches, validate sum = q-1/3
    print("## reduction check: U_c = Re(omega^s*c * alpha_reg_c) * E|z|^2 ,  Sum p_c U_c ?= q-1/3 = +4.1789e-4")
    for s in (+1, -1):
        recon = 0.0; rows = []
        for j in range(4):
            areg, Ezz = a_reg[j]
            Uc = (w2 ** (s * j) * areg).real * Ezz
            recon += p[j] * Uc
            rows.append(Uc)
        tag = "MATCH" if abs(recon - 4.1789e-4) / 4.1789e-4 < 0.05 else "no"
        print(f"   s={s:+d}: U_c = [{', '.join(f'{x:+.3e}' for x in rows)}]  Sum p_c U_c = {recon:+.4e}  [{tag}]")
    print()

    # common argument
    args = [np.degrees(np.angle(a_reg[j][0])) for j in range(4)]
    print(f"## arg alpha_reg per cell: {['%+.1f' % a for a in args]}")
    print(f"   named phase targets: 4pi/9 = {np.degrees(4*np.pi/9):.1f}deg,  2pi/9 = {np.degrees(2*np.pi/9):.1f}deg,  "
          f"pi/4 = 45deg,  pi/3 = 60deg,  2pi/5 = 72deg")
    # weighted-mean modulus/phase of z4/z overall
    good = np.abs(z) > 1e-9
    aall = complex((wgt[good] * (z4[good] / z[good])).sum() / wgt[good].sum())
    print(f"   overall E_mu[z4/z] = {aall:.4f}  |.|={abs(aall):.4f}  arg={np.degrees(np.angle(aall)):+.1f}deg")
    json.dump({'p': p.tolist(), 'skew': [skew.real, skew.imag],
               'alpha_reg': {str(j): [a_reg[j][0].real, a_reg[j][0].imag, a_reg[j][1]] for j in range(4)},
               'args': args},
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                'experiments_output', 'alpha.json'), 'w'))
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
