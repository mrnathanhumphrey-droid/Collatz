"""
probe_qx1_neg1_coset_2026_05_29.py

Refines the q=7 Plancherel null. q=7's exact 50/50 split was SYMMETRY-PROTECTED:
-1 not in <2> mod 7 (ord_2=3 odd) => the dark coset IS the conjugate mirror -<2>,
and reality of P (mu_hat(-xi)=conj mu_hat(xi)) forces equal mass. That is special.

KEY DISTINCTION:
  -1 in <2>  (e.g. q=17: ord_2 mod 17 = 8, 2^4 = -1):
      conjugation maps each <2>-coset to ITSELF => masses NOT symmetry-constrained
      => half-mixing dark-coset suppression CAN appear in Plancherel mass.
  -1 not in <2> (e.g. q=7,23: ord_2 odd):
      conjugation pairs <2> with its dark coset => masses LOCKED equal (50/50),
      half-mixing invisible to mass (need phase) AND the dark coset is a pure
      conjugate replica (no independent info to expose).

PREDICTION:
  q=11 full (idx 1): no dark coset.
  q=7, q=23 half, -1 not in <2>: EXACT 50/50 (symmetry-locked).
  q=17 half, -1 in <2>: GENUINE asymmetry, mass(<2>) > mass(dark).
"""
from __future__ import annotations
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
from probe_qx1_coset_plancherel_2026_05_29 import offset_distribution, subgroup_pow2

def main():
    print("=" * 104)
    print("qx+1 coset Plancherel mass, classified by whether -1 in <2> (the symmetry-protection switch)")
    print("=" * 104)
    print(f"{'q':>3} {'n':>2} {'N=q^n':>8} {'ord2':>6} {'idx':>4} {'-1 in<2>':>8} "
          f"{'mass_<2>':>10} {'mass_dark':>10} {'maxH/maxDark':>13} {'argmax in<2>':>12}")
    panel = [11, 7, 23, 17]
    for q in panel:
        for n in range(2, 5):
            N = q ** n
            P = offset_distribution(q, n)
            assert abs(P.sum() - 1.0) < 1e-9
            mu = np.fft.fft(P)
            xi = np.arange(N)
            units = (xi % q != 0)
            H, ord2 = subgroup_pow2(q, n)
            phi = q ** (n - 1) * (q - 1)
            idx = phi // ord2
            neg1_in_H = bool(H[(N - 1) % N])
            inH = units & H
            dark = units & (~H)
            a2 = np.abs(mu) ** 2
            mtot = a2[units].sum()
            mH = a2[inH].sum() / mtot
            mD = a2[dark].sum() / mtot if dark.any() else 0.0
            maxH = np.abs(mu[inH]).max() if inH.any() else 0.0
            maxD = np.abs(mu[dark]).max() if dark.any() else float('nan')
            ratio = (maxH / maxD) if (dark.any() and maxD > 0) else float('inf')
            amu = np.abs(mu).copy(); amu[~units] = -1
            argmax_inH = bool(H[int(np.argmax(amu))])
            print(f"{q:>3} {n:>2} {N:>8} {ord2:>6} {idx:>4} {str(neg1_in_H):>8} "
                  f"{mH:>10.6f} {mD:>10.6f} {ratio:>13.4g} {str(argmax_inH):>12}")
        print()

    print("Reading:")
    print("  q=11 (idx 1): single coset, mass_<2>=1 by construction (full mixing).")
    print("  q=7,23 (-1 not in <2>): expect EXACT 50/50, maxH/maxDark=1 -> symmetry-locked, mass-invisible.")
    print("  q=17 (-1 in <2>): expect mass_<2> > mass_dark -> half-mixing VISIBLE in Plancherel mass.")
    print("  => 'Plancherel mass cannot see ord_2(q)' was a q=7 artifact; the true switch is -1 in <2>.")

if __name__ == "__main__":
    main()
