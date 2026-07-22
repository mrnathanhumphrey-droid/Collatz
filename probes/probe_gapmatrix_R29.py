"""
PROBE R29 -- THE GAP MATRIX. The unblocked eigenproblem: the operator acts on {R(2m)}, a FIXED gap-index lattice
(odd gaps dead), same at every r -- so the finite matrix R26 couldn't build (growing spaces) IS buildable.

  M_r[m,m'] = 3 * P(2(m'-m)) * kappa_{2m}(r),  kappa_{2m}(r) = R_r(2m) / [3 sum_{m'} P(2(m'-m)) R_{r-1}(2m')].
  kappa_0 = 1 exactly (Diagonal Flatness: 3 sum_{m'} P(2m') R_{r-1}(2m') = 3 X_r = R_r(0)).
  M_r maps R_{r-1}-vector -> R_r-vector exactly (by construction of kappa). If kappa stabilizes in r, M_r -> M and
  its eigenvalues are the transfer rates: leading -> 1 (=rho at criticality, X grows linearly), subdominant = |lam2|.

C first (the real risk): does kappa_{2m}(r) stabilize? A: diagonalize. B: truncation D. D: deflation residual.
"""
import os, sys, math, json, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_gapop_R28 import build_nu, nu_hat, R_of_d

lam = 0.5
P = lambda d: (1 - lam) * lam ** abs(d) / (1 + lam)


def main():
    print("# PROBE R29 -- THE GAP MATRIX (lam=1/2).\n")
    hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                        'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in hist.items()}
    Sb = {k: float(F(7, 15) + EPS[k]) for k in EPS}
    Lam = {r: float((EPS[r + 1] - EPS[r]) / 2) for r in range(1, 8)}

    nus = build_nu(lam, 7)
    NH = {r: nu_hat(nus[r], 3 ** (r + 1)) for r in range(0, 8)}
    MM = 12                                             # compute R(2m) for |m|<=MM
    Rtab = {}                                           # Rtab[r][m] = R_r(2m)
    for r in range(0, 8):
        M1 = 3 ** (r + 1)
        Rtab[r] = {m: R_of_d(NH[r], M1, 2 * m).real for m in range(-MM, MM + 1)}

    def kappa(r, m, D):
        num = Rtab[r][m]
        den = 3 * sum(P(2 * (mp - m)) * Rtab[r - 1][mp] for mp in range(-D, D + 1))
        return num / den if den else float('nan')

    # ================= R29-C (run first) =================
    print("## R29-C  DOES kappa_{2m}(r) STABILIZE?  (the real risk; D=10)  [kappa_0 should be 1 exactly]")
    D = 10
    print(f"   {'r':>2} {'kappa_0':>9} {'kappa_2 (m=1)':>13} {'kappa_4 (m=2)':>13} {'kappa_6 (m=3)':>13} "
          f"{'chan1/S_r':>10} {'chan2/S_r':>10} {'chan3/S_r':>10}")
    for r in range(2, 8):
        k0 = kappa(r, 0, D); k1 = kappa(r, 1, D); k2 = kappa(r, 2, D); k3 = kappa(r, 3, D)
        # channel-m contribution to S_r = 2 P(2m) R_{r-1}(2m) / S_r
        c1 = 2 * P(2) * Rtab[r - 1][1] / Sb[r]
        c2 = 2 * P(4) * Rtab[r - 1][2] / Sb[r]
        c3 = 2 * P(6) * Rtab[r - 1][3] / Sb[r]
        print(f"   {r:>2} {k0:>9.5f} {k1:>13.5f} {k2:>13.5f} {k3:>13.5f} {c1:>10.5f} {c2:>10.5f} {c3:>10.5f}")
    print("   [stabilizes => M has a well-defined limit, A is legitimate. drifts => obstruction, report it.]\n")

    # ================= R29-A =================
    print("## R29-A  BUILD & DIAGONALIZE  M_r[m,m']=3 P(2(m'-m)) kappa_{2m}(r); spectrum per D, per r")
    print("   PRE-REG: leading eigenvalue -> 1 (=rho, sanity gate); |lam2| -> 1/2; arg -> period 2pi/arg")
    for D in (4, 6, 8, 10):
        idx = list(range(-D, D + 1))
        n = len(idx)
        print(f"   --- D={D} ({n} states) ---")
        for r in (5, 6, 7):
            Mr = np.zeros((n, n))
            for i, m in enumerate(idx):
                km = kappa(r, m, D)
                for j, mp in enumerate(idx):
                    Mr[i, j] = 3 * P(2 * (mp - m)) * km
            ev = np.linalg.eigvals(Mr)
            ev = sorted(ev, key=lambda z: -abs(z))
            lead = ev[0]
            lam2 = ev[1]
            arg2 = math.degrees(cmath.phase(lam2)) if abs(np.imag(lam2)) > 1e-9 else 0.0
            per = (360 / abs(arg2)) if abs(arg2) > 1e-9 else float('inf')
            top = ", ".join(f"{abs(z):.4f}" + (f"@{math.degrees(cmath.phase(z)):.0f}d" if abs(np.imag(z)) > 1e-9 else "") for z in ev[:5])
            print(f"     r={r}: leading={lead.real:+.5f}  |lam2|={abs(lam2):.5f}"
                  + (f" arg={arg2:.2f}d period={per:.3f}" if abs(np.imag(lam2)) > 1e-9 else " (real)")
                  + f"   top5|.|: {top}")
    print()

    # ================= R29-B =================
    print("## R29-B  TRUNCATION CONVERGENCE: does |lam2| stabilize as D grows? (r=7)")
    r = 7
    prev = None
    for D in (4, 6, 8, 10):
        idx = list(range(-D, D + 1)); n = len(idx)
        Mr = np.zeros((n, n))
        for i, m in enumerate(idx):
            km = kappa(r, m, D)
            for j, mp in enumerate(idx):
                Mr[i, j] = 3 * P(2 * (mp - m)) * km
        ev = sorted(np.linalg.eigvals(Mr), key=lambda z: -abs(z))
        l2 = abs(ev[1])
        dd = f"(D-step change {abs(l2-prev):.2e})" if prev else ""
        print(f"   D={D:>2}: leading={ev[0].real:.5f}  |lam2|={l2:.6f}  {dd}")
        prev = l2
    print("   [D=8,D=10 agree to several digits => legitimate. disagree => truncation illegitimate, A void.]\n")

    # ================= R29-D =================
    print("## R29-D  RELATIVE DEFLATION RESIDUAL: mu_r/|Lambda_r| for r=3,4 at mode z=0.49,0.5,0.503")
    print(f"   local ratios: Lambda_4/Lambda_3={Lam[4]/Lam[3]:.6f}  Lambda_5/Lambda_4={Lam[5]/Lam[4]:.6f} (bracket 0.5?)")
    for z in (0.49, 0.5, 0.503):
        row = []
        for r in (3, 4):
            mu = Lam[r + 1] - z * Lam[r]
            row.append(f"r{r}: mu/|Lam_r|={mu/abs(Lam[r]):+.2e}")
        print(f"   z={z}: " + "   ".join(row))
    print("   [which z minimizes |mu/Lam|? => direct read on whether |lam2|=1/2 exactly]")


if __name__ == "__main__":
    main()
