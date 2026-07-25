"""
PROBE POINCARE -- Wilson's pen reduction, verified + the dilation falsifier.

(1) POINCARE FORM (the solid deliverable, verify it):
    P = kill-3|m-modes multiplier. PC = q * q~ with q = P rho (fiber-fluctuation profile), because
    1_{3nmid m}|rho_hat|^2 = |1_{3nmid m} rho_hat|^2.  Target = (q*q~)(1) = sum_k q(k)q(k+1) > 0.
    Since sum q(k)q(k+1) = ||q||^2 - 1/2||Dq||^2  (Dq(k)=q(k+1)-q(k)):
        Re dhat_r(1) = 1 - ||Dq||^2 / (2||q||^2),   so  Re dhat_r(1)>0  <=>  ||Dq||^2 < 2||q||^2.
    Discrete Poincare on the fiber-fluctuation subspace; Rayleigh ||Dq||^2/||q||^2 in [0,4], threshold 2 = white
    (median-at-N/4). MARGINAL BY CONSTRUCTION: margin under 2 is exactly 2*d1. GATE: match banked d1 to ~1e-12.

(2) PARITY FALSIFIER (Wilson's cheap check) -- but m->N-m flips parity (N odd) and preserves Chat*cos, so
    Num_even = Num_odd = Num/2 EXACTLY, for ANY parity cut. Confirm numerically => the parity route is degenerate.

(3) DILATION-BAND split (the NON-degenerate version of Wilson's mechanism): m' = 2^{-1} m mod N is the coordinate
    the transport modulation |g_hat|^2 = 1/(5-4cos(2pi m'/N)) low-passes in. Split Num by |m'|<N/4 (modulation large)
    vs |m'|>N/4 (small). This cut is conjugation-INVARIANT (m'->-m'), so NOT symmetry-forced. Does the dilation
    organize the +-0.319 band cancellation?
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"
D1_BANK = {12: 2.963565168845e-3, 13: 2.696227117865e-3, 14: 2.440864723362e-3,
           15: 2.187132833772e-3, 16: 1.939224678822e-3}


def main():
    t0 = time.time()
    print("# PROBE POINCARE -- Poincare form (verify) + parity(degenerate) + dilation-band falsifier\n")
    print("## (1) POINCARE GATE:  Re dhat_r(1) = 1 - ||Dq||^2/(2||q||^2),  q = P rho")
    print(f"   {'r':>2} {'||q||^2':>13} {'||Dq||^2':>13} {'Rayleigh':>10} {'1-R/2':>14} {'banked d1':>14} {'rel':>9} {'2-R=2d1':>10}")
    for r in (12, 13, 14, 15, 16):
        N = 3 ** r
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
        H = np.fft.fft(rho); del rho
        H[0::3] = 0                                      # kill 3|m modes (P); index set symmetric => q real
        q = np.fft.ifft(H).real; del H
        nq2 = float(np.dot(q, q))
        dq = np.roll(q, -1) - q
        ndq2 = float(np.dot(dq, dq)); del dq
        lag1 = float(np.dot(q, np.roll(q, -1)))         # (q*q~)(1) directly
        R = ndq2 / nq2
        d1p = 1 - R / 2
        d1_direct = lag1 / nq2                           # must equal d1p
        rel = abs(d1p - D1_BANK[r]) / D1_BANK[r]
        print(f"   {r:>2} {nq2:>13.6e} {ndq2:>13.6e} {R:>10.6f} {d1p:>14.9e} {D1_BANK[r]:>14.9e} {rel:>9.1e} {2-R:>10.4e}"
              f"  [lag1/||q||^2={d1_direct:+.9e}]")
        del q
    print("   [Rayleigh threshold 2 = white; margin under 2 is exactly 2*d1 ~ 0.004 => MARGINAL BY CONSTRUCTION.]\n")

    print("## (2) PARITY FALSIFIER: Num_even vs Num_odd  (predicted EXACTLY Num/2 by m->N-m conjugate symmetry)")
    print(f"   {'r':>2} {'Num':>13} {'Num_even':>13} {'Num_odd':>13} {'even-odd':>11} {'(even-Num/2)/Num':>17}")
    for r in (12, 14, 16):
        N = 3 ** r
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
        prof = np.abs(np.fft.fft(rho)) ** 2; del rho
        m = np.arange(N)
        prim = (m % 3 != 0)
        cosw = np.cos(2 * np.pi * m / N)
        pm = prof[prim]; mm = m[prim]; cw = cosw[prim]
        num = float((pm * cw).sum())
        ev = (mm % 2 == 0)
        ne = float((pm[ev] * cw[ev]).sum()); no = float((pm[~ev] * cw[~ev]).sum())
        print(f"   {r:>2} {num:>13.5e} {ne:>13.5e} {no:>13.5e} {ne-no:>+11.2e} {(ne-num/2)/num:>17.2e}")
        del prof, pm, mm, cw, cosw
    print("   [even=odd=Num/2 => parity cannot organize the cancellation; Wilson's parity falsifier is DEGENERATE.]\n")

    print("## (3) DILATION-BAND: split Num by dilated coord m'=2^{-1}m mod N; |m'|<N/4 (mod large) vs |m'|>N/4 (small)")
    print("   (conjugation-invariant cut; the non-degenerate version of the dilation mechanism)")
    print(f"   {'r':>2} {'Num':>12} {'m<N/4 (big-mod)':>16} {'m>N/4 (small-mod)':>17} {'|resid|/|big|':>13} {'mass(m<N/4)':>12}")
    for r in (12, 14, 16):
        N = 3 ** r
        inv2 = pow(2, -1, N)
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
        prof = np.abs(np.fft.fft(rho)) ** 2; del rho
        m = np.arange(N)
        prim = (m % 3 != 0)
        pm = prof[prim]; mm = m[prim].astype(np.int64)
        cw = np.cos(2 * np.pi * mm / N)
        mp = (inv2 * mm) % N                              # dilated coord m'
        mpf = np.minimum(mp, N - mp)                      # |m'|
        big = mpf < (N // 4)                              # modulation large
        S = float(pm.sum())
        num = float((pm * cw).sum())
        nb = float((pm[big] * cw[big]).sum()); ns = float((pm[~big] * cw[~big]).sum())
        massb = float(pm[big].sum()) / S
        print(f"   {r:>2} {num:>12.5e} {nb:>+16.5e} {ns:>+17.5e} {abs(nb+ns)/max(abs(nb),abs(ns)):>13.4f} {massb:>12.5f}")
        del prof, pm, mm, cw, mp, mpf, big
    print("   [if big-mod & small-mod bands are large+opposing (|resid|/|big|<<1) => dilation organizes the cancellation,")
    print("    target becomes the dilation residue. if comparable to sum => dilation isn't the mechanism, route dead.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
