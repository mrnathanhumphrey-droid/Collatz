"""
PROBE MODREG (Wilson's cheap check, not a real probe) -- regress the fluctuation spectrum Chat(m) against the
pure transport modulation w(2^{-1} m mod N) over 3 nmid m, at r=12 and r=16.

w(u) = 1/(5 - 4 cos 2pi u),  u = m'/N,  m' = 2^{-1} m mod N   (the dilated coordinate the transport low-passes in).
Analytic: what = 2^{-|n|}/3, and since m=2m' the weight cos(2pi m/N)=cos(4pi u), <cos>_w = what(2)/what(0) = 1/4.
So the white part contributes 0 and the modulation contributes +1/4; d1's sign = sign of Chat's modulation-aligned part.

Report: OLS slope b (Chat ~ a + b*w over 3nmid m), sign, R^2 (=corr^2), and the implied modulation amplitude.
If b>0 and R^2 decent => d1>0 reduces to "the modulation component has positive amplitude" = a structural property
of the transport (provable). If R^2 poor => the ripple is not modulation-shaped and the 1/4 is a null coincidence.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def main():
    t0 = time.time()
    print("# PROBE MODREG -- regress Chat(m) vs w(2^{-1}m mod N) over 3 nmid m\n")
    print(f"   {'r':>2} {'slope b':>13} {'intercept a':>13} {'R^2':>9} {'corr':>8} {'b*Var(w)^.5/mean':>16} {'<cos>_check':>12}")
    for r in (12, 16):
        N = 3 ** r
        inv2 = pow(2, -1, N)
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
        prof = np.abs(np.fft.fft(rho)) ** 2; del rho
        m = np.arange(N)
        prim = (m % 3 != 0)
        y = prof[prim]                                  # Chat(m), fluctuation spectrum
        mm = m[prim].astype(np.int64)
        mp = (inv2 * mm) % N                            # dilated coord m'
        x = 1.0 / (5.0 - 4.0 * np.cos(2 * np.pi * mp / N))   # w(2^{-1}m)
        # centered OLS  y ~ a + b x
        xm = x.mean(); ym = y.mean()
        xc = x - xm; yc = y - ym
        b = float((xc * yc).sum() / (xc * xc).sum())
        a = float(ym - b * xm)
        ss_res = float(((yc - b * xc) ** 2).sum()); ss_tot = float((yc * yc).sum())
        R2 = 1 - ss_res / ss_tot
        corr = float((xc * yc).sum() / np.sqrt((xc * xc).sum() * (yc * yc).sum()))
        rel_ripple = b * np.sqrt(float((xc * xc).sum()) / len(x)) / ym    # ripple std / white level
        # independent <cos> reconstruction: d1 = sum Chat cos / sum Chat
        cosw = np.cos(2 * np.pi * mm / N)
        d1 = float((y * cosw).sum() / y.sum())
        print(f"   {r:>2} {b:>+13.5e} {a:>+13.5e} {R2:>9.4f} {corr:>+8.4f} {rel_ripple:>16.4e} {d1:>+12.5e}")
        del prof, y, mm, mp, x, cosw
    print("\n   [b>0 & R^2 decent => sign(d1)=sign(modulation amplitude)=structural(transport), provable.")
    print("    R^2 poor => ripple not modulation-shaped, 1/4 is a null coincidence. <cos>_check = banked d1.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
