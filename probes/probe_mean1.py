"""
PROBE MEAN1 -- verify Wilson's mean-1 constraint: mean_k gamma_r(k) = 1 exactly, and the M+/M- split (2026-07-26).

Since 4 generates G_r, as k runs 0..3^r-1 the ratio 4^k runs over the whole group => Sum_k C_r(k) = (Sum_s rho(s))^2 = 1
(rho normalized) => Sum_k gamma_r(k) = 3^r => mean_k gamma_r(k) = 1 EXACTLY at every r. White is the pinned MEAN.
Split: 1/3 of channels are 3|k (enriched, mean M+), 2/3 are 3-nmid (depleted, mean M-); (1/3)M+ + (2/3)M- = 1
=> M- ~ 0.75 forces M+ = 3 - 2 M- ~ 1.5. Wilson's small-k enriched sample averages 1.56.

Full autocorrelation via one FFT: C_r = ifft(|fft(rho_r)|^2). Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def main():
    t0 = time.time()
    print("# PROBE MEAN1 -- mean_k gamma_r(k) = 1 exactly; M+/M- split\n")
    for r in (14, 15, 16):
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho = rho / rho.sum()
        N = 3 ** r
        C = np.fft.ifft(np.abs(np.fft.fft(rho)) ** 2).real     # full autocorrelation, all lags
        gam = 3.0 ** r * C
        mean = gam.mean()
        k = np.arange(N)
        enr = (k % 3 == 0); dep = ~enr
        Mp = gam[enr].mean(); Mm = gam[dep].mean()
        # exclude k=0 (=X_r) from enriched to see the finite-channel enriched mean
        enr_no0 = enr.copy(); enr_no0[0] = False
        Mp_no0 = gam[enr_no0].mean()
        print(f"   r={r}: mean_k gamma = {mean:.10f}  (=1? {abs(mean-1)<1e-6})  Sum/3^r = {gam.sum()/N:.10f}")
        print(f"         M+ (3|k, incl k=0) = {Mp:.5f} ; M+ (excl k=0) = {Mp_no0:.5f} ; M- (3-nmid) = {Mm:.5f}")
        print(f"         identity (1/3)M+ + (2/3)M- = {Mp/3 + 2*Mm/3:.8f}  ; 3-2M- = {3-2*Mm:.5f} (predicted M+)")
        print(f"         gamma(0)=X_{r} = {gam[0]:.4f} ; #enriched={enr.sum()} #depleted={dep.sum()} "
              f"(1/3 vs 2/3: {enr.sum()/N:.4f})")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
