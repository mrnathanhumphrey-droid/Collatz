"""
OPTION 1 (diagnostic, read-only) -- is the stopping-time renewal residue eps_S(k)
log-periodic with the SAME log2/log3 incommensurate scales as the Lambda_i tower (LATTICE)?

LATTICE (tower, indexed by level i, modulus x=3^i): base-3 lattice mode ALIASED at integer i,
base-2 (div 2^v) VISIBLE at period ~9 = 2pi/log2; log3/log2=1.585 irrational => quasi-periodic.

eps_S is the exact MIRROR: indexed by k=log2(N). A log-periodic-base-b signal in log(N)
has period log(b)/log2 in k:
   base 2 -> period 1.000   (aliased to DC at INTEGER k = the old 2^32,34,36 data -> looked like noise)
   base 3 -> period 1.585   (= log3/log2, the same irrational; aliases to 2.71 at integer k)
Sampling k at 0.5 spacing (Nyquist period 1.0) resolves BOTH 1.000 and 1.585 directly.

TEST: periodogram of eps_S(k). Clean peak at 1.585 (or 1.000) => SAME mechanism (arcs bridged).
Flat/noise => eps_S oscillation is a finite-N convergence artifact, unrelated to the tower.

eps_S(k) sampled per-octave band [2^(k-1), 2^k): removes the huge cross-octave sigma_S variance
=> SE ~ 0.003 at a few M orbits/band instead of 50M.
"""
import time, numpy as np
from numba import njit, prange

LOG2 = np.log(2.0); LOG3 = np.log(3.0); LOG43 = 2*LOG2 - LOG3


@njit(parallel=True, cache=True)
def band_residuals(starts, max_value, max_steps):
    n = len(starts); nc = 64; chunk = (n + nc - 1)//nc
    s1 = np.zeros(nc); s2 = np.zeros(nc); cnt = np.zeros(nc, np.int64); ovf = np.zeros(nc, np.int64)
    for c in prange(nc):
        lo = c*chunk; hi = min(lo+chunk, n)
        for i in range(lo, hi):
            x = np.int64(starts[i]); logx = np.log(np.float64(x))
            sig = 0; steps = 0; failed = False
            while x != 1 and steps < max_steps:
                if x & 1:
                    if x > max_value//3: failed = True; break
                    x = 3*x + 1; sig += 1
                else:
                    x >>= 1
                steps += 1
            if failed or x != 1:
                ovf[c] += 1; continue
            r = sig - logx/LOG43            # per-orbit renewal residual
            s1[c] += r; s2[c] += r*r; cnt[c] += 1
    return s1.sum(), s2.sum(), cnt.sum(), ovf.sum()


def main():
    # IRRATIONAL k-grid: spacing sqrt(5)-2 = 0.236068 (incommensurate w/ 1.0 AND 1.585).
    # => any octave-binning systematic (period-1.0 in frac(k)) DECOHERES; true log-periods survive.
    delta = np.sqrt(5.0) - 2.0
    ks = 16.0 + delta*np.arange(0, int((32.0-16.0)/delta)+1)
    M = 10_000_000
    rng = np.random.default_rng(7)
    print(f"# OPTION 1  eps_S(k) per-octave sweep  ({len(ks)} bands, M={M:,}/band, spacing 0.5)\n")
    # warm up numba
    _ = band_residuals(np.array([3,5,7,9], np.int64), np.int64(1<<60), 100000)

    K = []; E = []; SE = []
    t0 = time.perf_counter()
    for k in ks:
        lo = int(2.0**(k-1)); hi = int(2.0**k)
        if lo % 2 == 0: lo += 1
        starts = rng.integers(lo//2, hi//2, size=M, dtype=np.int64)*2 + 1
        s1, s2, cnt, ovf = band_residuals(starts, np.int64(1<<60), 400000)
        mean = s1/cnt; var = s2/cnt - mean*mean; se = np.sqrt(var/cnt)
        K.append(k); E.append(mean); SE.append(se)
        print(f"   k={k:5.2f}  eps_S={mean:8.5f}  SE={se:.5f}  ovf={ovf}")
    print(f"\n   sweep {time.perf_counter()-t0:.1f}s")
    K = np.array(K); E = np.array(E); SE = np.array(SE); W = 1.0/SE**2
    np.savetxt("/c/tmp/eps_lattice.csv", np.column_stack([K,E,SE]), header="k eps_S SE", comments="")

    # flexible smooth trend a + b/k + c/k^2 (captures the 1/N transient; smooth => cannot fake short periods)
    A0 = np.column_stack([np.ones_like(K), 1.0/K, 1.0/K**2])
    cw = np.linalg.lstsq(A0*np.sqrt(W)[:,None], E*np.sqrt(W), rcond=None)[0]
    trend = A0 @ cw
    chi2_0 = np.sum(W*(E-trend)**2)
    print(f"\n   trend  eps_S = {cw[0]:.5f} + {cw[1]:+.3f}/k + {cw[2]:+.2f}/k^2   dof={len(K)-3}")
    print(f"   chi2(trend only) = {chi2_0:.2f}\n")

    def fit_at(P):
        c = np.cos(2*np.pi*K/P); s = np.sin(2*np.pi*K/P)
        Ad = np.column_stack([np.ones_like(K), 1.0/K, 1.0/K**2, c, s])
        cc = np.linalg.lstsq(Ad*np.sqrt(W)[:,None], E*np.sqrt(W), rcond=None)[0]
        chi2 = np.sum(W*(E - Ad@cc)**2)
        cov = np.linalg.inv((Ad*W[:,None]).T @ Ad)
        amp = np.hypot(cc[3], cc[4]); ampse = np.sqrt(cov[3,3]+cov[4,4])
        return amp, ampse, chi2_0 - chi2

    # SHORT-period periodogram (P in [0.9,4]; trend leakage lives at long P only)
    print("   ## short-period periodogram (P<=4; long-P trend-leakage excluded)")
    print(f"   {'P':>6} {'amp':>8} {'amp/SE':>7} {'dchi2':>8}   note")
    Ps = np.arange(0.90, 4.0001, 0.005); best = (0,0,0)
    for P in Ps:
        amp, ampse, d = fit_at(P)
        if d > best[0]: best = (d, P, amp)
    marks = {1.000:"base-2 (log2/log2)", 1.585:"base-3 (log3/log2) <-- HYP",
             2.710:"base-3 aliased@int-k", 3.170:"2*1.585 harmonic"}
    for Pm, note in marks.items():
        amp, ampse, d = fit_at(Pm)
        print(f"   {Pm:>6.3f} {amp:>8.5f} {amp/ampse:>7.1f} {d:>8.2f}   {note}")
    d,P,a = best; amp,ampse,_ = fit_at(P)
    print(f"\n   SHORT-BAND PEAK: P={P:.3f}  amp={amp:.5f}  ({amp/ampse:.1f} sigma)  dchi2={d:.2f}")
    print(f"   (predicted base-3 log-period = log3/log2 = {LOG3/LOG2:.4f})")


if __name__ == "__main__":
    main()
