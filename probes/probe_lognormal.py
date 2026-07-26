"""
PROBE LOGNORMAL (Wilson) -- is the pi-hat magnitude spectrum log-normal with sigma^2 ~ c*k? (2026-07-26)

VALPROFILE gave: gap ln(quad mean) - mean(ln) ~ (1/2)Var(ln|pi-hat|) grows linearly => Var(ln|pi-hat|) ~ c*k.
Linear-in-k variance = log-spectrum is a RANDOM WALK over levels (each single-level recursion factor multiplies) =>
model L := ln|pi-hat(a)| ~ Normal(mu, sigma^2), sigma^2 = c*k.

OVERDETERMINED CHECK (2 params, 3 moments): with L~N(mu,sigma^2),
  ln(gm)=mu ;  ln E|pi|^2 = 2mu+2sigma^2 ;  ln E|pi|^4 = 4mu+8sigma^2.
Fit (mu,sigma^2) from the DIRECT mean+variance of L, then PREDICT E2,E4 and compare to measured.
If E4 lands => 2-parameter closed form for the magnitude distribution (mu arithmetic, sigma^2 linear in k).
If it misses => log-normality wrong, the log-variance is just a fact.
Also: skew/excess-kurtosis of L (Gaussian => 0); is Sum_a|pi-hat|^4 = U^2 = 0.29754 (channel/rho object -- CAUTION,
pi-hat additive vs rho multiplicative; measure, don't assume).

All moments over the units U={3-nmid a}. Reuses fwd_hat. Cheap to k~13.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_singlerec import fwd_hat


def main():
    t0 = time.time()
    print("# PROBE LOGNORMAL -- L=ln|pi-hat| ~ N(mu,sigma^2), sigma^2~c*k ; predict E2,E4\n")
    print(f"   {'k':>2} {'mu':>8} {'sig2':>7} {'skew':>7} {'exkurt':>7} "
          f"{'E2rat':>7} {'E4rat':>7} {'Sum|pi|^4':>9}")
    rows = []
    for k in range(3, 14):
        N = 3 ** k
        ph, _ = fwd_hat(k)
        a = np.arange(N)
        U = a[a % 3 != 0]
        mag = np.abs(ph[U])
        L = np.log(mag)
        mu = L.mean(); s2 = L.var()
        sd = np.sqrt(s2)
        skew = ((L - mu) ** 3).mean() / sd ** 3
        exkurt = ((L - mu) ** 4).mean() / sd ** 4 - 3
        E2 = (mag ** 2).mean(); E4 = (mag ** 4).mean()
        # log-normal predictions from (mu, s2)
        E2pred = np.exp(2 * mu + 2 * s2)
        E4pred = np.exp(4 * mu + 8 * s2)
        sum4 = (mag ** 4).sum()
        rows.append((k, mu, s2, E2, E4, E2pred, E4pred))
        print(f"   {k:>2} {mu:>8.4f} {s2:>7.4f} {skew:>7.3f} {exkurt:>7.3f} "
              f"{E2/E2pred:>7.4f} {E4/E4pred:>7.4f} {sum4:>9.5f}")
    ks = np.array([r[0] for r in rows]); s2s = np.array([r[2] for r in rows])
    c, b = np.polyfit(ks, s2s, 1)
    r2 = 1 - np.sum((s2s - (b + c * ks)) ** 2) / np.sum((s2s - s2s.mean()) ** 2)
    print(f"\n   sigma^2 linear fit: sigma^2 = {c:.4f}*k + {b:.4f}  (R^2={r2:.5f})  [c*k random-walk model]")
    print("   [E2rat,E4rat -> 1 iff log-normal; E4 is the overdetermined test (fit uses only mu,sigma^2).")
    print("    skew/exkurt -> 0 iff Gaussian. Sum|pi|^4 vs U^2=0.29754 = pi-hat-ell4 vs rho-channel object.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
