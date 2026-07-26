"""
PROBE BRIDGE -- Wilson's 2 gates (2026-07-26): (1) additive recursion on the FORWARD measure; (2) the Gauss-sum bridge.

Resolution: build_nu = perpetuity X (intercept 1); Wilson's W_n=2^-a(3W+1) = forward Syracuse (intercept 2^-a), and
W = 2^-a X (the one-level offset). So test the recursion on the FORWARD chain measure = stationary_trunc(3,k) (=R66).

STEP 1: mu-hat_n(xi) =?= Sum_a 2^-a e(xi 2^-a/3^n) mu-hat_{n-1}(xi 2^-a mod 3^{n-1}),  mu-hat = fft(forward measure).
STEP 2 (the one Wilson cares about -- BRIDGE): rho-hat(a) = nu-hat(chi_a) = (1/tau(chi_a-bar)) Sum_t chi_a-bar(t) mu-hat(t),
   |tau|=sqrt(q). Pure Gauss-sum identity between the multiplicative (channel) and additive (dynamics) transforms.
   Verify to machine precision on the forward measure (holds for ANY measure if chi primitive & tau right).

Reuses probe_27 stationary_trunc + R10.dlog_table. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import gcd
import numpy as np
import probe_charledger_R10 as R10
from probe_27_high_k_rho_q5 import stationary_trunc

R66_MAX = {2: 0.14283, 3: 0.06360, 4: 0.03130, 5: 0.01670, 6: 0.00924}


def dense_fwd(k):
    N = 3 ** k
    pi, n = stationary_trunc(3, k)
    cp = np.array([r for r in range(N) if gcd(r, 3) == 1], dtype=np.int64)
    dense = np.zeros(N); dense[cp] = pi
    return dense / dense.sum()


def main():
    t0 = time.time()
    print("# PROBE BRIDGE -- additive recursion on forward measure + Gauss-sum bridge\n")

    fwd = {k: dense_fwd(k) for k in range(2, 7)}
    mh = {k: np.fft.fft(fwd[k]) for k in range(2, 7)}

    print("## sanity: max|fft(forward)|^2 over primitive vs R66 max")
    for k in range(2, 7):
        N = 3 ** k; a = np.arange(N); prim = a % 3 != 0
        v = (np.abs(mh[k]) ** 2)[prim].max()
        print(f"   k={k}: {v:.5f}  vs R66 {R66_MAX[k]:.5f}  {'MATCH' if abs(v-R66_MAX[k])<1e-4 else 'no'}")
    print()

    print("## STEP 1: additive recursion on FORWARD measure  mu_n(xi)=Sum_a 2^-a e(xi 2^-a/3^n) mu_{n-1}(xi 2^-a mod 3^{n-1})")
    for n in range(3, 7):
        N = 3 ** n; Nm = 3 ** (n - 1)
        inv2 = pow(2, -1, N)
        prim = [xi for xi in range(1, N) if xi % 3 != 0]
        errs = []
        for xi in prim:
            rhs = 0j; p = 1
            for a in range(1, 60):
                p = (p * inv2) % N
                wa = (xi * p) % N
                rhs += (0.5 ** a) * np.exp(2j * np.pi * wa / N) * mh[n - 1][wa % Nm]
            errs.append(abs(rhs - mh[n][xi]))
        rel = np.mean(errs) / np.mean([abs(mh[n][xi]) for xi in prim])
        print(f"   n={n}: rel = {rel:.3e}  [{'REPRODUCES' if rel<1e-9 else ('close' if rel<1e-3 else 'NO')}]")
    print()

    print("## STEP 2: GAUSS-SUM BRIDGE  rho-hat(a) = (1/tau(chi_a-bar)) Sum_t chi_a-bar(t) mu-hat(t)")
    for k in (3, 4, 5):
        N = 3 ** k
        nu = fwd[k]
        d = R10.dlog_table(k)
        # chi_a(x) = e(2pi i a * dlog(x)/N), dlog(x)=d[(x-1)//3] for x coprime to 3; 0 else
        x = np.arange(N)
        cop = (x % 3 != 0)
        dl = np.zeros(N, dtype=np.int64)
        dl[cop] = np.array([d[(xx - 1) // 3 % N] for xx in x[cop]])
        muhat = np.fft.fft(nu)                                   # additive
        # test a primitive (a coprime to 3)
        maxrel = 0.0
        for a in [1, 2, 4, 5, 7]:
            if a % 3 == 0:
                continue
            chi_a = np.where(cop, np.exp(2j * np.pi * a * dl / N), 0.0)      # chi_a(x) over x
            rho_hat_a = np.sum(nu * chi_a)                                    # = nu-hat(chi_a) = fft(rho)(a) analog
            # tau(chi_a-bar) = Sum_t chi_a-bar(t) e(2pi i t/N)
            chibar = np.where(cop, np.exp(-2j * np.pi * a * dl / N), 0.0)
            tau = np.sum(chibar * np.exp(2j * np.pi * x / N))
            rhs = (1.0 / tau) * np.sum(chibar * muhat)
            rel = abs(rhs - rho_hat_a) / (abs(rho_hat_a) + 1e-30)
            maxrel = max(maxrel, rel)
        # |tau| = sqrt(q)?
        chibar1 = np.where(cop, np.exp(-2j * np.pi * 1 * dl / N), 0.0)
        tau1 = np.sum(chibar1 * np.exp(2j * np.pi * x / N))
        print(f"   k={k}: bridge max rel (a=1,2,4,5,7) = {maxrel:.2e} [{'HOLDS' if maxrel<1e-9 else 'NO'}] ; "
              f"|tau|={abs(tau1):.4f} vs sqrt(3^k)={np.sqrt(N):.4f}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
