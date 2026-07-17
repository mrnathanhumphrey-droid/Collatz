"""
probe_coset_collision_op_2026_05_29.py

Decompose A_chi/A_triv into intrinsic objects via the collision-transfer-operator formula:

  A_chi / A_triv  ->  g(chi) * c_inf * (1 - lambda) / (q*lambda - 1)

where
  lambda  = limiting q-adic collision rate, lim Z_{n+1}/Z_n  for Z_n = sum_y P_n(y)^2
          = lim (1/q) * (sum_xi |mu_hat_n(xi)|^2) / level-(n-1) analog
  c_inf   = limiting chi-bar moment of (X-Y)/q^{n-1} mod q at deep collision
          = lim E[chi-bar((X-Y)/q^{n-1}) | v_q(X-Y) = n-1]
  g(chi)  = sqrt(q) for chi=Legendre mod 17 (q=17, conductor q)

Strategy: compute lambda(n), c(n), and asym(n) at n=4,5,6 from FFT and check
- consistency (asym(n) ~ g(chi) c(n) (1-lambda(n))/(q lambda(n)-1))
- whether c_inf is closer to a clean closed form than the raw asym
- Aitken/Shanks acceleration on each independently
- PSLQ at improved precision
"""
from __future__ import annotations
import sys, gc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
from probe_coset_closedform_2026_05_29 import offset_distribution

def legendre(x, q):
    x %= q
    return 0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)

q = 17
chi_q = np.array([legendre(x, q) for x in range(q)], dtype=np.float64)
g_chi = float(np.sqrt(q))   # Legendre mod 17 has Gauss sum +sqrt(17) (17 = 1 mod 4)

results = {}
for n in (3, 4, 5, 6):
    N = q ** n
    print(f"  computing n={n}, N={N}...", flush=True)
    P = offset_distribution(q, n, 100)
    mu = np.fft.fft(P)
    a2 = np.abs(mu) ** 2
    del mu
    xi = np.arange(N)
    units_mask = xi % q != 0
    # asymmetry A_chi / A_triv
    chi_vals_full = chi_q[xi % q]
    A_chi = float(np.sum(chi_vals_full * a2))
    A_triv = float(np.sum(a2[units_mask]))
    asym = A_chi / A_triv
    # energy sum Z_n via Parseval: sum_xi |mu|^2 = N * sum_y P^2; so sum P^2 = (sum |mu|^2)/N
    sumP2 = float(np.sum(a2)) / N
    # the conductor-q char-fn restricted to v_q=n-1 frequencies:
    # frequencies xi with v_q(xi) = 0 (units) carry the conductor-q tau_chi; this IS what asym uses.
    # the "intrinsic" c(n) = (asym(n) * A_triv) / (q^{n-1} * g(chi) * w(n-1))
    # where w(n-1) = P(v_q(X-Y)=n-1). Compute w(n-1) from P_D = autocorrelation.
    # P_D(d) = (1/N) IFFT(|mu|^2) at d
    PD = np.real(np.fft.ifft(a2))  # circular autocorrelation; sum_y P(y)P(y-d mod N)
    # v_q(d) for each d
    # P(v_q(D) = m) = sum_{d: v_q(d)=m} PD(d), for m in 0..n
    vq = np.zeros(N, dtype=np.int64)
    for m in range(1, n + 1):
        vq[(np.arange(N) % (q ** m) == 0)] = m
    # mass per depth
    w = {}
    for m in range(n + 1):
        w[m] = float(PD[vq == m].sum())
    # c(m) = E[chi-bar((X-Y)/q^m mod q) | v_q(D)=m]
    c = {}
    for m in range(n):
        mask = (vq == m)
        d_idx = np.nonzero(mask)[0]
        # (d / q^m) mod q for d_idx
        leading_digit = (d_idx // (q ** m)) % q
        chi_bar_vals = chi_q[leading_digit]
        num = float(np.sum(PD[mask] * chi_bar_vals))
        c[m] = num / w[m] if w[m] != 0 else float('nan')
    results[n] = dict(asym=asym, sumP2=sumP2, w=w, c=c, N=N)
    del a2, PD, P, vq
    gc.collect()

# table
print("\nn | asym(n)       sumP2(n)        w(n-1)         c(n-1)")
prev = None
for n in sorted(results):
    r = results[n]
    print(f"{n} | {r['asym']:.12f}  {r['sumP2']:.6e}   {r['w'][n-1]:.6e}   {r['c'][n-1]:.12f}")
# lambda(n) = sumP2(n) / sumP2(n-1)  ... but only consecutive levels comparable
print("\nlambda(n) = sumP2(n) / sumP2(n-1):")
ns = sorted(results)
for i in range(1, len(ns)):
    a, b = ns[i-1], ns[i]
    lam = results[b]['sumP2'] / results[a]['sumP2']
    print(f"  lambda from n={a}->{b}: {lam:.10f}")
# c(n-1) sequence and Shanks
print("\nc(n-1) sequence (deep-collision chi-bar moment):")
cs = [results[n]['c'][n-1] for n in ns]
for n, cv in zip(ns, cs):
    print(f"  n={n}: c(n-1)={cv:.12f}")

def shanks(a, b, c):
    den = c - 2*b + a
    return c - (c - b)**2 / den if abs(den) > 1e-20 else c

if len(cs) >= 3:
    print("\nShanks acceleration on c sequence:")
    for i in range(len(cs) - 2):
        s = shanks(cs[i], cs[i+1], cs[i+2])
        print(f"  Shanks(c[{ns[i]}..{ns[i+2]}]) = {s:.12f}")

# verify consistency: asym ~ g(chi) * c * (1-lambda)/(q*lambda - 1)
print("\nConsistency check (each n): asym vs g(chi)*c(n-1)*(1-lambda)/(q*lambda-1)")
for i in range(1, len(ns)):
    a, b = ns[i-1], ns[i]
    lam = results[b]['sumP2'] / results[a]['sumP2']
    c_val = results[b]['c'][b-1]
    pred = g_chi * c_val * (1 - lam) / (q * lam - 1)
    print(f"  n={b}: lambda={lam:.6f}  c={c_val:.10f}  pred={pred:.10f}  asym={results[b]['asym']:.10f}  diff={pred-results[b]['asym']:.2e}")

# Aitken/Shanks on asym
print("\nShanks on asym sequence:")
asyms = [results[n]['asym'] for n in ns]
for i in range(len(asyms) - 2):
    s = shanks(asyms[i], asyms[i+1], asyms[i+2])
    print(f"  Shanks(asym[{ns[i]}..{ns[i+2]}]) = {s:.12f}")

# PSLQ at improved precision
try:
    from mpmath import mp, mpf, pslq, identify, sqrt
    mp.dps = 30
    # use the best Shanks-accelerated asym
    if len(asyms) >= 3:
        best_asym = shanks(asyms[-3], asyms[-2], asyms[-1])
        best_c = shanks(cs[-3], cs[-2], cs[-1])
        # one more Shanks (Aitken iterated)
        if len(asyms) >= 5:
            s1 = [shanks(asyms[i], asyms[i+1], asyms[i+2]) for i in range(len(asyms)-2)]
            if len(s1) >= 3:
                best_asym = shanks(s1[-3], s1[-2], s1[-1])
            s1c = [shanks(cs[i], cs[i+1], cs[i+2]) for i in range(len(cs)-2)]
            if len(s1c) >= 3:
                best_c = shanks(s1c[-3], s1c[-2], s1c[-1])
        print(f"\nBest extrapolated asym = {best_asym:.12f}")
        print(f"Best extrapolated c    = {best_c:.12f}")
        x = mpf(repr(best_asym))
        xc = mpf(repr(best_c))
        print(f"  PSLQ asym  [1,x,x^2,x^3]   = {pslq([mpf(1), x, x**2, x**3], maxcoeff=10**6)}")
        print(f"  PSLQ asym  [1,x,sqrt17,x^2]= {pslq([mpf(1), x, sqrt(17), x**2], maxcoeff=10**5)}")
        print(f"  identify(asym) = {identify(x)}")
        print(f"  PSLQ c     [1,xc,xc^2,xc^3]= {pslq([mpf(1), xc, xc**2, xc**3], maxcoeff=10**6)}")
        print(f"  identify(c)    = {identify(xc)}")
except Exception as e:
    print(f"  [PSLQ step skipped: {e}]")
