"""
probe_c_inf_depth_extrap_2026_05_29.py

Chase c_inf for q=17 (Legendre chi). At a single level-n FFT we have access to c(m) for ALL
depths m=0..n-1 (via autocorrelation P_D and the leading-q-digit moments). Use the full depth
sequence at n=6 (deepest tractable) to extrapolate m -> infinity via iterated Shanks/Aitken,
then PSLQ at improved precision.
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

n = 6
N = q ** n
print(f"computing q={q}, n={n}, N={N}...", flush=True)
P = offset_distribution(q, n, 100)
mu = np.fft.fft(P)
a2 = np.abs(mu) ** 2
del mu
print("  FFT done, computing autocorrelation P_D...", flush=True)
PD = np.real(np.fft.ifft(a2))
del a2, P; gc.collect()

# v_q(d) for each d in [0, N)
print("  computing v_q...", flush=True)
vq = np.zeros(N, dtype=np.int8)
for m in range(1, n + 1):
    vq[(np.arange(N) % (q ** m) == 0)] = m

# c(m) = E[chi_bar(d/q^m mod q) | v_q(d) = m]
print(f"\n{'m':>3} {'w(m)=P(v_q=m)':>14} {'c(m)':>16}")
cm = []
for m in range(n):
    mask = (vq == m)
    d_idx = np.nonzero(mask)[0]
    leading = (d_idx // (q ** m)) % q
    chi_bar_vals = chi_q[leading]
    w_m = float(PD[mask].sum())
    num = float(np.sum(PD[mask] * chi_bar_vals))
    c_m = num / w_m if w_m != 0 else float('nan')
    cm.append(c_m)
    print(f"{m:>3} {w_m:>14.6e} {c_m:>16.12f}")

def shanks(a, b, c):
    den = c - 2*b + a
    return c - (c - b)**2 / den if abs(den) > 1e-20 else c

# iterated Shanks
print("\nIterated Shanks:")
seq = cm[:]
level = 0
print(f"  level {level}: " + " ".join(f"{x:.10f}" for x in seq))
while len(seq) >= 3:
    new = [shanks(seq[i], seq[i+1], seq[i+2]) for i in range(len(seq) - 2)]
    level += 1
    print(f"  level {level}: " + " ".join(f"{x:.10f}" for x in new))
    seq = new

c_inf_est = seq[-1] if seq else cm[-1]
print(f"\nbest c_inf estimate = {c_inf_est:.12f}")

# PSLQ at improved precision
try:
    from mpmath import mp, mpf, pslq, identify, sqrt, pi, exp, log, e as mpe
    mp.dps = 30
    x = mpf(repr(c_inf_est))
    print(f"\nCF: ", end="")
    t = x; cf = []
    for _ in range(18):
        ai = int(t); cf.append(ai); fr = t - ai
        if fr == 0: break
        t = 1/fr
    print(cf)
    print(f"identify({c_inf_est:.10f}) = {identify(x)}")
    print(f"PSLQ[1,x,x^2,x^3]      = {pslq([mpf(1), x, x**2, x**3], maxcoeff=10**7)}")
    print(f"PSLQ[1,x,x^2,x^3,x^4]  = {pslq([mpf(1), x, x**2, x**3, x**4], maxcoeff=10**6)}")
    # try against natural constants from the project
    s17 = sqrt(17)
    V_lead = mpf(76)/765   # leading-order V
    # c_leading from formula: asym_lead = (sqrt17/7) c_lead -> c_lead = 7 V / sqrt17
    c_lead = 7 * V_lead / s17
    print(f"c_leading (from V=76/765) = {float(c_lead):.12f}; ratio c/c_lead = {float(x/c_lead):.10f}")
    print(f"PSLQ[1,x,c_lead]       = {pslq([mpf(1), x, c_lead], maxcoeff=10**6)}")
    print(f"PSLQ[1,x,sqrt17,1/s17] = {pslq([mpf(1), x, s17, 1/s17], maxcoeff=10**5)}")
    # log identities
    print(f"PSLQ[1,x,log(2),log(3)]= {pslq([mpf(1), x, log(2), log(3)], maxcoeff=10**5)}")
    print(f"PSLQ[1,x,log(17)]      = {pslq([mpf(1), x, log(17)], maxcoeff=10**6)}")
except Exception as e:
    print(f"[PSLQ skipped: {e}]")
