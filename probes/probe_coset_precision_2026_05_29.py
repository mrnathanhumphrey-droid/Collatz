"""
probe_coset_precision_2026_05_29.py
High-precision push on the q=17 coset asymmetry + Richardson extrapolation + PSLQ/algebraic test.
asym(n) = A_chi/A_triv with chi=Legendre mod 17.  Goal: enough digits to accept/reject a closed form.
"""
from __future__ import annotations
import sys, gc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
from probe_coset_closedform_2026_05_29 import offset_distribution

def legendre(x, q):
    x %= q
    return 0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)

def asym(q, n, chi_q):
    N = q ** n
    P = offset_distribution(q, n, 100)
    mu = np.fft.fft(P); a2 = np.abs(mu) ** 2
    del mu, P; gc.collect()
    xi = np.arange(N); units = xi % q != 0
    chi_arr = np.zeros(N)
    cl = np.array([chi_q[r] for r in range(q)])
    chi_arr = cl[xi % q]
    A_chi = float(np.sum(chi_arr * a2))
    A_triv = float(np.sum(a2[units]))
    return A_chi / A_triv

def main():
    q = 17
    chi_q = {x: legendre(x, q) for x in range(q)}
    ns = [3, 4, 5]
    vals = {}
    for n in ns:
        vals[n] = asym(q, n, chi_q)
        print(f"  n={n}: asym = {vals[n]:.12f}")
    try:
        vals[6] = asym(q, 6, chi_q); ns.append(6)
        print(f"  n=6: asym = {vals[6]:.12f}")
    except MemoryError:
        print("  n=6: MemoryError (skipped)")

    # Richardson: asym(n) = L + C r^n. Use last three.
    last = ns[-3:]
    a, b, c = vals[last[0]], vals[last[1]], vals[last[2]]
    d1, d2 = b - a, c - b
    r = d2 / d1 if d1 != 0 else float('nan')
    L = c + d2 * r / (1 - r) if abs(1 - r) > 1e-12 else c
    print(f"\n  Richardson (n={last}): ratio r={r:.5f}, extrapolated limit L = {L:.12f}")

    try:
        from mpmath import mp, mpf, pslq, identify, sqrt, log
        mp.dps = 25
        x = mpf(repr(L))
        print(f"\n  CF of L: ", end="")
        t = x; cf = []
        for _ in range(15):
            ai = int(t); cf.append(ai); fr = t - ai
            if fr == 0: break
            t = 1 / fr
        print(cf)
        print(f"  identify(L) = {identify(x)}")
        # low-degree algebraic: pslq([1, x, x^2]) and ([1,x,x^2,x^3])
        print(f"  pslq[1,x,x^2]   = {pslq([mpf(1), x, x**2], maxcoeff=10**6)}")
        print(f"  pslq[1,x,x^2,x^3] = {pslq([mpf(1), x, x**2, x**3], maxcoeff=10**5)}")
        # against natural constants
        V = mpf(76)/765
        print(f"  V=76/765={float(V):.12f}; pslq[1,x,V] = {pslq([mpf(1), x, V], maxcoeff=10**6)}")
        print(f"  pslq[1,x,sqrt17] = {pslq([mpf(1), x, sqrt(17)], maxcoeff=10**6)}")
        print(f"  pslq[1,x,V,sqrt17/q] = {pslq([mpf(1), x, V, sqrt(17)/17], maxcoeff=10**5)}")
    except Exception as e:
        print(f"  [mpmath step skipped: {e}]")

if __name__ == "__main__":
    main()
