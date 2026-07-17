"""
probe_delta_pslq_2026_05_30.py

Compute c(2), c(3), c(4) at maximum float64 precision via the existing FFT machinery
(A_MAX = 200 ensures geometric-tail error << machine epsilon at every n).
Then Δ_m = c(m) - 19/127 and PSLQ against:
  - {1, log2, log(2), √2, √17, √(17-1), e, π}
  - Li_s(2^{-8}) = polylog at z=1/256, s=1..6
  - Lerch transcendents Φ(1/256, s, a) for s=1..3, a=1/8 .. 8/8

Also include Δ_1 EXACTLY (from c(1) exact rational - 19/127) at high precision.
"""
from __future__ import annotations
import sys, gc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 17

def offset_distribution(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a)
        p = (p * inv2) % N
    P_U /= P_U.sum()
    u_support = np.nonzero(P_U)[0]; u_weight = P_U[u_support]
    P_S = P_U.copy()
    for j in range(n - 1, 0, -1):
        v_idx = (1 + q * arange) % N
        P_V = np.zeros(N, dtype=np.float64)
        np.add.at(P_V, v_idx, P_S)
        P_new = np.zeros(N, dtype=np.float64)
        for u, w in zip(u_support, u_weight):
            P_new[(u * arange) % N] += w * P_V
        P_S = P_new
    return P_S

def chi(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq - 1) // 2, qq) == 1 else -1)

chi_q = np.array([chi(x, q) for x in range(q)], dtype=np.float64)

def c_at_n_m(n, m, A_MAX=200):
    """Compute c(m) = E[chi-bar((D)/q^m mod q) | v_q(D) = m] at FFT level n."""
    N = q ** n
    P = offset_distribution(q, n, A_MAX)
    # autocorrelation P_D via FFT
    mu = np.fft.fft(P)
    a2 = np.abs(mu) ** 2
    PD = np.real(np.fft.ifft(a2))
    # v_q for each d
    qm = q ** m
    qmp1 = q ** (m + 1)
    # we want d such that d % q^m == 0 and d % q^(m+1) != 0
    d_vals = np.arange(N)
    mask = (d_vals % qm == 0) & (d_vals % qmp1 != 0)
    leading_digit = (d_vals[mask] // qm) % q  # j = 1..q-1
    chi_vals = chi_q[leading_digit]
    num = (PD[mask] * chi_vals).sum()
    den = PD[mask].sum()
    return float(num / den)

print("Computing c(m) at high A_MAX=200, FFT level n=m+1...")
A_MAX = 200
results = {}
for m in (1, 2, 3, 4):
    n = m + 1
    print(f"  m={m}, n={n}, N=q^{n}={q**n}...", flush=True)
    cm = c_at_n_m(n=n, m=m, A_MAX=A_MAX)
    results[m] = cm
    print(f"    c({m}) = {cm:.16f}")
    gc.collect()

# Use the EXACT c(1) from the rational machinery as the high-precision Δ_1 input
from fractions import Fraction
c1_exact = Fraction(265011804960406635465672455997699, 1730087916969634762193659498034425)
c0_exact = Fraction(19, 127)

delta1_exact = c1_exact - c0_exact
print(f"\n=== EXACT Δ_1 = c(1) - 19/127 ===")
print(f"Δ_1 numerator   = {delta1_exact.numerator}")
print(f"Δ_1 denominator = {delta1_exact.denominator}")
print(f"Δ_1 (float)     = {float(delta1_exact):.16e}")

# Δ_m for m=2,3,4 from FFT (float64 precision)
print(f"\n=== Δ_m sequence ===")
c0_float = 19 / 127
delta = {}
delta[1] = float(delta1_exact)
for m in (2, 3, 4):
    delta[m] = results[m] - c0_float
print(f"  {'m':>2} {'c(m)':>20} {'Δ_m':>22}")
for m in (1, 2, 3, 4):
    cm = float(c1_exact) if m == 1 else results[m]
    print(f"  {m:>2} {cm:>20.15f} {delta[m]:>+22.15e}")

# Differences and ratios
print(f"\n  Δ_m / Δ_1:")
for m in (1, 2, 3, 4):
    print(f"    m={m}: {delta[m] / delta[1]:.12f}")
print(f"\n  c(m+1) - c(m) (alternating?):")
prev_cm = None
for m in (1, 2, 3, 4):
    cm = float(c1_exact) if m == 1 else results[m]
    if prev_cm is not None:
        print(f"    c({m}) - c({m-1}) = {cm - prev_cm:+.6e}")
    prev_cm = cm

# Now PSLQ on Δ_m
print(f"\n=== PSLQ against polylog / Lerch candidates ===")
try:
    from mpmath import mp, mpf, mpc, pslq, identify, polylog, sqrt, log, pi, e
    mp.dps = 30

    # Δ_1 as exact rational at high precision
    d1 = mpf(delta1_exact.numerator) / mpf(delta1_exact.denominator)
    d2 = mpf(delta[2])
    d3 = mpf(delta[3])
    d4 = mpf(delta[4])

    # polylogs Li_s(1/256) for s=1..6
    z = mpf(1) / mpf(256)
    polylogs = [polylog(s, z) for s in range(1, 7)]

    # Lerch transcendents Φ(z, s, a) for s=1..3, a=k/8 for k=1..8
    # Lerch Φ(z, s, a) = sum_{n>=0} z^n / (n+a)^s
    # mpmath has lerchphi(z, s, a)
    from mpmath import lerchphi
    lerchs = []
    lerch_names = []
    for s in (1, 2, 3):
        for k in range(1, 9):
            a = mpf(k) / mpf(8)
            lerchs.append(lerchphi(z, s, a))
            lerch_names.append(f"Φ(1/256, {s}, {k}/8)")

    print(f"\n--- Δ_1 against simple bases ---")
    basis_simple = [mpf(1), log(mpf(2)), log(mpf(17)), sqrt(mpf(17)), pi, e]
    names_simple = ["1", "log(2)", "log(17)", "√17", "π", "e"]
    for cand_name, cand_val in zip(["Δ_1", "Δ_2", "Δ_3", "Δ_4"], [d1, d2, d3, d4]):
        relation = pslq([cand_val] + basis_simple, maxcoeff=10**8, tol=mpf(10)**(-22))
        print(f"  {cand_name}: {relation}")

    print(f"\n--- Δ_m against polylogs Li_s(1/256) ---")
    for cand_name, cand_val in zip(["Δ_1", "Δ_2", "Δ_3", "Δ_4"], [d1, d2, d3, d4]):
        relation = pslq([cand_val, mpf(1)] + polylogs, maxcoeff=10**8, tol=mpf(10)**(-22))
        print(f"  {cand_name}: {relation}")

    print(f"\n--- Δ_m against Lerch transcendents ---")
    for cand_name, cand_val in zip(["Δ_1", "Δ_2", "Δ_3", "Δ_4"], [d1, d2, d3, d4]):
        relation = pslq([cand_val, mpf(1)] + lerchs, maxcoeff=10**6, tol=mpf(10)**(-15))
        print(f"  {cand_name}: {relation}")

    # also try identify on each
    print(f"\n--- mpmath.identify on each Δ_m ---")
    for cand_name, cand_val in zip(["Δ_1", "Δ_2", "Δ_3", "Δ_4"], [d1, d2, d3, d4]):
        print(f"  {cand_name} = {cand_val}")
        try:
            ident = identify(cand_val, ["log(2)", "log(17)", "sqrt(17)", "pi"])
            print(f"    identify: {ident}")
        except Exception as ex:
            print(f"    identify failed: {ex}")

    # Joint PSLQ: relations among the Δ_m's themselves
    print(f"\n--- PSLQ: relations AMONG (1, Δ_1, Δ_2, Δ_3, Δ_4) ---")
    rel = pslq([mpf(1), d1, d2, d3, d4], maxcoeff=10**10, tol=mpf(10)**(-12))
    print(f"  (1, Δ_1, Δ_2, Δ_3, Δ_4): {rel}")

    # Ratios
    print(f"\n--- PSLQ: (Δ_1, Δ_2, Δ_3, Δ_4) against {{1, log(2)}} ---")
    rel = pslq([d1, d2, d3, d4, mpf(1), log(mpf(2))], maxcoeff=10**8, tol=mpf(10)**(-12))
    print(f"  : {rel}")

except Exception as e:
    print(f"PSLQ section failed: {e}")
    import traceback; traceback.print_exc()
