"""
probe_test_A1a_q11_2026_05_30.py

T-A1a: For q=11 (≡ 3 mod 8, so Legendre(2/11) = -1, 2 is NQR mod 11),
the Prediction A claim is:
  - c_∞ = 0 exactly (Haar average of Legendre on (Z/11)*)
  - ρ = 1/3 (eigenvalue of dominant chain on χ_2 mode)

Sanity check: ord_11(2) = 10, so ⟨2⟩ = (Z/11)* (2 is primitive root).
But Legendre(2/11) = -1 since 11 ≡ 3 mod 8.
"""
from __future__ import annotations
import sys, gc, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 11

# Verify Legendre(2/q)
def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

print(f"q = {q}, q mod 8 = {q % 8}")
print(f"Legendre(2/{q}) = {legendre(2, q)} (predicted -1)")
ord_2 = 1; x = 2 % q
while x != 1: ord_2 += 1; x = (x * 2) % q
print(f"ord_{q}(2) = {ord_2}, index in (Z/{q})* = {(q-1)//ord_2}")
print(f"Legendre table for {q}: {[legendre(i, q) for i in range(q)]}")

chi_q = np.array([legendre(x, q) for x in range(q)], dtype=np.float64)

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

def c_at_n_m(n, m, A_MAX=200):
    N = q ** n
    P = offset_distribution(q, n, A_MAX)
    mu = np.fft.fft(P)
    a2 = np.abs(mu) ** 2
    PD = np.real(np.fft.ifft(a2))
    qm = q ** m
    qmp1 = q ** (m + 1)
    d_vals = np.arange(N)
    mask = (d_vals % qm == 0) & (d_vals % qmp1 != 0)
    leading_digit = (d_vals[mask] // qm) % q
    chi_vals = chi_q[leading_digit]
    num = (PD[mask] * chi_vals).sum()
    den = PD[mask].sum()
    return float(num / den)

print(f"\nComputing c(m) at q={q} for m=0..6...")
import time
results = {}
for m in range(0, 7):
    n = max(m + 1, 2)
    t0 = time.time()
    print(f"  m={m}, n={n}, N={q**n}...", flush=True)
    try:
        results[m] = c_at_n_m(n=n, m=m, A_MAX=200)
        print(f"    c({m}) = {results[m]:+.16f}  ({time.time()-t0:.1f}s)")
    except MemoryError:
        print(f"    OOM at n={n}, skipping")
        break
    gc.collect()

# Δ_m (relative to c(0))
c0 = results.get(0, 0)
print(f"\n=== Results for q={q} ===")
print(f"c(0) = {c0:+.12f}")
print(f"\nm  c(m)            c(m)-c(0)")
for m in sorted(results):
    print(f"{m}  {results[m]:+.12f}  {results[m]-c0:+.6e}")

# Check Prediction A1
print(f"\n=== Test T-A1a ===")
print(f"PREDICTION: c_∞ = 0 (Haar Legendre = 0 on (Z/{q})*)")
print(f"PREDICTION: ρ = 1/3 (dominant chain eigenvalue on χ_2 mode)")

# Did c(m) decay toward 0?
if 5 in results:
    cm = [results[m] for m in range(6) if m in results]
    print(f"\nc(m) sequence: {cm}")
    print(f"|c(5)| / |c(0)| = {abs(cm[5])/max(abs(cm[0]), 1e-20):.6e}")
    print(f"Predicted (1/3)^5 = {(1/3)**5:.6e}")
    # If c_∞ = 0 with ratio 1/3, expect |c(m+1)/c(m)| → 1/3.
    print(f"\nRatios |c(m+1)/c(m)|:")
    for m in range(5):
        if m+1 in results and abs(results[m]) > 1e-15:
            r = results[m+1] / results[m]
            print(f"  c({m+1})/c({m}) = {r:+.6f}")

# Verdict
if 5 in results:
    if abs(results[5]) < 1e-3 * abs(results[0]):
        print(f"\n*** T-A1a PASS: c(m) decaying toward 0 ***")
    else:
        print(f"\n*** T-A1a FAIL: c(m) does not decay to 0 ***")
