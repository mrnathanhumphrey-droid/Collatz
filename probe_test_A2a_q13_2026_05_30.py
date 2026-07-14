"""
T-A2a' (REVISED): q=13, ≡ 5 mod 8, splits in Z[i] as 13 = (3+2i)(3-2i).

Predict:
- c(0) ≠ 0, a clean rational (analog of 19/127 at q=17)
- c(m) → 0 with damping ρ close to 1/3 (NOT 0.076)
- Eigenvalue direction tied to Gaussian prime (3+2i), NOT (1+2i)
"""
from __future__ import annotations
import sys, gc, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 13
def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

print(f"q = {q}, q mod 4 = {q % 4}, q mod 8 = {q % 8}")
print(f"Legendre(-1/{q}) = {legendre(-1, q)} (predicted +1 for q ≡ 1 mod 4)")
print(f"Legendre(2/{q})  = {legendre(2, q)} (predicted -1 for q ≡ 5 mod 8)")
ord_2 = 1; x = 2 % q
while x != 1: ord_2 += 1; x = (x * 2) % q
print(f"ord_{q}(2) = {ord_2}, index in (Z/{q})* = {(q-1)//ord_2}")
print(f"13 splits in Z[i] as (3+2i)(3-2i), |3+2i|² = 13 ✓")

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
    PD = np.real(np.fft.ifft(np.abs(mu)**2))
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
results = {}
for m in range(0, 7):
    n = max(m + 1, 2)
    t0 = time.time()
    try:
        print(f"  m={m}, n={n}, N={q**n}...", flush=True)
        results[m] = c_at_n_m(n=n, m=m, A_MAX=200)
        print(f"    c({m}) = {results[m]:+.16f}  ({time.time()-t0:.1f}s)")
    except (MemoryError, ValueError) as e:
        print(f"    fail: {e}")
        break
    gc.collect()

print(f"\n=== q=13 results ===")
for m in sorted(results):
    print(f"  c({m}) = {results[m]:+.16f}")

# Try to identify c(0) as a clean rational
from fractions import Fraction
c0 = results.get(0, 0)
print(f"\n=== c(0) candidate rationals (denom ≤ 5000) ===")
f = Fraction(c0).limit_denominator(5000)
print(f"  c(0) ≈ {f} = {float(f):.16f}  diff = {c0 - float(f):.2e}")

# Sequence ratios → expected close to -1/3
print(f"\n=== Damping ratios c(m+1)/c(m) ===")
for m in range(5):
    if m+1 in results and abs(results[m]) > 1e-10:
        r = results[m+1] / results[m]
        print(f"  c({m+1})/c({m}) = {r:+.6f}  (predicted -1/3 = -0.3333)")

# Recurrence fit if possible (5 points)
if 5 in results:
    deltas = list(results.values())
    if all(abs(d) > 1e-10 for d in deltas[:5]):
        try:
            M = np.array([
                [deltas[1], -deltas[0], 1.0],
                [deltas[2], -deltas[1], 1.0],
                [deltas[3], -deltas[2], 1.0],
            ])
            rhs = np.array([deltas[2], deltas[3], deltas[4]])
            u, v, K = np.linalg.solve(M, rhs)
            print(f"\n=== Recurrence fit on c(0..4) ===")
            print(f"  u = {u}, v = {v}, K = {K}")
            if v > 0:
                rho = np.sqrt(v)
                print(f"  ρ = {rho:.6f}  (predicted ≈ 1/3)")
                if abs(u / (2*rho)) <= 1:
                    theta = np.arccos(u / (2*rho))
                    print(f"  θ = {theta:.6f} rad = {np.degrees(theta):.4f}°  (predicted π for real-negative eigenvalue)")
        except Exception as e:
            print(f"  recurrence fit failed: {e}")

print(f"\n=== Verdict ===")
if 0 in results and abs(results[0]) > 1e-10:
    print(f"  c(0) ≠ 0 ✓ (= {results[0]:.6f})")
    if 4 in results and abs(results[4]) < abs(results[0]) / 50:
        print(f"  c(m) decaying to 0 ✓")
        print(f"*** T-A2a' PASS ***")
    else:
        print(f"  c(m) NOT decaying as expected ✗")
        print(f"*** T-A2a' FAIL ***")
else:
    print(f"  c(0) = 0 unexpectedly")
    print(f"*** T-A2a' FAIL (swap symmetry should NOT apply) ***")
