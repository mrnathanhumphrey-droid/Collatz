"""
T-A1c (REVISED): q=23, ≡ 3 mod 4 (inert in Z[i]).
Predict: c(m) ≡ 0 for all m by swap-symmetry + Fermat-Euler.
"""
from __future__ import annotations
import sys, gc, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 23
def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

print(f"q = {q}, q mod 4 = {q % 4}, q mod 8 = {q % 8}")
print(f"Legendre(-1/{q}) = {legendre(-1, q)} (predicted -1 for q ≡ 3 mod 4)")
print(f"Legendre(2/{q})  = {legendre(2, q)}")
ord_2 = 1; x = 2 % q
while x != 1: ord_2 += 1; x = (x * 2) % q
print(f"ord_{q}(2) = {ord_2}")

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

# Limit to manageable n: 23^5 = 6.4M, 23^6 = 148M (probably OK in 36GB RAM)
print(f"\nComputing c(m) at q={q} for m=0..4...")
results = {}
for m in range(0, 5):
    n = max(m + 1, 2)
    t0 = time.time()
    print(f"  m={m}, n={n}, N={q**n}...", flush=True)
    try:
        results[m] = c_at_n_m(n=n, m=m, A_MAX=200)
        print(f"    c({m}) = {results[m]:+.16f}  ({time.time()-t0:.1f}s)")
    except (MemoryError, ValueError) as e:
        print(f"    fail: {e}")
        break
    gc.collect()

print(f"\n=== T-A1c verdict ===")
print(f"PREDICTION: c(m) ≡ 0 for all m (q ≡ 3 mod 4, inert in Z[i])")
max_c = max(abs(v) for v in results.values()) if results else 0
print(f"max |c(m)| = {max_c:.6e}")
if max_c < 1e-10:
    print(f"*** T-A1c PASS: c(m) ≡ 0 to machine precision ***")
else:
    print(f"*** T-A1c FAIL: c(m) non-zero, swap-symmetry argument broken ***")
