"""
T-A2c (NEW): q=29, ≡ 5 mod 8, splits in Z[i] as 29 = (5+2i)(5-2i).
Second confirmation of the q ≡ 5 mod 8 fast-decay prediction (after q=13).

Predict:
- c(0) ≠ 0 (clean rational)
- c(m+1)/c(m) → -3/5 (dominant Geom(4) Markov rate)
- c_∞ = 0
"""
from __future__ import annotations
import sys, gc, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 29
def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

print(f"q = {q}, q mod 4 = {q % 4}, q mod 8 = {q % 8}")
print(f"Legendre(-1/{q}) = {legendre(-1, q)} (predicted +1 for q ≡ 1 mod 4)")
print(f"Legendre(2/{q})  = {legendre(2, q)} (predicted -1 for q ≡ 5 mod 8)")
ord_2 = 1; x = 2 % q
while x != 1: ord_2 += 1; x = (x * 2) % q
print(f"ord_{q}(2) = {ord_2}, index in (Z/{q})* = {(q-1)//ord_2}")
print(f"29 = 25 + 4 = (5+2i)(5-2i) in Z[i]")

chi_q = np.array([legendre(x, q) for x in range(q)], dtype=np.float64)

def offset_distribution(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a); p = (p * inv2) % N
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
    qm = q ** m; qmp1 = q ** (m + 1)
    d_vals = np.arange(N)
    mask = (d_vals % qm == 0) & (d_vals % qmp1 != 0)
    leading_digit = (d_vals[mask] // qm) % q
    num = (PD[mask] * chi_q[leading_digit]).sum()
    den = PD[mask].sum()
    return float(num / den)

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
        print(f"    fail: {e}"); break
    gc.collect()

print(f"\n=== q={q} results ===")
for m in sorted(results):
    print(f"  c({m}) = {results[m]:+.16f}")

print(f"\n=== Damping ratios ===")
for m in range(4):
    if m+1 in results and abs(results[m]) > 1e-10:
        r = results[m+1] / results[m]
        print(f"  c({m+1})/c({m}) = {r:+.6f}  (predicted -3/5 = -0.6000)")

print(f"\n=== Verdict ===")
if 0 in results and abs(results[0]) > 1e-10:
    if 3 in results and abs(results[3] / results[0]) < 0.3:
        last_ratio = results[3] / results[2] if 2 in results else None
        if last_ratio and abs(last_ratio + 0.6) < 0.01:
            print(f"  *** T-A2c PASS: ρ ≈ 3/5 confirmed at q=29 ***")
        else:
            print(f"  *** T-A2c PARTIAL: c decays but ratio = {last_ratio}, not -3/5 ***")
    else:
        print(f"  *** T-A2c FAIL: c not decaying ***")
else:
    print(f"  *** T-A2c FAIL: c(0) = 0 unexpectedly ***")
