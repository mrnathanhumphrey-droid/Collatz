"""
probe_c5_highprec_2026_05_30.py

Compute c(5) at mpmath dps=30 via offset_distribution at n=6, N=17^6=24,137,569.

Estimated runtime: 4-8 hours.
Saves c(5) value to a file as a safety belt.
"""
from __future__ import annotations
import sys, gc, time
from mpmath import mp, mpf
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30
q = 17
n = 6
A_MAX = 110  # 110 * log2 ≈ 33 digits; cushion above dps=30

def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

chi_table = [legendre(x, q) for x in range(q)]

def offset_distribution_mp(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    P_U = [mpf(0)] * N
    p_val = inv2
    half = mpf(1) / mpf(2)
    cur_w = half
    for a in range(1, A_MAX + 1):
        P_U[p_val] += cur_w
        p_val = (p_val * inv2) % N
        cur_w *= half
    total = sum(P_U[i] for i in range(N))
    P_U = [x / total for x in P_U]

    u_sup = [i for i in range(N) if P_U[i] != 0]
    u_wt  = [P_U[i] for i in u_sup]
    print(f"    |u_sup| = {len(u_sup)}", flush=True)

    P_S = P_U
    for j in range(n - 1, 0, -1):
        t_j = time.time()
        P_V = [mpf(0)] * N
        for i in range(N):
            P_V[(1 + q * i) % N] += P_S[i]
        P_new = [mpf(0)] * N
        for u, w in zip(u_sup, u_wt):
            for i in range(N):
                P_new[(u * i) % N] += w * P_V[i]
        P_S = P_new
        del P_V
        gc.collect()
        print(f"    j={j} done in {time.time()-t_j:.1f}s  (elapsed total: {time.time()-t0:.1f}s)", flush=True)
    return P_S

print(f"Computing c(5) at n={n}, N={q**n}, A_MAX={A_MAX}, dps={mp.dps}")
t0 = time.time()
P_X = offset_distribution_mp(q, n, A_MAX)
print(f"P_X done in {time.time()-t0:.1f}s")

N = q ** n
qm = q ** 5
t1 = time.time()
num = mpf(0); den = mpf(0)
for j in range(1, q):
    d = j * qm
    Pd = mpf(0)
    for y in range(N):
        Pd += P_X[y] * P_X[(y - d) % N]
    num += Pd * mpf(chi_table[j])
    den += Pd
c5 = num / den
print(f"P_D extraction in {time.time()-t1:.1f}s")
print(f"\n*** c(5) = {c5}")

with open("c5_highprec_result.txt", "w") as f:
    f.write(f"c(5) at dps={mp.dps}, n={n}, A_MAX={A_MAX}\n")
    f.write(f"{c5}\n")
print("Saved to c5_highprec_result.txt")
