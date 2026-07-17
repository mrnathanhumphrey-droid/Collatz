"""
probe_c4_highprec_2026_05_30.py

Compute c(4) at mpmath dps=40 via offset_distribution at n=5, N=17^5=1419857.
Single deliverable: c(4) as a high-precision mpf, to feed into the damped-osc recurrence.

Expected runtime: 1-2 hours on this box.
"""
from __future__ import annotations
import sys, gc, time
from mpmath import mp, mpf
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 50  # extra headroom; effective precision ~ 45 digits
q = 17
n = 5
A_MAX = 170

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
    # normalize
    total = sum(P_U[i] for i in range(N))
    P_U = [x / total for x in P_U]

    u_sup = [i for i in range(N) if P_U[i] != 0]
    u_wt  = [P_U[i] for i in u_sup]
    print(f"    |u_sup| = {len(u_sup)}", flush=True)

    P_S = P_U
    for j in range(n - 1, 0, -1):
        t_j = time.time()
        # P_V[(1 + q*i) mod N] = P_S[i]
        P_V = [mpf(0)] * N
        for i in range(N):
            P_V[(1 + q * i) % N] += P_S[i]
        # P_new[(u*i) mod N] += w * P_V[i]
        P_new = [mpf(0)] * N
        for u, w in zip(u_sup, u_wt):
            for i in range(N):
                P_new[(u * i) % N] += w * P_V[i]
        P_S = P_new
        del P_V
        gc.collect()
        print(f"    j={j} done in {time.time()-t_j:.1f}s  (elapsed total: {time.time()-t0:.1f}s)", flush=True)
    return P_S

print(f"Computing c(4) at n={n}, N={q**n}, A_MAX={A_MAX}, dps={mp.dps}")
t0 = time.time()
P_X = offset_distribution_mp(q, n, A_MAX)
print(f"P_X done in {time.time()-t0:.1f}s; ||P||_1 = {sum(P_X[i] for i in range(q**n))}")

# Extract c(4) via direct autocorrelation at d = j * q^4
N = q ** n
qm = q ** 4
t1 = time.time()
num = mpf(0); den = mpf(0)
for j in range(1, q):
    d = j * qm
    Pd = mpf(0)
    for y in range(N):
        Pd += P_X[y] * P_X[(y - d) % N]
    num += Pd * mpf(chi_table[j])
    den += Pd
c4 = num / den
print(f"P_D extraction in {time.time()-t1:.1f}s")
print(f"\n*** c(4) = {c4}")

# Also save to a results file for safety
with open("c4_highprec_result.txt", "w") as f:
    f.write(f"c(4) at dps={mp.dps}, n={n}, A_MAX={A_MAX}\n")
    f.write(f"{c4}\n")
print("Saved to c4_highprec_result.txt")
