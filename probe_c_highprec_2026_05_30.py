"""
probe_c_highprec_2026_05_30.py

Compute c(2), c(3), c(4) at high precision via mpmath using the digit-by-digit
offset-distribution construction (no FFT — exact rational walk with mpf weights).
Then form Δ_m = c(m) - 19/127 at ~30 digit precision and PSLQ.

We use mpf with dps=40 throughout. c(m) only depends on P_X mod q^{m+1},
so we run offset_distribution at n=m+1.

For n=5 (c(4)), length = 1419857 → tractable but heavy. We may skip if too slow
and use the FFT-double-precision c(4) ≈ 0.152988709047 as a backup.
"""
from __future__ import annotations
import sys, time, gc
from mpmath import mp, mpf, mpc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40

q = 17
DPS = 40

def legendre(x, q):
    x %= q
    return 0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)

chi_table = [legendre(x, q) for x in range(q)]

def offset_distribution_mp(q, n, A_MAX):
    """Compute P_X (Syracuse-distributed q-adic at depth n) as mpf array length q^n."""
    N = q ** n
    inv2 = pow(2, -1, N)
    # Initial: P_U[2^-a mod N] = 2^-a for a >= 1 up to A_MAX
    P_U = [mpf(0)] * N
    p_val = inv2
    half = mpf("0.5")
    cur_w = half
    for a in range(1, A_MAX + 1):
        P_U[p_val] += cur_w
        p_val = (p_val * inv2) % N
        cur_w *= half
    # normalize so sum = 1 (geometric tail correction)
    total = sum(P_U)
    for i in range(N): P_U[i] /= total

    # support indices for sparse u-loop
    u_sup = [i for i in range(N) if P_U[i] != 0]
    u_wt = [P_U[i] for i in u_sup]
    P_S = list(P_U)

    arange_mod_N = list(range(N))
    for j in range(n - 1, 0, -1):
        t0 = time.time()
        # P_V[(1 + q*i) mod N] = P_S[i]
        P_V = [mpf(0)] * N
        for i in range(N):
            v_idx = (1 + q * i) % N
            P_V[v_idx] += P_S[i]
        # P_new[(u*i) mod N] += w * P_V[i]
        P_new = [mpf(0)] * N
        for u, w in zip(u_sup, u_wt):
            for i in range(N):
                idx = (u * i) % N
                P_new[idx] += w * P_V[i]
        P_S = P_new
        gc.collect()
        print(f"    j={j} done in {time.time()-t0:.1f}s")
    return P_S

def compute_c_m(n, m, A_MAX=300):
    """Compute c(m) = E[chi-bar((X-Y)/q^m mod q) | v_q(X-Y) = m] at FFT-level n.
    Returns c(m) as mpf at current dps.
    """
    N = q ** n
    print(f"  computing P_X at n={n}, N={N}, A_MAX={A_MAX}...")
    t0 = time.time()
    P_X = offset_distribution_mp(q, n, A_MAX)
    print(f"  P_X done in {time.time()-t0:.1f}s; sum = {sum(P_X)}")

    # Compute P_D(d) for d = j*q^m, j=1..q-1 (the v_q(d)=m values).
    # P_D(d) = sum_y P_X(y) P_X((y - d) mod N).
    # We only need q-1 values, so direct.
    qm = q ** m
    num = mpf(0)
    den = mpf(0)
    for j in range(1, q):
        d = j * qm
        # P_D(d) = sum_y P_X(y) * P_X((y - d) % N)
        Pd = mpf(0)
        for y in range(N):
            yp = (y - d) % N
            Pd += P_X[y] * P_X[yp]
        chi_j = chi_table[j]  # leading digit (d/q^m) % q = j
        num += Pd * chi_j
        den += Pd
        print(f"    j={j}: P_D({d}) = {Pd}")
    c_m = num / den
    return c_m

# Strategy: compute c(2) at n=3 first. If feasible time-wise, also c(3) at n=4.
# c(4) at n=5 is heavy.

print(f"=== High-precision c(m) computation (dps={DPS}) ===")
print()

print("--- c(2) at n=3 (length 4913) ---")
t_total = time.time()
c2 = compute_c_m(n=3, m=2, A_MAX=200)
print(f"\nc(2) = {c2}")
print(f"FFT-double  = 0.153247920779")
print(f"Total time: {time.time()-t_total:.1f}s")
print()

# Δ_2 = c(2) - 19/127
delta2 = c2 - mpf(19)/mpf(127)
print(f"Δ_2 = c(2) - 19/127 = {delta2}")
print(f"  approx = {float(delta2):.15e}")
