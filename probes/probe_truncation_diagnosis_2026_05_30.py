"""
probe_truncation_diagnosis_2026_05_30.py

Confirm that K=2 truncated kernel mixes to Haar (= washes out pi's structure).
Specifically: under Haar marginalization of higher digits, row sums of K depend only on
state mod q. Leading eigenvalue = average row sum = 1/q != true 1/3.

This is the "you can't see pi from finite K under Haar prior" diagnostic.
"""
from __future__ import annotations
import numpy as np, sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17; N2 = q**2; N3 = q**3
A_MAX = 60

def offset_distribution(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N); arange = np.arange(N); P_U = np.zeros(N)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a); p = (p * inv2) % N
    P_U /= P_U.sum()
    u_sup = np.nonzero(P_U)[0]; u_wt = P_U[u_sup]; P_S = P_U.copy()
    for j in range(n - 1, 0, -1):
        v_idx = (1 + q * arange) % N; P_V = np.zeros(N); np.add.at(P_V, v_idx, P_S)
        P_new = np.zeros(N)
        for u, w in zip(u_sup, u_wt): P_new[(u * arange) % N] += w * P_V
        P_S = P_new
    return P_S

P3 = offset_distribution(q, 3, A_MAX)
U = np.zeros(N3)
for a in range(1, A_MAX + 1):
    U += 2.0**(-a) * P3[(pow(2, a, N3) * np.arange(N3)) % N3]
U /= U.sum()
D = np.real(np.fft.ifft(np.fft.fft(U) * np.conj(np.fft.fft(U))))

# Row sums by mod-q residue of state s
rs_by_mod_q = {r: [] for r in range(q)}
for s in range(N2):
    deepen_mask = ((s + np.arange(N3)) % q == 0)
    rs_by_mod_q[s % q].append(float(D[deepen_mask].sum()))

print("Row sums of K=2 kernel, grouped by s mod q:")
print(f"  {'r':>3} {'min':>10} {'max':>10} {'mean':>10} {'std':>10}")
for r in range(q):
    vals = np.array(rs_by_mod_q[r])
    print(f"  {r:>3} {vals.min():>10.6f} {vals.max():>10.6f} {vals.mean():>10.6f} {vals.std():>10.2e}")

print(f"\n  -> if std ~ 0 within each row, then row sum depends only on s mod q (as predicted).")
print(f"  -> mean across all rows = {np.mean([np.mean(v) for v in rs_by_mod_q.values()]):.6f} = 1/q = {1/q:.6f}")
print(f"  -> the K=2 Haar-marginalized chain's leading eigenvalue = avg row sum = 1/q, NOT 1/3.")
print(f"\n  The true pi-weighted continuation rate of 1/3 requires the CORRECT (non-Haar) higher-digit")
print(f"  marginal — which is itself pi's own higher-digit distribution. Self-referential.")
print(f"  Resolution: lazy q-adic simulation OR self-consistent power iteration with on-demand digits.")
