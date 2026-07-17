"""
probe_coset_closedform_2026_05_29.py

Stage A of the coset-mass closed-form hunt. The converged per-coset Plancherel mass
fractions (q=17: ~0.5451/0.4549; q=31: ~0.228461/0.151084/0.120455) are limiting
spectral-energy fractions. Goal: get them to as many honest digits as float FFT allows,
check convergence in n and truncation A_MAX, then continued-fraction / rational tests.
"""
from __future__ import annotations
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

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

def subgroup_pow2(q, n):
    N = q ** n; H = np.zeros(N, dtype=bool); x = 1 % N; seen = set()
    while x not in seen:
        seen.add(x); H[x] = True; x = (x * 2) % N
    return H, len(seen)

def coset_labels(q, n, Helems):
    N = q ** n; lab = np.full(N, -1, dtype=np.int64); Helems = np.asarray(Helems, np.int64); nxt = 0
    for u in range(1, N):
        if u % q == 0 or lab[u] != -1: continue
        lab[(u * Helems) % N] = nxt; nxt += 1
    return lab, nxt

def coset_masses(q, n, A_MAX):
    N = q ** n
    P = offset_distribution(q, n, A_MAX)
    mu = np.fft.fft(P); a2 = np.abs(mu) ** 2
    H, ord2 = subgroup_pow2(q, n); Helems = np.nonzero(H)[0]
    lab, ncos = coset_labels(q, n, Helems)
    units = (np.arange(N) % q != 0); mtot = a2[units].sum()
    masses = np.array([a2[lab == c].sum() / mtot for c in range(ncos)])
    # collapse to distinct levels (within conj pair), sorted desc
    levels = sorted({round(float(m), 10) for m in masses}, reverse=True)
    return levels

def main():
    print("Convergence in n and truncation A_MAX (distinct per-coset mass levels):")
    for q, nmax in [(17, 5), (31, 4)]:
        print(f"\n q={q}")
        for A_MAX in (100, 300):
            print(f"  A_MAX={A_MAX}")
            for n in range(2, nmax + 1):
                lv = coset_masses(q, n, A_MAX)
                print(f"    n={n}: " + "  ".join(f"{x:.10f}" for x in lv))
    # continued-fraction / rational test on best estimates
    try:
        from mpmath import mp, mpf, identify
        mp.dps = 30
        print("\nContinued fractions of best (largest-n) estimates:")
        for q, n in [(17, 5), (31, 4)]:
            lv = coset_masses(q, n, 300)
            for x in lv:
                xm = mpf(repr(x))
                # simple CF
                cf = []
                t = xm
                for _ in range(12):
                    a = int(t); cf.append(a); fr = t - a
                    if fr == 0: break
                    t = 1 / fr
                print(f"  q={q} level {x:.10f}: CF={cf}  identify={identify(xm)}")
    except Exception as e:
        print(f"  [mpmath step skipped: {e}]")

if __name__ == "__main__":
    main()
