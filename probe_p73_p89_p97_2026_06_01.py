"""
Extend cross-prime c_∞(p) to p ∈ {73, 89, 97} — all p ≡ 1 mod 8.

Goal: 3 more nontrivial c_∞(p) data points → 5 total (17, 41, 73, 89, 97)
→ enough for universal-in-p PSLQ fit at low precision.

Memory budget: each p with m_max=3 needs N = p^4 ≈ p^4 entries.
  p=73, p^4 = 28.4M → ~1.4 GB
  p=89, p^4 = 62.7M → ~3 GB
  p=97, p^4 = 88.5M → ~4 GB
All feasible at 42 GB free.
"""
from __future__ import annotations
import sys, gc, time, json
from mpmath import mp, mpf
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30

def legendre(x, p):
    x %= p
    return 0 if x == 0 else (1 if pow(x, (p-1)//2, p) == 1 else -1)

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
    P_S = P_U
    for j in range(n - 1, 0, -1):
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
    return P_S

def compute_c(p, m, A_MAX=100):
    n = m + 1
    N = p ** n
    chi_L_table = [legendre(j, p) for j in range(p)]
    P_X = offset_distribution_mp(p, n, A_MAX)
    qm = p ** m
    num = mpf(0); den = mpf(0)
    for j in range(1, p):
        d = j * qm
        Pd = mpf(0)
        for y in range(N):
            Pd += P_X[y] * P_X[(y - d) % N]
        num += Pd * mpf(chi_L_table[j])
        den += Pd
    return num / den

def fit_c_inf(c_values_floats, m_indices):
    import numpy as np
    from scipy.optimize import least_squares
    c_arr = np.array(c_values_floats)
    m_arr = np.array(m_indices)
    def residual(params):
        c_inf, rho, theta, B_cos, B_sin = params
        return c_arr - c_inf - rho**m_arr * (B_cos*np.cos(m_arr*theta) - B_sin*np.sin(m_arr*theta))
    p0 = [c_arr[-1], 0.3, 0.5, 0.0, 0.0]
    result = least_squares(residual, p0, max_nfev=10000, ftol=1e-12, xtol=1e-12)
    return result.x

results = {}
for p in [73, 89, 97]:
    if p % 4 != 1:
        continue
    # Verify p mod 8
    p_mod8 = p % 8
    # ord_p(2)
    ord_p2 = 1
    x = 2
    while x != 1:
        x = (x * 2) % p
        ord_p2 += 1
    print(f"\n=== p = {p} (mod 8 = {p_mod8}, ord_p(2) = {ord_p2}) ===", flush=True)

    m_max = 3
    c_vals = {}
    t_p_start = time.time()
    for m in range(m_max + 1):
        n = m + 1
        N = p ** n
        t_m = time.time()
        print(f"  Computing c({m}) at N={N}...", flush=True)
        try:
            c_m = compute_c(p, m)
        except (MemoryError, Exception) as e:
            print(f"    FAILED: {type(e).__name__}: {str(e)[:80]}", flush=True)
            break
        c_vals[m] = c_m
        print(f"    c({m}) = {c_m}  (took {time.time()-t_m:.1f}s)", flush=True)

    if len(c_vals) < 3:
        print(f"  Insufficient data", flush=True)
        continue

    m_indices = sorted(c_vals.keys())
    c_floats = [float(c_vals[m]) for m in m_indices]
    try:
        c_inf_est, rho, theta, B_cos, B_sin = fit_c_inf(c_floats, m_indices)
        print(f"  Fit: c_∞={c_inf_est:.20f}, ρ={rho:.6f}, θ={theta:.6f}", flush=True)
    except Exception as e:
        c_inf_est = c_floats[-1]
        rho = theta = B_cos = B_sin = 0
        print(f"  Fit failed, using last c value: {c_inf_est}", flush=True)

    results[p] = {
        "p": p,
        "p_mod_8": p_mod8,
        "ord_p_2": ord_p2,
        "c_vals": {m: str(c_vals[m]) for m in c_vals},
        "c_inf_est": c_inf_est,
        "rho_fit": rho,
        "theta_fit": theta,
        "compute_time_s": time.time() - t_p_start,
    }
    print(f"  Total: {time.time()-t_p_start:.1f}s", flush=True)

with open("C:/Collatz/c_inf_p73_p89_p97_2026_06_01.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to C:/Collatz/c_inf_p73_p89_p97_2026_06_01.json", flush=True)

# Cross-prime PSLQ
print(f"\n=== Universal-formula PSLQ ===", flush=True)
print(f"  All p ≡ 1 mod 8 with nontrivial c_∞(p):", flush=True)
c_inf_17 = mpf("0.15298912060588517527891674877413229926086222622334")
c_inf_41 = mpf("0.08869012630946791764")
data_pts = [(17, c_inf_17), (41, c_inf_41)]
for p, r in sorted(results.items()):
    c_p = mpf(r['c_inf_est'])
    data_pts.append((p, c_p))
    print(f"    p={p:3d}, c_∞ = {float(c_p):.20f}", flush=True)

# Test if c_∞(p) follows: c_∞(p) = A/p + B/p^2 + ... or similar
print(f"\n  c_∞(p) * p:", flush=True)
for p, c_p in data_pts:
    print(f"    p={p:3d}, c_∞·p = {float(c_p * p):.10f}", flush=True)

print(f"\n  c_∞(p) * sqrt(p):", flush=True)
for p, c_p in data_pts:
    from mpmath import sqrt
    print(f"    p={p:3d}, c_∞·√p = {float(c_p * sqrt(mpf(p))):.10f}", flush=True)

print(f"\n  c_∞(p) * (p-1) / 16:", flush=True)
for p, c_p in data_pts:
    print(f"    p={p:3d}, c_∞·(p-1)/16 = {float(c_p * (p-1) / 16):.10f}", flush=True)

# Try PSLQ for 4+ data points
if len(data_pts) >= 4:
    from mpmath import pslq, log, sqrt, pi
    print(f"\n  PSLQ each c_∞(p) vs simple p-dependent constants:", flush=True)
    for p, c_p in data_pts:
        basis = [
            mpf(1),
            mpf(1) / p,
            mpf(1) / (p - 1),
            mpf(1) / (p * (p - 1)),
            log(mpf(p)) / p,
            1 / sqrt(mpf(p)),
            log(mpf(p)),
            pi / p,
            log(mpf(2)) / p,
        ]
        bnames = ["1", "1/p", "1/(p-1)", "1/(p(p-1))", "log(p)/p", "1/sqrt(p)", "log(p)", "pi/p", "log(2)/p"]
        for tol_exp in [4, 6, 8]:
            tol = mpf(10) ** (-tol_exp)
            rel = pslq([c_p] + basis, tol=tol, maxcoeff=200)
            if rel is not None and rel[0] != 0:
                terms = [f"({rel[0]:+d})*c_∞({p})"] + [f"({c:+d})*{n}" for c, n in zip(rel[1:], bnames) if c != 0]
                if len(terms) <= 6:
                    print(f"    p={p}, tol=10^-{tol_exp}: {' '.join(terms)}", flush=True)
                    break

print("\n=== Done ===", flush=True)
