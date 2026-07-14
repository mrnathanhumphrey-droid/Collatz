"""
Cross-prime c_∞(p) for p ≡ 1 mod 4: compute via numerical offset_distribution.

Goal: see if c_∞(p) has a universal formula in p, Legendre data, ord_p(2).

Method: same numerical offset_distribution as probe_p3_numerical, but parameterized
over p. Compute c(0)...c(m_max) for several m, extract c_∞ via damped-osc fit.

For each p ≡ 1 mod 4:
  - p = 5  (smallest, very fast)
  - p = 13 (small)
  - p = 17 (our reference, c_∞ = 0.15298912060588517...)
  - p = 29 (next)
  - p = 37 (skip — close to 41)
  - p = 41 (next, larger N)
"""
from __future__ import annotations
import sys, gc, time, json
from fractions import Fraction
from mpmath import mp, mpf, mpc, sqrt, log, pi, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30

def legendre(x, p):
    x %= p
    return 0 if x == 0 else (1 if pow(x, (p-1)//2, p) == 1 else -1)

def offset_distribution_mp(q, n, A_MAX):
    """X = Z(1+qX') self-similarity, distribution of X mod q^n."""
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

def compute_c(p, m, A_MAX=None):
    """c(m) = E[χ_L(σ_m) | v_p(D)=m] via numerical offset_distribution."""
    if A_MAX is None:
        A_MAX = max(50, int(mp.dps * 3.4))  # ~10·log10·dps
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

def fit_c_inf(c_values, m_indices):
    """Damped oscillation fit: c(m) = c_inf + B*ρ^m*cos(mθ+φ)."""
    import numpy as np
    from scipy.optimize import least_squares
    c_arr = np.array([float(c) for c in c_values])
    m_arr = np.array(m_indices)
    def residual(params):
        c_inf, rho, theta, B_cos, B_sin = params
        return c_arr - c_inf - rho**m_arr * (B_cos*np.cos(m_arr*theta) - B_sin*np.sin(m_arr*theta))
    # Initial guess
    p0 = [c_arr[-1], 0.3, 1.0, 0.0, 0.0]
    result = least_squares(residual, p0, max_nfev=10000, ftol=1e-12, xtol=1e-12)
    return result.x

# Process each prime
results = {}
for p in [5, 13, 17, 29, 41]:
    if p % 4 != 1:
        continue
    if p > 17:
        # Check memory
        for m_max in [3, 4]:
            if p ** (m_max+1) > 2 * 10**8:  # ~10 GB at dps=30
                break
        m_max -= 1  # one less to be safe
    else:
        m_max = 4
    print(f"\n=== p = {p}, ord_p(2) = ", end="")
    ord_p2 = 1
    x = 2
    while x != 1:
        x = (x * 2) % p
        ord_p2 += 1
    print(f"{ord_p2}, computing m=0..{m_max} ===")

    c_vals = {}
    t_p_start = time.time()
    for m in range(m_max + 1):
        n = m + 1
        N = p ** n
        t_m = time.time()
        try:
            c_m = compute_c(p, m)
        except (MemoryError, Exception) as e:
            print(f"  m={m} FAILED ({type(e).__name__}: {str(e)[:60]})")
            break
        c_vals[m] = c_m
        print(f"  c({m}) = {c_m}  (N={N}, took {time.time()-t_m:.1f}s)")

    if len(c_vals) < 3:
        print(f"  Insufficient data for p={p}")
        continue

    # Fit
    m_indices = sorted(c_vals.keys())
    c_values_list = [c_vals[m] for m in m_indices]
    try:
        c_inf_est, rho, theta, B_cos, B_sin = fit_c_inf(c_values_list, m_indices)
        print(f"  Fit: c_∞={c_inf_est:.20f}, ρ={rho:.6f}, θ={theta:.6f}")
        print(f"       B_cos={B_cos:.3e}, B_sin={B_sin:.3e}")
    except Exception as e:
        print(f"  Fit failed: {e}")
        c_inf_est = float(c_vals[max(m_indices)])
        rho, theta, B_cos, B_sin = 0, 0, 0, 0

    results[p] = {
        "p": p,
        "ord_p_2": ord_p2,
        "c_vals": {m: str(c_vals[m]) for m in c_vals},
        "c_inf_est": c_inf_est,
        "rho_fit": rho,
        "theta_fit": theta,
        "B_cos": B_cos,
        "B_sin": B_sin,
        "compute_time_s": time.time() - t_p_start,
    }
    print(f"  Total time for p={p}: {time.time()-t_p_start:.1f}s")

# Save and analyze
with open("C:/Collatz/c_inf_cross_prime_2026_06_01.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved results to C:/Collatz/c_inf_cross_prime_2026_06_01.json")

# === Cross-prime PSLQ analysis ===
print(f"\n=== Cross-prime c_∞(p) summary ===")
print(f"  p  | ord_p(2) | c_∞(p) estimate")
for p in sorted(results.keys()):
    r = results[p]
    print(f"  {p:2d} | {r['ord_p_2']:2d}       | {r['c_inf_est']:.20f}")

# Check for clean relationships
print(f"\n=== Test simple relationships ===")
# c_∞(p) vs 1/p, 1/(p-1), Legendre Gauss sum / p, etc.
for p in sorted(results.keys()):
    c = mpf(results[p]['c_inf_est'])
    print(f"\np = {p}, c_∞ = {float(c):.15f}")
    print(f"  c_∞ * p = {float(c * p):.10f}")
    print(f"  c_∞ * (p-1) = {float(c * (p-1)):.10f}")
    print(f"  c_∞ * (p^2-1) = {float(c * (p**2-1)):.10f}")
    print(f"  c_∞ * sqrt(p) = {float(c * sqrt(mpf(p))):.10f}")
    print(f"  c_∞ * p^(3/2) = {float(c * p * sqrt(mpf(p))):.10f}")
    if p >= 17:
        # Compare to c_∞(17) scaled
        c17 = mpf(results[17]['c_inf_est'])
        print(f"  c_∞(p) / c_∞(17) = {float(c / c17):.10f}")
        print(f"  c_∞(17) / c_∞(p) = {float(c17 / c):.10f}")
