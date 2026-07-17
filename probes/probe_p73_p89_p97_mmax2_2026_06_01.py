"""
Cross-prime c_∞(p) for p ∈ {73, 89, 97} at m_max=2.

Using c(2) as c_∞ approximation. With fast convergence observed at p=41,
c(2) is good to ~3-5 digits.

Memory: N = p^3 max.
  p=73:  389k entries → ~30 MB
  p=89:  705k entries → ~50 MB
  p=97:  913k entries → ~70 MB
Total: ~5 min compute.
"""
from __future__ import annotations
import sys, gc, time, json
from mpmath import mp, mpf, sqrt, log, pi, pslq
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

def compute_c(p, m, A_MAX=80):
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

def fit_c_inf(c_floats, m_indices):
    import numpy as np
    from scipy.optimize import least_squares
    c_arr = np.array(c_floats)
    m_arr = np.array(m_indices)
    def residual(params):
        c_inf, rho, theta, B = params
        return c_arr - c_inf - B * rho**m_arr * np.cos(m_arr*theta)
    p0 = [c_arr[-1], 0.3, 0.5, 0.0]
    try:
        result = least_squares(residual, p0, max_nfev=10000, ftol=1e-12, xtol=1e-12)
        return result.x
    except:
        return None

results = {}
for p in [73, 89, 97]:
    p_mod8 = p % 8
    ord_p2 = 1
    x = 2
    while x != 1:
        x = (x * 2) % p
        ord_p2 += 1
    print(f"\n=== p = {p} (mod 8 = {p_mod8}, ord_p(2) = {ord_p2}) ===", flush=True)
    m_max = 2
    c_vals = {}
    t_p_start = time.time()
    for m in range(m_max + 1):
        N = p ** (m + 1)
        t_m = time.time()
        print(f"  Computing c({m}) at N={N}...", flush=True)
        c_m = compute_c(p, m, A_MAX=80)
        c_vals[m] = c_m
        print(f"    c({m}) = {c_m}  (took {time.time()-t_m:.1f}s)", flush=True)
    # Fit (only 3 points → estimate via simple extrapolation, treat c(2) as c_∞)
    c_inf_est = float(c_vals[2])
    results[p] = {
        "p": p,
        "p_mod_8": p_mod8,
        "ord_p_2": ord_p2,
        "c_vals": {m: str(c_vals[m]) for m in c_vals},
        "c_inf_est_via_c2": c_inf_est,
        "c1_minus_c2": float(c_vals[1] - c_vals[2]),
        "compute_time_s": time.time() - t_p_start,
    }
    print(f"  Total: {time.time()-t_p_start:.1f}s, c_∞ ≈ c(2) = {c_inf_est:.10f}", flush=True)
    print(f"  c(1)-c(2) = {results[p]['c1_minus_c2']:+.3e} (gauge of precision)", flush=True)

with open("C:/Collatz/c_inf_p73_p89_p97_mmax2_2026_06_01.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved.", flush=True)

# Full family
print(f"\n=== Full family p ≡ 1 mod 8 with nontrivial c_∞ ===", flush=True)
c_inf_17 = mpf("0.15298912060588517527891674877413229926086222622334")
c_inf_41 = mpf("0.08869012630946791764")
family = [(17, c_inf_17), (41, c_inf_41)]
for p in sorted(results.keys()):
    family.append((p, mpf(results[p]['c_inf_est_via_c2'])))
print(f"  p   | c_∞(p)                | ord_p(2) | (p-1)/ord_p(2)")
for p, c_p in family:
    if p == 17:
        ord_p2 = 8
    elif p == 41:
        ord_p2 = 20
    else:
        ord_p2 = results[p]['ord_p_2']
    print(f"  {p:3d} | {float(c_p):.18f} | {ord_p2:2d}    | {(p-1)//ord_p2}", flush=True)

# Universal-formula PSLQ probes
print(f"\n=== Test simple universal formulas ===", flush=True)

# Probe 1: c_∞(p) = some function of p only
print(f"\n  c_∞(p) raw values:", flush=True)
for p, c_p in family:
    print(f"    p={p:3d}: {float(c_p):.10f}", flush=True)

print(f"\n  c_∞(p) · p:", flush=True)
for p, c_p in family:
    print(f"    p={p:3d}: {float(c_p * p):.10f}", flush=True)

print(f"\n  c_∞(p) · (p-1):", flush=True)
for p, c_p in family:
    print(f"    p={p:3d}: {float(c_p * (p-1)):.10f}", flush=True)

print(f"\n  c_∞(p) · sqrt(p):", flush=True)
for p, c_p in family:
    print(f"    p={p:3d}: {float(c_p * sqrt(mpf(p))):.10f}", flush=True)

print(f"\n  c_∞(p) · ord_p(2):", flush=True)
for p, c_p in family:
    if p == 17: o = 8
    elif p == 41: o = 20
    else: o = results[p]['ord_p_2']
    print(f"    p={p:3d}: c_∞ · ord = {float(c_p * o):.10f}", flush=True)

# Probe 2: PSLQ each c_∞(p) against {1, 1/p, 1/(p-1), ord_p(2)/p^2, ...}
print(f"\n=== PSLQ c_∞(p) against p-specific basis ===", flush=True)
for p, c_p in family:
    if p == 17: ord_p2 = 8
    elif p == 41: ord_p2 = 20
    else: ord_p2 = results[p]['ord_p_2']
    basis = [
        mpf(1),
        mpf(1) / p,
        mpf(1) / (p - 1),
        mpf(1) / p**2,
        mpf(ord_p2) / p,
        mpf(ord_p2) / (p**2),
        mpf(ord_p2) / (p * (p-1)),
        log(mpf(p)) / p,
        1 / sqrt(mpf(p)),
        pi / (p * sqrt(mpf(p))),
        mpf(p-1) / (p**2),
        mpf(ord_p2) / (p * ord_p2 - 1),  # related to Lambert sum
    ]
    bnames = ["1", "1/p", "1/(p-1)", "1/p^2", "ord/p", "ord/p^2", "ord/(p(p-1))",
              "log(p)/p", "1/sqrt(p)", "pi/(p sqrt(p))", "(p-1)/p^2", "ord/(p·ord-1)"]
    for tol_exp in [3, 5, 7]:
        tol = mpf(10) ** (-tol_exp)
        rel = pslq([c_p] + basis, tol=tol, maxcoeff=200)
        if rel is not None and rel[0] != 0:
            terms = [f"({rel[0]:+d})*c_∞({p})"]
            for c, n in zip(rel[1:], bnames):
                if c != 0:
                    terms.append(f"({c:+d})*{n}")
            if len(terms) <= 6:
                print(f"  p={p}, tol=10^-{tol_exp}: {' '.join(terms)}", flush=True)
                break

# Probe 3: Cross-prime LINEAR fit
# Hypothesize c_∞(p) = a/p + b/(p-1) + c·ord/p + d for small rational a,b,c,d
print(f"\n=== Cross-prime LINEAR fit: c_∞(p) = a/p + b/(p-1) + c·ord_p(2)/p ===", flush=True)
import numpy as np
M = []
y = []
for p, c_p in family:
    if p == 17: ord_p2 = 8
    elif p == 41: ord_p2 = 20
    else: ord_p2 = results[p]['ord_p_2']
    M.append([1/p, 1/(p-1), ord_p2/p, 1.0])
    y.append(float(c_p))
M = np.array(M); y = np.array(y)
coefs, res, rank, _ = np.linalg.lstsq(M, y, rcond=None)
print(f"  Fit: c_∞(p) ≈ {coefs[0]:.6f}/p + {coefs[1]:.6f}/(p-1) + {coefs[2]:.6f}·ord_p(2)/p + {coefs[3]:.6f}")
print(f"  Residual: {res if res.size > 0 else 'exact-fit'}")
pred = M @ coefs
print(f"  Predictions vs actual:")
for (p, c_p), pr in zip(family, pred):
    print(f"    p={p}: actual={float(c_p):.8f}, pred={pr:.8f}, diff={float(c_p)-pr:+.2e}")

print(f"\n=== Done ===", flush=True)
