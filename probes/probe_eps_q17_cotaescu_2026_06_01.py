"""
Cotaescu complex-eigenvalue analysis for c_inf(17) — the q=17 substrate analog
of the c=7/45 ε_k fit.

Goal:
  1. Compute c(m) for m=0..3 at q=17 at high precision (dps=80)
  2. Compute eps(m) = c(m) - c_inf(17) using the 50-digit reference value
  3. Fit Cotaescu Models 1, 3 (single complex; constant + complex)
  4. Compare fitted (rho, theta) to the 8 closed-form per-character eigenvalues
       lambda_dom(chi) = 3/(4 chi(2) - 1)   for chi(2) ∈ 8th roots of unity

Key advantage over q=3: the per-character formula gives 8 EXPLICIT candidate
eigenvalues in closed form. A direct match would close the spectrum question
algebraically, not just empirically.
"""
from __future__ import annotations
import sys, gc, time, json
from mpmath import mp, mpf, mpc, sqrt, log, pi, pslq
from scipy.optimize import differential_evolution
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 80
P = 17

# Reference c_inf(17) — 50 digits
C_INF_17 = mpf("0.15298912060588517527891674877413229926086222622334")

def legendre(x, p):
    x %= p
    return 0 if x == 0 else (1 if pow(x, (p-1)//2, p) == 1 else -1)

def offset_distribution_mp(q, n, A_MAX=120):
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

def compute_c(p, m, A_MAX=120):
    n = m + 1
    N = p ** n
    chi_L = [legendre(j, p) for j in range(p)]
    P_X = offset_distribution_mp(p, n, A_MAX)
    qm = p ** m
    num = mpf(0); den = mpf(0)
    for j in range(1, p):
        d = j * qm
        Pd = mpf(0)
        for y in range(N):
            Pd += P_X[y] * P_X[(y - d) % N]
        num += Pd * mpf(chi_L[j])
        den += Pd
    return num / den

# ============================================
# 1. Compute c(m) for m=0..M_MAX
# ============================================
M_MAX = 3   # m=3 -> N = 17^4 = 83521; m=4 would be 17^5 = 1.4M (much slower)

print(f"q={P}, mp.dps={mp.dps}, M_MAX={M_MAX}")
print(f"c_inf({P}) reference: {C_INF_17}")
print()

c_vals = {}
t0 = time.time()
for m in range(M_MAX + 1):
    N = P ** (m + 1)
    print(f"[c({m})]  N={N}  computing...", flush=True)
    tm = time.time()
    c_m = compute_c(P, m, A_MAX=120)
    c_vals[m] = c_m
    eps_m = c_m - C_INF_17
    print(f"  c({m}) = {c_m}")
    print(f"  eps({m}) = c({m}) - c_inf = {eps_m}")
    print(f"  |eps({m})| = {abs(eps_m)}  log10 = {mp.log10(abs(eps_m)) if abs(eps_m) > 0 else 'inf'}")
    print(f"  time = {time.time()-tm:.1f}s")
    print()
print(f"Total c(m) compute: {time.time()-t0:.1f}s")
print()

# ============================================
# 2. Build eps array (as Python floats for fitting)
# ============================================
ms = sorted(c_vals.keys())
eps_arr_mp = [c_vals[m] - C_INF_17 for m in ms]
print(f"eps(m) values:")
for m, e in zip(ms, eps_arr_mp):
    print(f"  m={m}:  eps = {float(e):+.6e}  |eps|·p^m = {float(abs(e) * P**m):+.6e}")
print()

# ============================================
# 3. Per-character eigenvalues  lambda_dom(chi) = 3 / (4 chi(2) - 1)
# ============================================
print(f"Per-character eigenvalues  lambda(chi) = 3/(4 chi(2) - 1):")
print(f"  chi(2) ∈ 8th roots of unity (since ord_{P}(2) = 8)")
print(f"  {'j':>3} {'chi(2)':>20} {'lambda':>26} {'|lambda|':>12} {'arg deg':>10}")
char_lams = []
for j in range(8):
    omega = mpc(0,0)
    omega = mpc(mp.cos(2*mp.pi*j/8), mp.sin(2*mp.pi*j/8))
    lam = mpc(3) / (mpc(4)*omega - mpc(1))
    char_lams.append(complex(float(lam.real), float(lam.imag)))
    abs_lam = abs(complex(float(lam.real), float(lam.imag)))
    arg_lam = np.angle(complex(float(lam.real), float(lam.imag))) * 180 / np.pi
    print(f"  {j:>3}  {complex(float(omega.real), float(omega.imag)):>20.4f}  {complex(float(lam.real), float(lam.imag)):>26.6f}  {abs_lam:>12.6f}  {arg_lam:>+10.3f}")
print()

# ============================================
# 4. Fit Cotaescu models if we have enough data
# ============================================
ms_arr = np.array(ms, dtype=np.float64)
eps_arr = np.array([float(e) for e in eps_arr_mp], dtype=np.float64)

print(f"Fit data: {len(ms_arr)} points")
if len(ms_arr) >= 3:
    # Model 1: eps(m) = A * rho^m * cos(m*theta + phi)
    def model1(m, A, rho, theta, phi):
        return A * rho**m * np.cos(m*theta + phi)
    def loss1(p, ms, eps):
        return float(np.sum((model1(ms, *p) - eps)**2))

    print("Model 1: eps(m) = A * rho^m * cos(m*theta + phi)")
    bounds1 = [(-2, 2), (0.05, 1.5), (0.0, np.pi), (-np.pi, np.pi)]
    best_loss1 = np.inf; best_x1 = None
    for seed in range(60):
        res = differential_evolution(loss1, bounds1, args=(ms_arr, eps_arr), tol=1e-15, seed=seed, maxiter=3000, polish=True)
        if res.fun < best_loss1:
            best_loss1 = res.fun; best_x1 = res.x
    A1, rho1, theta1, phi1 = best_x1
    pred1 = model1(ms_arr, *best_x1)
    resid1 = eps_arr - pred1
    if np.sum((eps_arr - eps_arr.mean())**2) > 0:
        r2_1 = 1 - np.sum(resid1**2) / np.sum((eps_arr - eps_arr.mean())**2)
    else:
        r2_1 = float('nan')
    print(f"  A={A1:+.6f}, rho={rho1:.6f}, theta={theta1*180/np.pi:.3f}°, phi={phi1*180/np.pi:.3f}°")
    print(f"  SSE={best_loss1:.4e}, R²={r2_1:.6f}")
    print(f"  Residuals: {resid1}")
    print()

    # Compare to per-character eigenvalues
    print(f"  Compare fitted rho={rho1:.4f}, |theta|={abs(theta1)*180/np.pi:.2f}° to per-character:")
    for j, lam in enumerate(char_lams):
        d_rho = abs(rho1 - abs(lam))
        d_theta_deg = abs(abs(theta1*180/np.pi) - abs(np.angle(lam)*180/np.pi))
        marker = " <-- MATCH" if (d_rho < 0.05 and d_theta_deg < 5) else ""
        print(f"    j={j}:  |λ|={abs(lam):.6f} (d_rho={d_rho:.4f}),  arg={np.angle(lam)*180/np.pi:+8.3f}° (d_theta={d_theta_deg:.2f}°){marker}")
    print()

# Model 3: eps(m) = C + A * rho^m * cos(m*theta + phi)
if len(ms_arr) >= 4:
    def model3(m, C, A, rho, theta, phi):
        return C + A * rho**m * np.cos(m*theta + phi)
    def loss3(p, ms, eps):
        return float(np.sum((model3(ms, *p) - eps)**2))

    print("Model 3: eps(m) = C + A * rho^m * cos(m*theta + phi)")
    eps_scale = max(abs(eps_arr))
    bounds3 = [(-eps_scale, eps_scale), (-2, 2), (0.05, 1.5), (0.0, np.pi), (-np.pi, np.pi)]
    best_loss3 = np.inf; best_x3 = None
    for seed in range(80):
        res = differential_evolution(loss3, bounds3, args=(ms_arr, eps_arr), tol=1e-15, seed=seed, maxiter=3000, polish=True)
        if res.fun < best_loss3:
            best_loss3 = res.fun; best_x3 = res.x
    C3, A3, rho3, theta3, phi3 = best_x3
    pred3 = model3(ms_arr, *best_x3)
    resid3 = eps_arr - pred3
    if np.sum((eps_arr - eps_arr.mean())**2) > 0:
        r2_3 = 1 - np.sum(resid3**2) / np.sum((eps_arr - eps_arr.mean())**2)
    else:
        r2_3 = float('nan')
    print(f"  C={C3:+.6e}, A={A3:+.6f}, rho={rho3:.6f}, theta={theta3*180/np.pi:.3f}°, phi={phi3*180/np.pi:.3f}°")
    print(f"  SSE={best_loss3:.4e}, R²={r2_3:.6f}")
    print(f"  Residuals: {resid3}")
    print()

    # Compare to per-character
    print(f"  Compare Model3 rho={rho3:.4f}, |theta|={abs(theta3)*180/np.pi:.2f}° to per-character:")
    for j, lam in enumerate(char_lams):
        d_rho = abs(rho3 - abs(lam))
        d_theta_deg = abs(abs(theta3*180/np.pi) - abs(np.angle(lam)*180/np.pi))
        marker = " <-- MATCH" if (d_rho < 0.05 and d_theta_deg < 5) else ""
        print(f"    j={j}:  |λ|={abs(lam):.6f} (d_rho={d_rho:.4f}),  arg={np.angle(lam)*180/np.pi:+8.3f}° (d_theta={d_theta_deg:.2f}°){marker}")

# Save
output = {
    "P": P,
    "mp_dps": mp.dps,
    "M_MAX": M_MAX,
    "c_inf_reference": str(C_INF_17),
    "c_vals_str": {str(m): str(c) for m, c in c_vals.items()},
    "eps_vals_str": {str(m): str(c - C_INF_17) for m, c in c_vals.items()},
    "char_eigenvalues": [
        {"j": j, "lambda_re": float(l.real), "lambda_im": float(l.imag),
         "abs": float(abs(l)), "arg_deg": float(np.angle(l)*180/np.pi)}
        for j, l in enumerate(char_lams)
    ],
}
if len(ms_arr) >= 3:
    output["model_1"] = {
        "A": float(A1), "rho": float(rho1), "theta_rad": float(theta1), "phi_rad": float(phi1),
        "SSE": float(best_loss1), "R2": float(r2_1) if not np.isnan(r2_1) else None,
    }
if len(ms_arr) >= 4:
    output["model_3"] = {
        "C": float(C3), "A": float(A3), "rho": float(rho3), "theta_rad": float(theta3), "phi_rad": float(phi3),
        "SSE": float(best_loss3), "R2": float(r2_3) if not np.isnan(r2_3) else None,
    }
with open("C:/Collatz/eps_q17_cotaescu_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nTotal wall time: {time.time()-t0:.1f}s")
print(f"Saved: C:/Collatz/eps_q17_cotaescu_2026_06_01.json")
