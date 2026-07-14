"""
Follow-up Model 3: explicit asymptote.

eps_k = C + A * rho^k * cos(k*theta + phi)

Tests whether eps_k -> nonzero constant as k -> inf (Cotaescu persistent residual)
or whether C compatible with 0 (pure decay).

Also runs constrained fits where rho is FIXED at each per-character |lambda|,
to test whether the per-character spectrum can be coerced to fit at all.
"""
from __future__ import annotations
import sys, json
import numpy as np
from scipy.optimize import differential_evolution
sys.stdout.reconfigure(encoding="utf-8")

eps_k = {
    1: +2.000000000000000e-01,
    2: +9.523809523809525e-03,
    3: -5.091986325893010e-03,
    4: -2.452258248318762e-03,
    5: -1.151746915130986e-03,
    6: -4.979056652203831e-04,
    7: -1.175236830374320e-03,
    8: -7.455463672855323e-04,
}
ks = np.array(sorted(eps_k.keys()), dtype=np.float64)
ev = np.array([eps_k[int(k)] for k in ks], dtype=np.float64)

# Per-character eigenvalue magnitudes (8th roots of unity for chi(2))
char_lam_abs = sorted(set(round(abs(3 / (4*np.exp(2j*np.pi*j/8) - 1)), 6) for j in range(8)))
print(f"Per-character |λ| values (deduped): {char_lam_abs}")
print()

# Model 3: C + A * rho^k * cos(k*theta + phi)  -- explicit asymptote
def model_const_cplx(k, C, A, rho, theta, phi):
    return C + A * rho**k * np.cos(k*theta + phi)

def loss_const_cplx(params, ks, ev):
    return float(np.sum((model_const_cplx(ks, *params) - ev)**2))

print("=" * 72)
print("Model 3: eps_k = C + A * rho^k * cos(k*theta + phi)  [fit k=1..8 full]")
print("=" * 72)
bounds3 = [(-0.01, 0.01), (-3, 3), (0.05, 0.95), (0.0, np.pi), (-np.pi, np.pi)]
best_loss3 = np.inf
best_x3 = None
for seed in range(80):
    res = differential_evolution(loss_const_cplx, bounds3, args=(ks, ev), tol=1e-15, seed=seed, maxiter=3000, polish=True)
    if res.fun < best_loss3:
        best_loss3 = res.fun
        best_x3 = res.x
C3, A3, rho3, theta3, phi3 = best_x3
pred3 = model_const_cplx(ks, *best_x3)
resid3 = ev - pred3
ss_tot = np.sum((ev - ev.mean())**2)
r2_3 = 1 - np.sum(resid3**2) / ss_tot
print(f"  C     = {C3:+.10e}  <-- asymptote epsilon_inf")
print(f"  A     = {A3:+.8f}")
print(f"  rho   = {rho3:.8f}")
print(f"  theta = {theta3:.8f} rad  ({theta3/np.pi*180:+.4f}°)")
print(f"  phi   = {phi3:+.8f} rad")
print(f"  SSE   = {best_loss3:.4e}")
print(f"  R²    = {r2_3:.10f}")
print(f"  residuals: {resid3}")
print()
print(f"  Predicted k=9:  {model_const_cplx(9.0, *best_x3):+.6e}")
print(f"  Predicted k=10: {model_const_cplx(10.0, *best_x3):+.6e}")
print(f"  Predicted k=15: {model_const_cplx(15.0, *best_x3):+.6e}")
print(f"  Predicted k=20: {model_const_cplx(20.0, *best_x3):+.6e}")
print(f"  Predicted k=50: {model_const_cplx(50.0, *best_x3):+.6e}")
print()

# Match rho to per-character |lambda|
print(f"  Per-character |λ| comparison:")
for L in char_lam_abs:
    print(f"    |λ|={L:.6f}  diff={abs(rho3-L):.6f}")
print()

# Constrained fits: fix rho at each per-character |lambda|, refit others
print("=" * 72)
print("Constrained fits: fix rho = per-character |λ|, refit A, theta, phi (k=2..8)")
print("=" * 72)
def model_cplx(k, A, rho, theta, phi):
    return A * rho**k * np.cos(k*theta + phi)
ks_tail = ks[1:]
ev_tail = ev[1:]
for rho_fixed in char_lam_abs:
    def loss_constr(p, ks, ev):
        return float(np.sum((model_cplx(ks, p[0], rho_fixed, p[1], p[2]) - ev)**2))
    best_lc = np.inf
    best_xc = None
    for seed in range(30):
        res = differential_evolution(loss_constr, [(-5, 5), (0.0, np.pi), (-np.pi, np.pi)],
                                     args=(ks_tail, ev_tail), tol=1e-14, seed=seed, maxiter=2000, polish=True)
        if res.fun < best_lc:
            best_lc = res.fun
            best_xc = res.x
    Ac, thetac, phic = best_xc
    pred = model_cplx(ks_tail, Ac, rho_fixed, thetac, phic)
    resid = ev_tail - pred
    r2 = 1 - np.sum(resid**2) / np.sum((ev_tail - ev_tail.mean())**2)
    print(f"  rho={rho_fixed:.6f}: A={Ac:+.4e}, theta={thetac/np.pi*180:+8.2f}°, phi={phic/np.pi*180:+8.2f}°  SSE={best_lc:.3e}  R²={r2:.6f}")
print()

# Also test: model with rho fixed=1 and unconstrained complex (decompose: persistent oscillation)
print("=" * 72)
print("Model 4: rho FIXED at 1.0  (pure non-decaying oscillation)  [fit k=1..8]")
print("=" * 72)
def loss_rho1(p, ks, ev):
    A, theta, phi = p
    return float(np.sum((A * np.cos(ks*theta + phi) - ev)**2))
best_4 = np.inf; bx4 = None
for seed in range(50):
    res = differential_evolution(loss_rho1, [(-1, 1), (0, np.pi), (-np.pi, np.pi)], args=(ks, ev),
                                 tol=1e-14, seed=seed, maxiter=2000, polish=True)
    if res.fun < best_4:
        best_4 = res.fun; bx4 = res.x
print(f"  A={bx4[0]:+.6f}, theta={bx4[1]/np.pi*180:+.4f}°, phi={bx4[2]/np.pi*180:+.4f}°")
print(f"  SSE: {best_4:.4e}, R²: {1 - best_4/ss_tot:.4f}  (expected bad — this is sanity check)")
print()

# Save
output = {
    "eps_values": {str(int(k)): float(v) for k, v in eps_k.items()},
    "model_const_plus_cplx": {
        "C_asymptote": float(C3),
        "A": float(A3), "rho": float(rho3), "theta_rad": float(theta3), "phi_rad": float(phi3),
        "SSE": float(best_loss3), "R2": float(r2_3),
    },
    "predictions_model3": {
        "k=9": float(model_const_cplx(9.0, *best_x3)),
        "k=10": float(model_const_cplx(10.0, *best_x3)),
        "k=15": float(model_const_cplx(15.0, *best_x3)),
        "k=20": float(model_const_cplx(20.0, *best_x3)),
        "k=50": float(model_const_cplx(50.0, *best_x3)),
        "k=inf (asymptote)": float(C3),
    },
}
with open("C:/Collatz/eps_asymptote_fit_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print("Saved: C:/Collatz/eps_asymptote_fit_2026_06_01.json")
