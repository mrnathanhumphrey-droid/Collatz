"""
R77 subdominant-rate reframe: fit eps_k as a real-projection of a COMPLEX
eigenvalue, not a real-exponential decay.

Tests:
  (1) eps_k = A * rho^k * cos(k*theta + phi)  -- single complex pair
  (2) Compare fitted (rho, theta) to per-character eigenvalues lambda(chi) = 3/(4*chi(2) - 1)
      where chi(2) runs over 8th roots of unity (ord_17(2) = 8).
  (3) eps_k = B*r^k + A*rho^k*cos(k*theta + phi)  -- one real + one complex pair
  (4) Predict eps_9, eps_10 from each model.

The conjecture being tested: R77.2 ("subdominant rate (1/2)^k") was wrong
about MODEL FORM, not just magnitude. Real-exponential decay is the wrong
shape; the empirical (1/2)^k fit for k=2..6 is the early-k tangent of a
slowly-oscillating cosine envelope whose zero falls near k=6, causing the
k=7 envelope break.
"""
from __future__ import annotations
import sys, json
import numpy as np
from scipy.optimize import differential_evolution
sys.stdout.reconfigure(encoding="utf-8")

# Exact eps values from C:/Collatz/experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json
eps_k = {
    1: +2.000000000000000e-01,  # = 1/5
    2: +9.523809523809525e-03,  # = 1/105
    3: -5.091986325893010e-03,
    4: -2.452258248318762e-03,
    5: -1.151746915130986e-03,
    6: -4.979056652203831e-04,
    7: -1.175236830374320e-03,
    8: -7.455463672855323e-04,
}
ks = np.array(sorted(eps_k.keys()), dtype=np.float64)
ev = np.array([eps_k[int(k)] for k in ks], dtype=np.float64)

print("Exact eps_k values:")
print(f"  {'k':>3} {'eps_k':>16} {'|eps|·2^k':>12} {'sign':>6}")
for k, v in zip(ks, ev):
    print(f"  {int(k):>3} {v:>+16.8e} {abs(v) * 2**int(k):>12.6f} {'+' if v >= 0 else '-':>6}")
print()

# Per-character eigenvalues for q=17
print("Per-character eigenvalues λ_dom(χ) = 3/(4χ(2)−1)  [χ(2) ∈ 8th roots of unity]:")
print(f"  {'idx':>3} {'χ(2)':>22} {'λ':>30} {'|λ|':>12} {'arg(λ)':>14}")
char_lams = []
for j in range(8):
    omega = np.exp(2j * np.pi * j / 8)
    lam = 3 / (4 * omega - 1)
    char_lams.append(lam)
    print(f"  {j:>3} {str(complex(np.round(omega, 4))):>22}  {str(complex(np.round(lam, 5))):>30} {abs(lam):>12.6f} {np.angle(lam):>+14.6f} ({np.angle(lam)/np.pi*180:+7.2f}°)")
print()

# ----- Model 1: single complex eigenvalue -----
def model_cplx(k, A, rho, theta, phi):
    return A * rho**k * np.cos(k*theta + phi)

def loss_cplx(params, ks, ev):
    return float(np.sum((model_cplx(ks, *params) - ev)**2))

# Fit on k=2..8 (exclude k=1 which is anomalously large -- likely initialization term)
ks_tail = ks[1:]
ev_tail = ev[1:]

print("=" * 72)
print("Model 1: eps_k = A * rho^k * cos(k*theta + phi)  [fit k=2..8]")
print("=" * 72)
bounds = [(-2, 2), (0.05, 0.99), (0.0, np.pi), (-np.pi, np.pi)]
best_loss = np.inf
best_x = None
for seed in range(50):
    res = differential_evolution(loss_cplx, bounds, args=(ks_tail, ev_tail), tol=1e-14, seed=seed, maxiter=2000, polish=True)
    if res.fun < best_loss:
        best_loss = res.fun
        best_x = res.x
A1, rho1, theta1, phi1 = best_x
pred_tail = model_cplx(ks_tail, A1, rho1, theta1, phi1)
resid_tail = ev_tail - pred_tail
ss_tot = np.sum((ev_tail - ev_tail.mean())**2)
r2 = 1 - np.sum(resid_tail**2) / ss_tot
print(f"  A     = {A1:+.8f}")
print(f"  rho   = {rho1:.8f}")
print(f"  theta = {theta1:.8f} rad  ({theta1/np.pi*180:+.4f}°)")
print(f"  phi   = {phi1:+.8f} rad")
print(f"  SSE   = {np.sum(resid_tail**2):.4e}")
print(f"  R²    = {r2:.10f}")
print(f"  residuals: {resid_tail}")
print()
print(f"  Predicted k=1: {model_cplx(1.0, A1, rho1, theta1, phi1):+.6e}  (actual {ev[0]:+.6e}  -- expect mismatch)")
print(f"  Predicted k=9: {model_cplx(9.0, A1, rho1, theta1, phi1):+.6e}")
print(f"  Predicted k=10: {model_cplx(10.0, A1, rho1, theta1, phi1):+.6e}")
print()

# Match fitted (rho, theta) against per-character eigenvalues
print("Match fitted (rho, theta) against per-character eigenvalues:")
print(f"  Fitted: rho={rho1:.6f}, |theta|={abs(theta1):.6f} ({abs(theta1)/np.pi*180:.2f}°)")
print(f"  {'idx':>3} {'|λ|':>12} {'|arg(λ)|':>14} {'|rho - |λ||':>14} {'|theta - |arg||':>16}")
for j, lam in enumerate(char_lams):
    abs_lam = abs(lam)
    abs_arg = abs(np.angle(lam))
    d_rho = abs(rho1 - abs_lam)
    d_theta = abs(abs(theta1) - abs_arg)
    print(f"  {j:>3} {abs_lam:>12.6f} {abs_arg:>14.6f} {d_rho:>14.6f} {d_theta:>16.6f}")
print()

# ----- Model 2: real exponential + complex pair (fit k=1..8, full) -----
def model_real_plus_cplx(k, B, r, A, rho, theta, phi):
    return B * r**k + A * rho**k * np.cos(k*theta + phi)

def loss_rpc(params, ks, ev):
    return float(np.sum((model_real_plus_cplx(ks, *params) - ev)**2))

print("=" * 72)
print("Model 2: eps_k = B*r^k + A*rho^k*cos(k*theta + phi)  [fit k=1..8 full]")
print("=" * 72)
bounds2 = [(-2, 2), (0.05, 0.99), (-2, 2), (0.05, 0.99), (0.0, np.pi), (-np.pi, np.pi)]
best_loss2 = np.inf
best_x2 = None
for seed in range(80):
    res = differential_evolution(loss_rpc, bounds2, args=(ks, ev), tol=1e-14, seed=seed, maxiter=3000, polish=True)
    if res.fun < best_loss2:
        best_loss2 = res.fun
        best_x2 = res.x
B2, r2_real, A2, rho2, theta2, phi2 = best_x2
pred2 = model_real_plus_cplx(ks, B2, r2_real, A2, rho2, theta2, phi2)
resid2 = ev - pred2
ss_tot2 = np.sum((ev - ev.mean())**2)
r2_score2 = 1 - np.sum(resid2**2) / ss_tot2
print(f"  Real:    B={B2:+.8f}, r={r2_real:.8f}")
print(f"  Complex: A={A2:+.8f}, rho={rho2:.8f}, theta={theta2:.6f} rad ({theta2/np.pi*180:+.4f}°), phi={phi2:+.6f}")
print(f"  SSE: {np.sum(resid2**2):.4e}")
print(f"  R²:  {r2_score2:.10f}")
print(f"  residuals: {resid2}")
print()
print(f"  Predicted k=9: {model_real_plus_cplx(9.0, B2, r2_real, A2, rho2, theta2, phi2):+.6e}")
print(f"  Predicted k=10: {model_real_plus_cplx(10.0, B2, r2_real, A2, rho2, theta2, phi2):+.6e}")
print()

# Match
print("Match Model-2 complex (rho, theta) against per-character eigenvalues:")
for j, lam in enumerate(char_lams):
    abs_lam = abs(lam)
    abs_arg = abs(np.angle(lam))
    d_rho = abs(rho2 - abs_lam)
    d_theta = abs(abs(theta2) - abs_arg)
    print(f"  {j:>3} {abs_lam:>12.6f} {abs_arg:>14.6f} {d_rho:>14.6f} {d_theta:>16.6f}")
print()
print("Match Model-2 real r against per-character eigenvalues:")
for j, lam in enumerate(char_lams):
    abs_lam = abs(lam)
    if abs(np.angle(lam)) < 1e-6 or abs(abs(np.angle(lam)) - np.pi) < 1e-6:
        print(f"  {j:>3} {abs_lam:>12.6f}  real  |r - |λ|| = {abs(r2_real - abs_lam):.6f}")
print()

# ----- Pure real-exponential model for comparison (the OLD R77 model) -----
print("=" * 72)
print("Model 0 (R77 baseline): eps_k = A * r^k  [pure real exponential]")
print("=" * 72)
def model_real(k, A, r):
    return A * r**k
def loss_real(params, ks, ev):
    return float(np.sum((model_real(ks, *params) - ev)**2))
best_loss0 = np.inf
best_x0 = None
for seed in range(30):
    res = differential_evolution(loss_real, [(-2, 2), (0.05, 0.99)], args=(ks_tail, ev_tail), tol=1e-14, seed=seed, maxiter=1000, polish=True)
    if res.fun < best_loss0:
        best_loss0 = res.fun
        best_x0 = res.x
A0, r0 = best_x0
pred0 = model_real(ks_tail, A0, r0)
resid0 = ev_tail - pred0
r2_0 = 1 - np.sum(resid0**2) / ss_tot
print(f"  A={A0:+.6f}, r={r0:.6f}")
print(f"  SSE: {np.sum(resid0**2):.4e}")
print(f"  R²:  {r2_0:.6f}")
print()

# Summary table
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"  Model           SSE             R²            Params")
print(f"  Real exp        {best_loss0:.4e}    {r2_0:.6f}     2  (R77 baseline)")
print(f"  Single complex  {best_loss:.4e}    {r2:.6f}     4  (Cotaescu frame)")
print(f"  Real + cplx     {best_loss2:.4e}    {r2_score2:.6f}     6")
print()

# Save
output = {
    "eps_values_used": {str(int(k)): float(v) for k, v in eps_k.items()},
    "char_eigenvalues_q17": [
        {
            "j": j,
            "chi_2": {"re": float(np.exp(2j*np.pi*j/8).real), "im": float(np.exp(2j*np.pi*j/8).imag)},
            "lambda": {"re": float(np.real(char_lams[j])), "im": float(np.imag(char_lams[j]))},
            "abs": float(abs(char_lams[j])),
            "arg_rad": float(np.angle(char_lams[j])),
            "arg_deg": float(np.angle(char_lams[j])/np.pi*180),
        }
        for j in range(8)
    ],
    "model_real": {"A": float(A0), "r": float(r0), "SSE": float(best_loss0), "R2": float(r2_0)},
    "model_complex_single": {
        "A": float(A1), "rho": float(rho1), "theta_rad": float(theta1), "phi_rad": float(phi1),
        "SSE": float(best_loss), "R2": float(r2), "fit_range": "k=2..8",
    },
    "model_real_plus_complex": {
        "B": float(B2), "r": float(r2_real),
        "A": float(A2), "rho": float(rho2), "theta_rad": float(theta2), "phi_rad": float(phi2),
        "SSE": float(best_loss2), "R2": float(r2_score2), "fit_range": "k=1..8",
    },
    "predictions": {
        "model_complex_single_k9": float(model_cplx(9.0, A1, rho1, theta1, phi1)),
        "model_complex_single_k10": float(model_cplx(10.0, A1, rho1, theta1, phi1)),
        "model_rpc_k9": float(model_real_plus_cplx(9.0, B2, r2_real, A2, rho2, theta2, phi2)),
        "model_rpc_k10": float(model_real_plus_cplx(10.0, B2, r2_real, A2, rho2, theta2, phi2)),
    },
}
with open("C:/Collatz/eps_complex_fit_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print("Saved: C:/Collatz/eps_complex_fit_2026_06_01.json")
