"""
watson_phase2_asymptotic.py — Phase 2 of WATSON probe.

Execute the Darboux multi-saddle asymptotic on ε_k k=2..13 data:
  ε_k ~ A · ρ_1^{-k} + 2|C| · ρ_2^{-k} · cos(k θ_2 + φ)

Fit (A, ρ_1, |C|, ρ_2, θ_2, φ) by nonlinear least squares to ε_k k=2..13 data.

Compare:
- PADE prediction: ρ_2 ≈ 1.57 (transient), θ_2 ≈ 0.68 rad (period 9.2), ρ_1 ≈ 1.016 (asymptotic).
- Faure prediction: ρ ≈ √3 ≈ 1.732 (semiclassical limit).

Also compute Hadamard radius lim sup |ε_k|^{1/k} on tail k=8..13 to confirm convergence radius
of f(z) = Σ ε_k z^k, and cross-check with the multi-saddle fit.

Output: WATSON_ASYMPTOTIC.md with the fitted parameters + residuals + comparison.
"""
import sys
import os
import math
import numpy as np
from scipy.optimize import least_squares, minimize

sys.stdout.reconfigure(encoding="utf-8")

# ε_k data k=1..13 from PADE_NUMERICAL_DATA.md
eps_k_data = {
    1:  +2.0000000000e-01,
    2:  +9.5238095238e-03,
    3:  -5.0919863259e-03,
    4:  -2.4522582483e-03,
    5:  -1.1517469151e-03,
    6:  -4.9790566522e-04,
    7:  -1.1752368304e-03,
    8:  -7.4554636729e-04,
    9:  -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
    12: +2.2747137206e-03,
    13: +2.9482473172e-03,
}

ks = np.array(sorted(eps_k_data.keys()))
es = np.array([eps_k_data[k] for k in ks])

# ============================================================
# Phase 2.1 — Hadamard radius from tail
# ============================================================
print("=" * 70)
print("PHASE 2.1 — Hadamard radius from tail ε_k^{1/k}")
print("=" * 70)
print()
print("  k     |eps_k|        |eps_k|^(1/k)   inferred rho = 1/that")
for k in [8, 9, 10, 11, 12, 13]:
    e = abs(eps_k_data[k])
    if e > 0:
        rad = e ** (1.0/k)
        rho = 1.0 / rad
    else:
        rad = float('nan')
        rho = float('nan')
    print(f"  {k:3d}   {e:.4e}     {rad:.4f}          {rho:.4f}")

print()
print("PADE_NUMERICAL_DISPOSITION reported Hadamard at n=13: |z|≈1.57.")
print()

# ============================================================
# Phase 2.2 — Pure complex-pair fit: ε_k = 2 |C| ρ^{-k} cos(k θ + φ)
# ============================================================
print("=" * 70)
print("PHASE 2.2 — Pure complex-pair fit on k=4..13 (skip near-zero k=9)")
print("=" * 70)
print()
print("Model: ε_k = R * (1/ρ)^k * cos(k θ + φ)")
print("  Free params: R (=2|C|), ρ, θ, φ")
print()

def model_pair(params, ks):
    R, rho, theta, phi = params
    return R * (1.0/rho)**ks * np.cos(ks*theta + phi)

def residuals_pair(params, ks, es):
    return model_pair(params, ks) - es

# Use k=4..13 — skip k=1,2 (outliers from initial transient) and k=3 (sign flip)
# Try several initial conditions including PADE prediction
best = None
for rho0 in [1.016, 1.3, 1.57, 1.732, 2.0]:
    for theta0 in [0.5, 0.68, 1.0, 1.5]:
        for R0 in [0.001, 0.01, 0.1, 1.0]:
            for phi0 in [0.0, math.pi/2, math.pi, -math.pi/2]:
                try:
                    fit_ks = np.array([k for k in ks if k >= 4])
                    fit_es = np.array([eps_k_data[k] for k in fit_ks])
                    result = least_squares(
                        residuals_pair, [R0, rho0, theta0, phi0],
                        args=(fit_ks, fit_es),
                        bounds=([1e-6, 0.5, 0.0, -2*math.pi], [10.0, 10.0, math.pi, 2*math.pi]),
                        max_nfev=2000,
                    )
                    rss = np.sum(result.fun**2)
                    if best is None or rss < best['rss']:
                        best = dict(rss=rss, params=result.x, fit_ks=fit_ks, fit_es=fit_es)
                except Exception:
                    pass

R_fit, rho_fit, theta_fit, phi_fit = best['params']
print(f"  Best pair fit (k=4..13):")
print(f"    R     (=2|C|) = {R_fit:.6f}")
print(f"    ρ            = {rho_fit:.6f}     (PADE: ~1.57; Faure: ~1.732)")
print(f"    θ            = {theta_fit:.6f} rad  (PADE: ~0.68; period 2π/θ = {2*math.pi/theta_fit:.3f})")
print(f"    φ            = {phi_fit:.6f} rad")
print(f"    RSS          = {best['rss']:.4e}")
print()
print("  Per-k fit:")
print(f"    {'k':>3} {'ε_k':>14} {'fit':>14} {'resid':>14} {'rel':>10}")
preds = model_pair(best['params'], ks)
for k, e, p in zip(ks, es, preds):
    rel = (p - e)/abs(e) if abs(e) > 0 else float('nan')
    marker = " " if k >= 4 else "*"
    print(f"  {marker}{k:>3} {e:>14.5e} {p:>14.5e} {p-e:>14.5e} {rel:>10.3f}")
print()
print("  * = not in fit window (k < 4)")

# ============================================================
# Phase 2.3 — Two-component fit: real slow-mode + complex pair
# ============================================================
print()
print("=" * 70)
print("PHASE 2.3 — Two-component fit: real ρ_1 + complex pair ρ_2 e^{±iθ_2}")
print("=" * 70)
print()
print("Model: ε_k = A * (1/ρ_1)^k + R * (1/ρ_2)^k * cos(k θ_2 + φ)")
print("  Free params: A, ρ_1, R, ρ_2, θ_2, φ")
print()

def model_two(params, ks):
    A, rho1, R, rho2, theta2, phi = params
    return A * (1.0/rho1)**ks + R * (1.0/rho2)**ks * np.cos(ks*theta2 + phi)

def residuals_two(params, ks, es):
    return model_two(params, ks) - es

# Multi-start
best2 = None
fit_ks = np.array([k for k in ks if k >= 2])
fit_es = np.array([eps_k_data[k] for k in fit_ks])

starts = []
for A0 in [-0.1, 0.0, 0.1, 1.0]:
    for rho1_0 in [1.016, 1.1, 1.3, 1.5, 1.732]:
        for R0 in [0.001, 0.01, 0.1, 1.0]:
            for rho2_0 in [1.3, 1.57, 1.732, 2.0]:
                for theta2_0 in [0.5, 0.68, 1.0]:
                    for phi0 in [0.0, math.pi/2, math.pi]:
                        starts.append([A0, rho1_0, R0, rho2_0, theta2_0, phi0])

print(f"  Trying {len(starts)} multi-start configurations...")
import warnings
warnings.filterwarnings("ignore")
for s in starts:
    try:
        result = least_squares(
            residuals_two, s, args=(fit_ks, fit_es),
            bounds=([-10.0, 0.5, -10.0, 0.5, 0.01, -2*math.pi],
                    [10.0, 10.0, 10.0, 10.0, math.pi, 2*math.pi]),
            max_nfev=5000,
        )
        rss = np.sum(result.fun**2)
        if best2 is None or rss < best2['rss']:
            best2 = dict(rss=rss, params=result.x)
    except Exception:
        pass

A_fit, rho1_fit, R_fit2, rho2_fit, theta2_fit, phi_fit2 = best2['params']
print(f"  Best two-component fit (k=2..13):")
print(f"    A          = {A_fit:.6f}")
print(f"    ρ_1 (real) = {rho1_fit:.6f}    (PADE asymp: ~1.016)")
print(f"    R   (=2|C|) = {R_fit2:.6f}")
print(f"    ρ_2 (pair) = {rho2_fit:.6f}    (PADE: ~1.57; Faure: ~1.732)")
print(f"    θ_2        = {theta2_fit:.6f} rad  (PADE: ~0.68; period {2*math.pi/theta2_fit:.3f})")
print(f"    φ          = {phi_fit2:.6f} rad")
print(f"    RSS        = {best2['rss']:.4e}")
print()
print("  Per-k fit (two-component model):")
print(f"    {'k':>3} {'ε_k':>14} {'fit':>14} {'resid':>14} {'rel':>10}")
preds2 = model_two(best2['params'], ks)
for k, e, p in zip(ks, es, preds2):
    rel = (p - e)/abs(e) if abs(e) > 0 else float('nan')
    print(f"   {k:>3} {e:>14.5e} {p:>14.5e} {p-e:>14.5e} {rel:>10.3f}")
print()

# ============================================================
# Phase 2.4 — Constrained fit: use PADE/Faure values for ρ
# ============================================================
print("=" * 70)
print("PHASE 2.4 — Predicted-values check: hold ρ_2 = √3 (Faure), let others fit")
print("=" * 70)
print()

def model_pair_rho_fixed(params, ks, rho_fixed):
    R, theta, phi = params
    return R * (1.0/rho_fixed)**ks * np.cos(ks*theta + phi)

for rho_test_name, rho_test in [
    ("Faure √3", math.sqrt(3.0)),
    ("PADE 1.57", 1.57),
    ("ρ=2 (R77.6 transient)", 2.0),
]:
    best_r = None
    fit_ks = np.array([k for k in ks if k >= 4])
    fit_es = np.array([eps_k_data[k] for k in fit_ks])
    for R0 in [0.01, 0.1, 1.0, 10.0]:
        for theta0 in [0.5, 0.68, 1.0, 1.5]:
            for phi0 in [0.0, math.pi/2, math.pi, -math.pi/2]:
                try:
                    result = least_squares(
                        lambda p, ks, es: model_pair_rho_fixed(p, ks, rho_test) - es,
                        [R0, theta0, phi0], args=(fit_ks, fit_es),
                        bounds=([-100.0, 0.01, -2*math.pi], [100.0, math.pi, 2*math.pi]),
                        max_nfev=2000,
                    )
                    rss = np.sum(result.fun**2)
                    if best_r is None or rss < best_r['rss']:
                        best_r = dict(rss=rss, params=result.x)
                except Exception:
                    pass
    R_, theta_, phi_ = best_r['params']
    print(f"  Held {rho_test_name} = {rho_test:.4f}:")
    print(f"    R   = {R_:.6f}")
    print(f"    θ   = {theta_:.6f} rad  (period 2π/θ = {2*math.pi/theta_:.3f})")
    print(f"    φ   = {phi_:.6f} rad")
    print(f"    RSS = {best_r['rss']:.4e}")
    print()

# ============================================================
# Phase 2.5 — Test the Darboux prediction at k=8..13 — does it match?
# ============================================================
print("=" * 70)
print("PHASE 2.5 — Numerical test of Darboux prediction at k=8..13")
print("=" * 70)
print()
print("Using best two-component fit, predict ε_k at k=8..13 and compute errors:")
print()
print(f"    {'k':>3} {'ε_k actual':>16} {'ε_k Darboux':>16} {'rel err':>12}")
for k in [8, 9, 10, 11, 12, 13]:
    actual = eps_k_data[k]
    pred = model_two(best2['params'], np.array([k]))[0]
    rel = (pred - actual)/abs(actual) if abs(actual) > 0 else float('nan')
    print(f"    {k:>3} {actual:>16.6e} {pred:>16.6e} {rel:>12.3%}")
print()

# ============================================================
# Phase 2.6 — Save outputs
# ============================================================
out = "C:/Collatz/watson_phase2_output.txt"
with open(out, 'w', encoding='utf-8') as f:
    f.write("WATSON Phase 2 — asymptotic fit outputs\n\n")
    f.write(f"Pair-only fit (k=4..13): rho={rho_fit:.4f}, theta={theta_fit:.4f}, period={2*math.pi/theta_fit:.3f}, RSS={best['rss']:.4e}\n")
    f.write(f"Two-component fit (k=2..13): rho_1={rho1_fit:.4f}, rho_2={rho2_fit:.4f}, theta_2={theta2_fit:.4f}, RSS={best2['rss']:.4e}\n")
print(f"  Saved {out}")
