"""
watson_phase2c_padeguided.py — PADE-guided fit + multi-saddle test

The k=2..13 data has only 1 sign change → the period inference is weak from data alone.
We use the PADE prediction (period 9.2) and Hadamard tail (ρ≈1.57 at k=13) as STRONG PRIORS
and test whether the multi-saddle Darboux model holds.

Key test: does the data agree with
  ε_k ≈ A · (1.016)^{-k} (slow-mode) + R · (1.57)^{-k} · cos(k * 2π/9.2 + φ)

with the slow-mode dominating eventually?

Also: examine SECOND-DIFFERENCE (Δ²ε_k) structure which damps slow modes.
"""
import sys
import math
import numpy as np
from scipy.optimize import least_squares

sys.stdout.reconfigure(encoding="utf-8")

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
# 2C.1 — Hold ρ_1=1.016, ρ_2=1.57, θ_2=2π/9.2; fit only amplitudes + phase
# ============================================================
print("=" * 70)
print("PHASE 2C.1 — Constrained PADE+Faure model, fit amplitudes only")
print("=" * 70)
print()

rho1 = 1.016
rho2 = 1.57
theta2 = 2*math.pi/9.2
print(f"  Held: ρ_1={rho1:.4f}, ρ_2={rho2:.4f}, θ_2=2π/9.2={theta2:.4f} rad")
print(f"  Free: A, R, φ")
print()

def model_constrained(params, ks):
    A, R, phi = params
    return A * (1.0/rho1)**ks + R * (1.0/rho2)**ks * np.cos(ks*theta2 + phi)

fit_ks = np.array([k for k in ks if k >= 2])
fit_es = np.array([eps_k_data[k] for k in fit_ks])

best = None
for A0 in [-1.0, -0.1, 0.0, 0.1, 1.0]:
    for R0 in [-1.0, -0.1, 0.01, 0.1, 1.0]:
        for phi0 in np.linspace(0, 2*math.pi, 10, endpoint=False):
            try:
                r = least_squares(
                    lambda p, ks, es: model_constrained(p, ks) - es,
                    [A0, R0, phi0], args=(fit_ks, fit_es),
                    bounds=([-100, -100, -2*math.pi], [100, 100, 2*math.pi]),
                    max_nfev=2000,
                )
                rss = np.sum(r.fun**2)
                if best is None or rss < best['rss']:
                    best = dict(rss=rss, params=r.x)
            except Exception:
                pass

A_, R_, phi_ = best['params']
print(f"  Best constrained fit:")
print(f"    A = {A_:.6e}")
print(f"    R = {R_:.6e}")
print(f"    φ = {phi_:.6f} rad")
print(f"    RSS = {best['rss']:.4e}")
print()
print("  Per-k:")
print(f"  {'k':>3} {'ε_k':>14} {'fit':>14} {'A·1.016^-k':>14} {'pair·1.57^-k':>14}")
for k in ks:
    e = eps_k_data[k]
    f_total = model_constrained([A_, R_, phi_], np.array([k]))[0]
    f_slow = A_ * (1.0/rho1)**k
    f_pair = R_ * (1.0/rho2)**k * math.cos(k*theta2 + phi_)
    print(f"  {k:>3} {e:>14.5e} {f_total:>14.5e} {f_slow:>14.5e} {f_pair:>14.5e}")
print()

# ============================================================
# 2C.2 — Fit period along with amplitude; check if period 9.2 emerges or other
# ============================================================
print("=" * 70)
print("PHASE 2C.2 — Fit period (let it find the optimal), constrain ρ_2 to Hadamard tail")
print("=" * 70)
print()

# Hold ρ_2 at the Hadamard tail estimate at k=13 (most accurate inference of leading singularity)
# Let period P, ρ_1, A, R, φ all vary
def model_5p(params, ks):
    A, rho1_v, R, theta_v, phi = params
    return A * (1.0/rho1_v)**ks + R * (1.0/1.57)**ks * np.cos(ks*theta_v + phi)

best5 = None
for A0 in [-0.1, 0.001, 0.1, 1.0]:
    for rho1_0 in [1.016, 1.1, 1.3, 1.5]:
        for R0 in [-1.0, -0.1, 0.01, 0.1, 1.0]:
            for theta0 in [0.3, 0.5, 0.68, 1.0, 1.5]:
                for phi0 in [0.0, math.pi/2, math.pi, -math.pi/2]:
                    try:
                        r = least_squares(
                            lambda p, ks, es: model_5p(p, ks) - es,
                            [A0, rho1_0, R0, theta0, phi0], args=(fit_ks, fit_es),
                            bounds=([-100, 0.5, -100, 0.01, -2*math.pi], [100, 10.0, 100, math.pi, 2*math.pi]),
                            max_nfev=2000,
                        )
                        rss = np.sum(r.fun**2)
                        if best5 is None or rss < best5['rss']:
                            best5 = dict(rss=rss, params=r.x)
                    except Exception:
                        pass

A_, rho1_, R_, theta_, phi_ = best5['params']
print(f"  Best 5-param fit (ρ_2 fixed at 1.57):")
print(f"    A    = {A_:.6e}")
print(f"    ρ_1  = {rho1_:.4f}    (PADE asymptotic: 1.016)")
print(f"    R    = {R_:.6e}")
print(f"    θ_2  = {theta_:.4f} rad, period {2*math.pi/theta_:.3f}  (PADE: 9.2)")
print(f"    φ    = {phi_:.4f}")
print(f"    RSS  = {best5['rss']:.4e}")
print()

# ============================================================
# 2C.3 — Test the slow-mode + transient picture: when does slow mode dominate?
# ============================================================
print("=" * 70)
print("PHASE 2C.3 — When does slow-mode (ρ_1=1.016) dominate over leading (ρ_2=1.57)?")
print("=" * 70)
print()

print("Ratio (ρ_1/ρ_2)^k = (1.016/1.57)^k for various k:")
for k in [13, 20, 30, 50, 100]:
    ratio = (1.016/1.57)**k
    print(f"  k = {k:3d}: (1.016/1.57)^k = {ratio:.4e}")
print()
print("If slow-mode (ρ_1) is the asymptotic, then ε_k ≈ A · 1.016^{-k} for large k.")
print("PADE predicts ρ_1 = 1.016 takes over around k ≈ 20..25, not at k=13.")
print()

# ============================================================
# 2C.4 — Detailed Hadamard tail behavior
# ============================================================
print("=" * 70)
print("PHASE 2C.4 — Hadamard radius trajectory and 'where is the asymptote?'")
print("=" * 70)
print()
print("  k    |ε_k|^{1/k}  → 1/|ε_k|^{1/k}=ρ_k")
for k in sorted(eps_k_data.keys()):
    e = abs(eps_k_data[k])
    if e > 0:
        rad = e ** (1.0/k)
        rho = 1.0/rad
    else:
        rad = float('nan'); rho = float('nan')
    print(f"  {k:3d}   {rad:.5f}    {rho:.5f}")
print()

# Extrapolate trend: linear fit log ρ_k vs 1/k → ρ_∞
ks_used = np.array([10, 11, 12, 13])
rho_used = np.array([1.0 / (abs(eps_k_data[k])**(1.0/k)) for k in ks_used])
inv_k = 1.0 / ks_used.astype(float)
A_mat = np.vstack([np.ones_like(inv_k), inv_k]).T
coef, *_ = np.linalg.lstsq(A_mat, rho_used, rcond=None)
rho_inf_intercept = coef[0]
print(f"  Linear extrapolation of ρ_k vs 1/k from k=10..13:")
print(f"    ρ_k ≈ {coef[0]:.4f} + ({coef[1]:.4f}) · (1/k)")
print(f"    ⟹  ρ_∞ = {rho_inf_intercept:.4f}")
print()
print(f"  PADE predicts ρ_∞ → 1.016 eventually; current k=10..13 trend gives ρ_∞ ≈ {rho_inf_intercept:.4f}")
print(f"  Faure 2009 predicts ρ_∞ → √3 = {math.sqrt(3):.4f}")
print()

# ============================================================
# 2C.5 — Ratio test for asymptotic singularity location
# ============================================================
print("=" * 70)
print("PHASE 2C.5 — Ratio test: |ε_{k+1}/ε_k| → 1/ρ_∞")
print("=" * 70)
print()
print("  k → k+1   ratio |ε_{k+1}/ε_k|    inferred 1/ρ")
for k in [10, 11, 12]:
    r = abs(eps_k_data[k+1]) / abs(eps_k_data[k])
    print(f"  {k} → {k+1}    {r:.4f}              {r:.4f}")
print()
print(f"  Latest ratios suggest 1/ρ ≈ 1.3-2.0, i.e. ρ ≈ 0.5-0.77")
print(f"  ⟹ |ε_k| is GROWING, not decaying! Series has radius of convergence < 1 at k=10..13.")
print()
print("WARNING: f(z) = Σ ε_k z^k may not converge for |z|≤1 in current k=10..13 window.")
print("This is the TRANSIENT regime per PADE_NUMERICAL_DISPOSITION; the asymptotic ρ_∞=1.016")
print("(if real) requires extrapolation BEYOND available data.")
