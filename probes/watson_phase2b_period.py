"""
watson_phase2b_period.py — diagnose the actual oscillation period in ε_k

The pure-pair fit on k=4..13 returned period 8766 (essentially DC + linear drift), which means
the data on 4..13 doesn't contain ENOUGH OSCILLATION CYCLES to identify a period like 9.2.

There's only ONE sign change in the data: between k=9 and k=10. Period > 2*(13-9) = 8 forced if
we want to fit one cycle in the data, but actually we need more cycles to fit. Period 9.2 would
match if we extrapolate backward to k=1 (positive) → k=10 (positive) ~= 9 step cycle.

Let's check WHICH period BEST matches the sign sequence + +  − − − − − − − + + + +.

Approach: for each candidate period P, compute the best-fit phase φ and amplitude envelope A
such that ε_k ≈ A * cos(2π k / P + φ) matches sign. Then plot RSS vs P.
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

print("=" * 70)
print("PHASE 2B — Period scan: ε_k = R * (1/ρ)^k * cos(k * 2π/P + φ)")
print("=" * 70)
print()

def model_period(params, ks, P):
    R, rho, phi = params
    return R * (1.0/rho)**ks * np.cos(ks * 2*math.pi/P + phi)

import warnings
warnings.filterwarnings("ignore")

# Use k=2..13 (skip k=1 which is huge transient)
fit_ks = np.array([k for k in ks if k >= 2])
fit_es = np.array([eps_k_data[k] for k in fit_ks])

results = []
for P in np.linspace(4.0, 30.0, 261):
    best_P = None
    for R0 in [-0.1, -0.01, 0.01, 0.1, 1.0]:
        for rho0 in [0.9, 1.016, 1.2, 1.4, 1.57, 1.732, 2.0]:
            for phi0 in np.linspace(0, 2*math.pi, 8, endpoint=False):
                try:
                    r = least_squares(
                        lambda p, ks, es: model_period(p, ks, P) - es,
                        [R0, rho0, phi0], args=(fit_ks, fit_es),
                        bounds=([-100, 0.5, -2*math.pi], [100, 10.0, 2*math.pi]),
                        max_nfev=1000,
                    )
                    rss = np.sum(r.fun**2)
                    if best_P is None or rss < best_P['rss']:
                        best_P = dict(rss=rss, params=r.x)
                except Exception:
                    pass
    results.append((P, best_P['rss'], best_P['params']))

# Find the best period
results.sort(key=lambda x: x[1])
print("Top 10 best periods (lowest RSS):")
print(f"  {'P':>7} {'rho':>9} {'R':>10} {'phi':>9} {'RSS':>12}")
for P, rss, p in results[:10]:
    R, rho, phi = p
    print(f"  {P:>7.3f} {rho:>9.4f} {R:>10.4f} {phi:>9.4f} {rss:>12.4e}")

# Also fit one specific PADE prediction
P_pade = 9.2
print()
print(f"  At PADE-predicted period P = {P_pade}:")
for r in results:
    if abs(r[0] - P_pade) < 0.06:
        P, rss, params = r
        R, rho, phi = params
        print(f"    P={P:.3f} rho={rho:.4f} R={R:.4f} phi={phi:.4f} RSS={rss:.4e}")
print()

# ============================================================
# Best fit detailed look
# ============================================================
P_best, rss_best, params_best = results[0]
print()
print(f"BEST period: P = {P_best:.4f}, ρ = {params_best[1]:.4f}, RSS = {rss_best:.4e}")
print(f"  θ = 2π/P = {2*math.pi/P_best:.4f} rad")
print()
print("Per-k fit:")
print(f"  {'k':>3} {'ε_k':>14} {'fit':>14} {'resid':>14}")
preds = model_period(params_best, ks, P_best)
for k, e, p in zip(ks, es, preds):
    print(f"  {k:>3} {e:>14.5e} {p:>14.5e} {p-e:>14.5e}")
print()

# ============================================================
# Now also test: maybe single-real-exponential decay fits the asymptotic tail?
# ============================================================
print("=" * 70)
print("PHASE 2B.2 — Single-exponential fit on tail k=8..13")
print("=" * 70)
print()
print("Model: ε_k = A * ρ^{-k} (single decay)")

tail_ks = np.array([k for k in [8, 9, 10, 11, 12, 13]])
tail_es = np.array([eps_k_data[k] for k in tail_ks])
# Note k=9 has sign-change zero crossing; near-zero, may be excluded
print("Using k=10..13 (after sign change):")
fit_ks2 = np.array([10, 11, 12, 13])
fit_es2 = np.array([eps_k_data[k] for k in fit_ks2])
# log linear regression
log_e = np.log(np.abs(fit_es2))
A_mat = np.vstack([np.ones_like(fit_ks2, dtype=float), fit_ks2.astype(float)]).T
coef, *_ = np.linalg.lstsq(A_mat, log_e, rcond=None)
log_A, neg_log_rho = coef
rho_exp = math.exp(-neg_log_rho)
A_exp = math.exp(log_A)
print(f"  Linear log-fit: log|ε_k| = {log_A:.4f} + {neg_log_rho:.4f}*k")
print(f"  ⟹ ε_k ≈ {A_exp:.4e} · ({rho_exp:.4f})^{{-k}}")
print(f"  ρ (singularity) = {1.0/rho_exp:.4f}   if interpreted as decay rate")
print(f"  Alternative: |ε_k|^{{1/k}} = {neg_log_rho:.4f}")
print()

# k=12 and k=13 ratio
r12_13 = eps_k_data[13] / eps_k_data[12]
print(f"  ε_13/ε_12 = {r12_13:.4f}    (period >> 13-12, so close to 1 means slow)")
print(f"  → tail ratio supports decay rate (1/ρ) = {r12_13:.4f}, i.e. ρ = {1.0/r12_13:.4f}")
print()
