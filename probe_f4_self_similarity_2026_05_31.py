"""
f.4 — Bernoulli-convolution / self-similar measure derivation.

Tao-Syracuse self-similarity: X = Z(1+qX') with Z = 2^(-Geom(1/2)) and X' iid copy.

Consequence for D = X-Y (X, Y iid):
  D = (Z-Z') + q(Z·X' - Z'·Y')

Case Z = Z' globally (P = 1/3 by Σ 2^(-2a) for a≥1):
  D = q·Z·D' where D' = X'-Y' iid copy of D.
  v_q(D) = 1 + v_q(D'); σ_m(D) = Z·σ_{m-1}(D').
  Since Z = 2^(-a) is in ⟨2⟩ ⊂ QR(17), χ_2(Z) = +1.
  → χ_2(σ_m(D)) | (Z=Z') = χ_2(σ_{m-1}(D'))
  → Contribution to N_m: (1/3) · N_{m-1}

Hence: N_m = (1/3)·N_{m-1} + N_m^{boundary}
       T_m = (1/3)·T_{m-1} + T_m^{boundary}

In the steady state c(m) → c_∞:
  c_∞ T_m = (1/3) c_∞ T_{m-1} + c_∞^{∂} T_m^{∂}
With T_m^{∂} = T_m - (1/3) T_{m-1}, this is c_∞ = c_∞ trivially (consistent, no constraint).

So self-similarity alone does NOT close c_∞. But it DOES tell us:
  - (1/3) is the natural "iteration constant" of the diagonal contribution
  - The boundary correction is everything else

Now extract z = ρ·e^(iθ) from the VECTOR p_m data (not just scalar c(m)).
With 3 vector samples p_0, p_1, p_2 (each 16-component), and σ↔−σ symmetry,
we have 24 effective equations in 17 (p_∞) + 1 (real ρ) + 1 (real θ) + 8 (B_+, B_-)
unknowns. Overdetermined.

Strategy:
  Joint nonlinear fit of all 16 component triples (p_0(σ), p_1(σ), p_2(σ))
  to model p_m(σ) = p_∞(σ) + B_+(σ)·z^m + B_-(σ)·z̄^m
  with shared z = ρ·e^(iθ) and B_- = conj(B_+) for real-valued p.

Output:
  z to higher precision than scalar c(m) fit
  PSLQ z against Q(i)-algebraic candidates
"""
from __future__ import annotations
import sys
import json
from fractions import Fraction
import numpy as np
from mpmath import mp, mpf, mpc, sqrt, pi, exp, log, pslq, atan2
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40

# Load p_0, p_1, p_2 from previous probe
with open("C:/Collatz/pm_distributions_2026_05_31.json") as f:
    data = json.load(f)
q = data["q"]

p0_frac = {int(s): Fraction(*v) for s, v in data["p_0"].items()}
p1_frac = {int(s): Fraction(*v) for s, v in data["p_1_rational"].items()}
p2_frac = {int(s): Fraction(*v) for s, v in data["p_2_rational"].items()}

# Convert to mpmath for numerical fitting
p0 = {s: mpf(p0_frac[s].numerator) / mpf(p0_frac[s].denominator) for s in range(q)}
p1 = {s: mpf(p1_frac[s].numerator) / mpf(p1_frac[s].denominator) for s in range(q)}
p2 = {s: mpf(p2_frac[s].numerator) / mpf(p2_frac[s].denominator) for s in range(q)}

print("=== p_m exact values loaded ===")
for s in range(1, q):
    print(f"  σ={s:2d}: p_0={float(p0[s]):.10f}  p_1={float(p1[s]):.10f}  p_2={float(p2[s]):.10f}")

# σ↔-σ pairs (8 distinct pair-sums)
pairs = [(1,16), (2,15), (3,14), (4,13), (5,12), (6,11), (7,10), (8,9)]
print("\n=== Pair-sums (verify σ↔-σ symmetry) ===")
for s1, s2 in pairs:
    print(f"  {s1:2d}/{s2:2d}: p_0 diff={float(p0[s1]-p0[s2]):+.2e}, p_1 diff={float(p1[s1]-p1[s2]):+.2e}, p_2 diff={float(p2[s1]-p2[s2]):+.2e}")

# Verify damped-osc model on c(m)
def chi2(x):
    x = x % q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

c0 = sum(chi2(s) * p0[s] for s in range(q))
c1 = sum(chi2(s) * p1[s] for s in range(q))
c2 = sum(chi2(s) * p2[s] for s in range(q))
print(f"\nc(0)={float(c0):.10f}, c(1)={float(c1):.10f}, c(2)={float(c2):.10f}")

# Joint fit: p_m(σ) = p_∞(σ) + 2*Re(B(σ)·z^m)
# Equivalently: p_m(σ) = p_∞(σ) + B_c(σ)·ρ^m·cos(mθ) - B_s(σ)·ρ^m·sin(mθ)
# Three unknowns per σ (p_∞, B_c, B_s), shared (ρ, θ).
#
# For fixed (ρ, θ), the per-σ system is linear:
#   M = [[1, 1, 0], [1, ρcosθ, -ρsinθ], [1, ρ²cos2θ, -ρ²sin2θ]]
# and we get p_∞(σ), B_c(σ), B_s(σ) from M^{-1} · [p_0, p_1, p_2].
#
# The total squared residual is 0 for this fit (we have 3 unknowns and 3 equations per σ),
# so there's no constraint that uniquely determines (ρ, θ).
#
# WAIT — that's wrong: with 3 unknowns per σ and 3 equations per σ, EVERY (ρ, θ) gives a perfect fit.
# So we need more depth (p_3) to constrain (ρ, θ).
#
# OR: we use the SHARED-z constraint differently. The B_c(σ), B_s(σ) for different σ should
# satisfy linear relations (e.g., be eigenvectors of the propagation operator at eigenvalue z).
# That's a softer constraint requiring depth structure.

# Test: at (ρ, θ) from scalar c(m) fit (ρ≈0.076, θ≈arctan(2)), extract per-σ p_∞ and check
# whether reconstructed c_∞ matches reference.

import numpy as np

print("\n=== Per-σ damped-osc decomposition at scalar-fit (ρ, θ) ===")
rho_scalar = 0.076
theta_scalar = float(np.arctan(2))
print(f"  Using ρ={rho_scalar}, θ={theta_scalar:.5f}")

M_fit = np.array([
    [1.0, 1.0,                       0.0                       ],
    [1.0, rho_scalar*np.cos(theta_scalar),       -rho_scalar*np.sin(theta_scalar)       ],
    [1.0, rho_scalar**2*np.cos(2*theta_scalar),  -rho_scalar**2*np.sin(2*theta_scalar)  ],
])

p_inf = {}
B_c = {}
B_s = {}
for s in range(q):
    if s == 0:
        p_inf[s] = 0.0; B_c[s] = 0.0; B_s[s] = 0.0
        continue
    pm_arr = np.array([float(p0[s]), float(p1[s]), float(p2[s])])
    sol = np.linalg.solve(M_fit, pm_arr)
    p_inf[s], B_c[s], B_s[s] = sol[0], sol[1], sol[2]

# Verify sum p_∞ = 1
sum_p_inf = sum(p_inf.values())
print(f"  Sum p_∞(σ) = {sum_p_inf:.10f}")
c_inf_pred = sum(chi2(s) * p_inf[s] for s in range(q))
print(f"  c_∞ predicted = {c_inf_pred:.10f}  vs ref 0.15298912...")

# Now compute z from c(m) fit alone to verify
# c(m) - c_∞ = (Σ χ_2(σ)·B_c(σ))·ρ^m·cos(mθ) - (Σ χ_2(σ)·B_s(σ))·ρ^m·sin(mθ)
B_c_total = sum(chi2(s) * B_c[s] for s in range(q))
B_s_total = sum(chi2(s) * B_s[s] for s in range(q))
print(f"  Σ χ_2·B_c = {B_c_total:.6e}")
print(f"  Σ χ_2·B_s = {B_s_total:.6e}")

# Now SEARCH for (ρ, θ) that best fits depth-3 prediction
# p_3(σ) = p_∞(σ) + ρ³·[B_c(σ)·cos(3θ) - B_s(σ)·sin(3θ)]
# We don't have p_3 yet, but we have c(3) = 0.15300533169151... from prior data
c3_ref = 0.1530053316915140267018082388619123587748
c4_ref = 0.15298870904682116522680466011500498985
c5_ref = 0.15298899941352190530916684536
c_inf_ref = 0.15298912060588517527891674877413229926086222622334

# Fit ρ, θ to scalar c(0..5) using least squares
c_data = [float(c0), float(c1), float(c2), c3_ref, c4_ref, c5_ref]
# Model: c(m) = c_inf + B·ρ^m·cos(mθ+φ)  (3 params: B, ρ, θ, φ → 4 actually... let's fit (B_c, B_s, ρ, θ))
# Treat c_inf as known at ref value, fit on Δ_m = c(m) - c_inf

from scipy.optimize import least_squares

def residual(params):
    rho, theta, B_c_total_p, B_s_total_p = params
    res = []
    for m in range(6):
        Delta_m = c_data[m] - c_inf_ref
        pred = rho**m * (B_c_total_p * np.cos(m*theta) - B_s_total_p * np.sin(m*theta))
        res.append(Delta_m - pred)
    return res

p0_guess = [0.076, np.arctan(2), -0.001, -0.002]
result = least_squares(residual, p0_guess, ftol=1e-15, xtol=1e-15)
rho_fit, theta_fit, Bc_fit, Bs_fit = result.x
print(f"\n=== Scalar c(0..5) fit (using c_inf at 50-digit ref) ===")
print(f"  ρ = {rho_fit:.10f}")
print(f"  θ = {theta_fit:.10f} rad = {np.degrees(theta_fit):.6f}°")
print(f"  B_c = {Bc_fit:.6e}, B_s = {Bs_fit:.6e}")
print(f"  arctan(2) = {np.arctan(2):.10f}")
print(f"  θ - arctan(2) = {theta_fit - np.arctan(2):+.6e}")
print(f"  Residual: {result.cost:.2e}")

# Now check z = ρ·e^(iθ) precision
z_re = rho_fit * np.cos(theta_fit)
z_im = rho_fit * np.sin(theta_fit)
print(f"  z = {z_re:+.10f} + {z_im:+.10f}i")
print(f"  |z|² = {rho_fit**2:.10f}")
print(f"  z·z̄ + z + z̄ = {rho_fit**2 + 2*z_re:.10f}")

# PSLQ z against Q(i)-algebraic basis
print(f"\n=== PSLQ z (re,im) against Q(i)-algebraic basis at fit precision ===")
mp.dps = 30
z_re_mp = mpf(z_re)
z_im_mp = mpf(z_im)

# Algebraic basis: 1, 1/n, 1/(a+bi) for small a,b
basis_names = ["1", "1/2", "1/3", "1/4", "1/5", "1/6", "1/8", "1/16", "1/17",
               "1/29", "1/30", "1/31", "1/sqrt(5)", "1/sqrt(17)", "1/sqrt(85)",
               "1/(1+2i).re", "1/(1+2i).im", "1/(2+i).re", "1/(2+i).im",
               "1/(4+i).re", "1/(4+i).im", "1/(1+4i).re", "1/(1+4i).im",
               "1/(2+3i).re", "1/(2+3i).im", "1/(3+2i).re", "1/(3+2i).im"]
basis_vals = []
for s in [mpf(1), mpf(1)/2, mpf(1)/3, mpf(1)/4, mpf(1)/5, mpf(1)/6, mpf(1)/8, mpf(1)/16, mpf(1)/17]:
    basis_vals.append(s)
for n in [29, 30, 31]:
    basis_vals.append(mpf(1)/n)
for n in [5, 17, 85]:
    basis_vals.append(mpf(1)/sqrt(mpf(n)))
# Gaussian integer reciprocals
for a, b in [(1,2), (2,1), (4,1), (1,4), (2,3), (3,2)]:
    denom_sq = mpf(a*a + b*b)
    re_part = mpf(a) / denom_sq
    im_part = -mpf(b) / denom_sq
    basis_vals.append(re_part)
    basis_vals.append(im_part)

# PSLQ z_re against basis
print(f"  z_re = {z_re_mp}")
for tol_exp in [4, 6, 8]:
    tol = mpf(10) ** (-tol_exp)
    rel = pslq([z_re_mp] + basis_vals, tol=tol, maxcoeff=50)
    if rel is not None and rel[0] != 0:
        terms = [f"({rel[0]:+d})·z_re"] + [f"({c:+d})·{n}" for c, n in zip(rel[1:], basis_names) if c != 0]
        print(f"  tol=10^-{tol_exp}: {' '.join(terms[:8])}")
        break
else:
    print(f"  No clean PSLQ relation for z_re")

print(f"  z_im = {z_im_mp}")
for tol_exp in [4, 6, 8]:
    tol = mpf(10) ** (-tol_exp)
    rel = pslq([z_im_mp] + basis_vals, tol=tol, maxcoeff=50)
    if rel is not None and rel[0] != 0:
        terms = [f"({rel[0]:+d})·z_im"] + [f"({c:+d})·{n}" for c, n in zip(rel[1:], basis_names) if c != 0]
        print(f"  tol=10^-{tol_exp}: {' '.join(terms[:8])}")
        break
else:
    print(f"  No clean PSLQ relation for z_im")

# Try z as (a+bi)/(c+di) for small (a,b,c,d)
print(f"\n=== Brute-force: z = (a+bi)/(c+di) over small integers ===")
hits = []
for a in range(-5, 6):
    for b in range(-5, 6):
        for c in range(1, 40):
            for d in range(-40, 41):
                if c == 0 and d == 0: continue
                denom = c*c + d*d
                if denom == 0: continue
                z_test_re = (a*c + b*d) / denom
                z_test_im = (b*c - a*d) / denom
                # Compare to fit
                err = abs(z_test_re - z_re) + abs(z_test_im - z_im)
                if err < 1e-4:
                    hits.append((a, b, c, d, z_test_re, z_test_im, err))
hits.sort(key=lambda x: x[6])
print(f"  Top 10 closest (a+bi)/(c+di) candidates:")
for a, b, c, d, zr, zi, err in hits[:10]:
    print(f"    ({a:+d}{b:+d}i)/({c:+d}{d:+d}i) = {zr:+.6e} + {zi:+.6e}i, err={err:.2e}")

# Try z = sqrt(p+qi)/N
print(f"\n=== Brute-force: z² candidate (since z²+z̄² and z·z̄ might be cleaner) ===")
z2_re = z_re**2 - z_im**2
z2_im = 2*z_re*z_im
print(f"  z² = {z2_re:+.10f} + {z2_im:+.10f}i")
print(f"  |z²| = {(z2_re**2 + z2_im**2)**0.5:.10f}")
print(f"  z·z̄ (|z|²) = {z_re**2 + z_im**2:.10f}")
print(f"  z + z̄ (2·Re(z)) = {2*z_re:.10f}")

# Look for clean rationals for |z|² and Re(z)
print(f"\n  |z|² candidates:")
for d in [16, 17, 64, 128, 256, 289, 17*256, 16*17, 5*17, 64*17, 169, 13*17]:
    n = round((z_re**2 + z_im**2) * d)
    err = abs(n/d - (z_re**2 + z_im**2))
    if err < 1e-5:
        print(f"    |z|² ≈ {n}/{d} = {n/d:.10f}, err={err:.2e}")
print(f"\n  2·Re(z) candidates:")
for d in [16, 17, 64, 128, 256, 289, 17*16, 5*17]:
    n = round(2 * z_re * d)
    err = abs(n/d - 2*z_re)
    if err < 1e-5:
        print(f"    2·Re(z) ≈ {n}/{d} = {n/d:.10f}, err={err:.2e}")

print("\n=== Done with f.4 ===")
