"""
probe_damped_oscillation_fit_2026_05_30.py

Treating (Δ_1, Δ_2, Δ_3, Δ_4, Δ_5) as a damped oscillation around an asymptote,
fit the model:
    Δ_m = A + B · ρ^m · cos(m·θ + φ)

Equivalently, with complex eigenvalue z = ρ·e^{iθ}:
    Δ_m = A + Re(C · z^m),   C = (B/2)·e^{iφ}

This satisfies the linear recurrence (after subtracting A):
    Δ_{m+2} = u · Δ_{m+1} - v · Δ_m + K
where u = 2·ρ·cos(θ), v = ρ², K = A·(1 - u + v).

5 data points → 3 linear equations in (u, v, K), exact solution. Then:
    ρ = √v
    θ = arccos(u / (2ρ))
    A = K / (1 - u + v)
    B, φ from Δ_1, Δ_2 initial conditions.

Checks:
  - ρ vs {1/√17, 1/16, 1/17, 1/√q²-q, 2^{-k}}
  - θ vs {π, π/2, 2π/q, π·(1-1/q)}
  - c_∞_model = A; compare to Shanks-extrapolated c_∞
"""
from __future__ import annotations
import sys, gc
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

q = 17

def offset_distribution(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a)
        p = (p * inv2) % N
    P_U /= P_U.sum()
    u_support = np.nonzero(P_U)[0]; u_weight = P_U[u_support]
    P_S = P_U.copy()
    for j in range(n - 1, 0, -1):
        v_idx = (1 + q * arange) % N
        P_V = np.zeros(N, dtype=np.float64)
        np.add.at(P_V, v_idx, P_S)
        P_new = np.zeros(N, dtype=np.float64)
        for u, w in zip(u_support, u_weight):
            P_new[(u * arange) % N] += w * P_V
        P_S = P_new
    return P_S

def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq - 1) // 2, qq) == 1 else -1)

chi_q = np.array([legendre(x, q) for x in range(q)], dtype=np.float64)

def c_at_n_m(n, m, A_MAX=200):
    N = q ** n
    P = offset_distribution(q, n, A_MAX)
    mu = np.fft.fft(P)
    a2 = np.abs(mu) ** 2
    PD = np.real(np.fft.ifft(a2))
    qm = q ** m
    qmp1 = q ** (m + 1)
    d_vals = np.arange(N)
    mask = (d_vals % qm == 0) & (d_vals % qmp1 != 0)
    leading_digit = (d_vals[mask] // qm) % q
    chi_vals = chi_q[leading_digit]
    num = (PD[mask] * chi_vals).sum()
    den = PD[mask].sum()
    return float(num / den)

# Use cached/precomputed c(m) for m=1..4; compute c(5), c(6) fresh
print("Recomputing c(m) at A_MAX=200 for high precision...")
import time
c_vals = {}
for m in (1, 2, 3, 4, 5):
    n = m + 1
    t0 = time.time()
    print(f"  m={m}, n={n}, N={q**n}, computing...", flush=True)
    c_vals[m] = c_at_n_m(n=n, m=m, A_MAX=200)
    print(f"    c({m}) = {c_vals[m]:.16f}   ({time.time()-t0:.1f}s)")
    gc.collect()

# Optional: c(6) for cross-check at n=7. N=17^7=410M — too big in float64. Skip.

c0 = 19 / 127
deltas = {m: c_vals[m] - c0 for m in c_vals}
print(f"\n=== Δ_m sequence ===")
for m in sorted(deltas):
    print(f"  Δ_{m} = {deltas[m]:+.16e}")

# Solve the recurrence: 3 equations in (u, v, K)
# Δ_{m+2} = u·Δ_{m+1} - v·Δ_m + K, for m=1,2,3
# i.e. Δ_3 = u·Δ_2 - v·Δ_1 + K
#      Δ_4 = u·Δ_3 - v·Δ_2 + K
#      Δ_5 = u·Δ_4 - v·Δ_3 + K
# Matrix form: [Δ_2, -Δ_1, 1; Δ_3, -Δ_2, 1; Δ_4, -Δ_3, 1] [u;v;K] = [Δ_3; Δ_4; Δ_5]
M = np.array([
    [deltas[2], -deltas[1], 1.0],
    [deltas[3], -deltas[2], 1.0],
    [deltas[4], -deltas[3], 1.0],
])
rhs = np.array([deltas[3], deltas[4], deltas[5]])
u, v, K = np.linalg.solve(M, rhs)
print(f"\n=== Linear recurrence fit ===")
print(f"  u = 2·ρ·cos(θ) = {u:.16f}")
print(f"  v = ρ²        = {v:.16e}")
print(f"  K = A·(1-u+v) = {K:.16e}")

# Derived quantities
if v < 0:
    print(f"  WARNING: v = {v} is NEGATIVE → not a damped oscillation (no real ρ).")
    print(f"  Model may not fit. Check for two real eigenvalues instead.")
    # Two real eigenvalues case: z^2 - u·z + v = 0, z = (u ± √(u²-4v))/2
    disc = u*u - 4*v
    print(f"  discriminant u²-4v = {disc:.6e}")
    if disc > 0:
        z1 = (u + np.sqrt(disc)) / 2
        z2 = (u - np.sqrt(disc)) / 2
        print(f"  z1 = {z1:.10f}   z2 = {z2:.10f}")
    rho = abs(v)**0.5
    theta = float('nan')
else:
    rho = np.sqrt(v)
    if abs(u / (2*rho)) > 1:
        print(f"  WARNING: u/(2ρ) = {u/(2*rho):.6f} > 1 → eigenvalues real.")
        z1 = (u + np.sqrt(u*u - 4*v)) / 2
        z2 = (u - np.sqrt(u*u - 4*v)) / 2
        print(f"  z1 = {z1:.16f}   z2 = {z2:.16f}")
        theta = 0.0 if u > 0 else np.pi
    else:
        theta = np.arccos(u / (2*rho))
        z = rho * np.exp(1j * theta)
        print(f"  ρ = {rho:.16f}   θ = {theta:.16f} rad = {np.degrees(theta):.6f}°")
        print(f"  complex eigenvalue z = ρ·e^{{iθ}} = {z}")

# Asymptote A
denom = 1.0 - u + v
if abs(denom) > 1e-30:
    A = K / denom
    print(f"  A (asymptote) = {A:.16e}")
    print(f"  c_∞ from model = 19/127 + A = {c0 + A:.16f}")
else:
    print(f"  WARNING: 1-u+v = {denom:.6e} ≈ 0 → A degenerate.")
    A = float('nan')

# Recover B, φ from Δ_1, Δ_2 (assuming complex eigenvalues)
if v > 0 and abs(u / (2*rho)) <= 1:
    # Δ_m - A = B·ρ^m·cos(mθ + φ)
    # → ε_m = B·ρ^m·cos(mθ + φ) = Re(B·e^{iφ} · z^m)
    # Set C = B·e^{iφ} / 2 (complex). Then Δ_m - A = Re((C·z^m) + (C̄·z̄^m)) = 2·Re(C·z^m) — wait.
    # Standard: Δ_m - A = α·cos(mθ)·ρ^m + β·sin(mθ)·ρ^m  with α = B·cos(φ), β = -B·sin(φ).
    e1 = deltas[1] - A
    e2 = deltas[2] - A
    # System: e_m = α·cos(mθ)·ρ^m + β·sin(mθ)·ρ^m
    # m=1: e1 = ρ·(α·cos(θ) - β·sin(-θ)) — let me be careful
    # cos(1·θ+φ) = cos(θ)cos(φ) - sin(θ)sin(φ); with α=B·cos(φ), β=-B·sin(φ):
    # B·cos(mθ+φ) = α·cos(mθ) + β·sin(mθ)
    # Yes.
    M2 = np.array([
        [rho * np.cos(theta),     rho * np.sin(theta)],
        [rho**2 * np.cos(2*theta), rho**2 * np.sin(2*theta)],
    ])
    rhs2 = np.array([e1, e2])
    alpha, beta = np.linalg.solve(M2, rhs2)
    B = np.sqrt(alpha**2 + beta**2)
    phi = np.arctan2(-beta, alpha)
    print(f"  B = {B:.10e}   φ = {phi:.10f} rad = {np.degrees(phi):.6f}°")

    # Verify the model on Δ_3, Δ_4, Δ_5
    print(f"\n=== Model verification ===")
    for m in (1, 2, 3, 4, 5):
        model_val = A + B * rho**m * np.cos(m*theta + phi)
        actual = deltas[m]
        print(f"  m={m}: actual={actual:+.12e}  model={model_val:+.12e}  diff={actual-model_val:+.6e}")

# Check ρ against algebraic candidates
print(f"\n=== ρ candidate check ===")
candidates_rho = {
    "1/√17":   1.0 / np.sqrt(17),
    "1/16":    1.0 / 16,
    "1/17":    1.0 / 17,
    "1/√255":  1.0 / np.sqrt(255),
    "1/√256":  1.0 / 16,
    "1/√(q²)":  1.0 / 17,
    "1/15":    1.0 / 15,
    "2^{-4}":  1.0 / 16,
    "2^{-5}":  1.0 / 32,
    "1/(q-1)":  1.0 / 16,
    "1/(q+1)":  1.0 / 18,
    "2/q":     2.0 / 17,
    "1/2^4·q": 1.0 / (16*17),
}
if not np.isnan(rho):
    for name, val in candidates_rho.items():
        diff = abs(rho - val)
        marker = "  ★" if diff < 1e-3 else ""
        print(f"  ρ={rho:.10f} vs {name}={val:.10f}  diff={diff:.6e}{marker}")

# Check θ against simple angle candidates
print(f"\n=== θ candidate check ===")
if not np.isnan(theta):
    candidates_theta = {
        "π":       np.pi,
        "π/2":     np.pi / 2,
        "2π/3":    2*np.pi / 3,
        "3π/4":    3*np.pi / 4,
        "5π/6":    5*np.pi / 6,
        "π·(1-1/q)": np.pi * (1 - 1.0/17),
        "π - 2π/q":  np.pi - 2*np.pi/17,
        "2π/q·8":    2*np.pi / 17 * 8,
        "2π · 7/15":  2*np.pi * 7 / 15,
        "arccos(-1/3)": np.arccos(-1.0/3),
    }
    for name, val in candidates_theta.items():
        diff = abs(theta - val)
        marker = "  ★" if diff < 1e-3 else ""
        print(f"  θ={theta:.10f} vs {name}={val:.10f}  diff={diff:.6e}{marker}")

# Algebraic candidates for v=ρ² and u=2ρcosθ
print(f"\n=== Algebraic candidates for v=ρ² ===")
cand_v = {
    "1/17":    1.0/17,
    "1/255":   1.0/255,
    "1/256":   1.0/256,
    "1/(15·17)": 1.0/(15*17),
    "1/(16·17)": 1.0/(16*17),
    "1/289":   1.0/289,
    "2/255":   2.0/255,
    "1/127":   1.0/127,
    "1/(2·127)": 1.0/(2*127),
    "1/40":    1.0/40,
}
for name, val in cand_v.items():
    diff = abs(v - val)
    marker = "  ★" if diff < 1e-4 else ""
    print(f"  v={v:.10e} vs {name}={val:.10e}  diff={diff:.6e}{marker}")

print(f"\n=== Algebraic candidates for u=2ρcosθ ===")
cand_u = {
    "-1/8":  -1.0/8,
    "-1/16": -1.0/16,
    "-2/17": -2.0/17,
    "0":      0.0,
    "-1/15": -1.0/15,
    "-1/(2·17)": -1.0/(2*17),
}
for name, val in cand_u.items():
    diff = abs(u - val)
    marker = "  ★" if diff < 1e-4 else ""
    print(f"  u={u:.10f} vs {name}={val:.10f}  diff={diff:.6e}{marker}")

# Compute c_∞ from model and compare to Shanks
print(f"\n=== c_∞ comparison ===")
print(f"  c_∞ from model   = {c0 + A:.16f}")
# Shanks on c_vals 2,3,4 and 3,4,5
def shanks(a, b, c):
    den = c - 2*b + a
    return c - (c - b)**2 / den if abs(den) > 1e-30 else c
c_inf_shanks_234 = shanks(c_vals[2], c_vals[3], c_vals[4])
c_inf_shanks_345 = shanks(c_vals[3], c_vals[4], c_vals[5])
print(f"  c_∞ Shanks 2,3,4 = {c_inf_shanks_234:.16f}")
print(f"  c_∞ Shanks 3,4,5 = {c_inf_shanks_345:.16f}")
print(f"  model - Shanks 345 = {(c0+A) - c_inf_shanks_345:+.6e}")
