"""
Build the full 8x8 class-resolved operator T from the level-jump sequence
v_1, v_2, ..., v_n (8-dim class-resolved bilinear-moment vectors).

Test whether the fitted complex eigenvalue (rho = 0.190, theta = 41°) from
eps_k data appears in T_8x8's spectrum.

v_n components:
  [P^++(1), P^++(2), P^--(1), P^--(2), Re P^+-(1), Im P^+-(1), Re P^+-(2), Im P^+-(2)]

If v_{n+1} = T · v_n holds (assuming T is n-independent in the asymptotic regime),
then with 7 transitions (n=1..7 -> n=2..8) we get 7*8 = 56 equations for 64
unknowns. Underdetermined by 8, so we add:
  - Eigenvalue 1 on the (1, 4, 1, 4, 0, 0, 0, 0) direction (mass conservation
    from R77's T_diag eigenvector (1, 4) on (P_+, P_-) extended trivially).
  - Class swap symmetry: components 0<->1, 2<->3, 4<->6, 5<->7 commute with T.
  - Real-valuedness on the v subspace (T has real entries).

With those constraints, T is determined and we can diagonalize.
"""
from __future__ import annotations
import sys, os, time, json, cmath
from fractions import Fraction
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:/Collatz")

from T_lead_operator import (
    build_markov_rational, stationary_rational,
    char_funcs_class_resolved, compute_P_vector,
)

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)

# Compute v_n for n = 1..N_MAX
N_MAX = 6  # n=6 => 486 coprime states, ~30 sec stationary

def P_to_vec(P):
    return np.array([
        float(P[('+', '+', 1)].real),
        float(P[('+', '+', 2)].real),
        float(P[('-', '-', 1)].real),
        float(P[('-', '-', 2)].real),
        float(P[('+', '-', 1)].real),
        float(P[('+', '-', 1)].imag),
        float(P[('+', '-', 2)].real),
        float(P[('+', '-', 2)].imag),
    ])

print(f"[1/3] Computing v_n at n=1..{N_MAX}...", flush=True)
t0 = time.time()
v_list = []
for k in range(1, N_MAX + 1):
    tk = time.time()
    K, coprime = build_markov_rational(k)
    pi_q = stationary_rational(K)
    mu_plus, mu_minus = char_funcs_class_resolved(pi_q, coprime, k)
    P = compute_P_vector(mu_plus, mu_minus, coprime, k)
    v = P_to_vec(P)
    v_list.append(v)
    print(f"  n={k}: |coprime|={len(coprime)}, v_n={[f'{x:+.6f}' for x in v]}, time={time.time()-tk:.1f}s", flush=True)
print(f"  total: {time.time()-t0:.1f}s", flush=True)
print()

# Stack as columns: V[:, i] = v_{i+1}, i=0..N_MAX-1
V = np.array(v_list).T  # shape (8, N_MAX)
print(f"V shape: {V.shape}")
print(f"v_n norms: {[np.linalg.norm(v) for v in v_list]}")
print()

# v_{n+1} = T · v_n  =>  V[:, 1:] = T · V[:, :-1]
# T_8x8 has 64 unknowns; we have 8*(N_MAX-1) equations.
# For N_MAX = 6: 5 transitions * 8 = 40 equations. Underdetermined.
# For N_MAX = 9: 8 transitions * 8 = 64 — exactly determined (if no noise).
print(f"[2/3] Solving for T via least squares: V[:, 1:] = T · V[:, :-1]", flush=True)
X = V[:, :-1]  # 8 x (N_MAX-1)
Y = V[:, 1:]   # 8 x (N_MAX-1)
print(f"  X shape {X.shape}, Y shape {Y.shape}")
# T · X = Y => T = Y · X^+
# Using pseudoinverse (will give min-norm T when underdetermined)
T_emp = Y @ np.linalg.pinv(X)
print(f"  T_emp condition: rank(X) = {np.linalg.matrix_rank(X)}, rank(T_emp) = {np.linalg.matrix_rank(T_emp)}")
print()

# Sanity check: how well does T_emp reproduce the transitions?
err = Y - T_emp @ X
print(f"  Residual ||Y - T·X||: {np.linalg.norm(err):.4e}")
print(f"  Residual per-transition norms:")
for i in range(X.shape[1]):
    print(f"    n={i+1} -> n={i+2}: err norm = {np.linalg.norm(err[:, i]):.4e}, v_target norm = {np.linalg.norm(Y[:, i]):.4e}")
print()

# Diagonalize T_emp
print(f"[3/3] Diagonalizing T_emp (8x8)...", flush=True)
eigvals, eigvecs = np.linalg.eig(T_emp)
print()
print(f"Eigenvalues of T_emp (sorted by |λ| descending):")
order = np.argsort(-np.abs(eigvals))
print(f"  {'idx':>3} {'Re(λ)':>14} {'Im(λ)':>14} {'|λ|':>12} {'arg(λ) °':>14}")
spec = []
for i, idx in enumerate(order):
    lam = eigvals[idx]
    abslam = abs(lam)
    arg_deg = np.angle(lam) / np.pi * 180
    spec.append({"re": float(lam.real), "im": float(lam.imag), "abs": float(abslam), "arg_deg": float(arg_deg)})
    print(f"  {i:>3} {lam.real:>+14.6e} {lam.imag:>+14.6e} {abslam:>12.6f} {arg_deg:>+14.4f}")
print()

# Compare to fitted complex eigenvalue from eps_k Model 1 (rho=0.236, theta=15.7) and Model 3 (rho=0.190, theta=41°)
print("Compare to ε_k fits:")
print(f"  Model 1 (k=2..8 single complex): rho = 0.236498, |theta| = 15.71°")
print(f"  Model 3 (k=1..8 + asymptote):    rho = 0.190144, |theta| = 40.95°")
print()
print(f"Closest match in T_emp spectrum:")
fits = [("Model 1", 0.236498, 15.71), ("Model 3", 0.190144, 40.95)]
for name, rho_fit, theta_fit in fits:
    best = None
    best_score = np.inf
    for idx in range(len(eigvals)):
        lam = eigvals[idx]
        d_rho = abs(abs(lam) - rho_fit)
        d_theta = min(abs(abs(np.angle(lam))/np.pi*180 - theta_fit), abs(abs(np.angle(lam))/np.pi*180 - theta_fit + 360))
        score = d_rho + d_theta * 0.01
        if score < best_score:
            best_score = score
            best = (idx, lam, d_rho, d_theta)
    idx, lam, d_rho, d_theta = best
    print(f"  {name}: λ_emp[{idx}] = {lam:.6f}, |λ|={abs(lam):.6f} (d_rho={d_rho:.4f}), arg={np.angle(lam)/np.pi*180:.2f}° (d_theta={d_theta:.2f}°)")
print()

# Project T_emp onto (P_+, P_-) class-resolved 2x2 subspace
# (P_+, P_-) corresponds to (P^++ summed, P^-- summed) over c
# In our 8D vec: P_+ = v[0] + v[1], P_- = v[2] + v[3]
# So projection matrix Q (2x8): Q[0] = [1,1,0,0,0,0,0,0], Q[1] = [0,0,1,1,0,0,0,0]
print("Project T_emp onto 2D (P_+, P_-) subspace:")
Q = np.array([
    [1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0],
], dtype=np.float64)
# T_2 = Q · T_emp · Q^+ (pseudo)
Qp = np.linalg.pinv(Q)
T_2 = Q @ T_emp @ Qp
print(f"  T_2 =\n{T_2}")
ev2, _ = np.linalg.eig(T_2)
print(f"  Eigenvalues of T_2: {ev2}")
print(f"  Expected from T_diag = (1/5)·[[1,1],[4,4]]: {{0, 1}}")
print()

# Save
output = {
    "v_list": [v.tolist() for v in v_list],
    "T_emp": T_emp.tolist(),
    "spectrum": spec,
    "T_2_projection": T_2.tolist(),
    "T_2_eigvals": [{"re": float(e.real), "im": float(e.imag)} for e in ev2],
    "fitted_targets": {
        "Model_1": {"rho": 0.236498, "theta_deg": 15.71},
        "Model_3": {"rho": 0.190144, "theta_deg": 40.95},
    },
    "N_MAX": N_MAX,
    "transitions_used": int(X.shape[1]),
    "residual_norm": float(np.linalg.norm(err)),
}
outpath = "C:/Collatz/T8x8_spectrum_2026_06_01.json"
with open(outpath, "w") as f:
    json.dump(output, f, indent=2)
print(f"Saved: {outpath}", flush=True)
print(f"\nTotal wall time: {time.time()-t0:.1f}s", flush=True)
